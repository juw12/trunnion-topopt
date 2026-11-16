"""
Rebuild the three print files in models/ and results/mesh_v4.json from the saved designs.

Run from the repository root:  python src/rebuild_print_files.py
Nothing is optimised or solved again: each surface is the saved density field plus the exact interface
features, extracted with marching tetrahedra. Every file is then read back the way a slicer reads it
(single precision, corners welded by position) and must be closed, or the script stops with an error.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from post5ax import (  # noqa: E402
    MIN_EDGE_FRACTION,
    build_field,
    marching_tetrahedra,
    mesh_checks,
    stl_checks,
    write_stl,
)
from topo5ax import Mesh, Params  # noqa: E402

PRINT_FILES = {
    "v0_naive.npz": "trunnion_upright_v0_naive_baseline.stl",
    "v2_tilt_belt_robust.npz": "trunnion_upright_v2_final.stl",
    "v4_tilt_gate.npz": "trunnion_upright_v4_final.stl",
}
REPORTED = {"v4_tilt_gate.npz": "mesh_v4.json"}


def interface_probes(params, field, origin, spacing):
    def probe(point):
        return float(map_coordinates(field, ((np.asarray([point]) - origin) / spacing).T, order=1)[0])

    axle_x, axle_z = params.axle_x, params.axle_z
    return {
        "seat front (void)": probe([axle_x + 13.5, 4, axle_z]) < 0,
        "seat back (void)": probe([axle_x, 36, axle_z + 13.5]) < 0,
        "bore (void)": probe([axle_x, 20, axle_z]) < 0,
        "shoulder lip (solid)": probe([axle_x + 12, 9, axle_z]) > 0,
        "boss (solid)": probe([axle_x + 18, 5, axle_z]) > 0,
        "bolt hole (void)": probe([10, 10, 3]) < 0,
        "pad (solid)": probe([16, 10, 3]) > 0,
        "key access (void)": probe([110, 30, 100]) < 0,
    }


def closed_for_a_slicer(checks):
    return (
        checks["collapsed_triangles"] == 0
        and checks["zero_area_triangles"] == 0
        and checks["edges_not_shared_by_2"] == 0
        and checks["consistently_oriented"]
        and checks["volume_cm3"] > 0
    )


def main():
    params = Params()
    grid = Mesh(params.nx, params.ny, params.nz, params.h)
    for design, stl_name in PRINT_FILES.items():
        density = np.load(ROOT / "results" / "designs" / design)["x_phys"]
        field, origin, spacing = build_field(params, grid.to3d(density), voxel=1.0)
        vertices, faces = marching_tetrahedra(field, origin, spacing)
        stl_path = ROOT / "models" / stl_name
        write_stl(stl_path, vertices, faces)
        as_read = stl_checks(stl_path)
        print(
            f"{stl_name}: {as_read['triangles_in_file']} triangles, "
            f"{as_read['collapsed_triangles']} collapsed, {as_read['zero_area_triangles']} zero-area, "
            f"{as_read['edges_not_shared_by_2']} edges not shared by 2, "
            f"volume {as_read['volume_cm3']:.4f} cm3"
        )
        if not closed_for_a_slicer(as_read):
            raise SystemExit(f"{stl_name} is not closed in single precision")
        if design in REPORTED:
            report = {
                "checks": mesh_checks(vertices, faces),
                "probes": interface_probes(params, field, origin, spacing),
                "as_read_from_the_stl_file": as_read,
                "min_edge_fraction": MIN_EDGE_FRACTION,
            }
            (ROOT / "results" / REPORTED[design]).write_text(
                json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8"
            )


if __name__ == "__main__":
    main()
