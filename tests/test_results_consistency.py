"""The committed results must follow from the committed designs and agree with each other.

What is recomputed here from `results/designs/*.npz`: the overhang and connectivity
audit, the printable surface and its interface probes. What is only cross-checked
(not recomputed, the solve has 1.16 M unknowns): the finite-element numbers in
`results/verification_*.json` against the README table and the results page data.
"""
import json
import struct
import sys
import types
import unittest
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from post5ax import build_field, marching_tetrahedra, mesh_checks, overhang_audit  # noqa: E402
from topo5ax import Mesh, Params, build_regions  # noqa: E402

RESULTS = ROOT / "results"
DENSITY_G_PER_CM3 = 1.25
TILT_LIMIT_ARCMIN = 5.0


def load_json(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def load_density(name):
    return np.load(RESULTS / "designs" / name)["x_phys"]


def fine_and_coarse():
    """{design: {"fine": block at 1.25 mm, "coarse": block at 2.5 mm}} for v0, v2 and v4."""
    corrected = load_json("verification_corrected.json")
    latest = load_json("verification_v4.json")
    return {
        "v0": {"fine": corrected["v0_naive@1.25"], "coarse": corrected["v0_naive@2.5"]},
        "v2": {
            "fine": corrected["v2_tilt_belt_robust@1.25"],
            "coarse": corrected["v2_tilt_belt_robust@2.5"],
        },
        "v4": {"fine": latest["v4_tilt_gate@1.25"], "coarse": latest["v4_tilt_gate@2.5"]},
    }


def peak_stress(block):
    return max(block["vm_belt_max_MPa"], block["vm_envelope_max_MPa"])


class DesignGridTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = Params()
        cls.mesh = Mesh(cls.params.nx, cls.params.ny, cls.params.nz, cls.params.h)
        cls.solid, cls.void = build_regions(cls.params, cls.mesh)
        cls.stand_in = types.SimpleNamespace(mesh=cls.mesh, solid=cls.solid, void=cls.void)

    def test_grid_matches_the_readme(self):
        self.assertEqual((self.params.nx, self.params.ny, self.params.nz), (48, 16, 60))
        self.assertEqual(self.mesh.n_el, 46080)

    def test_keep_in_and_keep_out_regions_do_not_overlap_and_are_mirror_symmetric(self):
        self.assertFalse(np.any(self.solid & self.void))
        for region in (self.solid, self.void):
            grid = self.mesh.to3d(region)
            self.assertTrue(np.array_equal(grid, grid[::-1, :, :]))

    def test_bore_is_open_and_bosses_and_pads_are_solid(self):
        centres = self.mesh.el_centers
        prm = self.params

        def nearest(point):
            return int(np.argmin(np.linalg.norm(centres - np.asarray(point), axis=1)))

        self.assertTrue(self.void[nearest([prm.axle_x, 20.0, prm.axle_z])])            # through bore
        self.assertTrue(self.void[nearest([prm.axle_x + 11.0, 3.0, prm.axle_z])])      # bearing seat
        self.assertTrue(self.solid[nearest([prm.axle_x + 18.0, 3.0, prm.axle_z])])     # boss ring
        self.assertTrue(self.solid[nearest([16.0, 20.0, 3.0])])                        # base pad
        self.assertTrue(self.void[nearest([10.0, 10.0, 3.0])])                         # M5 hole
        self.assertTrue(self.void[nearest([10.0, 10.0, 60.0])])                        # key access

    def test_every_saved_design_has_the_grid_size_and_density_range(self):
        for path in sorted((RESULTS / "designs").glob("*.npz")):
            density = np.load(path)["x_phys"]
            self.assertEqual(density.shape, (self.mesh.n_el,), path.name)
            self.assertGreaterEqual(density.min(), -1e-9, path.name)
            self.assertLessEqual(density.max(), 1.01, path.name)

    def test_v4_respects_the_interface_regions_and_is_symmetric(self):
        density = load_density("v4_tilt_gate.npz")
        self.assertTrue(np.all(density[self.solid] == 1.0))
        self.assertTrue(np.all(density[self.void] == 0.0))
        grid = self.mesh.to3d(density)
        self.assertLess(np.abs(grid - grid[::-1, :, :]).max(), 1e-9)

    def test_v4_uses_its_145_cm3_material_budget(self):
        density = load_density("v4_tilt_gate.npz")
        volume = density.sum() * self.params.h ** 3 / 1000.0
        self.assertAlmostEqual(volume, 145.0, delta=0.05)

    def test_overhang_audits_are_reproduced_from_the_saved_designs(self):
        for version, design in (("v3", "v3_belt_target.npz"), ("v4", "v4_tilt_gate.npz")):
            binary = (load_density(design) > 0.5).astype(float)
            recomputed = overhang_audit(self.stand_in, binary)
            self.assertEqual(recomputed, load_json(f"audit_{version}.json"), version)

    def test_v4_has_no_unsupported_overhang_and_is_one_body(self):
        audit = load_json("audit_v4.json")
        self.assertEqual(audit["unsupported_elements"], 0)
        self.assertEqual(audit["connected_components"], 1)
        self.assertGreater(audit["downfacing_elements"], 0)


class PrintableSurfaceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = Params()
        mesh = Mesh(cls.params.nx, cls.params.ny, cls.params.nz, cls.params.h)
        density = load_density("v4_tilt_gate.npz")
        cls.field, cls.origin, cls.spacing = build_field(cls.params, mesh.to3d(density), voxel=1.0)
        cls.vertices, cls.faces = marching_tetrahedra(cls.field, cls.origin, cls.spacing)
        cls.checks = mesh_checks(cls.vertices, cls.faces)
        cls.committed = load_json("mesh_v4.json")

    def probe(self, point):
        index = ((np.asarray([point]) - self.origin) / self.spacing).T
        return float(map_coordinates(self.field, index, order=1)[0])

    def test_rebuilt_surface_is_closed_and_oriented(self):
        self.assertTrue(self.checks["watertight"])
        self.assertTrue(self.checks["consistently_oriented"])
        self.assertEqual(self.checks["edges_not_shared_by_2"], 0)

    def test_rebuilt_surface_matches_the_committed_mesh_report(self):
        committed = self.committed["checks"]
        self.assertAlmostEqual(
            self.checks["volume_cm3"], committed["volume_cm3"], delta=1e-3 * committed["volume_cm3"]
        )
        self.assertAlmostEqual(
            self.checks["area_cm2"], committed["area_cm2"], delta=1e-3 * committed["area_cm2"]
        )
        self.assertAlmostEqual(
            self.checks["triangles"], committed["triangles"], delta=0.01 * committed["triangles"]
        )
        for rebuilt, reported in zip(self.checks["bbox_mm"], committed["bbox_mm"]):
            self.assertAlmostEqual(rebuilt, reported, delta=0.05)

    def test_part_fits_the_design_envelope_and_sits_flat_on_the_bed(self):
        length, depth, height = self.checks["bbox_mm"]
        self.assertAlmostEqual(length, self.params.length_x, delta=0.01)
        self.assertAlmostEqual(depth, self.params.depth_y, delta=0.01)
        self.assertLessEqual(height, self.params.height_z + 0.01)

    def test_interface_probes_are_reproduced(self):
        axle_x, axle_z = self.params.axle_x, self.params.axle_z
        probes = {
            "seat front (void)": self.probe([axle_x + 13.5, 4, axle_z]) < 0,
            "seat back (void)": self.probe([axle_x, 36, axle_z + 13.5]) < 0,
            "bore (void)": self.probe([axle_x, 20, axle_z]) < 0,
            "shoulder lip (solid)": self.probe([axle_x + 12, 9, axle_z]) > 0,
            "boss (solid)": self.probe([axle_x + 18, 5, axle_z]) > 0,
            "bolt hole (void)": self.probe([10, 10, 3]) < 0,
            "pad (solid)": self.probe([16, 10, 3]) > 0,
            "key access (void)": self.probe([110, 30, 100]) < 0,
        }
        self.assertEqual(probes, self.committed["probes"])
        self.assertTrue(all(probes.values()))

    def test_bearing_seat_and_bolt_hole_sizes_in_the_field(self):
        # Where an exact feature cuts kept-solid material the field is zero on the wall.
        # Each probe sits on a grid point inside a boss or a pad, so there is no rounding.
        axle_x, axle_z = self.params.axle_x, self.params.axle_z
        seat = self.params.seat_radius
        shoulder = self.params.shoulder_radius
        bolt_x, bolt_y = self.params.bolt_xy[0]
        self.assertAlmostEqual(self.probe([axle_x + seat, 4, axle_z]), 0.0, delta=1e-9)       # 28 mm seat, front
        self.assertAlmostEqual(self.probe([axle_x, 36, axle_z - seat]), 0.0, delta=1e-9)      # 28 mm seat, back
        self.assertAlmostEqual(self.probe([axle_x + shoulder, 9, axle_z]), 0.0, delta=1e-9)   # 16 mm bore, front lip
        self.assertAlmostEqual(self.probe([axle_x, 31, axle_z - shoulder]), 0.0, delta=1e-9)  # 16 mm bore, back lip
        self.assertAlmostEqual(
            self.probe([bolt_x + self.params.bolt_hole_radius, bolt_y, 3]), 0.0, delta=0.05   # M5 clearance
        )
        # one millimetre either side of each wall: open inside, material outside
        self.assertLess(self.probe([axle_x + seat - 1, 4, axle_z]), 0.0)
        self.assertGreater(self.probe([axle_x + seat + 1, 4, axle_z]), 0.0)
        self.assertLess(self.probe([axle_x + shoulder - 1, 9, axle_z]), 0.0)
        self.assertGreater(self.probe([axle_x + shoulder + 1, 9, axle_z]), 0.0)

    def test_shipped_print_file_is_this_surface(self):
        # models/trunnion_upright_v4_final.stl must be the surface rebuilt from the saved design,
        # triangle for triangle, in the single precision an STL file stores.
        data = (ROOT / "models" / "trunnion_upright_v4_final.stl").read_bytes()
        count = struct.unpack("<I", data[80:84])[0]
        self.assertEqual(count, self.committed["checks"]["triangles"])
        self.assertEqual(count, len(self.faces))
        self.assertEqual(len(data), 84 + 50 * count)
        record = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        shipped = np.frombuffer(data[84:], dtype=record)["v"]
        self.assertTrue(np.array_equal(shipped, self.vertices[self.faces].astype("<f4")))


class ShippedStlTest(unittest.TestCase):
    """What a slicer reads: the file's own single-precision coordinates, welded by position."""

    @classmethod
    def setUpClass(cls):
        data = (ROOT / "models" / "trunnion_upright_v4_final.stl").read_bytes()
        record = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        corners = np.frombuffer(data[84:], dtype=record)["v"].reshape(-1, 3)
        vertices, index = np.unique(corners, axis=0, return_inverse=True)
        faces = index.reshape(-1, 3)
        cls.collapsed = (
            (faces[:, 0] == faces[:, 1]) | (faces[:, 1] == faces[:, 2]) | (faces[:, 0] == faces[:, 2])
        )
        cls.checks = mesh_checks(vertices.astype(float), faces[~cls.collapsed])
        cls.committed = load_json("mesh_v4.json")["checks"]

    def test_volume_area_and_size_match_the_mesh_report(self):
        self.assertAlmostEqual(
            self.checks["volume_cm3"], self.committed["volume_cm3"], delta=1e-6 * self.committed["volume_cm3"]
        )
        self.assertAlmostEqual(
            self.checks["area_cm2"], self.committed["area_cm2"], delta=1e-6 * self.committed["area_cm2"]
        )
        for welded, reported in zip(self.checks["bbox_mm"], self.committed["bbox_mm"]):
            self.assertAlmostEqual(welded, reported, delta=1e-4)
        self.assertGreater(self.checks["volume_cm3"], 0.0)          # normals point outward

    def test_welding_leaves_only_a_few_open_edges(self):
        # Known limit, recorded so it cannot grow unnoticed. The mesh report says "watertight"
        # for the double-precision mesh. In the file's single precision about a fifth of the
        # triangles are slivers that collapse to a line, and welding leaves 83 edges that are
        # not shared by exactly two triangles. A slicer has to repair those.
        self.assertTrue(self.committed["watertight"])
        self.assertLess(self.collapsed.mean(), 0.25)
        edges_checked = 3 * int((~self.collapsed).sum()) // 2
        self.assertLessEqual(self.checks["edges_not_shared_by_2"], 83)
        self.assertLess(self.checks["edges_not_shared_by_2"] / edges_checked, 2e-4)


class ReportedNumbersTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.blocks = fine_and_coarse()
        cls.readme = (ROOT / "README.md").read_text(encoding="utf-8")
        cls.table = {}
        for line in cls.readme.splitlines():
            cells = [cell.strip().strip("*") for cell in line.strip().strip("|").split("|")]
            if line.startswith("|") and len(cells) == 4:
                cls.table[cells[0]] = dict(zip(("v0", "v2", "v4"), cells[1:]))

    def row(self, start):
        matches = [cells for label, cells in self.table.items() if label.startswith(start)]
        self.assertEqual(len(matches), 1, start)
        return matches[0]

    def test_readme_table_matches_the_fine_mesh_results(self):
        volume = self.row("Printed volume")
        compliance = self.row("Worst-case compliance")
        tilt = self.row("Axle tilt")
        deflection = self.row("Worst-case bearing-seat deflection")
        stress = self.row("Peak von Mises")
        for design, pair in self.blocks.items():
            fine = pair["fine"]
            mass = fine["volume_cm3"] * DENSITY_G_PER_CM3
            self.assertEqual(volume[design], f"{fine['volume_cm3']:.1f} cm³ / {mass:.0f} g", design)
            self.assertEqual(
                compliance[design], f"{fine['worst_case_compliance_Nmm']:.2f} N·mm", design
            )
            self.assertEqual(
                tilt[design], f"{fine['belt']['worst']['misalignment_arcmin']:.1f} arcmin", design
            )
            self.assertEqual(
                deflection[design],
                f"{fine['seat_deflection_um']['worst-case envelope']:.0f} µm",
                design,
            )
            self.assertEqual(stress[design], f"{peak_stress(fine):.1f} MPa", design)

    def test_readme_overhang_row_matches_the_audit_and_results_page_data(self):
        overhang = self.row("Unsupported overhang")
        designs = load_json("ui_summary.json")["designs"]
        self.assertEqual(overhang["v0"], f"{designs['v0_naive']['overhang_pct']:.1f} %")
        self.assertEqual(overhang["v2"], f"{designs['v2_tilt_belt_robust']['overhang_pct']:.1f} %")
        audit = load_json("audit_v4.json")
        self.assertEqual(
            overhang["v4"], f"{100 * audit['unsupported_fraction_of_downfacing']:.1f} %"
        )

    def test_results_page_data_matches_the_verification_files(self):
        designs = load_json("ui_summary.json")["designs"]
        names = {"v0": "v0_naive", "v2": "v2_tilt_belt_robust", "v4": "v4_tilt_gate"}
        for design, pair in self.blocks.items():
            fine = pair["fine"]
            shown = designs[names[design]]
            self.assertAlmostEqual(shown["volume_cm3"], fine["volume_cm3"], places=9)
            self.assertAlmostEqual(
                shown["worst_compliance"], fine["worst_case_compliance_Nmm"], places=9
            )
            self.assertAlmostEqual(
                shown["belt_arcmin"], fine["belt"]["worst"]["misalignment_arcmin"], places=9
            )
            self.assertAlmostEqual(shown["vm_max"], peak_stress(fine), places=9)
            self.assertAlmostEqual(
                shown["mass_g"], fine["volume_cm3"] * DENSITY_G_PER_CM3, places=6
            )

    def test_only_v4_is_marked_as_passing_every_gate(self):
        summary = load_json("ui_summary.json")
        self.assertEqual(summary["recommended"], "v4_tilt_gate")
        passing = [name for name, design in summary["designs"].items() if design["all_pass"]]
        self.assertEqual(passing, ["v4_tilt_gate"])

    def test_v4_meets_the_tilt_limit_on_both_meshes_and_v2_does_not(self):
        for mesh_size in ("fine", "coarse"):
            tilt = self.blocks["v4"][mesh_size]["belt"]["worst"]["misalignment_arcmin"]
            self.assertLess(tilt, TILT_LIMIT_ARCMIN, mesh_size)
        self.assertGreater(
            self.blocks["v2"]["fine"]["belt"]["worst"]["misalignment_arcmin"], TILT_LIMIT_ARCMIN
        )

    def test_v2_belt_pull_limit_quoted_in_the_readme(self):
        # linear FE: tilt scales with belt force, so the 5 arcmin limit is reached at
        # 100 N * 5 / tilt_at_100_N
        belt = self.blocks["v2"]["fine"]["belt"]
        limit = belt["force_N"] * TILT_LIMIT_ARCMIN / belt["worst"]["misalignment_arcmin"]
        self.assertIn(f"below {limit:.0f} N", self.readme)

    def test_stiffness_results_agree_between_coarse_and_fine_meshes_for_v2_and_v4(self):
        for design in ("v2", "v4"):
            fine = self.blocks[design]["fine"]
            coarse = self.blocks[design]["coarse"]
            pairs = {
                "compliance": (
                    coarse["worst_case_compliance_Nmm"],
                    fine["worst_case_compliance_Nmm"],
                ),
                "seat deflection": (
                    coarse["seat_deflection_um"]["worst-case envelope"],
                    fine["seat_deflection_um"]["worst-case envelope"],
                ),
                "belt tilt": (
                    coarse["belt"]["worst"]["misalignment_arcmin"],
                    fine["belt"]["worst"]["misalignment_arcmin"],
                ),
            }
            for quantity, (coarse_value, fine_value) in pairs.items():
                difference = abs(fine_value - coarse_value) / coarse_value
                self.assertLess(difference, 0.035, f"{design} {quantity}")

    def test_peak_stress_is_not_mesh_converged_and_the_readme_says_so(self):
        # the fine mesh finds 31 to 35 % more peak stress, so the README must not call it converged
        for design in ("v0", "v2", "v4"):
            rise = peak_stress(self.blocks[design]["fine"]) / peak_stress(self.blocks[design]["coarse"]) - 1.0
            self.assertGreater(rise, 0.30, design)
            self.assertLess(rise, 0.36, design)
        self.assertIn("Peak stress is not mesh-converged", self.readme)
        self.assertNotIn("fine meshes agree within", self.readme)

    def test_envelope_eigenvalues_are_sorted_and_the_largest_is_the_worst_case(self):
        for design, pair in self.blocks.items():
            for mesh_size, block in pair.items():
                eigenvalues = block["envelope_eig_Nmm"]
                self.assertEqual(eigenvalues, sorted(eigenvalues), f"{design} {mesh_size}")
                self.assertAlmostEqual(
                    eigenvalues[-1], block["worst_case_compliance_Nmm"], places=12
                )

    def test_belt_reactions_follow_from_statics(self):
        # pulley overhung outboard of the outer bearing: moments about the inner seat
        prm = Params()
        belt = self.blocks["v4"]["fine"]["belt"]
        span = prm.depth_y - prm.seat_depth          # distance between bearing centres
        outer = belt["force_N"] * (span + belt["pulley_offset_mm"] + prm.seat_depth / 2) / span
        self.assertAlmostEqual(belt["outer_reaction_N"] + belt["inner_reaction_N"], belt["force_N"])
        self.assertAlmostEqual(belt["outer_reaction_N"], outer, places=6)


if __name__ == "__main__":
    unittest.main()
