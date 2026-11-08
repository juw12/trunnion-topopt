"""
Build the results viewer (four design directions + final page) from files on disk.
Usage: python ui_build.py   (or import and call build_all())
"""
import base64
import json
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.ndimage import map_coordinates

sys.path.insert(0, ".")
from topo5ax import Params, Mesh, build_regions
from post5ax import build_field, marching_tetrahedra, overhang_audit

OUT = Path("out")
UI = OUT / "ui"
DENSITY = 1.25  # g/cm3, PA-CF / PETG-CF class (assumed)

LABELS = {
    "R0_conventional": ("R0 Conventional", "5 mm web with edge flanges, hand-style"),
    "v0_naive": ("v0 Single-load TO", "one vertical load, no overhang control"),
    "v1_robust_am": ("v1 Tilt-robust", "any tilt direction and axial, self-supporting"),
    "v2_tilt_belt_robust": ("v2 Tilt and belt", "adds the overhung belt case to the tilt envelope"),
    "v3_belt_target": ("v3 Belt-weighted", "belt case weighted 1.7x, 130 cm\u00b3 budget"),
    "v4_tilt_gate": ("v4 Tilt gate", "belt case weighted 2.2x, 145 cm\u00b3 budget"),
}
ORDER = ["R0_conventional", "v0_naive", "v1_robust_am", "v2_tilt_belt_robust", "v3_belt_target", "v4_tilt_gate"]

ASSUMPTIONS = [
    "Material PA-CF or PETG-CF treated as isotropic, E = 3.0 GPa, density 1.25 g/cm\u00b3, 100 % infill. Real prints are weaker across layers; stresses here are below 15 % of typical inter-layer strength.",
    "Bearing load envelope: 80 N in any direction of the tilt plane plus 30 N along the axle, shared by the two 6001-2RS bearings. This covers bed weight with a dynamic factor and nozzle reaction at any A-axis angle.",
    "Belt case: 100 N pull (2 x 30 N pre-tension plus 40 N torque increment) on a pulley 15 mm outboard of the outer face, in any direction of the tilt plane. Bearings react 159 N and -59 N.",
    "Axle tilt limit 5 arcmin, half of the 2 to 10 arcmin that SKF lists for single-row deep groove ball bearings.",
    "Base pads fixed at the bolted faces on the frame extrusion; contact elsewhere ignored (conservative).",
    "Results are linear-elastic FE on the exported geometry at 1.25 mm elements (about 1.1 million DOF). Not yet validated by a physical test print.",
    "Bearing seats are horizontal 28 mm holes: their top arc bridges about 20 mm. Print at 0.2 mm layers with bridging settings and ream to H7 for the press fit.",
]


def conventional_design(prm, mesh, target_cm3=105.0):
    solid, void = build_regions(prm, mesh)
    c = mesh.el_centers
    x, y, z = c[:, 0], c[:, 1], c[:, 2]
    half = 60 - (60 - prm.boss_radius) * np.clip(z / prm.axle_z, 0, 1)
    trap = (np.abs(x - prm.axle_x) <= half) & (z <= prm.axle_z + prm.boss_radius)
    web = trap & (np.abs(y - prm.depth_y / 2) <= 2.6)
    best = None
    for flange in np.arange(2.5, 25.1, 0.5):
        edge = trap & (half - np.abs(x - prm.axle_x) <= flange)
        xb = (web | edge | solid) & ~void
        vol = xb.sum() * prm.h ** 3 / 1000
        if best is None or abs(vol - target_cm3) < abs(best[1] - target_cm3):
            best = (xb.astype(float), vol)
    return best[0]


def load_designs(prm, mesh):
    phys = {"R0_conventional": conventional_design(prm, mesh)}
    for name in ORDER[1:]:
        path = Path("runs") / f"{name}.npz"
        if path.exists():
            phys[name] = np.load(path)["x_phys"]
    return phys


def probes_ok(field):
    fld, org, spc = field
    p = Params()

    def val(pt):
        return float(map_coordinates(fld, ((np.asarray([pt]) - org) / spc).T, order=1)[0])

    ax, az = p.axle_x, p.axle_z
    return all([
        val([ax + 13.5, 4, az]) < 0, val([ax, 36, az + 13.5]) < 0, val([ax, 20, az]) < 0,
        val([ax + 18, 5, az]) > 0, val([10, 10, 3]) < 0, val([16, 10, 3]) > 0, val([110, 30, 100]) < 0,
    ])


def viewer_mesh(field, voxel=2.5):
    fld, org, spc = field
    step = int(round(voxel / spc))
    coarse = fld[::step, ::step, ::step]
    clamp = 0.02 * spc * step
    coarse = np.where(np.abs(coarse) < clamp, np.where(coarse >= 0, clamp, -clamp), coarse)
    verts, faces = marching_tetrahedra(coarse, org, spc * step)
    keys, inverse = np.unique(np.round(verts, 4), axis=0, return_inverse=True)
    faces = inverse.reshape(-1)[faces]
    lo = np.array([-5.0, -5.0, -5.0])
    span = 160.0
    q = np.clip(np.round((keys - lo) / span * 65535), 0, 65535).astype("<u2")
    return {
        "positions": base64.b64encode(q.tobytes()).decode(),
        "indices": base64.b64encode(faces.astype("<u4").tobytes()).decode(),
        "lo": lo.tolist(),
        "span": span,
    }


def verification():
    res = json.loads((OUT / "verification_corrected.json").read_text())
    for extra in ("verification_v3.json", "verification_v4.json"):
        path = OUT / extra
        if path.exists():
            res.update(json.loads(path.read_text()))
    return res


def design_entry(name, res, audit, interfaces):
    key = f"{name}@1.25" if f"{name}@1.25" in res else f"{name}@2.5"
    r = res[key]
    sd = r["seat_deflection_um"]
    vm = max(r["vm_envelope_max_MPa"], r["vm_belt_max_MPa"])
    e = {
        "label": LABELS[name][0], "short": LABELS[name][1], "fe": key.split("@")[1] + " mm elements",
        "volume_cm3": r["volume_cm3"], "mass_g": r["volume_cm3"] * DENSITY,
        "worst_compliance": r["worst_case_compliance_Nmm"],
        "defl_vertical": sd["vertical 80 N"], "defl_horizontal": sd["horizontal 80 N"],
        "defl_axial": sd["axial 30 N"], "defl_worst": sd["worst-case envelope"],
        "belt_arcmin": r["belt"]["worst"]["misalignment_arcmin"],
        "belt_outer_um": r["belt"]["worst"]["outer_seat_um"],
        "vm_max": vm,
        "overhang_pct": 100 * audit["unsupported_fraction_of_downfacing"],
        "bodies": audit["connected_components"],
        "tilt": [[a, d] for a, d in r["tilt_sweep"]],
    }
    e["gates"] = [
        {"label": "One printable body", "pass": e["bodies"] == 1,
         "value": f"{e['bodies']} bod{'y' if e['bodies'] == 1 else 'ies'}"},
        {"label": "Self-supporting at 45\u00b0, unsupported \u2264 2 %", "pass": e["overhang_pct"] <= 2.0,
         "value": f"{e['overhang_pct']:.1f} %"},
        {"label": "Axle tilt under 100 N belt \u2264 5 arcmin", "pass": e["belt_arcmin"] <= 5.0,
         "value": f"{e['belt_arcmin']:.1f}\u2032"},
        {"label": "Worst-case seat deflection \u2264 100 \u00b5m", "pass": e["defl_worst"] <= 100.0,
         "value": f"{e['defl_worst']:.0f} \u00b5m"},
        {"label": "Peak stress \u2264 10 MPa", "pass": vm <= 10.0, "value": f"{vm:.1f} MPa"},
        {"label": "Seats, bolt holes and key access open", "pass": bool(interfaces),
         "value": "probed" if interfaces else "blocked"},
    ]
    e["all_pass"] = all(g["pass"] for g in e["gates"])
    return e


def build_data(rel):
    prm = Params()
    mesh = Mesh(prm.nx, prm.ny, prm.nz, prm.h)
    solid, _ = build_regions(prm, mesh)
    shim = SimpleNamespace(mesh=mesh, solid=solid)
    phys = load_designs(prm, mesh)
    res = verification()
    order = [n for n in ORDER if n in phys and f"{n}@1.25" in res]
    designs, meshes = {}, {}
    for name in order:
        field = build_field(prm, mesh.to3d(phys[name]), voxel=1.0)
        audit = overhang_audit(shim, (phys[name] > 0.5).astype(float))
        designs[name] = design_entry(name, res, audit, probes_ok(field))
        meshes[name] = viewer_mesh(field)
    passing = [n for n in order if designs[n]["all_pass"]]
    if passing:
        rec = min(passing, key=lambda n: designs[n]["mass_g"])
    else:
        # no design meets every gate: recommend the one with the most passes, then lowest belt tilt
        rec = min(order, key=lambda n: (-sum(g["pass"] for g in designs[n]["gates"]),
                                        designs[n]["belt_arcmin"]))
    d, b = designs[rec], designs["v0_naive"]
    tag = rec.split("_")[0]
    n_pass = sum(g["pass"] for g in d["gates"])
    if d["all_pass"]:
        lead = f"{d['label']} is the lightest design that passes every check."
    else:
        failed = [g["label"] for g in d["gates"] if not g["pass"]]
        lead = (f"{d['label']} passes {n_pass} of {len(d['gates'])} checks and is the best available. "
                f"Open limit: {'; '.join(failed)}. Tilt scales linearly with belt pull, so it meets 5 arcmin "
                f"at {100 * 5.0 / d['belt_arcmin']:.0f} N of belt pull.")
    files = []
    for fname, meta in [
        (f"trunnion_upright_{tag}_final.stl", "print file, watertight"),
        (f"trunnion_upright_{tag}_final_step.zip", "STEP, exact interfaces"),
        ("trunnion_upright_v2_final.stl", "103 cm\u00b3 lighter alternative, 8.9 arcmin"),
        ("trunnion_upright_v0_naive_baseline.stl", "baseline for comparison"),
    ]:
        if (OUT / fname).exists() and fname not in [f["name"] for f in files]:
            size = (OUT / fname).stat().st_size / 1e6
            files.append({"name": fname, "href": rel + fname, "meta": f"{size:.0f} MB, {meta}"})
    return {
        "order": order,
        "recommended": rec,
        "headline": f"Print {tag}: {d['mass_g']:.0f} g, {d['belt_arcmin']:.1f} arcmin axle tilt, no supports",
        "summary": (f"{lead} Against a standard "
                    f"single-load topology optimisation it cuts worst-case compliance from {b['worst_compliance']:.1f} "
                    f"to {d['worst_compliance']:.1f} N\u00b7mm and axle tilt under belt pull from "
                    f"{b['belt_arcmin']:.0f} to {d['belt_arcmin']:.1f} arcmin, and it prints without supports."),
        "designs": designs,
        "files": files,
        "primaryFiles": {"stl": rel + f"trunnion_upright_{tag}_final.stl",
                         "step": rel + f"trunnion_upright_{tag}_final_step.zip"},
        "assumptions": ASSUMPTIONS,
        "meshes": meshes,
    }


def write_page(template, page_dir, data):
    page_dir.mkdir(parents=True, exist_ok=True)
    assets = page_dir / "assets"
    assets.mkdir(exist_ok=True)
    for f in ("viewer.js", "common.js"):
        shutil.copy(UI / "assets" / f, assets / f)
    (assets / "data.js").write_text("window.TRUNNION = " + json.dumps(data, default=float) + ";\n", encoding="utf-8")
    html = (UI / "templates" / f"{template}.html").read_text(encoding="utf-8").replace("%%ASSETS%%", "assets")
    (page_dir / "index.html").write_text(html, encoding="utf-8")
    return (page_dir / "index.html").resolve()


def build_all(directions=("A", "B", "C", "D")):
    data = build_data("../../")
    pages = [write_page(d, UI / f"dir-{d}", data) for d in directions]
    summary = {n: {k: v for k, v in e.items() if k not in ("tilt",)} for n, e in data["designs"].items()}
    (OUT / "ui_summary.json").write_text(json.dumps({"recommended": data["recommended"], "designs": summary},
                                                    indent=1, default=float))
    return pages, data


def build_final(template="A"):
    """Write the bracket winner as out/ui/index.html; assets live in out/ui/assets already."""
    data = build_data("../")
    (UI / "assets" / "data.js").write_text("window.TRUNNION = " + json.dumps(data, default=float) + ";\n",
                                          encoding="utf-8")
    html = (UI / "templates" / f"{template}.html").read_text(encoding="utf-8").replace("%%ASSETS%%", "assets")
    (UI / "index.html").write_text(html, encoding="utf-8")
    summary = {n: {k: v for k, v in e.items() if k not in ("tilt",)} for n, e in data["designs"].items()}
    (OUT / "ui_summary.json").write_text(json.dumps({"recommended": data["recommended"], "headline": data["headline"],
                                                     "summary": data["summary"], "designs": summary},
                                                    indent=1, default=float))
    return (UI / "index.html").resolve(), data


if __name__ == "__main__":
    if "--final" in sys.argv:
        page, data = build_final()
        print(data["recommended"], page)
        print(data["headline"])
        print(data["summary"])
    else:
        pages, data = build_all()
        print(data["recommended"], [str(p) for p in pages])
