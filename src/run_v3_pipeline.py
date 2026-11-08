"""
Detached v3 pipeline (survives interactive-session restarts).
v3 = v2 objective with the belt envelope weighted 1.7x and a 130 cm3 budget,
warm-started from v2. Then: verification (2.5 and 1.25 mm), overhang audit,
printable STL with interface probes, CAD body STL, STEP via trunnion_cad.py.
"""
import json
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

sys.path.insert(0, ".")
from topo5ax import Params, Mesh
from topo5ax_v2 import TrunnionTopOptV2
from post5ax import build_field, build_cad_body_field, marching_tetrahedra, mesh_checks, write_stl, overhang_audit
import fine_fe

LOG = open("out/_v3_pipeline.log", "a", buffering=1)


def log(msg):
    LOG.write(f"{time.strftime('%H:%M:%S')} {msg}\n")


class TrunnionTopOptV3(TrunnionTopOptV2):
    def envelope_objective(self, x_phys, U):
        if self.scales is None:
            super().envelope_objective(x_phys, U)
            s_env, s_belt = self.scales
            self.scales = (s_env, 0.6 * s_belt)
        return super().envelope_objective(x_phys, U)


def main():
    out = Path("out")
    npz_path = Path("runs/v3_belt_target.npz")
    prm = Params(use_envelope=True, use_am_filter=True, target_volume_cm3=130.0,
                 beta_schedule=((0, 4.0), (12, 8.0), (25, 16.0)), max_iter=42, move=0.1)
    opt = TrunnionTopOptV3(prm, log=log)
    if npz_path.exists():
        x_phys = np.load(npz_path)["x_phys"]
        log("reusing saved v3 design")
    else:
        x0 = np.load("runs/v2_tilt_belt_robust.npz")["x"]
        x = opt.optimise(x0=x0)
        x_phys = opt.x_phys
        np.savez_compressed(npz_path, x=x, x_phys=x_phys, objective=np.array(opt.history.objective))
        log("v3 optimisation saved")

    base = Params()
    mesh = Mesh(base.nx, base.ny, base.nz, base.h)
    field = build_field(base, mesh.to3d(x_phys), voxel=1.0)

    results = {}
    for h in (2.5, 1.25):
        res = fine_fe.run_fine_verification(Params(h=h), *field, log=log)
        results[f"v3_belt_target@{h}"] = res
        (out / "verification_v3.json").write_text(json.dumps(results, indent=1, default=float))
        log(f"verified h={h}: belt {res['belt']['worst']['misalignment_arcmin']:.2f} arcmin, "
            f"worst C {res['worst_case_compliance_Nmm']:.3f}")

    audit = overhang_audit(opt, (x_phys > 0.5).astype(float))
    (out / "audit_v3.json").write_text(json.dumps(audit, indent=1, default=float))
    log(f"audit {audit}")

    verts, faces = marching_tetrahedra(*field)
    stl = out / "trunnion_upright_v3_final.stl"
    write_stl(stl, verts, faces)
    checks = mesh_checks(verts, faces)
    fld, org, spc = field

    def probe(p):
        return float(map_coordinates(fld, ((np.asarray([p]) - org) / spc).T, order=1)[0])

    ax, az = base.axle_x, base.axle_z
    probes = {
        "seat front (void)": probe([ax + 13.5, 4, az]) < 0,
        "seat back (void)": probe([ax, 36, az + 13.5]) < 0,
        "bore (void)": probe([ax, 20, az]) < 0,
        "shoulder lip (solid)": probe([ax + 12, 9, az]) > 0,
        "boss (solid)": probe([ax + 18, 5, az]) > 0,
        "bolt hole (void)": probe([10, 10, 3]) < 0,
        "pad (solid)": probe([16, 10, 3]) > 0,
        "key access (void)": probe([110, 30, 100]) < 0,
    }
    (out / "mesh_v3.json").write_text(json.dumps({"checks": checks, "probes": probes}, indent=1, default=float))
    log(f"stl {checks['watertight']} {checks['consistently_oriented']} probes {all(probes.values())}")

    body, borg, bspc = build_cad_body_field(base, mesh.to3d(x_phys), voxel=2.0)
    clamp = 0.02 * bspc
    body = np.where(np.abs(body) < clamp, np.where(body >= 0, clamp, -clamp), body)
    bv, bf = marching_tetrahedra(body, borg, bspc)
    body_stl = out / "_cad_body_v3.stl"
    write_stl(body_stl, bv, bf)
    step = out / "trunnion_upright_v3_final.step"
    cad = subprocess.run([sys.executable, "trunnion_cad.py", str(body_stl.resolve()), str(step.resolve())],
                         capture_output=True, text=True)
    log("cad: " + cad.stdout.strip()[-400:] + cad.stderr.strip()[-400:])
    with zipfile.ZipFile(out / "trunnion_upright_v3_final_step.zip", "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        zf.write(step, arcname=step.name)
    log("PIPELINE DONE")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback
        log("ERROR " + traceback.format_exc())
