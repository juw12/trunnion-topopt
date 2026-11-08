"""
Fine-mesh verification of the exported (printable) geometry.

The implicit field that produced the STL (optimised density + exact interface
features) is sampled at the centres of a finer hex mesh, and the resulting
black/white part is analysed with linear elasticity:

* bearing-load envelope (80 N any direction in the tilt plane, 30 N axial),
  shared equally by the two bearing seats,
* overhung belt pull on the axle pulley (bearings react unequally and the axle
  tilts), swept over every belt direction in the tilt plane,
* von Mises stress for the governing cases.
"""
from __future__ import annotations

import math
import time

import numpy as np
import scipy.sparse as sp
from scipy.ndimage import map_coordinates

from topo5ax import (
    Mesh,
    MultigridSolver,
    Params,
    element_stiffness,
    elasticity_matrix,
    strain_displacement,
    build_supports_and_loads,
)


BELT_FORCE_N = 100.0          # 2 x 30 N pre-tension + 40 N torque increment (assumed)
PULLEY_OFFSET_MM = 15.0       # pulley centre outboard of the outer face (assumed)


def sample_design(prm: Params, mesh: Mesh, field, origin, spacing):
    idx = ((mesh.el_centers - origin) / spacing).T
    values = map_coordinates(field, idx, order=1, mode="nearest")
    return (values > 0).astype(float)


def assemble_solid(mesh: Mesh, KE, x_bin, E0, free):
    solid = np.nonzero(x_bin > 0.5)[0]
    edof = mesh.edof[solid].astype(np.int64)
    rows = np.repeat(edof, 24, axis=1).ravel()
    cols = np.tile(edof, (1, 24)).ravel()
    vals = np.tile(KE.ravel() * E0, len(solid))
    K = sp.coo_matrix((vals, (rows, cols)), shape=(mesh.n_dofs, mesh.n_dofs)).tocsr()
    # tiny grounding spring keeps nodes that only touch void elements non-singular
    k_ground = 1e-6 * E0 * mesh.h
    K = K + sp.identity(mesh.n_dofs, format="csr") * k_ground
    return K[free][:, free].tocsr()


def run_fine_verification(prm: Params, field, origin, spacing, log=print):
    t0 = time.time()
    mesh = Mesh(prm.nx, prm.ny, prm.nz, prm.h)
    KE = element_stiffness(prm.h, prm.nu)
    fixed, free, _, _, _ = build_supports_and_loads(prm, mesh)
    x_bin = sample_design(prm, mesh, field, origin, spacing)
    # bearing contact = solid-attached nodes on the seat surface (innermost solid layer)
    attached = np.zeros(mesh.n_nodes, dtype=bool)
    attached[mesh.edof[x_bin > 0.5][:, ::3] // 3] = True
    nxyz = mesh.node_xyz
    r_axle = np.hypot(nxyz[:, 0] - prm.axle_x, nxyz[:, 2] - prm.axle_z)
    in_seat_zone = (nxyz[:, 1] <= prm.seat_depth + 1e-9) | (nxyz[:, 1] >= prm.depth_y - prm.seat_depth - 1e-9)
    seat = np.nonzero(attached & in_seat_zone & (r_axle >= prm.seat_radius - 0.75 * prm.h)
                      & (r_axle <= prm.seat_radius + 0.75 * prm.h))[0]
    log(f"  sampled design: {int(x_bin.sum())} solid elements, "
        f"{x_bin.sum() * prm.h ** 3 / 1000:.1f} cm3 ({time.time() - t0:.0f}s)")

    y_nodes = mesh.node_xyz[seat, 1]
    log(f"  bearing contact nodes: {len(seat)}")
    inner = seat[y_nodes <= prm.seat_depth + 1e-9]
    outer = seat[y_nodes >= prm.depth_y - prm.seat_depth - 1e-9]

    def unit_load(nodes, direction):
        f = np.zeros(mesh.n_dofs)
        f[3 * nodes + direction] = 1.0 / len(nodes)
        return f

    load_cases = {
        "inner_x": unit_load(inner, 0),
        "inner_z": unit_load(inner, 2),
        "outer_x": unit_load(outer, 0),
        "outer_z": unit_load(outer, 2),
        "both_y": 0.5 * (unit_load(inner, 1) + unit_load(outer, 1)),
    }

    K = assemble_solid(mesh, KE, x_bin, prm.E0, free)
    log(f"  assembled {K.shape[0]} free DOF, nnz {K.nnz / 1e6:.1f} M ({time.time() - t0:.0f}s)")
    mg = MultigridSolver(mesh, free, levels=4)
    mg.setup(K)
    log(f"  multigrid levels {mg.levels} ({time.time() - t0:.0f}s)")

    U = {}
    for name, f in load_cases.items():
        u = np.zeros(mesh.n_dofs)
        sol, its = mg.solve(f[free], tol=1e-7, maxit=400)
        u[free] = sol
        U[name] = u
        log(f"  solved {name}: {its} PCG its ({time.time() - t0:.0f}s)")

    def seat_disp(u, nodes):
        return u.reshape(-1, 3)[nodes].mean(axis=0)

    def combine(fx_in, fz_in, fx_out, fz_out, fy=0.0):
        """Displacement field for given per-seat forces (N); linear superposition."""
        return (fx_in * U["inner_x"] + fz_in * U["inner_z"]
                + fx_out * U["outer_x"] + fz_out * U["outer_z"] + fy * U["both_y"])

    # --- bearing-load envelope: load shared equally by both seats ---
    unit_uniform = {
        "x": 0.5 * (U["inner_x"] + U["outer_x"]),
        "y": U["both_y"],
        "z": 0.5 * (U["inner_z"] + U["outer_z"]),
    }
    f_unit = {
        "x": 0.5 * (load_cases["inner_x"] + load_cases["outer_x"]),
        "y": load_cases["both_y"],
        "z": 0.5 * (load_cases["inner_z"] + load_cases["outer_z"]),
    }
    axes = ["x", "y", "z"]
    C = np.array([[f_unit[a] @ unit_uniform[b] for b in axes] for a in axes])
    C = 0.5 * (C + C.T)
    R = np.diag([prm.radial_load, prm.axial_load, prm.radial_load])
    lam, vec = np.linalg.eigh(R @ C @ R)
    worst_load = R @ vec[:, -1]
    u_env = sum(worst_load[i] * unit_uniform[a] for i, a in enumerate(axes))

    def mean_seat(u):
        return 0.5 * (seat_disp(u, inner) + seat_disp(u, outer))

    seat_um = {
        "vertical 80 N": float(np.linalg.norm(mean_seat(prm.radial_load * unit_uniform["z"])) * 1e3),
        "horizontal 80 N": float(np.linalg.norm(mean_seat(prm.radial_load * unit_uniform["x"])) * 1e3),
        "axial 30 N": float(np.linalg.norm(mean_seat(prm.axial_load * unit_uniform["y"])) * 1e3),
        "worst-case envelope": float(np.linalg.norm(mean_seat(u_env)) * 1e3),
    }
    tilt = []
    for a_deg in range(-90, 91, 15):
        a = math.radians(a_deg)
        u = prm.radial_load * (math.sin(a) * unit_uniform["x"] - math.cos(a) * unit_uniform["z"])
        tilt.append((a_deg, float(np.linalg.norm(mean_seat(u)) * 1e3)))

    # --- overhung belt pull, every direction in the tilt plane ---
    y_in = 0.5 * prm.seat_depth
    y_out = prm.depth_y - 0.5 * prm.seat_depth
    y_pulley = prm.depth_y + PULLEY_OFFSET_MM
    span = y_out - y_in
    r_out = BELT_FORCE_N * (y_pulley - y_in) / span
    r_in = BELT_FORCE_N - r_out                       # negative: reacts the other way
    belt = []
    worst_belt = None
    for b_deg in range(0, 360, 15):
        b = math.radians(b_deg)
        cx, cz = math.sin(b), -math.cos(b)
        u = combine(r_in * cx, r_in * cz, r_out * cx, r_out * cz)
        d_in = seat_disp(u, inner)
        d_out = seat_disp(u, outer)
        rel = d_out - d_in
        misalign_rad = math.hypot(rel[0], rel[2]) / span
        entry = {
            "belt_dir_deg": b_deg,
            "outer_seat_um": float(np.linalg.norm(d_out) * 1e3),
            "inner_seat_um": float(np.linalg.norm(d_in) * 1e3),
            "misalignment_arcmin": float(math.degrees(misalign_rad) * 60.0),
        }
        belt.append(entry)
        if worst_belt is None or entry["misalignment_arcmin"] > worst_belt[0]["misalignment_arcmin"]:
            worst_belt = (entry, u)

    # --- stresses ---
    Dm = elasticity_matrix(prm.E0, prm.nu)
    B = strain_displacement((0.5, 0.5, 0.5), prm.h)
    solid = x_bin > 0.5

    def von_mises(u):
        ue = u[mesh.edof[solid]]
        s = (Dm @ (B @ ue.T)).T
        sx, sy, sz, txy, tyz, txz = s.T
        return np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)
                       + 3 * (txy ** 2 + tyz ** 2 + txz ** 2))

    vm_env = von_mises(u_env)
    vm_belt = von_mises(worst_belt[1])
    log(f"  done ({time.time() - t0:.0f}s)")
    return {
        "h_mm": prm.h,
        "dofs": int(mesh.n_dofs),
        "solid_elements": int(solid.sum()),
        "volume_cm3": float(solid.sum() * prm.h ** 3 / 1000),
        "envelope_eig_Nmm": [float(v) for v in lam],
        "worst_case_compliance_Nmm": float(lam[-1]),
        "worst_direction": [float(v) for v in vec[:, -1]],
        "seat_deflection_um": seat_um,
        "tilt_sweep": tilt,
        "belt": {
            "force_N": BELT_FORCE_N,
            "pulley_offset_mm": PULLEY_OFFSET_MM,
            "outer_reaction_N": r_out,
            "inner_reaction_N": r_in,
            "sweep": belt,
            "worst": worst_belt[0],
        },
        "vm_envelope_max_MPa": float(vm_env.max()),
        "vm_envelope_p99_MPa": float(np.percentile(vm_env, 99)),
        "vm_belt_max_MPa": float(vm_belt.max()),
        "vm_belt_p99_MPa": float(np.percentile(vm_belt, 99)),
        "seconds": time.time() - t0,
    }
