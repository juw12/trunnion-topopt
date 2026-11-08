"""
v2: tilt-robust + belt-robust, self-supporting trunnion upright.

Changes from v1
* Bearing loads act only on seat-surface nodes that touch non-void elements
  (the real contact between bearing outer ring and printed seat).
* Each seat is loaded separately, so the overhung belt pull on the axle pulley
  (bearings react unequally, axle tilts) becomes a design load case.
* Objective = p-norm over the eigenvalues of two worst-case matrices:
  bearing envelope (80 N any tilt-plane direction + 30 N axial, shared by the
  seats) and belt envelope (100 N in any tilt-plane direction on a pulley
  15 mm outboard), each normalised by its starting value.
"""
from __future__ import annotations

import numpy as np

from topo5ax import TrunnionTopOpt

BELT_FORCE_N = 100.0
PULLEY_OFFSET_MM = 15.0


class TrunnionTopOptV2(TrunnionTopOpt):
    def __init__(self, prm, log=print):
        super().__init__(prm, log=log)
        m = self.mesh
        nxyz = m.node_xyz
        attached = np.zeros(m.n_nodes, dtype=bool)
        attached[m.edof[~self.void][:, ::3] // 3] = True
        r_axle = np.hypot(nxyz[:, 0] - prm.axle_x, nxyz[:, 2] - prm.axle_z)
        y = nxyz[:, 1]
        band = (r_axle >= prm.seat_radius - 0.75 * prm.h) & (r_axle <= prm.seat_radius + 0.75 * prm.h)
        inner = np.nonzero(attached & band & (y <= prm.seat_depth + 1e-9))[0]
        outer = np.nonzero(attached & band & (y >= prm.depth_y - prm.seat_depth - 1e-9))[0]
        self.inner_nodes, self.outer_nodes = inner, outer

        def unit(nodes, d):
            f = np.zeros(m.n_dofs)
            f[3 * nodes + d] = 1.0 / len(nodes)
            return f

        # columns: inner_x, inner_z, outer_x, outer_z, both_y
        self.unit_loads = np.stack([
            unit(inner, 0),
            unit(inner, 2),
            unit(outer, 0),
            unit(outer, 2),
            0.5 * (unit(inner, 1) + unit(outer, 1)),
        ], axis=1)
        # uniform (shared) bearing load in x, y, z
        self.A_env = np.array([
            [0.5, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.5, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 1.0, 0.0],
        ])
        y_in = 0.5 * prm.seat_depth
        y_out = prm.depth_y - 0.5 * prm.seat_depth
        y_pulley = prm.depth_y + PULLEY_OFFSET_MM
        r_out = BELT_FORCE_N * (y_pulley - y_in) / (y_out - y_in)
        r_in = BELT_FORCE_N - r_out
        self.belt_reactions = (r_in, r_out)
        # belt pull of 1 (scaled below) in x and z
        self.A_belt = np.array([
            [r_in, 0.0],
            [0.0, r_in],
            [r_out, 0.0],
            [0.0, r_out],
            [0.0, 0.0],
        ])
        self.scales = None
        log(f"v2 loads: inner contact nodes {len(inner)}, outer {len(outer)}, "
            f"belt reactions inner {r_in:.1f} N, outer {r_out:.1f} N")

    def analyse(self, x_phys, tol=1e-6):
        prm = self.prm
        E = prm.E0 * (prm.Emin_ratio + (1 - prm.Emin_ratio) * x_phys ** prm.penal)
        K = self.assemble(E)
        self.mg.setup(K)
        n_cases = self.unit_loads.shape[1]
        U = np.zeros((self.mesh.n_dofs, n_cases))
        its = []
        for d in range(n_cases):
            b = self.unit_loads[self.free, d]
            x0 = None
            if self.u_prev is not None and self.u_prev.shape[1] == n_cases:
                x0 = self.u_prev[self.free, d]
            u, it = self.mg.solve(b, x0=x0, tol=tol)
            U[self.free, d] = u
            its.append(it)
        self.u_prev = U
        return U, its

    def envelope_objective(self, x_phys, U):
        prm = self.prm
        m = self.mesh
        G = self.unit_loads.T @ U
        G = 0.5 * (G + G.T)
        R = np.diag([prm.radial_load, prm.axial_load, prm.radial_load])
        M_env = R @ self.A_env.T @ G @ self.A_env @ R
        M_belt = self.A_belt.T @ G @ self.A_belt
        lam_env, vec_env = np.linalg.eigh(M_env)
        lam_belt, vec_belt = np.linalg.eigh(M_belt)
        if self.scales is None:
            self.scales = (lam_env.max(), lam_belt.max())
        s_env, s_belt = self.scales
        terms = []
        for lam, vec, s, A, pre in (
            (lam_env, vec_env, s_env, self.A_env, R),
            (lam_belt, vec_belt, s_belt, self.A_belt, np.eye(2)),
        ):
            for k in range(len(lam)):
                coeff = A @ (pre @ vec[:, k])
                terms.append((max(lam[k], 1e-30) / s, s, coeff))
        p = prm.pnorm
        phi = sum(t ** p for t, _, _ in terms) ** (1.0 / p)
        dE = prm.E0 * (1 - prm.Emin_ratio) * prm.penal * x_phys ** (prm.penal - 1)
        grad = np.zeros(m.n_el)
        for t, s, coeff in terms:
            w = (t / phi) ** (p - 1) / s
            u_c = U @ coeff
            ue = u_c[m.edof]
            energy = np.einsum("ij,jk,ik->i", ue, self.KE, ue)
            grad -= w * dE * energy
        self.last_eigs = {"env": lam_env.tolist(), "belt": lam_belt.tolist()}
        return phi, grad, np.concatenate([lam_env, lam_belt]), G
