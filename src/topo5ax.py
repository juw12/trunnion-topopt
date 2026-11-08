"""
Tilt-robust, self-supporting 3D topology optimisation of the A-axis trunnion
upright for an AC-trunnion 5-axis FDM printer.

Units: N, mm, MPa.  Build direction = +Z (part printed standing on its base).

Pipeline per iteration
    x (design) -> density filter -> Heaviside projection -> AM overhang filter
    -> passive overrides -> SIMP stiffness -> 3 unit-load solves (MG-PCG)
    -> worst-case compliance over the bearing-load envelope (eigen p-norm)
"""
from __future__ import annotations

import itertools
import math
import time
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy import ndimage


# ----------------------------------------------------------------------------
# Problem definition (all named, editable parameters)
# ----------------------------------------------------------------------------
@dataclass
class Params:
    # design domain
    length_x: float = 120.0      # along the frame rail
    depth_y: float = 40.0        # along the trunnion axle
    height_z: float = 150.0      # build direction
    h: float = 2.5               # element size

    # bearing interface (2x 6001-2RS: 12 x 28 x 8 mm, one seat per face)
    axle_x: float = 60.0
    axle_z: float = 120.0
    seat_radius: float = 14.0    # 28 mm bearing OD
    seat_depth: float = 8.0      # bearing width
    shoulder_radius: float = 8.0 # 16 mm through-hole between the seats
    boss_radius: float = 22.0    # kept-solid ring around each seat
    boss_depth: float = 10.0

    # base bolting to the frame extrusion (4x M5)
    pad_width: float = 20.0
    pad_height: float = 7.5
    bolt_hole_radius: float = 2.75
    bolt_access_radius: float = 5.5   # socket-head + hex key keep-out above pad
    bolt_xy: tuple = ((10.0, 10.0), (10.0, 30.0), (110.0, 10.0), (110.0, 30.0))

    # load envelope at the bearing (any direction, ellipsoidal bound)
    radial_load: float = 80.0    # N, any direction in the tilt (XZ) plane
    axial_load: float = 30.0     # N, along the axle

    # material (PETG-CF / PA-CF, effective FDM values, assumed)
    E0: float = 3000.0
    nu: float = 0.35
    Emin_ratio: float = 1e-5
    penal: float = 3.0

    # optimisation
    target_volume_cm3: float = 105.0
    rmin: float = 2.5            # filter radius in elements (~6 mm members)
    eta: float = 0.5
    beta_schedule: tuple = ((0, 1.0), (20, 2.0), (40, 4.0), (60, 8.0), (80, 16.0))
    max_iter: int = 120
    move: float = 0.1
    pnorm: float = 8.0           # eigen p-norm -> smooth worst case
    use_envelope: bool = True    # False -> naive single vertical load case
    use_am_filter: bool = True   # False -> no overhang control
    symmetric_x: bool = True

    @property
    def nx(self):
        return int(round(self.length_x / self.h))

    @property
    def ny(self):
        return int(round(self.depth_y / self.h))

    @property
    def nz(self):
        return int(round(self.height_z / self.h))


# ----------------------------------------------------------------------------
# Element stiffness (8-node hex, 2x2x2 Gauss), node order:
# (0,0,0) (1,0,0) (1,1,0) (0,1,0) (0,0,1) (1,0,1) (1,1,1) (0,1,1)
# ----------------------------------------------------------------------------
LOCAL_NODES = np.array(
    [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
     [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)


def elasticity_matrix(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] = lam + 2 * mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    return D


def strain_displacement(xi, h):
    """B matrix (6x24) at natural coords xi in [0,1]^3 for a cube of size h."""
    s = 2 * LOCAL_NODES - 1                     # +-1 signs
    q = 2 * np.asarray(xi) - 1                  # map to [-1,1]
    dN = np.zeros((8, 3))
    for a in range(8):
        dN[a, 0] = s[a, 0] * (1 + s[a, 1] * q[1]) * (1 + s[a, 2] * q[2]) / 8
        dN[a, 1] = s[a, 1] * (1 + s[a, 0] * q[0]) * (1 + s[a, 2] * q[2]) / 8
        dN[a, 2] = s[a, 2] * (1 + s[a, 0] * q[0]) * (1 + s[a, 1] * q[1]) / 8
    dN *= 2.0 / h                               # d/dq -> d/dx
    B = np.zeros((6, 24))
    for a in range(8):
        dx, dy, dz = dN[a]
        B[0, 3 * a] = dx
        B[1, 3 * a + 1] = dy
        B[2, 3 * a + 2] = dz
        B[3, 3 * a] = dy
        B[3, 3 * a + 1] = dx
        B[4, 3 * a + 1] = dz
        B[4, 3 * a + 2] = dy
        B[5, 3 * a] = dz
        B[5, 3 * a + 2] = dx
    return B


def element_stiffness(h, nu):
    D = elasticity_matrix(1.0, nu)
    g = 0.5 / math.sqrt(3)
    pts = [0.5 - g, 0.5 + g]
    KE = np.zeros((24, 24))
    detJ = (h ** 3)
    for xi in itertools.product(pts, pts, pts):
        B = strain_displacement(xi, h)
        KE += B.T @ D @ B * detJ / 8
    return KE


# ----------------------------------------------------------------------------
# Structured mesh
# ----------------------------------------------------------------------------
class Mesh:
    def __init__(self, nx, ny, nz, h):
        self.nx, self.ny, self.nz, self.h = nx, ny, nz, h
        self.nnx, self.nny, self.nnz = nx + 1, ny + 1, nz + 1
        self.n_nodes = self.nnx * self.nny * self.nnz
        self.n_dofs = 3 * self.n_nodes
        self.n_el = nx * ny * nz
        i, j, k = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
        i = i.ravel(order="F")
        j = j.ravel(order="F")
        k = k.ravel(order="F")
        corner_nodes = []
        for (a, b, c) in LOCAL_NODES.astype(int):
            corner_nodes.append(self.node_id(i + a, j + b, k + c))
        corner_nodes = np.stack(corner_nodes, axis=1)            # (ne, 8)
        self.edof = (3 * corner_nodes[:, :, None] + np.arange(3)).reshape(-1, 24)
        self.el_centers = np.stack([(i + 0.5) * h, (j + 0.5) * h, (k + 0.5) * h], axis=1)
        ni, nj, nk = np.meshgrid(np.arange(self.nnx), np.arange(self.nny), np.arange(self.nnz),
                                 indexing="ij")
        self.node_xyz = np.stack([ni.ravel(order="F") * h,
                                  nj.ravel(order="F") * h,
                                  nk.ravel(order="F") * h], axis=1)

    def node_id(self, i, j, k):
        return i + self.nnx * (j + self.nny * k)

    def to3d(self, v):
        return v.reshape((self.nx, self.ny, self.nz), order="F")

    @staticmethod
    def to1d(a):
        return a.ravel(order="F")


# ----------------------------------------------------------------------------
# Geometric multigrid preconditioned CG
# ----------------------------------------------------------------------------
def interp_1d(n_coarse_el):
    """Linear interpolation from n_coarse_el+1 coarse nodes to 2*n_coarse_el+1 fine nodes."""
    nf = 2 * n_coarse_el + 1
    nc = n_coarse_el + 1
    rows, cols, vals = [], [], []
    for f in range(nf):
        if f % 2 == 0:
            rows.append(f)
            cols.append(f // 2)
            vals.append(1.0)
        else:
            rows += [f, f]
            cols += [f // 2, f // 2 + 1]
            vals += [0.5, 0.5]
    return sp.csr_matrix((vals, (rows, cols)), shape=(nf, nc))


class MultigridSolver:
    def __init__(self, mesh: Mesh, free_dofs, levels=3, smooth_steps=3, omega=0.6):
        self.mesh = mesh
        self.free = free_dofs
        self.smooth_steps = smooth_steps
        self.omega = omega
        self.prolongations = []
        nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
        for level in range(levels - 1):
            if nx % 2 or ny % 2 or nz % 2:
                break
            Px = interp_1d(nx // 2)
            Py = interp_1d(ny // 2)
            Pz = interp_1d(nz // 2)
            P_scalar = sp.kron(Pz, sp.kron(Py, Px, format="csr"), format="csr")
            P = sp.kron(P_scalar, sp.identity(3, format="csr"), format="csr")
            self.prolongations.append(P)
            nx, ny, nz = nx // 2, ny // 2, nz // 2
        # restrict first prolongation to free fine DOFs
        if self.prolongations:
            self.prolongations[0] = self.prolongations[0][free_dofs, :].tocsr()
        self.levels = len(self.prolongations) + 1

    def setup(self, A):
        self.A = [A]
        for P in self.prolongations:
            Ac = (P.T @ self.A[-1] @ P).tocsr()
            d = Ac.diagonal()
            zero = np.abs(d) < 1e-12 * max(np.abs(d).max(), 1.0)
            if zero.any():
                Ac = Ac + sp.diags(zero.astype(float) * max(np.abs(d).max(), 1.0))
            self.A.append(Ac)
        self.Dinv = []
        self.cheb = []
        rng = np.random.default_rng(0)
        for Al in self.A[:-1]:
            dinv = 1.0 / Al.diagonal()
            v = rng.standard_normal(Al.shape[0])
            lam = 1.0
            for _ in range(15):
                v = dinv * (Al @ v)
                lam = np.linalg.norm(v)
                v /= lam
            upper = 1.1 * lam
            lower = upper / 30.0
            self.Dinv.append(dinv)
            self.cheb.append((lower, upper))
        self.coarse_lu = spla.splu(self.A[-1].tocsc())

    def _smooth(self, level, x, b):
        """Chebyshev polynomial smoother on D^-1 A over [lower, upper]."""
        A = self.A[level]
        dinv = self.Dinv[level]
        lower, upper = self.cheb[level]
        theta = 0.5 * (upper + lower)
        delta = 0.5 * (upper - lower)
        sigma = theta / delta
        rho = 1.0 / sigma
        r = b - A @ x
        d = dinv * r / theta
        for _ in range(self.smooth_steps):
            x = x + d
            r = r - A @ d
            rho_new = 1.0 / (2.0 * sigma - rho)
            d = rho_new * rho * d + (2.0 * rho_new / delta) * (dinv * r)
            rho = rho_new
        return x

    def vcycle(self, level, b):
        if level == self.levels - 1:
            return self.coarse_lu.solve(b)
        A = self.A[level]
        x = self._smooth(level, np.zeros_like(b), b)
        r = b - A @ x
        P = self.prolongations[level]
        x = x + P @ self.vcycle(level + 1, P.T @ r)
        x = self._smooth(level, x, b)
        return x

    def solve(self, b, x0=None, tol=1e-6, maxit=300):
        A = self.A[0]
        x = np.zeros_like(b) if x0 is None else x0.copy()
        r = b - A @ x
        bnorm = np.linalg.norm(b)
        if bnorm == 0:
            return x, 0
        z = self.vcycle(0, r)
        p = z.copy()
        rz = r @ z
        for it in range(1, maxit + 1):
            Ap = A @ p
            alpha = rz / (p @ Ap)
            x += alpha * p
            r -= alpha * Ap
            if np.linalg.norm(r) / bnorm < tol:
                return x, it
            z = self.vcycle(0, r)
            rz_new = r @ z
            p = z + (rz_new / rz) * p
            rz = rz_new
        return x, maxit


# ----------------------------------------------------------------------------
# Additive-manufacturing overhang filter (Langelaar 2017, 3D, 5-point support)
# ----------------------------------------------------------------------------
class OverhangFilter:
    def __init__(self, passive_solid3d, P=40.0, eps=1e-4, xi0=0.5):
        self.P = P
        self.eps = eps
        n_support = 5
        self.Q = P + math.log(n_support) / math.log(xi0)
        self.solid = passive_solid3d

    @staticmethod
    def _neighbours(layer):
        """Stack the 5 supporting values for every element of the next layer."""
        pad = np.pad(layer, 1)
        return np.stack([
            pad[1:-1, 1:-1],
            pad[:-2, 1:-1], pad[2:, 1:-1],
            pad[1:-1, :-2], pad[1:-1, 2:],
        ])

    @staticmethod
    def _neighbours_adjoint(values):
        """Transpose of _neighbours: scatter 5 stacked arrays back onto the layer."""
        nx, ny = values.shape[1:]
        pad = np.zeros((nx + 2, ny + 2))
        pad[1:-1, 1:-1] += values[0]
        pad[:-2, 1:-1] += values[1]
        pad[2:, 1:-1] += values[2]
        pad[1:-1, :-2] += values[3]
        pad[1:-1, 2:] += values[4]
        return pad[1:-1, 1:-1]

    def forward(self, x):
        P, Q, eps = self.P, self.Q, self.eps
        nz = x.shape[2]
        xp = np.empty_like(x)
        self.cache = []
        xp[:, :, 0] = np.where(self.solid[:, :, 0], 1.0, x[:, :, 0])
        for k in range(1, nz):
            nb = self._neighbours(xp[:, :, k - 1])
            s = np.maximum((nb ** P).sum(axis=0), 1e-250)
            Xi = s ** (1.0 / Q)
            diff = x[:, :, k] - Xi
            root = np.sqrt(diff ** 2 + eps)
            val = 0.5 * (x[:, :, k] + Xi - root + math.sqrt(eps))
            d_dx = 0.5 * (1 - diff / root)
            d_dXi = 0.5 * (1 + diff / root)
            dXi_dnb = (P / Q) * nb ** (P - 1) * s ** (1.0 / Q - 1.0)
            solid = self.solid[:, :, k]
            val = np.where(solid, 1.0, val)
            d_dx = np.where(solid, 0.0, d_dx)
            d_dXi = np.where(solid, 0.0, d_dXi)
            xp[:, :, k] = val
            self.cache.append((d_dx, d_dXi, dXi_dnb))
        self.base_free = ~self.solid[:, :, 0]
        return xp

    def backward(self, grad_xp):
        nz = grad_xp.shape[2]
        grad_x = np.zeros_like(grad_xp)
        lam = np.zeros(grad_xp.shape[:2])
        for k in range(nz - 1, 0, -1):
            d_dx, d_dXi, dXi_dnb = self.cache[k - 1]
            total = grad_xp[:, :, k] + lam
            grad_x[:, :, k] = total * d_dx
            mu = total * d_dXi
            lam = self._neighbours_adjoint(mu[None] * dXi_dnb)
        grad_x[:, :, 0] = (grad_xp[:, :, 0] + lam) * self.base_free
        return grad_x


# ----------------------------------------------------------------------------
# Problem setup: passive regions, supports, unit loads
# ----------------------------------------------------------------------------
def build_regions(prm: Params, mesh: Mesh):
    c = mesh.el_centers
    x, y, z = c[:, 0], c[:, 1], c[:, 2]
    r_axle = np.hypot(x - prm.axle_x, z - prm.axle_z)
    in_face_zone = (y <= prm.boss_depth) | (y >= prm.depth_y - prm.boss_depth)
    in_seat_zone = (y <= prm.seat_depth) | (y >= prm.depth_y - prm.seat_depth)

    boss = (r_axle <= prm.boss_radius) & in_face_zone
    bore = ((r_axle < prm.seat_radius) & in_seat_zone) | (r_axle < prm.shoulder_radius)
    pads = ((x <= prm.pad_width) | (x >= prm.length_x - prm.pad_width)) & (z <= prm.pad_height)

    bolt_hole = np.zeros_like(x, dtype=bool)
    bolt_access = np.zeros_like(x, dtype=bool)
    for bx, by in prm.bolt_xy:
        rb = np.hypot(x - bx, y - by)
        bolt_hole |= rb < prm.bolt_hole_radius
        bolt_access |= (rb < prm.bolt_access_radius) & (z > prm.pad_height)

    solid = (boss | pads) & ~bore & ~bolt_hole
    void = bore | bolt_hole | bolt_access
    return solid, void


def build_supports_and_loads(prm: Params, mesh: Mesh):
    nodes = mesh.node_xyz
    x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    on_base = np.isclose(z, 0.0)
    in_pad = (x <= prm.pad_width + 1e-9) | (x >= prm.length_x - prm.pad_width - 1e-9)
    fixed_nodes = np.where(on_base & in_pad)[0]
    fixed = np.concatenate([3 * fixed_nodes + d for d in range(3)])
    free = np.setdiff1d(np.arange(mesh.n_dofs), fixed)

    # bearing seat surface nodes: radius within half an element of the seat
    r_axle = np.hypot(x - prm.axle_x, z - prm.axle_z)
    in_seat_zone = (y <= prm.seat_depth + 1e-9) | (y >= prm.depth_y - prm.seat_depth - 1e-9)
    seat = np.where((np.abs(r_axle - prm.seat_radius) <= 0.75 * prm.h) & in_seat_zone)[0]

    unit_loads = np.zeros((mesh.n_dofs, 3))
    for d in range(3):
        unit_loads[3 * seat + d, d] = 1.0 / len(seat)
    envelope = np.diag([prm.radial_load, prm.axial_load, prm.radial_load])
    return fixed, free, seat, unit_loads, envelope


# ----------------------------------------------------------------------------
# Optimiser
# ----------------------------------------------------------------------------
@dataclass
class History:
    objective: list = field(default_factory=list)
    volume: list = field(default_factory=list)
    change: list = field(default_factory=list)
    beta: list = field(default_factory=list)
    seconds: list = field(default_factory=list)
    pcg_its: list = field(default_factory=list)


class TrunnionTopOpt:
    def __init__(self, prm: Params, log=print):
        self.prm = prm
        self.log = log
        self.mesh = Mesh(prm.nx, prm.ny, prm.nz, prm.h)
        m = self.mesh
        self.KE = element_stiffness(prm.h, prm.nu)
        self.solid, self.void = build_regions(prm, m)
        self.passive = self.solid | self.void
        self.active = ~self.passive
        (self.fixed, self.free, self.seat_nodes,
         self.unit_loads, self.envelope) = build_supports_and_loads(prm, m)
        self._build_assembly_map()
        self.mg = MultigridSolver(m, self.free)
        self.el_volume = prm.h ** 3 / 1000.0               # cm^3
        self.target_volume = prm.target_volume_cm3
        # density filter (cone kernel)
        r = prm.rmin
        R = int(math.ceil(r)) - 1 if r == int(r) else int(math.floor(r))
        g = np.arange(-R, R + 1)
        gx, gy, gz = np.meshgrid(g, g, g, indexing="ij")
        kern = np.maximum(0.0, r - np.sqrt(gx ** 2 + gy ** 2 + gz ** 2))
        self.kernel = kern
        self.kernel_norm = ndimage.convolve(np.ones((m.nx, m.ny, m.nz)), kern, mode="constant")
        self.am = OverhangFilter(m.to3d(self.solid)) if prm.use_am_filter else None
        self.u_prev = None
        self.history = History()

    # --- assembly with precomputed CSR mapping on free DOFs -----------------
    def _build_assembly_map(self):
        m = self.mesh
        dof_to_free = -np.ones(m.n_dofs, dtype=np.int64)
        dof_to_free[self.free] = np.arange(len(self.free))
        ii = np.repeat(m.edof, 24, axis=1)                 # (ne, 576)
        jj = np.tile(m.edof, (1, 24))
        fi = dof_to_free[ii]
        fj = dof_to_free[jj]
        keep = (fi >= 0) & (fj >= 0)
        self.keep_mask = keep
        rows = fi[keep]
        cols = fj[keep]
        n = len(self.free)
        key = rows * n + cols
        uniq, inverse = np.unique(key, return_inverse=True)
        self.csr_map = inverse
        self.csr_rows = uniq // n
        self.csr_cols = uniq % n
        self.n_free = n
        template = sp.csr_matrix((np.ones(len(uniq)), (self.csr_rows, self.csr_cols)), shape=(n, n))
        template.sort_indices()
        # map unique (row, col) order to CSR data order
        order = np.lexsort((self.csr_cols, self.csr_rows))
        self.unique_to_csr = np.empty_like(order)
        self.unique_to_csr[order] = np.arange(len(order))
        self.K_template = template
        self.ke_flat = self.KE.ravel()

    def assemble(self, E):
        vals = (E[:, None] * self.ke_flat[None, :])[self.keep_mask]
        data_unique = np.bincount(self.csr_map, weights=vals, minlength=len(self.csr_rows))
        K = self.K_template.copy()
        # unique keys are sorted row-major == canonical CSR data order
        K.data = data_unique
        return K

    # --- density chain ------------------------------------------------------
    def filt(self, x3):
        return ndimage.convolve(x3, self.kernel, mode="constant") / self.kernel_norm

    def filt_adj(self, g3):
        return ndimage.convolve(g3 / self.kernel_norm, self.kernel, mode="constant")

    def project(self, xt, beta):
        eta = self.prm.eta
        den = math.tanh(beta * eta) + math.tanh(beta * (1 - eta))
        xb = (math.tanh(beta * eta) + np.tanh(beta * (xt - eta))) / den
        dxb = beta * (1 - np.tanh(beta * (xt - eta)) ** 2) / den
        return xb, dxb

    def physical(self, x, beta, keep_cache=False):
        m = self.mesh
        x3 = m.to3d(x)
        xt = self.filt(x3)
        xb, dxb = self.project(xt, beta)
        xb = np.where(m.to3d(self.solid), 1.0, xb)
        xb = np.where(m.to3d(self.void), 0.0, xb)
        if self.am is not None:
            xp = self.am.forward(xb)
        else:
            xp = xb
        xp = np.where(m.to3d(self.solid), 1.0, xp)
        xp = np.where(m.to3d(self.void), 0.0, xp)
        if keep_cache:
            self._dxb = dxb
        return m.to1d(xp)

    def chain_backward(self, grad_phys):
        """d(.)/d x_phys -> d(.)/d x_design (needs physical(..., keep_cache=True) first)."""
        m = self.mesh
        g = m.to3d(grad_phys).copy()
        g[m.to3d(self.passive)] = 0.0
        if self.am is not None:
            g = self.am.backward(g)
        g[m.to3d(self.passive)] = 0.0
        g = g * self._dxb
        g = self.filt_adj(g)
        g = m.to1d(g)
        g[self.passive] = 0.0
        if self.prm.symmetric_x:
            g3 = m.to3d(g)
            g = m.to1d(0.5 * (g3 + g3[::-1, :, :]))
        return g

    # --- FE analysis --------------------------------------------------------
    def analyse(self, x_phys, tol=1e-6):
        prm = self.prm
        E = prm.E0 * (prm.Emin_ratio + (1 - prm.Emin_ratio) * x_phys ** prm.penal)
        K = self.assemble(E)
        self.mg.setup(K)
        U = np.zeros((self.mesh.n_dofs, 3))
        its = []
        for d in range(3):
            b = self.unit_loads[self.free, d]
            x0 = None if self.u_prev is None else self.u_prev[self.free, d]
            u, it = self.mg.solve(b, x0=x0, tol=tol)
            U[self.free, d] = u
            its.append(it)
        self.u_prev = U
        return U, its

    def envelope_objective(self, x_phys, U):
        """Worst-case compliance over the load envelope (eigen p-norm) + sensitivities."""
        prm = self.prm
        m = self.mesh
        C = self.unit_loads.T @ U                            # 3x3 unit-load compliance
        C = 0.5 * (C + C.T)
        R = self.envelope
        if prm.use_envelope:
            M = R @ C @ R
            lam, vec = np.linalg.eigh(M)
            lam = np.maximum(lam, 1e-30)
            p = prm.pnorm
            phi = (lam ** p).sum() ** (1 / p)
            weights = (lam / phi) ** (p - 1)                # dphi/dlam_i
        else:
            # naive: single vertical bearing load (gravity case) only
            lam = np.array([C[2, 2] * R[2, 2] ** 2])
            vec = np.array([[0.0], [0.0], [1.0]])
            phi = lam[0]
            weights = np.array([1.0])
        dE = prm.E0 * (1 - prm.Emin_ratio) * prm.penal * x_phys ** (prm.penal - 1)
        grad = np.zeros(m.n_el)
        for w, v in zip(weights, vec.T):
            u_v = U @ (R @ v)                                 # displacement for the eigen-load
            ue = u_v[m.edof]
            energy = np.einsum("ij,jk,ik->i", ue, self.KE, ue)
            grad -= w * dE * energy
        return phi, grad, lam, C

    # --- optimisation loop --------------------------------------------------
    def beta_at(self, it):
        beta = 1.0
        for start, b in self.prm.beta_schedule:
            if it >= start:
                beta = b
        return beta

    def volume(self, x_phys):
        return x_phys.sum() * self.el_volume

    def initial_design(self):
        prm = self.prm
        x = np.zeros(self.mesh.n_el)
        solid_vol = self.solid.sum() * self.el_volume
        n_active = self.active.sum()
        frac = (self.target_volume - solid_vol) / (n_active * self.el_volume)
        x[self.active] = frac
        x[self.solid] = 1.0
        return x

    def optimise(self, stop_flag=None, x0=None):
        prm = self.prm
        x = self.initial_design() if x0 is None else x0.copy()
        phi0 = None
        for it in range(prm.max_iter):
            t0 = time.time()
            beta = self.beta_at(it)
            x_phys = self.physical(x, beta, keep_cache=True)
            U, its = self.analyse(x_phys)
            phi, dphi_phys, lam, C = self.envelope_objective(x_phys, U)
            if phi0 is None:
                phi0 = phi
            dphi = self.chain_backward(dphi_phys) / phi0
            dvol = self.chain_backward(np.ones_like(x_phys))

            # optimality criteria with exact (nonlinear) volume via bisection
            act = self.active
            l1, l2 = 1e-12, 1e12
            move = prm.move if beta < 8 else 0.5 * prm.move
            ratio = np.zeros_like(x)
            ratio[act] = np.sqrt(np.maximum(-dphi[act], 1e-30) / np.maximum(dvol[act], 1e-12))
            while (l2 - l1) / (l1 + l2) > 1e-4:
                lmid = math.sqrt(l1 * l2)
                xnew = x.copy()
                xnew[act] = np.clip(x[act] * ratio[act] / math.sqrt(lmid),
                                    np.maximum(0.0, x[act] - move),
                                    np.minimum(1.0, x[act] + move))
                if prm.symmetric_x:
                    x3 = self.mesh.to3d(xnew)
                    xnew = self.mesh.to1d(0.5 * (x3 + x3[::-1, :, :]))
                vol = self.volume(self.physical(xnew, beta))
                if vol > self.target_volume:
                    l1 = lmid
                else:
                    l2 = lmid
            change = np.abs(xnew - x).max()
            x = xnew
            h = self.history
            h.objective.append(phi)
            h.volume.append(self.volume(x_phys))
            h.change.append(change)
            h.beta.append(beta)
            h.seconds.append(time.time() - t0)
            h.pcg_its.append(its)
            self.log(f"it {it:3d}  beta {beta:4.0f}  phi {phi:10.4f}  vol {h.volume[-1]:6.1f} cm3"
                     f"  change {change:5.3f}  pcg {its}  {h.seconds[-1]:5.1f}s")
            self.x = x
            self.x_phys = x_phys
            if stop_flag is not None and stop_flag():
                break
        self.x_phys = self.physical(x, self.beta_at(prm.max_iter - 1))
        return x
