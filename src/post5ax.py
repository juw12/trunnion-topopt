"""
Post-processing for the trunnion upright: exact interface features, watertight
surface extraction (marching tetrahedra), STL export, and verification.
"""
from __future__ import annotations

import itertools
import struct

import numpy as np
from scipy import ndimage

from topo5ax import Params, TrunnionTopOpt, elasticity_matrix, strain_displacement


# ----------------------------------------------------------------------------
# Implicit field: optimised density + exact CAD interface features
# ----------------------------------------------------------------------------
def build_field(prm: Params, x_phys3d: np.ndarray, voxel: float = 1.0):
    """Return (field, origin, spacing). field > 0 inside material, ~signed distance in mm."""
    pad = 2
    nx = int(round(prm.length_x / voxel)) + 1
    ny = int(round(prm.depth_y / voxel)) + 1
    nz = int(round(prm.height_z / voxel)) + 1
    gx = np.linspace(0, prm.length_x, nx)
    gy = np.linspace(0, prm.depth_y, ny)
    gz = np.linspace(0, prm.height_z, nz)
    gx = np.concatenate([gx[0] - voxel * np.arange(pad, 0, -1), gx, gx[-1] + voxel * np.arange(1, pad + 1)])
    gy = np.concatenate([gy[0] - voxel * np.arange(pad, 0, -1), gy, gy[-1] + voxel * np.arange(1, pad + 1)])
    gz = np.concatenate([gz[0] - voxel * np.arange(pad, 0, -1), gz, gz[-1] + voxel * np.arange(1, pad + 1)])
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")

    # density sampled at grid points: element-centred data, zero outside the domain
    rho = np.pad(x_phys3d, 1)                        # zero ring of ghost elements
    h = prm.h
    coords = [(X / h) + 0.5, (Y / h) + 0.5, (Z / h) + 0.5]   # index space of padded array
    rho_s = ndimage.map_coordinates(rho, coords, order=1, mode="constant", cval=0.0)
    rho_s = ndimage.gaussian_filter(rho_s, sigma=0.6 / voxel * 1.0)
    topo = (rho_s - 0.5) * h * 1.6                   # density -> approx. distance (mm)

    # base-pad bottom must be flat on the bed: clip anything below z = 0
    inside_domain = np.minimum.reduce([
        X, prm.length_x - X, Y, prm.depth_y - Y, Z, prm.height_z - Z,
    ])

    # exact features (inside-positive signed distances)
    r_axle = np.hypot(X - prm.axle_x, Z - prm.axle_z)
    face_band = np.maximum(-Y + 0.0, Y - prm.boss_depth)                   # >0 outside band
    face_band2 = np.maximum((prm.depth_y - prm.boss_depth) - Y, Y - prm.depth_y)
    boss = np.maximum(-np.minimum(face_band, face_band2) * 1.0 - 0.0, -1e9)
    boss = np.minimum(prm.boss_radius - r_axle, -np.minimum(face_band, face_band2))

    pad_x = np.maximum(
        np.minimum(prm.pad_width - X, X - (prm.length_x - prm.pad_width)) * -1.0,
        -1e9,
    )
    in_pad_x = np.maximum(prm.pad_width - X, X - (prm.length_x - prm.pad_width))
    pads = np.minimum.reduce([in_pad_x, prm.pad_height - Z, Z, Y, prm.depth_y - Y])

    # union of the two bearing-seat bands (inside-positive): front face OR back face
    seat_band = np.maximum(prm.seat_depth - Y, Y - (prm.depth_y - prm.seat_depth))
    seat_void = np.minimum(prm.seat_radius - r_axle, seat_band)
    shoulder_void = prm.shoulder_radius - r_axle
    bore = np.maximum(seat_void, shoulder_void)

    bolt_hole = np.full_like(X, -1e9)
    bolt_access = np.full_like(X, -1e9)
    for bx, by in prm.bolt_xy:
        rb = np.hypot(X - bx, Y - by)
        bolt_hole = np.maximum(bolt_hole, prm.bolt_hole_radius - rb)
        bolt_access = np.maximum(bolt_access, np.minimum(prm.bolt_access_radius - rb, Z - prm.pad_height))

    solid = np.maximum(topo, np.maximum(boss, pads))
    solid = np.minimum(solid, inside_domain)
    solid = np.minimum(solid, -bore)
    solid = np.minimum(solid, -bolt_hole)
    solid = np.minimum(solid, -bolt_access)
    origin = np.array([gx[0], gy[0], gz[0]])
    return solid, origin, voxel


# ----------------------------------------------------------------------------
# Marching tetrahedra (Freudenthal split -> consistent, closed 2-manifold)
# ----------------------------------------------------------------------------
def marching_tetrahedra(field, origin, spacing, level=0.0):
    f = field - level
    f = np.where(f == 0, 1e-9, f)
    nx, ny, nz = f.shape
    strides = np.array([1, nx, nx * ny])
    flat = f.ravel(order="F")

    # cubes with a sign change only
    cmin = np.minimum.reduce([f[a:nx - 1 + a, b:ny - 1 + b, c:nz - 1 + c]
                              for a, b, c in itertools.product((0, 1), repeat=3)])
    cmax = np.maximum.reduce([f[a:nx - 1 + a, b:ny - 1 + b, c:nz - 1 + c]
                              for a, b, c in itertools.product((0, 1), repeat=3)])
    ci, cj, ck = np.nonzero((cmin < 0) & (cmax > 0))
    base = ci + nx * (cj + ny * ck)

    unit = np.eye(3, dtype=int)
    tris_keys = []
    for perm in itertools.permutations(range(3)):
        offs = [np.zeros(3, dtype=int)]
        for axis in perm:
            offs.append(offs[-1] + unit[axis])
        vid = np.stack([base + int(o @ strides) for o in offs], axis=1)   # (n,4)
        val = flat[vid]
        inside = val > 0
        code = (inside * np.array([1, 2, 4, 8])).sum(axis=1)
        grad_dir = np.zeros((len(base), 3))
        for step, axis in enumerate(perm):
            grad_dir[:, axis] = val[:, step + 1] - val[:, step]
        for c in range(1, 15):
            sel = np.nonzero(code == c)[0]
            if len(sel) == 0:
                continue
            ins = [v for v in range(4) if c >> v & 1]
            outs = [v for v in range(4) if not c >> v & 1]
            if len(ins) == 1 or len(ins) == 3:
                lone = ins[0] if len(ins) == 1 else outs[0]
                others = [v for v in range(4) if v != lone]
                tri_list = [[(lone, others[0]), (lone, others[1]), (lone, others[2])]]
            else:
                a, b = ins
                c1, d1 = outs
                tri_list = [[(a, c1), (a, d1), (b, d1)], [(a, c1), (b, d1), (b, c1)]]
            for tri in tri_list:
                edges = np.stack([np.sort(np.stack([vid[sel, e0], vid[sel, e1]], axis=1), axis=1)
                                  for e0, e1 in tri], axis=1)   # (m,3,2)
                tris_keys.append((edges, grad_dir[sel]))

    edges = np.concatenate([e for e, _ in tris_keys])          # (T,3,2)
    grads = np.concatenate([g for _, g in tris_keys])          # (T,3) inward gradient
    keys = edges.reshape(-1, 2)
    uniq, inverse = np.unique(keys, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    # vertex positions on edges
    p0 = np.stack(np.unravel_index(uniq[:, 0], f.shape, order="F"), axis=1).astype(float)
    p1 = np.stack(np.unravel_index(uniq[:, 1], f.shape, order="F"), axis=1).astype(float)
    f0 = flat[uniq[:, 0]]
    f1 = flat[uniq[:, 1]]
    t = (f0 / (f0 - f1))[:, None]
    verts = origin + spacing * (p0 + t * (p1 - p0))
    faces = inverse.reshape(-1, 3)
    # orient outward: normal must oppose the (inward) field gradient
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    normals = np.cross(v1 - v0, v2 - v0)
    flip = (normals * grads).sum(axis=1) > 0
    faces[flip] = faces[flip][:, [0, 2, 1]]
    # drop degenerate triangles (collapsed vertices)
    good = (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 0] != faces[:, 2])
    return verts, faces[good]


def mesh_checks(verts, faces):
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    und = np.sort(e, axis=1)
    _, counts = np.unique(und, axis=0, return_counts=True)
    directed, dcounts = np.unique(e, axis=0, return_counts=True)
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    volume = np.einsum("ij,ij->i", v0, np.cross(v1, v2)).sum() / 6.0
    area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1).sum()
    return {
        "triangles": int(len(faces)),
        "vertices": int(len(verts)),
        "edges_not_shared_by_2": int((counts != 2).sum()),
        "directed_edge_duplicates": int((dcounts > 1).sum()),
        "watertight": bool((counts == 2).all()),
        "consistently_oriented": bool((dcounts == 1).all()),
        "volume_cm3": float(volume / 1000.0),
        "area_cm2": float(area / 100.0),
        "bbox_mm": [float(v) for v in (verts.max(0) - verts.min(0))],
    }


def write_stl(path, verts, faces, name="trunnion_upright"):
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    n = np.cross(v1 - v0, v2 - v0)
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    data = np.zeros(len(faces), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    data["n"] = n
    data["v"] = np.stack([v0, v1, v2], axis=1)
    with open(path, "wb") as fh:
        fh.write(name.encode()[:80].ljust(80, b" "))
        fh.write(struct.pack("<I", len(faces)))
        fh.write(data.tobytes())


# ----------------------------------------------------------------------------
# Verification on the thresholded (black/white) design
# ----------------------------------------------------------------------------
def verify_binary(opt: TrunnionTopOpt, x_bin: np.ndarray):
    prm = opt.prm
    m = opt.mesh
    opt.u_prev = None
    U, its = opt.analyse(x_bin.astype(float), tol=1e-9)
    C = opt.unit_loads.T @ U
    C = 0.5 * (C + C.T)
    R = opt.envelope
    M = R @ C @ R
    lam, vec = np.linalg.eigh(M)
    worst_dir = vec[:, -1]
    worst_load = R @ worst_dir
    u_worst = U @ worst_load

    seat = opt.seat_nodes
    seat_disp = {}
    for name, load in [
        ("vertical 80 N", np.array([0, 0, prm.radial_load])),
        ("horizontal 80 N", np.array([prm.radial_load, 0, 0])),
        ("axial 30 N", np.array([0, prm.axial_load, 0])),
        ("worst-case envelope", worst_load),
    ]:
        u = (U @ load).reshape(-1, 3)[seat]
        seat_disp[name] = float(np.linalg.norm(u.mean(axis=0)) * 1000.0)   # micron

    # tilt sweep: radial load direction in XZ plane for A = -90..+90 deg
    tilt = []
    for a_deg in range(-90, 91, 15):
        a = np.radians(a_deg)
        load = prm.radial_load * np.array([np.sin(a), 0.0, -np.cos(a)])
        u = (U @ load).reshape(-1, 3)[seat]
        tilt.append((a_deg, float(np.linalg.norm(u.mean(axis=0)) * 1000.0), float(load @ C @ load)))

    # von Mises at element centres for the worst-case load
    D = prm.E0 * 1.0
    Dm = elasticity_matrix(prm.E0, prm.nu)
    B = strain_displacement((0.5, 0.5, 0.5), prm.h)
    ue = u_worst[m.edof]
    sigma = (Dm @ (B @ ue.T)).T
    sx, sy, sz, txy, tyz, txz = sigma.T
    vm = np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3 * (txy ** 2 + tyz ** 2 + txz ** 2))
    vm_solid = vm[x_bin > 0.5]
    return {
        "eig_compliance_Nmm": [float(v) for v in lam],
        "worst_case_compliance_Nmm": float(lam[-1]),
        "worst_direction": [float(v) for v in worst_dir],
        "seat_deflection_um": seat_disp,
        "tilt_sweep": tilt,
        "vm_max_MPa": float(vm_solid.max()),
        "vm_p99_MPa": float(np.percentile(vm_solid, 99)),
        "pcg_its": its,
    }


def overhang_audit(opt: TrunnionTopOpt, x_bin: np.ndarray):
    """Solid, non-passive elements with none of the 5 supporting elements below (45 deg rule)."""
    m = opt.mesh
    s = m.to3d(x_bin > 0.5)
    passive = m.to3d(opt.solid)
    unsupported = np.zeros_like(s)
    below = np.pad(s, ((1, 1), (1, 1), (0, 0)))
    for k in range(1, s.shape[2]):
        b = below[:, :, k - 1]
        support = (b[1:-1, 1:-1] | b[:-2, 1:-1] | b[2:, 1:-1] | b[1:-1, :-2] | b[1:-1, 2:])
        unsupported[:, :, k] = s[:, :, k] & ~support & ~passive[:, :, k]
    n_down = 0
    for k in range(1, s.shape[2]):
        n_down += int((s[:, :, k] & ~s[:, :, k - 1]).sum())
    labels, n_comp = ndimage.label(s)
    sizes = np.bincount(labels.ravel())[1:]
    return {
        "unsupported_elements": int(unsupported.sum()),
        "downfacing_elements": int(n_down),
        "unsupported_fraction_of_downfacing": float(unsupported.sum() / max(n_down, 1)),
        "connected_components": int(n_comp),
        "largest_component_share": float(sizes.max() / sizes.sum()) if len(sizes) else 0.0,
    }


def build_cad_body_field(prm: Params, x_phys3d: np.ndarray, voxel: float = 1.5, overlap: float = 1.0):
    """Optimised body only: interface regions (bosses, pads) hollowed out except an
    `overlap` mm skin, so exact CAD features can be unioned with a clean volumetric
    overlap, and exact holes cut afterwards."""
    field, origin, spacing = build_field(prm, x_phys3d, voxel=voxel)
    pad = 2
    n = field.shape
    gx = origin[0] + spacing * np.arange(n[0])
    gy = origin[1] + spacing * np.arange(n[1])
    gz = origin[2] + spacing * np.arange(n[2])
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    r_axle = np.hypot(X - prm.axle_x, Z - prm.axle_z)
    in_face = np.maximum(prm.boss_depth - Y, Y - (prm.depth_y - prm.boss_depth))
    boss = np.minimum(prm.boss_radius - r_axle, np.minimum(in_face, np.minimum(Y, prm.depth_y - Y)))
    in_pad_x = np.maximum(prm.pad_width - X, X - (prm.length_x - prm.pad_width))
    pads = np.minimum.reduce([in_pad_x, prm.pad_height - Z, Z, Y, prm.depth_y - Y])
    body = np.minimum(field, -(boss - overlap))
    body = np.minimum(body, -(pads - overlap))
    return body, origin, spacing
