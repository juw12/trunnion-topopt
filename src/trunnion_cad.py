"""
Editable CAD rebuild of the 5-axis trunnion upright (CadQuery / OpenCascade).

Exact, parametric B-rep for every functional interface; the topology-optimised
body is brought in as a mesh-derived solid and unioned with a volumetric overlap.
Edit the parameters below and re-run to regenerate the STEP.

Usage:
    python trunnion_cad.py <optimised_body.stl> <output.step>
"""
from __future__ import annotations

import sys
import time

import cadquery as cq
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeSolid, BRepBuilderAPI_Sewing
from OCP.ShapeFix import ShapeFix_Solid
from OCP.StlAPI import StlAPI_Reader
from OCP.TopAbs import TopAbs_SHELL
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS, TopoDS_Shape

# ---------------------------------------------------------------- parameters (mm)
LENGTH_X = 120.0          # along the frame rail
DEPTH_Y = 40.0            # along the trunnion axle
AXLE_X = 60.0
AXLE_Z = 120.0
SEAT_DIAMETER = 28.0      # 6001-2RS outer ring, nominal (press fit: ream to H7 after printing)
SEAT_DEPTH = 8.0          # bearing width
SHOULDER_DIAMETER = 16.0  # through-bore between the seats (12 mm axle + clearance)
BOSS_DIAMETER = 44.0
BOSS_DEPTH = 10.0
PAD_WIDTH = 20.0
PAD_HEIGHT = 7.5
BOLT_CLEARANCE_DIAMETER = 5.5   # M5
HEX_KEY_ACCESS_DIAMETER = 11.0  # socket head + key, kept clear above each pad
BOLT_XY = [(10.0, 10.0), (10.0, 30.0), (110.0, 10.0), (110.0, 30.0)]
# --------------------------------------------------------------------------------


def mesh_solid(stl_path: str) -> cq.Solid:
    reader = StlAPI_Reader()
    shape = TopoDS_Shape()
    if not reader.Read(shape, stl_path):
        raise RuntimeError(f"cannot read {stl_path}")
    sewing = BRepBuilderAPI_Sewing(1e-3)
    sewing.Add(shape)
    sewing.Perform()
    sewn = sewing.SewedShape()
    explorer = TopExp_Explorer(sewn, TopAbs_SHELL)
    maker = BRepBuilderAPI_MakeSolid()
    while explorer.More():
        maker.Add(TopoDS.Shell_s(explorer.Current()))
        explorer.Next()
    fixer = ShapeFix_Solid(maker.Solid())
    fixer.Perform()
    return cq.Solid(fixer.Solid())


def interface_solids() -> cq.Workplane:
    # bearing bosses: cylinders along the axle (Y), one per face
    boss_front = (cq.Workplane("XZ", origin=(0, 0, 0))
                  .center(AXLE_X, AXLE_Z).circle(BOSS_DIAMETER / 2).extrude(-BOSS_DEPTH))
    boss_back = (cq.Workplane("XZ", origin=(0, DEPTH_Y - BOSS_DEPTH, 0))
                 .center(AXLE_X, AXLE_Z).circle(BOSS_DIAMETER / 2).extrude(-BOSS_DEPTH))
    pads = (cq.Workplane("XY")
            .box(PAD_WIDTH, DEPTH_Y, PAD_HEIGHT, centered=False)
            .union(cq.Workplane("XY", origin=(LENGTH_X - PAD_WIDTH, 0, 0))
                   .box(PAD_WIDTH, DEPTH_Y, PAD_HEIGHT, centered=False)))
    return boss_front.union(boss_back).union(pads)


def interface_cuts() -> cq.Workplane:
    seat_front = (cq.Workplane("XZ", origin=(0, 0, 0))
                  .center(AXLE_X, AXLE_Z).circle(SEAT_DIAMETER / 2).extrude(-SEAT_DEPTH))
    seat_back = (cq.Workplane("XZ", origin=(0, DEPTH_Y - SEAT_DEPTH, 0))
                 .center(AXLE_X, AXLE_Z).circle(SEAT_DIAMETER / 2).extrude(-SEAT_DEPTH))
    bore = (cq.Workplane("XZ", origin=(0, 1.0, 0))
            .center(AXLE_X, AXLE_Z).circle(SHOULDER_DIAMETER / 2).extrude(-(DEPTH_Y + 2.0)))
    bore = bore.translate((0, -2.0, 0))
    cuts = seat_front.union(seat_back).union(bore)
    for bx, by in BOLT_XY:
        hole = (cq.Workplane("XY", origin=(0, 0, -1.0))
                .center(bx, by).circle(BOLT_CLEARANCE_DIAMETER / 2).extrude(PAD_HEIGHT + 2.0))
        access = (cq.Workplane("XY", origin=(0, 0, PAD_HEIGHT))
                  .center(bx, by).circle(HEX_KEY_ACCESS_DIAMETER / 2).extrude(200.0))
        cuts = cuts.union(hole).union(access)
    return cuts


def build(stl_path: str, step_path: str) -> dict:
    t0 = time.time()
    body = cq.Workplane("XY").add(mesh_solid(stl_path))
    t_mesh = time.time() - t0
    print(f"mesh solid ready {t_mesh:.0f}s", flush=True)
    part = body.union(interface_solids()).cut(interface_cuts())
    solid = part.val()
    cq.exporters.export(part, step_path)
    bb = solid.BoundingBox()
    return {
        "valid": bool(solid.isValid()),
        "volume_cm3": solid.Volume() / 1000.0,
        "bbox_mm": [bb.xlen, bb.ylen, bb.zlen],
        "solids": len(part.solids().vals()),
        "faces": len(solid.Faces()),
        "mesh_import_s": t_mesh,
        "total_s": time.time() - t0,
    }


if __name__ == "__main__":
    info = build(sys.argv[1], sys.argv[2])
    print(info)
