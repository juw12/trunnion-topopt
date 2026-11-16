"""Overhang (45 degree) rule, surface extraction and STL export on small known cases."""
import struct
import sys
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from post5ax import marching_tetrahedra, mesh_checks, overhang_audit, write_stl  # noqa: E402
from topo5ax import Mesh, OverhangFilter  # noqa: E402

GRID = (7, 5, 6)


def cube_mesh():
    """Unit cube, 12 triangles, outward normals."""
    vertices = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float)
    faces = np.array([
        [0, 2, 1], [0, 3, 2],
        [4, 5, 6], [4, 6, 7],
        [0, 1, 5], [0, 5, 4],
        [1, 2, 6], [1, 6, 5],
        [2, 3, 7], [2, 7, 6],
        [3, 0, 4], [3, 4, 7],
    ])
    return vertices, faces


def audit_of(solid_grid):
    """Run the repository's audit on a hand-made boolean grid with no passive region."""
    mesh = Mesh(*solid_grid.shape, 1.0)
    stand_in = types.SimpleNamespace(mesh=mesh, solid=np.zeros(mesh.n_el, dtype=bool))
    return overhang_audit(stand_in, Mesh.to1d(solid_grid.astype(float)))


class OverhangFilterTest(unittest.TestCase):
    def setUp(self):
        self.filter = OverhangFilter(np.zeros(GRID, dtype=bool))

    def test_vertical_column_on_the_bed_is_kept(self):
        design = np.zeros(GRID)
        design[3, 2, :] = 1.0
        printed = self.filter.forward(design)
        self.assertTrue(np.allclose(printed[3, 2, :], 1.0, atol=1e-3))

    def test_45_degree_staircase_is_kept(self):
        design = np.zeros(GRID)
        for level in range(GRID[2]):
            design[level, 2, level] = 1.0
        printed = self.filter.forward(design)
        for level in range(GRID[2]):
            self.assertGreater(printed[level, 2, level], 0.99, level)

    def test_floating_element_is_removed(self):
        design = np.zeros(GRID)
        design[3, 2, 3] = 1.0
        self.assertLess(self.filter.forward(design).max(), 0.01)

    def test_overhang_steeper_than_45_degrees_is_removed(self):
        design = np.zeros(GRID)
        design[0, 2, 0] = 1.0
        design[2, 2, 1] = 1.0          # two elements sideways for one element up
        printed = self.filter.forward(design)
        self.assertGreater(printed[0, 2, 0], 0.99)
        self.assertLess(printed[2, 2, 1], 0.01)

    def test_printed_density_never_exceeds_the_design_by_more_than_the_smoothing(self):
        rng = np.random.default_rng(5)
        design = rng.uniform(0.0, 1.0, size=GRID)
        printed = self.filter.forward(design)
        self.assertLessEqual((printed - design).max(), 0.5 * np.sqrt(self.filter.eps) + 1e-12)

    def test_passive_solid_region_stays_solid(self):
        passive = np.zeros(GRID, dtype=bool)
        passive[3, 2, 4] = True
        printed = OverhangFilter(passive).forward(np.zeros(GRID))
        self.assertEqual(printed[3, 2, 4], 1.0)

    def test_neighbour_scatter_is_the_transpose_of_the_gather(self):
        rng = np.random.default_rng(2)
        layer = rng.standard_normal(GRID[:2])
        stacked = rng.standard_normal((5,) + GRID[:2])
        gathered = (OverhangFilter._neighbours(layer) * stacked).sum()
        scattered = (layer * OverhangFilter._neighbours_adjoint(stacked)).sum()
        self.assertAlmostEqual(gathered, scattered, places=10)

    def test_backward_pass_matches_finite_differences(self):
        rng = np.random.default_rng(3)
        design = rng.uniform(0.2, 0.9, size=GRID)
        weights = rng.standard_normal(GRID)
        self.filter.forward(design)
        gradient = self.filter.backward(weights)
        step = 1e-6
        for _ in range(25):
            index = tuple(int(rng.integers(0, size)) for size in GRID)
            plus = design.copy()
            plus[index] += step
            minus = design.copy()
            minus[index] -= step
            numeric = (
                (self.filter.forward(plus) * weights).sum()
                - (self.filter.forward(minus) * weights).sum()
            ) / (2 * step)
            self.assertAlmostEqual(gradient[index], numeric, delta=1e-5 * max(1.0, abs(numeric)))


class OverhangAuditTest(unittest.TestCase):
    def test_staircase_passes(self):
        solid = np.zeros(GRID, dtype=bool)
        for level in range(GRID[2]):
            solid[level, 2, level] = True
        audit = audit_of(solid)
        self.assertEqual(audit["unsupported_elements"], 0)

    def test_floating_element_is_counted_and_splits_the_part(self):
        solid = np.zeros(GRID, dtype=bool)
        solid[:, :, 0] = True
        solid[3, 2, 3] = True
        audit = audit_of(solid)
        self.assertEqual(audit["unsupported_elements"], 1)
        self.assertEqual(audit["downfacing_elements"], 1)
        self.assertEqual(audit["unsupported_fraction_of_downfacing"], 1.0)
        self.assertEqual(audit["connected_components"], 2)

    def test_flat_bridge_counts_only_its_unsupported_middle(self):
        solid = np.zeros(GRID, dtype=bool)
        solid[0, 2, 0:3] = True        # left post
        solid[6, 2, 0:3] = True        # right post
        solid[:, 2, 3] = True          # 7-element bridge across the posts
        audit = audit_of(solid)
        # elements 0, 1 and 5, 6 are held by a post directly or diagonally below
        self.assertEqual(audit["unsupported_elements"], 3)
        self.assertEqual(audit["downfacing_elements"], 5)
        self.assertEqual(audit["connected_components"], 1)


class MeshChecksTest(unittest.TestCase):
    def test_unit_cube(self):
        checks = mesh_checks(*cube_mesh())
        self.assertTrue(checks["watertight"])
        self.assertTrue(checks["consistently_oriented"])
        self.assertEqual(checks["triangles"], 12)
        self.assertAlmostEqual(checks["volume_cm3"], 1.0 / 1000.0, places=12)
        self.assertAlmostEqual(checks["area_cm2"], 6.0 / 100.0, places=12)
        self.assertEqual(checks["bbox_mm"], [1.0, 1.0, 1.0])

    def test_missing_triangle_is_not_watertight(self):
        vertices, faces = cube_mesh()
        checks = mesh_checks(vertices, faces[:-1])
        self.assertFalse(checks["watertight"])
        self.assertEqual(checks["edges_not_shared_by_2"], 3)

    def test_flipped_triangle_is_not_consistently_oriented(self):
        vertices, faces = cube_mesh()
        faces = faces.copy()
        faces[0] = faces[0][[0, 2, 1]]
        checks = mesh_checks(vertices, faces)
        self.assertTrue(checks["watertight"])
        self.assertFalse(checks["consistently_oriented"])


class MarchingTetrahedraTest(unittest.TestCase):
    def test_sphere_is_closed_oriented_and_has_the_right_volume(self):
        spacing = 0.5
        radius = 8.0
        axis = np.arange(-10.0, 10.0 + spacing / 2, spacing)
        x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
        field = radius - np.sqrt(x ** 2 + y ** 2 + z ** 2)      # positive inside
        origin = np.array([axis[0], axis[0], axis[0]])
        vertices, faces = marching_tetrahedra(field, origin, spacing)
        checks = mesh_checks(vertices, faces)
        self.assertTrue(checks["watertight"])
        self.assertTrue(checks["consistently_oriented"])
        exact_volume = 4.0 / 3.0 * np.pi * radius ** 3 / 1000.0
        exact_area = 4.0 * np.pi * radius ** 2 / 100.0
        # positive signed volume also proves the normals point outward
        self.assertAlmostEqual(checks["volume_cm3"], exact_volume, delta=0.01 * exact_volume)
        self.assertAlmostEqual(checks["area_cm2"], exact_area, delta=0.02 * exact_area)
        distance = np.linalg.norm(vertices, axis=1)
        self.assertLess(np.abs(distance - radius).max(), 0.05)


class StlExportTest(unittest.TestCase):
    def test_binary_stl_layout_and_contents(self):
        vertices, faces = cube_mesh()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "cube.stl"
            write_stl(path, vertices, faces, name="cube")
            data = path.read_bytes()
        self.assertEqual(len(data), 84 + 50 * len(faces))
        self.assertTrue(data[:80].startswith(b"cube"))
        self.assertEqual(struct.unpack("<I", data[80:84])[0], len(faces))
        record = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        triangles = np.frombuffer(data[84:], dtype=record)
        self.assertTrue(np.allclose(triangles["v"], vertices[faces]))
        self.assertTrue(np.allclose(triangles["n"][0], [0.0, 0.0, -1.0]))
        self.assertTrue(np.allclose(np.linalg.norm(triangles["n"], axis=1), 1.0))


if __name__ == "__main__":
    unittest.main()
