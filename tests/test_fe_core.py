"""Finite-element building blocks checked against closed-form answers.

These tests need only numpy and scipy and run in a few seconds. They do not
re-run the optimisation or the 1.16 M DOF verification solve.
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from topo5ax import (  # noqa: E402
    LOCAL_NODES,
    Mesh,
    MultigridSolver,
    elasticity_matrix,
    element_stiffness,
    interp_1d,
    strain_displacement,
)


def assemble_uniform(mesh, poisson):
    """Global stiffness of a uniform block with Young's modulus 1."""
    element = element_stiffness(mesh.h, poisson)
    rows = np.repeat(mesh.edof, 24, axis=1).ravel()
    cols = np.tile(mesh.edof, (1, 24)).ravel()
    values = np.tile(element.ravel(), mesh.n_el)
    return sp.csr_matrix((values, (rows, cols)), shape=(mesh.n_dofs, mesh.n_dofs))


def dofs_of(nodes):
    return np.concatenate([3 * nodes + direction for direction in range(3)])


class ElasticityMatrixTest(unittest.TestCase):
    def test_uniaxial_stress_gives_youngs_modulus_and_poisson_ratio(self):
        youngs, poisson = 3000.0, 0.35
        stress = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        strain = np.linalg.solve(elasticity_matrix(youngs, poisson), stress)
        self.assertAlmostEqual(strain[0], 1.0 / youngs, places=12)
        self.assertAlmostEqual(strain[1], -poisson / youngs, places=12)
        self.assertAlmostEqual(strain[2], -poisson / youngs, places=12)

    def test_shear_terms_equal_the_shear_modulus(self):
        youngs, poisson = 3000.0, 0.35
        matrix = elasticity_matrix(youngs, poisson)
        shear_modulus = youngs / (2 * (1 + poisson))
        for index in (3, 4, 5):
            self.assertAlmostEqual(matrix[index, index], shear_modulus, places=9)

    def test_matrix_is_symmetric_positive_definite(self):
        matrix = elasticity_matrix(1.0, 0.35)
        self.assertTrue(np.allclose(matrix, matrix.T))
        self.assertGreater(np.linalg.eigvalsh(matrix).min(), 0.0)


class StrainDisplacementTest(unittest.TestCase):
    def test_linear_displacement_field_gives_its_exact_strain(self):
        size = 2.5
        rng = np.random.default_rng(1)
        gradient = rng.standard_normal((3, 3))
        nodes = LOCAL_NODES * size
        displacement = (nodes @ gradient.T).ravel()
        expected = np.array([
            gradient[0, 0],
            gradient[1, 1],
            gradient[2, 2],
            gradient[0, 1] + gradient[1, 0],
            gradient[1, 2] + gradient[2, 1],
            gradient[0, 2] + gradient[2, 0],
        ])
        for point in [(0.5, 0.5, 0.5), (0.1, 0.9, 0.3), (0.0, 0.0, 1.0)]:
            strain = strain_displacement(point, size) @ displacement
            self.assertTrue(np.allclose(strain, expected, atol=1e-12), point)


class ElementStiffnessTest(unittest.TestCase):
    def setUp(self):
        self.size = 2.5
        self.element = element_stiffness(self.size, 0.35)

    def test_symmetric(self):
        self.assertTrue(np.allclose(self.element, self.element.T, atol=1e-12))

    def test_exactly_six_rigid_body_modes(self):
        eigenvalues = np.linalg.eigvalsh(self.element)
        scale = eigenvalues.max()
        self.assertEqual(int((np.abs(eigenvalues) < 1e-9 * scale).sum()), 6)
        self.assertGreater(eigenvalues[6], 1e-6 * scale)

    def test_translations_and_rotations_produce_no_force(self):
        nodes = LOCAL_NODES * self.size
        for direction in range(3):
            translation = np.zeros((8, 3))
            translation[:, direction] = 1.0
            self.assertTrue(np.allclose(self.element @ translation.ravel(), 0.0, atol=1e-10))
        for axis in np.eye(3):
            rotation = np.cross(axis, nodes)
            self.assertTrue(np.allclose(self.element @ rotation.ravel(), 0.0, atol=1e-10))

    def test_stiffness_scales_linearly_with_element_size(self):
        unit = element_stiffness(1.0, 0.35)
        self.assertTrue(np.allclose(self.element, self.size * unit, rtol=1e-12))


class MeshTest(unittest.TestCase):
    def test_counts_and_element_connectivity(self):
        mesh = Mesh(4, 3, 2, 2.5)
        self.assertEqual(mesh.n_el, 24)
        self.assertEqual(mesh.n_nodes, 5 * 4 * 3)
        self.assertEqual(mesh.edof.shape, (24, 24))
        # every element's 8 corner nodes sit at its centre +- half an element
        corner_nodes = mesh.edof[:, ::3] // 3
        corners = mesh.node_xyz[corner_nodes]
        self.assertTrue(np.allclose(corners.mean(axis=1), mesh.el_centers))
        self.assertTrue(np.allclose(np.ptp(corners, axis=1), mesh.h))

    def test_flat_and_grid_views_round_trip(self):
        mesh = Mesh(4, 3, 2, 1.0)
        values = np.arange(mesh.n_el, dtype=float)
        self.assertTrue(np.array_equal(Mesh.to1d(mesh.to3d(values)), values))
        grid = mesh.to3d(mesh.el_centers[:, 0])
        self.assertTrue(np.allclose(grid[:, 0, 0], [0.5, 1.5, 2.5, 3.5]))


class AxialBarTest(unittest.TestCase):
    def test_bar_extension_matches_force_length_over_area_modulus(self):
        # Poisson's ratio 0: a clamped bar under uniform end traction is in pure
        # uniaxial stress, so the exact answer is u = F L / (E A) with E = 1.
        mesh = Mesh(6, 2, 2, 2.0)
        stiffness = assemble_uniform(mesh, 0.0)
        x = mesh.node_xyz[:, 0]
        length = 6 * 2.0
        area = (2 * 2.0) * (2 * 2.0)
        clamped = np.where(np.isclose(x, 0.0))[0]
        free = np.setdiff1d(np.arange(mesh.n_dofs), dofs_of(clamped))

        force = 7.0
        end_face = np.isclose(x, length)
        load = np.zeros(mesh.n_dofs)
        # consistent nodal loads of a uniform traction: weight = share of face area
        for element_index in range(mesh.n_el):
            nodes = mesh.edof[element_index, ::3] // 3
            on_face = nodes[end_face[nodes]]
            if len(on_face) == 4:
                load[3 * on_face] += force * (mesh.h ** 2 / area) / 4
        self.assertAlmostEqual(load.sum(), force, places=12)

        displacement = np.zeros(mesh.n_dofs)
        displacement[free] = spla.spsolve(stiffness[free][:, free].tocsc(), load[free])
        end_extension = displacement[3 * np.where(end_face)[0]]
        self.assertTrue(np.allclose(end_extension, force * length / area, rtol=1e-9))
        # the bar stretches linearly along its length and does not move sideways
        self.assertTrue(np.allclose(displacement[0::3], force * x / area, atol=1e-9))
        self.assertTrue(np.allclose(displacement[1::3], 0.0, atol=1e-9))
        self.assertTrue(np.allclose(displacement[2::3], 0.0, atol=1e-9))


class MultigridTest(unittest.TestCase):
    def test_interpolation_reproduces_linear_functions(self):
        prolongation = interp_1d(4)
        self.assertEqual(prolongation.shape, (9, 5))
        self.assertTrue(np.allclose(prolongation.sum(axis=1), 1.0))
        coarse = np.linspace(0.0, 4.0, 5)
        self.assertTrue(np.allclose(prolongation @ coarse, np.linspace(0.0, 4.0, 9)))

    def test_multigrid_solution_matches_direct_solve_on_a_cantilever(self):
        mesh = Mesh(8, 4, 4, 1.0)
        stiffness = assemble_uniform(mesh, 0.3)
        x = mesh.node_xyz[:, 0]
        clamped = np.where(np.isclose(x, 0.0))[0]
        free = np.setdiff1d(np.arange(mesh.n_dofs), dofs_of(clamped))
        tip = np.where(np.isclose(x, 8.0))[0]
        load = np.zeros(mesh.n_dofs)
        load[3 * tip + 2] = -1.0 / len(tip)

        reduced = stiffness[free][:, free].tocsr()
        direct = spla.spsolve(reduced.tocsc(), load[free])

        solver = MultigridSolver(mesh, free)
        self.assertEqual(solver.levels, 3)
        solver.setup(reduced)
        solution, iterations = solver.solve(load[free], tol=1e-10)
        self.assertLess(iterations, 60)
        error = np.linalg.norm(solution - direct) / np.linalg.norm(direct)
        self.assertLess(error, 1e-8)

        # sanity of the physics: the tip deflection lies between the Euler-Bernoulli
        # beam (no shear) and the Timoshenko beam (shear factor 5/6) for E = 1
        tip_deflection = -direct[np.searchsorted(free, 3 * tip + 2)].mean()
        inertia = 4 * 4 ** 3 / 12
        euler = 8.0 ** 3 / (3 * inertia)
        timoshenko = euler + 8.0 / ((1 / (2 * 1.3)) * 16 * 5 / 6)
        self.assertGreater(tip_deflection, euler)
        self.assertLess(tip_deflection, timoshenko)

    def test_zero_load_returns_zero_without_iterating(self):
        mesh = Mesh(4, 2, 2, 1.0)
        stiffness = assemble_uniform(mesh, 0.3)
        clamped = np.where(np.isclose(mesh.node_xyz[:, 0], 0.0))[0]
        free = np.setdiff1d(np.arange(mesh.n_dofs), dofs_of(clamped))
        solver = MultigridSolver(mesh, free)
        solver.setup(stiffness[free][:, free].tocsr())
        solution, iterations = solver.solve(np.zeros(len(free)))
        self.assertEqual(iterations, 0)
        self.assertTrue(np.array_equal(solution, np.zeros(len(free))))


if __name__ == "__main__":
    unittest.main()
