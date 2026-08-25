"""
Seeded integration tests for the two existing noiseless VQE workflows.
This file tests the two noiseless VQE workflows. It contains helper functions and two actual pytest tests.
"""

import warnings

import numpy as np
import pytest
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import Statevector
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import SLSQP
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit_nature.second_q.mappers import JordanWignerMapper
from scipy.optimize import minimize
from scipy.sparse import SparseEfficiencyWarning

from quantum_hubbard import (
    N_SITES,
    double_occupancy_operator,
    expectation,
    hubbard_hamiltonian,
    particle_number_operator,
    solve_particle_sector,
    spin_correlation_operator,
)
from quantum_hubbard.operators import (
    pauli_decomposition,
    pauli_matrix,
    qiskit_hubbard_fermionic_op,
)

T = 1.0
U = 4.0
ATOL = 1e-10
VQE_ATOL = 1e-8
EXPECTED_GROUND_ENERGY = -0.8284271247461901
EXPECTED_DOUBLE_OCCUPANCY = 0.1464466094067262
EXPECTED_SPIN_CORRELATION = -0.640165042944955
SYMMETRY_BASIS = (5, 9, 6, 10)


def _two_level_rotation(state, a, b, theta):
    """
    This helper performs a rotation only between two computational basis states, \(|a\rangle\) and \(|b\rangle\).
    The leading underscore means this is an internal helper, not a pytest test. Pytest only automatically runs functions whose names begin with test_.
    """
    result = state.copy()
    amplitude_a, amplitude_b = state[a], state[b]
    result[a] = np.cos(theta) * amplitude_a - np.sin(theta) * amplitude_b
    result[b] = np.sin(theta) * amplitude_a + np.cos(theta) * amplitude_b
    return result


def _manual_variational_state(theta):
    """
    This helper constructs the manual three-parameter VQE ansatz.
    """
    state = np.zeros(16, dtype=complex)
    state[SYMMETRY_BASIS[0]] = 1.0 # state 5, |1010⟩
    ## three two-level rotations
    for target, angle in zip(SYMMETRY_BASIS[1:], theta, strict=True):
        state = _two_level_rotation(state, SYMMETRY_BASIS[0], target, angle)
    return state


@pytest.mark.slow
def test_manual_vqe_reproduces_reference_ground_state_and_observables():
    """
    To test that the manual ansatz and BFGS optimization can still reproduce the trusted ground-state energy, wavefunction, and observables.
    
    ---
    The decorator
    @pytest.mark.slow
    labels this as a slower integration test because it runs a numerical optimizer.
    It can be excluded with:
    pytest -m "not slow"
    """
    hamiltonian = hubbard_hamiltonian(t=T, u=U)
    pauli_terms = pauli_decomposition(hamiltonian)
    exact = solve_particle_sector(hamiltonian)

    def energy(theta):
        state = _manual_variational_state(theta)
        value = sum(
            coefficient * np.vdot(state, pauli_matrix(label) @ state)
            for label, coefficient in pauli_terms.items()
        )
        return float(np.real(value))

    result = minimize(
        energy,
        np.array([0.35, -0.25, 0.15]),
        method="BFGS",
        options={"gtol": 1e-7, "maxiter": 500},
    )
    vqe_state = _manual_variational_state(result.x)
    vqe_energy = energy(result.x)
    fidelity = abs(np.vdot(exact.ground_state, vqe_state)) ** 2

    assert vqe_energy >= exact.ground_energy - ATOL
    assert abs(vqe_energy - EXPECTED_GROUND_ENERGY) < VQE_ATOL
    assert fidelity > 1 - VQE_ATOL
    assert np.isclose(
        expectation(vqe_state, double_occupancy_operator()),
        EXPECTED_DOUBLE_OCCUPANCY,
        atol=VQE_ATOL,
    )
    assert np.isclose(
        expectation(vqe_state, spin_correlation_operator()),
        EXPECTED_SPIN_CORRELATION,
        atol=VQE_ATOL,
    )


@pytest.mark.slow
def test_seeded_qiskit_vqe_reproduces_reference_ground_state_and_observables():
    """
    
    """
    warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
    mapper = JordanWignerMapper()
    qubit_hamiltonian = mapper.map(qiskit_hubbard_fermionic_op(t=T, u=U)).simplify()
    exact = solve_particle_sector(hubbard_hamiltonian(t=T, u=U))

    num_particles = (1, 1)
    initial_state = HartreeFock(
        num_spatial_orbitals=N_SITES,
        num_particles=num_particles,
        qubit_mapper=mapper,
    )
    ansatz = UCCSD(
        num_spatial_orbitals=N_SITES,
        num_particles=num_particles,
        qubit_mapper=mapper,
        reps=2,
        initial_state=initial_state,
    )
    ## The random seed is fixed at 7. Without a fixed seed, every test run would start from different parameters and could follow a different optimization path. A fixed seed makes the regression reproducible. The initial parameters are random but deterministic
    initial_point = 0.2 * np.random.default_rng(7).standard_normal(
        ansatz.num_parameters
    )
    vqe = VQE(
        estimator=StatevectorEstimator(),
        ansatz=ansatz,
        optimizer=SLSQP(maxiter=1000, ftol=1e-12),
        initial_point=initial_point,
    )
    result = vqe.compute_minimum_eigenvalue(qubit_hamiltonian)

    optimal_circuit = ansatz.assign_parameters(result.optimal_parameters)
    vqe_state = Statevector.from_instruction(optimal_circuit).data
    vqe_energy = float(np.real(result.eigenvalue))
    fidelity = abs(np.vdot(exact.ground_state, vqe_state)) ** 2

    assert vqe_energy >= exact.ground_energy - ATOL
    assert abs(vqe_energy - EXPECTED_GROUND_ENERGY) < VQE_ATOL
    assert fidelity > 1 - VQE_ATOL
    ## This confirms that the Hartree–Fock plus UCCSD construction remained in the intended sector.
    assert np.isclose(
        expectation(vqe_state, particle_number_operator()), 2.0, atol=ATOL
    )
    assert np.isclose(
        expectation(vqe_state, double_occupancy_operator()),
        EXPECTED_DOUBLE_OCCUPANCY,
        atol=VQE_ATOL,
    )
    assert np.isclose(
        expectation(vqe_state, spin_correlation_operator()),
        EXPECTED_SPIN_CORRELATION,
        atol=VQE_ATOL,
    )
