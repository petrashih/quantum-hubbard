"""Deterministic tests for exact quench dynamics and dimer references."""

import numpy as np
import pytest
from scipy.linalg import expm

from quantum_hubbard import (
    analytic_singlet_energies,
    analytic_singlet_state,
    average_double_occupancy_operator,
    dimer_quench_frequency,
    double_occupancy_operator,
    exact_trajectory,
    hubbard_hamiltonian,
    local_double_occupancy_operator,
    particle_number_operator,
    solve_particle_sector,
    spin_z_correlation_operator,
    spin_z_operator,
    total_spin_z_operator,
    trajectory_expectation,
)

T = 1.0
U_INITIAL = 100.0
U_FINAL = 5.0
ATOL = 1e-10
EXPECTED_INITIAL_DOUBLE_OCCUPANCY = 0.0001997603195526407
EXPECTED_MAX_DOUBLE_OCCUPANCY = 0.1854174990094452


def _literature_quench():
    initial_hamiltonian = hubbard_hamiltonian(t=T, u=U_INITIAL)
    final_hamiltonian = hubbard_hamiltonian(t=T, u=U_FINAL)
    initial_state = solve_particle_sector(initial_hamiltonian).ground_state
    return final_hamiltonian, initial_state


def test_exact_trajectory_matches_independent_matrix_exponential():
    hamiltonian, initial_state = _literature_quench()
    times = np.array([0.0, 0.13, 0.47, 0.91])

    actual = exact_trajectory(hamiltonian, initial_state, times)
    expected = np.asarray(
        [expm(-1j * hamiltonian * time) @ initial_state for time in times]
    )

    assert np.allclose(actual, expected, atol=ATOL)


def test_exact_trajectory_preserves_initial_state_and_conserved_quantities():
    hamiltonian, initial_state = _literature_quench()
    times = np.linspace(0.0, 2.0, 101)
    states = exact_trajectory(hamiltonian, initial_state, times)

    norms = np.sum(np.abs(states) ** 2, axis=1)
    energies = trajectory_expectation(states, hamiltonian)
    particle_numbers = trajectory_expectation(states, particle_number_operator())
    total_spin_z = trajectory_expectation(states, total_spin_z_operator())

    assert np.allclose(states[0], initial_state, atol=ATOL)
    assert np.allclose(norms, 1.0, atol=ATOL)
    assert np.allclose(energies, energies[0], atol=ATOL)
    assert np.allclose(particle_numbers, 2.0, atol=ATOL)
    assert np.allclose(total_spin_z, 0.0, atol=ATOL)


def test_no_quench_ground_state_only_accumulates_a_global_phase():
    hamiltonian = hubbard_hamiltonian(t=T, u=4.0)
    initial_state = solve_particle_sector(hamiltonian).ground_state
    times = np.linspace(0.0, 3.0, 31)
    states = exact_trajectory(hamiltonian, initial_state, times)

    fidelities = np.abs(states.conj() @ initial_state) ** 2
    double_occupancy = trajectory_expectation(
        states, average_double_occupancy_operator()
    )

    assert np.allclose(fidelities, 1.0, atol=ATOL)
    assert np.allclose(double_occupancy, double_occupancy[0], atol=ATOL)


@pytest.mark.parametrize("u", [-5.0, 0.0, 5.0, 100.0])
def test_analytic_singlets_are_orthonormal_eigenstates(u):
    hamiltonian = hubbard_hamiltonian(t=T, u=u)
    energy_minus, energy_plus = analytic_singlet_energies(t=T, u=u)
    ground = analytic_singlet_state(t=T, u=u, branch="ground")
    excited = analytic_singlet_state(t=T, u=u, branch="excited")
    exact_ground = solve_particle_sector(hamiltonian).ground_state

    assert np.isclose(np.linalg.norm(ground), 1.0, atol=ATOL)
    assert np.isclose(np.linalg.norm(excited), 1.0, atol=ATOL)
    assert np.isclose(np.vdot(ground, excited), 0.0, atol=ATOL)
    assert np.allclose(hamiltonian @ ground, energy_minus * ground, atol=ATOL)
    assert np.allclose(hamiltonian @ excited, energy_plus * excited, atol=ATOL)
    assert abs(np.vdot(ground, exact_ground)) ** 2 > 1 - ATOL


def test_literature_quench_matches_analytic_solution_and_period():
    hamiltonian, initial_state = _literature_quench()
    omega = dimer_quench_frequency(t=T, u_final=U_FINAL)
    period = 2 * np.pi / omega
    times = np.array([0.0, period / 2, period])
    exact_states = exact_trajectory(hamiltonian, initial_state, times)

    analytic_initial = analytic_singlet_state(t=T, u=U_INITIAL, branch="ground")
    final_ground = analytic_singlet_state(t=T, u=U_FINAL, branch="ground")
    final_excited = analytic_singlet_state(t=T, u=U_FINAL, branch="excited")
    energy_minus, energy_plus = analytic_singlet_energies(t=T, u=U_FINAL)
    alpha = np.vdot(final_ground, analytic_initial)
    beta = np.vdot(final_excited, analytic_initial)
    analytic_states = np.asarray(
        [
            alpha * np.exp(-1j * energy_minus * time) * final_ground
            + beta * np.exp(-1j * energy_plus * time) * final_excited
            for time in times
        ]
    )

    fidelities = (
        np.abs(np.einsum("ti,ti->t", analytic_states.conj(), exact_states)) ** 2
    )
    double_occupancy = trajectory_expectation(
        exact_states, average_double_occupancy_operator()
    )

    assert np.isclose(omega, np.sqrt(41.0), atol=ATOL)
    assert np.allclose(fidelities, 1.0, atol=ATOL)
    assert np.allclose(
        double_occupancy,
        [
            EXPECTED_INITIAL_DOUBLE_OCCUPANCY,
            EXPECTED_MAX_DOUBLE_OCCUPANCY,
            EXPECTED_INITIAL_DOUBLE_OCCUPANCY,
        ],
        atol=ATOL,
    )


def test_local_and_average_observable_definitions_are_consistent():
    double_0 = local_double_occupancy_operator(0)
    double_1 = local_double_occupancy_operator(1)
    total_double = double_occupancy_operator()
    average_double = average_double_occupancy_operator()
    spin_z_0 = spin_z_operator(0)
    spin_z_1 = spin_z_operator(1)

    assert np.allclose(total_double, double_0 + double_1, atol=ATOL)
    assert np.allclose(average_double, 0.5 * total_double, atol=ATOL)
    assert np.allclose(
        spin_z_correlation_operator(0, 1), spin_z_0 @ spin_z_1, atol=ATOL
    )
    assert np.allclose(total_spin_z_operator(), spin_z_0 + spin_z_1, atol=ATOL)


@pytest.mark.parametrize("site", [-1, 2, 0.5])
def test_local_observables_reject_invalid_site_indices(site):
    with pytest.raises(ValueError, match="site must be"):
        local_double_occupancy_operator(site)
    with pytest.raises(ValueError, match="site must be"):
        spin_z_operator(site)


def test_dynamics_helpers_reject_incompatible_inputs():
    hamiltonian = np.eye(2, dtype=complex)
    state = np.array([1.0, 0.0], dtype=complex)

    with pytest.raises(ValueError, match="square"):
        exact_trajectory(np.zeros((2, 3)), state, np.array([0.0]))
    with pytest.raises(ValueError, match="dimensions"):
        exact_trajectory(hamiltonian, np.ones(3), np.array([0.0]))
    with pytest.raises(ValueError, match="non-empty"):
        exact_trajectory(hamiltonian, state, np.array([]))
    with pytest.raises(ValueError, match="Hermitian"):
        exact_trajectory(
            np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex),
            state,
            np.array([0.0]),
        )
    with pytest.raises(ValueError, match="states must have shape"):
        trajectory_expectation(state, hamiltonian)
