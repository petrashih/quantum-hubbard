"""Tests for the exact two-time correlator reference."""

import numpy as np
import pytest
from scipy.linalg import expm

from quantum_hubbard import (
    exact_trajectory,
    hubbard_hamiltonian,
    solve_particle_sector,
    spin_z_correlation_operator,
    spin_z_operator,
    trajectory_expectation,
    two_time_correlator,
)

T = 1.0
U_INITIAL = 100.0
U_FINAL = 5.0
ATOL = 1e-10


def _literature_quench():
    final_hamiltonian = hubbard_hamiltonian(t=T, u=U_FINAL)
    initial_state = solve_particle_sector(
        hubbard_hamiltonian(t=T, u=U_INITIAL)
    ).ground_state
    return final_hamiltonian, initial_state


def _heisenberg(hamiltonian, operator, time):
    propagator = expm(-1j * hamiltonian * time)
    return propagator.conj().T @ operator @ propagator


@pytest.mark.parametrize("reference_time", [0.0, 0.37])
def test_correlator_matches_heisenberg_picture_definition(reference_time):
    hamiltonian, initial_state = _literature_quench()
    spin_0, spin_1 = spin_z_operator(0), spin_z_operator(1)
    times = np.array([0.0, 0.21, 0.58, 1.3])

    actual = two_time_correlator(
        hamiltonian, initial_state, spin_0, spin_1, times, reference_time
    )
    spin_1_heisenberg = _heisenberg(hamiltonian, spin_1, reference_time)
    expected = [
        np.vdot(
            initial_state,
            _heisenberg(hamiltonian, spin_0, time) @ spin_1_heisenberg @ initial_state,
        )
        for time in times
    ]

    assert np.allclose(actual, expected, atol=ATOL)


def test_correlator_matches_lehmann_representation():
    # C(tau, 0) = sum_mn <psi|m><m|A|n><n|B|psi> exp(i (E_m - E_n) tau)
    hamiltonian, initial_state = _literature_quench()
    spin_0, spin_1 = spin_z_operator(0), spin_z_operator(1)
    times = np.linspace(0.0, 2.0, 17)
    energies, vectors = np.linalg.eigh(hamiltonian)
    left = vectors.conj().T @ initial_state
    right = vectors.conj().T @ spin_1 @ initial_state
    a_matrix = vectors.conj().T @ spin_0 @ vectors
    phases = np.exp(
        1j * np.subtract.outer(energies, energies)[None] * times[:, None, None]
    )

    expected = np.einsum("m,mn,n,tmn->t", left.conj(), a_matrix, right, phases)
    actual = two_time_correlator(hamiltonian, initial_state, spin_0, spin_1, times)

    assert np.allclose(actual, expected, atol=ATOL)


def test_equal_times_reduce_to_equal_time_spin_correlation():
    hamiltonian, initial_state = _literature_quench()
    spin_0, spin_1 = spin_z_operator(0), spin_z_operator(1)
    times = np.linspace(0.0, 2.0, 9)
    states = exact_trajectory(hamiltonian, initial_state, times)
    equal_time = trajectory_expectation(states, spin_z_correlation_operator(0, 1))

    diagonal = [
        two_time_correlator(
            hamiltonian, initial_state, spin_0, spin_1, np.array([time]), time
        )[0]
        for time in times
    ]

    assert np.allclose(diagonal, equal_time, atol=ATOL)


def test_hermitian_conjugate_swaps_operators_and_times():
    # C_AB(tau, tau')^* = C_BA(tau', tau) for Hermitian A and B.
    hamiltonian, initial_state = _literature_quench()
    spin_0, spin_1 = spin_z_operator(0), spin_z_operator(1)
    tau, tau_prime = 0.83, 0.29

    c_ab = two_time_correlator(
        hamiltonian, initial_state, spin_0, spin_1, np.array([tau]), tau_prime
    )[0]
    c_ba = two_time_correlator(
        hamiltonian, initial_state, spin_1, spin_0, np.array([tau_prime]), tau
    )[0]

    assert np.isclose(np.conj(c_ab), c_ba, atol=ATOL)
    assert abs(c_ab.imag) > 1e-3  # genuinely complex, not a real observable


def test_stationary_state_depends_only_on_time_difference():
    hamiltonian = hubbard_hamiltonian(t=T, u=U_FINAL)
    ground_state = solve_particle_sector(hamiltonian).ground_state
    spin_0, spin_1 = spin_z_operator(0), spin_z_operator(1)
    lags = np.array([0.0, 0.4, 1.1])

    from_zero = two_time_correlator(hamiltonian, ground_state, spin_0, spin_1, lags)
    shifted = two_time_correlator(
        hamiltonian, ground_state, spin_0, spin_1, lags + 0.7, reference_time=0.7
    )

    assert np.allclose(from_zero, shifted, atol=ATOL)


def test_quenched_state_is_not_time_translation_invariant():
    hamiltonian, initial_state = _literature_quench()
    spin_0, spin_1 = spin_z_operator(0), spin_z_operator(1)
    lags = np.array([0.0, 0.4, 1.1])

    from_zero = two_time_correlator(hamiltonian, initial_state, spin_0, spin_1, lags)
    shifted = two_time_correlator(
        hamiltonian, initial_state, spin_0, spin_1, lags + 0.25, reference_time=0.25
    )

    assert not np.allclose(from_zero, shifted, atol=1e-3)


def test_correlator_rejects_incompatible_inputs():
    hamiltonian, initial_state = _literature_quench()
    spin_0 = spin_z_operator(0)
    times = np.array([0.0])

    with pytest.raises(ValueError, match="operator_a"):
        two_time_correlator(hamiltonian, initial_state, np.eye(2), spin_0, times)
    with pytest.raises(ValueError, match="operator_b"):
        two_time_correlator(hamiltonian, initial_state, spin_0, np.eye(2), times)
    with pytest.raises(ValueError, match="reference_time"):
        two_time_correlator(hamiltonian, initial_state, spin_0, spin_0, times, np.inf)
