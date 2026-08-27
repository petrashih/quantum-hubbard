"""Regression tests for decomposed Suzuki-Trotter quench evolution."""

import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp

from quantum_hubbard import (
    average_double_occupancy_operator,
    build_suzuki_trotter_step,
    dimer_quench_frequency,
    exact_trajectory,
    hubbard_hamiltonian,
    particle_number_operator,
    solve_particle_sector,
    total_spin_z_operator,
    trajectory_expectation,
    trotter_trajectory,
)
from quantum_hubbard.operators import pauli_decomposition

T = 1.0
U_INITIAL = 100.0
U_FINAL = 5.0
ATOL = 1e-10


def _quench_components():
    initial_hamiltonian = hubbard_hamiltonian(t=T, u=U_INITIAL)
    final_hamiltonian = hubbard_hamiltonian(t=T, u=U_FINAL)
    initial_state = solve_particle_sector(initial_hamiltonian).ground_state
    pauli_terms = pauli_decomposition(final_hamiltonian)
    qubit_hamiltonian = SparsePauliOp.from_list(list(pauli_terms.items()))
    return final_hamiltonian, qubit_hamiltonian, initial_state


def test_suzuki_step_is_a_decomposed_product_formula_circuit():
    _, qubit_hamiltonian, _ = _quench_components()
    period = 2 * np.pi / dimer_quench_frequency(t=T, u_final=U_FINAL)
    step = build_suzuki_trotter_step(qubit_hamiltonian, time_step=period / 20, order=2)
    operation_names = set(step.count_ops())

    assert {"rz", "rxx", "ryy", "rzz"} <= operation_names
    assert all("PauliEvolution" not in name for name in operation_names)


def test_second_order_trotter_trajectory_is_unitary_and_convergent():
    final_hamiltonian, qubit_hamiltonian, initial_state = _quench_components()
    period = 2 * np.pi / dimer_quench_frequency(t=T, u_final=U_FINAL)
    double_occupancy = average_double_occupancy_operator()
    errors = []
    finest_states = None

    for steps_per_period in (20, 40, 80):
        time_step = period / steps_per_period
        num_steps = 2 * steps_per_period
        times = np.arange(num_steps + 1) * time_step
        exact_states = exact_trajectory(final_hamiltonian, initial_state, times)
        step = build_suzuki_trotter_step(
            qubit_hamiltonian, time_step=time_step, order=2
        )
        trotter_states = trotter_trajectory(initial_state, step, num_steps)
        exact_values = trajectory_expectation(exact_states, double_occupancy)
        trotter_values = trajectory_expectation(trotter_states, double_occupancy)

        errors.append(np.max(np.abs(trotter_values - exact_values)))
        assert np.allclose(np.sum(np.abs(trotter_states) ** 2, axis=1), 1.0, atol=ATOL)
        finest_states = trotter_states

    errors = np.asarray(errors)
    error_reduction = errors[:-1] / errors[1:]

    assert errors[0] > 1e-4
    assert np.all((3.5 < error_reduction) & (error_reduction < 4.5))
    assert errors[-1] < 1.1e-4
    assert finest_states is not None
    assert np.allclose(
        trajectory_expectation(finest_states, particle_number_operator()),
        2.0,
        atol=ATOL,
    )
    assert np.allclose(
        trajectory_expectation(finest_states, total_spin_z_operator()),
        0.0,
        atol=ATOL,
    )


def test_trotter_helpers_reject_invalid_inputs():
    _, qubit_hamiltonian, initial_state = _quench_components()
    valid_step = build_suzuki_trotter_step(qubit_hamiltonian, time_step=0.01)

    with pytest.raises(ValueError, match="time_step"):
        build_suzuki_trotter_step(qubit_hamiltonian, time_step=0.0)
    with pytest.raises(ValueError, match="order"):
        build_suzuki_trotter_step(qubit_hamiltonian, time_step=0.01, order=3)
    with pytest.raises(ValueError, match="reps"):
        build_suzuki_trotter_step(qubit_hamiltonian, time_step=0.01, reps=0)
    with pytest.raises(ValueError, match="num_steps"):
        trotter_trajectory(initial_state, valid_step, num_steps=-1)
    with pytest.raises(ValueError, match="dimensions"):
        trotter_trajectory(np.ones(8, dtype=complex), valid_step, num_steps=1)
