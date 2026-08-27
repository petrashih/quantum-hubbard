"""Exact and product-formula real-time evolution utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]
ComplexVector = NDArray[np.complex128]


def exact_trajectory(
    hamiltonian: ComplexMatrix,
    initial_state: ComplexVector,
    times: NDArray[np.float64],
) -> ComplexMatrix:
    """Evolve a state exactly under a time-independent Hermitian Hamiltonian."""

    if hamiltonian.ndim != 2 or hamiltonian.shape[0] != hamiltonian.shape[1]:
        raise ValueError("hamiltonian must be a square matrix")
    dimension = hamiltonian.shape[0]
    if initial_state.ndim != 1 or initial_state.shape != (dimension,):
        raise ValueError("initial_state and hamiltonian dimensions do not match")
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or times.size == 0:
        raise ValueError("times must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(times)):
        raise ValueError("times must contain only finite values")
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise ValueError("hamiltonian must be Hermitian")

    energies, eigenvectors = np.linalg.eigh(hamiltonian)
    coefficients = eigenvectors.conj().T @ initial_state
    phases = np.exp(-1j * np.outer(times, energies))
    return (phases * coefficients) @ eigenvectors.T


def build_suzuki_trotter_step(
    qubit_hamiltonian,
    time_step: float,
    order: int = 2,
    reps: int = 1,
):
    """Build one decomposed Suzuki product-formula circuit for ``time_step``."""

    from qiskit import QuantumCircuit
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.synthesis import SuzukiTrotter

    if not np.isfinite(time_step) or time_step <= 0:
        raise ValueError("time_step must be finite and positive")
    if not isinstance(order, int) or order < 2 or order % 2:
        raise ValueError("order must be a positive even integer")
    if not isinstance(reps, int) or reps < 1:
        raise ValueError("reps must be a positive integer")
    if not hasattr(qubit_hamiltonian, "num_qubits"):
        raise TypeError("qubit_hamiltonian must expose num_qubits")

    product_formula = SuzukiTrotter(order=order, reps=reps)
    evolution_gate = PauliEvolutionGate(
        qubit_hamiltonian,
        time=time_step,
        synthesis=product_formula,
    )
    step_circuit = QuantumCircuit(qubit_hamiltonian.num_qubits)
    step_circuit.append(evolution_gate, range(qubit_hamiltonian.num_qubits))
    return step_circuit.decompose()


def trotter_trajectory(
    initial_state: ComplexVector,
    step_circuit,
    num_steps: int,
) -> ComplexMatrix:
    """Apply a decomposed product-formula step circuit repeatedly."""

    from qiskit.quantum_info import Statevector

    if not isinstance(num_steps, int) or num_steps < 0:
        raise ValueError("num_steps must be a non-negative integer")
    expected_dimension = 2**step_circuit.num_qubits
    if initial_state.ndim != 1 or initial_state.shape != (expected_dimension,):
        raise ValueError("initial_state and step_circuit dimensions do not match")

    state = Statevector(initial_state)
    states = [state.data.copy()]
    for _ in range(num_steps):
        state = state.evolve(step_circuit)
        states.append(state.data.copy())
    return np.asarray(states)
