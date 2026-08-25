"""Qubit and Qiskit operator utilities for the two-site Hubbard model."""

from __future__ import annotations

import itertools
from collections.abc import Mapping

import numpy as np
from numpy.typing import NDArray

from quantum_hubbard.model import N_SITES, N_SPIN_ORBITALS

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}

ComplexMatrix = NDArray[np.complex128]


def operator_on_qubit(
    single_qubit_operator: ComplexMatrix,
    qubit: int,
    n_qubits: int = N_SPIN_ORBITALS,
) -> ComplexMatrix:
    """Embed a single-qubit operator using little-endian qubit indexing."""

    if single_qubit_operator.shape != (2, 2):
        raise ValueError("single_qubit_operator must be a 2 x 2 matrix")
    if not 0 <= qubit < n_qubits:
        raise ValueError("qubit is outside the register")

    factors = [
        single_qubit_operator if index == qubit else I2
        for index in reversed(range(n_qubits))
    ]
    result = factors[0]
    for factor in factors[1:]:
        result = np.kron(result, factor)
    return result


def jordan_wigner_annihilator(
    orbital: int, n_qubits: int = N_SPIN_ORBITALS
) -> ComplexMatrix:
    """Return a Jordan-Wigner annihilator as a dense matrix."""

    if not 0 <= orbital < n_qubits:
        raise ValueError("orbital is outside the register")
    dimension = 2**n_qubits
    x_operator = operator_on_qubit(X, orbital, n_qubits)
    y_operator = operator_on_qubit(Y, orbital, n_qubits)
    parity = np.eye(dimension, dtype=complex)
    for preceding_orbital in range(orbital):
        parity = parity @ operator_on_qubit(Z, preceding_orbital, n_qubits)
    return parity @ (x_operator + 1j * y_operator) / 2


def pauli_matrix(label: str) -> ComplexMatrix:
    """Convert a big-endian Pauli label such as ``IIXX`` to a matrix."""

    if not label or any(symbol not in PAULI for symbol in label):
        raise ValueError("label must be a non-empty string over I, X, Y, and Z")
    result = PAULI[label[0]]
    for symbol in label[1:]:
        result = np.kron(result, PAULI[symbol])
    return result


def pauli_decomposition(
    operator: ComplexMatrix, atol: float = 1e-12
) -> dict[str, complex]:
    """Decompose a square power-of-two operator into Pauli strings."""

    if operator.ndim != 2 or operator.shape[0] != operator.shape[1]:
        raise ValueError("operator must be square")
    dimension = operator.shape[0]
    n_qubits = int(np.log2(dimension))
    if 2**n_qubits != dimension:
        raise ValueError("operator dimension must be a power of two")

    terms: dict[str, complex] = {}
    for symbols in itertools.product("IXYZ", repeat=n_qubits):
        label = "".join(symbols)
        matrix = pauli_matrix(label)
        coefficient = np.trace(matrix.conj().T @ operator) / dimension
        if abs(coefficient) > atol:
            terms[label] = complex(np.real_if_close(coefficient))
    return terms


def reconstruct_pauli_operator(terms: Mapping[str, complex]) -> ComplexMatrix:
    """Reconstruct a dense operator from Pauli coefficients."""

    if not terms:
        raise ValueError("terms cannot be empty")
    labels = tuple(terms)
    if len({len(label) for label in labels}) != 1:
        raise ValueError("all Pauli labels must have equal length")
    dimension = 2 ** len(labels[0])
    result = np.zeros((dimension, dimension), dtype=complex)
    for label, coefficient in terms.items():
        result += coefficient * pauli_matrix(label)
    return result


def qiskit_hubbard_fermionic_op(t: float = 1.0, u: float = 4.0):
    """Construct Qiskit Nature's Hubbard operator in spin-blocked order."""

    from qiskit_nature.second_q.hamiltonians import FermiHubbardModel
    from qiskit_nature.second_q.hamiltonians.lattices import (
        BoundaryCondition,
        LineLattice,
    )

    lattice = LineLattice(
        num_nodes=N_SITES,
        edge_parameter=-t,
        onsite_parameter=0.0,
        boundary_condition=BoundaryCondition.OPEN,
    )
    site_interleaved = FermiHubbardModel(
        lattice=lattice, onsite_interaction=u
    ).second_q_op()
    return site_interleaved.permute_indices([0, 2, 1, 3]).simplify()


def qiskit_fermionic_observables() -> dict[str, object]:
    """Return Qiskit Nature observables in the canonical spin-blocked order."""

    from qiskit_nature.second_q.operators import FermionicOp

    def term(label: str) -> object:
        return FermionicOp({label: 1.0}, num_spin_orbitals=N_SPIN_ORBITALS)

    number = [term(f"+_{p} -_{p}") for p in range(N_SPIN_ORBITALS)]
    spin_plus_0, spin_minus_0 = term("+_0 -_2"), term("+_2 -_0")
    spin_plus_1, spin_minus_1 = term("+_1 -_3"), term("+_3 -_1")

    double_occupancy = (number[0] @ number[2] + number[1] @ number[3]).simplify()
    spin_correlation = (
        0.25 * ((number[0] - number[2]) @ (number[1] - number[3]))
        + 0.5 * (spin_plus_0 @ spin_minus_1 + spin_minus_0 @ spin_plus_1)
    ).simplify()
    return {
        "double_occupancy": double_occupancy,
        "spin_correlation": spin_correlation,
    }


def qiskit_pauli_dict(operator, atol: float = 1e-12) -> dict[str, complex]:
    """Return nonzero coefficients from a Qiskit ``SparsePauliOp``."""

    return {
        label: complex(coefficient)
        for label, coefficient in operator.simplify(atol=atol).to_list()
        if abs(coefficient) > atol
    }
