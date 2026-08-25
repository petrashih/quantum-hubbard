"""
Regression tests for the Qiskit Nature representation and mapping.
This file checks that Qiskit Nature translates the fermionic problem into the same qubit problem as the independent manual implementation.

Qiskit fermionic operator
        ↓
Jordan–Wigner mapping              test_qiskit_mapping.py
        ↓
Qubit Hamiltonian


"""

import numpy as np
from qiskit_nature.second_q.mappers import JordanWignerMapper

from quantum_hubbard import (
    double_occupancy_operator,
    hubbard_hamiltonian,
    spin_correlation_operator,
)
from quantum_hubbard.operators import (
    pauli_decomposition,
    qiskit_fermionic_observables,
    qiskit_hubbard_fermionic_op,
    qiskit_pauli_dict,
)

T = 1.0
U = 4.0
ATOL = 1e-10

EXPECTED_PAULI_TERMS = {
    "IIII": 2.0,
    "IIIZ": -1.0,
    "IIXX": -0.5,
    "IIYY": -0.5,
    "IIZI": -1.0,
    "IZII": -1.0,
    "IZIZ": 1.0,
    "XXII": -0.5,
    "YYII": -0.5,
    "ZIII": -1.0,
    "ZIZI": 1.0,
}


def test_qiskit_jordan_wigner_mapping_matches_fixed_reference():
    ## constructing the trusted manual Hamiltonian
    reference_matrix = hubbard_hamiltonian(t=T, u=U)
    ## Create Qiskit Jordan-Wigner mapper
    mapper = JordanWignerMapper()
    ## Qiskit Nature constructs the fermionic Hamiltonian
    fermionic_operator = qiskit_hubbard_fermionic_op(t=T, u=U)
    ## The fermionic operator is then mapped to qubits
    qubit_operator = mapper.map(fermionic_operator).simplify()
    ## extracts Pauli dictionaries from two independent paths
    qiskit_terms = qiskit_pauli_dict(qubit_operator)
    manual_terms = pauli_decomposition(reference_matrix)

    assert qiskit_terms.keys() == EXPECTED_PAULI_TERMS.keys()
    assert manual_terms.keys() == EXPECTED_PAULI_TERMS.keys()
    for label, expected_coefficient in EXPECTED_PAULI_TERMS.items():
        assert np.isclose(qiskit_terms[label], expected_coefficient, atol=ATOL)
        assert np.isclose(manual_terms[label], expected_coefficient, atol=ATOL)

    qiskit_matrix = qubit_operator.to_matrix()
    assert np.allclose(qiskit_matrix, reference_matrix, atol=ATOL)
    assert np.allclose(
        np.linalg.eigvalsh(qiskit_matrix),
        np.linalg.eigvalsh(reference_matrix),
        atol=ATOL,
    )


def test_qiskit_observable_mappings_match_first_principles_matrices():
    """
    This function checks physical observables rather than the Hamiltonian.
    """
    mapper = JordanWignerMapper()
    qiskit_observables = qiskit_fermionic_observables()

    mapped_double_occupancy = mapper.map(
        qiskit_observables["double_occupancy"]
    ).simplify()
    mapped_spin_correlation = mapper.map(
        qiskit_observables["spin_correlation"]
    ).simplify()

    assert np.allclose(
        mapped_double_occupancy.to_matrix(),
        double_occupancy_operator(),
        atol=ATOL,
    )
    assert np.allclose(
        mapped_spin_correlation.to_matrix(),
        spin_correlation_operator(),
        atol=ATOL,
    )
