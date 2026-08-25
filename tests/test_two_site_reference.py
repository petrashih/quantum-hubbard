"""
Deterministic regression tests for the trusted t=1, U=4 reference.
This is to test algebra, spectrum, energy, observables, and fixed Pauli coefficients for the trusted t=1, U=4 reference.
"""

import numpy as np

from quantum_hubbard import (
    HILBERT_DIMENSION,
    N_SPIN_ORBITALS,
    analytic_ground_energy,
    basis_indices,
    double_occupancy_operator,
    expectation,
    hubbard_hamiltonian,
    particle_number_operator,
    solve_particle_sector,
    spin_correlation_operator,
)
from quantum_hubbard.model import fermionic_operators
from quantum_hubbard.operators import (
    jordan_wigner_annihilator,
    pauli_decomposition,
    reconstruct_pauli_operator,
)

T = 1.0
U = 4.0
ATOL = 1e-10

EXPECTED_GROUND_ENERGY = -0.8284271247461901
EXPECTED_HALF_FILLED_BASIS = (3, 5, 6, 9, 10, 12)
EXPECTED_HALF_FILLED_SPECTRUM = np.array(
    [EXPECTED_GROUND_ENERGY, 0.0, 0.0, 0.0, 4.0, 4.82842712474619]
)
EXPECTED_DOUBLE_OCCUPANCY = 0.1464466094067262
EXPECTED_SPIN_CORRELATION = -0.640165042944955
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


def test_fermionic_operators_obey_canonical_anticommutation_relations():
    annihilators, creators, _ = fermionic_operators()
    identity = np.eye(HILBERT_DIMENSION)
    zero = np.zeros((HILBERT_DIMENSION, HILBERT_DIMENSION))

    for p in range(N_SPIN_ORBITALS):
        for q in range(N_SPIN_ORBITALS):
            expected = identity if p == q else zero
            anticommutator_creation = (
                annihilators[p] @ creators[q] + creators[q] @ annihilators[p]
            )
            anticommutator_annihilation = (
                annihilators[p] @ annihilators[q] + annihilators[q] @ annihilators[p]
            )
            assert np.allclose(anticommutator_creation, expected, atol=ATOL)
            assert np.allclose(anticommutator_annihilation, zero, atol=ATOL)


def test_half_filled_reference_energy_and_spectrum():
    """
    This does not test VQE. It protects the exact-diagonalization reference that VQE will later be compared against.
    """
    hamiltonian = hubbard_hamiltonian(t=T, u=U)
    number = particle_number_operator()
    solution = solve_particle_sector(hamiltonian, num_particles=2)

    ## check hamiltonian is hermitian
    assert np.allclose(hamiltonian, hamiltonian.conj().T, atol=ATOL)
    ## check particle number is conserved
    assert np.allclose(hamiltonian @ number - number @ hamiltonian, 0, atol=ATOL)
    ## check the half-filling basis convention
    assert basis_indices(2) == EXPECTED_HALF_FILLED_BASIS
    assert solution.basis == EXPECTED_HALF_FILLED_BASIS
    ## check the spectrum (six eigen values), the spectrum is solved by ED
    assert np.allclose(solution.energies, EXPECTED_HALF_FILLED_SPECTRUM, atol=ATOL)
    ## check that the numerical ED result matches the saved hard-coded value
    assert np.isclose(solution.ground_energy, EXPECTED_GROUND_ENERGY, atol=ATOL)
    ## check that the numerical ED results mateches the independent analytic formula
    assert np.isclose(
        solution.ground_energy, analytic_ground_energy(t=T, u=U), atol=ATOL
    )
    ## Confirm that the resulting state has two particles
    assert np.isclose(expectation(solution.ground_state, number), 2.0, atol=ATOL)


def test_half_filled_reference_observables():
    """
    checks physical properties of ground-state vector through double occupancy and spin correlation
    """
    solution = solve_particle_sector(hubbard_hamiltonian(t=T, u=U))

    double_occupancy = expectation(solution.ground_state, double_occupancy_operator())
    spin_correlation = expectation(solution.ground_state, spin_correlation_operator())

    assert np.isclose(double_occupancy, EXPECTED_DOUBLE_OCCUPANCY, atol=ATOL)
    assert np.isclose(spin_correlation, EXPECTED_SPIN_CORRELATION, atol=ATOL)


def test_manual_jordan_wigner_and_pauli_representations():
    """
    This test verifies encoding and representation correctness. It makes sure the Pauli Hamiltonian correctly represent the fermionic model.
    """
    annihilators, _, _ = fermionic_operators()
    hamiltonian = hubbard_hamiltonian(t=T, u=U)

    for orbital, annihilator in enumerate(annihilators):
        assert np.allclose(jordan_wigner_annihilator(orbital), annihilator, atol=ATOL)

    pauli_terms = pauli_decomposition(hamiltonian)
    assert pauli_terms.keys() == EXPECTED_PAULI_TERMS.keys()
    for label, expected_coefficient in EXPECTED_PAULI_TERMS.items():
        assert np.isclose(pauli_terms[label], expected_coefficient, atol=ATOL)
    assert np.allclose(reconstruct_pauli_operator(pauli_terms), hamiltonian, atol=ATOL)
