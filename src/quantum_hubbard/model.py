"""First-principles model construction for the two-site Hubbard problem.

The canonical convention is the spin-blocked orbital order
``[0 up, 1 up, 0 down, 1 down]``. Orbital 0 is the least-significant bit in
the NumPy basis index.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

N_SITES = 2
N_SPIN_ORBITALS = 2 * N_SITES
HILBERT_DIMENSION = 2**N_SPIN_ORBITALS

SPIN_BLOCKED_ORBITALS = (
    (0, "up"),
    (1, "up"),
    (0, "down"),
    (1, "down"),
)
HOPPING_PAIRS = ((0, 1), (2, 3))
ONSITE_PAIRS = ((0, 2), (1, 3))

ComplexMatrix = NDArray[np.complex128]
ComplexVector = NDArray[np.complex128]


@dataclass(frozen=True)
class ParticleSectorSolution:
    """Eigenproblem solution embedded back into the full Fock space."""

    basis: tuple[int, ...]
    energies: NDArray[np.float64]
    eigenvectors: ComplexMatrix
    ground_state: ComplexVector

    @property
    def ground_energy(self) -> float:
        return float(self.energies[0])


def occupations(index: int, n_spin_orbitals: int = N_SPIN_ORBITALS) -> tuple[int, ...]:
    """Return ``(n_0, ..., n_(M-1))`` for a Fock-basis index."""

    if not 0 <= index < 2**n_spin_orbitals:
        raise ValueError("basis index is outside the Hilbert space")
    return tuple((index >> orbital) & 1 for orbital in range(n_spin_orbitals))


def ket_label(index: int, n_spin_orbitals: int = N_SPIN_ORBITALS) -> str:
    """Format a basis index in ascending orbital order."""

    return "|" + "".join(map(str, occupations(index, n_spin_orbitals))) + "⟩"


def basis_indices(
    num_particles: int, n_spin_orbitals: int = N_SPIN_ORBITALS
) -> tuple[int, ...]:
    """Return all basis indices in a fixed-particle-number sector."""

    if not 0 <= num_particles <= n_spin_orbitals:
        raise ValueError("num_particles must be between 0 and n_spin_orbitals")
    return tuple(
        index
        for index in range(2**n_spin_orbitals)
        if sum(occupations(index, n_spin_orbitals)) == num_particles
    )


def apply_annihilation(
    index: int,
    orbital: int,
    n_spin_orbitals: int = N_SPIN_ORBITALS,
) -> tuple[int, int] | None:
    """Apply a fermionic annihilator to a basis index.

    Returns ``(new_index, phase)`` or ``None`` when the orbital is empty.
    """

    if not 0 <= orbital < n_spin_orbitals:
        raise ValueError("orbital is outside the spin-orbital range")
    occupation = occupations(index, n_spin_orbitals)
    if occupation[orbital] == 0:
        return None
    phase = (-1) ** sum(occupation[:orbital])
    return index ^ (1 << orbital), phase


def annihilation_matrix(
    orbital: int, n_spin_orbitals: int = N_SPIN_ORBITALS
) -> ComplexMatrix:
    """Construct a fermionic annihilation matrix from the Fock-space rule."""

    dimension = 2**n_spin_orbitals
    matrix = np.zeros((dimension, dimension), dtype=complex)
    for column in range(dimension):
        result = apply_annihilation(column, orbital, n_spin_orbitals)
        if result is not None:
            row, phase = result
            matrix[row, column] = phase
    return matrix


def fermionic_operators(
    n_spin_orbitals: int = N_SPIN_ORBITALS,
) -> tuple[
    tuple[ComplexMatrix, ...], tuple[ComplexMatrix, ...], tuple[ComplexMatrix, ...]
]:
    """Return annihilation, creation, and number operators."""

    annihilators = tuple(
        annihilation_matrix(orbital, n_spin_orbitals)
        for orbital in range(n_spin_orbitals)
    )
    creators = tuple(operator.conj().T for operator in annihilators)
    number_operators = tuple(
        creators[orbital] @ annihilators[orbital] for orbital in range(n_spin_orbitals)
    )
    return annihilators, creators, number_operators


def hubbard_hamiltonian(t: float = 1.0, u: float = 4.0) -> ComplexMatrix:
    """Construct the open two-site Hubbard Hamiltonian in spin-blocked order."""

    annihilators, creators, number_operators = fermionic_operators()
    zero = np.zeros((HILBERT_DIMENSION, HILBERT_DIMENSION), dtype=complex)
    hopping = sum(
        (
            creators[left] @ annihilators[right] + creators[right] @ annihilators[left]
            for left, right in HOPPING_PAIRS
        ),
        start=zero,
    )
    interaction = sum(
        (number_operators[up] @ number_operators[down] for up, down in ONSITE_PAIRS),
        start=zero,
    )
    return -t * hopping + u * interaction


def particle_number_operator() -> ComplexMatrix:
    """Return the total particle-number operator."""

    _, _, number_operators = fermionic_operators()
    zero = np.zeros((HILBERT_DIMENSION, HILBERT_DIMENSION), dtype=complex)
    return sum(number_operators, start=zero)


def solve_particle_sector(
    hamiltonian: ComplexMatrix, num_particles: int = 2
) -> ParticleSectorSolution:
    """Diagonalize a Hamiltonian in a fixed-particle-number sector."""

    expected_shape = (HILBERT_DIMENSION, HILBERT_DIMENSION)
    if hamiltonian.shape != expected_shape:
        raise ValueError(f"hamiltonian must have shape {expected_shape}")

    basis = basis_indices(num_particles)
    sector_hamiltonian = hamiltonian[np.ix_(basis, basis)]
    energies, eigenvectors = np.linalg.eigh(sector_hamiltonian)
    ground_state = np.zeros(HILBERT_DIMENSION, dtype=complex)
    ground_state[list(basis)] = eigenvectors[:, 0]
    return ParticleSectorSolution(basis, energies, eigenvectors, ground_state)


def analytic_ground_energy(t: float = 1.0, u: float = 4.0) -> float:
    """Return the half-filled two-site Hubbard ground-state energy."""

    return float((u - np.sqrt(u**2 + 16 * t**2)) / 2)
