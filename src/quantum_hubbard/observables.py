"""
Physical observables for the canonical two-site Hubbard convention.
expectation values, double occupancy, and spin correlation
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from quantum_hubbard.model import N_SITES, fermionic_operators

ComplexMatrix = NDArray[np.complex128]
ComplexVector = NDArray[np.complex128]


def expectation(state: ComplexVector, operator: ComplexMatrix) -> float:
    """Return the real expectation value of a Hermitian operator."""

    if state.ndim != 1 or operator.shape != (state.size, state.size):
        raise ValueError("state and operator dimensions do not match")
    value = np.vdot(state, operator @ state)
    if not np.isclose(value.imag, 0.0, atol=1e-12):
        raise ValueError("expectation value has a non-negligible imaginary part")
    return float(value.real)


def trajectory_expectation(
    states: ComplexMatrix, operator: ComplexMatrix
) -> NDArray[np.float64]:
    """Return an observable expectation value for every state in a trajectory."""

    if states.ndim != 2:
        raise ValueError("states must have shape (num_times, state_dimension)")
    state_dimension = states.shape[1]
    if operator.shape != (state_dimension, state_dimension):
        raise ValueError("states and operator dimensions do not match")
    values = np.einsum("ti,ij,tj->t", states.conj(), operator, states)
    if not np.allclose(values.imag, 0.0, atol=1e-12):
        raise ValueError("expectation values have a non-negligible imaginary part")
    return values.real


def _validate_site(site: int) -> None:
    if not isinstance(site, int) or not 0 <= site < N_SITES:
        raise ValueError(f"site must be an integer between 0 and {N_SITES - 1}")


def local_double_occupancy_operator(site: int) -> ComplexMatrix:
    """Return ``n_(site,up) n_(site,down)``."""

    _validate_site(site)
    _, _, number = fermionic_operators()
    return number[site] @ number[N_SITES + site]


def double_occupancy_operator() -> ComplexMatrix:
    """Return total on-site double occupancy, summed over both sites."""

    return sum(
        (local_double_occupancy_operator(site) for site in range(N_SITES)),
        start=np.zeros((2 ** (2 * N_SITES), 2 ** (2 * N_SITES)), dtype=complex),
    )


def average_double_occupancy_operator() -> ComplexMatrix:
    """Return double occupancy averaged over lattice sites."""

    return double_occupancy_operator() / N_SITES


def spin_z_operator(site: int) -> ComplexMatrix:
    """Return the local spin projection ``S_site^z``."""

    _validate_site(site)
    _, _, number = fermionic_operators()
    return 0.5 * (number[site] - number[N_SITES + site])


def total_spin_z_operator() -> ComplexMatrix:
    """Return total spin projection summed over all sites."""

    zero = np.zeros((2 ** (2 * N_SITES), 2 ** (2 * N_SITES)), dtype=complex)
    return sum((spin_z_operator(site) for site in range(N_SITES)), start=zero)


def spin_z_correlation_operator(site_i: int, site_j: int) -> ComplexMatrix:
    """Return the longitudinal equal-time operator ``S_i^z S_j^z``."""

    _validate_site(site_i)
    _validate_site(site_j)
    return spin_z_operator(site_i) @ spin_z_operator(site_j)


def spin_correlation_operator() -> ComplexMatrix:
    """Return the inter-site spin correlation ``S_0 dot S_1``."""

    annihilators, creators, number = fermionic_operators()
    spin_z_0 = 0.5 * (number[0] - number[2])
    spin_z_1 = 0.5 * (number[1] - number[3])
    spin_plus_0 = creators[0] @ annihilators[2]
    spin_minus_0 = creators[2] @ annihilators[0]
    spin_plus_1 = creators[1] @ annihilators[3]
    spin_minus_1 = creators[3] @ annihilators[1]
    return spin_z_0 @ spin_z_1 + 0.5 * (
        spin_plus_0 @ spin_minus_1 + spin_minus_0 @ spin_plus_1
    )


def observable_operators() -> dict[str, ComplexMatrix]:
    """Return the ground-state observables used by the benchmark."""

    return {
        "double_occupancy": double_occupancy_operator(),
        "spin_correlation": spin_correlation_operator(),
    }
