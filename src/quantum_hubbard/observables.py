"""
Physical observables for the canonical two-site Hubbard convention.
expectation values, double occupancy, and spin correlation
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from quantum_hubbard.model import fermionic_operators

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


def double_occupancy_operator() -> ComplexMatrix:
    """Return total on-site double occupancy, summed over both sites."""

    _, _, number = fermionic_operators()
    return number[0] @ number[2] + number[1] @ number[3]


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
