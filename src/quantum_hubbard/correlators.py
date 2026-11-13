"""Exact classical references for unequal-time (two-time) correlators.

For a state ``|psi_0>`` prepared at ``tau = 0`` and evolved by a
time-independent Hamiltonian ``H`` (the post-quench Hamiltonian), the
two-time correlator is

    C_AB(tau, tau') = <psi_0| A_H(tau) B_H(tau') |psi_0>,
    A_H(tau) = U^dagger(tau) A U(tau),  U(tau) = exp(-i H tau).

Inserting ``U U^dagger = 1`` gives the form used here,

    C_AB(tau, tau') = <psi(tau)| A U(tau - tau') B |psi(tau')>,

so only state propagation is needed (no operator exponentials). After a
quench ``|psi_0>`` is not an eigenstate of ``H``, so ``C_AB`` depends on both
times, not only on ``tau - tau'``. ``A_H(tau) B_H(tau')`` is not Hermitian in
general, so ``C_AB`` is complex.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from quantum_hubbard.dynamics import exact_trajectory

ComplexMatrix = NDArray[np.complex128]
ComplexVector = NDArray[np.complex128]


def two_time_correlator(
    hamiltonian: ComplexMatrix,
    initial_state: ComplexVector,
    operator_a: ComplexMatrix,
    operator_b: ComplexMatrix,
    times: NDArray[np.float64],
    reference_time: float = 0.0,
) -> NDArray[np.complex128]:
    """Return ``C_AB(tau, reference_time)`` for every ``tau`` in ``times``.

    ``operator_a`` acts at the later time ``tau`` and ``operator_b`` at
    ``reference_time``. Neither time ordering nor a connected part is applied.
    """

    dimension = initial_state.shape[0] if initial_state.ndim == 1 else -1
    for name, operator in (("operator_a", operator_a), ("operator_b", operator_b)):
        if operator.shape != (dimension, dimension):
            raise ValueError(f"{name} and initial_state dimensions do not match")
    if not np.isfinite(reference_time):
        raise ValueError("reference_time must be finite")
    times = np.asarray(times, dtype=float)

    # |psi(tau)> for every tau.
    bra_states = exact_trajectory(hamiltonian, initial_state, times)
    # B |psi(tau')>, then propagated by U(tau - tau').
    psi_reference = exact_trajectory(
        hamiltonian, initial_state, np.array([reference_time])
    )[0]
    ket_states = exact_trajectory(
        hamiltonian, operator_b @ psi_reference, times - reference_time
    )
    return np.einsum("ti,ij,tj->t", bra_states.conj(), operator_a, ket_states)
