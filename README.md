# Ground-State Simulation of the Two-Site Hubbard Model

An educational quantum-computing project that solves the half-filled,
two-site Fermi-Hubbard model in three independent ways:

1. exact diagonalization (ED),
2. a manual variational quantum eigensolver (VQE), and
3. Qiskit Nature's VQE workflow.

The goal is not to claim a quantum advantage for a four-qubit problem. It is to
build and validate the complete path from a fermionic many-body Hamiltonian to
a qubit representation and a hybrid quantum-classical ground-state algorithm.

> **Project status:** The Phase I ground-state core is implemented and protected
> by automated $U/t=4$ regression tests. Broader ansatz, optimizer, finite-shot,
> and resource comparisons remain planned validation work and do not block the
> start of Phase II, real-time dynamics.

## Physical problem

The project studies the open two-site Hubbard Hamiltonian

$$
H = -t\sum_{\sigma \in \{\uparrow,\downarrow\}}
\left(c^\dagger_{0\sigma}c_{1\sigma}
+c^\dagger_{1\sigma}c_{0\sigma}\right)
+U\sum_{i=0}^{1}n_{i\uparrow}n_{i\downarrow}.
$$

There are four spin orbitals, mapped to four qubits with the Jordan-Wigner
transformation. The calculations target the half-filled sector with two
electrons and use $t=1$ as the energy unit. The worked example uses $U=4$.

## Three solution paths

| Method | What is implemented | Why it is included |
| --- | --- | --- |
| Exact diagonalization | Fermionic creation and annihilation matrices, explicit Hamiltonian construction, and diagonalization of the six-dimensional $N=2$ sector | Provides an analytic and numerical reference |
| Manual VQE | Manual Jordan-Wigner operators, Pauli decomposition, a three-parameter symmetry-preserving ansatz, exact Pauli expectation values, and SciPy optimization | Exposes every layer hidden by a packaged VQE |
| Qiskit Nature VQE | `FermiHubbardModel`, mode-order alignment, `JordanWignerMapper`, Hartree-Fock reference, six-parameter UCCSD ansatz, and Qiskit's `VQE` | Reproduces the result with a standard quantum simulation workflow |

Both VQE implementations use exact statevector expectation values. They do not
yet include finite-shot sampling, a noise model, or execution on quantum
hardware.

## Current result

For $t=1$, $U=4$, and $N=2$, the analytic ground-state energy is

$$
E_0=\frac{U-\sqrt{U^2+16t^2}}{2}=-0.828427124746.
$$

| Quantity | ED | Manual VQE | Qiskit Nature VQE |
| --- | ---: | ---: | ---: |
| Ground-state energy | -0.8284271247 | -0.8284271247 | -0.8284271247 |
| State fidelity with ED | 1.0 | 1.0 | 1.0 |
| Total double occupancy | 0.14644661 | 0.14644661 | 0.14644661 |
| $\langle\mathbf S_0\cdot\mathbf S_1\rangle$ | -0.64016504 | -0.64016504 | -0.64016504 |

Agreement is checked with executable assertions, not only by comparing the
final energy. The notebooks verify:

- the canonical anticommutation relations and particle-number conservation;
- equality between the explicit fermionic and Jordan-Wigner matrices;
- every nonzero Pauli coefficient and the full $16\times16$ Hamiltonian;
- the full spectrum and the half-filled spectrum;
- the variational bound, ground-state fidelity, and particle number; and
- double occupancy and inter-site spin correlation.

## Notebooks

The recommended reading order is:

1. [`01_two_site_hubbard_ps.ipynb`](notebook/01_two_site_hubbard_ps.ipynb) —
   first-principles ED, Jordan-Wigner mapping, Pauli decomposition, and manual
   VQE using the spin-blocked orbital order
   $[0\uparrow,1\uparrow,0\downarrow,1\downarrow]$.
2. [`02_two_site_hubbard_qiskit_nature.ipynb`](notebook/02_two_site_hubbard_qiskit_nature.ipynb) —
   the equivalent Qiskit Nature construction and UCCSD-VQE calculation,
   cross-validated against the first notebook.
3. [`03_quench_dynamics.ipynb`](notebook/03_quench_dynamics.ipynb) — exact and
   Trotterized real-time evolution after an interaction quench, including the
   $U_i/t=100\rightarrow U_f/t=5$ Hubbard-dimer literature regression.
4. [`04_reusable_quench_experiment.ipynb`](notebook/04_reusable_quench_experiment.ipynb) —
   the same validated quench expressed as a thin experiment that imports
   analytic states, propagation, and observables from `src/quantum_hubbard`.

[`01_two_site_hubbard_codex.ipynb`](notebook/01_two_site_hubbard_codex.ipynb)
is an alternate first-principles version using a site-interleaved orbital
ordering. Keeping both versions makes the effect of basis and qubit-ordering
conventions explicit.

## Run locally

This project requires Python 3.10 or newer.

```bash
git clone https://github.com/petrashih/quantum-hubbard.git
cd quantum-hubbard
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
jupyter lab
```

Open the notebooks in the order above and select the environment's Python
kernel. A clean **Run All** in each primary notebook serves as the acceptance
test; a failed physical or numerical comparison raises an assertion.

The exact top-level versions used for the latest successful regression run are
recorded in [`requirements-tested.txt`](requirements-tested.txt). The version
ranges in `pyproject.toml` remain the normal installation source.

## Automated regression tests

The reusable implementation is organized as:

```text
src/quantum_hubbard/
├── dynamics.py       # exact and Suzuki-Trotter real-time propagation
├── model.py          # Fock basis, Hamiltonian, ED, and analytic dimer states
├── operators.py      # Jordan-Wigner, Pauli, and Qiskit operator utilities
└── observables.py    # local, averaged, and trajectory observables

tests/
├── test_two_site_reference.py
├── test_qiskit_mapping.py
└── test_vqe_regression.py
```

Run the deterministic physics and mapping tests during normal development:

```bash
pytest -m "not slow"
```

Run the complete suite, including the seeded manual and Qiskit VQE integration
tests, before committing a physics or algorithm change:

```bash
pytest
```

The fixed $t=1$, $U=4$, $N=2$ tests preserve the trusted spectrum, ground-state
energy, Pauli coefficients, double occupancy, spin correlation, particle
number, and VQE fidelity. Numerical comparisons use explicit tolerances rather
than exact floating-point equality.

## Skills demonstrated

- second-quantized fermionic Hamiltonians and Fock-space sign conventions;
- exact diagonalization in a conserved-particle-number sector;
- Jordan-Wigner fermion-to-qubit mapping and Pauli decomposition;
- symmetry-preserving variational ansatz design;
- hybrid quantum-classical optimization with VQE;
- Qiskit and Qiskit Nature model, mapper, ansatz, and estimator APIs; and
- scientific validation through independent implementations and observables.

## Next steps

Alongside the initial dynamics work, Phase I can be extended with:

- a sweep over interaction strength $U/t$;
- ansatz, optimizer, initialization, and convergence comparisons;
- finite-shot energy and observable estimation;
- circuit depth, parameter count, and measurement-cost reporting; and
- automated plots and broader regression coverage across parameter sweeps.

The prioritized experiments, metrics, and completion criteria are tracked in
the [Phase I roadmap](TODO.md).

## Phase II: quench dynamics

Phase II will prepare the ground state $|\psi_0\rangle$ of an initial
Hamiltonian $H_0$, suddenly change the interaction from $U_0$ to $U_1$, and
evolve under the post-quench Hamiltonian,

$$
|\psi(t)\rangle=e^{-iH_1t}|\psi_0\rangle.
$$

The first target is a noiseless two-site interaction quench that reports the
site-resolved double occupancy
$D_i(t)=\langle n_{i\uparrow}n_{i\downarrow}\rangle_t$ and equal-time spin
correlation $C_{ij}^{zz}(t)=\langle S_i^zS_j^z\rangle_t$. Validation has two
levels:

1. Exact classical propagation provides the pointwise reference for the
   Trotterized quantum simulation. Acceptance checks cover the $t=0$ values,
   state norm, post-quench energy, particle number, total $S^z$, and convergence
   as the Trotter time step is reduced.
2. A literature regression reproduces the closed-system ($\Gamma=0$)
   $U_i/t=100\rightarrow U_f/t=5$ quench in
   [Zavatti, Bellomia, and Capone](https://arxiv.org/abs/2605.18494), including
   their per-site double occupancy
   $\langle d(t)\rangle=\frac12\sum_i D_i(t)$. Its oscillation frequency must
   agree with $E_+-E_-=\sqrt{U_f^2+16t^2}$ in units with $\hbar=1$.

Implementation references include the
[Qiskit real-time Trotterization tutorial](https://qiskit-community.github.io/qiskit-algorithms/tutorials/13_trotterQRTE.html)
and the
[OpenFermion/FQE Fermi-Hubbard tutorial](https://quantumai.google/openfermion/fqe/tutorials/fermi_hubbard).
