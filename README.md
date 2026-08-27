# Two-Site Hubbard Model: Ground States and Quench Dynamics

An educational quantum-computing project that develops a validated simulation
workflow for the half-filled, two-site Fermi-Hubbard model. The project begins
with ground-state preparation and now extends to real-time evolution after an
interaction quench.

The goal is not to claim quantum advantage for a four-qubit system. The dimer
is small enough to solve analytically and by exact diagonalization, making it a
useful testbed for fermion-to-qubit mappings, variational algorithms, quantum
time evolution, observables, and scientific validation.

## Project phases

| Phase | Scope | Status |
| --- | --- | --- |
| Phase I | Ground states from exact diagonalization, manual VQE, and Qiskit Nature VQE | Implemented and regression-tested |
| Phase IIa | Sudden interaction quenches, exact propagation, analytic dimer dynamics, and Suzuki-Trotter circuits | Implemented and regression-tested |
| Phase IIb | Unequal-time dynamical correlation and response functions | Next |

Both VQE and real-time quantum calculations currently use noiseless statevector
simulation. Finite-shot sampling, noise models, and quantum-hardware execution
remain future work.

## Physical model and conventions

The project studies the open two-site Hubbard Hamiltonian

$$
H(U) = -t\sum_{\sigma \in \{\uparrow,\downarrow\}}
\left(c^\dagger_{0\sigma}c_{1\sigma}
+c^\dagger_{1\sigma}c_{0\sigma}\right)
+U\sum_{i=0}^{1}n_{i\uparrow}n_{i\downarrow}.
$$

There are four spin orbitals in the canonical spin-blocked order

$$
[0\uparrow,1\uparrow,0\downarrow,1\downarrow],
$$

mapped to four qubits with the Jordan-Wigner transformation. Calculations
target half filling with $N_\uparrow=N_\downarrow=1$ and use $t=1$ as the
energy unit, so time is expressed in units of $1/t$.

## Phase I: ground-state simulation

Phase I implements three independent solution paths:

| Method | What is implemented | Role |
| --- | --- | --- |
| Exact diagonalization | Explicit fermionic matrices, Hamiltonian construction, and diagonalization of the six-dimensional $N=2$ sector | Trusted classical reference |
| Manual VQE | Manual Jordan-Wigner operators, Pauli decomposition, a symmetry-preserving ansatz, exact Pauli expectations, and SciPy optimization | Exposes the complete variational workflow |
| Qiskit Nature VQE | `FermiHubbardModel`, mode-order alignment, `JordanWignerMapper`, Hartree-Fock preparation, UCCSD, and Qiskit's `VQE` | Standard-library reproduction |

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

## Phase IIa: quench dynamics

Phase IIa prepares the ground state $|\psi_0\rangle$ of an initial Hamiltonian
$H(U_i)$, suddenly changes the interaction to $U_f$ at time zero, and evolves
under the post-quench Hamiltonian:

$$
|\psi(\tau)\rangle=e^{-iH(U_f)\tau}|\psi_0\rangle.
$$

The implemented observables are the site-resolved and averaged double
occupancy,

$$
D_i(\tau)=\langle n_{i\uparrow}n_{i\downarrow}\rangle_\tau,
\qquad
\langle d(\tau)\rangle=\frac12\sum_iD_i(\tau),
$$

and the equal-time longitudinal spin correlation

$$
C_{ij}^{zz}(\tau)=\langle\psi(\tau)|S_i^zS_j^z|\psi(\tau)\rangle.
$$

The Phase IIa backbone includes:

- exact classical propagation under a time-independent Hamiltonian;
- the analytic two-singlet solution of the Hubbard dimer;
- reusable local and trajectory observable operators;
- decomposed second-order Suzuki-Trotter step circuits;
- repeated statevector circuit evolution;
- norm, energy, particle-number, and total-$S^z$ checks; and
- time-step convergence measurements against the exact trajectory.

The literature regression uses the closed-system $U_i/t=100\rightarrow U_f/t=5$
quench studied by
[Zavatti, Bellomia, and Capone](https://arxiv.org/abs/2605.18494). The dynamics
reproduce the analytical oscillation frequency and period

$$
\omega=E_+-E_-=\sqrt{U_f^2+16t^2}=\sqrt{41},
\qquad
T=\frac{2\pi}{\omega}=0.98126869.
$$

The analytic and exact trajectories have unit fidelity. The measured
Suzuki-Trotter convergence order is $1.9949$, consistent with the expected
global $O(\Delta\tau^2)$ error of a second-order formula.

## Phase IIb: dynamical correlation functions

Phase IIa measures an operator or equal-time product in the evolving state.
Phase IIb will instead introduce unequal-time correlations such as

$$
C_{AB}(\tau,\tau')=
\langle\psi_0|A_H(\tau)B_H(\tau')|\psi_0\rangle,
\qquad
A_H(\tau)=e^{iH_f\tau}Ae^{-iH_f\tau},
$$

and response functions such as

$$
\chi_{AB}^{R}(\tau,\tau')=
-i\,\Theta(\tau-\tau')
\langle[A_H(\tau),B_H(\tau')]\rangle.
$$

The first Phase IIb target will be a spin or density correlator with an exact
classical reference, followed by a quantum-circuit measurement strategy and,
where appropriate, a frequency-domain spectrum. Because a quenched initial
state is generally not stationary under $H_f$, the full two-time dependence
must be retained unless one time is fixed by the protocol.

## Notebooks

Recommended reading order:

### Phase I

1. [`01_two_site_hubbard_ps.ipynb`](notebook/01_two_site_hubbard_ps.ipynb) —
   first-principles exact diagonalization, Jordan-Wigner mapping, Pauli
   decomposition, and manual VQE.
2. [`02_two_site_hubbard_qiskit_nature.ipynb`](notebook/02_two_site_hubbard_qiskit_nature.ipynb) —
   the equivalent Qiskit Nature construction and UCCSD-VQE workflow.

[`01_two_site_hubbard_codex.ipynb`](notebook/01_two_site_hubbard_codex.ipynb)
is an alternate first-principles notebook using a site-interleaved orbital
ordering. It is retained to make basis and qubit-ordering conventions explicit.

### Phase IIa

3. [`03_quench_dynamics.ipynb`](notebook/03_quench_dynamics.ipynb) — derivation
   and first implementation of exact and Trotterized quench dynamics, including
   the literature regression and detailed Suzuki-Trotter explanation.
4. [`04_reusable_quench_experiment.ipynb`](notebook/04_reusable_quench_experiment.ipynb) —
   the same validated quench written as a thin experiment that imports the
   reusable physics operations from `src/quantum_hubbard`.

## Reusable implementation

```text
src/quantum_hubbard/
├── dynamics.py       # exact and Suzuki-Trotter real-time propagation
├── model.py          # Fock basis, Hamiltonian, ED, and analytic dimer states
├── operators.py      # Jordan-Wigner, Pauli, and Qiskit operator utilities
└── observables.py    # local, averaged, and trajectory observables

tests/
├── test_dynamics.py
├── test_qiskit_mapping.py
├── test_trotter_dynamics.py
├── test_two_site_reference.py
└── test_vqe_regression.py
```

The automated suite currently contains 24 tests. It protects the fermionic
algebra, fixed ground-state reference, qubit mappings, VQE results, exact
propagation, analytic singlet states, observable conventions, conservation
laws, the literature quench, decomposed Suzuki circuits, and second-order
Trotter convergence.

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
kernel. A clean **Run All** in each primary notebook is an executable acceptance
check.

Run the deterministic tests during normal development:

```bash
pytest -m "not slow"
```

Run the complete suite, including the seeded manual and Qiskit VQE integration
tests, before committing a physics or algorithm change:

```bash
pytest
```

The exact top-level versions used for the latest successful regression run are
recorded in [`requirements-tested.txt`](requirements-tested.txt). Package
version ranges remain in `pyproject.toml`.

## Remaining work

The immediate priority is Phase IIb. Phase I and IIa can also be extended with:

- sweeps over $U_i/t$, $U_f/t$, and evolution time;
- ansatz, optimizer, and initialization comparisons;
- finite-shot energy and observable estimation;
- transpiled circuit depth and two-qubit gate counts;
- noise models, mitigation, and hardware execution; and
- broader regression coverage across parameter sweeps.

The original Phase I experiments and completion criteria remain in the
[Phase I roadmap](TODO.md).

## References and implementation guides

- [Exact Hubbard-dimer quench study](https://arxiv.org/abs/2605.18494)
- [Qiskit real-time Trotterization tutorial](https://qiskit-community.github.io/qiskit-algorithms/tutorials/13_trotterQRTE.html)
- [OpenFermion/FQE Fermi-Hubbard tutorial](https://quantumai.google/openfermion/fqe/tutorials/fermi_hubbard)

## Skills demonstrated

- fermionic Fock-space construction and sign conventions;
- exact diagonalization in conserved-particle sectors;
- Jordan-Wigner mapping and Pauli decomposition;
- variational ground-state preparation;
- analytic and numerical real-time propagation;
- Suzuki-Trotter circuit synthesis and convergence analysis;
- local, equal-time, and trajectory observables; and
- independent scientific references and executable regression tests.
