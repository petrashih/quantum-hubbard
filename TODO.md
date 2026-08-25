# Phase I Validation and Benchmarking Backlog

This roadmap strengthens the two-site Hubbard ground-state study alongside
Phase II (dynamics). The core $U/t=4$ calculation is already validated; this
backlog turns it into a broader reproducible benchmark that explains not only
whether VQE works, but also when and why it becomes difficult.

## Definition of done for the Phase I benchmark

The broader Phase I benchmark is complete when the repository contains:

- [ ] an exact-diagonalization reference across a documented $U/t$ grid;
- [ ] a reproducible comparison of several ansätze and optimizers;
- [ ] finite-shot error results separated from optimization error;
- [ ] consistently transpiled circuit-resource measurements;
- [ ] documented qubit mappings and any symmetry reductions;
- [ ] plots and tables generated from saved benchmark data; and
- [ ] a short conclusions section stating what becomes harder as $U/t$ grows.

The core benchmark grid will be
$U/t \in \{0, 0.5, 1, 2, 4, 8, 16\}$ with $t=1$. Additional points can be
added near interesting transitions in behavior.

## 1. Establish the benchmark and metrics — highest priority

- [ ] Refactor the shared Hamiltonian, ED, observable, and validation code into
  reusable functions so all experiments use identical conventions.
- [ ] Fix and record random seeds, package versions, stopping tolerances, and
  maximum function evaluations.
- [ ] For every $U/t$, save the ED ground-state energy and state, total double
  occupancy, and $\langle\mathbf S_0\cdot\mathbf S_1\rangle$.
- [ ] For every VQE run, report the absolute energy error
  $|E_{\mathrm{VQE}}-E_{\mathrm{ED}}|$.
- [ ] Also record fidelity with the ED state, observable errors, optimizer
  evaluations, iterations, termination status, and wall-clock time.
- [ ] Save machine-readable results (CSV or JSON) rather than copying values
  manually into plots.
- [x] Add regression assertions for the existing $U/t=4$ result.

**Deliverable:** one baseline table and plots of energy, absolute energy error,
fidelity, double occupancy, and spin correlation versus $U/t$.

## 2. Understand and compare ansätze

### Short literature and documentation review

- [ ] Write a one-paragraph explanation of what an ansatz must provide:
  expressibility, symmetry preservation, shallow circuits, and a trainable
  energy landscape.
- [ ] Make a comparison table covering at least:
  - the current three-angle symmetry-preserving ansatz;
  - UCCSD with one and two repetitions;
  - a Hamiltonian variational ansatz (alternating hopping and interaction
    evolution);
  - a particle-number-preserving Givens-rotation ansatz; and
  - one adaptive method, such as ADAPT-VQE or qubit-ADAPT-VQE, as a stretch
    study.
- [ ] Briefly review hardware-efficient circuits as a contrast, including
  their shallow construction and their risk of leaving the target symmetry
  sector.
- [ ] For any method described as "state of the art," cite a recent primary
  source and state the problem setting in which it is competitive. Avoid
  implying that one ansatz is universally best.

### Controlled ansatz benchmark

- [ ] Use the same reference state and optimizer budget where the ansätze allow
  it.
- [ ] Run each ansatz from at least 20 recorded initial parameter sets at every
  $U/t$ value.
- [ ] Include a warm-start experiment, using the optimized parameters at one
  $U/t$ as the starting point for the next value.
- [ ] Report median and worst-case absolute energy error, median fidelity,
  success rate, and evaluation count—not only the best run.
- [ ] Define success before running the benchmark (for example,
  $|E_{\mathrm{VQE}}-E_{\mathrm{ED}}|<10^{-6}$ in the noiseless study).
- [ ] Separate **expressibility failure** (the ansatz cannot represent the
  target within the budget) from **optimization failure** (some starts reach
  the answer and others do not).

**Deliverable:** an ansatz comparison table and distributions of error and
fidelity versus $U/t$.

## 3. Compare classical optimizers

- [ ] First benchmark deterministic optimizers with exact statevector
  expectations: SLSQP, L-BFGS-B or BFGS, and COBYLA.
- [ ] Then benchmark an optimizer designed for stochastic objectives, such as
  SPSA, with finite-shot estimates.
- [ ] Use identical initial points and comparable evaluation budgets for fair
  pairwise comparisons.
- [ ] Repeat each ansatz–optimizer–$U/t$ combination over the same set of
  initializations.
- [ ] Plot convergence as energy error versus objective evaluations, since an
  "iteration" can cost different numbers of evaluations across optimizers.
- [ ] Report success rate, median final error, spread across initializations,
  and total circuit/shot cost.

**Deliverable:** a recommendation that is explicitly conditional on the
setting, such as noiseless simulation versus finite-shot estimation.

## 4. Quantify finite-shot error

### Measurement error at fixed parameters

- [ ] Freeze the ED-quality or noiseless-VQE parameters so optimizer behavior
  cannot contaminate the measurement study.
- [ ] Estimate energy and observables with
  $N_{\mathrm{shots}} \in \{100, 1{,}000, 10{,}000, 100{,}000\}$.
- [ ] Repeat each estimate at least 30 times with recorded seeds.
- [ ] Report bias, standard deviation, root-mean-square error, and confidence
  intervals relative to exact statevector values.
- [ ] Check the expected approximately $1/\sqrt{N_{\mathrm{shots}}}$
  statistical scaling.
- [ ] Compare term-by-term measurement with a documented commuting-group
  strategy and record the total shot budget.

### Optimization with finite shots

- [ ] Only after the fixed-parameter study, rerun VQE with shot-based energy
  estimates.
- [ ] Compare a deterministic optimizer with SPSA under the same total shot
  budget.
- [ ] Report the distribution of final errors, not a single favorable run.
- [ ] Distinguish sampling noise from device noise; hardware noise models are a
  separate stretch goal.

**Deliverable:** error-versus-shots plots for fixed-state estimation and for
full shot-based VQE.

## 5. Document mappings and symmetry reduction

- [ ] Document the spin-orbital order, bit significance, Pauli-label order,
  particle-number sector, and sign convention in one canonical location.
- [ ] Compare Jordan-Wigner, parity, and Bravyi-Kitaev mappings for the same
  fermionic operator.
- [ ] Verify matrix or spectrum agreement after accounting for each mapping's
  basis convention.
- [ ] Investigate two-qubit reduction or $\mathbb Z_2$ symmetry tapering and
  state exactly which fixed symmetries permit the reduction.
- [ ] Confirm that the reduced problem reproduces the target-sector ED energy
  and observables.
- [ ] Record qubit count and Pauli-term count before and after reduction.

**Deliverable:** a mapping/symmetry table with correctness checks and resource
savings.

## 6. Measure circuit resources consistently

- [ ] Record logical parameter count, undecomposed circuit depth, and logical
  two-qubit operations for every ansatz.
- [ ] Choose and document a common transpilation target: basis gates, coupling
  map or backend, optimization level, layout method, and transpiler seed.
- [ ] After transpilation, record total depth, two-qubit depth, two-qubit gate
  count, total gate count, and qubit count.
- [ ] Keep state-preparation cost separate from measurement-basis rotations.
- [ ] Report the number of Pauli terms, commuting measurement groups, circuits
  per energy evaluation, and total shots.
- [ ] Plot resources versus $U/t$ only when the selected ansatz or adaptive
  circuit actually changes with $U/t$. For a fixed circuit template, explain
  that depth and gate count remain constant even if optimization becomes
  harder.

**Deliverable:** a resource table that distinguishes logical, transpiled,
measurement, and optimization costs.

## 7. Answer the central physics/algorithm question

For each $U/t$, diagnose difficulty in four separate layers:

- [ ] **Representation:** What is the best error achievable by the ansatz?
- [ ] **Optimization:** How often does the optimizer find that solution, and
  how many energy evaluations are required?
- [ ] **Measurement:** How many shots are required for a target precision?
- [ ] **Circuit:** Does the required depth or two-qubit gate count change?

Relate these results to physical indicators of correlation, including double
occupancy, spin correlation, and the overlap of the chosen reference state
with the ED ground state. Do not infer that stronger correlation directly
causes every resource metric to increase; use the separate measurements above
to identify the actual bottleneck.

**Deliverable:** a concise conclusion answering: *As the system becomes more
correlated, which part of the VQE workflow becomes harder for each ansatz?*

## 8. Finalize Phase I

- [ ] Create publication-quality figures with labeled axes, units, legends,
  seeds, and uncertainty bars where appropriate.
- [ ] Update the README results and limitations from generated benchmark data.
- [x] Add an environment-lock file or otherwise record exact dependency
  versions used for the final results.
- [ ] Ensure all notebooks run from a fresh environment without hidden state.
- [x] Add automated tests for Hamiltonians, mappings, conserved quantities,
  and benchmark reference values.
- [ ] Write a short Phase I summary suitable for the project entry on a CV.

## Minimum scope for the broader Phase I benchmark

To keep Phase I bounded, the minimum milestone is Sections 1, 2 (three
ansätze), 3 (three deterministic optimizers plus SPSA), 4, and 6. Mapping
comparisons, symmetry tapering, adaptive ansätze, and hardware noise are useful
extensions. None of these backlog items blocks exploratory work on dynamics.
