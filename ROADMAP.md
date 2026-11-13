# Roadmap

Current focus: **Phase IIb, circuit-based unequal-time correlators**, on a
four-week schedule. Phases I and IIa are complete and regression-tested (see
the [README](README.md)). The earlier Phase I benchmark backlog is deferred and
lives in [`docs/phase1-backlog.md`](docs/phase1-backlog.md).

Status key: `[ ]` not started, `[~]` in progress, `[x]` done.

## Week 1: theory and exact reference

- [ ] Choose the first correlator, e.g. $\langle S^z_0(\tau)\,S^z_1(0)\rangle$
      after the $U_i/t=100\rightarrow U_f/t=5$ quench.
- [ ] Exact classical reference for the two-time correlator and its
      frequency-domain spectrum (`src/quantum_hubbard/correlators.py`).
- [ ] Derivation note: Hadamard-test circuit for
      $\langle\psi_0|U^\dagger(\tau)\,A\,U(\tau)\,B|\psi_0\rangle$ with Pauli
      $A,B$ controlled and the evolution $U$ uncontrolled; real/imaginary part
      extraction; measurement cost.
- [ ] Decide initial-state preparation on the circuit (`StatePreparation` from
      the ED vector first; Phase I VQE ansatz later if time allows).
- [ ] Tests for the exact reference (symmetries, $\tau=0$ limit, analytic
      dimer check).

## Week 2: circuit implementation

- [ ] Hadamard-test circuits built from the existing Suzuki-Trotter step.
- [ ] Noiseless statevector agreement with the exact correlator.
- [ ] Trotter-step convergence test for the correlator.
- [ ] Notebook `05_dynamical_correlators.ipynb` importing from `src/`.

## Week 3: finite shots and headline result

- [ ] Shot-based estimation; error vs. shots (expect roughly $1/\sqrt{N}$),
      repeated runs with recorded seeds.
- [ ] Headline figure: exact vs. circuit-measured correlator (and spectrum if
      time permits).
- [ ] Transpiled depth and two-qubit gate count for the dimer circuits;
      compare term ordering/grouping in the Trotter step.

## Week 4: packaging and applications

- [ ] README: short summary and headline figure at the top.
- [ ] Repository description, topics, and license.
- [ ] Short technical note (2-3 pages, PDF in `docs/`) with a clearly labeled
      "path to scale" discussion.
- [ ] CV bullets, cover-letter paragraph, interview preparation.

## After applying

- [ ] Four-site chain or 2x2 lattice: Trotter steps for fixed error,
      transpiled depth, and gate count vs. system size.
- [ ] Classical baseline with a tensor-network method (DMRG/TEBD).
- [ ] Noise model with one mitigation technique; optional hardware run.
- [ ] Selected items from the deferred Phase I backlog (ansatz and optimizer
      comparisons, mapping comparisons, symmetry tapering).

## Scope rules

- A small, fully validated result beats a broad, unvalidated one. If Week 2
  slips, narrow to one correlator at a few time points.
- Every new quantity gets an exact classical reference and a regression test
  before it is reported.
- All results are noiseless or finite-shot simulation unless stated otherwise.
  No quantum-advantage claims are made for this system size.
