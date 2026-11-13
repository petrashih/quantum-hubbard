# quantum-hubbard: handoff summary for Claude Code
The owner is learning and must be able to explain every derivation. For physics derivations, ask the owner for their derivation first and check it rather than supplying one. Write tests before implementation. Keep commits small.

Written 2026-11-12 from a planning conversation with the repo owner. Repo:
https://github.com/petrashih/quantum-hubbard (local: `~/Documents/project-QuantumComputing/quantum-hubbard`).

**Labels used throughout**
- **[IMPLEMENTED]**: exists in the repo, written by the owner. I (the planning assistant) read it but did not run it.
- **[PREPARED]**: files I generated in the conversation; the owner applied them locally; merge on GitHub was still pending.
- **[PROPOSED]**: discussed only. Nothing exists yet. Treat as a plan to confirm, not as a fact.
- **[UNVERIFIED]**: I did not check it.

## 1. Project goal and constraints

- Owner: PhD + postdoc in theoretical computational chemical physics (classical algorithms for materials simulation). New to quantum algorithms.
- Goals: (1) learn quantum algorithms by building something real; (2) have a CV project for **algorithm/application research roles** (chemistry, materials, optimization) in quantum computing industry.
- Timeline: **apply within about one month** of 2026-11-12. Available time: about 2 hours/day.
- Positioning: the PhD/postdoc is the headline credential; this repo is evidence of the transition. The project must **not claim quantum advantage** (4 qubits, exactly solvable).
- Working style: the owner wants to understand and be able to explain every part in interviews. Prefer small steps, explain design choices, and ask before large rewrites. One notebook is named `..._codex.ipynb`, so part of the repo was built with AI assistance; the owner should be able to explain all of it.

## 2. Physics and code conventions (from README and notebook 03)

- Model: open two-site Fermi-Hubbard dimer, H(U) = -t Σ_σ (c†_{0σ}c_{1σ} + h.c.) + U Σ_i n_{i↑}n_{i↓}. t = 1 is the energy unit; ħ = 1; time is in units of 1/t. Notebooks write time as τ to avoid clashing with t.
- Four spin orbitals in **spin-blocked order [0↑, 1↑, 0↓, 1↓]**, mapped to four qubits with **Jordan-Wigner**. Orbital 0 is the least-significant bit of the statevector index (Qiskit little-endian).
- Half filling: N↑ = N↓ = 1 (N = 2, six-dimensional sector inside the 16-dimensional Fock space).
- Reference energy: E₀ = (U − √(U²+16t²))/2 = −0.828427124746 at U = 4, t = 1.
- Observables: the repo's `double_occupancy_operator()` is **total** double occupancy D₀ + D₁. The paper plots **per-site** ⟨d⟩ = ½ΣD_i; the notebook keeps both. S^z_i = ½(n_{i↑} − n_{i↓}). Equal-time C^zz_ij = ⟨S^z_i S^z_j⟩.
- Quench protocol: prepare ground state of H(U_i), change U to U_f at τ = 0, evolve under H(U_f). Literature regression: U_i/t = 100 → U_f/t = 5 (Zavatti, Bellomia, Capone, arXiv:2605.18494 as cited in the README; I did not verify the reference). ω = √(U_f²+16t²) = √41, T = 2π/ω = 0.98126869.
- Analytic dimer states in notebook 03 use basis indices 5, 6, 9, 10 in this mode order.

## 3. Completed work [IMPLEMENTED]

Per the README, `src/quantum_hubbard/`, and notebooks. Tests: README says 24 tests (`pytest -m "not slow"` for fast; plain `pytest` includes seeded VQE integration tests).

- **Phase I (ground state)**: exact diagonalization, manual VQE (JW, Pauli decomposition, symmetry-preserving three-angle ansatz, SciPy optimizer), and Qiskit Nature VQE (UCCSD + Hartree-Fock). All three agree on energy, fidelity, double occupancy, and ⟨S₀·S₁⟩ at U = 4 (notebooks 01, 02).
- **Phase IIa (quench dynamics)**:
  - `dynamics.py` (read in full, 90 lines): `exact_trajectory` (eigendecomposition with input validation), `build_suzuki_trotter_step` (Qiskit `PauliEvolutionGate` + `SuzukiTrotter(order, reps)`, decomposed), `trotter_trajectory` (repeated `Statevector.evolve`).
  - Notebook 03: exact vs analytic fidelity 1, period and frequency reproduced, conservation checks (norm, energy, N, total S^z), Trotter vs exact, convergence at three step sizes with measured order ≈ 1.995 (maximum double-occupancy error over two periods).
  - Notebook 04: the same quench as a thin experiment importing from `src`.
  - Other modules: `model.py`, `operators.py`, `observables.py`.
- Packaging: `pyproject.toml` (Python ≥ 3.10, qiskit ~= 2.5.1), `requirements-tested.txt` (recorded on Python 3.14.6; not a full lock file).

## 4. Findings from reading the code [IMPLEMENTED facts, my assessment]

1. **Not end to end**: the initial state in notebook 03 comes from classical ED and is loaded directly into `Statevector`. The Phase I VQE state is not fed into the dynamics, and there is no state-preparation circuit.
2. **No transpilation or resource counts for dynamics**: notebook 03 prints only the decomposed one-step depth and operation counts. No transpiled depth or two-qubit gate count.
3. **Pauli-term ordering in the Trotter step** follows the dictionary order from `pauli_decomposition`. No grouping of commuting terms (hopping XX+YY, interaction ZZ).
4. **Convergence study checks observable error only** (max double-occupancy error, three step sizes).
5. **Duplication**: notebook 03 redefines `exact_trajectory` and `trotter_trajectory` locally while `dynamics.py` has the package versions (notebook 04 uses the package).
6. Notebook 03 closes by listing next steps: finite-shot estimation, post-transpilation resources, noisy simulation, then the paper's dephasing channel Γ > 0.

## 5. Repo cleanup [PREPARED, merge pending]

I generated these; the owner reported applying them on a local branch and was merging on GitHub. Confirm they landed.
- `ROADMAP.md`: four-week Phase IIb plan (see section 7).
- `docs/phase1-backlog.md`: old `TODO.md` moved with a "deferred" banner. The Phase I benchmark (about 70 unchecked boxes: ansatz/optimizer sweeps, mappings, tapering) is **deferred until after applying**.
- `README.md`: "Remaining work" replaced with a "Roadmap" section linking both files.
- `.github/workflows/tests.yml`: ruff + `pytest -m "not slow"` on Python 3.12. **[UNVERIFIED]**: never run. Tested environment was 3.14; the first CI run may fail on dependencies or on existing ruff warnings.
- Not done (GitHub settings and owner decisions): repo description and topics, license (MIT suggested; the owner should check postdoc/employer IP policy first), test badge in the README.

## 6. Technical proposals [PROPOSED only, nothing implemented]

**Phase IIb target**: unequal-time correlator C_AB(τ, 0) = ⟨ψ₀| U†(τ) A U(τ) B |ψ₀⟩ with U = e^{-iH_f τ}. The README defines the general C_AB(τ, τ′); because the quenched state is not stationary under H_f, the full two-time dependence must be kept unless a time is fixed by protocol. Suggested first correlator: ⟨S^z₀(τ) S^z₁(0)⟩, with an exact classical reference.

**Hadamard-test idea (to be derived and checked on paper first)**:
- If A and B are Pauli strings (or sums of them), only A and B need controlled versions; the evolution U can stay **uncontrolled**: ancilla in |+⟩, controlled-B, then U, then controlled-A, then measure the ancilla in X (real part) and Y (imaginary part).
- Under this JW convention n_i = (I − Z)/2, so S^z_i = (Z_{i↓} − Z_{i↑})/4 and S^z₀S^z₁ is a sum of four ZZ strings. Watch Qiskit's reversed Pauli-label order.
- The product A_H(τ)B is not Hermitian in general, so the correlator is complex.
- This is my derivation from memory; the owner should re-derive it and validate against the exact reference before relying on it.

**Other proposals**
- Initial-state preparation on the circuit: first `StatePreparation` from the ED vector (clearly labeled a stand-in), later the Phase I VQE ansatz with optimized parameters.
- New module `src/quantum_hubbard/correlators.py` with tests; notebook `05_dynamical_correlators.ipynb` importing from `src` (do not redefine functions in the notebook).
- Finite-shot error vs shots (expect roughly 1/√N), with recorded seeds and repeated runs.
- Transpiled depth and two-qubit gate count for the dimer; compare term ordering/grouping in the Trotter step; consider state-fidelity or unitary-error convergence in addition to observable error.
- After applying: 4-site chain or 2×2 lattice resource scaling (Trotter steps for fixed error, depth, gate count); DMRG/TEBD classical baseline (TeNPy or quimb); noise model plus one mitigation technique; optional hardware run; DMFT/impurity-solver context as discussion text only, clearly labeled as such.
- Technical note (2–3 page PDF in `docs/`), CV bullets, cover-letter paragraph, interview prep topics (JW locality, Trotter error scaling, Hadamard test, VQE limits and where quantum might beat DMRG/QMC). Resource-estimation papers (e.g., Wecker et al., Cade et al., Kivlichan et al.) were mentioned from memory; verify before citing.

## 7. Four-week plan [PROPOSED; written into ROADMAP.md]

1. **Week 1**: choose the first correlator; exact classical two-time reference and spectrum; paper derivation of the Hadamard circuit; decide state preparation; tests for the reference.
2. **Week 2**: Hadamard circuits from the existing Trotter step; noiseless agreement with exact; correlator convergence test; notebook 05.
3. **Week 3**: shot-based estimation and error scaling; headline figure (exact vs circuit-measured); transpiled depth and gate counts.
4. **Week 4**: README summary and headline figure, repo description/topics/license, technical note, CV bullets, start applying.
- Risk rule: if Week 2 slips, narrow to one correlator at a few time points. A small fully validated result beats a broad unvalidated one.

## 8. Unresolved questions

- Did the cleanup branch merge, and does CI pass on Python 3.12? Should the CI Python version match the tested 3.14?
- Which correlator first (spin ⟨S^z S^z⟩ vs density), and is a frequency-domain spectrum in scope (finite time window limits resolution; windowing needed)?
- How to prepare |ψ₀⟩ on the circuit for the first version (`StatePreparation` vs VQE ansatz)?
- License choice and IP policy check.
- Does the README's reference (arXiv:2605.18494) check out? Not verified.
- Coverage of the existing 24 tests and the contents of `model.py`, `operators.py`, `observables.py`, notebooks 01/02/04, and notebook output cells were **not reviewed** by me. Only the first part of each notebook-03 cell was printed, and the end of the Trotter explanation cell was cut off.
- Whether any of the Phase I backlog should be pulled forward (I suggested only its finite-shot and resource-count sections overlap with Weeks 3–4).

## 9. Suggested first actions for Claude Code

1. Confirm the cleanup branch merged; run `pytest -m "not slow"` and `ruff check src tests` locally and report results.
2. Read `model.py`, `operators.py`, `observables.py`, and the tests before adding code, to match existing conventions.
3. Start Week 1: write `correlators.py` with an exact two-time correlator reference and tests (τ = 0 limit equals the equal-time value, symmetries, analytic dimer check) **before** any circuit code.
4. Keep commits small, keep results validated against exact references, and keep notebooks as thin consumers of `src/`.
