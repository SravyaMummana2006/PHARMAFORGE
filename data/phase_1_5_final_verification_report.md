# PHARMAFORGE — Phase 1.5: Final Research Integrity & Numerical Consistency Report

**Date:** November 21, 2026  
**Status:** AUDITED, RECONCILED, RECALCULATED & SCIENTIFICALLY CONSTRAINED  
**Project:** PHARMAFORGE Research Prototype
**Research Focus:** *Self-Improving Multi-Agent Framework for Scientific Reasoning and Pharmaceutical Discovery*  

---

## 1. Hypervolume Methodological Reconciliation & QED Normalization Fix

### A. Correction of QED Normalization
In [`src/pharmaforge/benchmarks/metrics.py`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/src/pharmaforge/benchmarks/metrics.py#L50-L65), the multi-objective hypervolume normalizes each objective value into the unit interval $[0, 1]$.

For maximization objectives, the previous implementation computed:
$$\tilde{y} = \max\left(0.0, \, \min\left(1.0, \, \frac{\text{val} - \text{ref\_val}}{10.0 - \text{ref\_val}}\right)\right)$$

While the divisor $(10.0 - \text{ref\_val})$ is suitable for $p\text{IC}_{50}$ (where $\text{ref} = 4.0, \, \text{max} \approx 10.0$), applying it uniformly to QED ($\text{val} \in [0, 1], \, \text{ref\_val} = 0.0$) compressed QED into $[0.0, 0.1]$. Because hypervolume in 4 dimensions is a product of normalized spans, this compressed QED scale artificially shrank all 4D hypervolumes by a factor of 10.

**Code Correction in `metrics.py`:**
```python
if obj.name.lower() == "qed":
    # QED is naturally bounded in [0.0, 1.0]
    span = max(1.0 - ref_val, 1e-6)
else:
    span = max(10.0 - ref_val, 1e-6)
norm_val = (val - ref_val) / span
```
This correctly scales QED to $[0.0, 1.0]$.

### B. Monte Carlo Numerical Integration Protocol
Hypervolume is computed using **Monte Carlo numerical integration** with $N = 20,000$ pseudo-random samples drawn from $\mathcal{U}(0, 1)^4$ using a deterministic generator (`np.random.default_rng(42)`). 

Hypervolume values are reported as Monte Carlo estimates ($\text{HV}_{\text{MC}}$) with their standard errors:
$$\text{SE} = \sqrt{\frac{\text{HV}_{\text{MC}} \cdot (1 - \text{HV}_{\text{MC}})}{N}}$$

### C. Recalculated 4-Way Ablation Hypervolumes
All configurations were recalculated using the corrected QED normalization, the fixed 4D reference point $r = [\text{affinity}: 4.0, \, \text{qed}: 0.0, \, \text{logp}: 6.0, \, \text{sascore}: 10.0]$, and the exact same $N = 20,000$ Monte Carlo sampling protocol. Previous values are preserved in [`artifacts/benchmark_runs/benchmark_optimization_results.json`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/artifacts/benchmark_runs/benchmark_optimization_results.json) under `previous_pareto_hypervolume_scaled10`.

| Benchmark Lead | Configuration | Pareto Front Size | Corrected $\text{HV}_{\text{MC}}$ | Monte Carlo $\text{SE}$ | 95% Confidence Interval | Previous Scaled $\text{HV}_{\text{scaled10}}$ |
|---|---|:---:|:---:|:---:|:---:|:---:|
| **`CHEMBL289959`** | **Config A** (No Memory, No Repair) | 7 | **0.1578** | $\pm 0.0026$ | $[0.1528, \, 0.1629]$ | 0.0159 |
| ($R = \text{H}$, unsubstituted) | **Config B** (Unverified Memory) | 1 | **0.1070** | $\pm 0.0022$ | $[0.1027, \, 0.1113]$ | 0.0109 |
| | **Config C** (Verified Memory, No Repair) | 1 | **0.1070** | $\pm 0.0022$ | $[0.1027, \, 0.1113]$ | 0.0109 |
| | **Config D** (Full PHARMAFORGE) | 7 | **0.1578** | $\pm 0.0026$ | $[0.1528, \, 0.1629]$ | 0.0159 |
| **`CHEMBL290096`** | **Config A** (No Memory, No Repair) | 5 | **0.1604** | $\pm 0.0026$ | $[0.1553, \, 0.1654]$ | 0.0162 |
| ($R = 3-\text{Br}$, 3-bromophenyl) | **Config B** (Unverified Memory) | 1 | **0.1215** | $\pm 0.0023$ | $[0.1170, \, 0.1260]$ | 0.0124 |
| | **Config C** (Verified Memory, No Repair) | 1 | **0.1215** | $\pm 0.0023$ | $[0.1170, \, 0.1260]$ | 0.0124 |
| | **Config D** (Full PHARMAFORGE) | 5 | **0.1604** | $\pm 0.0026$ | $[0.1553, \, 0.1654]$ | 0.0162 |

*Analytical Validation:* For the single-point baseline front (Configs B & C, 1 point: $p\text{IC}_{50}=6.419, \, \text{QED}=0.613, \, \text{LogP}=3.373, \, \text{SAScore}=1.815$ on `CHEMBL289959`), the analytical 4D box volume is $0.1084$. For Configs A & D, the exact analytical volume via principle of inclusion-exclusion is $0.1569$ (for `CHEMBL289959`) and $0.1591$ (for `CHEMBL290096`). The Monte Carlo numerical estimates ($0.1578$ and $0.1604$) match the analytical values within $1\text{ SE}$ ($\pm 0.0026$).

---

## 2. LogP Constraint Non-Compliance & Objective Definition Reconciliation

### A. Identification of Near-Threshold Candidate Non-Compliance
The stated drug-likeness constraint for lead optimization is:
$$\text{LogP} \le 3.500$$

The generated pyrido-pyrimidine isosteres have the following properties:
- `Brc1cccc(Nc2ncnc3cnccc23)c1` (Pyrido[4,3-d]pyrimidine): $\text{LogP} = \mathbf{3.531}$, Pred $p\text{IC}_{50} = 7.389$, $\text{QED} = 0.786$, $\text{SAScore} = 2.115$.
- `Brc1cccc(Nc2ncnc3ccncc23)c1` (Pyrido[3,4-d]pyrimidine): $\text{LogP} = \mathbf{3.531}$, Pred $p\text{IC}_{50} = 7.458$, $\text{QED} = 0.786$, $\text{SAScore} = 2.126$.
- `Brc1cccc(Nc2ncnc3cccnc23)c1` (Pyrido[3,2-d]pyrimidine): $\text{LogP} = \mathbf{3.531}$, Pred $p\text{IC}_{50} = 7.334$, $\text{QED} = 0.786$, $\text{SAScore} = 2.062$.

**Strict Finding:**
Because $3.531 > 3.500$, these candidates **fail the stated $\text{LogP} \le 3.500$ constraint** by $+0.031\log$ units.

**Scientific Clarification:**
1. **Why they appear on the Pareto front:** In multi-objective optimization, the non-dominated Pareto front consists of solutions that are not dominated by any other observed candidate across the continuous objective vectors $(\text{affinity}, \text{qed}, -\text{logp}, -\text{sascore})$. They occupy non-dominated trade-off positions because their predicted affinity ($7.334 - 7.458$) and QED ($0.786$) are substantially higher than other low-LogP proposals.
2. **Epistemic Constraint:** These candidates are strictly identified as **near-threshold candidates that fail the constraint**. They are **not** claimed to satisfy all 4 objectives or to be feasible drug candidates.
3. **Dashboard Update:** The Streamlit dashboard ([`ui/app.py`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/ui/app.py#L260-L280)) now includes an explicit `"LogP <= 3.5"` status column marking these candidates as `"FAIL (3.531 > 3.5)"`, accompanied by a prominent warning caption.

### B. Objective Function Discrepancy
- **Earlier Design Specification (Phase 1.3.2)**: Proposed center-point deviation ($|\text{LogP} - 2.5|$).
- **Implemented Benchmark Code**: Defined direct monotonic minimization ($\min \text{LogP}$ with upper threshold $3.5$, normalized via $\tilde{y} = (6.0 - \text{LogP})/6.0$).
- **Implication**: Direct minimization drove the agent to introduce polar heteroatoms (lowering LogP from $4.136$ down to $3.531$, $\Delta = -0.605\log$). Because runs were bounded to 5 iterations on quinazoline substrates, the molecules remained in the drug-like range ($3.3 \le \text{LogP} \le 4.1$) without drifting into extreme hydrophilicity. All benchmark records are preserved as executed, with direct minimization documented as the active benchmark objective.

---

## 3. Disentangled QSAR Cross-Validation Metrics

Previous report drafts inadvertently conflated metrics across distinct validation protocols (e.g. citing the random-split $R^2 = 0.7050$ alongside the scaffold-split Spearman $\rho = 0.7640$).

The validation metrics are derived from 5,226 core biochemical compounds (one sample per standardized compound, median experimental $p\text{IC}_{50}$, equal weighting) evaluated on the selected primary surrogate model (**Model M2: Random Forest Hybrid** with 2,048 Morgan bits + 10 RDKit physicochemical descriptors).

### Disentangled Cross-Validation Table:
| Validation Protocol | Evaluation Purpose | $R^2$ Score | Spearman $\rho$ | MAE ($\log$ units) | RMSE ($\log$ units) |
|---|---|:---:|:---:|:---:|:---:|
| **Random 5-Fold CV** | In-distribution interpolation across chemical space | **0.7050** | **0.8397** | 0.5791 | 0.7737 |
| **Bemis–Murcko Scaffold 5-Fold CV** | Generalization to unseen chemical core scaffolds | **0.5849** | **0.7640** | 0.7074 | 0.9177 |
| **Document (Publication) 5-Fold CV** | Generalization across independent assay laboratories | **0.3199** | **0.5816** | 0.9276 | 1.1747 |

*Methodological Rule:* Random-split metrics ($R^2 = 0.7050, \, \rho = 0.8397$) and scaffold-split metrics ($R^2 = 0.5849, \, \rho = 0.7640$) characterize distinct generalization regimes and must never be mixed.

### Comparative Baseline Summary:
- **M0 (Dummy Mean Baseline)**: Random $R^2 = -0.0009$, $\rho = -0.0359$; Scaffold $R^2 = -0.0078$, $\rho = -0.1048$.
- **M1 (Random Forest Morgan Only)**: Random $R^2 = 0.7042$, $\rho = 0.8377$; Scaffold $R^2 = 0.5878$, $\rho = 0.7665$.
- **M2 (Random Forest Hybrid)**: Random $R^2 = 0.7050$, $\rho = 0.8397$; Scaffold $R^2 = 0.5849$, $\rho = 0.7640$.
- **M3 (HistGradientBoosting Descriptors Only)**: Random $R^2 = 0.3995$, $\rho = 0.6214$; Scaffold $R^2 = 0.2587$, $\rho = 0.5005$.

---

## 4. Narrowed Claim on Execution-Trace Self-Repair (ET-SRE)

Previous claims of "100% functional validity restoration" are formally narrowed to the specific empirical scope tested:

### A. Specific Fault Tested and Verified
ET-SRE was evaluated against **syntactic namespace faults injected during procedural tool synthesis**:
- **Injected Fault**: Omission of `from rdkit import Chem` during tool generation, producing `NameError: name 'Chem' is not defined` when `Chem.MolFromSmiles()` or `Chem.SanitizeMol()` was called.
- **Verification Trigger**: During skill registration, the VPSL sandbox executed the skill's unit test harness (`run_tests()`). The sandbox intercepted the runtime exception, and `VerificationFailureError` captured the full execution traceback.
- **Diagnosis and Patch**: ET-SRE inspected the traceback, diagnosed the missing AST symbol `Chem`, prepended `from rdkit import Chem\n`, re-executed the unit tests in the sandbox, verified test passing, and registered the repaired skill.
- **Subsequent Reuse**: In iterations 2–5, the planner retrieved the repaired skill from VPSL (`Skill Reuse = 4` on `CHEMBL289959`), successfully producing 50/50 valid candidate proposals.

### B. Explicit Non-Claims (Scope Restrictions)
- **No Claim of Arbitrary Fault Recovery**: ET-SRE has **not** been tested on, and does not claim to repair, complex semantic chemical errors (e.g. invalid SMARTS reaction transforms, stereochemical inversion bugs, or incorrect retrosynthetic step logic).
- **Narrowed Claim**: The claim is strictly restricted to **successful recovery of tested missing-import faults, sandbox re-verification, and subsequent procedural skill reuse in VPSL**.

---

## 5. Benchmark Scaffold Representation & Epistemic Boundaries

### A. Scaffold Familiarity vs Unseen Generalization
- **Strict Benchmark Holdout**: Neither `CHEMBL289959` nor `CHEMBL290096` was present in the core biochemical training dataset (`egfr_core_biochemical_compounds.csv`, 5,226 compounds).
- **Shared Scaffold Presence**: Exactly **276 distinct training compounds** share the 4-anilinoquinazoline Bemis–Murcko scaffold ($\texttt{c1ccc(Nc2ncnc3ccccc23)cc1}$).
- **Scientific Characterization**:
  - The autonomous rediscovery of `CHEMBL290096` ($R = 3-\text{Br}$) starting from `CHEMBL289959` ($R = \text{H}$) demonstrates that the surrogate model has learned substituent SAR around this specific hinge-binding scaffold.
  - It constitutes **in-domain / familiar-scaffold analogue recovery**.
  - It must **not** be presented as evidence of generalization to novel, unseen chemical scaffolds.

### B. Surrogate Model Epistemic Boundary
- Model predictions ($p\text{IC}_{50}$) are computational surrogate estimates from an empirical Random Forest model.
- Predicted affinity serves as an in silico ranking heuristic; it does **not** constitute experimental biochemical activity, binding confirmation, or clinical efficacy.

---

## 6. Dashboard & Runtime Verification Audit Summary

### A. Programmatic & Service Verification
- `ui/app.py` runs under Python 3.11 with Streamlit 1.28.0 compatibility (standard `st.container()`).
- Streamlit bound to port 8508, served `HTTP 200 OK` (892 bytes), and returned `ok` on `_stcore/health`.
- FastAPI backend (`api/main.py`) bound to port 8008, served `HTTP 200 OK` on `/health` and `/`.
- Subprocess execution bridge (`scripts/run_optimization_task.py`) was verified live, executing 1 closed-loop iteration and returning 5 non-dominated Pareto candidates.

### B. Stated Limitation on Visual Browser Checks
- Under Windows sandbox isolation, sandboxed processes cannot bind to external host-visible network ports without elevation (`BypassSandbox: true`).
- Consequently, **no visual browser screenshot check was performed**. We report only verified programmatic HTTP responses.

---

## 7. What PHARMAFORGE Actually Demonstrates & Remaining Limitations

### What PHARMAFORGE Demonstrates:
1. **Closed-Loop Multi-Agent Integration**: Seamless coordination of LLM-guided scientific planning, dynamic procedural tool synthesis in RDKit, deterministic sandboxed execution, unit-test verification (STO-IV), and surrogate-guided Pareto evaluation.
2. **Mitigation of Memory Poisoning**:
   - Without verification (Config B), buggy procedural tools poison persistent memory, resulting in repeated runtime crashes across subsequent iterations.
   - Deterministic verification (Config C) prevents memory poisoning by rejecting unverified code.
   - Traceback-driven self-repair (Config D) successfully repairs missing-import faults, re-verifies them, and enables persistent tool reuse across optimization iterations.
3. **Autonomous SAR Exploitation on Known Scaffolds**: Successfully navigates multi-objective trade-offs on the 4-anilinoquinazoline series, rediscovering known potent leads (`CHEMBL290096`) and proposing polar bioisosteric modifications that lower lipophilicity.

### Remaining Scientific Limitations:
1. **Narrow Repair Scope**: ET-SRE was evaluated against missing imports (`NameError`). Automated repair of arbitrary chemical logic bugs remains unverified.
2. **Limited Chemotype Diversity**: Evaluated exclusively on the 4-anilinoquinazoline kinase series with 276 training scaffold analogs. Cross-family generalization remains unproven.
3. **LogP Feasibility Boundary**: Pyrido-pyrimidine candidates reduce lipophilicity from $4.136$ to $3.531$, but strictly remain near-threshold non-compliant solutions (+0.031 above $\text{LogP} \le 3.500$).
4. **Surrogate Approximation**: Surrogate predictions ($p\text{IC}_{50}$) are statistical approximations and require downstream wet-lab experimental confirmation.
5. **Monte Carlo Volume Uncertainty**: Multi-objective hypervolumes are numerical Monte Carlo estimates ($N = 20,000$, $\text{SE} \approx 0.0026$), not exact continuous analytical integrals.

---

## 8. Artifact Inventory & Exact Locations

| Artifact | File Path | Status |
|---|---|:---:|
| **Final Phase 1.5 Report** | [`data/phase_1_5_final_verification_report.md`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/data/phase_1_5_final_verification_report.md) | Audited & Corrected |
| **Root Copy of Phase 1.5 Report** | [`phase_1_5_final_verification_report.md`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/phase_1_5_final_verification_report.md) | Audited & Corrected |
| **Phase 1.4 Benchmark Report** | [`data/phase_1_4_benchmark_optimization_report.md`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/data/phase_1_4_benchmark_optimization_report.md) | Reconciled |
| **Benchmark Optimization Runs (JSON)** | [`artifacts/benchmark_runs/benchmark_optimization_results.json`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/artifacts/benchmark_runs/benchmark_optimization_results.json) | Preserved & Recalculated |
| **Trained QSAR Surrogate Model** | [`artifacts/qsar/egfr_qsar_best_model.joblib`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/artifacts/qsar/egfr_qsar_best_model.joblib) | Preserved |
| **Streamlit Dashboard** | [`ui/app.py`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/ui/app.py) | Verified (Constraint Feasibility Column Added) |
| **FastAPI Backend** | [`api/main.py`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/api/main.py) | Verified (Live HTTP 200) |
| **Optimization CLI Runner** | [`scripts/run_optimization_task.py`](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/scripts/run_optimization_task.py) | Verified |
| **Walkthrough Artifact** | [`walkthrough.md`](file:///C:/Users/Sravya%20mummana/.gemini/antigravity/brain/d8ade012-07b4-4f00-984d-3f56c9323abc/walkthrough.md) | Updated |
