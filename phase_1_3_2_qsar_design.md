# PHARMAFORGE — Phase 1.3.2: QSAR Architecture & Experimental Design

**Document ID:** `data/phase_1_3_2_qsar_design.md`  
**Date:** 2026-11-20 (Reconciled and Finalized in Phase 1.3.2.2)  
**Target:** Human Epidermal Growth Factor Receptor (EGFR; `CHEMBL203`)  
**Phase Status:** Architecture and Experimental Design Only — Zero Model Training / Zero Package Installation  

---

## 1. Executive Summary & Objective

The objective of Phase 1.3.2 is to specify a scientifically rigorous, leakage-resistant Quantitative Structure-Activity Relationship (QSAR) modeling architecture for human EGFR binding affinity.

This predictive engine will replace earlier empirical proxy heuristics in PHARMAFORGE, providing the generative multi-agent framework with a calibrated surrogate estimator for target affinity ($p\text{IC}_{50}$) alongside established chemoinformatics properties (QED, LogP, SAScore).

**Core Design Tenet:** In accordance with computational medicinal chemistry standards, the QSAR model is formulated as an **imperfect statistical estimator** bounded by an explicit structural applicability domain and predictive dispersion, rather than an infallible experimental test-oracle.

---

## 2. Environment & Dependency Audit

Targeted inspection of the existing environment and workspace confirmed:

### 2.1 Python Environment & Package Versions
- **Python Version:** `3.13.16` (64-bit AMD64)
- **RDKit:** `2026.09.2` (Verified; includes `rdFingerprintGenerator`, `Chem.Descriptors`, `Chem.QED`, `SaltRemover`, `FilterCatalog`, and `MurckoScaffold`)
- **scikit-learn:** `1.9.1` (Verified; includes `DummyRegressor`, `RandomForestRegressor`, `HistGradientBoostingRegressor`, `StandardScaler`, and validation splitters)
- **NumPy:** `2.5.4` (Verified)
- **SciPy:** `1.18.1` (Verified)
- **PyTorch:** `2.15.0+cpu` (Available; intentionally deferred in favor of tabular tree baselines)
- **Standard Library:** Python's native `csv` module provides robust, zero-dependency streaming of datasets without requiring `pandas`.

### 2.2 Reusable Existing Components in PHARMAFORGE
To prevent code duplication, the QSAR architecture directly leverages:
1. **`ScientificTestOracle.generate_morgan_fingerprint()`** ([src/pharmaforge/verifier/oracle.py](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/src/pharmaforge/verifier/oracle.py#L172-L195)): Implements canonical Morgan circular fingerprints (ECFP4 equivalent) via `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)`.
2. **`ScientificTestOracle.verify_smiles()`** ([src/pharmaforge/verifier/oracle.py](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/src/pharmaforge/verifier/oracle.py#L53-L156)): Extracts canonical sanitized SMILES and exact physicochemical descriptors (`MolWt`, `MolLogP`, `TPSA`, `NumHDonors`, `NumHAcceptors`, `NumRotatableBonds`, `QED`, `SAScore`).
3. **`MoleculeCandidate` Schema** ([src/pharmaforge/core/types.py](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/src/pharmaforge/core/types.py#L33-L56)): Contains the native field `affinity_score: float | None`, which will store predicted $p\text{IC}_{50}$ alongside confidence and applicability domain metadata.

---

## 3. Prediction Task Formulation & Epistemic Boundaries

### 3.1 Task Definition
The modeling task is framed as a supervised non-linear regression problem:
$$\widehat{y} = f_{\boldsymbol{\theta}}(\mathbf{x})$$
where:
- $\mathbf{x} \in \mathbb{R}^D$ is the molecular representation derived from standardized SMILES.
- $\widehat{y} \in \mathbb{R}$ is the predicted EGFR in vitro biochemical potency in $p\text{IC}_{50}$ units:
  $$p\text{IC}_{50} = -\log_{10}(\text{IC}_{50}\text{ [M]}) = 9.0 - \log_{10}(\text{IC}_{50}\text{ [nM]})$$

### 3.2 Epistemic Distinction of Scientific Quantities
To prevent hallucination and scientific over-claiming within the agent system, PHARMAFORGE strictly enforces four distinct ontological categories:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                       EPISTEMIC ONTOLOGY OF POTENCY                         │
├───────────────────────────────────┬─────────────────────────────────────────┤
│ Concept                           │ Definition & Operational Boundary       │
├───────────────────────────────────┼─────────────────────────────────────────┤
│ 1. Experimental IC50              │ Direct laboratory observation (nM) of   │
│                                   │ 50% enzyme inhibition under specified   │
│                                   │ in vitro assay conditions.              │
├───────────────────────────────────┼─────────────────────────────────────────┤
│ 2. Derived pIC50                  │ Exact mathematical transformation       │
│                                   │ (9.0 - log10(IC50)) of experimental     │
│                                   │ measurement; ground-truth training label│
├───────────────────────────────────┼─────────────────────────────────────────┤
│ 3. Model-Predicted pIC50          │ Statistical estimate produced by QSAR   │
│                                   │ model; carries error bounds and         │
│                                   │ applicability domain limits.            │
├───────────────────────────────────┼─────────────────────────────────────────┤
│ 4. Clinical Drug Efficacy         │ In vivo therapeutic efficacy depending   │
│                                   │ on cellular uptake, ADMET, clearance,   │
│                                   │ and off-target safety. NOT modeled here.│
└───────────────────────────────────┴─────────────────────────────────────────┘
```
**Strict System Boundary:** The multi-agent orchestrator, planner, and verifier must **never** report model-predicted $p\text{IC}_{50}$ as experimentally validated activity or clinical efficacy.

### 3.3 Training Dataset Selection & Stratification Strategy
From the Phase 1.3.1.3 audit, 8,944 unified records were stratified into 5 assay contexts. Among the 7,053 high-confidence records, 7,044 are annotated with `eligible_for_qsar_training == 'True'` (5,970 unique compounds).

#### Stratified Context Breakdown:
1. **Core Homogeneous Cell-Free Biochemical Subset (Primary Training Corpus):**
   - Composed of the union of:
     - `purified_recombinant_egfr.csv`: 1,791 eligible records (1,638 unique compounds)
     - `other_cell_free_biochemical.csv`: 4,202 eligible records (3,752 unique compounds)
   - **Total Core Records:** **5,993 records** (85.08% of eligible high-confidence records).
   - **Total Core Unique Compounds:** **5,226 compounds** (87.54% of eligible high-confidence compounds).
   - **Biochemical Justification:** Both strata measure direct in vitro enzymatic kinase inhibition on synthetic or peptide substrates in homogeneous cell-free solutions. They operate without intact cells, membrane partitioning barriers, or organelle compartmentalization. These represent **confirmed cell-free biochemical formats**.
2. **`cell_free_membrane_derived.csv` (169 eligible records, 168 unique compounds):**
   - Assays using cell-membrane vesicle homogenates (e.g. A431 membranes), including benchmark assay `CHEMBL680017`.
   - **Protocol:** Held out as an **External Protocol-Transfer Test Set** to evaluate how models trained on soluble enzyme generalize to membrane-bound receptor inhibition, accounting for lipophilicity-driven lipid bilayer partitioning.
3. **`single_protein_terse_metadata.csv` (882 eligible records, 743 unique compounds):**
   - Single-protein format (`BAO_0000357`) with brief textual descriptions.
   - **Protocol:** Reserved for secondary sensitivity analysis and domain-expansion evaluation.
4. **Exclusion of Severe Conflicts (`egfr_conflict_records.csv`, 1,870 records; 421 compounds):**
   - Records with $> 1.0\log$ divergence across publications are excluded to prevent conflicting gradient updates.
5. **Exclusion of Non-Standard Formats (`unresolved_non_standard.csv`, 22 records):**
   - Excluded due to unstandardized assay formats.

---

## 4. Repeated-Measurement Policy

### 4.1 Empirical Distribution
Inspection of the 7,044 eligible records in `data/egfr_high_confidence.csv` reveals:
- **5,128 compounds (85.9%)** possess a single observation (5,128 records; 72.8%).
- **842 compounds (14.1%)** possess repeated observations (1,916 records; 27.2%), ranging from 2 up to 9 measurements per compound.

### 4.2 Reconciled Strategy: Compound-Level Median Target ($w_i = 1.0$)
In Phase 1.3.2.1, a concordance weight formula $w_i = 1 + \ln(M_i) \cdot \exp(-\Delta p\text{IC}_{50, i})$ was proposed with the claim that weights reached a maximum of approximately 1.5. 
- **Numerical Audit:** Re-calculation shows weights reach **2.7563** (for `CHEMBL1230609`, $M_i = 6$, spread $= 0.020\log$). 674 compounds exceed $w_i = 1.5$.
- **Scientific Rejection of Concordance Weighting:** Publication frequency ($M_i$) in ChEMBL is an artifact of chemical popularity, patent history, and commercial availability—not intrinsic molecular energetics. Up-weighting popular series by up to $2.75\times$ severely distorts loss gradients, biasing the model against novel, single-observation chemotypes.

**Adopted Primary Specification (Strategy A):**
1. **One Compound, One Sample:** Each unique standardized molecule represents exactly **one training sample** ($N = 5,226$ in the core homogeneous corpus; $N = 5,970$ across all high-confidence strata).
2. **Median Potency Target:** The ground-truth target is the median of concordant measurements:
   $$\bar{y}_i = \text{median}\left(\{y_{ij}\}_{j=1}^{M_i}\right)$$
   Taking the median eliminates single-assay noise and robustly estimates central potency within the $\Delta \le 1.0\log$ concordant envelope.
3. **Equal Compound Weights:** Every compound receives identical weight $w_i = 1.0000$, ensuring loss gradients reflect structural diversity rather than literature frequency.
4. **Provenance Preservation:** All 7,044 individual activity records, assay IDs, document IDs, and raw measurements remain 100% intact in source CSV files.

---

## 5. Molecular Feature Representations

A robust, low-complexity feature engineering strategy is specified:

### 5.1 Feature Set 1: Circular Morgan Fingerprints (ECFP4 Equivalent)
- **Algorithm:** Morgan circular fingerprint via RDKit `rdFingerprintGenerator.GetMorganGenerator`.
- **Parameters:** Radius $r = 2$ (diameter 4), bit length $N = 2,048$ bits.
- **Biochemical Rationale:** Captures circular atom environments up to 4 bonds away, encoding specific pharmacophoric features (hydrogen-bond donor/acceptor patterns, aromatic rings, halogen substitutions, kinase hinge-binding motifs).
- **Representation:** Sparse binary vector $\mathbf{x}_{\text{fp}} \in \{0, 1\}^{2048}$.

### 5.2 Feature Set 2: Compact Physicochemical Descriptors
To complement local circular fingerprints with global molecular properties, a compact 10-dimensional vector of interpretable RDKit descriptors is computed:
1. `MolWt`: Molecular Weight (size/bulk)
2. `MolLogP`: Wildman-Crippen octanol-water partition coefficient (lipophilicity)
3. `TPSA`: Topological Polar Surface Area (polar surface available for hydrogen bonding)
4. `NumHDonors`: Number of hydrogen bond donors (critical for kinase hinge interactions at Met793)
5. `NumHAcceptors`: Number of hydrogen bond acceptors (hinge and gatekeeper interactions)
6. `NumRotatableBonds`: Molecular flexibility and conformational entropy
7. `NumAromaticRings`: $\pi$-$\pi$ stacking propensity in the kinase hydrophobic pocket
8. `HeavyAtomCount`: Total non-hydrogen atom count
9. `FractionCSP3`: Fraction of $sp^3$-hybridized carbons (three-dimensionality)
10. `QED`: Quantitative Estimate of Drug-Likeness (Bickerton et al.)

### 5.3 Feature Set 3: Hybrid Concatenation
$$\mathbf{x}_{\text{hybrid}} = [\mathbf{x}_{\text{fp}} \,\|\, \mathbf{z}_{\text{desc}}] \in \mathbb{R}^{2058}$$
where $\mathbf{z}_{\text{desc}} = \text{StandardScaler}(\mathbf{x}_{\text{desc}})$ is z-score normalized strictly using the training partition's mean and variance.

---

## 6. Finalized Baseline Model-Feature Comparison Matrix

The 4-way baseline matrix is finalized as follows:

| Model ID | Model Family | Feature Representation | Input Dimensions | Primary Scientific Hypothesis Tested |
| :---: | :--- | :--- | :---: | :--- |
| **M0** | `DummyRegressor(strategy='mean')` | None (null model) | 0 | **Baseline Mean Floor:** Predicts training set mean $\bar{y}_{\text{train}}$. On test sets, $R^2_{\text{test}} \le 0$ depending on train vs test mean deviation. Any valid QSAR model must achieve substantially lower MAE/RMSE and positive $R^2$. |
| **M1** | `RandomForestRegressor` | Morgan Fingerprints (ECFP4; $r=2$) | 2,048 bits | **H1:** Local 2D circular substructures and pharmacophoric bits alone provide sufficient signal for target affinity. |
| **M2** | `RandomForestRegressor` | Hybrid (Morgan 2048 + 10 Descriptors) | 2,058 dims | **H2:** Augmenting local circular bits with global physical properties (MW, LogP, TPSA, HBD, HBA) improves generalization. |
| **M3** | `HistGradientBoostingRegressor` | Compact Descriptors Only | 10 continuous | **H3:** Low-dimensional interpretable physical properties with non-linear boosting thresholds can approximate affinity rapidly. |

---

## 7. Experimental Split Protocols & Leakage Prevention

PHARMAFORGE defines **four distinct evaluation protocols**:

### 7.1 Protocol A: Molecule-Grouped Random 5-Fold Cross-Validation
- **Structure:** Grouped 5-Fold CV using `standardized_smiles` as grouping key.
- **Rule:** All activity records for the exact same standardized molecule are strictly confined to the same fold.
- **Role:** Diagnostic baseline measuring pure in-distribution interpolation.

### 7.2 Protocol B: Bemis–Murcko Scaffold 5-Fold Cross-Validation
- **Structure:** `GroupKFold` using the Bemis–Murcko scaffold string as grouping key.
- **Rule:** All molecules sharing the same core ring-and-linker scaffold are assigned together to either the training set or test set.
- **Role:** Primary evaluation metric. Tests whether the model learns transferable SAR rules that generalize to novel chemical scaffolds.

### 7.3 Protocol C: Leakage-Safe Document Evaluation (Connected-Component Split)
- **Structure:** Bipartite graph partitioning of the compound-document network (538 connected components; maximum component 15.4% / 920 compounds / 58 documents).
- **Rule:** Greedy packing of connected components into 5 folds guarantees **ZERO shared compounds AND ZERO shared documents** between training and test sets simultaneously.
- **Alternative:** Disjoint single-document series holdout (evaluating on complete lead-optimization series published by an unseen laboratory).

### 7.4 Protocol D: Exact Benchmark Holdout Evaluation
- **Benchmark Leads:** `CHEMBL289959` ($344.0\text{ nM}$, $p\text{IC}_{50} = 6.46$) and `CHEMBL290096` ($27.0\text{ nM}$, $p\text{IC}_{50} = 7.57$).
- **Rules:**
  1. Both benchmark starting compounds (`is_benchmark_compound == 'True'`, 9 records) are strictly held out from the training partition.
  2. For the scaffold-generalization experiment, all **354 unique compounds (466 eligible records)** sharing the 4-anilinoquinazoline scaffold (`shares_benchmark_scaffold == 'True'`) are completely purged from training ($6,578$ records remain).
- **Scientific Limitation Caveat:** Evaluating on only two benchmark starting points evaluates local matched-molecular-pair sensitivity ($R = \text{H}$ vs $R = \text{Br}$ near gatekeeper Thr790), but cannot prove broad multi-target or multi-chemotype generative capability.

---

## 8. Evaluation Metrics, Uncertainty & Applicability Domain Strategy

### 8.1 Quantitative Performance Metrics
Every evaluation split will report:
1. **MAE (Mean Absolute Error):** Average absolute prediction error in $p\text{IC}_{50}$ units.
2. **RMSE (Root Mean Squared Error):** Penalizes catastrophic mispredictions.
3. **$R^2$ (Coefficient of Determination):** Variance explained relative to mean.
4. **Spearman Rank Correlation ($\rho$):** Non-parametric rank ordering fidelity.

### 8.2 Uncertainty Evaluation Protocol
- **Provisional Status:** Random Forest tree dispersion ($\sigma_{\text{trees}}$) is formally categorized as an **uncalibrated heuristic measure of ensemble disagreement**, reflecting decision path divergence. It is not Bayesian epistemic uncertainty.
- **Error-versus-Retained-Coverage Curves:**
  - Evaluated on held-out test predictions across broad retained-coverage tiers: **100%, 90%, 80%, 70%, 50%**.
  - Assessed by the Spearman rank correlation $\rho(|y - \widehat{y}|, \sigma_{\text{trees}})$ and overall downward error trend across major deciles. Strict step-by-step monotonicity is not asserted due to finite fold sample noise.
- **Prediction Intervals & Calibration:**
  - Nominal intervals ($\pm 1.645\sigma$, $\pm 1.96\sigma$) will not be reported as calibrated 90%/95% confidence intervals without empirical coverage testing.
  - Where coverage guarantees are required, **Split Conformal Prediction** will be applied on an independent calibration set to guarantee distribution-free coverage $\mathbb{P}(Y \in \widehat{C}(X)) \ge 1 - \alpha$.

### 8.3 Dual-Axis Telemetry: Applicability Domain vs Predictive Dispersion
Structural similarity and predictive dispersion are reported as **independent, orthogonal diagnostics**:
1. **Structural Applicability Domain ($T_{\max}$):** Maximum Morgan fingerprint Tanimoto similarity to the nearest training set neighbor:
   - High: $T_{\max} \ge 0.60$
   - Moderate: $0.45 \le T_{\max} < 0.60$
   - Low: $0.35 \le T_{\max} < 0.45$
   - Novel / Out-of-Domain: $T_{\max} < 0.35$
2. **Predictive Dispersion ($\sigma_{\text{trees}}$):** Ensemble standard deviation across decision trees.
3. **Agent Rule:** Never combine $T_{\max}$ and $\sigma_{\text{trees}}$ into a single uncalibrated score. Report both distinctly to enable the planner to differentiate activity cliffs (high $T_{\max}$, high $\sigma$) from novel extrapolated space (low $T_{\max}$, low $\sigma$).

---

## 9. Multi-Agent System Integration

The QSAR predictor integrates with PHARMAFORGE's existing components:
1. **Molecular Optimization Planner (`orchestrator.py`):** Uses $\widehat{p\text{IC}}_{50}$ to guide multi-objective Pareto optimization alongside QED, LogP, and SAScore.
2. **Scientific Test-Oracle (`oracle.py`):** Validates chemical validity (valence, sanitization, PAINS, Brenk) first. Chemically invalid molecules are discarded immediately before scoring.
3. **Self-Repair & Reflexion Engine:** If the predictor indicates $T_{\max} < 0.35$ or high dispersion $\sigma > 1.0$, the self-repair loop prompts the planner to propose smaller bioisosteric modifications rather than drastic scaffold hops.
4. **Procedural Skill Library (`skill_library.py`):** Featurization, inference, and applicability domain checks are encapsulated into a reusable procedural skill.

---

## 10. Final Verdict

### **DESIGN APPROVED FOR IMPLEMENTATION**

All architectural components, feature designs, baseline model choices, split protocols, repeated-measurement policies, assay-aggregation scopes, and uncertainty standards are fully reconciled and methodologically verified.

**Implementation is strictly stopped. Awaiting explicit user approval before proceeding to Phase 1.3.3.**
