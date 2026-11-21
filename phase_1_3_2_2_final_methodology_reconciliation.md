# PHARMAFORGE — Phase 1.3.2.2: Final QSAR Methodology Reconciliation

**Document ID:** `data/phase_1_3_2_2_final_methodology_reconciliation.md`  
**Date:** 2026-11-20  
**Target:** Human Epidermal Growth Factor Receptor (EGFR; `CHEMBL203`)  
**Phase Status:** Methodological Reconciliation & Design Finalization — Zero Model Training / Zero Package Installation  

---

## 1. Executive Summary

Phase 1.3.2.2 resolves the remaining methodological and numerical questions identified in the Phase 1.3.2 and 1.3.2.1 QSAR experimental designs before code implementation.

This audit establishes four foundational methodological decisions:
1. **Repeated-Measurement Weighting:** The previous claim that $w_i = 1 + \ln(M_i) \cdot \exp(-\Delta p\text{IC}_{50, i})$ yields a maximum weight of $\approx 1.5$ is mathematically refuted; actual weights reach **2.7563**. More importantly, concordance weighting is discarded on scientific grounds: publication frequency is an artifact of chemical popularity rather than molecular energetics. The primary model adopts **Strategy A: compound-level median target with equal compound weights ($w_i = 1.0$)**, preserving all 7,044 raw observation records and provenance.
2. **Assay-Compatible Aggregation:** The primary QSAR model is restricted to the **Core Homogeneous Cell-Free Biochemical Subset** (union of `purified_recombinant_egfr` and `other_cell_free_biochemical`), encompassing **5,993 records (85.08%) and 5,226 unique standardized compounds (87.54%)**. The `cell_free_membrane_derived` stratum (169 eligible records, 168 compounds) is held out as an external protocol-transfer evaluation partition to prevent lipid partitioning artifacts.
3. **Uncertainty & Applicability Domain:** Random Forest tree dispersion ($\sigma_{\text{trees}}$) is formally categorized as an uncalibrated heuristic measure of ensemble disagreement, not Bayesian epistemic uncertainty. Nominal intervals ($\pm 1.645\sigma$, $\pm 1.96\sigma$) will not be claimed as calibrated confidence intervals without empirical coverage validation; split conformal prediction is specified if valid coverage guarantees are required. Error-versus-coverage curves are evaluated via rank correlation and broad percentile deciles rather than requiring strict step-by-step monotonicity. Structural domain similarity ($T_{\max}$) and predictive dispersion ($\sigma_{\text{trees}}$) are decoupled into orthogonal diagnostics.
4. **Baseline & Quality Claims:** `DummyRegressor(strategy='mean')` is accurately formulated as a mean-prediction baseline with test $R^2 \le 0$ rather than exactly zero. The 85.1% assay proportion is cautiously phrased as "confirmed cell-free biochemical formats" rather than asserting universal verification of catalytic domain construct boundaries.

---

## 2. Repeated-Measurement Policy Reconciliation

### 2.1 Mathematical Audit of the Concordance Weight Formula

The formula previously proposed in Phase 1.3.2.1 was:
$$w_i = 1.0 + \ln(M_i) \cdot \exp\left(-\Delta p\text{IC}_{50, i}\right)$$
where $M_i$ is the number of eligible measurements for standardized compound $i$, and $\Delta p\text{IC}_{50, i} = \max_j(y_{ij}) - \min_j(y_{ij})$ is the measurement spread.

An empirical re-calculation across all 5,970 unique compounds in the 7,044 eligible records reveals:
- **Minimum weight:** $1.0000$ (for all 5,128 compounds with $M_i = 1$, where $\ln(1) = 0$).
- **Maximum weight:** **$2.7563$** (observed for `CHEMBL1230609`, with $M_i = 6$, spread $= 0.020\log$, where $1 + \ln(6) \cdot \exp(-0.020) = 1 + 1.7918 \cdot 0.9802 = 2.7563$).
- **Distribution of weights above 1.0:**
  - $w_i > 1.5$: **674 compounds** (11.29% of all compounds; 80.05% of repeated-measurement compounds).
  - $w_i > 2.0$: **48 compounds** (0.80% of all compounds; 5.70% of repeated-measurement compounds).
  - Maximum observation count: $M_i = 9$ (observed for `CHEMBL296407`, spread $= 0.700\log$, weight $= 2.0911$).

**Error Rectification:** The prior report's statement that the formula produces a maximum weight of "approximately 1.5" was mathematically incorrect. The weight reaches nearly $2.76$, granting certain replicated compounds nearly triple the gradient influence of singleton compounds.

### 2.2 Scientific Critique of Concordance-Based Weighting

Weighting training samples proportionally to concordance ($\ln(M_i) \cdot \exp(-\Delta)$) is scientifically invalid for QSAR modeling:
1. **Publication Bias vs Biophysical Reality:** The number of times a compound appears in the scientific literature ($M_i$) reflects historical popularity, commercial reagent availability, or extensive patent follow-up (e.g. gefitinib or erlotinib derivatives). It is not an intrinsic physical or thermodynamic property of the ligand-receptor binding interface.
2. **Gradient Distortion:** Up-weighting well-characterized series by up to $2.75\times$ causes supervised algorithms (decision tree splits, gradient boosts, or neural loss functions) to over-fit to saturated chemical clusters at the expense of novel, single-observation chemotypes. A novel compound synthesized and tested once by an academic group should not have its gradient suppressed to one-third the importance of a commercial kinase inhibitor.
3. **Redundancy with Quality Curation:** All compounds in `data/egfr_high_confidence.csv` have already passed the severe-conflict filter ($\Delta p\text{IC}_{50} \le 1.0\log$). Discordant data have already been excluded. Within this concordant envelope, single-assay noise is effectively attenuated by taking the median without requiring artificial sample weight scaling.

### 2.3 Systematic Strategy Comparison

| Evaluation Dimension | Strategy A: Compound-Level Median ($w_i = 1.0$) | Strategy B: Bounded Concordance Weighting | Strategy C: Observation-Level Inverse Weights ($w_{ij} = 1/M_i$) |
| :--- | :--- | :--- | :--- |
| **Sample Count ($N$)** | $5,970$ unique molecules | $5,970$ unique molecules | $7,044$ observation rows |
| **Target Formulation** | $\bar{y}_i = \text{median}(\{y_{ij}\}_{j=1}^{M_i})$ | $\bar{y}_i = \text{median}(\{y_{ij}\}_{j=1}^{M_i})$ | Individual $y_{ij}$ per row |
| **Weight Formulation** | $w_i = 1.0$ (equal weight) | $w_i = 1 + \alpha \frac{\ln(M_i)}{\ln(M_{\max})} e^{-\Delta} \le 1.25$ | $w_{ij} = 1 / M_i$ (compound sum $= 1.0$) |
| **Publication Bias Risk** | **Zero** (equal per compound) | Moderate (favors popular series) | Low (per-compound sum bounded) |
| **Tree-Model Compatibility** | **Optimal** (unique $\mathbf{x}_i \to \bar{y}_i$ mapping) | Good (sample weights supported) | **Sub-optimal** (identical $\mathbf{x}_i$ mapped to differing $y_{ij}$ creates artificial split noise) |
| **Robustness to Assay Noise**| High (median filters outliers) | High (median filters outliers) | Low (each individual assay is a target) |
| **Provenance Preservation**| Full (original rows preserved) | Full (original rows preserved) | Full (original rows preserved) |

### 2.4 Adopted Primary Policy: Strategy A

**Decision:** PHARMAFORGE designates **Strategy A (Compound-level median target with equal compound weights $w_i = 1.0$)** as the primary QSAR training specification.
- Every unique standardized compound represents exactly **one training sample** ($N = 5,970$ in high-confidence; $N = 5,226$ in core homogeneous).
- The scalar target is the median derived potency:
  $$\bar{y}_i = \text{median}\left(\{y_{ij}\}_{j=1}^{M_i}\right)$$
- Every compound receives uniform weight $w_i = 1.0000$.
- All 7,044 individual activity records, assay IDs, and publication sources remain 100% preserved in `data/egfr_high_confidence.csv` and stratum CSV files.

---

## 3. Assay-Compatible Aggregation Policy

### 3.1 Stratified Distribution in Eligible QSAR Records

Inspection of the 7,044 eligible records across the context strata establishes:

| Assay Stratum | Eligible Records | Unique Compounds | Percentage of Records | Percentage of Compounds | Assay Format Characteristics |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **`purified_recombinant_egfr`** | 1,791 | 1,638 | 25.43% | 27.44% | Soluble recombinant/purified enzyme; synthetic/peptide substrate phosphorylation. |
| **`other_cell_free_biochemical`** | 4,202 | 3,752 | 59.65% | 62.85% | Cell-free enzymatic inhibition (HTRF, FRET, radiometric, Caliper mobility shift). |
| **`single_protein_terse_metadata`** | 882 | 743 | 12.52% | 12.45% | Single-protein format (`BAO_0000357`) with minimal textual assay descriptions. |
| **`cell_free_membrane_derived`** | 169 | 168 | 2.40% | 2.81% | Cell membrane vesicles (e.g. A431 membranes); includes benchmark assay `CHEMBL680017`. |
| **`unresolved_non_standard`** | 0 | 0 | 0.00% | 0.00% | Excluded non-standard formats (22 records total, 0 eligible). |
| **Total High-Confidence** | **7,044** | **5,970** | **100.0%** | **100.0%** | (Union across all 4 eligible strata) |

### 3.2 Cross-Strata Compound Overlap

Only **316 compounds (5.29%)** possess measurements that span multiple assay strata:
- `('other_cell_free_biochemical', 'purified_recombinant_egfr')`: 160 compounds
- `('other_cell_free_biochemical', 'single_protein_terse_metadata')`: 78 compounds
- `('cell_free_membrane_derived', 'other_cell_free_biochemical')`: 27 compounds
- `('cell_free_membrane_derived', 'single_protein_terse_metadata')`: 25 compounds
- `('cell_free_membrane_derived', 'other_cell_free_biochemical', 'single_protein_terse_metadata')`: 11 compounds
- `('purified_recombinant_egfr', 'single_protein_terse_metadata')`: 10 compounds
- `('cell_free_membrane_derived', 'other_cell_free_biochemical', 'purified_recombinant_egfr')`: 2 compounds
- `('other_cell_free_biochemical', 'purified_recombinant_egfr', 'single_protein_terse_metadata')`: 2 compounds
- `('cell_free_membrane_derived', 'purified_recombinant_egfr')`: 1 compound

### 3.3 Core Homogeneous Cell-Free Biochemical Subset

To ensure biophysical consistency and eliminate assay-format artifacts, the primary QSAR model will be trained on the **Core Homogeneous Cell-Free Biochemical Subset**:
$$\text{Core Dataset} = \text{purified\_recombinant\_egfr} \cup \text{other\_cell\_free\_biochemical}$$
- **Total Records:** **5,993** (85.08% of all eligible records).
- **Unique Standardized Compounds:** **5,226** (87.54% of all eligible compounds).
- **Biochemical Rationale:** Both strata evaluate direct in vitro catalytic phosphorylation of peptide substrates by soluble EGFR kinase without living cells or lipid bilayers. Because discordant measurements across laboratories have already been excluded ($\Delta \le 1.0\log$), pooling these two strata creates a large, statistically powerful, and biophysically homogeneous training corpus.

### 3.4 Separate Partitions & Protocol-Transfer Evaluation

The remaining strata are handled under explicit scientific protocols:
1. **`cell_free_membrane_derived` (169 records, 168 compounds):**
   - *Biophysical Distinction:* These assays utilize membrane vesicle preparations (e.g. A431 epidermoid carcinoma membrane fragments). Molecules must partition into or traverse lipid bilayers to reach the receptor kinase domain. Potency measurements in membrane vesicles can systematically diverge from purified enzyme assays due to compound lipophilicity and non-specific membrane partitioning.
   - *Protocol:* Held out as an **External Protocol-Transfer Test Set**. The primary QSAR model (trained on soluble biochemical assays) will predict these membrane-derived potencies to quantify protocol-transfer shift ($\Delta = y_{\text{membrane}} - \widehat{y}_{\text{purified}}$) as a function of $c\text{LogP}$.
2. **`single_protein_terse_metadata` (882 records, 743 compounds):**
   - *Biophysical Distinction:* Classified under `BAO_0000357` (single-protein format) but lacking detailed assay protocol text in ChEMBL annotations.
   - *Protocol:* Evaluated as an **External Domain Expansion Split** or pooled only in secondary sensitivity runs.
3. **Cross-Publication & Cross-Format Aggregation Protocol:**
   - For any compound with repeated measurements within the core homogeneous pool across different papers or assay IDs (164 shared compounds), the target label is the median $p\text{IC}_{50}$ of those measurements.
   - Observations from `cell_free_membrane_derived` are **never averaged** into the core biochemical target; they remain isolated in their stratum for transfer testing.

---

## 4. Uncertainty & Applicability-Domain Evaluation Reconciliation

### 4.1 Epistemic Boundaries of Random Forest Tree Dispersion

In the proposed Random Forest Regressors (M1 and M2), the standard deviation across individual tree predictions is defined as:
$$\sigma_{\text{trees}}(\mathbf{x}) = \sqrt{\frac{1}{B} \sum_{b=1}^B \left(f_b(\mathbf{x}) - \bar{f}(\mathbf{x})\right)^2}$$
- **Scientific Clarification:** $\sigma_{\text{trees}}$ is an **empirical heuristic measuring ensemble disagreement**, reflecting the variance of bagged estimators.
- **Explicit Caveats:**
  1. It is **not** calibrated Bayesian epistemic uncertainty.
  2. It does **not** capture intrinsic aleatoric assay measurement noise.
  3. It does **not** detect out-of-distribution extrapolation failure reliably: decision trees predict constant values outside the convex hull of training features, which can lead to low tree variance despite complete model failure.

### 4.2 Error-versus-Retained-Coverage Curve Evaluation

The previous design claimed that pruning uncertain predictions would yield a strictly monotonic decline in test error at every step. This claim is corrected:
- On finite held-out test splits ($N \approx 1,000$ to $1,200$ molecules per fold), removing small increments of compounds (e.g. 1% to 5%) is subject to stochastic sampling noise, causing minor non-monotonic fluctuations.
- **Reconciled Evaluation Standard:**
  1. Evaluate test MAE and RMSE across broad retained-coverage tiers: **100%, 90%, 80%, 70%, 50%**. A successful uncertainty heuristic must demonstrate a clear general downward trend across these major coverage deciles.
  2. Compute the **Spearman rank correlation** between the absolute prediction error $|y_i - \widehat{y}_i|$ and the uncertainty score $\sigma_{\text{trees}, i}$:
     $$\rho_{\text{error}} = \text{SpearmanRank}\left(|y_i - \widehat{y}_i|,\, \sigma_{\text{trees}, i}\right) > 0$$
  3. Strict step-by-step monotonicity is not asserted or required.

### 4.3 Prediction Intervals & Split Conformal Prediction

- **Rejection of Nominal Gaussian Intervals:** Intervals computed as $[\widehat{y} \pm 1.645\sigma_{\text{trees}}]$ or $[\widehat{y} \pm 1.96\sigma_{\text{trees}}]$ assume normality and unbiased dispersion. They **will not be claimed as calibrated 90% or 95% confidence intervals** without explicit empirical validation.
- **Specification of Split Conformal Prediction:** If statistically rigorous prediction intervals with guaranteed finite-sample marginal coverage are required:
  1. Reserve an independent calibration split $(X_{\text{cal}}, y_{\text{cal}})$ of $n_{\text{cal}}$ compounds unseen during model training.
  2. Compute normalized non-conformity scores:
     $$s_i = \frac{|y_i - \widehat{y}_i|}{\sigma_{\text{trees}, i}}$$
  3. Determine the empirical quantile:
     $$\widehat{q} = \text{Quantile}\left(\{s_i\}_{i=1}^{n_{\text{cal}}},\, \frac{\lceil(n_{\text{cal}}+1)(1-\alpha)\rceil}{n_{\text{cal}}}\right)$$
  4. Construct the calibrated prediction interval for any unseen molecule:
     $$\widehat{C}(x_{\text{test}}) = \left[\widehat{y}(x_{\text{test}}) - \widehat{q} \cdot \sigma_{\text{trees}}(x_{\text{test}}),\, \widehat{y}(x_{\text{test}}) + \widehat{q} \cdot \sigma_{\text{trees}}(x_{\text{test}})\right]$$
  5. Under exchangeability, this guarantees $\mathbb{P}\left(Y \in \widehat{C}(X)\right) \ge 1 - \alpha$ without distributional assumptions.

### 4.4 Decoupling Structural Applicability Domain from Predictive Dispersion

To provide clean, interpretable diagnostics to the multi-agent system, structural distance and predictive variance are reported as **independent, orthogonal signals**:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                   DUAL-AXIS DIAGNOSTIC TELEMETRY MATRIX                     │
├────────────────────────────────┬────────────────────────────────────────────┤
│ Diagnostic Dimension           │ Definition & Operational Boundary          │
├────────────────────────────────┼────────────────────────────────────────────┤
│ 1. Structural Applicability    │ Nearest-Neighbor Tanimoto Similarity:      │
│    Domain ($T_{\max}$)         │ $T_{\max}(\mathbf{x}) = \max_{\mathbf{x}_k \in \mathcal{D}_{\text{train}}} Tanimoto(\mathbf{x}, \mathbf{x}_k)$ │
│                                │ Binned into: High ($\ge 0.60$), Medium     │
│                                │ ($0.45\text{--}0.60$), Low ($0.35\text{--}0.45$), Novel ($< 0.35$).│
├────────────────────────────────┼────────────────────────────────────────────┤
│ 2. Predictive Dispersion       │ Standard deviation of ensemble trees:      │
│    ($\sigma_{\text{trees}}$)   │ Measures local ambiguity in tree splits;   │
│                                │ high near activity cliffs.                 │
└────────────────────────────────┴────────────────────────────────────────────┘
```

**Diagnostic Decoupling Rule:** The agent planner and self-repair engine must never collapse $T_{\max}$ and $\sigma_{\text{trees}}$ into a single uncalibrated score. 
- A candidate molecule with $T_{\max} = 0.82$ and $\sigma_{\text{trees}} = 1.15$ signals an **activity cliff** (close structural analog with conflicting SAR paths).
- A candidate molecule with $T_{\max} = 0.28$ and $\sigma_{\text{trees}} = 0.35$ signals **extrapolated flatland** (novel scaffold falling into an unbranched leaf where all trees predict the global mean).

---

## 5. Correction of Baseline and Quality Claims

### 5.1 Formulation of Baseline Model M0 (`DummyRegressor`)

- **Prior Statement:** "Baseline statistical floor ($R^2 = 0$). Any predictive QSAR model must substantially beat this."
- **Mathematical Correction:** `DummyRegressor(strategy='mean')` always outputs the training partition mean $\bar{y}_{\text{train}}$. On an independent test partition with mean $\bar{y}_{\text{test}}$ and variance $\sigma^2_{\text{test}}$, the test $R^2$ is:
  $$R^2_{\text{test}} = 1 - \frac{\frac{1}{N_{\text{test}}}\sum_{i=1}^{N_{\text{test}}} (y_i - \bar{y}_{\text{train}})^2}{\frac{1}{N_{\text{test}}}\sum_{i=1}^{N_{\text{test}}} (y_i - \bar{y}_{\text{test}})^2} = -\frac{(\bar{y}_{\text{test}} - \bar{y}_{\text{train}})^2}{\sigma^2_{\text{test}}} \le 0$$
- If the test fold mean matches the training fold mean exactly, $R^2_{\text{test}} = 0$. If the test fold mean deviates, $R^2_{\text{test}}$ is strictly negative. 
- **Adopted Description:** `DummyRegressor(strategy='mean')` is designated as a **standard mean-prediction baseline**. Its test MAE and RMSE represent the uninformative statistical floor; meaningful QSAR models must achieve substantially lower MAE/RMSE and positive $R^2$.

### 5.2 Cautious Phrasing of Cell-Free Assay Formats

- **Prior Statement:** "85.1% of assays are confirmed EGFR catalytic domain assays."
- **Evidence-Based Correction:** In the eligible high-confidence dataset, 85.08% (5,993 / 7,044) of records belong to the `purified_recombinant_egfr` (1,791 records) and `other_cell_free_biochemical` (4,202 records) strata. These are confirmed **cell-free biochemical formats** (measuring in vitro enzymatic kinase inhibition of EGFR on synthetic/peptide substrates, distinct from membrane-bound vesicles or intact cells).
- Because primary literature publications vary in the granularity with which constructs are annotated (some specifying exact catalytic domain residues 668–1210 or 696–1022, while others specify recombinant full-length intracellular domain or purified enzyme), this proportion is strictly described as **"confirmed cell-free biochemical formats"**, rather than asserting universal verification of isolated catalytic domain boundaries.

### 5.3 Dataset Partition Terminology Reference

All dataset partitions and paths strictly adhere to the established Phase 1.3.1.3 and 1.3.1.3.1 nomenclature:

```
data/
├── egfr_unified_dataset.csv                   # 8,944 records (6,414 unique compounds)
├── egfr_conflict_records.csv                  # 1,870 records (421 unique compounds; spread > 1.0 log)
├── egfr_high_confidence.csv                   # 7,053 records (5,970 unique compounds; spread <= 1.0 log)
│   ├── eligible_for_qsar_training == 'True'   # 7,044 records (5,970 unique compounds)
│   ├── is_benchmark_compound == 'True'        # 9 records (2 benchmark leads; excluded from training)
│   └── shares_benchmark_scaffold == 'True'    # 475 records (held out in Protocol D)
└── egfr_assays_by_context/
    ├── purified_recombinant_egfr.csv          # 2,075 records (1,791 eligible; 1,638 compounds)
    ├── other_cell_free_biochemical.csv        # 5,572 records (4,202 eligible; 3,752 compounds)
    ├── cell_free_membrane_derived.csv         # 231 records (169 eligible; 168 compounds; benchmark assay)
    ├── single_protein_terse_metadata.csv      # 1,044 records (882 eligible; 743 compounds)
    └── unresolved_non_standard.csv            # 22 records (0 eligible)
```

---

## 6. Summary of Reconciled Design Specifications

| Design Element | Previous Formulation | Reconciled Formulation (Phase 1.3.2.2) |
| :--- | :--- | :--- |
| **Repeated Measurements** | $w_i = 1 + \ln(M_i) e^{-\Delta}$ (claimed max $\approx 1.5$) | **Strategy A:** Compound-level median target with equal weights ($w_i = 1.0$). Actual formula max of $2.7563$ rejected as scientifically unjustified publication bias. |
| **Primary Training Scope** | All 7,044 eligible records | **Core Homogeneous Subset:** `purified_recombinant_egfr` + `other_cell_free_biochemical` (5,993 records; 5,226 compounds). |
| **Membrane Vesicle Assays** | Included in general training | **Isolated:** Held out as external protocol-transfer evaluation partition (169 records; 168 compounds). |
| **RF Tree Dispersion** | Described alongside epistemic uncertainty | Formally designated as an **uncalibrated heuristic measure of ensemble disagreement**. |
| **Prediction Intervals** | Nominal $\pm 1.645\sigma$, $\pm 1.96\sigma$ intervals | **Split Conformal Prediction** specified for distribution-free coverage guarantees; nominal intervals not claimed as calibrated. |
| **Coverage Curves** | Monotonic decline asserted | Evaluated via Spearman rank correlation and broad percentile deciles (100% to 50%); strict step monotonicity dropped. |
| **Applicability Domain** | Combined uncertainty score | **Decoupled:** Nearest-neighbor Morgan $T_{\max}$ and RF dispersion $\sigma_{\text{trees}}$ reported as orthogonal diagnostics. |
| **M0 Baseline** | Floor where test $R^2 = 0$ | Standard mean-predictor baseline where test $R^2 \le 0$ depending on split mean shift. |
| **Assay Homogeneity Claim** | "85.1% confirmed catalytic domain" | **"85.1% confirmed cell-free biochemical formats"** reflecting actual curation boundaries. |

---

## 7. Sign-Off & Implementation Gate

### **ALL METHODOLOGICAL INCONSISTENCIES RESOLVED**

With this reconciliation, the QSAR experimental design is mathematically verified, biophysically sound, and fully aligned with computational medicinal chemistry standards.

**Phase Gate Status:**
- Zero machine learning models were trained.
- Zero external packages were installed.
- Zero source CSV datasets were modified.
- Zero agent architecture code was altered.

**Awaiting explicit user authorization before proceeding to Phase 1.3.3 implementation.**
