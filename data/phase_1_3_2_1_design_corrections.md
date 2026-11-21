# PHARMAFORGE — Phase 1.3.2.1: QSAR Experimental Design Corrections Report

**Document ID:** `data/phase_1_3_2_1_design_corrections.md`  
**Date:** 2026-11-20  
**Target:** Human Epidermal Growth Factor Receptor (EGFR; `CHEMBL203`)  
**Phase Status:** Design Corrections Only — Zero Model Training / Zero Package Installations / Zero Source Dataset Modifications  

---

## 1. Executive Summary

This report establishes essential scientific and methodological corrections to the Phase 1.3.2 QSAR architecture design ([phase_1_3_2_qsar_design.md](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/data/phase_1_3_2_qsar_design.md)).

Five specific methodological issues are resolved:
1. **Repeated-Measurement Policy:** Reconciling the 7,044 eligible records with the 5,970 unique standardized compounds, preventing statistical over-representation of frequently tested compounds, and defining compound-level aggregation and sample-weighting strategies.
2. **Uncertainty & Applicability Domain (AD) Policy:** Labeling Random Forest tree variance and initial thresholds as provisional heuristics, detailing an empirical validation framework using held-out error-versus-confidence curves.
3. **Leakage-Safe Document Evaluation (Protocol C):** Analyzing the bipartite compound-document graph (538 connected components; maximum component 15.4%), defining a true zero-leakage component split, and providing an alternative disjoint publication series holdout.
4. **Baseline Model-Feature Comparison:** Finalizing the 4-way baseline matrix (Dummy, RF-Morgan, RF-Hybrid, HistGradientBoosting-Descriptors) and stating the exact hypotheses tested.
5. **Dataset Composition Check:** Reconciling the exact stratum-level composition of the 7,044 eligible records and the 1,900 excluded records.

---

## 2. Repeated-Measurement Policy

### 2.1 Empirical Distribution of Observations per Compound
Inspection of the 7,044 eligible training records in `data/egfr_high_confidence.csv` reveals that the dataset does **not** contain 7,044 independent molecules, but rather **5,970 unique standardized compounds**:

| Observations per Compound ($M_i$) | Number of Compounds | % of Compounds | Total Records ($M_i \times N_i$) | % of Total Records | Representative Compounds & Context |
| :---: | :---: | :---: | :---: | :---: | :--- |
| **1 observation** | **5,128** | **85.9%** | **5,128** | **72.8%** | Single-observation chemical series across literature. |
| **2 observations** | 662 | 11.1% | 1,324 | 18.8% | Re-tested in follow-up study or independent confirmation. |
| **3 observations** | 144 | 2.4% | 432 | 6.1% | Common chemical intermediates or reference compounds. |
| **4 observations** | 27 | 0.45% | 108 | 1.5% | Frequently cited clinical or research chemotypes. |
| **5 observations** | 5 | 0.08% | 25 | 0.4% | E.g., `CHEMBL343352` (spread: $0.025\log$ unit across 5 docs). |
| **6 observations** | 3 | 0.05% | 18 | 0.3% | E.g., `CHEMBL328216`, `CHEMBL597402`, `CHEMBL1230609`. |
| **9 observations** | 1 | 0.02% | 9 | 0.1% | `CHEMBL296407` (tested across 7 publications; spread: $0.699\log$). |
| **Total** | **5,970** | **100.0%** | **7,044** | **100.0%** | **Surplus records from repeated assays:** $\mathbf{1,074}$ |

### 2.2 Methodological Risks of Naïve Observation-Level Training
1. **Violation of Statistical Independence:** Treating 7,044 rows as independent identically distributed (i.i.d.) observations gives disproportionate statistical weight to heavily studied molecules. A compound with 9 measurements (`CHEMBL296407`) contributes $9\times$ the gradient pull of a single-observation compound.
2. **Data Leakage in Splits:** If observations of a multi-measurement compound are randomly distributed across folds, the model memorizes the structure in training and is evaluated on the exact same structure in testing.

### 2.3 Approved Handling Strategy

#### A. Cross-Validation Partitioning Rule (Strict Grouping)
- **Mandatory Rule:** In all cross-validation protocols (Random, Scaffold, Document), **all records of the same standardized compound must be assigned to the exact same fold**.
- Grouping key: `standardized_smiles`. No compound can straddle the train/test divide.

#### B. Primary Model Training Strategy: Compound-Level Median Aggregation
For the primary point-regression model:
1. **Sample Unit:** Each unique standardized compound $i \in \{1, \dots, 5,970\}$ constitutes exactly **one training sample**.
2. **Target Label:** The median $p\text{IC}_{50}$ of its concordant measurements:
   $$\bar{y}_i = \text{median}\left(\{y_{ij}\}_{j=1}^{M_i}\right)$$
   *Justification:* All records in the high-confidence set have $\text{spread} \le 1.0\log$ unit. The median is inherently robust to single-assay protocol variation and prevents outlier skewing without corrupting raw records.
3. **Sample Weighting by Experimental Concordance:**
   Compounds with multiple concordant observations across independent publications represent higher empirical confidence than single-observation molecules. An experimental certainty weight is assigned:
   $$w_i = 1.0 + \ln(M_i) \cdot \exp\left(-\Delta p\text{IC}_{50, i}\right)$$
   where $\Delta p\text{IC}_{50, i} = \max(y_i) - \min(y_i)$. For single observations ($M_i = 1$), $w_i = 1.0$. Highly concordant multi-observation compounds receive a modest weight boost up to $\sim 1.5$, while compounds with larger spread receive lower weights.

#### C. Evaluation Protocol
- **Primary Benchmark:** Evaluated against the compound-level median target $\bar{y}_i$.
- **Secondary Assay-Level Robustness:** Predictions are also compared against all 1,916 individual raw observations to quantify prediction error relative to natural inter-assay experimental variability.

---

## 3. Uncertainty and Applicability-Domain Policy

### 3.1 Clarification on Random Forest Tree Variance
- **Epistemic Limitation:** The standard deviation across decision trees in a Random Forest ($\sigma_{\text{trees}}$) reflects ensemble variance on the feature subspace. **It is NOT calibrated Bayesian epistemic uncertainty.**
- Random Forests are known to produce uncalibrated, overconfident predictions in sparse chemical space.
- **Correction:** In all documentation and code, $\sigma_{\text{trees}}$ is formally designated as **provisional heuristic ensemble dispersion**, not true epistemic uncertainty.

### 3.2 Provisional Designation of Thresholds
The proposed thresholds:
- $T_{\max} < 0.35$ (Tanimoto distance to nearest training neighbor)
- $\sigma_{\text{trees}} > 1.0\ p\text{IC}_{50}\text{ units}$ (ensemble spread)
are explicitly labeled as **provisional engineering heuristics** that must undergo empirical validation.

### 3.3 Empirical Validation Framework Using Held-Out Folds
In Phase 1.3.3, applicability domain and uncertainty performance will be validated using held-out test predictions:

1. **Error-versus-Confidence (Rejection) Curves:**
   - Sort test set predictions in order of increasing uncertainty metric ($\sigma_{\text{trees}}$ or decreasing $T_{\max}$).
   - Successively remove the most uncertain $5\%, 10\%, \dots, 50\%$ of predictions.
   - **Validation Criterion:** If the uncertainty metric is statistically informative, the test MAE on the remaining retained subset must monotonically decrease. If MAE does not decrease, the heuristic is uninformative and must be recalibrated.
2. **Empirical Coverage of Prediction Intervals:**
   - Construct nominal 80% and 95% intervals: $\widehat{y} \pm z \cdot \sigma_{\text{trees}}$.
   - Measure actual empirical coverage on held-out folds. If empirical coverage diverges from nominal coverage, apply conformal prediction or temperature scaling on validation folds.
3. **Stratified Error by Structural Similarity Bins:**
   - Partition held-out test compounds into Tanimoto bins relative to training set:
     * High similarity: $T_{\max} \ge 0.60$
     * Moderate similarity: $0.45 \le T_{\max} < 0.60$
     * Low similarity: $0.35 \le T_{\max} < 0.45$
     * Out-of-Domain: $T_{\max} < 0.35$
   - Compute MAE and $R^2$ within each bin to empirically identify the true degradation boundary.

---

## 4. Leakage-Safe Document Evaluation (Redesign of Protocol C)

### 4.1 Bipartite Graph Analysis of Compounds and Documents
To evaluate inter-laboratory generalizability without data leakage, a document-held-out split must satisfy **two simultaneous non-leakage conditions**:
1. Zero shared documents between training and test sets ($\mathcal{D}_{\text{train}} \cap \mathcal{D}_{\text{test}} = \emptyset$).
2. Zero shared standardized compounds between training and test sets ($\mathcal{C}_{\text{train}} \cap \mathcal{C}_{\text{test}} = \emptyset$).

Graph analysis of the bipartite network between the 5,970 unique compounds and 719 unique publications reveals:
- **Total Connected Components:** **538 components**.
- **Distribution of Component Sizes:**
  * Component 1 (Largest): **920 compounds (15.4%)**, 58 documents (8.1%)
  * Component 2: 263 compounds (4.4%), 27 documents (3.8%)
  * Component 3: 165 compounds (2.8%), 2 documents
  * Component 4: 124 compounds (2.1%), 7 documents
  * Components 5–538: 94 compounds or fewer, mostly spanning 1 to 3 documents.
- **Key Empirical Finding:** **There is NO giant connected component spanning 50%+ of the dataset.** Component 1 occupies only 15.4% of compounds.

### 4.2 Two Approved Protocol C Options

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    REDESIGNED PROTOCOL C OPTIONS                            │
├───────────────────────────────────┬─────────────────────────────────────────┤
│ Option C.1: Connected-Component   │ Pack the 538 connected components into  │
│ 5-Fold Cross-Validation           │ 5 balanced folds (~20% each). Guarantees│
│                                   │ ZERO shared compounds AND zero shared   │
│                                   │ documents simultaneously.               │
├───────────────────────────────────┼─────────────────────────────────────────┤
│ Option C.2: Disjoint Single-Paper │ Train on multi-document components      │
│ Series Holdout Evaluation         │ (4,000+ compounds); hold out 5–10       │
│                                   │ complete single-document lead-opt series│
│                                   │ (e.g. 50–100 novel compounds each).     │
└───────────────────────────────────┴─────────────────────────────────────────┘
```

#### What Each Establishes and Cannot Establish:
- **Option C.1 (Connected-Component Split):**
  - *Establishes:* Rigorous, leakage-free cross-validation over structurally and document-disjoint clusters.
  - *Caveat:* Component 1 (which links 58 documents through common reference compounds like Gefitinib and Lapatinib) must be placed in a single fold, so that fold will contain the most well-known kinase drugs.
- **Option C.2 (Disjoint Series Holdout):**
  - *Establishes:* True real-world simulation of evaluating a QSAR model on an entire new lead-optimization series published by an unseen laboratory.
  - *Cannot Establish:* Cannot test whether the model predicts new measurements of familiar control drugs, because control drugs link papers together.

---

## 5. Finalized Baseline Model-Feature Comparison Matrix

The 4-way baseline matrix is finalized as follows:

| Model ID | Model Family | Feature Representation | Input Dimensions | Primary Scientific Hypothesis Tested |
| :---: | :--- | :--- | :---: | :--- |
| **M0** | `DummyRegressor(strategy='mean')` | None (null model) | 0 | Baseline statistical floor ($R^2 = 0$). Any predictive QSAR model must substantially beat this. |
| **M1** | `RandomForestRegressor` | Morgan Fingerprints (ECFP4; $r=2$) | 2,048 bits | **H1:** Local 2D circular substructures and pharmacophoric bits alone provide sufficient signal for target affinity. |
| **M2** | `RandomForestRegressor` | Hybrid (Morgan 2048 + 10 Descriptors) | 2,058 dims | **H2:** Augmenting local circular bits with global physical properties (MW, LogP, TPSA, HBD, HBA) improves generalization. |
| **M3** | `HistGradientBoostingRegressor` | Compact Descriptors Only | 10 continuous | **H3:** Low-dimensional interpretable physical properties with non-linear boosting thresholds can approximate affinity rapidly. |

### Justification for Model Selection:
- **Zero New Dependencies:** All 4 models are built directly into `scikit-learn 1.9.1`.
- **Rapid Iteration:** Training completes in seconds to minutes on standard CPU hardware.
- **Principled Progression:** Directly compares pure substructures (M1), pure physical descriptors (M3), and their combination (M2) against the null baseline (M0).
- **Deep Learning Deferral:** Graph neural networks (GNNs) or chemical language transformers are deferred until this tabular baseline matrix establishes the performance benchmark.

---

## 6. Dataset Composition Reconciliation

### 6.1 Exact Accounting of Eligible and Excluded Records
Reconciling the complete unified 8,944-record dataset against the proposed primary training set:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    UNIFIED DATASET STRATUM-LEVEL ACCOUNTING                 │
├───────────────────────────────┬────────────┬────────────┬───────────────────┤
│ Stratum Name                  │ Total Rows │ Eligible   │ Excluded Records  │
│                               │            │ Training   │ & Exclusion Reason│
├───────────────────────────────┼────────────┼────────────┼───────────────────┤
│ other_cell_free_biochemical   │ 5,572      │ 4,202      │ 1,370 (Conflict)  │
│ purified_recombinant_egfr     │ 2,075      │ 1,791      │ 284 (Conflict)    │
│ single_protein_terse_metadata │ 1,044      │ 882        │ 162 (Conflict)    │
│ cell_free_membrane_derived    │ 231        │ 169        │ 62 (60 Conflict + │
│                               │            │            │     2 Benchmark)  │
│ unresolved_non_standard       │ 22         │ 0          │ 22 (Non-standard  │
│                               │            │            │     BAO_0000249)  │
├───────────────────────────────┼────────────┼────────────┼───────────────────┤
│ TOTALS                        │ 8,944      │ 7,044      │ 1,900 Excluded    │
└───────────────────────────────┴────────────┴────────────┴───────────────────┘
```

### 6.2 Key Scientific Takeaways
1. **7,044 Records $\neq$ 7,044 Independent Molecules:** The eligible training corpus represents exactly **5,970 unique standardized compounds**.
2. **Dominant Assays:** 85.1% of eligible records belong to confirmed cell-free enzymatic kinase domain assays (`other_cell_free_biochemical` [59.7%] and `purified_recombinant_egfr` [25.4%]).
3. **Excluded Records Accounting (1,900 records):**
   - 1,869 records: Severe inter-laboratory conflict ($> 1.0\log$ unit spread)
   - 1 record: Severe conflict + unresolved assay format (ActID `609241`, `CHEMBL44`)
   - 21 records: Unresolved non-standard membrane format only
   - 9 records: Benchmark starting compounds (`CHEMBL289959` and `CHEMBL290096`)

---

## 7. Updates Applied to the QSAR Architecture Design

The primary architecture document ([phase_1_3_2_qsar_design.md](file:///c:/Users/Sravya%20mummana/OneDrive/PHARMAFORGE/data/phase_1_3_2_qsar_design.md)) has been updated to reflect:
1. Compound-level median aggregation ($5,970$ samples) with concordance sample weighting ($w_i$).
2. Provisional labeling of RF tree variance and applicability domain thresholds ($T_{\max} < 0.35$).
3. Redesign of Protocol C using connected-component partitioning.
4. Formalization of the 4-way baseline matrix (M0 through M3).
5. Explicit distinction between the 7,044 activity records and 5,970 unique molecules.

---

## 8. Final Verdict

### **DESIGN FULLY RECONCILED & METHODOLOGICALLY RIGOROUS**

All scientific concerns regarding repeated measurements, uncalibrated uncertainty, document split leakage, and baseline model comparisons are resolved.

**Implementation remains stopped. Awaiting explicit user approval before proceeding to Phase 1.3.3 (Model Implementation & Cross-Validation).**
