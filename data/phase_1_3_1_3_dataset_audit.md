# PHARMAFORGE — Phase 1.3.1.3: Training-Ready Dataset Preparation Audit

**Target:** Human Epidermal Growth Factor Receptor (EGFR; `CHEMBL203`)  
**Phase Status:** Data Preparation Only — Zero Training / Zero Packages Installed  

---

## 1. Executive Summary & Stratified Architecture

In accordance with Phase 1.3.1.3 specifications, biochemical EGFR assay data have been stratified into **five mutually exclusive, scientifically distinct context datasets** under `data/egfr_assays_by_context/`.
Heterogeneous experimental formats are **not** combined into an unstratified pool. In addition, transparent conflict partitioning has generated dedicated high-confidence and conflict inspection datasets.

### Primary Dataset Accounting

| Dataset Identifier | File Path | Record Count | Unique Compounds | Primary Assay Context & Ontology |
| :--- | :--- | :---: | :---: | :--- |
| **`purified_recombinant_egfr`** | `data/egfr_assays_by_context/purified_recombinant_egfr.csv` | 2,075 | 1,756 | Confirmed purified / recombinant human EGFR catalytic domain (`BAO_0000357`, `BAO_0000219`). |
| **`other_cell_free_biochemical`** | `data/egfr_assays_by_context/other_cell_free_biochemical.csv` | 5,572 | 4,149 | Established cell-free enzymatic kinase assays (radiometric, FRET, HTRF, Caliper, ELISA). |
| **`cell_free_membrane_derived`** | `data/egfr_assays_by_context/cell_free_membrane_derived.csv` | 231 | 227 | Cell-free phosphorylation using membrane preparations/vesicles (e.g. A431 vesicles including `CHEMBL680017`). |
| **`single_protein_terse_metadata`** | `data/egfr_assays_by_context/single_protein_terse_metadata.csv` | 1,044 | 809 | Single-protein format (`BAO_0000357`) with generic literature description (e.g. 'Inhibition of EGFR'). |
| **`unresolved_non_standard`** | `data/egfr_assays_by_context/unresolved_non_standard.csv` | 22 | 22 | Non-standard BAO ontology (`BAO_0000249` membrane fraction) retained via enzymatic regex. |
| **Total Unified Dataset** | *Sum of all 5 strata* | **8,944** | **6,414** | Zero overlap across activity IDs. |
| **High-Confidence Subset** | `data/egfr_high_confidence.csv` | **7,053** | **5,972** | Excludes severe conflicts (> 1.0 log unit) and unresolved contexts. |
| **Conflict Records Subset** | `data/egfr_conflict_records.csv` | **1,870** | **421** | All records for compounds with inter-laboratory spread > 1.0 log unit. |

---

## 2. Correction of Benchmark Assay Classification (CHEMBL680017)

### 2.1 Investigation and Biochemical Context
- **Assay ID:** `CHEMBL680017`
- **Primary Document:** Bridges et al., *J. Med. Chem.* 1996, 39, 267–276 (DOI: 10.1021/jm9503613; `CHEMBL1129291`)
- **Primary Assay Description:** *Concentration required to inhibit the phosphorylation of a 14-residue fragment of phospholipase C gamma 1 by Epidermal growth factor receptor from A431 cell vesicles*
- **BioAssay Ontology:** `BAO_0000219` (*cell-free format*)

### 2.2 Classification Decision & Evidence
1. **Why It Is NOT a Live-Cell Assay:** The experimental procedure in Bridges et al. uses shed plasma membrane vesicles harvested from A431 cells. The vesicles are suspended in cell-free buffer with $\gamma$-[$^{32}\text{P}$]ATP and synthetic peptide substrate (PLC-$\gamma 1$ peptide fragment). There are no intact living cells, no cell proliferation, and no intracellular transport barriers.
2. **Why It Was Previously Excluded:** In Phase 1.3.1, `curate_chembl_egfr.py` used regex `CELLULAR_REGEX` which matched the substring `'A431 cell'`. Because it also contained `'phosphorylation'`, the script assigned it to 'Cellular Autophosphorylation'.
3. **Correction:** `CHEMBL680017` is now correctly classified under **`cell_free_membrane_derived`**. All 64 records from `CHEMBL680017` (including benchmark activity `996053` for `CHEMBL289959` and `1016106` for `CHEMBL290096`) are fully integrated into `data/egfr_assays_by_context/cell_free_membrane_derived.csv`.

---

## 3. Transparent Handling of Conflicting Measurements

In accordance with Phase 1.3.1.3 guidelines:
- **Zero Data Destruction:** Conflicting measurements are **not** replaced with medians or means. Every raw observation is preserved with its exact `activity_id`, `document_chembl_id`, and `standard_value`.
- **Explicit Conflict Metadata:** Columns `conflict_status`, `measurement_count`, and `spread_log10` are embedded in every record.
- **Inter-Laboratory Variance Recognition:** Spread $> 1.0 \log$ unit is recognized as protocol diversity (differences in ATP concentration relative to $K_m$, construct lengths, incubation buffers) across independent academic and industrial laboratories, rather than simple experimental error.

### Conflict Distribution Across Unique Compounds

- **Single-Observation Compounds:** 5,149 compounds (80.3%)
- **Concordant Multi-Observation Compounds (Spread <= 1.0 log unit):** 844 compounds (13.2%)
- **Severe-Conflict Compounds (Spread > 1.0 log unit):** 421 compounds (6.6%)

---

## 4. Preservation of Uncertainty Regarding Wild-Type Status

- **Wild-Type Field:** Every record in the stratified datasets contains `wild_type_status = 'inferred_unannotated'`.
- **Scientific Rationale:** All 7,815 records with explicit mutation annotations (`assay_variant_mutation`) or mutant mentions in descriptions were purged during Phase 1.3.1. However, the absence of mutation annotations does not constitute formal proof of wild-type construct. Designating these as `inferred_unannotated` prevents misleading future modeling phases.

---

## 5. Independent Verification of Benchmark Scaffold Identity

RDKit scaffold analysis was independently executed on `CHEMBL289959` and `CHEMBL290096`:

| Compound Identifier | Original Canonical SMILES | Standardized SMILES | Bemis–Murcko Scaffold | Generic Scaffold of Bemis–Murcko Scaffold |
| :--- | :--- | :--- | :--- | :--- |
| **CHEMBL289959** | `c1ccc(Nc2ncnc3ccccc23)cc1` | `c1ccc(Nc2ncnc3ccccc23)cc1` | `c1ccc(Nc2ncnc3ccccc23)cc1` | `C1CCC(CC2CCCC3CCCCC32)CC1` |
| **CHEMBL290096** | `Brc1cccc(Nc2ncnc3ccccc23)c1` | `Brc1cccc(Nc2ncnc3ccccc23)c1` | `c1ccc(Nc2ncnc3ccccc23)cc1` | `C1CCC(CC2CCCC3CCCCC32)CC1` |

### Verification Findings
1. **Scaffold Identity:** The calculated Bemis–Murcko scaffolds are **identical** (`c1ccc(Nc2ncnc3ccccc23)cc1`; 4-anilinoquinazoline).
2. **Generic Scaffold Verification:** When RDKit's `MakeScaffoldGeneric` is applied to the calculated Bemis–Murcko scaffold, the resulting generic framework is **`C1CCC(CC2CCCC3CCCCC32)CC1`** for both compounds. (Note: Running `MakeScaffoldGeneric` directly on the halogenated molecule `CHEMBL290096` without prior scaffold extraction preserves the bromine as a methyl carbon `C`, demonstrating why Bemis–Murcko stripping must precede generic framework conversion).

---

## 6. Benchmark Leakage Prevention & Exclusion Rules

### 6.1 Compound-Level Benchmark Exclusion
- **Rule 1:** `CHEMBL289959` and `CHEMBL290096` (and standardized SMILES `c1ccc(Nc2ncnc3ccccc23)cc1` and `Brc1cccc(Nc2ncnc3ccccc23)c1`) must be excluded from the QSAR training partition (`is_benchmark_compound == 'True'`).

### 6.2 Scaffold-Level Generalization Exclusion
- **Rule 2:** To evaluate true zero-shot out-of-scaffold generalization, all **372 unique compounds (696 records)** in the dataset that share the 4-anilinoquinazoline Bemis–Murcko scaffold (`shares_benchmark_scaffold == 'True'`) must be held out in the test partition.
- **Scope Constraint:** This two-compound benchmark evaluates single-point halogen optimization ($R = \text{H}$ vs $R = \text{Br}$) targeting the hydrophobic pouch near gatekeeper Thr790. It does **not** demonstrate broad scaffold generalization across alternative kinase hinge-binding motifs.

---

## 7. Known Limitations & Recommendations for Phase 1.3.2

1. **Stratum-Specific QSAR Modeling:** QSAR models should be benchmarked primarily on the `purified_recombinant_egfr` and `other_cell_free_biochemical` strata, or on `egfr_high_confidence.csv`.
2. **Inferred Wild-Type Uncertainty:** Wild-type status remains inferred. Future validation against prospective assays should verify target construct sequencing.
3. **Protocol Variation:** The 421 severe-conflict compounds should be held out of benchmark evaluation or modeled using distribution/uncertainty intervals rather than scalar regression targets.

---

## 8. Final Verdict

### **PASS — DATASETS PREPARED AND AUDITED**

All required datasets, stratified directories, high-confidence subsets, conflict partitions, and metadata annotations are generated and verified. Phase 1.3.1.3 is complete.
