# PHARMAFORGE — Phase 1.3.1.2: Numerical and Scientific Reconciliation Report

**Date Generated:** 2026-11-20 20:07:00 UTC  
**Target:** Human Epidermal Growth Factor Receptor (EGFR; `CHEMBL203`)  
**Pipeline Script Audited:** `scripts/curate_chembl_egfr.py`  
**Datasets Audited:** `data/egfr_raw_activities.csv`, `data/egfr_curated_activities.csv`, `data/benchmark_lead_provenance.csv`  

---

## 1. Exact Numerical Reconciliation of the Curation Funnel

The curation pipeline in `scripts/curate_chembl_egfr.py` executes a strictly sequential single-pass funnel. The table below presents the exact stage-by-stage accounting where **every stage's input count strictly equals the previous stage's output count**, and every exclusion is mutually exclusive within the sequence.

### 1.1 Stage-by-Stage Accounting Table

| Stage | Filter Name | Input Records | Excluded Records | Output Records | Specific Exclusion Reason & Category Breakdown |
| :---: | :--- | :---: | :---: | :---: | :--- |
| **0** | **Raw API Extraction** | — | — | **18,381** | Query: `target_chembl_id=CHEMBL203`, `standard_type=IC50`, `units=nM`, `relation='='`, `assay_type=B`. |
| **1** | **Valid Numeric Endpoint** | 18,381 | 1 | **18,380** | Excluded ActID `740166` (`standard_value=''` empty string; unparseable float). |
| **2** | **Target Organism Specificity** | 18,380 | 0 | **18,380** | 100% of records confirmed `target_organism == 'Homo sapiens'`. |
| **3** | **Strict Variant / Mutant Exclusion** | 18,380 | 7,815 | **10,565** | Excluded mutant constructs: <br>• Both `assay_variant_mutation` & `MUTATION_REGEX`: 4,252 <br>• `MUTATION_REGEX` alone (free text): 3,549 <br>• `assay_variant_mutation` alone: 14 |
| **4** | **Assay Context Validation** | 10,565 | 1,677 | **8,888** | Excluded non-biochemical / cellular formats: <br>• Cellular Autophosphorylation: 1,490 <br>• Cell Viability / Intact Proliferation: 170 <br>• Ambiguous / Unclear Context: 17 |
| **5** | **Chemical Standardization & Desalting** | 8,888 | 8 | **8,880** | Excluded records with missing `canonical_smiles`: 8 records (ActIDs: `7601338`, `7621398`, `7666254`, `7679918`, `7690831`, `7693413`, `7753173`, `7790477`). |

$$\text{Final Curated Records} = 18,381 - (1 + 0 + 7,815 + 1,677 + 8) = 18,381 - 9,501 = \mathbf{8,880}$$

### 1.2 Reconciliation of 10,566 vs 10,565 vs 10,557 Discrepancy

The discrepancy between these numbers arises from evaluating mutation annotations at different points relative to endpoint validity and chemical structure completeness:

1. **10,566 Records:** The total number of records in the **entire raw dataset** (18,381 records) that have no mutation annotation (`variant_id` null, `assay_variant_mutation` empty, and no regex match in `assay_description`):
   $$18,381 - 7,815\text{ (mutants)} = \mathbf{10,566}$$
2. **10,565 Records:** The number of unmutated records **after Stage 1 valid numeric endpoint filtering**. Exactly 1 record (ActID `740166`, standard_value = `''`) had no mutation annotation but failed Stage 1:
   $$10,566 - 1\text{ (invalid endpoint)} = \mathbf{10,565}$$
3. **10,557 Records:** The number of unmutated records that have both valid numeric endpoints **and** valid canonical SMILES, evaluated prior to assay-context filtering:
   $$10,565 - 8\text{ (missing SMILES)} = \mathbf{10,557}$$

All three numbers are mathematically reconciled and account for every record in the dataset.

---

## 2. Severe Assay Conflict Reconciliation (420 vs 303)

### 2.1 Origin and Clarification of the Counts
- **Official Pipeline Count (420 Compounds):** In `scripts/curate_chembl_egfr.py` and `data/curation_report.md` (lines 72–74), the grouping key is `standardized_smiles` (desalted, canonical RDKit SMILES). Exactly **420 unique compounds** possess measurements spanning $> 1.0 \log$ unit across different assays/publications. These 420 compounds account for **1,861 activity rows** labeled `SEVERE_ASSAY_CONFLICT_GT_1LOG` in `egfr_curated_activities.csv`.
- **Origin of 303:** The count of 303 in the previous narrative summary was an **accidental drafting error** in the prose text of Section 5. The underlying database files, curation scripts, and raw execution logs have consistently contained 420 compounds. 

### 2.2 Conflict Recalculation Across Defined Grouping Keys

| Grouping Key | Unique Entities | Single-Observation Entities | Concordant Entities ($\text{Spread} \le 1.0 \log$) | Severe Conflict Entities ($\text{Spread} > 1.0 \log$) | Activity Rows Flagged as Severe Conflict |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`standardized_smiles` (Official Pipeline)** | **6,398** | **5,158 (80.6%)** | **820 (12.8%)** | **420 (6.6%)** | **1,861 rows** |
| `canonical_smiles` (Pre-desalting) | 6,435 | 5,213 (81.0%) | 801 (12.4%) | 421 (6.5%) | 1,862 rows |
| `molecule_chembl_id` | 6,435 | 5,213 (81.0%) | 801 (12.4%) | 421 (6.5%) | 1,862 rows |
| `standardized_smiles` within identical BAO format | 6,398 | 5,158 (80.6%) | 916 (14.3%) | 324 (5.1%) | 1,448 rows |

*Note:* When conflict analysis is restricted strictly within identical BAO ontology terms (e.g., comparing `BAO_0000357` only against `BAO_0000357`), 324 compounds remain in severe conflict ($> 1.0 \log$ unit), indicating genuine inter-laboratory protocol variance (varying ATP concentrations, incubation times, enzyme constructs) rather than ontology mismatch.

---

## 3. Refined Assay-Context Interpretation

It is scientifically inaccurate to assume all 8,880 retained records represent homogeneous recombinant kinase domain assays. A granular audit of BAO ontology terms and primary assay descriptions reveals four distinct experimental categories:

### 3.1 Granular Breakdown of Retained Activities (8,880 Records)

| Assay Classification Category | Record Count | % of Curated Set | BioAssay Ontology (BAO) Terms | What This Establishes | What This Does NOT Establish |
| :--- | :---: | :---: | :--- | :--- | :--- |
| **Other Cell-Free Biochemical Assays** | **5,572** | **62.75%** | `BAO_0000357`, `BAO_0000219` | Established cell-free enzymatic kinase activity (radiometric, FRET, HTRF, Caliper mobility shift, or ELISA on exogenous substrates). | Specific construct length, phosphorylation state, or ATP concentration relative to $K_{m,\text{ATP}}$. |
| **Purified / Recombinant Protein Assays** | **2,075** | **23.37%** | `BAO_0000357`, `BAO_0000219` | Explicit confirmation of purified or recombinant human EGFR (e.g. Sf9 insect cell expression, GST/His-tagged catalytic domain). | Homogeneous ATP concentration across literature sources. |
| **Single-Protein Ontology with Terse Metadata** | **1,044** | **11.76%** | `BAO_0000357` | Classified as single-protein format by ChEMBL curators; confirms non-cellular format. | Description is generic (e.g., `"Inhibition of EGFR"`); detailed experimental protocol cannot be reconstructed from text alone. |
| **Cell-Free Membrane-Derived Assays** | **167** | **1.88%** | `BAO_0000357`, `BAO_0000219` | Cell-free kinase phosphorylation using membrane-bound EGFR preparations (e.g. isolated from A431 membranes). | Not purified soluble enzyme; contains membrane lipids and full-length EGFR. |
| **Unresolved / Non-Standard BAO** | **22** | **0.25%** | `BAO_0000249` (membrane format) | Cell-free phosphorylation confirmed by text regex, but carrying legacy membrane ontology tag. | Standardized single-protein format. |

### 3.2 Deep Investigation of Benchmark Assay CHEMBL680017
The selected benchmark assay `CHEMBL680017` is described in Bridges et al. 1996 (*J. Med. Chem.*, DOI: 10.1021/jm9503613) as:
> *"Concentration required to inhibit the phosphorylation of a 14-residue fragment of phospholipase C gamma 1 by Epidermal growth factor receptor from A431 cell vesicles"*

**Biochemical and Engineering Audit:**
1. **Biochemical Reality:** This is a **cell-free, membrane-derived in vitro kinase assay**. Plasma membrane vesicles shed from human A431 carcinoma cells overexpress wild-type EGFR. The vesicles are incubated with exogenous synthetic peptide substrate (PLC-$\gamma 1$ fragment) and $\gamma$-[$^{32}\text{P}$]ATP. There are no living cells, no cell permeability barriers, and no cellular metabolic clearance. In ChEMBL, it is classified under `BAO_0000219` (*cell-free format*).
2. **Exclusion from Curated Set:** In `scripts/curate_chembl_egfr.py`, `CELLULAR_REGEX` matched the substring `"A431 cell"`, while `BIOCHEMICAL_REGEX` matched `"phosphorylation"`. Under the script's decision tree, any record matching both cellular and biochemical patterns was flagged as `"Cellular Autophosphorylation"` and excluded. Consequently, all 64 records from `CHEMBL680017` were excluded from `egfr_curated_activities.csv`.
3. **Preservation:** The benchmark measurements for `CHEMBL289959` and `CHEMBL290096` from `CHEMBL680017` are fully preserved in `data/benchmark_lead_provenance.csv` and `data/egfr_raw_activities.csv`. Both compounds also have alternative concordant cell-free measurements present directly inside `egfr_curated_activities.csv` (e.g., ActIDs `894026` [344 nM] and `894021` [27 nM] from Ward et al. 1995, `CHEMBL677204`).

---

## 4. Benchmark Comparability Verification

### 4.1 Verification of Primary Benchmark Records
Both benchmark compound records have been verified against `data/benchmark_lead_provenance.csv` and `data/egfr_raw_activities.csv`:

| Parameter | Compound 1 (CHEMBL289959) | Compound 2 (CHEMBL290096) | Status |
| :--- | :---: | :---: | :---: |
| **Activity ID** | `996053` | `1016106` | **Verified Present** |
| **Assay ID** | `CHEMBL680017` | `CHEMBL680017` | **Identical** |
| **Document ID** | `CHEMBL1129291` | `CHEMBL1129291` | **Identical (Bridges et al. 1996)** |
| **Standard Type** | `IC50` | `IC50` | **Identical** |
| **Standard Relation** | `=` | `=` | **Identical** |
| **Standard Value** | `344.0` | `27.0` | **Exact Raw Values Preserved** |
| **Standard Units** | `nM` | `nM` | **Identical** |
| **pChEMBL Value** | `6.46` | `7.57` | **Exact Derived Values Preserved** |

### 4.2 Limitations of Shared Bemis-Murcko Scaffolds
RDKit Murcko scaffold analysis confirms that `CHEMBL289959` and `CHEMBL290096` share the **exact identical Bemis-Murcko scaffold**:
$$\text{Scaffold: } \mathbf{c1ccc(Nc2ncnc3ccccc23)cc1}\text{ (N-phenylquinazolin-4-amine)}$$
$$\text{Generic Framework: } \mathbf{C1CCC(CC2CCCC3CCCCC32)CC1}$$

**Scientific Limitations:**
1. **Lack of Chemotype Diversity:** Using both compounds as separate benchmark tasks evaluates optimization solely within the 4-anilinoquinazoline chemotype. It does not test whether the multi-agent system can optimize alternative kinase hinge-binding motifs (e.g., pyrrolopyrimidines, indoles, pyrazolopyrimidines).
2. **Scaffold-Split Data Leakage:** In a strict Bemis-Murcko scaffold split for QSAR validation, all molecules containing `c1ccc(Nc2ncnc3ccccc23)cc1` must reside in the same data split. If this scaffold is in the training set, the QSAR model has memorized the core; if it is held out in the test set, the model must predict both strictly by zero-shot out-of-scaffold generalization.
3. **Matched Molecular Pair (MMP) Utility:** The single structural difference is a bromine atom at the meta position of the aniline ring ($R = \text{H}$ vs $R = \text{Br}$), which fills the hydrophobic pocket adjacent to gatekeeper Thr790. This serves as an excellent test of whether the agent understands subtle halogenation SAR, but it represents a single chemical series rather than two independent starting points.

---

## 5. Remaining Scientific Limitations

1. **Inferred Wild-Type Status:** 100% of curated records are designated wild-type based on the absence of mutant annotations. This is an unavoidable limitation of public chemoinformatics databases.
2. **Inter-Laboratory Protocol Variance:** 420 compounds (6.6%) have conflicting IC50 values spanning $> 1.0 \log$ unit across different publications. These records are explicitly flagged in `conflict_status` and must be handled conservatively during QSAR model development.
3. **Filter 4 Regex False Exclusions:** Keyword matching on `"A431 cell"` excluded membrane vesicle assays like `CHEMBL680017` that are cell-free in execution. While this conservative filter preserved 8,880 high-confidence biochemical records, future pipeline refinements could allow `BAO_0000219` to override cell mentions when `"vesicle"` or `"cell-free"` is explicitly annotated.

---

## 6. Final Verdict

### **PASS WITH LIMITATIONS**

**Justification:**
1. **Mathematical Reconciliation:** The stage-by-stage funnel accounting is strictly verified. Every stage input matches the preceding output, and the final count reconciles exactly to 8,880.
2. **Conflict Reconciliation:** The 420 vs 303 discrepancy is resolved: 420 is the true, verified compound count in the dataset and report, while 303 was an acknowledged narrative drafting typo.
3. **Benchmark Verification:** Primary records for `CHEMBL289959` (ActID `996053`) and `CHEMBL290096` (ActID `1016106`) are verified with 100% provenance integrity, identical assays, identical publications, and preserved raw values.
4. **Limitations Documented:** Scaffold redundancy and inferred wild-type classifications are fully documented with clear recommendations for Phase 1.3.2.

**Implementation is paused. Awaiting explicit user approval before proceeding to Phase 1.3.2.**
