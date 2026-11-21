# PHARMAFORGE — Phase 1.3.1: ChEMBL EGFR Data Curation & Provenance Report

**Date Generated:** 2026-11-20 14:24:00 UTC  
**Target:** Human Epidermal Growth Factor Receptor (EGFR; `CHEMBL203`)  
**Primary Objective:** Audit biochemical assay provenance and establish a leakage-resistant dataset.

---

## 1. Executive Summary & Provenance Metrics

- **Total Raw Activities Retrieved:** 18,381
- **Curated Activities Retained:** 8,880
- **Unique Standardized Compounds:** 6,398
- **Unique Assays in Curated Dataset:** 1,070
- **Unique Documents in Curated Dataset:** 804
- **Total Variant/Mutant Records Excluded:** 7,815

---

## 2. Multi-Stage Curation Funnel

| Stage | Filter Condition | Retained Records | Cumulative Retention |
| :--- | :--- | :--- | :--- |
| **Raw API Query** | `target_chembl_id=CHEMBL203`, `standard_type=IC50`, `units=nM`, `relation='='`, `assay_type=B` | 18,381 | 100.0% |
| **Filter 1: Valid Endpoint** | Positive numerical `standard_value > 0.0` | 18,380 | 100.0% |
| **Filter 2: Organism Specificity** | Confirmed `Homo sapiens` | 18,380 | 100.0% |
| **Filter 3: Wild-Type Verification** | Excluded `variant_id IS NOT NULL`, `assay_variant_mutation`, and mutant regex | 10,565 | 57.5% |
| **Filter 4: Biochemical Context** | Excluded intact-cell / functional assays; retained purified kinase domain | 8,888 | 48.4% |
| **Filter 5: Chemical Standardization** | RDKit sanitization, salt removal, canonical SMILES | 8,880 | 48.3% |

---

## 3. Variant and Mutant Exclusions

A total of **7,815 records** were identified as involving mutant EGFR constructs (e.g. T790M, L858R, C797S, G719S, Exon 19 deletions) and were purged to maintain wild-type kinase domain specificity.

Representative excluded mutant records:

| Activity ID | Molecule ID | Variant Annotation | Reason for Exclusion | Assay Description Preview |
| :--- | :--- | :--- | :--- | :--- |
| `2136929` | `CHEMBL387187` | `None` | mutation=L858R; mutant_text='L858R' | Inhibition of EGFR L858R mutant |
| `2136930` | `CHEMBL387187` | `None` | mutation=L861Q; mutant_text='L861Q' | Inhibition of EGFR L861Q mutant |
| `2271252` | `CHEMBL249496` | `None` | mutation=L858R; mutant_text='L858R' | Inhibition of EGFR L858R mutant by time resolved fluorescence assay |
| `2271255` | `CHEMBL248877` | `None` | mutation=L858R; mutant_text='L858R' | Inhibition of EGFR L858R mutant by time resolved fluorescence assay |
| `2271258` | `CHEMBL251732` | `None` | mutation=L858R; mutant_text='L858R' | Inhibition of EGFR L858R mutant by time resolved fluorescence assay |
| `2271261` | `CHEMBL437584` | `None` | mutation=L858R; mutant_text='L858R' | Inhibition of EGFR L858R mutant by time resolved fluorescence assay |

---

## 4. Assay-Context Distribution

| Assay Classification Context | Record Count | Disposition |
| :--- | :--- | :--- |
| Cell-Free / Recombinant Biochemical Kinase Domain (Retained) | 8,888 | **Retained in Curated Set** |
| Cellular Autophosphorylation (Excluded) | 1,490 | *Excluded from Curated Set* |
| Cell-Based / Intact Cell Viability (Excluded) | 170 | *Excluded from Curated Set* |
| Unclear / Ambiguous Context (Excluded) | 17 | *Excluded from Curated Set* |

---

## 5. Missing-Value Statistics in Raw Dataset

| Field Name | Missing Count | Impact |
| :--- | :--- | :--- |
| `standard_value` | 1 | Dropped during Filter 1 |
| `canonical_smiles` | 17 | Dropped during Filter 5 |

---

## 6. Conflict Analysis and Provenance Preservation

- **Compounds with Single Unique Measurement:** 5,158 (80.6%)
- **Compounds with Concordant Multiple Measurements (Spread <= 1.0 log unit):** 820 (12.8%)
- **Compounds with Severe Protocol Conflicts (Spread > 1.0 log unit):** 420 (6.6%)

### Strict Handling Policy for Conflicts
In accordance with Phase 1.3.1 requirements, measurements from differing assay protocols are **not** blindly averaged or medians computed across incompatible assay contexts. Every measurement is preserved individually with its exact `activity_id`, `assay_chembl_id`, `document_chembl_id`, and `conflict_status`. Compounds flagged with `SEVERE_ASSAY_CONFLICT_GT_1LOG` are preserved for auditability but will be partitioned or pruned during QSAR training to avoid corrupting ground-truth SAR.

---

## 7. Exact Benchmark Lead Provenance

### 7.1 Compound 1: CHEMBL289959 (N-phenylquinazolin-4-amine)
- **Canonical SMILES:** `c1ccc(Nc2ncnc3ccccc23)cc1`
- **Selected Benchmark Record:** Activity ID **`996053`**
- **Assay ID:** `CHEMBL680017`
- **Document ID:** `CHEMBL1129291` (Bridges et al., *J. Med. Chem.* 1996, 39, 267; DOI: 10.1021/jm9503613)
- **Assay Description:** *Concentration required to inhibit the phosphorylation of a 14-residue fragment of phospholipase C gamma 1 by Epidermal growth factor receptor enzyme*
- **Standard Measurement:** `344.0 nM` (Relation: `=`; pChEMBL: `6.46`)
- **Format:** Recombinant purified human EGFR kinase catalytic domain; cell-free in vitro radiometric assay

### 7.2 Compound 2: CHEMBL290096 (N-(3-bromophenyl)quinazolin-4-amine)
- **Canonical SMILES:** `Brc1cccc(Nc2ncnc3ccccc23)c1`
- **Investigation of 27 nM vs 28 nM Discrepancy:**
  * **1996 Bridges et al. Record:** Activity ID **`1016106`**, Assay **`CHEMBL680017`**, Document **`CHEMBL1129291`** (Bridges et al., 1996). Standard Value: **`27.0 nM`** (Relation: `=`; pChEMBL: `7.57`). Description: *Concentration required to inhibit the phosphorylation of a 14-residue fragment of phospholipase C gamma 1 by Epidermal growth factor receptor enzyme*.
  * **2008 Gontarewitz et al. Record:** Activity ID **`2111246`**, Assay **`CHEMBL952131`**, Document **`CHEMBL1143196`** (Gontarewitz et al., 2008). Standard Value: **`27.54 nM`** (Relation: `=`; pChEMBL: `7.56`). Description: *Inhibition of EGFR*.
- **Definitive Provenance Attribution:**
  The benchmark lead selection strictly uses the **1996 primary measurement (Activity ID `1016106`, Assay `CHEMBL680017`, `27.0 nM`, pIC50 = 7.57)**, which was generated in the **exact same paper, same laboratory, and identical assay** as `CHEMBL289959` (Activity ID `996053`). The 2008 value (`27.54 nM`, rounded in literature as ~28 nM) is a secondary compilation record and is kept separate.

---

## 8. Preserved Raw Values Guarantee

All raw measurements (`standard_value`, `standard_units`, `standard_relation`) have been strictly preserved in both `egfr_raw_activities.csv` and `egfr_curated_activities.csv`. No raw values have been overwritten with derived `pIC50` values.
