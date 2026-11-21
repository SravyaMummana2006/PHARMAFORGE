# PHARMAFORGE: Live Demonstration & Research Evaluation Guide

**Project:** PHARMAFORGE (*A Self-Improving Multi-Agent Framework for Scientific Reasoning and Pharmaceutical Discovery*)  
**Research Alignment:** Professor Lin Tan (ASSET Lab, Purdue University)  
**Target:** Undergraduate Research Evaluation, Committee Review & Live Demonstration  

---

## 1. System Requirements & Prerequisites

- **Operating System:** Windows 10/11 (or Linux / macOS)
- **Python Runtime:** Python 3.10+ (Recommended: Python 3.11 for UI/API services; Python 3.11 or 3.13 with RDKit and scikit-learn for scientific computation)
- **Key Scientific Packages:** `rdkit`, `scikit-learn`, `joblib`, `numpy`, `scipy`, `pandas`
- **Key Web & API Packages:** `streamlit`, `fastapi`, `uvicorn`, `pydantic`

> **Note on Environment Architecture:**  
> On Windows, if your scientific stack (RDKit) and web stack (Streamlit) are installed in separate Python interpreters, PHARMAFORGE automatically bridges them via its subprocess execution runner. You can also explicitly specify the RDKit Python interpreter path by setting:
> ```powershell
> $env:PHARMAFORGE_PYTHON = "C:\Path\To\Python\python.exe"
> ```

---

## 2. Setup & Installation Commands

### Option A: Standard Virtual Environment (Recommended)
```powershell
# 1. Clone and navigate to the repository
cd PHARMAFORGE

# 2. Create and activate a virtual environment
python -m venv venv
.\venv\Scripts\Activate.ps1

# 3. Install dependencies in editable mode
pip install -r requirements.txt
pip install -e .
```

---

## 3. Launching Services

### 3.1 Launching the FastAPI Backend Server
The FastAPI backend exposes REST endpoints for deterministic molecular verification (STO-IV), skill retrieval (VPSL), and multi-objective optimization.

```powershell
python -m uvicorn api.main:app --host 127.0.0.1 --port 8000
```
- **Service Root:** `http://127.0.0.1:8000/`
- **Interactive OpenAPI Documentation:** `http://127.0.0.1:8000/docs`
- **Health Check:** `http://127.0.0.1:8000/health`

### 3.2 Launching the Streamlit Scientific Dashboard
The interactive dashboard provides instant inspection of pre-computed benchmark ablation studies, live single-step lead optimization, semantic skill retrieval, and STO-IV molecular validation.

```powershell
python -m streamlit run ui/app.py --server.port 8501
```
- **Local URL:** `http://localhost:8501`

---

## 4. Suggested 5–7 Minute Live Demonstration Sequence

### Step 1: Scientific Test-Oracle (STO-IV) [1 Minute]
1. Navigate to **Tab 4: Scientific Test-Oracle (STO-IV)**.
2. Enter the clinical lead molecule `Brc1cccc(Nc2ncnc3ccccc23)c1` (`CHEMBL290096`).
3. Click **"Validate Molecule"**.
4. **Key Talking Point:** Show that STO-IV deterministically executes RDKit sanitization, explicit atomic valence checks, and PAINS/Brenk toxicophore filters in an isolated sandbox, returning property estimates ($\text{QED} = 0.773, \, \text{LogP} = 4.136, \, \text{SAScore} = 1.802$).

### Step 2: The Core Ablation Study & Research Findings [2–3 Minutes]
1. Navigate to **Tab 2: Research Contribution & Ablation Study**.
2. Select target lead **`CHEMBL289959`** ($R = \text{H}$, unsubstituted 4-anilinoquinazoline).
3. Walk through the 4-way comparative table:
   - **Configuration A (No Memory, No Repair):** Stateless baseline. Clean skills synthesized each iteration; zero persistent reuse.
   - **Configuration B (Unverified Memory):** *Demonstrates Memory Poisoning.* Unverified buggy skills containing missing imports (`NameError`) are admitted directly into VPSL (`Admitted = 1`). Across iterations 2–5, the planner retrieves this faulty skill (`Reuse = 4`), triggering 5 sandbox crashes and 0 valid candidates (validity rate is **N/A**, not 0%).
   - **Configuration C (Verified Memory, No Repair):** *Demonstrates Pure Memory Without Recovery.* STO-IV unit-test verification rejects buggy skills (`Rejected = 5`). Memory is protected from poisoning, but the agent is left with no working generative tools (0 valid candidates).
   - **Configuration D (Full PHARMAFORGE):** *Demonstrates Self-Repair & Re-Verification.* STO-IV intercepts the failure; ET-SRE diagnoses the runtime traceback, patches the missing import, re-verifies in the sandbox, admits the repaired skill, and reuses it across subsequent iterations (`Reuse = 4`), discovering 50/50 valid candidates (100% validity) and expanding the Pareto frontier to 7 non-dominated solutions.
4. **Key Talking Point:** Point out the **Monte Carlo Hypervolume** ($N = 20,000, \, \text{SE} \approx 0.0026$) values ($0.1578$ for Configs A & D vs $0.1070$ for Configs B & C), explaining why Config A and D achieve identical hypervolume (repaired skills execute identical chemical transforms as clean skills).

### Step 3: Lead Optimization Studio & Constraint Non-Compliance [1–2 Minutes]
1. Navigate to **Tab 1: Lead Optimization Studio**.
2. Select task **`EGFR-OPT-CHEMBL289959`** or **`EGFR-OPT-CHEMBL290096`**.
3. View the **Non-Dominated Pareto Frontier** table.
4. Highlight the **"LogP <= 3.5"** constraint column:
   - Point out that pyrido-pyrimidine isosteres (`Brc1cccc(Nc2ncnc3cnccc23)c1`, $\text{LogP} = 3.531$) significantly lower lipophilicity ($\Delta = -0.605\log$ relative to starting lead $4.136$).
   - **Epistemic Honesty Point:** Emphasize that because $3.531 > 3.500$, this is a near-threshold candidate that strictly **fails** the feasibility constraint. Explain that Pareto selection is based on mathematical non-domination across continuous objectives, not constraint satisfaction.

### Step 4: Verified Procedural Skill Library (VPSL) [1 Minute]
1. Navigate to **Tab 3: Verified Procedural Skill Library (VPSL)**.
2. Search for `"aromatic affinity optimization"` or `"polar solubilizing"`.
3. Show the retrieved Python procedure, its RDKit Reaction SMARTS, and its unit-test harness (`run_tests()`).

---

## 5. Automated CLI Demonstration Script

For a complete programmatic demonstration without opening a browser, execute the verified 7-step test script:

```powershell
python scripts/verify_demo_pipeline.py
```

### Expected Output Summary:
```text
======================================================================
PHARMAFORGE: LIVE REPRODUCIBLE END-TO-END DEMO VERIFICATION
======================================================================

[STEP 1] Benchmark Lead Molecule: CHEMBL289959
  SMILES: c1ccc(Nc2ncnc3ccccc23)cc1
  Description: 4-anilinoquinazoline with unsubstituted phenyl (Rewcastle 1996)
  Experimental Activity: IC50 = 344 nM, pChEMBL = 6.463
  Lead STO-IV Validity: True
  Lead Baseline Properties: QED=0.721, LogP=3.373, SAScore=1.588
  Lead Baseline Predicted pIC50: 6.419 (IC50 = 380.85 nM)
  Applicability Domain: T_max = 0.686

[STEP 2 & 3] Executing Procedural Skill on Lead Molecule
  Execution Success: True
  Candidates Generated: 8

[STEP 4] Scientific Verification & Multi-Objective Evaluation (STO-IV + QSAR)
  Verified Valid Candidates: 8 / 8
    Brc1cccc(Nc2ncnc3ccccc23)c1 | Pred pIC50=7.683 | QED=0.773 | LogP=4.136 | SAScore=1.802 <-- AUTONOMOUS REDISCOVERY OF CHEMBL290096!

[STEP 5] Testing Fault Interception & ET-SRE Automated Program Repair
  Testing registration in VPSL (should fail unit-test execution due to missing import)...
  Registration Blocked: True
  Interception Reason: Skill 'rdkit_qed_bioisostere_optimization' failed unit test verification...
  Triggering ET-SRE Traceback Diagnosis and Patching...
  Repair Success: True
    Attempt 1: Applied patch for 'undefined_name'
    Success: Unit tests passed on attempt 1.
  Patched Code Preview (Line 1):
    from rdkit import Chem

[STEP 6] Re-Verification and Subsequent Skill Reuse in VPSL
  Repaired Skill Admitted to VPSL: True
  VPSL Retrieval Match: rdkit_qed_bioisostere_optimization (Score: 0.21)
  Reused Skill Executed Successfully: True
  Candidates Produced on Reuse: 10

[STEP 7] Multi-Objective Pareto Comparison & Constraint Feasibility
  Lead CHEMBL289959: pIC50=6.419, QED=0.721, LogP=3.373
  Top Brominated Candidate: SMILES=Brc1ccc(Nc2ncnc3ccccc23)cc1
    pIC50=6.85 (Delta: +0.431 log units)
    QED=0.773, LogP=4.136, SAScore=1.734
    LogP Constraint Feasibility: FAIL (4.136 > 3.5)

======================================================================
ALL 7 DEMONSTRATION WORKFLOW STEPS SUCCESSFULLY VERIFIED PROGRAMMATICALLY
======================================================================
```

---

## 6. Common Errors & Troubleshooting

| Error | Root Cause | Verified Solution |
| :--- | :--- | :--- |
| `Address already in use` (Port 8501 or 8000) | Another service or terminal is occupying the default port. | Launch on a different port:<br>`python -m streamlit run ui/app.py --server.port 8502`<br>`python -m uvicorn api.main:app --port 8001` |
| `ModuleNotFoundError: No module named 'rdkit'` | Active Python environment lacks RDKit. | Ensure RDKit is installed in your active environment (`pip install rdkit`) or set `$env:PHARMAFORGE_PYTHON` pointing to your RDKit interpreter. |
| `ModuleNotFoundError: No module named 'pharmaforge'` | Package not installed in editable mode. | Run `pip install -e .` from repository root, or verify that `src` is in `sys.path` (handled automatically by `ui/app.py` and `api/main.py`). |
| `LayoutsMixin.container() got unexpected keyword argument 'border'` | Using Streamlit < 1.30.0. | Dashboard code has been updated to use standard `st.container()`, compatible with Streamlit 1.28.0+. |
