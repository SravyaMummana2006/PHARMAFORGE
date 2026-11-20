# PHARMAFORGE: A Self-Improving Multi-Agent Framework for Scientific Reasoning and Pharmaceutical Discovery

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Research Alignment](https://img.shields.io/badge/Research-Purdue%20ASSET%20Lab-orange.svg)](https://www.cs.purdue.edu/homes/lintan/)

Aligned with the research initiative:
**"Self-Improving AI Agents for Pharmaceutical Applications and Beyond"**  
Associated with **Professor Lin Tan**, Purdue University (ASSET Lab).

---

## 1. Executive Overview

Current scientific LLM agents (such as ChemCrow and Coscientist) suffer from high hallucination rates, lack independent verification, and exhibit **episodic amnesia**—every discovery task starts from scratch, or errors accumulate via unverified prompt reflections.

Drawing upon foundational principles in **software dependability, automated program repair (APR), and test-driven synthesis** developed by Professor Lin Tan's research group (*USEagent*, *TENET*, *CoRe*, *CURE*), **PHARMAFORGE** conceptualizes chemical workflows as **executable, testable, and self-repairing computational programs**.

```mermaid
flowchart TD
    Lead["Candidate Lead Molecule + Multi-Objective Goals"] --> Planner
    
    subgraph "PHARMAFORGE Core"
        Planner["Scientific Discovery Planner"]
        VPSL["Verified Procedural Skill Library (VPSL)"]
        Synthesizer["Tool & Skill Synthesizer"]
        Sandbox["Isolated Execution Sandbox"]
        STO_IV["Scientific Test-Oracle (STO-IV)"]
        Repair["Execution-Trace Self-Repair Engine (APR)"]
    end
    
    Planner <--> VPSL
    Planner --> Synthesizer
    Synthesizer --> Sandbox
    Sandbox --> STO_IV
    STO_IV -- "Fail / Exception Trace" --> Repair
    Repair --> Sandbox
    STO_IV -- "Passed Unit Tests" --> VPSL
    STO_IV --> Pareto["Validated Non-Dominated Pareto Frontier"]
```

### Key Technical Contributions:
1. **Verified Procedural Skill Library (VPSL):** Dynamic, vector-indexed memory of executable, parameterized Python procedures gated by automated test oracles.
2. **Scientific Test-Oracle (STO-IV):** Deterministic, independent verification engine validating RDKit sanitization, explicit atomic valence, and PAINS/Brenk toxicophore filters.
3. **Execution-Trace Self-Repair Engine (ET-SRE):** Adapting automated program repair to scientific skills, diagnosing runtime exceptions, and iteratively synthesizing code patches.
4. **Pareto Hypervolume Optimization:** Multi-objective balancing across conflicting pharmacological criteria (affinity, QED, logP, SAScore, hERG).

---

## 2. Directory Structure

```text
PHARMAFORGE/
├── .gitignore
├── README.md
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
├── pyproject.toml
├── config/
│   └── default_config.yaml           # System settings, timeouts, and thresholds
├── src/
│   └── pharmaforge/
│       ├── core/
│       │   ├── types.py              # Pydantic data schemas
│       │   └── exceptions.py         # Domain-specific exceptions
│       ├── verifier/
│       │   ├── oracle.py             # Scientific Test-Oracle (STO-IV)
│       │   └── filters.py            # PAINS and Brenk structural alert catalogs
│       ├── sandbox/
│       │   └── executor.py           # Subprocess execution sandbox with AST checks
│       ├── memory/
│       │   ├── episodic.py           # Short-term trajectory buffers
│       │   └── skill_library.py      # Verified Procedural Skill Library (VPSL)
│       ├── agents/
│       │   ├── base.py               # Base agent interface
│       │   ├── planner.py            # Scientific Discovery Planner
│       │   ├── synthesizer.py        # Tool & Skill Synthesizer
│       │   ├── repair.py             # Execution-Trace Self-Repair Engine
│       │   └── orchestrator.py       # Closed-loop multi-agent coordinator
│       └── benchmarks/
│           ├── metrics.py            # Pareto Hypervolume, chemical validity metrics
│           └── tdc_loader.py         # Therapeutics Data Commons benchmark suites
├── tests/                            # Comprehensive Pytest suite
│   ├── test_verifier.py
│   ├── test_sandbox.py
│   ├── test_skill_library.py
│   ├── test_repair_engine.py
│   ├── test_metrics.py
│   └── test_orchestrator.py
├── api/
│   └── main.py                       # FastAPI REST API server
└── ui/
    └── app.py                        # Streamlit scientific dashboard
```

---

## 3. Quickstart & Installation

### Option A: Local Python Environment
```bash
# Clone the repository and navigate to root
cd PHARMAFORGE

# Install dependencies in editable mode
pip install -e ".[science,agent,ui,dev]"

# Run comprehensive test suite
pytest
```

### Option B: Docker Compose Deployment
```bash
# Build and launch MLflow, FastAPI backend, and Streamlit dashboard
docker-compose up --build
```
- **Streamlit Scientific Dashboard:** `http://localhost:8501`
- **FastAPI Documentation:** `http://localhost:8000/docs`
- **MLflow Experiment Tracking:** `http://localhost:5000`

---

## 4. Benchmark Tasks (Therapeutics Data Commons)

| Task ID | Target / Focus | Target Class | Key Objectives |
| :--- | :--- | :--- | :--- |
| `TDC-EGFR-01` | EGFR Kinase | Kinase | Maximize Affinity ($\ge 7.0$), Maximize QED, Minimize logP, Minimize SAScore |
| `TDC-DRD2-02` | Dopamine D2 | GPCR | Maximize Affinity ($\ge 7.2$), Eliminate hERG liability, Maximize QED |
| `TDC-AqSol-03` | Aqueous Solubility | ADMET | Minimize logP ($\le 2.5$), Maximize QED, Minimize SAScore |
| `TDC-hERG-04` | Cardiotoxicity Mitigation | Safety | De-risk basic amine center, maintain target affinity |

---

## 5. Research Integrity & Disclaimer

> **IMPORTANT:**
> All molecular properties, binding affinities, and ADMET profiles evaluated by PHARMAFORGE are **in silico computational approximations**. Computational predictions do not constitute medical, biological, or pharmacological proof of efficacy or safety. This platform is an experimental scientific reasoning prototype for computational research and methodological advancement.
