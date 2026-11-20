"""Streamlit dashboard for PHARMAFORGE scientific reasoning and lead optimization."""
from __future__ import annotations

import streamlit as st
import pandas as pd
import json

from pharmaforge import __version__
from pharmaforge.agents.orchestrator import PharmaForgeOrchestrator
from pharmaforge.benchmarks.tdc_loader import get_standard_benchmark_suite
from pharmaforge.core.types import LeadOptimizationTask, OptimizationObjective, ObjectiveDirection
from pharmaforge.memory.skill_library import VerifiedProceduralSkillLibrary
from pharmaforge.verifier.oracle import ScientificTestOracle

st.set_page_config(
    page_title="PHARMAFORGE Dashboard",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Initialize core instances
@st.cache_resource
def get_system_components():
    oracle = ScientificTestOracle()
    library = VerifiedProceduralSkillLibrary(oracle=oracle)
    orchestrator = PharmaForgeOrchestrator(oracle=oracle, skill_library=library)
    return oracle, library, orchestrator

oracle, library, orchestrator = get_system_components()

# --- SIDEBAR ---
with st.sidebar:
    st.title("🧪 PHARMAFORGE")
    st.caption(f"Research Prototype v{__version__}")
    st.markdown("---")
    st.markdown("**Lab Alignment:**")
    st.markdown("Prof. Lin Tan (ASSET Lab, Purdue)")
    st.markdown("**Research Focus:**")
    st.markdown("*Self-Improving AI Agents for Pharmaceutical Applications and Beyond*")
    st.markdown("---")
    st.metric("Verified Skills in VPSL", len(library.skills))
    st.metric("Verifier Backend", "RDKit" if oracle.is_rdkit_backed else "Heuristic Engine")
    st.markdown("---")
    st.info("⚠️ **Research Disclaimer**: All properties are computational approximations. In silico predictions do not establish biological efficacy or clinical safety.")

# --- MAIN TABS ---
tab1, tab2, tab3, tab4 = st.tabs([
    "🚀 Lead Optimization Studio",
    "📚 Verified Skill Library (VPSL)",
    "🔬 Scientific Test-Oracle (STO-IV)",
    "📊 Benchmarks & Architecture",
])

# TAB 1: LEAD OPTIMIZATION STUDIO
with tab1:
    st.header("Closed-Loop Multi-Objective Lead Optimization")
    st.markdown(
        "Execute automated, program-guided lead optimization. The agent plans transformations, "
        "recalls or synthesizes verified procedural skills, executes them in sandbox, and validates candidates via STO-IV."
    )

    benchmark_tasks = get_standard_benchmark_suite()
    task_options = {f"{t.task_id} ({t.target_name} - {t.target_class})": t for t in benchmark_tasks}
    
    col_sel, col_btn = st.columns([3, 1])
    with col_sel:
        selected_task_name = st.selectbox("Select Standard Benchmark Task:", list(task_options.keys()))
        selected_task = task_options[selected_task_name]
    with col_btn:
        st.write("")
        st.write("")
        run_btn = st.button("▶ Run Optimization", type="primary", use_container_width=True)

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Task Specification")
        st.text_input("Target Protein / Endpoint:", selected_task.target_name, disabled=True)
        st.text_input("Initial Lead SMILES:", selected_task.initial_smiles, disabled=True)
    with col2:
        st.subheader("Multi-Objective Constraints")
        obj_df = pd.DataFrame([
            {
                "Property": obj.name,
                "Goal": obj.direction.value,
                "Target Threshold": obj.threshold,
                "Weight": obj.weight,
            }
            for obj in selected_task.objectives
        ])
        st.dataframe(obj_df, use_container_width=True)

    if run_btn:
        with st.spinner("Executing closed-loop optimization through VPSL, Sandbox, and STO-IV..."):
            result = orchestrator.run_optimization(selected_task)

        if result.get("success"):
            st.success(f"Optimization completed in {result['wall_clock_seconds']}s! Discovered {result['total_candidates_generated']} candidates ({result['pareto_front_size']} on Pareto front).")

            if result["newly_learned_skills"]:
                st.info(f"✨ **New Skills Synthesized & Verified**: {', '.join(result['newly_learned_skills'])}")

            st.subheader("Non-Dominated Pareto Frontier")
            pareto_data = result["pareto_front"]
            if pareto_data:
                rows = []
                for p in pareto_data:
                    rows.append({
                        "SMILES": p["smiles"],
                        "Affinity Score": p.get("affinity_score"),
                        "QED": p.get("qed"),
                        "LogP": p.get("logp"),
                        "SAScore": p.get("sascore"),
                        "Derived From": p.get("parent_smiles"),
                        "Skill Applied": p.get("transformation_applied"),
                    })
                df_pareto = pd.DataFrame(rows)
                st.dataframe(df_pareto, use_container_width=True)

            with st.expander("View Full Trajectory & Deliberation Steps"):
                st.json(result["trajectory_summary"])

# TAB 2: VERIFIED PROCEDURAL SKILL LIBRARY (VPSL)
with tab2:
    st.header("Verified Procedural Skill Library (VPSL)")
    st.markdown(
        "Skills in VPSL are **parameterized, executable Python procedures** that passed automated unit test harnesses. "
        "Unlike text reflections, these procedural assets persist across research episodes."
    )

    query = st.text_input("🔍 Semantic Skill Retrieval Query:", "reduce logP and improve aqueous solubility")
    if query:
        retrieved = library.retrieve_skills(query, top_k=5, threshold=0.05)
        st.caption(f"Found {len(retrieved)} relevant skills matching query.")
        for skill, score in retrieved:
            with st.container(border=True):
                st.markdown(f"### `{skill.name}` (Relevance: {score:.2f})")
                st.markdown(f"**Description:** {skill.description}")
                col_m1, col_m2, col_m3 = st.columns(3)
                col_m1.metric("Reliability Score", f"{skill.reliability_score * 100:.1f}%")
                col_m2.metric("Success Executions", skill.success_count)
                col_m3.metric("Failures", skill.failure_count)
                st.code(skill.code, language="python")
                with st.expander("Unit Test Harness"):
                    st.code(skill.unit_tests, language="python")

# TAB 3: SCIENTIFIC TEST-ORACLE (STO-IV)
with tab3:
    st.header("Scientific Test-Oracle (STO-IV)")
    st.markdown(
        "Deterministic verification gate that intercepts chemical hallucinations, explicit valence violations, "
        "and toxicophore alerts (PAINS, Brenk) prior to admitting candidates or skills."
    )

    test_smiles = st.text_input("Enter Candidate Molecule (SMILES):", "CC(=O)Oc1ccccc1C(=O)O")
    if st.button("Validate Molecule"):
        report = oracle.verify_smiles(test_smiles)
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Chemical Validity", "PASSED ✅" if report.is_valid else "FAILED ❌")
        c2.metric("RDKit Sanitization", "PASSED" if report.passed_sanitization else "FAILED")
        c3.metric("PAINS Filter", "PASSED" if report.passed_pains else "ALERT ⚠️")
        c4.metric("Brenk Filter", "PASSED" if report.passed_brenk else "ALERT ⚠️")

        if report.error_reasons:
            st.error(f"Rejection Reasons: {', '.join(report.error_reasons)}")

        st.subheader("Physicochemical Property Profile")
        st.json(report.property_scores)

# TAB 4: BENCHMARKS & ARCHITECTURE
with tab4:
    st.header("System Architecture & Research Foundations")
    st.markdown("""
    ### Intellectual Synergy with Professor Lin Tan's ASSET Lab
    1. **Automated Program Repair (APR):** Skills with runtime exceptions are repaired via localized AST fault diagnosis.
    2. **Test-Harness Selection (TENET):** Procedural skills are gated by formal unit test suites rather than subjective LLM self-evaluation.
    3. **Unified Agent Capabilities (USEagent):** Seamless state switching between planning, tool synthesis, sandboxed execution, and verification.
    
    ### Formal Verification Pipeline
    - **Step 1:** Scientific Discovery Planner breaks multi-objective goal into prioritized deficiencies.
    - **Step 2:** VPSL semantic retrieval matches procedural skills (or triggers dynamic skill synthesis).
    - **Step 3:** Isolated Docker Execution Sandbox runs Python functions under strict time/memory quotas.
    - **Step 4:** Independent Verifier (STO-IV) deterministically checks validity, PAINS, and property constraints.
    - **Step 5:** Pareto frontier update and skill reliability recording.
    """)
