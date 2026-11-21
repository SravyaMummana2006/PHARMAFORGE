"""
Script to verify that PHARMAFORGE Docker Compose deployment is healthy and operational.
Tests:
  1. FastAPI health endpoint (/health) including QSAR model and RDKit readiness.
  2. FastAPI documentation endpoint (/docs).
  3. FastAPI benchmark tasks endpoint (/api/v1/benchmarks).
  4. Deterministic molecular verification endpoint (/api/v1/verify).
  5. Multi-objective lead optimization endpoint (/api/v1/optimize).
  6. Streamlit health endpoint (/_stcore/health).
  7. Streamlit main application page (/).
"""

import sys
import time
import json
import urllib.request
import urllib.error

FASTAPI_URL = "http://localhost:8000"
STREAMLIT_URL = "http://localhost:8501"

def wait_for_endpoint(url, timeout=60, step=2):
    start = time.time()
    while time.time() - start < timeout:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PharmaForgeVerifier/1.0"})
            with urllib.request.urlopen(req, timeout=5) as resp:
                if resp.status == 200:
                    return True
        except Exception:
            pass
        time.sleep(step)
    return False

def http_get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "PharmaForgeVerifier/1.0"})
    with urllib.request.urlopen(req, timeout=10) as resp:
        return resp.status, resp.read().decode("utf-8")

def http_post(url, data_dict):
    data_bytes = json.dumps(data_dict).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data_bytes,
        headers={"Content-Type": "application/json", "User-Agent": "PharmaForgeVerifier/1.0"},
        method="POST"
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.status, json.loads(resp.read().decode("utf-8"))

def main():
    print("=" * 70)
    print("PHARMAFORGE DOCKER COMPOSE DEPLOYMENT VERIFICATION")
    print("=" * 70)

    # 1. Wait for FastAPI
    print("\n[1/7] Probing FastAPI backend health endpoint...")
    if not wait_for_endpoint(f"{FASTAPI_URL}/health", timeout=30):
        print(f"FAILED: FastAPI backend not reachable at {FASTAPI_URL}/health after 30s.")
        sys.exit(1)
    status, body = http_get(f"{FASTAPI_URL}/health")
    health = json.loads(body)
    print(f"  HTTP Status: {status}")
    print(f"  Service: {health.get('service')} | Version: {health.get('version')}")
    print(f"  RDKit Available: {health.get('rdkit_available')}")
    print(f"  QSAR Model Loaded: {health.get('qsar_model_loaded')}")
    assert health.get("status") == "healthy", f"FastAPI health check failed: {health}"
    assert health.get("rdkit_available") is True, "RDKit is not available in FastAPI container"
    assert health.get("qsar_model_loaded") is True, "QSAR model failed to load in FastAPI container"
    print("  --> PASS: FastAPI is healthy, RDKit ready, QSAR model loaded.")

    # 2. Check OpenAPI Docs
    print("\n[2/7] Probing FastAPI OpenAPI documentation (/docs)...")
    status, body = http_get(f"{FASTAPI_URL}/docs")
    print(f"  HTTP Status: {status}")
    assert status == 200, f"FastAPI /docs returned non-200: {status}"
    assert "swagger-ui" in body.lower() or "pharmaforge" in body.lower(), "Docs page content unexpected"
    print("  --> PASS: Interactive OpenAPI documentation accessible.")

    # 3. Check Benchmark Endpoints
    print("\n[3/7] Probing FastAPI benchmark tasks (/api/v1/benchmarks)...")
    status, body = http_get(f"{FASTAPI_URL}/api/v1/benchmarks")
    benchmarks = json.loads(body)
    print(f"  HTTP Status: {status}")
    tasks = benchmarks if isinstance(benchmarks, list) else benchmarks.get("tasks", [])
    print(f"  Available Benchmark Tasks: {len(tasks)}")
    task_ids = [t.get("task_id") for t in tasks]
    print(f"  Task IDs: {task_ids}")
    assert any("CHEMBL289959" in tid for tid in task_ids), "CHEMBL289959 missing from benchmarks"
    assert any("CHEMBL290096" in tid for tid in task_ids), "CHEMBL290096 missing from benchmarks"
    print("  --> PASS: Standard benchmark suite correctly loaded.")

    # 4. Check Molecular Verification Endpoint (STO-IV)
    print("\n[4/7] Testing STO-IV verification (/api/v1/verify)...")
    verify_payload = {"smiles": "c1ccc(Nc2ncnc3ccccc23)cc1"}
    status, verify_resp = http_post(f"{FASTAPI_URL}/api/v1/verify", verify_payload)
    print(f"  HTTP Status: {status}")
    is_valid = verify_resp.get("is_valid", verify_resp.get("valid"))
    prop_scores = verify_resp.get("property_scores", verify_resp.get("properties", {}))
    print(f"  Chemical Validity (is_valid): {is_valid}")
    print(f"  Property Scores: {prop_scores}")
    assert is_valid is True, f"Lead molecule failed STO-IV verification: {verify_resp}"
    assert "qed" in prop_scores and "logp" in prop_scores, f"Missing properties: {prop_scores}"
    print("  --> PASS: STO-IV molecular verification operational.")

    # 5. Check Lead Optimization Endpoint
    print("\n[5/7] Testing 1-iteration lead optimization (/api/v1/optimize)...")
    opt_payload = {"task_id": "EGFR-BENCH-CHEMBL289959", "max_iterations": 1}
    status, opt_resp = http_post(f"{FASTAPI_URL}/api/v1/optimize", opt_payload)
    print(f"  HTTP Status: {status}")
    print(f"  Task ID: {opt_resp.get('task_id')}")
    pareto_front = opt_resp.get("pareto_front", [])
    print(f"  Pareto Front Size: {len(pareto_front)}")
    for i, cand in enumerate(pareto_front[:3]):
        pIC50 = cand.get('affinity_score') or cand.get('predicted_pic50')
        print(f"    Candidate {i+1}: {cand.get('smiles')} | Pred pIC50={pIC50} | QED={cand.get('qed')}")
    assert len(pareto_front) > 0, "Optimization returned empty Pareto front"
    print("  --> PASS: Closed-loop multi-objective optimization operational.")

    # 6. Check Streamlit Health
    print("\n[6/7] Probing Streamlit health endpoint (/_stcore/health)...")
    if not wait_for_endpoint(f"{STREAMLIT_URL}/_stcore/health", timeout=30):
        print(f"FAILED: Streamlit health endpoint not reachable at {STREAMLIT_URL}/_stcore/health after 30s.")
        sys.exit(1)
    status, body = http_get(f"{STREAMLIT_URL}/_stcore/health")
    print(f"  HTTP Status: {status} | Body: {body.strip()}")
    assert status == 200 and body.strip() == "ok", f"Streamlit health check unexpected: {status}, {body}"
    print("  --> PASS: Streamlit server is healthy.")

    # 7. Check Streamlit Main Page
    print("\n[7/7] Probing Streamlit main dashboard page (/)...")
    status, body = http_get(f"{STREAMLIT_URL}/")
    print(f"  HTTP Status: {status}")
    assert status == 200, f"Streamlit returned non-200 status: {status}"
    assert "PHARMAFORGE" in body or "streamlit" in body.lower(), "Streamlit page content unexpected"
    print("  --> PASS: Streamlit scientific dashboard UI accessible.")

    print("\n" + "=" * 70)
    print("ALL 7 DOCKER DEPLOYMENT VERIFICATION CHECKS PASSED SUCCESSFULLY!")
    print("=" * 70)

if __name__ == "__main__":
    main()
