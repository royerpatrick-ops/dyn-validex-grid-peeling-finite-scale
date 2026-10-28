"""Read-only provenance audit of the completed v0.5 replication.

Only replays the saved PCG64 seed to verify the existing ordered list. This
does not create a seed or a new collection, run project code, or analyze the
scientific outcomes. All outputs are confined to the audit output path.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", type=Path, required=True)
    parser.add_argument("--execute", type=Path, required=True)
    parser.add_argument("--prior-prepare", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("provenance_results_v05.json"))
    args = parser.parse_args()
    P, E = args.prepare.resolve(), args.execute.resolve()
    prefixes = {"_prepare_replication_v0_4_20261027": P, "_execute_replication_v0_5_20261027": E}

    def remap(text):
        normalized = str(text).replace("\\", "/")
        parts = normalized.split("/")
        for prefix, root in prefixes.items():
            if prefix in parts:
                start = parts.index(prefix)
                return root.joinpath(*parts[start + 1:])
        return None

    def verify_records(records):
        ok, bad, unavailable = [], [], []
        for record in records:
            path = remap(record["path"])
            if path is None or not path.is_file():
                unavailable.append(record["path"])
                continue
            valid = digest(path) == record["sha256"] and ("bytes" not in record or path.stat().st_size == record["bytes"])
            (ok if valid else bad).append(record["path"])
        return {"matched_count": len(ok), "mismatches": bad, "unavailable": unavailable}

    frozen_path = P / "frozen_manifest_v04.json"
    frozen = read(frozen_path)
    frozen_sha = digest(frozen_path)
    state = read(P / "launcher_state_v04.json")
    translations = read(P / "translations_v04.json")
    plan = read(P / "plan/PLAN_REPLICATION_INDEPENDANTE_v0_2.json")
    execution = read(E / "execution_manifest_v05.json")
    receipt = read(E / "control/control_creation_receipt_v05.json")
    complete = read(P / "campaign_v04/campaign_complete_manifest_v04.json")
    costs = read(P / "campaign_v04/campaign_costs_v04.json")
    budget = read(P / "pilot/cost_budget_v04.json")
    delivery = read(E / "delivery_manifest_v05.json")
    analysis_manifest = read(P / "analysis_v04/analysis_manifest_v04.json")
    frozen_check = verify_records(frozen["frozen_files"])
    delivery_check = verify_records(delivery["records"])
    campaign_check = verify_records(complete["outputs"])
    execution_check = verify_records(execution["artifacts"].values())
    analysis_check = verify_records([analysis_manifest["result"], analysis_manifest["cell_contrasts"]])

    seed_input_path = E / "control/seed_input_v05.txt"
    stored_seed_path = P / "master_seed_v04.txt"
    seed = seed_input_path.read_text(encoding="ascii").strip()
    stored_seed = stored_seed_path.read_text(encoding="ascii").strip()
    seed_format_ok = re.fullmatch(r"[0-9a-f]{64}", seed) is not None
    blacklist_groups = plan["randomization"]["blacklist"]
    blacklist = {tuple(map(int, pair)) for group in blacklist_groups.values() for pair in group}
    old16 = {tuple(map(int, pair)) for pair in blacklist_groups["P06_P19"]}
    saved_rows = translations["translations"]
    saved_pairs = [(int(row["tx"]), int(row["ty"])) for row in saved_rows]
    # Verification replay of the saved stream only, never a new seed/list.
    generator = np.random.Generator(np.random.PCG64(int(seed, 16)))
    replayed, seen, calls = [], set(), 0
    while len(replayed) < 64:
        candidate = generator.integers(0, 2**20, size=2, dtype=np.uint64)
        calls += 1
        pair = tuple(map(int, candidate))
        if pair in blacklist or pair in seen:
            continue
        replayed.append(pair)
        seen.add(pair)
    auth_checks = {}
    for action in ("DRAW", "RUN"):
        expected = {"authorize_action": action, "frozen_manifest_sha256": frozen_sha, "acknowledge_external_state_change": True}
        auth_checks[action] = read(E / f"control/authorization_{action.lower()}_v05.json") == expected

    quickhull = []
    for index in (0, 8, 16, 24, 32, 40, 48, 56):
        q = read(P / f"campaign_v04/rows/row_{index:02d}_quickhull_check.json")
        quickhull.append({
            "row": index, "layers_compared": q["layers_compared"],
            "status_ok": q["status"] == "PASS_COMPLETE_LAYER_BY_LAYER_SECOND_ENGINE",
            "exact_checks": all(q[k] is True for k in ("exact_layer_indices", "exact_remaining_and_vertex_counts", "exact_vertex_sets", "exact_target_brackets", "cumulative_layer_record_sha256_equal")),
            "digests_equal": q["cumulative_layer_record_sha256_primary"] == q["cumulative_layer_record_sha256_secondary"],
            "minimum_containment_cross": q["minimum_exact_adjacent_containment_cross"],
        })
    row_costs = costs["rows"]
    main_wall = sum(row["main_wall_seconds"] for row in row_costs)
    second_wall = sum(row["second_wall_seconds"] or 0 for row in row_costs)
    main_cpu = sum(row["main_cpu_seconds"] for row in row_costs)
    second_cpu = sum(row["second_cpu_seconds"] or 0 for row in row_costs)
    campaign_bytes = sum(path.stat().st_size for path in (P / "campaign_v04").rglob("*") if path.is_file())
    analysis_bytes = sum(path.stat().st_size for path in (P / "analysis_v04").rglob("*") if path.is_file())
    upper_time = budget["full_project_including_already_executed_reference_qualification"]["high_wall_seconds_factor_1_5"]
    upper_storage = budget["full_project_including_already_executed_reference_qualification"]["high_storage_bytes_factor_1_25"]
    prior_checks = None
    if args.prior_prepare:
        prior = args.prior_prepare.resolve()
        paths = [record["path"].split("/", 1)[1] for record in frozen["frozen_files"] if record["path"].startswith("_prepare_replication_v0_4_20261027/")]
        prior_checks = {
            "frozen_manifest_byte_identical": (prior / "frozen_manifest_v04.json").read_bytes() == frozen_path.read_bytes(),
            "frozen_files_byte_identical": all((prior / relative).read_bytes() == (P / relative).read_bytes() for relative in paths),
            "prior_state_hash_matches_pre_draw_receipt": digest(prior / "launcher_state_v04.json") == receipt["pre_draw_state_sha256"],
            "files_compared": len(paths),
        }
    checks = {
        "frozen_digest_matches_preapproved": frozen_sha == "023b69c57501ed27ccba892b0e2c2afee88fa1621f711d6c057e3dc0a58f2518",
        "state_complete": state["last_completed_action"] == "ANALYZE" and state["next_action"] == "COMPLETE",
        "state_frozen_manifest_anchor": state["frozen_manifest_sha256"] == frozen_sha,
        "state_translation_anchor": state["translation_manifest_sha256"] == digest(P / "translations_v04.json"),
        "state_seed_anchor": state["master_seed_sha256"] == digest(stored_seed_path),
        "state_completion_anchor": state["campaign_complete_manifest_sha256"] == digest(P / "campaign_v04/campaign_complete_manifest_v04.json"),
        "seed_input_and_stored_same_value": seed == stored_seed == execution["randomization"]["master_seed_hex"],
        "seed_format_ok": seed_format_ok,
        "seed_receipt_hash_match": receipt["seed_sha256"] == digest(seed_input_path),
        "auth_payloads_match": all(auth_checks.values()),
        "blacklist_48_unique": len(blacklist) == 48,
        "64_unique_ordered_rows": len(saved_pairs) == len(set(saved_pairs)) == 64 and [row["row_index"] for row in saved_rows] == list(range(64)),
        "no_historical_48_overlap": not (set(saved_pairs) & blacklist),
        "no_old16_overlap": not (set(saved_pairs) & old16),
        "prng_replay_matches_exact_order": replayed == saved_pairs,
        "prng_candidate_count_matches": calls == translations["candidate_calls"] == execution["randomization"]["candidate_calls"],
        "execution_lists_same_rows": execution["randomization"]["translations_in_frozen_order"] == saved_rows,
        "execution_frozen_files_unchanged": execution["frozen_files"] == frozen["frozen_files"],
        "64_main_artifacts_present": len(list((P / "campaign_v04/rows").glob("row_*_main.npz"))) == 64,
        "eight_quickhull_checks_pass": all(q["status_ok"] and q["exact_checks"] and q["digests_equal"] and q["minimum_containment_cross"] >= 0 for q in quickhull),
        "costs_have_64_rows_8_second": len(row_costs) == 64 and sum(row["second_wall_seconds"] is not None for row in row_costs) == 8,
        "run_wall_below_project_upper_forecast": costs["elapsed_wall_seconds"] <= upper_time,
        "campaign_plus_analysis_storage_below_upper_forecast": campaign_bytes + analysis_bytes <= upper_storage,
        "all_available_records_hash_match": not any(section["mismatches"] for section in (frozen_check, delivery_check, campaign_check, execution_check, analysis_check)),
    }
    if prior_checks:
        checks["preapproved_static_package_identical"] = all(v for k, v in prior_checks.items() if k != "files_compared")
    result = {
        "status": "PROVENANCE_AND_SAVED_STREAM_REPLAY_CONSISTENT" if all(checks.values()) else "PROVENANCE_MISMATCH_FOUND",
        "checks": checks,
        "frozen_records": frozen_check,
        "delivery_records": delivery_check,
        "campaign_records": campaign_check,
        "execution_artifact_records": execution_check,
        "analysis_artifacts": analysis_check,
        "prior_package_comparison": prior_checks,
        "saved_stream_replay": {"candidate_calls": calls, "historical_overlap": len(set(saved_pairs) & blacklist), "old16_overlap": len(set(saved_pairs) & old16), "numpy_audit_version": np.__version__, "new_draw_or_collection": False},
        "quickhull_reports": quickhull,
        "costs": {"run_wall_seconds": costs["elapsed_wall_seconds"], "main_wall_sum": main_wall, "second_wall_sum": second_wall, "run_residual_wall": costs["elapsed_wall_seconds"] - main_wall - second_wall, "main_cpu_sum": main_cpu, "second_cpu_sum": second_cpu, "campaign_bytes": campaign_bytes, "analysis_bytes": analysis_bytes, "upper_project_wall_forecast": upper_time, "upper_storage_forecast": upper_storage, "production_peak_memory": None},
        "prospective_documentary_times": {"freeze_created_utc": frozen["created_utc"], "control_created_utc": receipt["created_utc"], "execution_manifest_created_utc": execution["created_utc"]},
        "limits": ["Four original-workspace frozen dependencies are not included, so their hashes remain reported rather than reverified here.", "Local receipts, state and hashes support a single prospective execution; they are not an external timestamp or a forensic proof that no unrecorded trial ever existed.", "QuickHull JSON reports are verified as delivered and anchored; no hull or peeling engine is rerun.", "No scientific outcome is recomputed by this script; no new seed, DRAW, RUN, oracle or ANALYZE occurs.", "Production peak memory is not serialized and cannot be certified after exit."],
        "script_sha256": digest(Path(__file__)),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "checks": checks, "frozen_local_matches": frozen_check["matched_count"], "frozen_external_unavailable": len(frozen_check["unavailable"]), "delivery_matches": delivery_check["matched_count"], "delivery_unavailable": len(delivery_check["unavailable"]), "costs": result["costs"]}))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
