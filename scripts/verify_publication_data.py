#!/usr/bin/env python3
"""Lightweight, portable arithmetic verification of the published campaign CSV.

Dependencies: Python 3.10+, NumPy, SciPy.

Example:
  python verify_publication_data.py --campaign-csv campaign_cells_v04.csv \
      --output verification.json --analysis-json analysis_results_v04.json

This script reads the supplied CSV only (plus the optional analysis receipt).
It performs no hull construction, PDE solve, random draw, or new scientific test.
It recomputes the originally specified finite-comparison rule. Geometry, seed
chronology and qualification of the empirical reference errors are out of scope.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import t as student_t

N = 64
RADII = (256, 512)
ORIENTATIONS = ("theta0", "atan_1_8", "pi_over_4")
Q_VALUES = (0.95, 0.90, 0.85)
BETA_LABELS = ("0.30", "0.31", "0.32", "1/3", "0.34", "0.35", "0.36")
BETA_VALUES = (0.30, 0.31, 0.32, 1.0 / 3.0, 0.34, 0.35, 0.36)
DENOMINATOR = 1 << 20
PROVENANCE = "replication_v04_new_draw"
FIELDS = (
    "row_index", "tx", "ty", "denominator", "provenance", "R", "orientation",
    "q_historical", "q_absolute", "beta_label", "beta", "peel_a3_shape_real",
    "peel_a3_shape_imag", "reference_a3_shape_real", "reference_a3_shape_imag",
    "epsilon_interp_a3", "epsilon_ref_a3_qualification_grid_max",
)
FLOAT_FIELDS = (
    "q_historical", "q_absolute", "beta", "peel_a3_shape_real",
    "peel_a3_shape_imag", "reference_a3_shape_real", "reference_a3_shape_imag",
    "epsilon_interp_a3", "epsilon_ref_a3_qualification_grid_max",
)


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def finite_number(text: str | None, name: str, line: int) -> float:
    if text is None or not text.strip():
        raise ValueError(f"CSV line {line}: missing {name}")
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f"CSV line {line}: nonfinite {name}")
    return value


def description(values) -> dict:
    array = np.asarray(values, dtype=float)
    return {
        "n": int(array.size), "mean": float(array.mean()),
        "sd_ddof1": float(array.std(ddof=1)),
        "min": float(array.min()), "max": float(array.max()),
    }


def read_campaign(path: Path):
    groups = defaultdict(dict)
    identities = defaultdict(set)
    eps_by_beta = defaultdict(set)
    full_keys = set()
    count = 0
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError("CSV header differs from the frozen 17-field schema")
        for line, source in enumerate(reader, start=2):
            count += 1
            if None in source:
                raise ValueError(f"CSV line {line}: extra unlabelled fields")
            value = {name: finite_number(source[name], name, line) for name in FLOAT_FIELDS}
            row = int(source["row_index"])
            radius = int(source["R"])
            tx, ty, denominator = (int(source[name]) for name in ("tx", "ty", "denominator"))
            provenance = source["provenance"]
            orientation, label = source["orientation"], source["beta_label"]
            q, q_absolute = value["q_historical"], value["q_absolute"]
            if row not in range(N) or radius not in RADII or orientation not in ORIENTATIONS or q not in Q_VALUES:
                raise ValueError(f"CSV line {line}: row/R/orientation/q outside the frozen grid")
            if label not in BETA_LABELS or value["beta"] != BETA_VALUES[BETA_LABELS.index(label)]:
                raise ValueError(f"CSV line {line}: beta label/value mismatch")
            if denominator != DENOMINATOR or not (0 <= tx < denominator and 0 <= ty < denominator):
                raise ValueError(f"CSV line {line}: invalid dyadic translation")
            if provenance != PROVENANCE:
                raise ValueError(f"CSV line {line}: unexpected translation provenance")
            if not 0.75 <= q_absolute <= 1.0:
                raise ValueError(f"CSV line {line}: q_absolute outside the frozen oracle domain")
            epsilon_interp = value["epsilon_interp_a3"]
            epsilon_ref = value["epsilon_ref_a3_qualification_grid_max"]
            if epsilon_interp < 0 or epsilon_ref < 0:
                raise ValueError(f"CSV line {line}: negative sensitivity")
            key = (row, radius, orientation, q, label)
            if key in full_keys:
                raise ValueError(f"CSV line {line}: duplicate scientific key {key}")
            full_keys.add(key)
            identities[row].add((tx, ty, denominator, provenance))
            eps_by_beta[label].add(epsilon_ref)
            groups[(row, radius, orientation, q)][label] = {
                "peel": complex(value["peel_a3_shape_real"], value["peel_a3_shape_imag"]),
                "reference": complex(value["reference_a3_shape_real"], value["reference_a3_shape_imag"]),
                "q_absolute": q_absolute, "epsilon_interp": epsilon_interp,
                "epsilon_ref": epsilon_ref,
            }
    expected = {
        (row, radius, orientation, q, label)
        for row in range(N) for radius in RADII for orientation in ORIENTATIONS
        for q in Q_VALUES for label in BETA_LABELS
    }
    if count != 8064 or full_keys != expected or len(groups) != 1152:
        raise ValueError("Expected exactly 8064 distinct beta cells and 1152 complete target groups")
    if set(identities) != set(range(N)) or any(len(values) != 1 for values in identities.values()):
        raise ValueError("Each of the 64 row indices must have one consistent translation identity")
    identity_by_row = {row: next(iter(values)) for row, values in identities.items()}
    if len(set(identity_by_row.values())) != N or len({value[:3] for value in identity_by_row.values()}) != N:
        raise ValueError("The 64 translation coordinates must be distinct")
    if any(len(values) != 1 for values in eps_by_beta.values()):
        raise ValueError("Reference sensitivity must be constant within each beta")
    for key, group in groups.items():
        if set(group) != set(BETA_LABELS):
            raise ValueError(f"Incomplete beta group {key}")
        first = group["0.32"]
        for record in group.values():
            if any(record[name] != first[name] for name in ("peel", "q_absolute", "epsilon_interp")):
                raise ValueError(f"The paired target differs across beta in group {key}")
    return groups, identity_by_row, {label: next(iter(values)) for label, values in eps_by_beta.items()}


def recompute(path: Path) -> dict:
    groups, identities, epsilons = read_campaign(path)
    cells = {}
    losses = defaultdict(list)
    for key, group in groups.items():
        p = group["0.32"]["peel"]
        r032, r13 = group["0.32"]["reference"], group["1/3"]["reference"]
        ei = group["0.32"]["epsilon_interp"]
        e032, e13 = group["0.32"]["epsilon_ref"], group["1/3"]["epsilon_ref"]
        l032, l13 = abs(p - r032) ** 2, abs(p - r13) ** 2
        b_interp = 2 * ei * abs(r13 - r032)
        b_ref13 = 2 * abs(p - r13) * e13 + e13**2
        b_ref032 = 2 * abs(p - r032) * e032 + e032**2
        b_mixed = 2 * ei * (e13 + e032)
        bound = b_interp + b_ref13 + b_ref032 + b_mixed
        cells[key] = {
            "D": float(l13 - l032), "B": float(bound),
            "D_minus_B": float(l13 - l032 - bound),
            "D_plus_B": float(l13 - l032 + bound),
        }
        for label, record in group.items():
            loss = abs(p - record["reference"]) ** 2
            if not math.isfinite(loss):
                raise ValueError(f"Loss overflow in group {key}")
            losses[(key[1], label)].append(loss)
        if not all(math.isfinite(value) for value in cells[key].values()):
            raise ValueError(f"Contrast/sensitivity overflow in group {key}")
    w = np.asarray([
        math.fsum(cells[(row, 512, orientation, q)]["D_minus_B"] for orientation in ORIENTATIONS for q in Q_VALUES) / 9
        for row in range(N)
    ])
    mean, sd = float(w.mean()), float(w.std(ddof=1))
    critical = float(student_t.ppf(0.975, N - 1))
    lcb = mean - critical * sd / math.sqrt(N)
    loo = []
    for omitted in range(N):
        retained = np.delete(w, omitted)
        mean_loo = float(retained.mean())
        lcb_loo = float(mean_loo - student_t.ppf(0.975, N - 2) * retained.std(ddof=1) / math.sqrt(N - 1))
        loo.append({"omitted_row_index": omitted, "mean_W": mean_loo, "lcb": lcb_loo})
    pi4 = {}
    for radius in RADII:
        pi4[str(radius)] = {}
        for q in Q_VALUES:
            selected = [cells[(row, radius, "pi_over_4", q)] for row in range(N)]
            pi4[str(radius)][f"{q:.2f}"] = {
                "D_raw": description([cell["D"] for cell in selected]),
                "B_empirical": description([cell["B"] for cell in selected]),
                "D_minus_B": description([cell["D_minus_B"] for cell in selected]),
                "D_plus_B": description([cell["D_plus_B"] for cell in selected]),
                "D_plus_B_negative_count": sum(cell["D_plus_B"] < 0 for cell in selected),
                "role": "descriptive secondary arithmetic; no additional inferential decision",
            }
    identity_rows = [
        {"row_index": row, "tx": identities[row][0], "ty": identities[row][1],
         "denominator": identities[row][2], "provenance": identities[row][3], "W": float(w[row])}
        for row in range(N)
    ]
    return {
        "status": "VERIFIED_DATA_RECOMPOSITION",
        "scientific_decision_from_original_rule": "FINITE_R512_PREFERENCE_0_32" if lcb > 0 else "NON_DISCRIMINE",
        "counts": {"translation_rows": N, "beta_cells": 8064, "paired_target_groups": 1152, "cells_per_primary_row": 9},
        "primary": {
            "n_independent_translation_rows": N, "mean_W": mean, "sd_W_ddof1": sd,
            "student_t_quantile_0_975_df63": critical, "lower_confidence_bound_nominal": lcb,
            "row_score": "equal mean over 3 orientations x 3 q at R=512 of D-B",
            "W_min": float(w.min()), "W_max": float(w.max()),
            "W_positive_count": int(np.sum(w > 0)), "row_scores_with_exact_identity": identity_rows,
            "coverage_label": "nominal approximate Student bound",
        },
        "leave_one_translation_out": {
            "minimum_mean_W": min(row["mean_W"] for row in loo),
            "maximum_mean_W": max(row["mean_W"] for row in loo),
            "minimum_lcb": min(row["lcb"] for row in loo),
            "maximum_lcb": max(row["lcb"] for row in loo),
            "positive_count_of_64": sum(row["lcb"] > 0 for row in loo), "rows": loo,
        },
        "pi_over_4_by_radius_and_q": pi4,
        "mean_losses_by_radius_and_beta": {
            str(radius): {label: float(np.mean(losses[radius, label])) for label in BETA_LABELS}
            for radius in RADII
        },
        "epsilon_ref_values_from_CSV": epsilons,
        "sensitivity_formula": "B=2*ei*|r13-r032|+2*|p-r13|*e13+e13^2+2*|p-r032|*e032+e032^2+2*ei*(e13+e032)",
        "scope": [
            "CSV-only arithmetic verification of the original finite primary comparison",
            "no new statistical test, hull construction, PDE solve or random draw",
            "reference sensitivities are read from CSV; their qualification is not independently established here",
            "no geometric, chronological or freeze-integrity audit",
            "empirical sensitivities are not certified enclosures; Student coverage is approximate",
            "no uniformity, exact exponent, continuous optimum or asymptotic conclusion",
        ],
    }


def compare_receipt(result: dict, path: Path) -> dict:
    receipt = json.loads(path.read_text(encoding="utf-8"))
    comparisons = []
    def check(name, observed, expected):
        if isinstance(expected, (int, float)):
            if not math.isfinite(float(expected)):
                raise ValueError(f"Nonfinite receipt value: {name}")
            equal = math.isclose(float(observed), float(expected), rel_tol=1e-12, abs_tol=1e-25)
            difference = abs(float(observed) - float(expected))
        else:
            equal = observed == expected
            difference = None
        comparisons.append({"field": name, "matches": equal, "absolute_difference": difference})
    check("status", result["scientific_decision_from_original_rule"], receipt["status"])
    for key in ("n_independent_translation_rows", "mean_W", "sd_W_ddof1", "student_t_quantile_0_975_df63", "lower_confidence_bound_nominal"):
        check("primary." + key, result["primary"][key], receipt["primary"][key])
    received_scores = {int(row["row_index"]): row for row in receipt["primary"]["row_scores_with_exact_identity"]}
    if len(received_scores) != N:
        raise ValueError("Optional analysis receipt does not contain 64 row scores")
    for row in result["primary"]["row_scores_with_exact_identity"]:
        other = received_scores[row["row_index"]]
        for key in ("tx", "ty", "denominator", "provenance", "W"):
            check(f"row[{row['row_index']}].{key}", row[key], other[key])
    received_loo = {int(row["omitted_row_index"]): row for row in receipt["leave_one_translation_out"]["rows"]}
    if len(received_loo) != N:
        raise ValueError("Optional analysis receipt does not contain 64 leave-one-out rows")
    for row in result["leave_one_translation_out"]["rows"]:
        for key in ("mean_W", "lcb"):
            check(f"leave_one_out[{row['omitted_row_index']}].{key}", row[key], received_loo[row["omitted_row_index"]][key])
    if not all(item["matches"] for item in comparisons):
        bad = [item["field"] for item in comparisons if not item["matches"]]
        raise ValueError("Optional analysis receipt mismatch: " + ", ".join(bad))
    return {
        "status": "MATCHES_WITH_FLOATING_POINT_TOLERANCE", "sha256": digest(path),
        "comparison_count": len(comparisons),
        "maximum_absolute_numerical_difference": max(item["absolute_difference"] or 0.0 for item in comparisons),
        "relative_tolerance": 1e-12, "absolute_tolerance": 1e-25,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--campaign-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--analysis-json", type=Path, default=None)
    args = parser.parse_args()
    output = args.output.resolve()
    inputs = [args.campaign_csv.resolve()] + ([args.analysis_json.resolve()] if args.analysis_json else [])
    if output in inputs:
        parser.error("--output must differ from the input files")
    try:
        result = recompute(args.campaign_csv)
        result["campaign_csv_sha256"] = digest(args.campaign_csv)
        result["verifier_source_sha256"] = digest(Path(__file__))
        result["environment"] = {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}
        if args.analysis_json:
            result["optional_analysis_receipt"] = compare_receipt(result, args.analysis_json)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(json.dumps({"status": result["status"], "decision": result["scientific_decision_from_original_rule"],
                          "mean_W": result["primary"]["mean_W"], "LCB": result["primary"]["lower_confidence_bound_nominal"],
                          "output": str(output)}))
        return 0
    except Exception as exc:
        failure = {"status": "VERIFICATION_ERROR", "error_type": type(exc).__name__, "error": str(exc)}
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(failure, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(failure), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
