"""Fail-closed confirmatory analyzer for the future v0.4 replication.

This file is prepared and schema-checked in the current pass, but no
confirmatory analysis is executed.  Scientific analysis requires one complete
frozen campaign containing exactly 64 unique translation rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PREP_ROOT = HERE.parent
WORKSPACE = PREP_ROOT.parent
RUNNER_PATH = PREP_ROOT / "launcher" / "campaign_runner_v04.py"
QUALIFICATION_PATH = PREP_ROOT / "pilot" / "reference_qualification_v04.json"
CONFIG_PATH = PREP_ROOT / "protocol_prepare_v04.json"
STATE_PATH = PREP_ROOT / "launcher_state_v04.json"

ORIENTATIONS = ("theta0", "atan_1_8", "pi_over_4")
RADII = (256, 512)
Q_VALUES = (0.95, 0.90, 0.85)
BETA_LABELS = ("0.30", "0.31", "0.32", "1/3", "0.34", "0.35", "0.36")
BETA_VALUES = (0.30, 0.31, 0.32, 1.0 / 3.0, 0.34, 0.35, 0.36)
EXPECTED_FIELDS = (
    "row_index", "tx", "ty", "denominator", "provenance", "R", "orientation",
    "q_historical", "q_absolute", "beta_label", "beta", "peel_a3_shape_real",
    "peel_a3_shape_imag", "reference_a3_shape_real", "reference_a3_shape_imag",
    "epsilon_interp_a3", "epsilon_ref_a3_qualification_grid_max",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def finite_float(value: str, name: str) -> float:
    if value is None or value.strip() == "":
        raise ValueError(f"missing {name}")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"nonfinite {name}")
    return result


def validate_campaign(campaign_dir: Path):
    campaign_dir = campaign_dir.resolve()
    complete_path = campaign_dir / "campaign_complete_manifest_v04.json"
    if not complete_path.is_file():
        raise RuntimeError("REPLICATION_INCOMPLETE: completion manifest is absent")
    if not STATE_PATH.is_file():
        raise RuntimeError("analysis requires the launcher's anchored state")
    launcher_state = read_json(STATE_PATH)
    if launcher_state.get("next_action") not in ("ANALYZE", "COMPLETE"):
        raise RuntimeError("ANALYZE is impossible before the launcher's complete RUN transition")
    if (
        Path(launcher_state.get("campaign_complete_manifest_path", "")).resolve() != complete_path.resolve()
        or launcher_state.get("campaign_complete_manifest_sha256") != sha256(complete_path)
    ):
        raise RuntimeError("completion manifest is not anchored by launcher_state_v04.json")
    complete = read_json(complete_path)
    if complete.get("status") != "COMPLETE_64_ROWS" or complete.get("row_count") != 64:
        raise RuntimeError("REPLICATION_INCOMPLETE: exactly 64 complete rows are required")
    if complete.get("geometry_divergence_count") != 0:
        raise RuntimeError("NO_PRIMARY_DECISION: recorded geometry divergence")
    if complete.get("second_engine_row_indices_zero_based") != [0, 8, 16, 24, 32, 40, 48, 56]:
        raise RuntimeError("second-engine row schedule mismatch")

    runner = load_module("v04_analysis_runner_guards", RUNNER_PATH)
    frozen_sha = complete["frozen_manifest_sha256"]
    runner.verify_frozen_manifest(frozen_sha)
    translation_path = Path(complete["translation_manifest_path"])
    if not translation_path.is_file() or sha256(translation_path) != complete["translation_manifest_sha256"]:
        raise RuntimeError("translation manifest missing or changed")
    translations = runner.validate_translations(translation_path, frozen_sha)

    for record in complete.get("outputs", []):
        path = Path(record["path"]).resolve()
        try:
            path.relative_to(campaign_dir)
        except ValueError as exc:
            raise RuntimeError("campaign manifest references an output outside campaign directory") from exc
        if not path.is_file() or sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise RuntimeError(f"campaign output missing or changed: {path.name}")
    cells_path = Path(complete["cells_path"])
    if not cells_path.is_file() or sha256(cells_path) != complete["cells_sha256"]:
        raise RuntimeError("campaign cell table missing or changed")
    qualification = read_json(QUALIFICATION_PATH)
    if qualification.get("status") != "QUALIFIED_EMPIRICAL_ORACLE_NOT_INTERVAL_CERTIFIED":
        raise RuntimeError("NO_PRIMARY_DECISION: oracle is not qualified")
    if sha256(QUALIFICATION_PATH) != complete.get("oracle_qualification_sha256"):
        raise RuntimeError("NO_PRIMARY_DECISION: oracle qualification hash mismatch")
    if qualification.get("config_sha256") != sha256(CONFIG_PATH):
        raise RuntimeError("NO_PRIMARY_DECISION: oracle config hash mismatch")

    # Every one of the eight prespecified checks must be present and successful.
    for row_index in complete["second_engine_row_indices_zero_based"]:
        check_path = campaign_dir / "rows" / f"row_{row_index:02d}_quickhull_check.json"
        check = read_json(check_path)
        if (
            check.get("status") != "PASS_COMPLETE_LAYER_BY_LAYER_SECOND_ENGINE"
            or not check.get("cumulative_layer_record_sha256_equal")
        ):
            raise RuntimeError(f"NO_PRIMARY_DECISION: geometry divergence at row {row_index}")

    expected_rows = {int(row["row_index"]): row for row in translations}
    qualification_eps = qualification["epsilon_ref_qualification_grid_max_final"]
    expected_eps = {
        label: float(qualification_eps[f"{value:.17g}"])
        for label, value in zip(BETA_LABELS, BETA_VALUES)
    }
    records = []
    full_keys = set()
    base_groups = defaultdict(list)
    row_keys = defaultdict(set)
    with cells_path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != EXPECTED_FIELDS:
            raise RuntimeError("campaign cell schema mismatch")
        for source in reader:
            row_index = int(source["row_index"])
            if row_index not in expected_rows:
                raise RuntimeError("cell references an unknown row")
            expected = expected_rows[row_index]
            tx, ty, denominator = int(source["tx"]), int(source["ty"]), int(source["denominator"])
            provenance = source["provenance"]
            if (tx, ty, denominator, provenance) != (
                int(expected["tx"]), int(expected["ty"]), int(expected["denominator"]), expected["provenance"]
            ):
                raise RuntimeError("join failure on (tx,ty,D,provenance)")
            radius = int(source["R"])
            orientation = source["orientation"]
            q_historical = finite_float(source["q_historical"], "q_historical")
            q_absolute = finite_float(source["q_absolute"], "q_absolute")
            label = source["beta_label"]
            beta = finite_float(source["beta"], "beta")
            if radius not in RADII or orientation not in ORIENTATIONS or q_historical not in Q_VALUES:
                raise RuntimeError("cell outside the frozen R/orientation/q grid")
            if label not in BETA_LABELS or beta != BETA_VALUES[BETA_LABELS.index(label)]:
                raise RuntimeError("cell outside the frozen beta grid")
            if not 0.75 <= q_absolute <= 1.0:
                raise RuntimeError("q_absolute outside qualified oracle domain")
            full_key = (row_index, tx, ty, denominator, provenance, radius, orientation, q_historical, label)
            if full_key in full_keys:
                raise RuntimeError("duplicate scientific cell")
            full_keys.add(full_key)
            peel = complex(
                finite_float(source["peel_a3_shape_real"], "peel real"),
                finite_float(source["peel_a3_shape_imag"], "peel imag"),
            )
            reference = complex(
                finite_float(source["reference_a3_shape_real"], "reference real"),
                finite_float(source["reference_a3_shape_imag"], "reference imag"),
            )
            epsilon_interp = finite_float(source["epsilon_interp_a3"], "epsilon_interp_a3")
            epsilon_ref = finite_float(
                source["epsilon_ref_a3_qualification_grid_max"],
                "epsilon_ref_a3_qualification_grid_max",
            )
            if epsilon_interp < 0 or epsilon_ref < 0 or epsilon_ref != expected_eps[label]:
                raise RuntimeError("cell sensitivity value mismatch")
            record = {
                "row_index": row_index, "tx": tx, "ty": ty, "denominator": denominator,
                "provenance": provenance, "R": radius, "orientation": orientation,
                "q_historical": q_historical, "q_absolute": q_absolute, "beta_label": label,
                "beta": beta, "peel": peel, "reference": reference,
                "epsilon_interp": epsilon_interp, "epsilon_ref": epsilon_ref,
            }
            records.append(record)
            base_key = (row_index, radius, orientation, q_historical)
            base_groups[base_key].append(record)
            row_keys[row_index].add((radius, orientation, q_historical, label))

    if len(records) != 64 * 2 * 3 * 3 * 7 or len(base_groups) != 64 * 2 * 3 * 3:
        raise RuntimeError("REPLICATION_INCOMPLETE: expected exactly 8064 beta cells")
    expected_per_row = {
        (radius, orientation, q_value, label)
        for radius in RADII for orientation in ORIENTATIONS for q_value in Q_VALUES for label in BETA_LABELS
    }
    if set(row_keys) != set(range(64)) or any(keys != expected_per_row for keys in row_keys.values()):
        raise RuntimeError("REPLICATION_INCOMPLETE: 2R x 3theta x 3q x 7beta grid is incomplete")
    for base_key, group in base_groups.items():
        if len(group) != 7 or {row["beta_label"] for row in group} != set(BETA_LABELS):
            raise RuntimeError(f"incomplete beta group {base_key}")
        first = group[0]
        if any(
            row["peel"] != first["peel"]
            or row["q_absolute"] != first["q_absolute"]
            or row["epsilon_interp"] != first["epsilon_interp"]
            for row in group[1:]
        ):
            raise RuntimeError("paired target fields differ across beta values")

    expected_order = [
        (row_index, radius, orientation, q_value, label)
        for row_index in range(64)
        for radius in RADII
        for orientation in ORIENTATIONS
        for q_value in Q_VALUES
        for label in BETA_LABELS
    ]
    actual_order = [
        (row["row_index"], row["R"], row["orientation"], row["q_historical"], row["beta_label"])
        for row in records
    ]
    if actual_order != expected_order:
        raise RuntimeError("campaign CSV order differs from the frozen row/R/orientation/q/beta order")

    # Reconcile every CSV target with the raw main-engine NPZ.  The CSV is not
    # treated as an independent source of geometry or interpolation bounds.
    orientation_pq = np.asarray(((1, 0), (8, 1), (1, 1)), dtype=np.int64)
    for row_index in range(64):
        raw_path = campaign_dir / "rows" / f"row_{row_index:02d}_main.npz"
        if not raw_path.is_file():
            raise RuntimeError(f"REPLICATION_INCOMPLETE: raw main row {row_index} is absent")
        expected_translation = expected_rows[row_index]
        with np.load(raw_path, allow_pickle=False) as raw:
            if raw["radii"].tolist() != [256, 256, 256, 512, 512, 512]:
                raise RuntimeError("raw row radius order mismatch")
            if not np.array_equal(raw["orientation_pq"], np.vstack((orientation_pq, orientation_pq))):
                raise RuntimeError("raw row orientation order mismatch")
            expected_translation_array = np.tile(
                np.asarray([[expected_translation["tx"], expected_translation["ty"], 1048576]], dtype=np.int64),
                (6, 1),
            )
            if not np.array_equal(raw["translations"], expected_translation_array):
                raise RuntimeError("raw row coordinates differ from (tx,ty,D,provenance) join")
            if raw["target_q_historical"].tolist() != list(Q_VALUES):
                raise RuntimeError("raw row q target order mismatch")
            modes = raw["fourier_modes"].tolist()
            if modes != [0, 1, 3, 6, 9, 12]:
                raise RuntimeError("raw row Fourier mode order mismatch")
            mode3 = modes.index(3)
            if (
                raw["a_shape_exp_plus_i_m_theta"].shape != (6, 3, 6)
                or raw["epsilon_interpolation_a3"].shape != (6, 3)
                or raw["target_q_absolute"].shape != (6, 3)
            ):
                raise RuntimeError("raw row target array shape mismatch")
            if max(
                float(np.max(raw["before_c0_complex_residual"])),
                float(np.max(raw["after_c0_complex_residual"])),
                float(np.max(raw["interpolated_c0_complex_residual_before_correction"])),
            ) > 5e-14:
                raise RuntimeError("raw row c0 identity/target residual exceeds frozen tolerance")
            sentinel = int(raw["layer0_containment_sentinel"])
            containment = np.asarray(raw["layer_minimum_containment_cross_from_previous"], dtype=np.int64)
            if np.any(containment[containment != sentinel] < 0) or np.any(
                np.asarray(raw["minimum_all_adjacent_layer_containment_cross"], dtype=np.int64) < 0
            ):
                raise RuntimeError("raw row adjacent-layer nesting failure")
            reconstructed_global = raw["interpolated_coefficients_corrected_c0"] / raw[
                "target_q_absolute"
            ][:, :, None]
            if not np.array_equal(reconstructed_global, raw["a_global"]):
                raise RuntimeError("raw row a_global normalization mismatch")
            reconstructed_shape = raw["a_global"] * np.exp(
                1j * raw["fourier_modes"][None, None, :] * raw["orientation_theta"][:, None, None]
            )
            if not np.array_equal(reconstructed_shape, raw["a_shape_exp_plus_i_m_theta"]):
                raise RuntimeError("raw row shape-frame convention mismatch")
            for radius_index, radius in enumerate(RADII):
                for orientation_index, orientation in enumerate(ORIENTATIONS):
                    trajectory_index = radius_index * 3 + orientation_index
                    for target_index, q_value in enumerate(Q_VALUES):
                        base = base_groups[(row_index, radius, orientation, q_value)][0]
                        q_formula = q_value * float(raw["initial_perimeter"][trajectory_index]) / (
                            2 * math.pi * radius
                        )
                        q_raw = float(raw["target_q_absolute"][trajectory_index, target_index])
                        if not (
                            math.isclose(q_formula, q_raw, rel_tol=0.0, abs_tol=2e-15)
                            and base["q_absolute"] == q_raw
                        ):
                            raise RuntimeError("q_abs is not q_hist*P_initial/(2*pi*R)")
                        peel_raw = complex(raw["a_shape_exp_plus_i_m_theta"][trajectory_index, target_index, mode3])
                        epsilon_raw = float(raw["epsilon_interpolation_a3"][trajectory_index, target_index])
                        if base["peel"] != peel_raw or base["epsilon_interp"] != epsilon_raw:
                            raise RuntimeError("CSV peel a3/epsilon_interp differs from raw main NPZ")

    # Reconcile the flattened beta-expanded CSV with the frozen oracle NPZ,
    # including the exact q order and per-beta empirical-grid sensitivity.
    reference_path = campaign_dir / "campaign_reference_v04.npz"
    metadata_path = reference_path.with_suffix(".json")
    if not reference_path.is_file() or not metadata_path.is_file():
        raise RuntimeError("REPLICATION_INCOMPLETE: frozen oracle outputs are absent")
    base_sequence = [records[index] for index in range(0, len(records), 7)]
    with np.load(reference_path, allow_pickle=False) as reference_raw:
        if reference_raw["beta_grid"].tolist() != list(BETA_VALUES):
            raise RuntimeError("oracle beta order mismatch")
        reference_modes = reference_raw["modes"].tolist()
        if reference_modes != [0, 3, 6, 9, 12]:
            raise RuntimeError("oracle mode order mismatch")
        if not np.array_equal(
            reference_raw["q_absolute"],
            np.asarray([row["q_absolute"] for row in base_sequence], dtype=float),
        ):
            raise RuntimeError("oracle q_abs order differs from finite-grid target order")
        eps_array = np.asarray(reference_raw["epsilon_ref_qualification_grid_max"], dtype=float)
        if not np.array_equal(eps_array, np.asarray([expected_eps[label] for label in BETA_LABELS])):
            raise RuntimeError("oracle NPZ empirical sensitivities differ from qualification")
        reference_coefficients = np.asarray(reference_raw["coefficients"])
        reference_a = np.asarray(reference_raw["normalized_coefficients"])
        if reference_a.shape != (7, 64 * 2 * 3 * 3, 5):
            raise RuntimeError("oracle normalized coefficient shape mismatch")
        q_array = np.asarray(reference_raw["q_absolute"], dtype=float)
        if float(np.max(np.abs(reference_coefficients[:, :, 0] - q_array[None, :]))) > 2e-12:
            raise RuntimeError("oracle complex c0 differs from q_abs beyond qualified tolerance")
        reconstructed_reference_a = reference_coefficients / reference_coefficients[:, :, [0]].real
        if not np.array_equal(reconstructed_reference_a, reference_a):
            raise RuntimeError("oracle a_m normalization differs from coefficients/Re(c0)")
        mode3 = reference_modes.index(3)
        for base_index, base in enumerate(base_sequence):
            for beta_index, label in enumerate(BETA_LABELS):
                csv_record = records[base_index * 7 + beta_index]
                if (
                    csv_record["beta_label"] != label
                    or csv_record["reference"] != complex(reference_a[beta_index, base_index, mode3])
                    or csv_record["epsilon_ref"] != float(eps_array[beta_index])
                ):
                    raise RuntimeError("CSV reference/epsilon/order differs from oracle NPZ")
    reference_metadata = read_json(metadata_path)
    if (
        reference_metadata.get("status") != "EVALUATED_BY_FROZEN_ORACLE"
        or reference_metadata.get("q_count") != 64 * 2 * 3 * 3
        or reference_metadata.get("output_sha256") != sha256(reference_path)
        or reference_metadata.get("qualification_sha256") != sha256(QUALIFICATION_PATH)
    ):
        raise RuntimeError("oracle evaluation metadata/hash mismatch")
    return {
        "complete": complete,
        "translations": translations,
        "records": records,
        "base_groups": base_groups,
        "cells_path": cells_path,
        "frozen_manifest_sha256": frozen_sha,
        "raw_npz_rows_reconciled": 64,
        "oracle_targets_reconciled": 64 * 2 * 3 * 3,
    }


def descriptive(values) -> dict:
    array = np.asarray(values, dtype=float)
    return {
        "n": int(len(array)),
        "mean": float(np.mean(array)),
        "sd_ddof1": float(np.std(array, ddof=1)) if len(array) > 1 else None,
        "min": float(np.min(array)),
        "max": float(np.max(array)),
    }


def analyze(validated, output_dir: Path) -> dict:
    from scipy.stats import t as student_t

    if output_dir.exists():
        raise FileExistsError("refusing to overwrite an analysis directory")
    output_dir.mkdir(parents=True)
    contrast_rows = []
    losses = defaultdict(list)
    for base_key, group in validated["base_groups"].items():
        values = {row["beta_label"]: row for row in group}
        peel = values["0.32"]["peel"]
        r032 = values["0.32"]["reference"]
        r13 = values["1/3"]["reference"]
        ei = values["0.32"]["epsilon_interp"]
        e032 = values["0.32"]["epsilon_ref"]
        e13 = values["1/3"]["epsilon_ref"]
        loss032 = abs(peel - r032) ** 2
        loss13 = abs(peel - r13) ** 2
        b_interp = 2 * ei * abs(r13 - r032)
        b_ref13 = 2 * abs(peel - r13) * e13 + e13**2
        b_ref032 = 2 * abs(peel - r032) * e032 + e032**2
        b_mixed = 2 * ei * (e13 + e032)
        bound = b_interp + b_ref13 + b_ref032 + b_mixed
        row_index, radius, orientation, q_value = base_key
        identity = values["0.32"]
        contrast_rows.append({
            "row_index": row_index, "tx": identity["tx"], "ty": identity["ty"],
            "denominator": identity["denominator"], "provenance": identity["provenance"],
            "R": radius, "orientation": orientation, "q_historical": q_value,
            "q_absolute": identity["q_absolute"], "epsilon_interp_a3": ei,
            "epsilon_ref_1_over_3": e13, "epsilon_ref_0_32": e032,
            "loss_1_over_3": float(loss13),
            "loss_0_32": float(loss032), "d_raw": float(loss13 - loss032),
            "B_common_interpolation": float(b_interp),
            "B_reference_1_over_3": float(b_ref13),
            "B_reference_0_32": float(b_ref032),
            "B_mixed": float(b_mixed),
            "B_empirical": float(bound), "d_minus_B": float(loss13 - loss032 - bound),
        })
        for label, row in values.items():
            losses[(radius, orientation, q_value, label)].append(abs(peel - row["reference"]) ** 2)

    by_row_radius = defaultdict(list)
    for row in contrast_rows:
        by_row_radius[(row["row_index"], row["R"])].append(row)
    if any(len(rows) != 9 for rows in by_row_radius.values()):
        raise AssertionError("row-radius pairing failure")
    w = np.asarray([
        math.fsum(row["d_minus_B"] for row in by_row_radius[(index, 512)]) / 9
        for index in range(64)
    ])
    mean_w = float(np.mean(w))
    sd_w = float(np.std(w, ddof=1))
    critical = float(student_t.ppf(0.975, 63))
    lcb = mean_w - critical * sd_w / math.sqrt(64)
    leave_one_out = []
    for omitted in range(64):
        retained = np.delete(w, omitted)
        loo_lcb = float(np.mean(retained) - student_t.ppf(0.975, 62) * np.std(retained, ddof=1) / math.sqrt(63))
        identity = validated["translations"][omitted]
        leave_one_out.append({
            "omitted_row_index": omitted,
            "omitted_tx": int(identity["tx"]),
            "omitted_ty": int(identity["ty"]),
            "omitted_denominator": int(identity["denominator"]),
            "omitted_provenance": identity["provenance"],
            "mean_W": float(np.mean(retained)),
            "lcb": loo_lcb,
        })

    row_raw_by_radius = {
        (index, radius): math.fsum(row["d_raw"] for row in by_row_radius[(index, radius)]) / 9
        for index in range(64) for radius in RADII
    }
    paired_radius = [row_raw_by_radius[(index, 512)] - row_raw_by_radius[(index, 256)] for index in range(64)]
    by_angle = {
        orientation: descriptive([
            math.fsum(
                row["d_raw"] for row in contrast_rows
                if row["row_index"] == index and row["R"] == 512 and row["orientation"] == orientation
            ) / 3
            for index in range(64)
        ])
        for orientation in ORIENTATIONS
    }
    by_q = {
        f"{q_value:.2f}": descriptive([
            math.fsum(
                row["d_raw"] for row in contrast_rows
                if row["row_index"] == index and row["R"] == 512 and row["q_historical"] == q_value
            ) / 3
            for index in range(64)
        ])
        for q_value in Q_VALUES
    }
    dedicated = {
        str(radius): descriptive([
            row["d_raw"] for row in contrast_rows
            if row["R"] == radius and row["orientation"] == "pi_over_4" and row["q_historical"] == 0.95
        ])
        for radius in RADII
    }
    loss_summary = [
        {"R": key[0], "orientation": key[1], "q_historical": key[2], "beta_label": key[3], **descriptive(values)}
        for key, values in sorted(losses.items(), key=lambda item: (
            item[0][0], ORIENTATIONS.index(item[0][1]), Q_VALUES.index(item[0][2]), BETA_LABELS.index(item[0][3])
        ))
    ]
    result = {
        "status": "FINITE_R512_PREFERENCE_0_32" if lcb > 0 else "NON_DISCRIMINE",
        "primary": {
            "n_independent_translation_rows": 64,
            "row_score": "equal mean over 3 orientations x 3 q at R=512 of d-B",
            "mean_W": mean_w,
            "sd_W_ddof1": sd_w,
            "student_t_quantile_0_975_df63": critical,
            "lower_confidence_bound_nominal": lcb,
            "decision_rule": "LCB>0",
            "coverage_label": "nominal approximate Student bound",
            "sensitivity_label": "B uses empirical qualification-grid maxima, not certified enclosures",
            "uniformity_label": "a positive mean result does not imply uniformity over angle or contraction",
            "row_scores_with_exact_identity": [
                {
                    "row_index": index,
                    "tx": int(validated["translations"][index]["tx"]),
                    "ty": int(validated["translations"][index]["ty"]),
                    "denominator": int(validated["translations"][index]["denominator"]),
                    "provenance": validated["translations"][index]["provenance"],
                    "W": float(w[index]),
                }
                for index in range(64)
            ],
        },
        "leave_one_translation_out": {
            "minimum_mean_W": min(row["mean_W"] for row in leave_one_out),
            "maximum_mean_W": max(row["mean_W"] for row in leave_one_out),
            "minimum_lcb": min(row["lcb"] for row in leave_one_out),
            "maximum_lcb": max(row["lcb"] for row in leave_one_out),
            "positive_count_of_64": sum(row["lcb"] > 0 for row in leave_one_out),
            "rows": leave_one_out,
        },
        "predeclared_secondary": {
            "raw_paired_D512_minus_D256_by_translation": descriptive(paired_radius),
            "raw_D_at_R512_by_orientation": by_angle,
            "raw_D_at_R512_by_q": by_q,
            "raw_D_pi_over_4_q_0_95_by_R": dedicated,
            "all_seven_beta_losses": loss_summary,
            "cannot_rescue_negative_primary": True,
        },
        "sensitivity_formula": (
            "B=2*ei*|r13-r032|+2*|p-r13|*e13+e13^2+"
            "2*|p-r032|*e032+e032^2+2*ei*(e13+e032)"
        ),
        "frozen_manifest_sha256": validated["frozen_manifest_sha256"],
        "campaign_cells_sha256": sha256(validated["cells_path"]),
        "scope_limits": [
            "finite-scale paired comparison only",
            "no continuous beta optimum and no exact exponent",
            "no proof or refutation of an affine asymptotic law",
            "no H_rot or mode-4 claim",
        ],
    }

    contrast_path = output_dir / "cell_contrasts_v04.csv"
    with contrast_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(contrast_rows[0]))
        writer.writeheader()
        writer.writerows(contrast_rows)
    result_path = output_dir / "analysis_results_v04.json"
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    manifest = {
        "status": "ANALYZED_COMPLETE_64_ROWS",
        "result": {"path": str(result_path.resolve()), "sha256": sha256(result_path)},
        "cell_contrasts": {"path": str(contrast_path.resolve()), "sha256": sha256(contrast_path)},
        "source_sha256": sha256(Path(__file__).resolve()),
    }
    (output_dir / "analysis_manifest_v04.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return result


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--schema-check", action="store_true")
    group.add_argument("--validate-only", action="store_true")
    group.add_argument("--analyze", action="store_true")
    parser.add_argument("--campaign-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.schema_check:
        print(json.dumps({
            "status": "SCHEMA_ONLY_NO_SCIENTIFIC_ANALYSIS",
            "expected_rows": 64,
            "expected_cells": 64 * 2 * 3 * 3 * 7,
            "fields": EXPECTED_FIELDS,
        }))
        return
    if args.campaign_dir is None:
        parser.error("--validate-only and --analyze require --campaign-dir")
    validated = validate_campaign(args.campaign_dir)
    if args.validate_only:
        print(json.dumps({
            "status": "VALIDATED_COMPLETE_INPUT_NO_SCIENTIFIC_ANALYSIS",
            "rows": 64,
            "cells": len(validated["records"]),
        }))
        return
    if args.output_dir is None:
        parser.error("--analyze requires --output-dir")
    print(json.dumps(analyze(validated, args.output_dir)))


if __name__ == "__main__":
    main()
