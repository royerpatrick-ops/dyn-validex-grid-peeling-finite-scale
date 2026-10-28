"""Prepared v0.4 campaign collector; never invoked during preparation.

The launcher may call :func:`run_campaign` only after an externally supplied
seed has produced a frozen 64-row translation manifest.  Geometry is collected
for every prespecified cell.  The continuous oracle is evaluated only after
all finite-grid rows exist, using the frozen mechanical rule and no loss
feedback.  There is no efficacy/futility stop.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PREP_ROOT = HERE.parent
WORKSPACE = PREP_ROOT.parent
CONFIG_PATH = PREP_ROOT / "protocol_prepare_v04.json"
PLAN_PATH = PREP_ROOT / "plan" / "PLAN_REPLICATION_INDEPENDANTE_v0_2.json"
PILOT_MODULE_PATH = PREP_ROOT / "pilot" / "run_historical_cost_pilot_v04.py"
ORACLE_PATH = PREP_ROOT / "pilot" / "reference_oracle_v04.py"
QUALIFICATION_PATH = PREP_ROOT / "pilot" / "reference_qualification_v04.json"
FROZEN_MANIFEST_PATH = PREP_ROOT / "frozen_manifest_v04.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def json_ready(value):
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    return value


def write_json_new(path: Path, value) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(json_ready(value), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def verify_frozen_manifest(expected_sha256: str) -> dict:
    if not FROZEN_MANIFEST_PATH.is_file() or sha256(FROZEN_MANIFEST_PATH) != expected_sha256:
        raise RuntimeError("frozen manifest missing or SHA-256 mismatch")
    manifest = read_json(FROZEN_MANIFEST_PATH)
    if manifest.get("status") != "FROZEN_BEFORE_DRAW":
        raise RuntimeError("manifest is not in the frozen pre-draw state")
    current_environment = {
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "mpmath", "threadpoolctl", "psutil")
        },
        "threadpool_limits": 1,
    }
    if manifest.get("environment") != current_environment:
        raise RuntimeError("runtime environment differs from the frozen manifest")
    for record in manifest["frozen_files"]:
        path = WORKSPACE / record["path"]
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise RuntimeError(f"frozen file changed: {record['path']}")
    return manifest


def validate_translations(path: Path, frozen_sha256: str) -> list[dict]:
    payload = read_json(path)
    if payload.get("status") != "DRAWN_64_TRANSLATIONS_FROZEN_ORDER":
        raise RuntimeError("translation manifest status mismatch")
    if payload.get("frozen_manifest_sha256") != frozen_sha256:
        raise RuntimeError("translation manifest was drawn against another frozen manifest")
    rows = payload.get("translations", [])
    if len(rows) != 64 or [row.get("row_index") for row in rows] != list(range(64)):
        raise RuntimeError("exactly 64 ordered translation rows are required")
    plan = read_json(PLAN_PATH)
    blacklist = {
        (int(tx), int(ty))
        for group in plan["randomization"]["blacklist"].values()
        for tx, ty in group
    }
    seen = set()
    for row in rows:
        if row.get("denominator") != 1048576 or row.get("provenance") != "replication_v04_new_draw":
            raise RuntimeError("translation denominator/provenance mismatch")
        key = (int(row["tx"]), int(row["ty"]))
        if key in blacklist or key in seen:
            raise RuntimeError("historical or duplicate translation in new manifest")
        seen.add(key)
    return rows


def run_geometry_line(geometry, engine, fourier, row: dict, *, second: bool):
    runner = geometry.run_second_quickhull_trajectory if second else geometry.run_primary_trajectory
    trajectories = []
    for radius in (256, 512):
        for p, q in ((1, 0), (8, 1), (1, 1)):
            trajectories.append(
                runner(
                    radius=radius,
                    p=p,
                    q=q,
                    tx=int(row["tx"]),
                    ty=int(row["ty"]),
                    denominator=int(row["denominator"]),
                    engine=engine,
                    fourier=fourier,
                )
            )
    if len(trajectories) != 6 or any(len(item["targets"]) != 3 for item in trajectories):
        raise AssertionError("incomplete geometry line")
    return trajectories


def run_campaign(translation_manifest_path: Path, campaign_dir: Path, frozen_sha256: str) -> dict:
    """Collect all 64 frozen rows; called only by the authorized RUN transition."""
    frozen = verify_frozen_manifest(frozen_sha256)
    rows = validate_translations(translation_manifest_path, frozen_sha256)
    qualification = read_json(QUALIFICATION_PATH)
    if qualification.get("status") != "QUALIFIED_EMPIRICAL_ORACLE_NOT_INTERVAL_CERTIFIED":
        raise RuntimeError("continuous oracle qualification failed")
    if qualification.get("source_sha256") != sha256(ORACLE_PATH):
        raise RuntimeError("oracle source changed after qualification")
    if qualification.get("config_sha256") != sha256(CONFIG_PATH):
        raise RuntimeError("oracle config changed after qualification")
    qualified_environment = qualification.get("environment", {})
    if any(
        qualified_environment.get(name) != current_value
        for name, current_value in (
            ("python", platform.python_version()),
            ("numpy", np.__version__),
            ("scipy", importlib.metadata.version("scipy")),
            ("mpmath", importlib.metadata.version("mpmath")),
            ("threadpoolctl", importlib.metadata.version("threadpoolctl")),
            ("threadpool_limits", 1),
            ("platform", platform.platform()),
            ("psutil", importlib.metadata.version("psutil")),
        )
    ):
        raise RuntimeError("runtime environment differs from oracle qualification")
    if campaign_dir.exists():
        raise FileExistsError("campaign directory already exists; automatic overwrite/resume is forbidden")
    campaign_dir.mkdir(parents=True)
    rows_dir = campaign_dir / "rows"
    rows_dir.mkdir()

    pilot = load_module("v04_campaign_pilot_helpers", PILOT_MODULE_PATH)
    geometry = pilot.geometry
    engine, fourier = geometry.load_historical_modules()
    q_records = []
    row_costs = []
    outputs = []
    second_indices = set(read_json(CONFIG_PATH)["campaign"]["second_engine_row_indices_zero_based"])
    started_campaign = time.perf_counter()

    for row in rows:
        row_index = int(row["row_index"])
        started_wall = time.perf_counter()
        started_cpu = time.process_time()
        primary = run_geometry_line(geometry, engine, fourier, row, second=False)
        primary_wall = time.perf_counter() - started_wall
        primary_cpu = time.process_time() - started_cpu
        primary_path = rows_dir / f"row_{row_index:02d}_main.npz"
        pilot.save_line_npz(primary_path, primary)
        outputs.append(primary_path)

        second_check = None
        second_wall = second_cpu = None
        if row_index in second_indices:
            started_wall = time.perf_counter()
            started_cpu = time.process_time()
            secondary = run_geometry_line(geometry, engine, fourier, row, second=True)
            second_wall = time.perf_counter() - started_wall
            second_cpu = time.process_time() - started_cpu
            second_check = pilot.compare_engines(primary, secondary)
            if second_check["status"] != "PASS_COMPLETE_LAYER_BY_LAYER_SECOND_ENGINE":
                raise RuntimeError(f"second-engine divergence at frozen row {row_index}")
            check_path = rows_dir / f"row_{row_index:02d}_quickhull_check.json"
            write_json_new(check_path, second_check)
            outputs.append(check_path)

        for trajectory in primary:
            orientation = {(1, 0): "theta0", (8, 1): "atan_1_8", (1, 1): "pi_over_4"}[
                (trajectory["p"], trajectory["q"])
            ]
            for target in trajectory["targets"]:
                q_abs = float(target["q_absolute"])
                if not 0.75 <= q_abs <= 1.0:
                    raise RuntimeError("future q_abs outside frozen oracle domain")
                mode3 = geometry.MODES.index(3)
                q_records.append(
                    {
                        "row_index": row_index,
                        "tx": int(row["tx"]),
                        "ty": int(row["ty"]),
                        "denominator": int(row["denominator"]),
                        "provenance": row["provenance"],
                        "R": int(trajectory["R"]),
                        "orientation": orientation,
                        "q_historical": float(target["q_historical"]),
                        "q_absolute": q_abs,
                        "peel_a3_shape": complex(target["a_shape_exp_plus_i_m_theta"][mode3]),
                        "epsilon_interp_a3": float(target["epsilon_interpolation_a3"]),
                    }
                )
        row_costs.append(
            {
                "row_index": row_index,
                "main_wall_seconds": primary_wall,
                "main_cpu_seconds": primary_cpu,
                "second_wall_seconds": second_wall,
                "second_cpu_seconds": second_cpu,
                "main_artifact_sha256": sha256(primary_path),
                "second_check_status": None if second_check is None else second_check["status"],
            }
        )

    if len(q_records) != 64 * 2 * 3 * 3:
        raise AssertionError("finite-grid campaign is not complete")
    q_input_path = campaign_dir / "oracle_q_input_v04.json"
    write_json_new(q_input_path, {"q_absolute": [record["q_absolute"] for record in q_records]})
    outputs.append(q_input_path)
    reference_path = campaign_dir / "campaign_reference_v04.npz"
    oracle = load_module("v04_campaign_oracle", ORACLE_PATH)
    oracle.evaluate(q_input_path, reference_path)
    outputs.extend((reference_path, reference_path.with_suffix(".json")))
    with np.load(reference_path, allow_pickle=False) as archive:
        beta_grid = np.asarray(archive["beta_grid"], dtype=float)
        modes = np.asarray(archive["modes"], dtype=np.int64)
        reference_a = np.asarray(archive["normalized_coefficients"])
        epsilon_ref = np.asarray(archive["epsilon_ref_qualification_grid_max"], dtype=float)
    if beta_grid.tolist() != [0.30, 0.31, 0.32, 1.0 / 3.0, 0.34, 0.35, 0.36]:
        raise AssertionError("oracle beta grid changed")
    mode3 = int(np.where(modes == 3)[0][0])

    cells_path = campaign_dir / "campaign_cells_v04.csv"
    fieldnames = [
        "row_index", "tx", "ty", "denominator", "provenance", "R", "orientation",
        "q_historical", "q_absolute", "beta_label", "beta", "peel_a3_shape_real",
        "peel_a3_shape_imag", "reference_a3_shape_real", "reference_a3_shape_imag",
        "epsilon_interp_a3", "epsilon_ref_a3_qualification_grid_max",
    ]
    beta_labels = ("0.30", "0.31", "0.32", "1/3", "0.34", "0.35", "0.36")
    with cells_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for q_index, record in enumerate(q_records):
            for beta_index, (beta, label) in enumerate(zip(beta_grid, beta_labels)):
                reference = complex(reference_a[beta_index, q_index, mode3])
                peel = record["peel_a3_shape"]
                writer.writerow(
                    {
                        **{key: record[key] for key in (
                            "row_index", "tx", "ty", "denominator", "provenance", "R",
                            "orientation", "q_historical", "q_absolute"
                        )},
                        "beta_label": label,
                        "beta": repr(float(beta)),
                        "peel_a3_shape_real": repr(float(peel.real)),
                        "peel_a3_shape_imag": repr(float(peel.imag)),
                        "reference_a3_shape_real": repr(float(reference.real)),
                        "reference_a3_shape_imag": repr(float(reference.imag)),
                        "epsilon_interp_a3": repr(float(record["epsilon_interp_a3"])),
                        "epsilon_ref_a3_qualification_grid_max": repr(float(epsilon_ref[beta_index])),
                    }
                )
    outputs.append(cells_path)

    costs_path = campaign_dir / "campaign_costs_v04.json"
    write_json_new(costs_path, {"rows": row_costs, "elapsed_wall_seconds": time.perf_counter() - started_campaign})
    outputs.append(costs_path)
    complete_path = campaign_dir / "campaign_complete_manifest_v04.json"
    complete = {
        "status": "COMPLETE_64_ROWS",
        "scientific_collection_executed": True,
        "confirmatory_analysis_executed": False,
        "translation_manifest_path": str(translation_manifest_path.resolve()),
        "translation_manifest_sha256": sha256(translation_manifest_path),
        "frozen_manifest_sha256": frozen_sha256,
        "row_count": 64,
        "finite_grid_cell_count_before_beta_expansion": len(q_records),
        "expanded_cell_count": len(q_records) * 7,
        "second_engine_row_indices_zero_based": sorted(second_indices),
        "geometry_divergence_count": 0,
        "oracle_qualification_sha256": sha256(QUALIFICATION_PATH),
        "cells_path": str(cells_path.resolve()),
        "cells_sha256": sha256(cells_path),
        "outputs": [
            {"path": str(path.resolve()), "sha256": sha256(path), "bytes": path.stat().st_size}
            for path in outputs
        ],
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
    }
    write_json_new(complete_path, complete)
    return complete
