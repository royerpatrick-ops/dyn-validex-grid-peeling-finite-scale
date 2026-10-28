"""Frozen continuous-reference oracle and empirical qualification for v0.4.

The oracle is fixed independently of every peeling result.  It implements
the seven prescribed beta values on q_abs in [0.75, 1].  Qualification uses
separate spatial and time refinements, a distinct support/FD4/Radau method,
disk controls for every beta, and an ellipse control for beta=1/3.

All reported errors are empirical numerical discrepancies, not interval
enclosures.  No exponent, phase, speed, or tolerance is fitted here.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import threading
import time
import tracemalloc
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

try:
    import psutil
except ImportError:  # pragma: no cover
    psutil = None

from threadpoolctl import threadpool_limits


HERE = Path(__file__).resolve().parent
PREP_ROOT = HERE.parent
CONFIG_PATH = PREP_ROOT / "protocol_prepare_v04.json"
QUALIFICATION_PATH = HERE / "reference_qualification_v04.json"
REFERENCE_GRID_PATH = HERE / "reference_grid_v04.npz"
MODES = (0, 3, 6, 9, 12)
RECEIVED_FLOORS = {
    1.0 / 3.0: 2.083506908667765e-9,
    0.32: 3.999023335258609e-13,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


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


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(json_ready(value), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def load_config():
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    oracle = config["continuous_oracle"]
    if oracle["beta_grid"] != [0.30, 0.31, 0.32, 1.0 / 3.0, 0.34, 0.35, 0.36]:
        raise RuntimeError("beta grid changed")
    if oracle["q_absolute_domain"] != [0.75, 1.0]:
        raise RuntimeError("q_abs domain changed")
    if not oracle["no_peeling_feedback"]:
        raise RuntimeError("oracle must be independent of peeling results")
    if oracle.get("threadpool_limits") != 1:
        raise RuntimeError("oracle threadpool limit must be exactly one")
    return config


CONFIG = load_config()
ORACLE = CONFIG["continuous_oracle"]
BETAS = np.asarray(ORACLE["beta_grid"], dtype=float)
T_MAX = float(ORACLE["t_max"])
PRIMARY = ORACLE["primary"]
SPATIAL = ORACLE["spatial_refinement"]
TEMPORAL = ORACLE["temporal_refinement"]
ALTERNATE = ORACLE["distinct_method"]
THRESHOLDS = ORACLE["qualification_thresholds"]


class MemorySampler:
    def __init__(self, interval=0.005):
        self.process = psutil.Process() if psutil else None
        self.interval = interval
        self.baseline = self.process.memory_info().rss if self.process else None
        self.peak = self.baseline
        self.stop_event = threading.Event()
        self.thread = None

    def start(self):
        if self.process is None:
            return

        def sample():
            while not self.stop_event.wait(self.interval):
                try:
                    self.peak = max(self.peak, self.process.memory_info().rss)
                except psutil.Error:
                    return

        self.thread = threading.Thread(target=sample, daemon=True)
        self.thread.start()

    def stop(self):
        if self.process is None:
            return
        self.stop_event.set()
        if self.thread is not None:
            self.thread.join(timeout=1)
        try:
            self.peak = max(self.peak, self.process.memory_info().rss)
        except psutil.Error:
            pass


def full_spectral_operator(n: int) -> np.ndarray:
    frequencies = np.fft.fftfreq(n) * n
    identity = np.eye(n)
    return np.fft.ifft(
        (1 - frequencies[:, None] ** 2) * np.fft.fft(identity, axis=0), axis=0
    ).real


def reduced_operator(n: int, factor: int, method: str):
    if method == "spectral":
        frequencies = np.fft.fftfreq(n) * n
        identity = np.eye(n)
        second = np.fft.ifft(
            -(frequencies[:, None] ** 2) * np.fft.fft(identity, axis=0), axis=0
        ).real
        return identity + factor * factor * second
    if method != "fd4":
        raise ValueError(method)
    dx = 2 * math.pi / n
    rows, columns, data = [], [], []
    for index in range(n):
        for offset, value in ((-2, -1), (-1, 16), (0, -30), (1, 16), (2, -1)):
            rows.append(index)
            columns.append((index + offset) % n)
            data.append(factor * factor * value / (12 * dx * dx))
        rows.append(index)
        columns.append(index)
        data.append(1.0)
    return sparse.csr_matrix((data, (rows, columns)), shape=(n, n))


def solve_rho(beta: float, specification: dict) -> dict:
    n = int(specification["N"])
    phi = 2 * math.pi * np.arange(n) / n
    rho0 = 1 + 0.6 * np.cos(3 * phi)
    multiplier = 1 - (np.fft.fftfreq(n) * n) ** 2
    operator = full_spectral_operator(n)
    minimum = [float("inf")]

    def rhs(_time, rho):
        low = float(np.min(rho))
        minimum[0] = min(minimum[0], low)
        if low <= 0:
            raise ArithmeticError("nonpositive curvature radius")
        return -np.fft.ifft(multiplier * np.fft.fft(rho ** (-beta))).real

    def jac(_time, rho):
        return beta * operator * (rho ** (-beta - 1))[None, :]

    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    solution = solve_ivp(
        rhs,
        (0.0, T_MAX),
        rho0,
        method=specification["time"],
        jac=jac,
        dense_output=True,
        rtol=float(specification["rtol"]),
        atol=float(specification["atol"]),
    )
    elapsed_wall = time.perf_counter() - started_wall
    elapsed_cpu = time.process_time() - started_cpu
    if not solution.success:
        raise ArithmeticError(solution.message)
    return {
        "variable": "rho",
        "beta": beta,
        "N": n,
        "space_method": specification["space"],
        "time_method": specification["time"],
        "rtol": float(specification["rtol"]),
        "atol": float(specification["atol"]),
        "solution": solution,
        "minimum_rho_rhs": minimum[0],
        "elapsed_wall_seconds": elapsed_wall,
        "elapsed_cpu_seconds": elapsed_cpu,
    }


def solve_support(beta: float, specification: dict, *, factor=None, initial=None) -> dict:
    n = int(specification["N"])
    factor = int(specification.get("symmetry_factor", factor))
    coordinate = 2 * math.pi * np.arange(n) / n
    if initial is None:
        initial = lambda phi: 1 - 0.075 * np.cos(3 * phi)
    h0 = np.asarray(initial(coordinate / factor), dtype=float)
    space_method = "fd4" if "FD4" in specification["space"] or "fd4" in specification["space"] else "spectral"
    operator = reduced_operator(n, factor, space_method)
    frequencies = np.fft.rfftfreq(n) * n
    minimum = [float("inf")]

    def rho(values):
        if space_method == "spectral":
            result = values + factor * factor * np.fft.irfft(
                -(frequencies**2) * np.fft.rfft(values), n=n
            )
        else:
            result = operator @ values
        low = float(np.min(result))
        minimum[0] = min(minimum[0], low)
        if low <= 0:
            raise ArithmeticError("nonpositive curvature radius")
        return result

    def rhs(_time, values):
        return -(rho(values) ** (-beta))

    def jac(_time, values):
        weights = beta * rho(values) ** (-beta - 1)
        if sparse.issparse(operator):
            return sparse.diags(weights) @ operator
        return weights[:, None] * operator

    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    solution = solve_ivp(
        rhs,
        (0.0, T_MAX),
        h0,
        method=specification["time"],
        jac=jac,
        dense_output=True,
        rtol=float(specification["rtol"]),
        atol=float(specification["atol"]),
    )
    elapsed_wall = time.perf_counter() - started_wall
    elapsed_cpu = time.process_time() - started_cpu
    if not solution.success:
        raise ArithmeticError(solution.message)
    return {
        "variable": "h",
        "beta": beta,
        "N": n,
        "factor": factor,
        "space_method": space_method,
        "time_method": specification["time"],
        "rtol": float(specification["rtol"]),
        "atol": float(specification["atol"]),
        "solution": solution,
        "initial": h0,
        "minimum_rho_rhs": minimum[0],
        "elapsed_wall_seconds": elapsed_wall,
        "elapsed_cpu_seconds": elapsed_cpu,
    }


def description(run: dict):
    result = {
        key: run[key]
        for key in (
            "variable",
            "beta",
            "N",
            "space_method",
            "time_method",
            "rtol",
            "atol",
            "minimum_rho_rhs",
            "elapsed_wall_seconds",
            "elapsed_cpu_seconds",
        )
    }
    if "factor" in run:
        result["factor"] = run["factor"]
    solution = run["solution"]
    result.update(nfev=solution.nfev, njev=solution.njev, nlu=solution.nlu)
    return result


def target_times(run: dict, q_values: np.ndarray) -> np.ndarray:
    q_values = np.asarray(q_values, dtype=float)
    if np.min(q_values) < 0.75 or np.max(q_values) > 1.0:
        raise ValueError("q_abs outside frozen [0.75,1] domain")
    solution = run["solution"]
    if float(np.mean(solution.y[:, -1])) >= float(np.min(q_values)):
        raise ArithmeticError("reference window does not reach q_abs=0.75")
    times = []
    for q in q_values:
        if q == 1.0:
            times.append(0.0)
        else:
            times.append(
                brentq(
                    lambda t, target=float(q): float(np.mean(solution.sol(t))) - target,
                    0.0,
                    float(solution.t[-1]),
                    xtol=2e-14,
                    rtol=8 * np.finfo(float).eps,
                )
            )
    return np.asarray(times)


def rho_observables(run: dict, q_values: np.ndarray, modes=MODES):
    times = target_times(run, q_values)
    states = run["solution"].sol(times).T
    n = states.shape[1]
    frequencies = np.rint(np.fft.fftfreq(n) * n).astype(int)
    rho_coefficients = np.fft.fft(states, axis=1) / n
    coefficients = np.empty((len(q_values), len(modes)), dtype=complex)
    for mode_index, mode in enumerate(modes):
        index = int(np.flatnonzero(frequencies == mode)[0])
        coefficients[:, mode_index] = rho_coefficients[:, index] / (1 - mode * mode)
    return {
        "times": times,
        "states": states,
        "coefficients": coefficients,
        "a": coefficients / coefficients[:, [0]].real,
        "c0_residual": np.abs(coefficients[:, 0].real - q_values),
        "rho_min": np.min(states, axis=1),
    }


def support_observables(run: dict, q_values: np.ndarray, modes=MODES):
    times = target_times(run, q_values)
    states = run["solution"].sol(times).T
    n = states.shape[1]
    frequencies = np.rint(np.fft.fftfreq(n) * n).astype(int)
    factor = int(run["factor"])
    h_coefficients = np.fft.fft(states, axis=1) / n
    coefficients = np.empty((len(q_values), len(modes)), dtype=complex)
    for mode_index, mode in enumerate(modes):
        if mode % factor:
            coefficients[:, mode_index] = 0
        else:
            index = int(np.flatnonzero(frequencies == mode // factor)[0])
            coefficients[:, mode_index] = h_coefficients[:, index]
    return {
        "times": times,
        "states": states,
        "coefficients": coefficients,
        "a": coefficients / coefficients[:, [0]].real,
        "c0_residual": np.abs(coefficients[:, 0].real - q_values),
    }


def evaluate_rho_support(states: np.ndarray, angles: np.ndarray) -> np.ndarray:
    states = np.atleast_2d(states)
    n = states.shape[1]
    frequencies = np.rint(np.fft.fftfreq(n) * n).astype(int)
    rho_coefficients = np.fft.fft(states, axis=1) / n
    denominator = 1 - frequencies * frequencies
    support_coefficients = np.zeros_like(rho_coefficients)
    mask = denominator != 0
    support_coefficients[:, mask] = rho_coefficients[:, mask] / denominator[mask]
    return (support_coefficients @ np.exp(1j * np.outer(frequencies, angles))).real


def evaluate_reduced_support(states: np.ndarray, factor: int, angles: np.ndarray) -> np.ndarray:
    states = np.atleast_2d(states)
    n = states.shape[1]
    frequencies = np.rint(np.fft.fftfreq(n) * n).astype(int) * factor
    coefficients = np.fft.fft(states, axis=1) / n
    return (coefficients @ np.exp(1j * np.outer(frequencies, angles))).real


def compare_rho(primary: dict, candidate: dict, q_values: np.ndarray):
    left = rho_observables(primary, q_values)
    right = rho_observables(candidate, q_values)
    angles = 2 * math.pi * np.arange(1024) / 1024
    support_left = evaluate_rho_support(left["states"], angles)
    support_right = evaluate_rho_support(right["states"], angles)
    return {
        "support_max_absolute": float(np.max(np.abs(support_left - support_right))),
        "a3_max_absolute": float(np.max(np.abs(left["a"][:, 1] - right["a"][:, 1]))),
        "target_time_max_absolute": float(np.max(np.abs(left["times"] - right["times"]))),
        "target_c0_residual": float(max(np.max(left["c0_residual"]), np.max(right["c0_residual"]))),
        "rho_min": float(min(np.min(left["rho_min"]), np.min(right["rho_min"]), primary["minimum_rho_rhs"], candidate["minimum_rho_rhs"])),
    }


def compare_alternate(primary: dict, alternate: dict, q_values: np.ndarray):
    left = rho_observables(primary, q_values)
    right = support_observables(alternate, q_values)
    angles = 2 * math.pi * np.arange(1024) / 1024
    support_left = evaluate_rho_support(left["states"], angles)
    support_right = evaluate_reduced_support(right["states"], int(alternate["factor"]), angles)
    return {
        "support_max_absolute": float(np.max(np.abs(support_left - support_right))),
        "a3_max_absolute": float(np.max(np.abs(left["a"][:, 1] - right["a"][:, 1]))),
        "target_time_max_absolute": float(np.max(np.abs(left["times"] - right["times"]))),
        "target_c0_residual": float(max(np.max(left["c0_residual"]), np.max(right["c0_residual"]))),
        "rho_min": float(min(primary["minimum_rho_rhs"], alternate["minimum_rho_rhs"], np.min(left["rho_min"]))),
    }


def disk_control(beta: float, q_values: np.ndarray):
    n = 32
    multiplier = 1 - (np.fft.fftfreq(n) * n) ** 2
    operator = full_spectral_operator(n)

    def rhs(_time, rho):
        return -np.fft.ifft(multiplier * np.fft.fft(rho ** (-beta))).real

    def jac(_time, rho):
        return beta * operator * (rho ** (-beta - 1))[None, :]

    solution = solve_ivp(
        rhs,
        (0.0, T_MAX),
        np.ones(n),
        method="BDF",
        jac=jac,
        dense_output=True,
        rtol=float(SPATIAL["rtol"]),
        atol=float(SPATIAL["atol"]),
    )
    if not solution.success:
        raise ArithmeticError(solution.message)
    times = []
    for q in q_values:
        times.append(0.0 if q == 1.0 else brentq(lambda t, target=q: float(np.mean(solution.sol(t))) - target, 0, T_MAX))
    times = np.asarray(times)
    exact = (1 - q_values ** (beta + 1)) / (beta + 1)
    states = solution.sol(times).T
    return {
        "time_max_absolute": float(np.max(np.abs(times - exact))),
        "state_max_absolute": float(np.max(np.abs(states - q_values[:, None]))),
    }


def ellipse_control(q_values: np.ndarray):
    initial = lambda phi: np.sqrt(2 * np.cos(phi) ** 2 + 0.5 * np.sin(phi) ** 2)
    specification = {
        "N": 128,
        "space": "spectral",
        "time": "BDF",
        "rtol": float(SPATIAL["rtol"]),
        "atol": float(SPATIAL["atol"]),
    }
    run = solve_support(1 / 3, specification, factor=2, initial=initial)
    initial_c0 = float(np.mean(run["initial"]))
    targets = q_values * initial_c0
    # Generic target_times enforces the nonellipse [0.75,1] scale.  Ellipse
    # targets are absolute supports, so evaluate roots directly here.
    solution = run["solution"]
    times = np.asarray(
        [
            0.0
            if q == 1.0
            else brentq(lambda t, target=float(target): float(np.mean(solution.sol(t))) - target, 0, T_MAX)
            for q, target in zip(q_values, targets)
        ]
    )
    states = solution.sol(times).T
    exact_times = 0.75 * (1 - q_values ** (4 / 3))
    exact_states = q_values[:, None] * run["initial"][None, :]
    return {
        "time_max_absolute": float(np.max(np.abs(times - exact_times))),
        "support_max_absolute": float(np.max(np.abs(states - exact_states))),
        "description": description(run),
    }


def solve_primary_for_targets(q_values: np.ndarray):
    """Production rule: fixed primary solver, no access to peeling losses."""
    q_values = np.asarray(q_values, dtype=float)
    if np.min(q_values) < 0.75 or np.max(q_values) > 1.0:
        raise ValueError("q_abs outside frozen oracle domain")
    values = {}
    with threadpool_limits(limits=1):
        for beta in BETAS:
            run = solve_rho(float(beta), PRIMARY)
            values[float(beta)] = rho_observables(run, q_values)
    return values


def qualify():
    q_values = np.linspace(0.75, 1.0, int(ORACLE["qualification_grid_count"]))
    sampler = MemorySampler()
    tracemalloc.start()
    sampler.start()
    total_wall_start = time.perf_counter()
    total_cpu_start = time.process_time()
    by_beta = {}
    reference_values = {}
    disk_controls = {}
    comparisons = []
    with threadpool_limits(limits=1):
        for beta_value in BETAS:
            beta = float(beta_value)
            primary = solve_rho(beta, PRIMARY)
            spatial = solve_rho(beta, SPATIAL)
            temporal = solve_rho(beta, TEMPORAL)
            alternate = solve_support(beta, ALTERNATE)
            primary_observed = rho_observables(primary, q_values)
            # Isolate space from time: N=192-tight versus N=384-tight.
            # The separate time check remains N=192-primary versus N=192-tight.
            spatial_check = compare_rho(temporal, spatial, q_values)
            temporal_check = compare_rho(primary, temporal, q_values)
            alternate_check = compare_alternate(primary, alternate, q_values)
            comparisons.extend((spatial_check, temporal_check, alternate_check))
            empirical = max(
                spatial_check["a3_max_absolute"],
                temporal_check["a3_max_absolute"],
                alternate_check["a3_max_absolute"],
            )
            received_floor = RECEIVED_FLOORS.get(beta, 0.0)
            qualification_grid_epsilon = max(empirical, received_floor)
            by_beta[f"{beta:.17g}"] = {
                "primary": description(primary),
                "spatial_refinement": {
                    "baseline_description": description(temporal),
                    "description": description(spatial),
                    "comparison_N192_tight_to_N384_tight": spatial_check,
                },
                "time_refinement": {"description": description(temporal), "comparison_to_primary": temporal_check},
                "distinct_support_FD4_Radau": {"description": description(alternate), "comparison_to_primary": alternate_check},
                "empirical_a3_reference_error_max_on_qualification_grid": empirical,
                "received_v03_floor": received_floor,
                "epsilon_ref_qualification_grid_max_used": qualification_grid_epsilon,
                "epsilon_rule": "max observed on the finite qualification grid across space, time and distinct-method a3 comparisons, then max with the received v0.3 floor when prescribed",
                "primary_target_c0_residual": float(np.max(primary_observed["c0_residual"])),
                "primary_rho_min_on_q_grid": float(np.min(primary_observed["rho_min"])),
            }
            reference_values[beta] = primary_observed
            disk_controls[f"{beta:.17g}"] = disk_control(beta, q_values)
        ellipse = ellipse_control(q_values)
    total_cpu = time.process_time() - total_cpu_start
    total_wall = time.perf_counter() - total_wall_start
    trace_current, trace_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    sampler.stop()

    maxima = {
        key: float(max(check[key] for check in comparisons))
        for key in ("support_max_absolute", "a3_max_absolute", "target_time_max_absolute", "target_c0_residual")
    }
    analytic_max = max(
        [max(control.values()) for control in disk_controls.values()]
        + [ellipse["time_max_absolute"], ellipse["support_max_absolute"]]
    )
    positivity = min(check["rho_min"] for check in comparisons) > 0
    qualified = (
        positivity
        and maxima["support_max_absolute"] <= float(THRESHOLDS["support_max_absolute"])
        and maxima["a3_max_absolute"] <= float(THRESHOLDS["a3_max_absolute"])
        and maxima["target_time_max_absolute"] <= float(THRESHOLDS["target_time_max_absolute"])
        and maxima["target_c0_residual"] <= float(THRESHOLDS["target_c0_residual"])
        and analytic_max <= float(THRESHOLDS["analytic_control_max_absolute"])
        and max(
            row["epsilon_ref_qualification_grid_max_used"] for row in by_beta.values()
        ) <= float(THRESHOLDS["a3_max_absolute"])
    )
    np.savez_compressed(
        REFERENCE_GRID_PATH,
        beta_grid=BETAS,
        q_absolute=q_values,
        modes=np.asarray(MODES, dtype=np.int64),
        times=np.vstack([reference_values[float(beta)]["times"] for beta in BETAS]),
        coefficients=np.stack([reference_values[float(beta)]["coefficients"] for beta in BETAS]),
        normalized_coefficients=np.stack([reference_values[float(beta)]["a"] for beta in BETAS]),
        epsilon_ref_qualification_grid_max=np.asarray(
            [by_beta[f"{float(beta):.17g}"]["epsilon_ref_qualification_grid_max_used"] for beta in BETAS]
        ),
    )
    primary_wall = sum(by_beta[key]["primary"]["elapsed_wall_seconds"] for key in by_beta)
    primary_cpu = sum(by_beta[key]["primary"]["elapsed_cpu_seconds"] for key in by_beta)
    result = {
        "status": "QUALIFIED_EMPIRICAL_ORACLE_NOT_INTERVAL_CERTIFIED" if qualified else "ORACLE_QUALIFICATION_FAILED",
        "source_sha256": sha256(Path(__file__).resolve()),
        "config_sha256": sha256(CONFIG_PATH),
        "domain": {"q_abs_min": 0.75, "q_abs_max": 1.0, "grid_count": len(q_values), "grid_spacing": float(q_values[1] - q_values[0])},
        "beta_grid": BETAS.tolist(),
        "no_peeling_result_used": True,
        "parameters": ORACLE,
        "by_beta": by_beta,
        "epsilon_ref_rule": (
            "per beta maximum observed on the finite qualification grid across separate space refinement, "
            "separate time refinement and distinct support formulation; beta=1/3 and .32 also use max with received v0.3 floors; not a domain envelope"
        ),
        "received_v03_floors": {"1/3": RECEIVED_FLOORS[1 / 3], "0.32": RECEIVED_FLOORS[0.32]},
        "epsilon_ref_qualification_grid_max_final": {
            f"{float(beta):.17g}": by_beta[f"{float(beta):.17g}"]["epsilon_ref_qualification_grid_max_used"]
            for beta in BETAS
        },
        "disk_controls_every_beta": disk_controls,
        "ellipse_beta_1_over_3": ellipse,
        "max_refinement_discrepancies": maxima,
        "max_analytic_control_discrepancy": analytic_max,
        "positivity_all_beta_and_methods": positivity,
        "thresholds": THRESHOLDS,
        "cost": {
            "qualification_wall_seconds": total_wall,
            "qualification_cpu_seconds": total_cpu,
            "production_primary_seven_beta_wall_seconds_observed_within_qualification": primary_wall,
            "production_primary_seven_beta_cpu_seconds_observed_within_qualification": primary_cpu,
            "peak_rss_bytes": sampler.peak,
            "baseline_rss_bytes": sampler.baseline,
            "incremental_peak_rss_bytes": None if sampler.peak is None else max(0, sampler.peak - sampler.baseline),
            "peak_tracemalloc_bytes": trace_peak,
            "tracemalloc_current_at_stop_bytes": trace_current,
            "reference_grid_storage_bytes": REFERENCE_GRID_PATH.stat().st_size,
        },
        "reference_grid": {"path": str(REFERENCE_GRID_PATH), "sha256": sha256(REFERENCE_GRID_PATH), "bytes": REFERENCE_GRID_PATH.stat().st_size},
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "mpmath": importlib.metadata.version("mpmath"),
            "threadpoolctl": importlib.metadata.version("threadpoolctl"),
            "threadpool_limits": 1,
            "platform": platform.platform(),
            "psutil": None if psutil is None else psutil.__version__,
        },
        "scope": "finite 1025-node-grid empirical qualification over the frozen q_abs domain; neither interval enclosure nor domain-wide envelope",
    }
    write_json(QUALIFICATION_PATH, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "qualification_wall_seconds": total_wall,
                "max_a3_discrepancy": maxima["a3_max_absolute"],
                "minimum_rho_positive": positivity,
                "epsilon_ref_qualification_grid_max_final": result["epsilon_ref_qualification_grid_max_final"],
            }
        ),
        flush=True,
    )
    if not qualified:
        raise SystemExit(2)


def evaluate(input_path: Path, output_path: Path):
    qualification = json.loads(QUALIFICATION_PATH.read_text(encoding="utf-8"))
    if qualification["status"] != "QUALIFIED_EMPIRICAL_ORACLE_NOT_INTERVAL_CERTIFIED":
        raise RuntimeError("frozen oracle is not qualified")
    if qualification["source_sha256"] != sha256(Path(__file__).resolve()):
        raise RuntimeError("oracle code changed after qualification")
    if qualification["config_sha256"] != sha256(CONFIG_PATH):
        raise RuntimeError("oracle protocol config changed after qualification")
    qualified_environment = qualification["environment"]
    current_environment = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "mpmath": importlib.metadata.version("mpmath"),
        "threadpoolctl": importlib.metadata.version("threadpoolctl"),
        "threadpool_limits": 1,
        "platform": platform.platform(),
        "psutil": None if psutil is None else psutil.__version__,
    }
    if current_environment != qualified_environment:
        raise RuntimeError("runtime environment differs from the qualified oracle environment")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    q_values = np.asarray(payload["q_absolute"], dtype=float)
    if q_values.ndim != 1 or len(q_values) == 0 or not np.all(np.isfinite(q_values)):
        raise ValueError("q_absolute must be a nonempty finite one-dimensional array")
    q_min, q_max = map(float, ORACLE["q_absolute_domain"])
    if float(np.min(q_values)) < q_min or float(np.max(q_values)) > q_max:
        raise ValueError("q_absolute is outside the frozen qualification domain")
    if output_path.exists() or output_path.with_suffix(".json").exists():
        raise FileExistsError("refusing to overwrite an oracle evaluation")
    values = solve_primary_for_targets(q_values)
    np.savez_compressed(
        output_path,
        q_absolute=q_values,
        beta_grid=BETAS,
        modes=np.asarray(MODES, dtype=np.int64),
        coefficients=np.stack([values[float(beta)]["coefficients"] for beta in BETAS]),
        normalized_coefficients=np.stack([values[float(beta)]["a"] for beta in BETAS]),
        epsilon_ref_qualification_grid_max=np.asarray(
            [qualification["epsilon_ref_qualification_grid_max_final"][f"{float(beta):.17g}"] for beta in BETAS]
        ),
    )
    metadata = {
        "status": "EVALUATED_BY_FROZEN_ORACLE",
        "q_count": len(q_values),
        "q_min": float(np.min(q_values)),
        "q_max": float(np.max(q_values)),
        "input_sha256": sha256(input_path),
        "output_sha256": sha256(output_path),
        "oracle_source_sha256": sha256(Path(__file__).resolve()),
        "qualification_sha256": sha256(QUALIFICATION_PATH),
        "no_parameter_fit": True,
    }
    write_json(output_path.with_suffix(".json"), metadata)
    print(json.dumps(metadata))


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--qualify", action="store_true")
    group.add_argument("--evaluate", action="store_true")
    parser.add_argument("--q-input", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.qualify:
        qualify()
    else:
        if args.q_input is None or args.output is None:
            parser.error("--evaluate requires --q-input and --output")
        evaluate(args.q_input, args.output)


if __name__ == "__main__":
    main()
