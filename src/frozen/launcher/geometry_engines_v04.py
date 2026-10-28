"""Frozen finite-grid geometry used by the v0.4 preparation.

The primary engine is loaded read-only from the qualified P06 source after
checking its SHA-256.  The independent hull route is an exact-integer
QuickHull over row endpoints.  It has an independent layer driver and
vertex-removal routine, but deliberately shares the qualified initial row
intervals; this limitation is reported rather than hidden.

No random number generator is used in this module.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
from pathlib import Path

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[2]
ENGINE_PATH = (
    WORKSPACE
    / "_audit_codex_pi_morphodynamique_20261027"
    / "working_copies"
    / "P06"
    / "DYN-Transfert_NonElliptique_Horloge-Gelee_v0_1"
    / "support_shape_engine.py"
)
FOURIER_PATH = (
    WORKSPACE
    / "_audit_codex_pi_morphodynamique_20261027"
    / "working_copies"
    / "P07"
    / "DYN-Comparaison_Contraction-Identique_v0_1"
    / "baseline"
    / "polygon_fourier.py"
)
EXPECTED_ENGINE_SHA256 = "f1769af06d914f5ff078a3c22c5c259b8777656d5c3596321362766a6e51a191"
EXPECTED_FOURIER_SHA256 = "ec5bcec9b6e5c52bc136f39ee9c27a82a2fbc7c03b2493b24de321c3b4446245"
MODES = (0, 1, 3, 6, 9, 12)
TARGETS = (0.95, 0.90, 0.85)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_historical_modules():
    if sha256(ENGINE_PATH) != EXPECTED_ENGINE_SHA256:
        raise RuntimeError("Historical geometry engine SHA-256 mismatch")
    if sha256(FOURIER_PATH) != EXPECTED_FOURIER_SHA256:
        raise RuntimeError("Historical polygon-Fourier SHA-256 mismatch")
    return (
        _load_module("v04_historical_geometry", ENGINE_PATH),
        _load_module("v04_historical_fourier", FOURIER_PATH),
    )


def row_endpoints(rows: dict[int, list[int]]) -> np.ndarray:
    points = []
    for y, (lo, hi) in rows.items():
        points.append((int(lo), int(y)))
        if hi != lo:
            points.append((int(hi), int(y)))
    return np.asarray(points, dtype=np.int64)


def canonical_vertices(vertices) -> np.ndarray:
    points = np.asarray(vertices, dtype=np.int64)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("vertices must be an Nx2 array")
    points = np.unique(points, axis=0)
    order = np.lexsort((points[:, 1], points[:, 0]))
    return np.ascontiguousarray(points[order].astype("<i8", copy=False))


def vertex_digest(vertices) -> str:
    return hashlib.sha256(canonical_vertices(vertices).tobytes()).hexdigest()


def cross_int(a, b, c) -> int:
    return (int(b[0]) - int(a[0])) * (int(c[1]) - int(a[1])) - (
        int(b[1]) - int(a[1])
    ) * (int(c[0]) - int(a[0]))


def _cross_values(a: np.ndarray, b: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Exact int64 cross products after a proved no-overflow guard."""
    points = np.asarray(points, dtype=np.int64)
    if len(points) == 0:
        return np.empty(0, dtype=np.int64)
    # Every bound is formed as a Python arbitrary-precision integer before
    # any NumPy subtraction; opposite int64 extremes therefore cannot wrap.
    ax, ay = int(a[0]), int(a[1])
    bx, by = int(b[0]), int(b[1])
    dx, dy = bx - ax, by - ay
    min_x, max_x = int(np.min(points[:, 0])), int(np.max(points[:, 0]))
    min_y, max_y = int(np.min(points[:, 1])), int(np.max(points[:, 1]))
    span_x = max(abs(dx), abs(min_x - ax), abs(max_x - ax))
    span_y = max(abs(dy), abs(min_y - ay), abs(max_y - ay))
    worst = 2 * span_x * span_y
    if worst > np.iinfo(np.int64).max:
        raise OverflowError("int64 cross-product safety guard failed")
    # The proved spans fit int64, so these subtractions and products are exact.
    return np.int64(dx) * (points[:, 1] - np.int64(ay)) - np.int64(dy) * (
        points[:, 0] - np.int64(ax)
    )


def _signed_twice_area(vertices: np.ndarray) -> int:
    total = 0
    for a, b in zip(vertices, np.vstack((vertices[1:], vertices[:1]))):
        total += int(a[0]) * int(b[1]) - int(a[1]) * int(b[0])
    return total


def quickhull_integer_endpoints(points) -> list[tuple[int, int]]:
    """Exact QuickHull using integer cross products; collinear points drop.

    For a row-convex lattice set every point strictly between the two row
    endpoints is a convex combination of those endpoints and cannot be a
    strict extreme point.  Therefore the endpoint candidate set has exactly
    the same extreme vertices as exhaustive all-point hull construction.
    """
    pts = canonical_vertices(points).astype(np.int64, copy=False)
    if len(pts) < 3:
        raise ValueError("Degenerate point set")
    a = pts[0]
    b = pts[-1]
    if np.array_equal(a, b):
        raise ValueError("Degenerate point set")

    def chain(left: np.ndarray, right: np.ndarray, candidates: np.ndarray):
        if len(candidates) == 0:
            return [tuple(map(int, right))]
        distances = _cross_values(left, right, candidates)
        positive = distances > 0
        candidates = candidates[positive]
        distances = distances[positive]
        if len(candidates) == 0:
            return [tuple(map(int, right))]
        maximum = int(np.max(distances))
        farthest = candidates[distances == maximum]
        # Deterministic tie rule; all tied points are exact arithmetic ties.
        order = np.lexsort((farthest[:, 1], farthest[:, 0]))
        pivot = farthest[order[0]]
        left_set = candidates[_cross_values(left, pivot, candidates) > 0]
        right_set = candidates[_cross_values(pivot, right, candidates) > 0]
        return chain(left, pivot, left_set) + chain(pivot, right, right_set)

    upper = pts[_cross_values(a, b, pts) > 0]
    lower = pts[_cross_values(b, a, pts) > 0]
    ordered = [tuple(map(int, a))] + chain(a, b, upper) + chain(b, a, lower)
    if ordered[-1] == ordered[0]:
        ordered.pop()
    hull = np.asarray(ordered, dtype=np.int64)
    if len(hull) < 3 or _signed_twice_area(hull) == 0:
        raise ValueError("Collinear or degenerate hull")
    if _signed_twice_area(hull) < 0:
        hull = hull[::-1]
    if any(cross_int(hull[i - 1], hull[i], hull[(i + 1) % len(hull)]) <= 0 for i in range(len(hull))):
        raise AssertionError("QuickHull retained a collinear point or wrong turn")
    return [tuple(map(int, point)) for point in hull]


def remove_vertices_independent(rows: dict[int, list[int]], vertices) -> None:
    """Independent mutation: remove endpoints that are strict hull vertices."""
    for x_value, y_value in vertices:
        x, y = int(x_value), int(y_value)
        lo, hi = rows[y]
        if x == lo:
            lo += 1
        elif x == hi:
            hi -= 1
        else:
            raise AssertionError("QuickHull vertex is not a row endpoint")
        if lo > hi:
            del rows[y]
        else:
            rows[y] = [lo, hi]


def perimeter(vertices) -> float:
    values = list(vertices)
    return math.fsum(
        math.hypot(int(a[0]) - int(b[0]), int(a[1]) - int(b[1]))
        for a, b in zip(values, values[1:] + values[:1])
    )


def row_point_count(rows: dict[int, list[int]]) -> int:
    return sum(int(hi) - int(lo) + 1 for lo, hi in rows.values())


def exact_containment_min_cross(outer, inner) -> int:
    outside = np.asarray(outer, dtype=np.int64)
    inside = np.asarray(inner, dtype=np.int64)
    if len(outside) < 3 or len(inside) < 3:
        raise ValueError("degenerate containment input")
    if _signed_twice_area(outside) <= 0 or _signed_twice_area(inside) <= 0:
        raise AssertionError("containment requires strict CCW polygons")
    min_x, max_x = int(np.min(inside[:, 0])), int(np.max(inside[:, 0]))
    min_y, max_y = int(np.min(inside[:, 1])), int(np.max(inside[:, 1]))
    global_min = None
    for index, a in enumerate(outside):
        b = outside[(index + 1) % len(outside)]
        ax, ay = int(a[0]), int(a[1])
        dx, dy = int(b[0]) - ax, int(b[1]) - ay
        span_x = max(abs(dx), abs(min_x - ax), abs(max_x - ax))
        span_y = max(abs(dy), abs(min_y - ay), abs(max_y - ay))
        if 2 * span_x * span_y > np.iinfo(np.int64).max:
            raise OverflowError("containment cross-product safety guard failed")
        values = np.int64(dx) * (inside[:, 1] - np.int64(ay)) - np.int64(dy) * (
            inside[:, 0] - np.int64(ax)
        )
        local = int(np.min(values))
        global_min = local if global_min is None else min(global_min, local)
    return int(global_min)


def qualify_quickhull() -> dict:
    tests = {}
    square_points = np.asarray([(x, y) for y in range(-1, 2) for x in range(-1, 2)], dtype=np.int64)
    square_hull = quickhull_integer_endpoints(square_points)
    expected = {(-1, -1), (1, -1), (1, 1), (-1, 1)}
    if set(square_hull) != expected:
        raise AssertionError("3x3 square hull mismatch")
    square_rows = {-1: [-1, 1], 0: [-1, 1], 1: [-1, 1]}
    remove_vertices_independent(square_rows, square_hull)
    remaining = {
        (x, y)
        for y, (lo, hi) in square_rows.items()
        for x in range(int(lo), int(hi) + 1)
    }
    expected_remaining = {(0, 0), (-1, 0), (1, 0), (0, -1), (0, 1)}
    if remaining != expected_remaining:
        raise AssertionError("aligned non-vertices were removed")
    tests["square_3x3"] = {
        "hull": sorted(map(list, expected)),
        "remaining_after_one_layer": sorted(map(list, remaining)),
        "aligned_nonvertices_retained": True,
    }
    rejected = 0
    for points in (
        [(0, 0)],
        [(0, 0), (1, 0)],
        [(0, 0), (1, 0), (2, 0), (3, 0)],
    ):
        try:
            quickhull_integer_endpoints(points)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("degenerate input accepted")
    tests["degenerate_inputs_rejected"] = rejected
    # A large-coordinate test exercises the explicit overflow guard while
    # remaining well inside signed int64 arithmetic.
    large = np.asarray([(-10**8, -10**8), (10**8, -10**8), (10**8, 10**8), (-10**8, 10**8), (0, 0)], dtype=np.int64)
    large_hull = quickhull_integer_endpoints(large)
    if len(large_hull) != 4:
        raise AssertionError("large-coordinate hull mismatch")
    tests["large_integer_predicates"] = {
        "maximum_coordinate_absolute": 10**8,
        "hull_vertices": len(large_hull),
        "overflow_guard": "passed",
    }
    try:
        _cross_values(
            np.asarray([np.iinfo(np.int64).min, 0], dtype=np.int64),
            np.asarray([np.iinfo(np.int64).max, 1], dtype=np.int64),
            np.asarray([[0, np.iinfo(np.int64).max]], dtype=np.int64),
        )
    except OverflowError:
        tests["overflow_rejection_at_int64_extremes"] = True
    else:
        raise AssertionError("overflow guard failed to reject an unsafe predicate")
    return {
        "status": "PASS_EXACT_INTEGER_QUICKHULL_ENDPOINTS",
        "collinear_candidates_eliminated": True,
        "row_convex_endpoint_equivalence": (
            "each omitted row-interior point is a convex combination of its two retained row endpoints"
        ),
        "initial_membership_independence_limit": (
            "QuickHull shares the qualified P06 initial row intervals; only hull construction, layer driver and removal are independent"
        ),
        "tests": tests,
    }


def _assert_strict_ccw(vertices) -> np.ndarray:
    ordered = np.asarray(vertices, dtype=np.int64)
    if len(ordered) < 3 or _signed_twice_area(ordered) <= 0:
        raise AssertionError("polygon is not nondegenerate CCW")
    if any(
        cross_int(ordered[index - 1], ordered[index], ordered[(index + 1) % len(ordered)]) <= 0
        for index in range(len(ordered))
    ):
        raise AssertionError("polygon is not strictly convex CCW")
    return ordered


def _target_record(
    *, target_index, target, n, p0, radius, p, q, tx, ty, denominator,
    previous_vertices, vertices, previous_perimeter, current_perimeter,
    previous_count, count, fourier,
):
    before_ordered = _assert_strict_ccw(previous_vertices)
    after_ordered = _assert_strict_ccw(vertices)
    if not previous_perimeter > target * p0 >= current_perimeter:
        raise AssertionError("target bracket mismatch")
    containment = exact_containment_min_cross(before_ordered, after_ordered)
    if containment < 0:
        raise AssertionError("adjacent target polygons are not nested")
    weight = (previous_perimeter - target * p0) / (previous_perimeter - current_perimeter)
    if not 0 <= weight <= 1:
        raise AssertionError("Minkowski weight outside [0,1]")
    translation = (tx / denominator, ty / denominator)
    before_map = fourier.coefficients(
        [tuple(map(int, value)) for value in before_ordered],
        translation=translation,
        scale=radius,
        modes=MODES,
    )
    after_map = fourier.coefficients(
        [tuple(map(int, value)) for value in after_ordered],
        translation=translation,
        scale=radius,
        modes=MODES,
    )
    before_raw = np.asarray([before_map[mode] for mode in MODES], dtype=np.complex128)
    after_raw = np.asarray([after_map[mode] for mode in MODES], dtype=np.complex128)
    before_expected_c0 = previous_perimeter / (2 * math.pi * radius)
    after_expected_c0 = current_perimeter / (2 * math.pi * radius)
    before_c0_complex_residual = abs(before_raw[0] - before_expected_c0)
    after_c0_complex_residual = abs(after_raw[0] - after_expected_c0)
    if max(before_c0_complex_residual, after_c0_complex_residual) > 5e-14:
        raise AssertionError("raw endpoint c0/perimeter identity failed")
    interpolated_raw = (1 - weight) * before_raw + weight * after_raw
    target_c0 = target * p0 / (2 * math.pi * radius)
    interpolated_c0_complex_residual = abs(interpolated_raw[0] - target_c0)
    if interpolated_c0_complex_residual > 5e-14:
        raise AssertionError("raw Minkowski c0 target mismatch")
    interpolated_corrected = interpolated_raw.copy()
    interpolated_corrected[0] = complex(target_c0, 0.0)
    a_global = interpolated_corrected / target_c0
    theta = math.atan2(q, p)
    a_shape = a_global * np.exp(1j * np.asarray(MODES) * theta)
    interpolation_coefficient_bound = 2 * weight * (1 - weight) * (
        previous_perimeter - current_perimeter
    ) / (2 * math.pi * radius)
    return {
        "target_index": target_index,
        "q_historical": target,
        "n_before": n - 1,
        "n_after": n,
        "perimeter_before": previous_perimeter,
        "perimeter_after": current_perimeter,
        "remaining_before": previous_count,
        "remaining_after": count,
        "weight_after": weight,
        "q_absolute": target_c0,
        "before_c0_complex_residual": before_c0_complex_residual,
        "after_c0_complex_residual": after_c0_complex_residual,
        "interpolated_c0_complex_residual_before_correction": interpolated_c0_complex_residual,
        "minimum_exact_containment_cross": containment,
        "before_vertices_ordered": before_ordered,
        "after_vertices_ordered": after_ordered,
        "before_vertices_canonical": canonical_vertices(before_ordered),
        "after_vertices_canonical": canonical_vertices(after_ordered),
        "before_vertex_sha256": vertex_digest(before_ordered),
        "after_vertex_sha256": vertex_digest(after_ordered),
        "before_coefficients_raw": before_raw,
        "after_coefficients_raw": after_raw,
        "interpolated_coefficients_raw": interpolated_raw,
        "interpolated_coefficients_corrected_c0": interpolated_corrected,
        "a_global": a_global,
        "a_shape_exp_plus_i_m_theta": a_shape,
        "interpolation_coefficient_bound": interpolation_coefficient_bound,
        "epsilon_interpolation_a3": interpolation_coefficient_bound / target_c0,
    }


def run_primary_trajectory(
    *, radius: int, p: int, q: int, tx: int, ty: int, denominator: int,
    engine, fourier, targets=TARGETS,
):
    """Historical main driver: Andrew hull and historical removal routine."""
    if tuple(map(float, targets)) != TARGETS:
        raise ValueError("target allowlist mismatch")
    rows = engine.initial_rows(radius, p, q, tx, ty, denominator)
    initial_count = row_point_count(rows)
    vertices = engine.hull(engine.endpoints(rows))
    initial_ordered = _assert_strict_ccw(vertices)
    p0 = perimeter(vertices)
    layers, target_records = [], []
    previous_vertices = previous_perimeter = previous_count = None
    global_containment = None
    n = removed_count = 0
    while True:
        ordered = _assert_strict_ccw(vertices)
        count = row_point_count(rows)
        if count + removed_count != initial_count:
            raise AssertionError("primary point-count conservation failed")
        current_perimeter = perimeter(vertices)
        containment = None
        if n:
            if not current_perimeter < previous_perimeter:
                raise AssertionError("primary perimeter failed strict decrease")
            containment = exact_containment_min_cross(previous_vertices, ordered)
            if containment < 0:
                raise AssertionError("primary adjacent layers are not nested")
            global_containment = containment if global_containment is None else min(global_containment, containment)
        layers.append({
            "n": n,
            "remaining_points": count,
            "vertex_count": len(ordered),
            "perimeter": current_perimeter,
            "vertex_sha256": vertex_digest(ordered),
            "vertices_ordered": ordered.copy(),
            "minimum_containment_cross_from_previous": containment,
        })
        while len(target_records) < len(targets) and current_perimeter <= targets[len(target_records)] * p0:
            if n == 0 or previous_vertices is None:
                raise AssertionError("target crossed at primary initial layer")
            target_index = len(target_records)
            target_records.append(_target_record(
                target_index=target_index, target=float(targets[target_index]), n=n,
                p0=p0, radius=radius, p=p, q=q, tx=tx, ty=ty, denominator=denominator,
                previous_vertices=previous_vertices, vertices=ordered,
                previous_perimeter=previous_perimeter, current_perimeter=current_perimeter,
                previous_count=previous_count, count=count, fourier=fourier,
            ))
        if len(target_records) == len(targets):
            break
        previous_vertices = ordered.copy()
        previous_perimeter, previous_count = current_perimeter, count
        removed_count += len(vertices)
        engine.remove_vertices(rows, vertices)
        vertices = engine.hull(engine.endpoints(rows))
        n += 1
        if n > 10000:
            raise RuntimeError("primary guard: q=.85 not reached")
    return {
        "R": int(radius), "p": int(p), "q": int(q), "tx": int(tx), "ty": int(ty),
        "denominator": int(denominator), "method": "main_monotone_historical_driver",
        "initial_point_count": initial_count, "initial_vertices_ordered": initial_ordered,
        "initial_vertex_sha256": vertex_digest(initial_ordered), "initial_perimeter": p0,
        "minimum_all_adjacent_layer_containment_cross": global_containment,
        "layers": layers, "targets": target_records,
    }


def run_second_quickhull_trajectory(
    *, radius: int, p: int, q: int, tx: int, ty: int, denominator: int,
    engine, fourier, targets=TARGETS,
):
    """Separate QuickHull driver; shares only qualified initial row intervals."""
    if tuple(map(float, targets)) != TARGETS:
        raise ValueError("target allowlist mismatch")
    rows = {int(y): [int(lo), int(hi)] for y, (lo, hi) in engine.initial_rows(radius, p, q, tx, ty, denominator).items()}
    initial_count = row_point_count(rows)
    vertices = quickhull_integer_endpoints(row_endpoints(rows))
    initial_ordered = _assert_strict_ccw(vertices)
    p0 = perimeter(vertices)
    layers, target_records = [], []
    previous_vertices = previous_perimeter = previous_count = None
    global_containment = None
    n = removed_count = 0
    while True:
        ordered = _assert_strict_ccw(vertices)
        count = row_point_count(rows)
        if count + removed_count != initial_count:
            raise AssertionError("secondary point-count conservation failed")
        current_perimeter = perimeter(vertices)
        containment = None
        if n:
            if not current_perimeter < previous_perimeter:
                raise AssertionError("secondary perimeter failed strict decrease")
            containment = exact_containment_min_cross(previous_vertices, ordered)
            if containment < 0:
                raise AssertionError("secondary adjacent layers are not nested")
            global_containment = containment if global_containment is None else min(global_containment, containment)
        layers.append({
            "n": n,
            "remaining_points": count,
            "vertex_count": len(ordered),
            "perimeter": current_perimeter,
            "vertex_sha256": vertex_digest(ordered),
            "vertices_ordered": ordered.copy(),
            "minimum_containment_cross_from_previous": containment,
        })
        while len(target_records) < len(targets) and current_perimeter <= targets[len(target_records)] * p0:
            if n == 0 or previous_vertices is None:
                raise AssertionError("target crossed at secondary initial layer")
            target_index = len(target_records)
            target_records.append(_target_record(
                target_index=target_index, target=float(targets[target_index]), n=n,
                p0=p0, radius=radius, p=p, q=q, tx=tx, ty=ty, denominator=denominator,
                previous_vertices=previous_vertices, vertices=ordered,
                previous_perimeter=previous_perimeter, current_perimeter=current_perimeter,
                previous_count=previous_count, count=count, fourier=fourier,
            ))
        if len(target_records) == len(targets):
            break
        previous_vertices = ordered.copy()
        previous_perimeter, previous_count = current_perimeter, count
        removed_count += len(vertices)
        remove_vertices_independent(rows, vertices)
        vertices = quickhull_integer_endpoints(row_endpoints(rows))
        n += 1
        if n > 10000:
            raise RuntimeError("secondary guard: q=.85 not reached")
    return {
        "R": int(radius), "p": int(p), "q": int(q), "tx": int(tx), "ty": int(ty),
        "denominator": int(denominator), "method": "second_quickhull_independent_driver",
        "initial_point_count": initial_count, "initial_vertices_ordered": initial_ordered,
        "initial_vertex_sha256": vertex_digest(initial_ordered), "initial_perimeter": p0,
        "minimum_all_adjacent_layer_containment_cross": global_containment,
        "layers": layers, "targets": target_records,
    }
