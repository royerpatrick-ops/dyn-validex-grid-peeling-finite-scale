"""Exact-cone sector decomposition of a convex polygon's third support mode.

The constructor takes unshifted lattice vertices in strict CCW order. Physical
body-frame coordinates are ((vertex + translation) / radiusR) * exp(-i*theta).
Nothing is centered or normalized by mean support here. weighted_integrals
returns integral h_body(phi)*cos(3*phi) dphi / (2*pi), on arbitrary intervals.
All finite elementary-function evaluations use float64, not intervals.
"""
from __future__ import annotations

import cmath
import math

import numpy as np


TAU = 2 * math.pi


def _intervals(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError('intervals must have shape (N,2)')
    if not np.all(np.isfinite(values)) or np.any(values[:, 1] < values[:, 0]):
        raise ValueError('finite intervals with upper >= lower are required')
    return values


def _weighted_piece(x, y, lower, width):
    """Primitive difference for (x cos(phi)+y sin(phi))*cos(3 phi).

    Midpoint/sine differences avoid cancellation on very narrow normal cones.
    Broadcasting is intentional, for arrays of requested interval endpoints.
    """
    mid = lower + width / 2
    int_cos2 = np.sin(width) * np.cos(2 * mid)
    int_cos4 = np.sin(2 * width) * np.cos(4 * mid) / 2
    int_sin2 = np.sin(width) * np.sin(2 * mid)
    int_sin4 = np.sin(2 * width) * np.sin(4 * mid) / 2
    return (x * (int_cos2 + int_cos4) + y * (int_sin4 - int_sin2)) / 2


def dipole_support(c1_shape, angles):
    """The first support harmonic 2 Re(c1*exp(i*phi))."""
    c1_shape = complex(c1_shape)
    phi = np.asarray(angles, dtype=float)
    return 2 * (c1_shape.real * np.cos(phi) - c1_shape.imag * np.sin(phi))


def dipole_weighted_integrals(c1_shape, intervals):
    """Integral of the first harmonic times cos(3 phi), divided by 2*pi."""
    intervals = _intervals(intervals)
    c1_shape = complex(c1_shape)
    return _weighted_piece(2 * c1_shape.real, -2 * c1_shape.imag,
                           intervals[:, 0], intervals[:, 1] - intervals[:, 0]) / TAU


class PolygonSupport:
    """Piecewise analytic polygon support in the shape's body frame."""

    dipole_support = staticmethod(dipole_support)
    dipole_weighted_integrals = staticmethod(dipole_weighted_integrals)

    def __init__(self, vertices, theta, translation=(0.0, 0.0), radiusR=1.0):
        raw = np.asarray(vertices)
        if raw.ndim != 2 or raw.shape[1] != 2 or len(raw) < 3:
            raise ValueError('At least three 2D vertices are required')
        if not np.all(np.isfinite(raw)):
            raise ValueError('Vertices must be finite')
        if not math.isfinite(float(radiusR)) or radiusR <= 0:
            raise ValueError('radiusR must be positive and finite')
        if not math.isfinite(float(theta)):
            raise ValueError('theta must be finite')
        translation = np.asarray(translation, dtype=float)
        if translation.shape != (2,) or not np.all(np.isfinite(translation)):
            raise ValueError('translation must contain two finite coordinates')
        # Convert individual scalars before products: integer input does not
        # inherit NumPy int64 overflow in the strict-convexity predicate.
        points = [tuple(v.item() if hasattr(v, 'item') else v for v in p) for p in raw]
        for j, b in enumerate(points):
            a, c = points[j - 1], points[(j + 1) % len(points)]
            cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
            if cross <= 0:
                raise ValueError('Vertices must be strictly convex CCW')
        xy = (np.asarray(raw, dtype=float) + translation) / float(radiusR)
        z = (xy[:, 0] + 1j * xy[:, 1]) * cmath.exp(-1j * float(theta))
        self.positions = np.column_stack((z.real, z.imag))
        self.theta = float(theta)
        self.translation = tuple(map(float, translation))
        self.radiusR = float(radiusR)
        edge = np.roll(self.positions, -1, axis=0) - self.positions
        # For a CCW edge (dx,dy), outward normal is (dy,-dx).
        normals = np.mod(np.arctan2(-edge[:, 0], edge[:, 1]), TAU)
        pieces = []
        for j, pos in enumerate(self.positions):
            lower = float(normals[j - 1])
            upper = float(normals[j])
            if upper <= lower:
                upper += TAU
            if not 0 < upper - lower < math.pi:
                raise ArithmeticError('Invalid normal cone width')
            if upper > TAU:
                pieces.append((lower, TAU, float(pos[0]), float(pos[1])))
                pieces.append((0.0, upper - TAU, float(pos[0]), float(pos[1])))
            else:
                pieces.append((lower, upper, float(pos[0]), float(pos[1])))
        pieces.sort(key=lambda p: p[0])
        a = np.asarray(pieces, dtype=float)
        if a[0, 0] != 0 or abs(a[-1, 1] - TAU) > 4e-15:
            raise ArithmeticError('Normal cones do not span one full turn')
        if len(a) > 1 and np.max(np.abs(a[:-1, 1] - a[1:, 0])) > 4e-15:
            raise ArithmeticError('Normal cone partition is not contiguous')
        self._lower, self._upper = a[:, 0], a[:, 1]
        self._x, self._y = a[:, 2], a[:, 3]
        pieces_integral = _weighted_piece(self._x, self._y, self._lower, self._upper - self._lower)
        self._prefix = np.concatenate(([0.0], np.cumsum(pieces_integral)))
        self._period_integral = math.fsum(map(float, pieces_integral))
        self._prefix[-1] = self._period_integral

    def support(self, angles):
        """Evaluate the uncentered, unnormalized body-frame support."""
        angles = np.asarray(angles, dtype=float)
        if not np.all(np.isfinite(angles)):
            raise ValueError('angles must be finite')
        phi = np.mod(angles, TAU)
        idx = np.searchsorted(self._upper, phi, side='right')
        idx = np.minimum(idx, len(self._upper) - 1)
        return self._x[idx] * np.cos(phi) + self._y[idx] * np.sin(phi)

    def _weighted_antiderivative(self, angles):
        angles = np.asarray(angles, dtype=float)
        turns = np.floor(angles / TAU)
        phi = angles - turns * TAU
        idx = np.searchsorted(self._upper, phi, side='right')
        idx = np.minimum(idx, len(self._upper) - 1)
        partial = _weighted_piece(self._x[idx], self._y[idx], self._lower[idx], phi - self._lower[idx])
        return turns * self._period_integral + self._prefix[idx] + partial

    def weighted_integrals(self, intervals):
        """Return raw signed mode-3 contributions on arbitrary intervals.

        Endpoints may be outside [0,2*pi], and interval length may exceed a turn.
        The return values already include 1/(2*pi), but not division by c0.
        """
        intervals = _intervals(intervals)
        return (self._weighted_antiderivative(intervals[:, 1])
                - self._weighted_antiderivative(intervals[:, 0])) / TAU

    def coefficients(self, modes=(0, 1, 3)):
        """Raw body-frame Fourier coefficients from the normal-cone pieces."""
        modes = tuple(modes)
        if any(not isinstance(m, int) for m in modes) or len(set(modes)) != len(modes):
            raise ValueError('modes must be distinct Python integers')
        result = {}
        lower, width = self._lower, self._upper - self._lower
        z = self._x + 1j * self._y

        def integral_exp(k):
            if k == 0:
                return width.astype(complex)
            return (2 * np.sin(k * width / 2) / k) * np.exp(-1j * k * (lower + width / 2))

        for m in modes:
            values = (z.conjugate() * integral_exp(m - 1) + z * integral_exp(m + 1)) / (2 * TAU)
            result[m] = complex(math.fsum(map(float, values.real)), math.fsum(map(float, values.imag)))
        return result
