"""Evaluate a saved affine-flow support reference, without solving an ODE/PDE.

Only NumPy and scipy.optimize.brentq are used. The archive stores BDF dense
polynomials for the support h(u), u=3*phi. Fourier coefficients are RAW
coefficients c_k=FFT(h)/N; physical angular frequencies are 3*k.

``at_q`` and ``static`` return dictionaries containing raw coefficients,
physical frequencies, c0, normalized a3 and time. ``support`` and
``weighted_integrals`` divide the raw series by Re(c0), so their support is
normalized to mean one. No first harmonic is present in this triple-symmetric
reference.

The single even-grid FFT Nyquist coefficient is real. It is retained once,
without doubling, and interpreted through the real part of the Fourier
series, hence as one real cosine. This applies both to support evaluation
and analytic interval integration.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.optimize import brentq


class SavedAffineReference:
    """Read-only evaluator for ``reference_contraction.npz``.

    The archive's dense polynomials are the qualified numerical reference;
    this class does not turn their numerical qualification into interval
    certification. It never imports an integration routine or source module.
    """

    def __init__(self, path):
        self.path = Path(path)
        required = (
            'dense_breakpoints', 'dense_order', 'dense_D', 'dense_t_shift',
            'dense_denom', 'factor', 'fourier_modes_u',
        )
        with np.load(self.path, allow_pickle=False) as archive:
            missing = [key for key in required if key not in archive.files]
            if missing:
                raise ValueError('Missing dense-reference fields: ' + ', '.join(missing))
            self.breakpoints = np.asarray(archive['dense_breakpoints'], dtype=float).copy()
            self.orders = np.asarray(archive['dense_order'], dtype=np.int64).copy()
            self.dense_D = np.asarray(archive['dense_D'], dtype=float).copy()
            self.shifts = np.asarray(archive['dense_t_shift'], dtype=float).copy()
            self.denominators = np.asarray(archive['dense_denom'], dtype=float).copy()
            self.factor = int(archive['factor'])
            self.modes_u = np.asarray(archive['fourier_modes_u'], dtype=np.int64).copy()
        if self.factor != 3:
            raise ValueError('This evaluator requires the triple-symmetric factor=3 archive')
        if self.breakpoints.ndim != 1 or len(self.breakpoints) < 2 or np.any(np.diff(self.breakpoints) <= 0):
            raise ValueError('Dense breakpoints must be strictly increasing')
        n_intervals = len(self.breakpoints) - 1
        if self.dense_D.ndim != 3 or self.dense_D.shape[:2] != (n_intervals, 6):
            raise ValueError('Unexpected dense-polynomial shape')
        self.n = int(self.dense_D.shape[2])
        if self.orders.shape != (n_intervals,) or np.any((self.orders < 0) | (self.orders > 5)):
            raise ValueError('Unexpected BDF orders')
        if self.shifts.shape != (n_intervals, 5) or self.denominators.shape != (n_intervals, 5):
            raise ValueError('Unexpected dense shift/denominator shapes')
        if self.modes_u.shape != (self.n,):
            raise ValueError('Unexpected Fourier-mode count')
        fft_modes = np.rint(np.fft.fftfreq(self.n) * self.n).astype(np.int64)
        if not np.array_equal(self.modes_u, fft_modes):
            raise ValueError('Saved frequencies do not use NumPy full-FFT ordering')
        for array in (self.breakpoints, self.dense_D, self.shifts, self.denominators):
            if not np.all(np.isfinite(array)):
                raise ValueError('Non-finite dense-reference value')
        for index, order in enumerate(self.orders):
            if np.any(self.denominators[index, :int(order)] == 0):
                raise ValueError('Zero active dense-polynomial denominator')
        self.frequencies = self.factor * self.modes_u
        self.index3 = int(np.flatnonzero(self.frequencies == 3)[0])
        self.index_minus3 = int(np.flatnonzero(self.frequencies == -3)[0])
        self.time_min = float(self.breakpoints[0])
        self.time_max = float(self.breakpoints[-1])
        self.q_max = float(np.mean(self.state_at_time(self.time_min)))
        self.q_min = float(np.mean(self.state_at_time(self.time_max)))
        if not self.q_max > self.q_min > 0:
            raise ValueError('Reference has no positive decreasing mean-support range')
        self._cache = {}
        self._integral_cache = {}

    def state_at_time(self, time):
        """Reconstruct one saved BDF polynomial at a prescribed window time."""
        time = float(time)
        if not np.isfinite(time) or time < self.time_min - 2e-15 or time > self.time_max + 2e-15:
            raise ValueError('Time outside the saved-reference window')
        time = min(max(time, self.time_min), self.time_max)
        index = min(int(np.searchsorted(self.breakpoints[1:], time, side='left')),
                    len(self.orders) - 1)
        order = int(self.orders[index])
        x = (time - self.shifts[index, :order]) / self.denominators[index, :order]
        products = np.cumprod(x)
        values = self.dense_D[index]
        return values[0] + products @ values[1:order + 1]

    def coefficients_from_state(self, state):
        """Return RAW support coefficients FFT(state)/N, not normalized ones."""
        state = np.asarray(state, dtype=float)
        if state.shape != (self.n,) or not np.all(np.isfinite(state)):
            raise ValueError('State must contain N finite real support values')
        coefficients = np.fft.fft(state) / self.n
        coefficients[0] = complex(float(np.mean(state)), 0.0)
        if self.n % 2 == 0:
            coefficients[self.n // 2] = complex(coefficients[self.n // 2].real, 0.0)
        return coefficients

    def _record(self, coefficients, time):
        coefficients = self._check_coefficients(coefficients)
        c0 = float(coefficients[0].real)
        return {
            'coefficients': coefficients.copy(),
            'frequencies': self.frequencies.copy(),
            'c0': c0,
            'a3': complex(coefficients[self.index3] / c0),
            'time': float(time),
        }

    def at_q(self, q_abs):
        """Evaluate the saved support at its actual absolute mean q_abs.

        The only root search evaluates the existing polynomial mean. No
        evolution equation is integrated. ``c0`` records the reconstructed
        numerical mean, rather than silently imposing q_abs after the root.
        """
        q_abs = float(q_abs)
        if not np.isfinite(q_abs) or not self.q_min <= q_abs <= self.q_max:
            raise ValueError(f'q_abs outside [{self.q_min:.17g}, {self.q_max:.17g}]')
        if q_abs not in self._cache:
            if q_abs == self.q_max:
                time = self.time_min
            elif q_abs == self.q_min:
                time = self.time_max
            else:
                time = brentq(
                    lambda t: float(np.mean(self.state_at_time(t))) - q_abs,
                    self.time_min, self.time_max,
                    xtol=2e-15, rtol=8 * np.finfo(float).eps,
                )
            coefficients = self.coefficients_from_state(self.state_at_time(time))
            self._cache[q_abs] = (coefficients, time)
        coefficients, time = self._cache[q_abs]
        return self._record(coefficients, time)

    def static(self):
        """Exact initial support h0=1-3/40*cos(3*phi), in the same interface."""
        coefficients = np.zeros(self.n, dtype=complex)
        coefficients[0] = 1.0
        coefficients[self.index3] = -3 / 80
        coefficients[self.index_minus3] = -3 / 80
        return self._record(coefficients, self.time_min)

    def _check_coefficients(self, coefficients):
        coefficients = np.asarray(coefficients, dtype=complex)
        if coefficients.shape != (self.n,) or not np.all(np.isfinite(coefficients)):
            raise ValueError('Expected N finite raw complex coefficients')
        c0 = coefficients[0]
        if c0.real <= 0 or abs(c0.imag) > 5e-14:
            raise ValueError('Mean support must be positive and real')
        if self.n % 2 == 0 and abs(coefficients[self.n // 2].imag) > 5e-14:
            raise ValueError('The single FFT Nyquist coefficient must be real')
        return coefficients

    def support(self, coefficients, angles):
        """Evaluate mean-one H(phi)=Re(sum c_k exp(i*f_k*phi))/Re(c0).

        Angles can have any finite real values and any array shape. They are
        not rounded to an angular grid. A scalar angle returns a scalar.
        """
        coefficients = self._check_coefficients(coefficients)
        angles = np.asarray(angles, dtype=float)
        if not np.all(np.isfinite(angles)):
            raise ValueError('Angles must be finite')
        values = np.real(coefficients @ np.exp(1j * np.outer(self.frequencies, angles.reshape(-1))))
        result = (values / coefficients[0].real).reshape(angles.shape)
        return float(result) if result.ndim == 0 else result

    def weighted_integrals(self, coefficients, intervals):
        """Integrate normalized H(phi)*cos(3*phi)/(2*pi) on each interval.

        Each row is [left,right], with arbitrary finite real endpoints,
        including endpoints outside [0,2*pi]. Integration is oriented if
        right < left. Exact exponential antiderivatives are evaluated with
        a stable sinc expression; no angular quadrature is used.
        """
        coefficients = self._check_coefficients(coefficients)
        intervals = np.asarray(intervals, dtype=float)
        if intervals.ndim != 2 or intervals.shape[1] != 2 or not np.all(np.isfinite(intervals)):
            raise ValueError('Intervals must be a finite N-by-2 array')
        # A localization reuses the same angular bins for many translations.
        # Cache the coefficient-independent antiderivative kernel, with a
        # bounded number of keys so arbitrary user interval sets stay modest.
        key = (intervals.shape, intervals.tobytes())
        weights = self._integral_cache.get(key)
        if weights is None:
            widths = intervals[:, 1] - intervals[:, 0]
            midpoints = (intervals[:, 1] + intervals[:, 0]) / 2

            def exponential_integrals(frequencies):
                product = np.outer(widths, frequencies)
                return (widths[:, None] * np.sinc(product / (2 * np.pi))
                        * np.exp(1j * np.outer(midpoints, frequencies)))

            weights = (exponential_integrals(self.frequencies + 3)
                       + exponential_integrals(self.frequencies - 3)) / 2
            if len(self._integral_cache) >= 32:
                del self._integral_cache[next(iter(self._integral_cache))]
            self._integral_cache[key] = weights
        result = np.real(weights @ coefficients) / (2 * np.pi * coefficients[0].real)
        return result
