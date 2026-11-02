"""Raw fixed-frame evaluation of existing BDF support polynomials.

No ODE/PDE integration is imported or performed. Time is the unchanged
physical frozen clock t=alpha*n/R**(4/3); alpha is read from the archive.
The inherited SavedAffineReference is kept unchanged in its own source file.

Unlike the inherited ``support`` method, ``raw_support`` never divides by c0
and accepts derivative coefficient arrays whose mean is negative. The primary
reference displacement over [a,b] is h(b)-h(a), evaluated directly from the
saved polynomials. Instantaneous dense-polynomial derivatives are diagnostics;
their startup error and their jumps at BDF knots are reported separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np

from spectral_reference_base import SavedAffineReference


class FrozenTimeReference(SavedAffineReference):
    """Evaluate raw support and its time derivative from saved BDF output."""

    def __init__(self, path):
        super().__init__(path)
        with np.load(self.path, allow_pickle=False) as archive:
            if 'alpha_clock' not in archive.files:
                raise ValueError('Frozen reference has no recorded clock')
            self.alpha_clock = float(archive['alpha_clock'])
        if not np.isfinite(self.alpha_clock) or self.alpha_clock <= 0:
            raise ValueError('Frozen clock must be finite and positive')
        self._time_cache = {}

    def state_and_derivative_at_time(self, time, *, side='left'):
        """Return raw h grid, raw h_t grid, interval index and polynomial order.

        Product derivatives are evaluated by recurrence, without division by
        t-shift. Evaluation is therefore defined at every interpolation knot.
        The default left choice matches the unchanged base evaluator. ``side``
        exists only for the one-sided knot-continuity diagnostics.
        """
        time = float(time)
        if side not in ('left', 'right'):
            raise ValueError('side must be left or right')
        if not np.isfinite(time) or time < self.time_min - 2e-15 or time > self.time_max + 2e-15:
            raise ValueError('Time outside the saved-reference window')
        time = min(max(time, self.time_min), self.time_max)
        index = min(int(np.searchsorted(self.breakpoints[1:], time, side=side)),
                    len(self.orders) - 1)
        order = int(self.orders[index])
        values = self.dense_D[index]
        state = values[0].copy()
        derivative = np.zeros(self.n)
        product, product_derivative = 1.0, 0.0
        for j in range(order):
            x = (time - self.shifts[index, j]) / self.denominators[index, j]
            product_derivative, product = (
                product_derivative * x + product / self.denominators[index, j],
                product * x,
            )
            state += product * values[j + 1]
            derivative += product_derivative * values[j + 1]
        return state, derivative, index, order

    def at_time(self, time):
        """Raw support record plus raw time-derivative Fourier coefficients.

        Coefficients are FFT(h)/N and FFT(h_t)/N, with physical frequencies
        factor*modes_u. No perimeter normalization or mode-1 subtraction is
        applied. The existing ``a3`` field is only descriptive compatibility
        with the base record; primary physical comparisons use raw arrays.
        """
        time = float(time)
        if time not in self._time_cache:
            state, derivative, index, order = self.state_and_derivative_at_time(time)
            coefficients = self.coefficients_from_state(state)
            derivative_coefficients = self.coefficients_from_state(derivative)
            self._time_cache[time] = (coefficients, derivative_coefficients, index, order)
        coefficients, derivative_coefficients, index, order = self._time_cache[time]
        result = self._record(coefficients, time)
        result.update(derivative_coefficients=derivative_coefficients.copy(),
                      dense_interval=index, dense_order=order)
        return result

    def raw_coeffs_at_times(self, times):
        """Return an Ntime-by-Nfrequency array of raw support coefficients."""
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or not np.all(np.isfinite(times)):
            raise ValueError('times must be a finite one-dimensional array')
        if len(times) == 0:
            return np.empty((0, self.n), dtype=complex)
        return np.stack([self.at_time(float(value))['coefficients'] for value in times])

    def _check_raw_coefficients(self, coefficients, ndim):
        coefficients = np.asarray(coefficients, dtype=complex)
        if coefficients.ndim != ndim or coefficients.shape[-1] != self.n:
            raise ValueError('Unexpected raw-coefficient array shape')
        if not np.all(np.isfinite(coefficients)):
            raise ValueError('Non-finite raw Fourier coefficient')
        if self.n % 2 == 0 and np.any(np.abs(coefficients[..., self.n // 2].imag) > 5e-14):
            raise ValueError('Unique FFT Nyquist coefficient must be real')
        return coefficients

    def raw_support(self, coefficients, angles):
        """Raw Re(sum c_k exp(i*f_k*angle)), including negative-mean h_t."""
        coefficients = self._check_raw_coefficients(coefficients, 1)
        angles = np.asarray(angles, dtype=float)
        if not np.all(np.isfinite(angles)):
            raise ValueError('Angles must be finite')
        values = np.real(coefficients @ np.exp(1j * np.outer(self.frequencies, angles.reshape(-1))))
        values = values.reshape(angles.shape)
        return float(values) if values.ndim == 0 else values

    def raw_support_many(self, coefficients, angles):
        """Batch raw supports with shape (Ntime,*angles.shape)."""
        coefficients = self._check_raw_coefficients(coefficients, 2)
        angles = np.asarray(angles, dtype=float)
        if not np.all(np.isfinite(angles)):
            raise ValueError('Angles must be finite')
        values = np.real(coefficients @ np.exp(1j * np.outer(self.frequencies, angles.reshape(-1))))
        return values.reshape((len(coefficients),) + angles.shape)


def qualify(reference_path):
    """Check saved-output reconstruction and report empirical h_t accuracy.

    The PDE residual is a diagnostic of the already-saved interpolant, not a
    new solution or an interval enclosure. Structural checks and exact
    polynomial derivative identities determine PASS/FAIL. Instantaneous
    derivative accuracy is reported explicitly rather than assumed equal to
    the much smaller state-reconstruction accuracy.
    """
    started = time.perf_counter()
    ref = FrozenTimeReference(reference_path)
    root = Path(__file__).resolve().parent
    alpha_expected = 0.6450983601211226
    structural = []
    angles = np.linspace(-math.pi, 3 * math.pi, 517)
    initial = ref.at_time(0.0)
    initial_error = float(np.max(np.abs(ref.raw_support(initial['coefficients'], angles)
                                        - (1 - 3 / 40 * np.cos(3 * angles)))))
    structural.append({'check': 'initial_analytic_raw_support', 'error': initial_error,
                       'pass': initial_error < 5e-14})
    structural.append({'check': 'archive_clock_matches_frozen_value',
                       'value': ref.alpha_clock, 'expected': alpha_expected,
                       'pass': ref.alpha_clock == alpha_expected})

    # All layer times through the worst-case q=.95 upper brackets, together
    # with a uniform diagnostic sample reaching the declared .051 boundary.
    times = np.unique(np.concatenate((
        ref.alpha_clock * np.arange(127) / 256 ** (4 / 3),
        ref.alpha_clock * np.arange(314) / 512 ** (4 / 3),
        np.linspace(0.0, .051, 511),
    )))
    pde_errors = []
    state_base_errors = []
    minimum_rho = float('inf')
    for t in times:
        h, ht, index, order = ref.state_and_derivative_at_time(float(t))
        rho = np.fft.ifft((1 - ref.frequencies ** 2) * np.fft.fft(h)).real
        minimum_rho = min(minimum_rho, float(np.min(rho)))
        if np.any(rho <= 0):
            raise ValueError('Nonpositive saved-reference curvature radius')
        error = float(np.max(np.abs(ht + rho ** (-1 / 3))))
        pde_errors.append({'time': float(t), 'error': error, 'order': order})
        state_base_errors.append(float(np.max(np.abs(h - ref.state_at_time(float(t))))))
    max_base_error = max(state_base_errors)
    structural.append({'check': 'state_matches_unchanged_base', 'error': max_base_error,
                       'pass': max_base_error < 5e-14})

    # Independent algebra: convert each local product into a polynomial in
    # x=(t-left)/(right-left), differentiate its power coefficients, evaluate
    # at x=.37, and compare against the product-derivative recurrence.
    derivative_errors = []
    integral_errors = []
    gauss_x, gauss_w = np.polynomial.legendre.leggauss(4)
    for i, left in enumerate(ref.breakpoints[:-1]):
        if left > .051:
            break
        right = float(ref.breakpoints[i + 1])
        width = right - float(left)
        order = int(ref.orders[i])
        powers = np.zeros((order + 1, ref.n))
        powers[0] = ref.dense_D[i, 0]
        product = np.array([1.0])
        for j in range(order):
            factor = np.array([(float(left) - ref.shifts[i, j]) / ref.denominators[i, j],
                               width / ref.denominators[i, j]])
            product = np.polynomial.polynomial.polymul(product, factor)
            powers[:len(product)] += product[:, None] * ref.dense_D[i, j + 1]
        x = .37
        if order:
            derivative_power = np.arange(1, order + 1)[:, None] * powers[1:] / width
            independent = np.sum(derivative_power * x ** np.arange(order)[:, None], axis=0)
        else:
            independent = np.zeros(ref.n)
        t = float(left) + width * x
        _, actual, _, _ = ref.state_and_derivative_at_time(t)
        derivative_errors.append(float(np.max(np.abs(independent - actual))))
        # A 4-node Gauss rule integrates any degree<=7 exactly; BDF h_t is
        # degree<=4. This quadrature is over an existing polynomial only.
        nodes = (float(left) + right) / 2 + width / 2 * gauss_x
        derivatives = np.stack([ref.state_and_derivative_at_time(float(t))[1] for t in nodes])
        integrated = width / 2 * (gauss_w @ derivatives)
        before = ref.state_at_time(float(left))
        after = ref.state_at_time(right)
        integral_errors.append(float(np.max(np.abs(integrated - (after - before)))))
    max_derivative_identity = max(derivative_errors)
    max_integral_identity = max(integral_errors)
    structural.append({'check': 'independent_power_polynomial_derivative',
                       'error': max_derivative_identity, 'pass': max_derivative_identity < 5e-12})
    structural.append({'check': 'integral_of_derivative_equals_saved_support_difference',
                       'error': max_integral_identity, 'pass': max_integral_identity < 5e-14})

    state_jumps, derivative_jumps = [], []
    for knot in ref.breakpoints[1:-1]:
        if knot > .051:
            break
        left, dleft, _, _ = ref.state_and_derivative_at_time(float(knot), side='left')
        right, dright, _, _ = ref.state_and_derivative_at_time(float(knot), side='right')
        state_jumps.append(float(np.max(np.abs(left - right))))
        derivative_jumps.append(float(np.max(np.abs(dleft - dright))))
    structural.append({'check': 'state_continuity_at_knots', 'error': max(state_jumps),
                       'pass': max(state_jumps) < 5e-14})
    batch_times = np.array([0.0, .001, .02, .04, .051])
    batch_coefficients = ref.raw_coeffs_at_times(batch_times)
    raw_many = ref.raw_support_many(batch_coefficients, angles[:13])
    repeated = np.stack([ref.raw_support(ref.at_time(float(t))['coefficients'], angles[:13]) for t in batch_times])
    batch_error = float(np.max(np.abs(raw_many - repeated)))
    structural.append({'check': 'batch_support_matches_individual', 'error': batch_error,
                       'pass': batch_error < 5e-14})
    negative_mean = ref.at_time(.02)['derivative_coefficients']
    derivative_value = ref.raw_support(negative_mean, .3)
    structural.append({'check': 'raw_support_accepts_negative_mean_derivative',
                       'derivative_c0': float(negative_mean[0].real),
                       'example_value': derivative_value,
                       'pass': bool(negative_mean[0].real < 0 and np.isfinite(derivative_value))})

    worst_pde = max(pde_errors, key=lambda row: row['error'])
    worst_nonzero = max((row for row in pde_errors if row['time'] > 0), key=lambda row: row['error'])
    worst_after_startup = max((row for row in pde_errors if row['time'] >= .001), key=lambda row: row['error'])
    return {
        'status': 'PASS_SAVED_TIME_REFERENCE_STRUCTURAL_CHECKS' if all(row['pass'] for row in structural) else 'FAIL_SAVED_TIME_REFERENCE_STRUCTURAL_CHECKS',
        'date': '2026-10-30',
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'base_source_sha256': hashlib.sha256((root / 'spectral_reference_base.py').read_bytes()).hexdigest(),
        'reference_sha256': hashlib.sha256(Path(reference_path).read_bytes()).hexdigest(),
        'environment': {'python': platform.python_version(), 'numpy': np.__version__},
        'clock_alpha': ref.alpha_clock, 'saved_time_domain': [ref.time_min, ref.time_max],
        'checked_time_domain': [0.0, .051], 'sampled_time_count': len(times),
        'checks': structural, 'minimum_sampled_rho': minimum_rho,
        'max_sampled_ht_plus_rho_minus_one_third': worst_pde,
        'max_sampled_nonzero_time_pde_residual': worst_nonzero,
        'max_sampled_pde_residual_at_or_after_001': worst_after_startup,
        'knot_count_checked': len(state_jumps),
        'max_state_knot_jump': max(state_jumps),
        'max_derivative_knot_jump': max(derivative_jumps),
        'primary_displacement_method': 'Direct difference of saved raw supports h(t_b)-h(t_a); not h_t(t_a)*(t_b-t_a).',
        'derivative_accuracy_limit': 'Dense state accuracy does not certify equal derivative accuracy; BDF startup and one-sided derivative jumps are explicitly retained.',
        'new_pde_solve': False, 'new_peeling': False, 'clock_refitted': False,
        'wall_seconds': time.perf_counter() - started,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', default=str(Path(__file__).resolve().parent / 'inputs/reference_contraction.npz'))
    parser.add_argument('--output', default=str(Path(__file__).resolve().parent / 'qualification_time_reference.json'))
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise FileExistsError('Choose a fresh qualification output')
    result = qualify(args.reference)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf8')
    print(result['status'], 'sampled_times', result['sampled_time_count'],
          'max_pde_residual', result['max_sampled_ht_plus_rho_minus_one_third']['error'])
    raise SystemExit(0 if result['status'].startswith('PASS') else 1)
