"""Independent review regressions: truthful intervals and retained batch cohorts."""
import json
import math
from fractions import Fraction
from math import comb

import pytest

from panel_mic import Reference, ValidationError, design_upper, frontier_upper
from panel_mic.cli import run_batch
from panel_mic.numerics import binomial_cutoff


def test_subnormal_drift_does_not_drop_true_exceedances():
    s = math.ulp(0.0)
    reference = Reference.from_intervals(
        [{'lower': s, 'upper': 2*s, 'count': 2000}], iid=True)
    # A degenerate iid latent population at 2*s truthfully fits this interval.
    # The admissible shifted population is at 4*s, ABOVE the candidate 3*s.
    result = design_upper(reference, [3*s], future_n=1000, drift=1)
    assert result['status'] != 'FOUND'
    assert result['candidate_diagnostics'][0]['possible_exceedances'] == 2000


def test_subnormal_drift_frontier_does_not_certify_zero_power():
    s = math.ulp(0.0)
    reference = Reference.from_intervals(
        [{'lower': s, 'upper': 2*s, 'count': 2000}], iid=True)
    result = frontier_upper(reference, [3*s], future_ns=[1000], drifts=[1],
                            grid_prespecified=True)
    assert result['rows'][0]['status'] != 'FOUND'
    assert result['rows'][0]['candidate_bounds'][0]['power_lower_bound'] == 0


@pytest.mark.parametrize('metadata', [
    {'source_score': math.inf}, {'nested': [math.nan]}, {'unsupported': object()},
])
def test_metadata_cannot_break_promised_json_output(metadata):
    with pytest.raises(ValidationError):
        Reference([1], [2000, 0], iid=True, metadata=metadata)


def test_batch_retains_valid_cohort_beside_nonfinite_metadata():
    # Valid JSON decimal notation overflows binary64 without parse_constant.
    doc = json.loads('''{"cohorts":[
      {"id":"ok","operation":"assess-upper","reference":
       {"thresholds":[1],"counts":[2000,0],"iid":true},"upper_limit":1},
      {"id":"overflow","operation":"assess-upper","reference":
       {"thresholds":[1],"counts":[2000,0],"iid":true,
        "metadata":{"source_score":1e999}},"upper_limit":1}]}''')
    result = run_batch(doc)
    assert result['cohorts'][0]['status'] == 'PANEL_REACHES_QUANTILE'
    assert result['cohorts'][1]['status'] == 'INVALID_INPUT'
    json.dumps(result, allow_nan=False)


def test_small_cutoffs_against_independent_rational_enumeration():
    alpha = Fraction(1, 16)
    for n in range(1, 25):
        for q in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)):
            cumulative = Fraction(0)
            expected = None
            for k in range(n + 1):
                cumulative += comb(n, k) * q**k * (1-q)**(n-k)
                if cumulative <= alpha:
                    expected = k
            assert binomial_cutoff(n, float(q), float(alpha)) == expected


def test_shift_threshold_is_exact_or_directed_down():
    from panel_mic.upper import _shift_threshold
    s = math.ulp(0.0)
    candidates = [i*s for i in range(1, 129)]
    candidates += [math.ldexp(1.23456789, e) for e in (-1022, -1021, -500, 0, 500, 1023)]
    for c in candidates:
        for drift in (0, 1, 2, 3, 10, 100, 1023):
            exact = Fraction(c) / 2**drift
            transformed = _shift_threshold(c, drift)
            assert Fraction(transformed) <= exact
            following = math.nextafter(transformed, math.inf)
            if math.isfinite(following):
                assert Fraction(following) > exact


def test_effective_target_rounding_is_conservative():
    reference = Reference([1], [20000, 0], iid=True)
    eta = 0.00013062592196620468
    result = design_upper(reference, [1], future_n=368, unmodelled_share=eta)
    exact = (Fraction(result['q_star']) - Fraction(eta)) / (1 - Fraction(eta))
    assert Fraction(result['q_eff']) <= exact


def test_nested_metadata_keys_cannot_silently_collide():
    with pytest.raises(ValidationError):
        Reference([1], [1, 0], metadata={'source': {1: 'discarded', '1': 'retained'}})


def test_metadata_nested_input_is_copied():
    original = {'source': {'labels': ['initial']}}
    reference = Reference([1], [1, 0], metadata=original)
    original['source']['labels'].append('later')
    assert reference.metadata == {'source': {'labels': ['initial']}}


def test_small_exact_cdf_and_strict_sf_ties():
    from panel_mic.numerics import excess_tail_significant, _cdf_at_most
    for n in range(2, 10):
        for q in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)):
            cumulative = Fraction(0)
            for k in range(n + 1):
                cumulative += comb(n, k) * q**k * (1-q)**(n-k)
                if 0 < cumulative < Fraction(1, 2):
                    a = float(cumulative)
                    assert _cdf_at_most(k, n, float(q), a)
                    assert not _cdf_at_most(k, n, float(q), math.nextafter(a, 0))
                survival = 1 - cumulative
                if 0 < survival < Fraction(1, 2):
                    a = float(survival)
                    assert not excess_tail_significant(k+1, n, float(q), a)
                    assert excess_tail_significant(k+1, n, float(q), math.nextafter(a, 1))
