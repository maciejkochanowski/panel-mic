"""Independent exact-oracle regressions for the 2.0 release audit findings."""
import math
from fractions import Fraction

import pytest

from panel_mic import Reference, diagnose_upper, frontier_upper, protocol_assurance
from panel_mic._version import __version__
from panel_mic.numerics import (
    binomial_cutoff, certification_power_lower, cp_lower, cp_upper,
    minimum_reference,
)


def exact_cdf(n, k, q):
    """Closed polynomial, independent of the production positive recurrence."""
    if k < 0:
        return Fraction(0)
    if k >= n:
        return Fraction(1)
    a, b = Fraction(q).as_integer_ratio()
    return Fraction(sum(math.comb(n, j)*a**j*(b-a)**(n-j)
                        for j in range(k+1)), b**n)


@pytest.mark.parametrize('n,k,alpha', [
    (4, 1, .05), (8, 0, .05), (10, 1, .025), (874, 0, .01),
    (1250, 0, .05639926893788604), (1303, 0, .01),
    (1499, 0, .001), (3, 1, .125), (3, 3, .05),
])
def test_cp_upper_is_on_safe_side_of_exact_inverse(n, k, alpha):
    upper = cp_upper(k, n, alpha)
    if k == n:
        assert upper == 1.
    else:
        assert exact_cdf(n, k, upper) <= Fraction(alpha)
    lower = cp_lower(k, n, alpha)
    if k == 0:
        assert lower == 0.
    else:
        assert 1-exact_cdf(n, k-1, lower) <= Fraction(alpha)


def test_frontier_does_not_certify_the_audited_unsafe_cp_boundary():
    result = frontier_upper(Reference([1], [1250, 0], iid=True), [1],
        future_ns=[1000], target_power=.8, delta=.05639926893788604,
        grid_prespecified=True)
    row = result['rows'][0]
    assert row['status'] == 'NO_CERTIFIED_BOUNDARY'
    assert row['certified_power_floor'] is None
    candidate = row['candidate_bounds'][0]
    exact = exact_cdf(1000, 3, candidate['worst_tail_upper_bound'])
    assert Fraction(candidate['power_lower_bound']) <= exact < Fraction(.8)


@pytest.mark.parametrize('nr,delta,target', [
    (874, .01, .14384498882876628),
    (1303, .01, .2723639692336389),
    (1499, .001, .1834463208152981),
])
def test_diagnostic_cutoff_windows_and_meets_target_reject_same_boundary(nr, delta, target):
    result = diagnose_upper(Reference([1], [nr, 0], iid=True), 1,
        future_n=368, max_future_n=368, target_power=target, delta=delta)
    analysis = result['conditional_analysis']
    current = analysis['current_scenario']
    assert current['reference_cutoff'] is None
    assert not current['meets_target']
    assert analysis['future_n_intervals'] == []
    assert Fraction(current['power_lower_bound']) <= exact_cdf(
        368, 0, analysis['worst_tail_upper']) < Fraction(target)


def test_minimum_reference_proves_both_sides_of_logarithmic_tie():
    q = .002991249545095296
    n = minimum_reference(q, .05)
    assert n == 1001
    assert exact_cdf(n, 0, q) <= Fraction(.05)
    assert exact_cdf(n-1, 0, q) > Fraction(.05)
    assert binomial_cutoff(n, q, .05) == 0
    assert binomial_cutoff(n-1, q, .05) is None


def test_minimum_reference_exact_equality_and_large_zero_count_regime():
    assert minimum_reference(.5, .25) == 2
    q = .0001
    n = minimum_reference(q, .05)
    assert n > 10000
    assert exact_cdf(n, 0, q) <= Fraction(.05)
    assert exact_cdf(n-1, 0, q) > Fraction(.05)
    assert binomial_cutoff(n, q, .05) == 0
    assert binomial_cutoff(n-1, q, .05) is None


@pytest.mark.parametrize('nr,nc,feasible_stage', [(1139, 1401, 0), (1138, 1402, 1)])
@pytest.mark.parametrize('independent', [False, True])
def test_infeasible_other_stage_preserves_feasible_marginal(nr, nc, feasible_stage, independent):
    result = protocol_assurance(nr, nc, future_n=874,
        possible_tail_reference=.001, possible_tail_control=.001,
        independent_stages=independent)
    expected = float((1-Fraction(.001))**(nr if feasible_stage == 0 else nc))
    assert result['stage_success_probabilities'][feasible_stage] == pytest.approx(expected)
    assert result['stage_success_probabilities'][1-feasible_stage] == 0.
    assert result['joint_success_probability'] == 0.
    assert result['joint_success_interval'] == [0., 0.]
    assert result['status'] == 'REFERENCE_OR_CONTROL_INFORMATION_LIMIT'
    assert result['software_version'] == __version__
    assert result['schema_version'] == 'panel_mic.protocol/v2'


@pytest.mark.parametrize('q', [0., .001, .0022975957430846745, .5, 1.])
def test_reported_power_lower_is_an_exact_lower_enclosure(q):
    assert Fraction(certification_power_lower(1000, q)) <= exact_cdf(1000, 3, q)


def test_success_floor_is_the_declared_target():
    result = frontier_upper(Reference([1], [2000, 0], iid=True), [1],
        future_ns=[1000], grid_prespecified=True)
    row = result['rows'][0]
    assert row['status'] == 'FOUND'
    assert row['certified_power_floor'] == result['target_power'] == .8


def test_unresolved_numerics_do_not_invent_a_bound_or_an_exact_minimum(monkeypatch):
    from panel_mic import numerics
    numerics._cp_upper.cache_clear()
    monkeypatch.setattr(numerics, '_cdf_relation', lambda *args: None)
    assert numerics.cp_upper(1, 12, .05) == 1.
    assert numerics.minimum_reference(.01, .05) is None
    numerics._cp_upper.cache_clear()
