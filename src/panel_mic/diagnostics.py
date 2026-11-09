"""Fixed-upper diagnostic certificates and complete, nonmonotone future-N windows.

Reference sample size is fixed. This is not optional-stopping inference for a
growing reference cohort. A pointwise bound is not a simultaneous c/D band.
"""
import math
from functools import lru_cache
import numpy as np
from scipy.stats import binom
from .errors import ValidationError
from .models import integer, real
from .numerics import (probability, error_probability, percentile_tail, cp_upper,
                       binomial_cutoff, tail_target, minimum_reference, _cdf_at_most,
                       _decimal_power_at_least, certification_power_lower)
from .upper import (_base, _drift, _shift_threshold, _effective_target,
                    _worst_tail_upper)


def _null_accepts(ks, ns, epsilon, alpha):
    values = np.asarray(binom.cdf(ks, ns, epsilon))
    answer = values <= alpha
    near = np.flatnonzero(np.abs(values-alpha) <= 1e-11*alpha)
    for i in near:
        answer[i] = _cdf_at_most(int(ks[i]), int(ns[i]), epsilon, alpha)
    return answer


def _power_at_least(ks, ns, q, target):
    values = np.asarray(binom.cdf(ks, ns, q))
    answer = values >= target
    for i in np.flatnonzero(np.abs(values-target) <= 1e-12*target):
        answer[i] = _decimal_power_at_least(int(ns[i]), int(ks[i]), q, target)
    return answer


@lru_cache(maxsize=8)
def _critical_blocks(epsilon, alpha, cap):
    last = binomial_cutoff(cap, epsilon, alpha)
    if last is None:
        empty = np.array([], dtype=np.int64)
        empty.setflags(write=False)
        return empty, empty, empty
    ks = np.arange(last+1, dtype=np.int64)
    lo = ks+1
    hi = np.full_like(lo, cap)
    while np.any(lo < hi):
        active = np.flatnonzero(lo < hi)
        mid = (lo[active]+hi[active])//2
        good = _null_accepts(ks[active], mid, epsilon, alpha)
        hi[active] = np.where(good, mid, hi[active])
        lo[active] = np.where(good, lo[active], mid+1)
    starts = lo
    ends = np.r_[starts[1:]-1, cap]
    for array in (ks, starts, ends):
        array.setflags(write=False)
    return ks, starts, ends


@lru_cache(maxsize=256)
def _windows(q, percentile, target_power, alpha, cap, budget):
    epsilon = percentile_tail(percentile)
    if q >= epsilon:
        return ()
    ks, starts, ends = _critical_blocks(epsilon, alpha, cap)
    if not len(ks):
        return ()
    good = (ks >= budget) & _power_at_least(ks-budget, starts, q, target_power)
    cuts = ks[good]-budget
    beginnings = starts[good]
    lo, hi = beginnings.copy(), ends[good].copy()
    while np.any(lo < hi):
        active = np.flatnonzero(lo < hi)
        mid = (lo[active]+hi[active]+1)//2
        enough = _power_at_least(cuts[active], mid, q, target_power)
        lo[active] = np.where(enough, mid, lo[active])
        hi[active] = np.where(enough, hi[active], mid-1)
    merged = []
    for start, end in zip(beginnings, lo):
        a, b = int(start), int(end)
        if merged and merged[-1][1]+1 == a:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return tuple(merged)


def future_n_intervals(q, *, percentile=.99, target_power=.8, alpha=.025,
                       max_future_n=1000000, future_ambiguous_budget=0):
    """Numerical power windows for an input tail bound; not a confidence claim."""
    q = real(q, 'q', minimum=0, maximum=1)
    percentile_tail(percentile)
    target_power = probability(target_power, 'target_power')
    alpha = error_probability(alpha)
    if target_power <= alpha:
        raise ValidationError('diagnostic target_power must exceed alpha')
    cap = integer(max_future_n, 'max_future_n', 1)
    budget = integer(future_ambiguous_budget, 'future_ambiguous_budget')
    if cap > 10000000:
        raise ValidationError('max_future_n exceeds diagnostic limit 10000000')
    if q == 0:
        ks, starts, _ = _critical_blocks(percentile_tail(percentile), alpha, cap)
        return [[int(starts[budget]), cap]] if len(ks)>budget else []
    return [list(row) for row in _windows(q, percentile, target_power, alpha, cap, budget)]


def maximum_future_ambiguous_budget(n, q, *, percentile=.99, target_power=.8, alpha=.025):
    """Numerical budget from a supplied tail bound; its protection is caller-scoped."""
    n=integer(n,'n',1);q=real(q,'q',minimum=0,maximum=1)
    pi=probability(target_power,'target_power');alpha=error_probability(alpha)
    k=binomial_cutoff(n,percentile_tail(percentile),alpha)
    if k is None or not _power_at_least(np.array([k]),np.array([n]),q,pi)[0]:
        return None
    lo,hi=0,k
    while lo<hi:
        mid=(lo+hi+1)//2
        if _power_at_least(np.array([k-mid]),np.array([n]),q,pi)[0]:lo=mid
        else:hi=mid-1
    return lo


def diagnose_upper(reference, upper_limit, *, future_n=None, percentile=.99,
                   target_power=.8, alpha=.025, delta=.05, drift=0,
                   unmodelled_share=0., max_future_n=1000000, future_ambiguous_budget=0):
    result = _base(reference, 'diagnose-upper')
    c = real(upper_limit, 'upper_limit', minimum=0, open_min=True)
    epsilon = percentile_tail(percentile)
    alpha = error_probability(alpha)
    delta = error_probability(delta, 'delta')
    pi = probability(target_power, 'target_power')
    if pi <= alpha:
        raise ValidationError('diagnostic target_power must exceed alpha')
    cap = integer(max_future_n, 'max_future_n', 1)
    if cap > 10000000:
        raise ValidationError('max_future_n exceeds diagnostic limit 10000000')
    n = None if future_n is None else integer(future_n, 'future_n', 1)
    budget = integer(future_ambiguous_budget, 'future_ambiguous_budget')
    d = _drift(drift)
    eta = real(unmodelled_share, 'unmodelled_share', minimum=0, maximum=1, open_max=True)
    t = _shift_threshold(c, d)
    a, m = reference.tail_counts(t)
    # The shared CP routine already verifies its upward direction.
    u = cp_upper(m, reference.n, delta)
    w = _worst_tail_upper(eta, u)
    windows = future_n_intervals(w, percentile=percentile, target_power=pi,
                                 alpha=alpha, max_future_n=cap, future_ambiguous_budget=budget)
    blockers = [] if reference.iid else ['IID_NOT_DECLARED']
    if eta >= epsilon:
        feasibility = 'NO_N_POSSIBLE_IN_UNCERTAINTY_CLASS'
    elif w >= epsilon:
        feasibility = 'NO_N_CERTIFIABLE_FROM_CURRENT_REFERENCE'
    else:
        feasibility = 'EVENTUALLY_CERTIFIABLE'
    if feasibility != 'EVENTUALLY_CERTIFIABLE':
        blockers.append(feasibility)
    elif not windows:
        blockers.append('NO_ADMISSIBLE_N_WITHIN_SEARCH_CAP')
    eventual_min = None
    if eta < epsilon:
        qe_inf = _effective_target(epsilon, eta)
        eventual_min = minimum_reference(qe_inf, delta)
        while eventual_min is not None:
            uu=cp_upper(0,eventual_min,delta)
            if _worst_tail_upper(eta,uu)<epsilon:break
            eventual_min=eventual_min+1 if eventual_min<2**53-1 else None
    current = None
    if n is not None:
        qs = tail_target(n, percentile, pi, alpha, future_ambiguous_budget=budget)
        cut = binomial_cutoff(n, epsilon, alpha)
        meets = cut is not None and cut>=budget and bool(_power_at_least(np.array([cut-budget]), np.array([n]), w, pi)[0])
        current = {'future_n':n, 'q_star':qs, 'q_eff':None,
                   'power_lower_bound':certification_power_lower(n,w,percentile,alpha,future_ambiguous_budget=budget),
                   'maximum_future_ambiguous_budget':maximum_future_ambiguous_budget(n,w,percentile=percentile,target_power=pi,alpha=alpha),
                   'reference_cutoff':None, 'minimum_reference_zero_tail':None,
                   'meets_target':meets, 'repair_status':'NO_REFERENCE_CUTOFF',
                   'minimum_ambiguous_resolutions_below_threshold':None}
        if qs is None:
            blockers.append('FUTURE_N_TOO_SMALL' if cut is None else 'FUTURE_OBSERVABILITY_BUDGET_EXCEEDS_CUTOFF')
        elif eta > qs:
            blockers.append('STRUCTURALLY_IMPOSSIBLE_AT_CURRENT_N')
        elif eta == qs:
            blockers.append('ZERO_TAIL_NOT_CERTIFIABLE_AT_CURRENT_N')
        else:
            qe = _effective_target(qs, eta)
            k = binomial_cutoff(reference.n, qe, delta) if qe > 0 else None
            current.update(q_eff=qe, reference_cutoff=k,
                           minimum_reference_zero_tail=minimum_reference(qe, delta))
            if k is None:
                blockers.append('REFERENCE_N_TOO_SMALL')
            else:
                if m <= k:
                    current['repair_status'] = 'ALREADY_MEETS_COUNT_TARGET'
                elif a > k:
                    current['repair_status'] = 'DEFINITE_EXCEEDANCES_EXCEED_CUTOFF'
                    blockers.append('DEFINITE_EXCEEDANCES_EXCEED_CUTOFF')
                else:
                    current['repair_status'] = 'COMPATIBLE_REFINEMENT_COULD_MEET_TARGET'
                    current['minimum_ambiguous_resolutions_below_threshold'] = m-k
                    blockers.append('AMBIGUOUS_OUTCOMES_EXCEED_CUTOFF')
                unbounded = sum(row.count for row in reference.intervals if math.isinf(row.upper))
                if unbounded > k:
                    blockers.append('UNBOUNDED_REFERENCE_BLOCKS_ALL_FINITE_CANDIDATES_AT_CURRENT_N')
        if not meets:
            blockers.append('CURRENT_POWER_TARGET_NOT_CERTIFIED')
    result.update(status='DIAGNOSED' if reference.iid else 'IID_NOT_DECLARED',
                  conditional_only=not reference.iid, upper_limit=c, percentile=percentile,
                  alpha=alpha, delta=delta, drift=d, unmodelled_share=eta, target_power=pi,
                  future_ambiguous_budget=budget,
                  future_observation_contract='At most b unresolved exceedance indicators; exceeding b is outside the guaranteed observation class.',
                  guarantee_scope='fixed_upper_diagnostic',
                  obstructions=list(dict.fromkeys(blockers)),
                  conditional_analysis={'reference_threshold':t, 'definite_exceedances':a,
                    'possible_exceedances':m, 'reference_tail_upper':u, 'worst_tail_upper':w,
                    'feasibility':feasibility, 'population_inadequacy_proved':False,
                    'minimum_reference_zero_tail_eventual':eventual_min,
                    'future_n_intervals':windows,
                    'first_admissible_future_n':windows[0][0] if windows else None,
                    'complete_through_n':cap, 'beyond_cap_evaluated':False,
                    'current_scenario':current},
                  limitations=['Pointwise diagnostic for an independently specified upper limit and drift. Broader post-data scenario selection needs simultaneous protection; fixed-scenario nested threshold selection is provided by design_upper.',
                               'Reference n is fixed; collecting until significance requires a different inference protocol.',
                               'Compatible interval refinement is not a guarantee of a repeat assay outcome.',
                               'Future ambiguity budget is prespecified; it is not an estimated missingness rate.',
                               'No bound on error conditional on a recommendation being issued.'])
    return result
