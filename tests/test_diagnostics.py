import json
import math
import pytest
from scipy.stats import binom
import panel_mic
from panel_mic import Reference, ValidationError


def diagnose(*args, **kwargs):
    assert hasattr(panel_mic, 'diagnose_upper'), 'approved diagnostic API is not implemented'
    return panel_mic.diagnose_upper(*args, **kwargs)


def zero_ref(n=445, iid=True):
    return Reference([1, 2], [0, n, 0], iid=iid)


def test_complete_future_windows_keep_sawtooth_gaps():
    r = diagnose(zero_ref(), 2, future_n=6190, max_future_n=20000)
    d = r['conditional_analysis']
    assert r['status'] == 'DIAGNOSED'
    assert d['future_n_intervals'][:2] == [[6129, 6136], [6243, 6276]]
    assert d['future_n_intervals'][-1] == [6698, 20000]
    assert all(b+1<c for (a,b),(c,d) in zip(d['future_n_intervals'],d['future_n_intervals'][1:]))
    assert d['first_admissible_future_n'] == 6129
    assert d['complete_through_n'] == 20000
    assert d['beyond_cap_evaluated'] is False
    assert d['current_scenario']['meets_target'] is False
    assert 'CURRENT_POWER_TARGET_NOT_CERTIFIED' in r['obstructions']


def test_minimum_N_does_not_claim_all_larger_N():
    r = diagnose(zero_ref(), 2, max_future_n=6200)
    windows = r['conditional_analysis']['future_n_intervals']
    assert windows == [[6129, 6136]]
    assert not any(a <= 6190 <= b for a, b in windows)


def test_all_N_classification_is_about_evidence_not_true_population():
    r = diagnose(zero_ref(812), 2, percentile=.999, future_n=10000, max_future_n=1000000)
    d = r['conditional_analysis']
    assert d['future_n_intervals'] == []
    assert d['feasibility'] == 'NO_N_CERTIFIABLE_FROM_CURRENT_REFERENCE'
    assert d['minimum_reference_zero_tail_eventual'] == 2995
    assert d['current_scenario']['minimum_reference_zero_tail'] == 13042
    assert 'REFERENCE_N_TOO_SMALL' in r['obstructions']
    assert d['population_inadequacy_proved'] is False


def test_no_admissible_N_within_cap_is_not_global_impossibility():
    r = diagnose(zero_ref(), 2, max_future_n=5000)
    d = r['conditional_analysis']
    assert d['feasibility'] == 'EVENTUALLY_CERTIFIABLE'
    assert d['future_n_intervals'] == []
    assert 'NO_ADMISSIBLE_N_WITHIN_SEARCH_CAP' in r['obstructions']


def test_IID_gate_keeps_conditional_results_without_certifying():
    r = diagnose(zero_ref(iid=False), 2, max_future_n=20000)
    assert r['status'] == 'IID_NOT_DECLARED'
    assert r['conditional_only'] is True
    assert r['conditional_analysis']['first_admissible_future_n'] == 6129
    assert 'IID_NOT_DECLARED' in r['obstructions']
    json.dumps(r, allow_nan=False)


def test_repair_count_and_definite_failure_are_distinct():
    ambiguous = Reference.from_intervals([{'lower':0,'upper':1,'count':1990},
                                         {'lower':1,'upper':4,'count':10}], iid=True)
    a = diagnose(ambiguous, 2, future_n=1000, max_future_n=1000)
    s = a['conditional_analysis']['current_scenario']
    assert s['reference_cutoff'] is not None
    assert s['minimum_ambiguous_resolutions_below_threshold'] == 10-s['reference_cutoff']
    assert s['repair_status'] == 'COMPATIBLE_REFINEMENT_COULD_MEET_TARGET'
    definite = Reference.from_intervals([{'lower':0,'upper':1,'count':1990},
                                        {'lower':2,'upper':None,'count':10}], iid=True)
    b = diagnose(definite, 2, future_n=1000, max_future_n=1000)
    s = b['conditional_analysis']['current_scenario']
    assert s['minimum_ambiguous_resolutions_below_threshold'] is None
    assert s['repair_status'] == 'DEFINITE_EXCEEDANCES_EXCEED_CUTOFF'
    assert 'UNBOUNDED_REFERENCE_BLOCKS_ALL_FINITE_CANDIDATES_AT_CURRENT_N' in b['obstructions']


def test_all_blockers_retained_and_structural_limit_scoped():
    r = diagnose(zero_ref(100, False), 2, percentile=.999, future_n=10000,
                 unmodelled_share=.0005, max_future_n=50000)
    assert 'IID_NOT_DECLARED' in r['obstructions']
    assert 'STRUCTURALLY_IMPOSSIBLE_AT_CURRENT_N' in r['obstructions']
    assert 'NO_N_CERTIFIABLE_FROM_CURRENT_REFERENCE' in r['obstructions']
    assert r['conditional_analysis']['feasibility'] != 'NO_N_POSSIBLE_IN_UNCERTAINTY_CLASS'
    r = diagnose(zero_ref(5000), 2, unmodelled_share=.01, max_future_n=1000)
    assert r['conditional_analysis']['feasibility'] == 'NO_N_POSSIBLE_IN_UNCERTAINTY_CLASS'


def test_small_cap_and_future_N_are_separately_reported():
    r = diagnose(zero_ref(2000), 2, future_n=100, max_future_n=100)
    assert r['conditional_analysis']['future_n_intervals'] == []
    assert 'FUTURE_N_TOO_SMALL' in r['obstructions']


@pytest.mark.parametrize('kw', [{'target_power':.01},{'max_future_n':True},
                             {'max_future_n':10000001},{'max_future_n':0},
                             {'future_n':False},{'drift':-1},{'unmodelled_share':1}])
def test_diagnostic_validation(kw):
    with pytest.raises(ValidationError):
        diagnose(zero_ref(), 2, **kw)


def test_diagnostic_windows_match_independent_direct_scan():
    r = diagnose(zero_ref(), 2, max_future_n=20000)
    d = r['conditional_analysis'];q = d['worst_tail_upper']
    windows = d['future_n_intervals']
    for n in range(6000, 7001):
        k = int(binom.ppf(.025, n, .01))
        if binom.cdf(k, n, .01) > .025:k -= 1
        expected = k >= 0 and binom.cdf(k, n, q) >= .8
        assert any(a <= n <= b for a, b in windows) == expected


def test_drift_and_underflow_do_not_lose_exceedances():
    a = diagnose(zero_ref(2000), 2, future_n=1000, drift=1, max_future_n=1000)
    assert a['conditional_analysis']['possible_exceedances'] == 2000
    assert a['conditional_analysis']['future_n_intervals'] == []
    tiny = Reference([math.ulp(0.)], [1,0], iid=True)
    b = diagnose(tiny, math.ulp(0.), drift=1, max_future_n=1000)
    assert b['conditional_analysis']['definite_exceedances'] == 1


def test_cli_operation_and_schema(tmp_path):
    from panel_mic.cli import run_operation, explain, run_batch
    data={'reference':{'thresholds':[1,2],'counts':[0,445,0],'iid':True},
          'upper_limit':2,'max_future_n':6200}
    r=run_operation('diagnose-upper',data)
    assert r['conditional_analysis']['future_n_intervals']==[[6129,6136]]
    assert '6129' in explain(r) and '6136' in explain(r)
    batch=run_batch({'cohorts':[{'id':'a','operation':'diagnose-upper',**data}]})
    assert batch['cohorts_invalid']==0
    import jsonschema
    from pathlib import Path
    schema=json.loads((Path(panel_mic.__file__).parent/'schemas/input.json').read_text())
    jsonschema.validate(data,schema)


def test_numeric_boundary_guards_and_zero_tail():
    import numpy as np
    from panel_mic.diagnostics import _null_accepts, _power_at_least, future_n_intervals
    assert _null_accepts(np.array([1]),np.array([2]),.75,.4375).tolist()==[True]
    # The exact dyadic equality is now resolved, rather than conservatively withheld.
    assert _power_at_least(np.array([0]),np.array([2]),.5,.25).tolist()==[True]
    assert _power_at_least(np.array([0]),np.array([2]),.5,.25-1e-13).tolist()==[True]
    assert _power_at_least(np.array([0]),np.array([2]),.5,.25+1e-13).tolist()==[False]
    assert future_n_intervals(0,max_future_n=400)==[[368,400]]
    assert future_n_intervals(0,max_future_n=100)==[]
    with pytest.raises(ValidationError):future_n_intervals(.001,target_power=.01)
    with pytest.raises(ValidationError):future_n_intervals(.001,max_future_n=10000001)


def test_current_count_target_and_zero_tail_equality():
    from panel_mic.numerics import tail_target
    r=diagnose(zero_ref(2000),2,future_n=1000,max_future_n=1000)
    assert r['conditional_analysis']['current_scenario']['meets_target']
    assert r['conditional_analysis']['current_scenario']['repair_status']=='ALREADY_MEETS_COUNT_TARGET'
    e=diagnose(zero_ref(2000),2,future_n=1000,unmodelled_share=tail_target(1000),max_future_n=1000)
    assert 'ZERO_TAIL_NOT_CERTIFIABLE_AT_CURRENT_N' in e['obstructions']
    # Extreme valid eta leaves a mathematical reference minimum outside the
    # supported exact-integer range; this is retained as unknown, never zero.
    z=diagnose(zero_ref(),2,unmodelled_share=math.nextafter(.01,0),max_future_n=1000)
    assert z['conditional_analysis']['minimum_reference_zero_tail_eventual'] is None


def test_zero_tail_window_respects_conservative_large_n_cutoff_guard():
    from panel_mic.diagnostics import future_n_intervals
    from panel_mic.numerics import certification_power
    n=10001;alpha=math.exp(n*math.log1p(-.01))
    intervals=future_n_intervals(0,alpha=alpha,max_future_n=10004)
    assert intervals==[[10002,10004]]
    assert certification_power(10001,0,alpha=alpha)==0


def test_eventual_minimum_uses_same_directed_bound_as_diagnostic():
    delta=.4948386596002073
    r=diagnose(zero_ref(70),2,delta=delta,max_future_n=400)
    d=r['conditional_analysis']
    assert d['feasibility']=='NO_N_CERTIFIABLE_FROM_CURRENT_REFERENCE'
    assert d['minimum_reference_zero_tail_eventual']==71
    from panel_mic.diagnostics import _power_at_least
    import numpy as np
    assert _power_at_least(np.array([1]),np.array([3]),.5,.5).tolist()==[True]
