"""Independent count laws catch spurious assurance and missingness independence."""
from fractions import Fraction
import math

import pytest
import panel_mic
from panel_mic.cli import run_operation, explain


def cdf(n, k, q):
    q=Fraction(q)
    return sum((math.comb(n,j)*q**j*(1-q)**(n-j) for j in range(k+1)),Fraction(0))


def scenario():
    return dict(reference_n=1139,control_n=1402,future_n=874,
                possible_tail_reference=.001,possible_tail_control=.001,
                independent_stages=True)


def test_protocol_checks_joint_success_rather_than_zero_tail_minimum():
    r=panel_mic.protocol_assurance(**scenario())
    assert r['reference_acceptance_count']==r['control_acceptance_count']==0
    assert r['joint_success_probability']==pytest.approx(.07868747006270206)
    assert r['joint_success_interval'][0]<=float((1-Fraction(.001))**2541)
    assert r['future_n']==874


def test_unknown_stage_dependence_returns_frechet_interval():
    args=scenario();args['independent_stages']=False
    r=panel_mic.protocol_assurance(**args)
    assert r['joint_success_probability'] is None
    assert r['joint_success_interval']==pytest.approx([0.,.24593168922618386])


def test_larger_fixed_protocol_has_high_assurance_without_optimizing_costs():
    args=scenario();args.update(reference_n=10000,control_n=10000)
    assert panel_mic.protocol_assurance(**args)['joint_success_probability']==pytest.approx(.9591744821003543)


def test_qc_overrun_bound_handles_dependence_on_exceedance_count():
    p=panel_mic.qc_power_lower(874,.001,[(2,.05)])
    exact=cdf(874,2,.001)-Fraction(.05)
    assert Fraction(p)<=exact
    assert float(exact)-p<2e-15
    assert panel_mic.qc_power_lower(874,0.,[(1,.05)])==pytest.approx(.95)
    assert panel_mic.qc_power_lower(874,.001,[])==0.
    assert panel_mic.qc_power_lower(1000,0.,[(5,0.)])==0.


def test_joint_qc_profile_uses_best_protected_budget():
    caps=[(1,.2),(2,.06),(3,.02),(4,.001)]
    exact=max(Fraction(0),*(cdf(874,4-r,.001)-Fraction(g) for r,g in caps))
    bound=panel_mic.qc_power_lower(874,.001,caps)
    assert Fraction(bound)<=exact
    assert float(exact)-bound<2e-15


def test_uncertain_qc_profile_cannot_authorize_complete_data_assurance():
    args=scenario();args['ambiguity_tail_caps']=[(1,.25)]
    r=panel_mic.protocol_assurance(**args)
    assert r['status']=='POWER_TARGET_UNATTAINABLE_IN_QC_CLASS'
    assert r['joint_success_probability']==0.


def test_qc_calibration_requires_independent_same_process_prespecified_batches():
    args=dict(future_n=874,thresholds=[1,2],error=.05)
    withheld=panel_mic.qc_profile_from_counts([0]*30,**args)
    assert withheld['status']=='QC_ASSUMPTIONS_NOT_DECLARED'
    assert withheld['uniform_tail_caps'] is None
    r=panel_mic.qc_profile_from_counts([0]*30,**args,iid_batches=True,
        same_target_process=True,thresholds_prespecified=True)
    expected=1-.025**(1/30)
    assert r['status']=='QC_PROFILE_CONDITIONAL_ON_DECLARATIONS'
    assert [r for r,_ in r['uniform_tail_caps']]==[1,2]
    assert [p for _,p in r['uniform_tail_caps']]==pytest.approx([expected,expected])
    assert r['batches']==30


@pytest.mark.parametrize('bad',[None,True,-1,875,1.5])
def test_qc_never_drops_unknown_or_invalid_batches(bad):
    with pytest.raises(panel_mic.ValidationError):
        panel_mic.qc_profile_from_counts([0,bad],future_n=874,thresholds=[1])


@pytest.mark.parametrize('caps',[[(0,.1)],[(1,-.1)],[(1,True)],[(1,1.1)],[(True,.1)],[(1,.1),(1,.2)]])
def test_invalid_qc_caps_are_rejected(caps):
    with pytest.raises(panel_mic.ValidationError):
        panel_mic.qc_power_lower(874,.001,caps)


def test_protocol_cli_explains_assumptions_and_preserves_fractional_counts():
    r=run_operation('protocol-assurance',scenario())
    text=explain(r)
    assert '7.9%' in text
    assert '1139' in text and '1402' in text
    bad=scenario();bad['reference_n']=1139.5
    with pytest.raises(panel_mic.ValidationError):run_operation('protocol-assurance',bad)


def test_protocol_scenario_has_an_installed_machine_readable_contract():
    import json
    import jsonschema
    from importlib.resources import files
    schema=json.loads(files('panel_mic').joinpath('schemas/protocol.json').read_text())
    jsonschema.validate(scenario(),schema)
    bad=scenario();bad['reference_n']=True
    with pytest.raises(jsonschema.ValidationError):jsonschema.validate(bad,schema)


def test_protocol_grid_preserves_a_real_downward_step_with_more_future_data():
    args=scenario();args.pop('future_n');args.update(reference_n=1499,control_n=10000,
        possible_tail_reference=0.,possible_tail_control=0.)
    r=panel_mic.protocol_assurance_grid([769,770],**args)
    assert [x['joint_success_probability'] for x in r['rows']]==[1.,0.]
    assert r['admissible_future_n_intervals']==[[769,769]]
    assert r['complete_integer_scan'] is True


def test_protocol_grid_never_fills_unexamined_gaps():
    args=scenario();args.pop('future_n');args.update(reference_n=10000,control_n=10000,
        possible_tail_reference=0.,possible_tail_control=0.)
    r=panel_mic.protocol_assurance_grid([719,720,770,874],**args)
    assert r['admissible_future_n_intervals']==[[719,720],[770,770],[874,874]]
    assert r['complete_integer_scan'] is False


def test_huge_valid_budget_remains_an_explicit_zero_power_class():
    args=scenario();args['future_ambiguous_budget']=2**53-1
    r=panel_mic.protocol_assurance(**args)
    assert r['future_ambiguous_budget']==2**53-1
    assert r['joint_success_probability']==0.


def test_invalid_qc_threshold_sequence_is_a_validation_error():
    with pytest.raises(panel_mic.ValidationError):
        panel_mic.qc_profile_from_counts([0],future_n=10,thresholds=None)


@pytest.mark.parametrize('batches,events,error',[(4,1,.05),(8,0,.05),(10,1,.025)])
def test_qc_profile_upper_is_certified_beyond_one_ulp_of_beta_quantile(batches,events,error):
    from fractions import Fraction as F
    import math
    result=panel_mic.qc_profile_from_counts([1]*events+[0]*(batches-events),
        future_n=874,thresholds=[1],error=error,iid_batches=True,
        same_target_process=True,thresholds_prespecified=True)
    q=F(result['uniform_tail_caps'][0][1])
    exact=sum(F(math.comb(batches,j))*q**j*(1-q)**(batches-j) for j in range(events+1))
    assert exact<=F(error)


def test_rounded_unit_power_does_not_invent_a_feasible_count_threshold():
    import math
    args=scenario();args.update(future_n=100,target_power=.75,
        ambiguity_tail_caps=[(1,math.nextafter(.25,0.))])
    r=panel_mic.protocol_assurance(**args)
    assert r['status']=='POWER_TARGET_UNATTAINABLE_IN_QC_CLASS'
    assert r['q_star'] is None


def test_population_contamination_uses_a_conservative_transformed_target():
    from fractions import Fraction as F
    args=scenario();args.update(reference_n=10000,control_n=10000,unmodelled_share=.001)
    r=panel_mic.protocol_assurance(**args)
    expected=(F(r['q_star'])-F(.001))/(1-F(.001))
    assert 0<F(r['q_effective'])<=expected
    assert r['q_effective']<r['q_star']
    args['unmodelled_share']=.1
    assert panel_mic.protocol_assurance(**args)['status']=='POPULATION_CLASS_UNATTAINABLE'


def test_QC_event_observed_in_every_batch_gives_no_small_tail_guarantee():
    r=panel_mic.qc_profile_from_counts([1]*12,future_n=874,thresholds=[1],
        iid_batches=True,same_target_process=True,thresholds_prespecified=True)
    assert r['uniform_tail_caps']==[(1,1.)]


@pytest.mark.parametrize('change',[
    {'future_n':10000001},{'independent_stages':1},
    {'future_ambiguous_budget':1,'ambiguity_tail_caps':[(1,.05)]},
    {'reference_error':5e-324,'threshold_family_size':2},
    {'ambiguity_tail_caps':None,'target_power':True},
])
def test_invalid_planning_contract_is_rejected_before_results(change):
    args=scenario();args.update(change)
    with pytest.raises(panel_mic.ValidationError):panel_mic.protocol_assurance(**args)


@pytest.mark.parametrize('ns,extra',[(None,{}),([],{}),([1,1],{}),([1],{'future_n':1})])
def test_invalid_grid_never_returns_a_partial_success(ns,extra):
    args=scenario();args.pop('future_n');args.update(extra)
    with pytest.raises(panel_mic.ValidationError):panel_mic.protocol_assurance_grid(ns,**args)


def test_qc_input_and_numerical_regimes_are_explicit():
    for caps in (None,[(1,)],[(1,.2,3)]):
        with pytest.raises(panel_mic.ValidationError):panel_mic.qc_power_lower(874,.001,caps)
    with pytest.raises(panel_mic.ValidationError):panel_mic.qc_power_lower(10000001,.001,[(1,.1)])
    for counts,kwargs in [([],{}),([0],{'thresholds':[]}),([0],{'thresholds':[875]}),([0],{'iid_batches':1})]:
        args={'future_n':874,'thresholds':[1]};args.update(kwargs)
        with pytest.raises(panel_mic.ValidationError):panel_mic.qc_profile_from_counts(counts,**args)
    args=scenario();args.update(reference_n=100000000,control_n=100000000)
    with pytest.raises(panel_mic.ValidationError,match='directed numerical regime'):panel_mic.protocol_assurance(**args)
