import math
from fractions import Fraction
from math import comb
import pytest

def test_reference_rejects_silent_count_casts():
    from panel_mic import Reference, ValidationError
    for counts in ([.9,200.9,0],[True,2,0],[-1,3,0],[1.0,2,0],['1',2,0]):
        with pytest.raises(ValidationError): Reference([1,2],counts,iid=True)

def test_reference_endpoints():
    from panel_mic import Reference
    r=Reference([1,2],[5,3,2],iid=True)
    assert r.n==10
    assert r.tail_counts(1)==(5,5)
    assert r.tail_counts(2)==(2,2)
    assert r.tail_counts(1.5)==(2,5)
    assert r.tail_counts(4)==(0,2)
    z=Reference([.1,.3,2.7],[1,2,3,4],iid=True)
    assert z.thresholds==(.1,.3,2.7)
    assert z.tail_counts(.3)==(7,7)

@pytest.mark.parametrize('t,c',[([1,1],[1,2,3]),([2,1],[1,2,3]),([0],[1,2]),([math.nan],[1,2]),([1],[0,0]),([1],[1])])
def test_bad_reference(t,c):
    from panel_mic import Reference,ValidationError
    with pytest.raises(ValidationError): Reference(t,c)

def test_missing_outcomes_retained():
    from panel_mic import Reference
    r=Reference.from_intervals([{'lower':None,'upper':None,'count':3},{'lower':1,'upper':2,'count':7}],iid=True)
    assert r.n==10 and r.tail_counts(2)==(0,3)

def test_exact_cutoffs():
    from panel_mic.numerics import binomial_cutoff
    for n,q,k in [(367,.01,None),(368,.01,0),(1000,.01,3),(3687,.001,None),(3688,.001,0)]:
        assert binomial_cutoff(n,q,.025)==k
    assert binomial_cutoff(1,.75,.25)==0
    assert binomial_cutoff(2,.75,.4375)==1
    for k,expected in [(3,True),(4,False)]:
        v=sum(Fraction(comb(1000,j)*99**(1000-j),100**1000) for j in range(k+1))
        assert (v<=Fraction(1,40))==expected

def test_inverse_power():
    from panel_mic.numerics import tail_target,minimum_reference,certification_power
    q=tail_target(1000)
    assert q==pytest.approx(.002297595743084674,abs=1e-15)
    assert minimum_reference(q)==1303
    assert minimum_reference((q-.001)/.999)==2305
    assert minimum_reference((q-.002)/.998)==10045
    assert certification_power(1000,.005)==pytest.approx(.2643224034591385)
    assert certification_power(100,0)==0 and tail_target(100) is None

def test_design_and_simultaneous_differ():
    from panel_mic import Reference,design_upper,frontier_upper
    r=Reference([1,2,4],[1290,8,5,0],iid=True)
    d=design_upper(r,[1,2,4,8],future_n=1000)
    assert d['status']=='FOUND' and d['upper_limit']==4
    assert d['certified_power_floor']==.8 and d['guarantee_scope']=='fixed_scenario'
    assert 'selected_pointwise_upper_bound' not in d
    f=frontier_upper(r,[1,2,4,8],future_ns=[1000],grid_prespecified=True)
    assert f['guarantee_scope']=='simultaneous_frontier' and f['threshold_count']==4
    assert f['rows'][0]['status']!='FOUND'

def test_assess_small_n_and_location():
    from panel_mic import Reference,assess_upper
    a=assess_upper(Reference([1,2],[0,368,0],iid=True),2)
    assert a['status']=='PANEL_REACHES_QUANTILE'
    assert a['location']['upper']==2 and a['location']['lower']==1
    assert assess_upper(Reference([1,2],[0,18,2],iid=True),2)['status']=='PANEL_TOO_NARROW'

def test_iid_required():
    from panel_mic import Reference,assess_upper,design_upper
    r=Reference([1,2],[0,2000,0])
    assert assess_upper(r,2)['status']=='IID_NOT_DECLARED'
    assert design_upper(r,[1,2],future_n=1000)['status']=='IID_NOT_DECLARED'

def test_structural_failure_and_censor_plateau():
    from panel_mic import Reference,design_upper
    from panel_mic.numerics import tail_target
    r=Reference([1,2],[0,2000,0],iid=True)
    assert design_upper(r,[1,2],future_n=100)['status']=='FUTURE_N_TOO_SMALL'
    assert design_upper(r,[1,2],future_n=1000,unmodelled_share=.01)['status']=='STRUCTURALLY_IMPOSSIBLE'
    assert design_upper(r,[1,2],future_n=1000,unmodelled_share=tail_target(1000))['status']=='ZERO_TAIL_NOT_CERTIFIABLE'
    c=Reference([1,2],[0,1990,10],iid=True)
    assert design_upper(c,[2,4,8,16,1024],future_n=1000)['status']=='REFERENCE_RANGE_UNRESOLVED'
    assert design_upper(r,[1,2,4],future_n=1000,drift=1)['upper_limit']==4
    assert design_upper(r,[1,2,4],future_n=1000,unmodelled_share=.001)['status']=='REFERENCE_N_TOO_SMALL'

def test_grid_declaration_required():
    from panel_mic import Reference,frontier_upper
    r=Reference([1,2],[0,2000,0],iid=True)
    assert frontier_upper(r,[1,2],future_ns=[1000])['status']=='GRID_NOT_PRESPECIFIED'

@pytest.mark.parametrize('kwargs',[{'percentile':.9},{'alpha':True},{'delta':math.nan},{'target_power':1},{'drift':-1},{'drift':1.5},{'unmodelled_share':-1},{'future_n':3.5}])
def test_invalid_design(kwargs):
    from panel_mic import Reference,design_upper,ValidationError
    r=Reference([1,2],[0,2000,0],iid=True)
    with pytest.raises(ValidationError): design_upper(r,[1,2],**{'future_n':1000,**kwargs})

def test_cache_cannot_bypass_boolean_validation():
    from panel_mic.numerics import binomial_cutoff,tail_target
    from panel_mic import ValidationError
    binomial_cutoff(1,.01,.025); tail_target(1)
    with pytest.raises(ValidationError): binomial_cutoff(True,.01,.025)
    with pytest.raises(ValidationError): tail_target(True)

def test_one_ulp_probability_boundary_uses_exact_small_sample_comparison():
    from panel_mic.numerics import binomial_cutoff
    assert binomial_cutoff(2,.75,math.nextafter(.4375,0))==0
    assert binomial_cutoff(2,.75,.4375)==1

@pytest.mark.parametrize('kwargs',[{'iid':'yes'},{'unit':''},{'label':False},{'metadata':[]},{'metadata':{'original_n':5}}])
def test_reference_declarations_and_denominator(kwargs):
    from panel_mic import Reference,ValidationError
    with pytest.raises(ValidationError):Reference([1],[2,0],**kwargs)

@pytest.mark.parametrize('rows',[None,[{}],[{'lower':2,'upper':1,'count':2}],[{'lower':0,'upper':1,'count':-1}]])
def test_bad_interval_rows(rows):
    from panel_mic import Reference,ValidationError
    with pytest.raises(ValidationError):Reference.from_intervals(rows)

def test_nonsequence_inputs_and_parameter_guards():
    from panel_mic import Reference,ValidationError,design_upper
    from panel_mic.numerics import cp_upper,minimum_reference
    for t,c in [(None,[1,0]),([1],None)]:
        with pytest.raises(ValidationError):Reference(t,c)
    with pytest.raises(ValidationError):cp_upper(3,2,.05)
    with pytest.raises(ValidationError):design_upper(None,[1],future_n=1000)
    with pytest.raises(ValidationError):design_upper(Reference([1],[2000,0]),[1],future_n=1000,drift=1024)
    assert minimum_reference(0) is None

def test_low_power_success_does_not_imply_quantile_adequacy():
    from panel_mic import Reference,design_upper
    from panel_mic.cli import explain
    r=design_upper(Reference([1],[25,0],iid=True),[1],future_n=400,target_power=1e-30)
    assert r['status']=='FOUND' and not r['population_quantile_adequacy_implied']
    assert 'nie implikuje' in explain(r)
