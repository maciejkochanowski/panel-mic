"""Upper range decisions and explicitly scoped design guarantees."""
import math
from fractions import Fraction
from scipy.stats import binom
from .errors import ValidationError
from .models import Reference, integer, real, grid
from .numerics import (percentile_tail,error_probability,probability,binomial_cutoff,
                       certification_power_lower,tail_target,minimum_reference,cp_upper,cp_lower,excess_tail_significant,_power_ge)

from ._version import __version__ as VERSION

def _reference(reference):
    if not isinstance(reference,Reference): raise ValidationError('reference must be a Reference')

def _drift(drift):
    drift=integer(drift,'drift')
    if drift>1023: raise ValidationError('drift exceeds the supported floating-point exponent range')
    return drift

def _shift_threshold(concentration,drift):
    """Conservative division including rounded subnormal results."""
    threshold=math.ldexp(concentration,-drift)
    try: rounded_up=math.ldexp(threshold,drift)>concentration
    except OverflowError: rounded_up=True
    return math.nextafter(threshold,0.) if rounded_up else threshold

def _effective_target(qstar,eta):
    exact=(Fraction(qstar)-Fraction(eta))/(1-Fraction(eta))
    value=float(exact)
    return math.nextafter(value,0.) if Fraction(value)>exact else value

def _worst_tail_upper(eta,upper):
    exact=Fraction(eta)+(1-Fraction(eta))*Fraction(upper)
    value=float(exact)
    return math.nextafter(value,1.) if Fraction(value)<exact else value

def _base(reference,operation):
    _reference(reference)
    return {'schema_version':'panel_mic.upper/v2','software_version':VERSION,'operation':operation,
            'reference_n':reference.n,'reference_label':reference.label,'unit':reference.unit,
            'metadata':reference.metadata,'iid_declared':reference.iid}

def assess_upper(reference,upper_limit,*,percentile=.99,alpha=.025):
    result=_base(reference,'assess-upper'); c=real(upper_limit,'upper_limit',minimum=0,open_min=True)
    epsilon=percentile_tail(percentile); alpha=error_probability(alpha)
    a,b=reference.tail_counts(c)
    result.update(upper_limit=c,percentile=percentile,alpha=alpha,definite_exceedances=a,possible_exceedances=b,
                  empirical_tail_interval=[a/reference.n,b/reference.n],guarantee_scope='fixed_upper_assessment')
    if not reference.iid: return {**result,'status':'IID_NOT_DECLARED','location':None,'tail_bounds':None}
    cutoff=binomial_cutoff(reference.n,epsilon,alpha)
    good=cutoff is not None and b<=cutoff
    bad=excess_tail_significant(a,reference.n,epsilon,alpha)
    lower=0.; upper=None
    for t in reference.thresholds:
        aa,bb=reference.tail_counts(t)
        if excess_tail_significant(aa,reference.n,epsilon,alpha): lower=t
        if upper is None and cutoff is not None and bb<=cutoff: upper=t
    result.update(status='PANEL_REACHES_QUANTILE' if good else 'PANEL_TOO_NARROW' if bad else 'UNDECIDED',
                  tail_bounds=[cp_lower(a,reference.n,alpha),cp_upper(b,reference.n,alpha)],
                  location={'lower':lower,'upper':upper,'lower_open':True,'upper_closed':True,'confidence':1-2*alpha,
                            'note':'Population quantile interval; null upper is unbounded. No exact continuous MIC is imputed.'})
    return result

def design_upper(reference,candidates,*,future_n,percentile=.99,target_power=.8,alpha=.025,delta=.05,drift=0,unmodelled_share=0.,future_ambiguous_budget=0):
    result=_base(reference,'design-upper'); candidates=grid(candidates)
    future_n=integer(future_n,'future_n',1); percentile_tail(percentile)
    pi=probability(target_power,'target_power'); alpha=error_probability(alpha); delta=error_probability(delta,'delta')
    drift=_drift(drift); eta=real(unmodelled_share,'unmodelled_share',minimum=0,maximum=1,open_max=True)
    budget=integer(future_ambiguous_budget,'future_ambiguous_budget')
    qs=tail_target(future_n,percentile,pi,alpha,future_ambiguous_budget=budget)
    result.update(future_n=future_n,percentile=percentile,target_power=pi,alpha=alpha,delta=delta,drift=drift,
                  unmodelled_share=eta,candidates=list(candidates),q_star=qs,q_eff=None,upper_limit=None,
                  reference_cutoff=None,minimum_reference_zero_tail=None,certified_power_floor=None,
                  guarantee_scope='fixed_scenario',candidate_diagnostics=[],population_quantile_adequacy_implied=False,
                  future_ambiguous_budget=budget)
    if not reference.iid: return {**result,'status':'IID_NOT_DECLARED'}
    if qs is None:
        cutoff=binomial_cutoff(future_n,percentile_tail(percentile),alpha)
        return {**result,'status':'FUTURE_N_TOO_SMALL' if cutoff is None else 'FUTURE_OBSERVABILITY_BUDGET_EXCEEDS_CUTOFF'}
    if eta>qs: return {**result,'status':'STRUCTURALLY_IMPOSSIBLE'}
    if eta==qs: return {**result,'status':'ZERO_TAIL_NOT_CERTIFIABLE'}
    qe=_effective_target(qs,eta)
    result['q_eff']=qe; result['minimum_reference_zero_tail']=minimum_reference(qe,delta)
    if qe<=0: return {**result,'status':'ZERO_TAIL_NOT_CERTIFIABLE'}
    k=binomial_cutoff(reference.n,qe,delta); result['reference_cutoff']=k
    if k is None: return {**result,'status':'REFERENCE_N_TOO_SMALL'}
    for c in candidates:
        t=_shift_threshold(c,drift); a,b=reference.tail_counts(t)
        result['candidate_diagnostics'].append({'upper_limit':c,'reference_threshold':t,'definite_exceedances':a,'possible_exceedances':b})
        if result['upper_limit'] is None and b<=k: result['upper_limit']=c
    if result['upper_limit'] is not None:
        return {**result,'status':'FOUND','certified_power_floor':pi,
                'population_quantile_adequacy_implied':qs<=percentile_tail(percentile)}
    unresolved=sum(r.count for r in reference.intervals if math.isinf(r.upper))
    return {**result,'status':'REFERENCE_RANGE_UNRESOLVED' if unresolved>k else 'NO_PERMITTED_BOUNDARY'}

def frontier_upper(reference,candidates,*,future_ns,drifts=(0,),unmodelled_shares=(0.,),percentile=.99,target_power=.8,alpha=.025,delta=.05,grid_prespecified=False,future_ambiguous_budget=0):
    result=_base(reference,'frontier-upper'); candidates=grid(candidates)
    try:
        ns=tuple(integer(n,'future_n',1) for n in future_ns)
        ds=tuple(_drift(d) for d in drifts)
        es=tuple(real(e,'unmodelled_share',minimum=0,maximum=1,open_max=True) for e in unmodelled_shares)
    except TypeError as exc: raise ValidationError('frontier parameter grids must be sequences') from exc
    if not ns or not ds or not es: raise ValidationError('frontier grids must be nonempty')
    if type(grid_prespecified) is not bool: raise ValidationError('grid_prespecified must be Boolean')
    percentile_tail(percentile); pi=probability(target_power,'target_power')
    alpha=error_probability(alpha); delta=error_probability(delta,'delta')
    budget=integer(future_ambiguous_budget,'future_ambiguous_budget')
    thresholds=sorted({_shift_threshold(c,d) for c in candidates for d in ds})
    result.update(guarantee_scope='simultaneous_frontier',delta=delta,percentile=percentile,alpha=alpha,target_power=pi,
                  threshold_count=len(thresholds),simultaneous_tail_bounds=[],rows=[],grid_prespecified=grid_prespecified,
                  future_ambiguous_budget=budget)
    if not reference.iid: return {**result,'status':'IID_NOT_DECLARED'}
    if not grid_prespecified: return {**result,'status':'GRID_NOT_PRESPECIFIED'}
    per=math.nextafter(delta/len(thresholds),0.)
    bounds={}
    for t in thresholds:
        a,b=reference.tail_counts(t); u=cp_upper(b,reference.n,per); bounds[t]=u
        result['simultaneous_tail_bounds'].append({'reference_threshold':t,'possible_exceedances':b,'upper_tail_bound':u,'error_allocation':per})
    for n in ns:
        qs=tail_target(n,percentile,pi,alpha,future_ambiguous_budget=budget)
        for d in ds:
            for e in es:
                row={'future_n':n,'drift':d,'unmodelled_share':e,'q_star':qs,'upper_limit':None,'certified_power_floor':None,'candidate_bounds':[],'status':'NO_CERTIFIED_BOUNDARY'}
                row['population_quantile_adequacy_implied']=False
                if qs is None:
                    row['status']='FUTURE_N_TOO_SMALL' if binomial_cutoff(n,percentile_tail(percentile),alpha) is None else 'FUTURE_OBSERVABILITY_BUDGET_EXCEEDS_CUTOFF'
                elif e>qs: row['status']='STRUCTURALLY_IMPOSSIBLE'
                elif e==qs: row['status']='ZERO_TAIL_NOT_CERTIFIABLE'
                for c in candidates:
                    u=bounds[_shift_threshold(c,d)]; worst=_worst_tail_upper(e,u)
                    power=certification_power_lower(n,worst,percentile,alpha,future_ambiguous_budget=budget)
                    from .diagnostics import maximum_future_ambiguous_budget
                    maxb=maximum_future_ambiguous_budget(n,worst,percentile=percentile,target_power=pi,alpha=alpha)
                    row['candidate_bounds'].append({'upper_limit':c,'worst_tail_upper_bound':worst,'power_lower_bound':power,
                                                    'maximum_future_ambiguous_budget':maxb})
                    cutoff=binomial_cutoff(n,percentile_tail(percentile),alpha)
                    if row['upper_limit'] is None and cutoff is not None and _power_ge(n,cutoff-budget,worst,pi):
                        row.update(upper_limit=c,certified_power_floor=pi,status='FOUND',
                                   population_quantile_adequacy_implied=worst<=percentile_tail(percentile))
                result['rows'].append(row)
    return {**result,'status':'COMPUTED'}
