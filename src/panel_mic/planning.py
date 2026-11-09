"""Fixed-protocol assurance and QC-aware future power; no cost optimization."""
from fractions import Fraction

from .errors import ValidationError
from .models import integer, real
from .numerics import (probability, error_probability, percentile_tail,
    binomial_cutoff, tail_target, cp_upper, _rounded, _cdf_bounds)
from ._version import __version__ as VERSION


def _caps(values):
    result=[];seen=set()
    try:
        for item in values:
            if not isinstance(item,(list,tuple)) or len(item)!=2:
                raise ValidationError('each QC cap requires [unresolved_count_threshold, probability_upper]')
            r=integer(item[0],'QC count threshold',1)
            gamma=real(item[1],'QC probability upper',minimum=0,maximum=1)
            if r in seen:raise ValidationError('duplicate QC count threshold')
            seen.add(r);result.append((r,gamma))
    except TypeError as exc:raise ValidationError('QC caps must be a sequence') from exc
    return tuple(sorted(result))


def qc_power_lower(future_n,q,ambiguity_tail_caps,*,percentile=.99,alpha=.025):
    """Sharp floor under Pr(R>=r)<=gamma_r, without K/R independence.

    Bounds must hold uniformly over the declared future population/process.
    A calibrated profile needs its own simultaneous confidence protection.
    """
    n=integer(future_n,'future_n',1)
    if n>10000000:raise ValidationError('future_n exceeds the QC numerical limit')
    q=real(q,'q',minimum=0,maximum=1);caps=_caps(ambiguity_tail_caps)
    k=binomial_cutoff(n,percentile_tail(percentile),error_probability(alpha))
    if k is None:return 0.
    best=Fraction(0)
    for r,gamma in caps:
        lower,_=_cdf_bounds(n,k-r+1,q)
        best=max(best,Fraction(lower)-Fraction(gamma))
    return _rounded(best)


def _allocated(error,m):
    allocated=_rounded(Fraction(error)/m)
    if allocated<=0:raise ValidationError('error allocation is below the supported numerical range')
    return allocated


def _certified_cp_upper(events,batches,error):
    """Compatibility alias for the shared, directed CP implementation."""
    return cp_upper(events,batches,error)


def qc_profile_from_counts(unresolved_counts,*,future_n,thresholds,error=.05,
                          iid_batches=False,same_target_process=False,
                          thresholds_prespecified=False):
    """Simultaneous CP upper bounds from fixed-size independent QC batches.

    Counts describe unresolved indicators for the SAME physical threshold,
    assay protocol, population/process and future batch size. Those properties
    are declarations that require source evidence, not facts inferred here.
    """
    n=integer(future_n,'future_n',1);error=error_probability(error,'QC error')
    try:counts=tuple(integer(x,'unresolved QC count') for x in unresolved_counts)
    except TypeError as exc:raise ValidationError('QC counts must be a sequence') from exc
    if not counts:raise ValidationError('QC batches must be nonempty; unknown batches cannot be dropped')
    if any(x>n for x in counts):raise ValidationError('QC unresolved count exceeds its batch size')
    try:caps=_caps([(r,0.) for r in thresholds])
    except TypeError as exc:raise ValidationError('QC thresholds must be a sequence') from exc
    if not caps or any(r>n for r,_ in caps):raise ValidationError('QC thresholds must lie in [1, future_n]')
    flags=(iid_batches,same_target_process,thresholds_prespecified)
    if any(type(x) is not bool for x in flags):raise ValidationError('QC declarations must be Boolean')
    result={'future_n':n,'batches':len(counts),'QC_error':error,
        'iid_batches_declared':iid_batches,'same_target_process_declared':same_target_process,
        'thresholds_prespecified':thresholds_prespecified,
        'scope':'Conditional on independent fixed-size batches and the same prespecified target process; not proof of transport or biological independence.'}
    if not all(flags):return {**result,'status':'QC_ASSUMPTIONS_NOT_DECLARED','uniform_tail_caps':None}
    level=_allocated(error,len(caps));profile=[]
    for r,_ in caps:
        events=sum(x>=r for x in counts)
        upper=_certified_cp_upper(events,len(counts),level)
        profile.append((r,upper))
    return {**result,'status':'QC_PROFILE_CONDITIONAL_ON_DECLARATIONS',
        'uniform_tail_caps':profile,'per_threshold_error':level,'simultaneous_confidence_at_least':1-error}


def protocol_assurance(reference_n,control_n,*,future_n,possible_tail_reference,
                       possible_tail_control,independent_stages=False,
                       percentile=.99,target_power=.8,alpha=.025,
                       reference_error=.05,control_error=.025,threshold_family_size=1,
                       future_ambiguous_budget=0,ambiguity_tail_caps=None,
                       unmodelled_share=0.):
    """Probability of both count certificates passing at given sample sizes.

    Possible-tail probabilities are prespecified hypothetical observation-law
    values at the transformed reference/control threshold, not estimates here.
    QC caps are treated as given; any calibration confidence is separate.
    """
    nr=integer(reference_n,'reference_n',1);nc=integer(control_n,'control_n',1)
    n=integer(future_n,'future_n',1);b=integer(future_ambiguous_budget,'future_ambiguous_budget')
    if n>10000000:raise ValidationError('future_n exceeds the protocol numerical limit')
    family=integer(threshold_family_size,'threshold_family_size',1)
    pr=real(possible_tail_reference,'possible_tail_reference',minimum=0,maximum=1)
    pc=real(possible_tail_control,'possible_tail_control',minimum=0,maximum=1)
    pi=probability(target_power,'target_power');alpha=error_probability(alpha)
    dr=error_probability(reference_error,'reference_error');dc=error_probability(control_error,'control_error')
    eta=real(unmodelled_share,'unmodelled_share',minimum=0,maximum=1,open_max=True)
    percentile_tail(percentile)
    if type(independent_stages) is not bool:raise ValidationError('independent_stages must be Boolean')
    if ambiguity_tail_caps is not None and b:
        raise ValidationError('supply either a hard ambiguity budget or QC tail caps')
    caps=_caps(([] if b>=n else [(b+1,0.)]) if ambiguity_tail_caps is None else ambiguity_tail_caps)
    result={'schema_version':'panel_mic.protocol/v2','software_version':VERSION,
        'operation':'protocol-assurance','reference_n':nr,'control_n':nc,'future_n':n,
        'percentile':percentile,'target_power':pi,'alpha':alpha,'reference_error':dr,'control_error':dc,
        'threshold_family_size':family,'independent_stages_declared':independent_stages,
        'possible_tail_reference':pr,'possible_tail_control':pc,'ambiguity_tail_caps':caps,
        'future_ambiguous_budget':b if ambiguity_tail_caps is None else None,
        'unmodelled_share':eta,'q_star':None,'q_effective':None,
        'reference_acceptance_count':None,'control_acceptance_count':None,
        'stage_success_probabilities':[0.,0.],'joint_success_probability':0.,'joint_success_interval':[0.,0.],
        'scope':'Fixed hypothetical count-law planning, conditional on stated observations and QC bounds; not biological validation or a sample-allocation optimizer.'}
    targets=[]
    for r,gamma in caps:
        required=Fraction(pi)+Fraction(gamma)
        if required>1:continue
        if required==1 or _rounded(required,True)==1.:
            k=binomial_cutoff(n,percentile_tail(percentile),alpha)
            candidate=0. if k is not None and k>=r-1 else None
        else:
            rounded=_rounded(required,True)
            candidate=tail_target(n,percentile,rounded,alpha,future_ambiguous_budget=r-1)
        if candidate is not None:targets.append(candidate)
    if not targets:return {**result,'status':'POWER_TARGET_UNATTAINABLE_IN_QC_CLASS'}
    qs=max(targets);result['q_star']=qs
    if eta>qs:return {**result,'status':'POPULATION_CLASS_UNATTAINABLE'}
    qe=max(0.,_rounded((Fraction(qs)-Fraction(eta))/(1-Fraction(eta))))
    result['q_effective']=qe
    kr=binomial_cutoff(nr,qe,_allocated(dr,family)) if qe else None
    kc=binomial_cutoff(nc,qe,_allocated(dc,family)) if qe else None
    result.update(reference_acceptance_count=kr,control_acceptance_count=kc)
    # An infeasible control certificate does not erase the reference marginal
    # (and conversely). Only joint success must then be zero.
    lr,ur=(0.,0.) if kr is None else _cdf_bounds(nr,kr,pr)
    lc,uc=(0.,0.) if kc is None else _cdf_bounds(nc,kc,pc)
    mr=(lr+ur)/2;mc=(lc+uc)/2
    result['stage_success_probabilities']=[mr,mc]
    if kr is None or kc is None:return {**result,'status':'REFERENCE_OR_CONTROL_INFORMATION_LIMIT'}
    if independent_stages:
        lo=Fraction(lr)*Fraction(lc);hi=Fraction(ur)*Fraction(uc)
        result['joint_success_probability']=mr*mc
    else:
        lo=max(Fraction(0),Fraction(lr)+Fraction(lc)-1);hi=min(Fraction(ur),Fraction(uc))
        result['joint_success_probability']=None
    result['joint_success_interval']=[_rounded(lo),_rounded(hi,True)]
    return {**result,'status':'FIXED_PROTOCOL_SCENARIO_EVALUATED'}


def protocol_assurance_grid(future_ns,*,minimum_joint_success=.8,**scenario):
    """Evaluate declared Ns; neither allocate samples nor fill unevaluated gaps."""
    goal=probability(minimum_joint_success,'minimum_joint_success')
    try:ns=tuple(integer(n,'future_n',1) for n in future_ns)
    except TypeError as exc:raise ValidationError('future_ns must be a sequence') from exc
    if not ns or len(set(ns))!=len(ns):raise ValidationError('future_ns must be nonempty and distinct')
    if 'future_n' in scenario:raise ValidationError('supply future_ns, not a second future_n')
    ns=sorted(ns);rows=[];windows=[]
    for n in ns:
        row=protocol_assurance(future_n=n,**scenario);rows.append(row)
        if row['joint_success_interval'][0]<goal:continue
        if windows and windows[-1][1]+1==n:windows[-1][1]=n
        else:windows.append([n,n])
    return {'rows':rows,'minimum_joint_success':goal,'evaluated_ns':ns,
        'admissible_future_n_intervals':windows,'complete_integer_scan':len(ns)==ns[-1]-ns[0]+1,
        'scope':'Only the declared Ns and fixed stage sample sizes were evaluated; no monotone minimum or sample-allocation optimum is asserted.'}
