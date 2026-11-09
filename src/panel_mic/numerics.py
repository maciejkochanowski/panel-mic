"""Exact-distribution formulas with directed checks for reported guarantees.

SciPy supplies inverse seeds and ordinary power evaluations. Confidence bounds,
reported power floors and zero-count sample minima use directed arithmetic.
"""
import math
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from functools import lru_cache
from scipy.stats import binom, beta
from .models import integer, real
from .errors import ValidationError

def probability(value,name):
    return real(value,name,minimum=0,maximum=1,open_min=True,open_max=True)

def percentile_tail(percentile):
    percentile=real(percentile,'percentile')
    if percentile not in (.99,.999):
        raise ValidationError('public percentiles are .99 and .999 only')
    return .01 if percentile == .99 else .001

def error_probability(value,name='alpha'):
    return real(value,name,minimum=0,maximum=.5,open_min=True,open_max=True)

def _exact_cdf(k,n,q):
    """Exact arithmetic on the supplied binary float, only near small-sample ties."""
    if k<0: return Fraction(0)
    if k>=n: return Fraction(1)
    a,b=q.as_integer_ratio(); z=b-a
    term=z**n; total=term
    for j in range(1,k+1):
        term=term*(n-j+1)*a//(j*z)
        total+=term
    return Fraction(total,b**n)

def _cdf_at_most(k,n,q,alpha):
    if k==0:
        return _cdf_relation(n,k,q,alpha) == -1
    value=float(binom.cdf(k,n,q))
    guard=1e-11*alpha
    if abs(value-alpha)>guard: return value<=alpha
    if n<=10000: return _exact_cdf(k,n,q)<=Fraction(alpha)
    return value<alpha-guard  # abstain on an unresolved numerical boundary

def excess_tail_significant(k,n,q,alpha):
    """Strict upper-tail test, with exact comparison near small-sample ties."""
    value=float(binom.sf(k-1,n,q));guard=1e-11*alpha
    if abs(value-alpha)>guard: return value<alpha
    if n<=10000: return 1-_exact_cdf(k-1,n,q)<Fraction(alpha)
    return value<alpha-guard

def binomial_cutoff(n,q,alpha):
    n=integer(n,'n',1); q=probability(q,'q'); alpha=error_probability(alpha)
    return _cutoff(n,q,alpha)

@lru_cache(maxsize=8192)
def _cutoff(n,q,alpha):
    if not _cdf_at_most(0,n,q,alpha):
        return None
    lo,hi=0,n
    while lo<hi:
        mid=(lo+hi+1)//2
        if _cdf_at_most(mid,n,q,alpha): lo=mid
        else: hi=mid-1
    return lo

def certification_power(n,q,percentile=.99,alpha=.025,*,future_ambiguous_budget=0):
    n=integer(n,'n',1); q=real(q,'q',minimum=0,maximum=1)
    b=integer(future_ambiguous_budget,'future_ambiguous_budget')
    k=binomial_cutoff(n,percentile_tail(percentile),alpha)
    return 0. if k is None or k<b else float(binom.cdf(k-b,n,q))


def _decimal_integer_power(base, exponent):
    """Positive integer power using only the active directed multiplications."""
    value=Decimal(1)
    while exponent:
        if exponent & 1:value*=base
        exponent//=2
        if exponent:base*=base
    return value


def _decimal_cdf_bounds(n,k,q):
    """Enclose a positive binomial sum without an absolute probability margin."""
    dq=Decimal.from_float(q)
    with localcontext() as ctx:
        ctx.prec=100;ctx.rounding=ROUND_FLOOR;complement_lower=1-dq
    with localcontext() as ctx:
        ctx.prec=100;ctx.rounding=ROUND_CEILING;complement_upper=1-dq
    bounds=[]
    for rounding,base,denominator in (
        (ROUND_FLOOR,complement_lower,complement_upper),
        (ROUND_CEILING,complement_upper,complement_lower)):
        with localcontext() as ctx:
            ctx.prec=100;ctx.rounding=rounding
            term=_decimal_integer_power(base,n);value=term
            for j in range(1,k+1):
                term=term*Decimal(n-j+1)*dq/Decimal(j)/denominator
                value+=term
            bounds.append(value)
    return bounds[0],min(Decimal(1),bounds[1])


def _rounded(value, upward=False):
    """Convert a Decimal or Fraction to a float in the specified direction."""
    result=float(value)
    actual=Decimal.from_float(result) if isinstance(value,Decimal) else Fraction(result)
    if (upward and actual<value) or (not upward and actual>value):
        result=math.nextafter(result,1. if upward else 0.)
    return result


def _cdf_bounds(n,k,q):
    """Binary64 enclosure, including endpoints; reject unsupported sum sizes."""
    if k<0 or (q==1 and k<n):return 0.,0.
    if k>=n or q==0:return 1.,1.
    if k>100000:raise ValidationError('count cutoff exceeds the directed numerical regime')
    lo,hi=_decimal_cdf_bounds(n,k,q)
    return _rounded(lo),_rounded(hi,True)


def _cdf_relation(n,k,q,level):
    """Return -1 for CDF <= level, +1 for CDF > level, or None if unresolved.

    Equality belongs to the acceptance set. The zero-count version is shared
    by minimum_reference and binomial_cutoff, including n above 10000.
    """
    if k<0 or (q==1 and k<n):return -1
    if k>=n or q==0:return 1 if level<1 else -1
    if k>100000:return None
    lower,upper=_decimal_cdf_bounds(n,k,q)
    exact_level=Decimal.from_float(level)
    if upper<=exact_level:return -1
    if lower>exact_level:return 1
    if n<=10000:
        return -1 if _exact_cdf(k,n,q)<=Fraction(level) else 1
    return None


def certification_power_lower(n,q,percentile=.99,alpha=.025,*,future_ambiguous_budget=0):
    """Directed lower enclosure of certification power, not a rounded estimate."""
    n=integer(n,'n',1); q=real(q,'q',minimum=0,maximum=1)
    b=integer(future_ambiguous_budget,'future_ambiguous_budget')
    k=binomial_cutoff(n,percentile_tail(percentile),alpha)
    if k is None or k<b:return 0.
    return _cdf_bounds(n,k-b,q)[0]


def _decimal_power_at_least(n,k,q,target):
    """Directed bounds, then rational equality checks; unresolved cases defer."""
    if k<0:return False
    if q==0 or k>=n:return True
    if q==1:return False
    if k>100000:return float(binom.cdf(k,n,q))>target+1e-12*target
    lower,upper=_decimal_cdf_bounds(n,k,q);exact_target=Decimal.from_float(target)
    if lower>=exact_target:return True
    if upper<exact_target:return False
    if n<=10000:return _exact_cdf(k,n,q)>=Fraction(target)
    return False


def _power_ge(n,k,q,target):
    value=float(binom.cdf(k,n,q))
    return _decimal_power_at_least(n,k,q,target) if abs(value-target)<=1e-12*target else value>=target

def tail_target(n,percentile=.99,target_power=.8,alpha=.025,*,future_ambiguous_budget=0):
    n=integer(n,'n',1); epsilon=percentile_tail(percentile)
    pi=probability(target_power,'target_power'); alpha=error_probability(alpha)
    b=integer(future_ambiguous_budget,'future_ambiguous_budget')
    return _target(n,epsilon,pi,alpha,b)

@lru_cache(maxsize=4096)
def _target(n,epsilon,pi,alpha,b):
    k=binomial_cutoff(n,epsilon,alpha)
    if k is None or k<b: return None
    k-=b
    lo,hi=0.,1.
    for _ in range(90):
        mid=(lo+hi)/2
        if mid==lo or mid==hi:break
        if binom.cdf(k,n,mid)>=pi: lo=mid
        else: hi=mid
    # A rounded CDF equal to pi can conceal a slightly unsafe inverse root.
    for _ in range(128):
        if _decimal_power_at_least(n,k,lo,pi):return lo
        lo=math.nextafter(lo,0.)
    lo*=1-1e-8
    if not _decimal_power_at_least(n,k,lo,pi):return 0.
    return lo

def minimum_reference(q,delta=.05):
    q=real(q,'q',minimum=0,maximum=1,open_max=True); delta=error_probability(delta,'delta')
    if q==0: return None
    value=math.log(delta)/math.log1p(-q)
    if value>2**53-1: return None
    n=max(1,math.ceil(value))
    # Logs are only a seed: their rounded equality can be off by one.
    # Never report an exact minimum without proving both sides of the boundary.
    for _ in range(128):
        relation=_cdf_relation(n,0,q,delta)
        if relation is None:return None
        if relation==1:
            if n==2**53-1:return None
            n+=1
            continue
        if n==1:return n
        previous=_cdf_relation(n-1,0,q,delta)
        if previous is None:return None
        if previous==1:return n
        n-=1
    return None

def _count_pair(k,n,alpha):
    n=integer(n,'n',1); k=integer(k,'k'); alpha=error_probability(alpha)
    if k>n: raise ValidationError('k cannot exceed n')
    return k,n,alpha

def _cp_upper_seed(k,n,alpha):
    """Unverified inverse seed; never expose it as a confidence bound."""
    if k==0: return -math.expm1(math.log(alpha)/n)
    return float(beta.isf(alpha,k+1,n-k))


def cp_upper(k,n,alpha):
    """Conservative one-sided CP bound verified by directed CDF inversion.

    A rounded beta quantile, even with one added ULP, is not a certificate.
    Outside the directed regime or after an unresolved inversion, return the
    valid uninformative upper bound 1 rather than an optimistic estimate.
    """
    k,n,alpha=_count_pair(k,n,alpha)
    return _cp_upper(k,n,alpha)


@lru_cache(maxsize=8192)
def _cp_upper(k,n,alpha):
    if k==n or k>100000:return 1.
    candidate=_cp_upper_seed(k,n,alpha)
    if not math.isfinite(candidate):return 1.
    for _ in range(128):
        if candidate>=1.:return 1.
        if _cdf_relation(n,k,candidate,alpha)==-1:return candidate
        candidate=math.nextafter(candidate,1.)
    return 1.


def cp_lower(k,n,alpha):
    """Conservative lower CP bound by complemented upper-bound symmetry."""
    k,n,alpha=_count_pair(k,n,alpha)
    if k==0: return 0.
    return _rounded(1-Fraction(cp_upper(n-k,n,alpha)))
