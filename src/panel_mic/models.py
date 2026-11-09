"""Validated measured-cohort intervals; outcomes have positive support."""
from dataclasses import dataclass
from numbers import Integral, Real
import math
import json
from .errors import ValidationError

MAX_INTEGER = 2**53 - 1

def integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or not minimum <= value <= MAX_INTEGER:
        raise ValidationError(f'{name} must be an exact integer in [{minimum}, {MAX_INTEGER}]')
    return int(value)

def real(value, name, *, minimum=None, maximum=None, open_min=False, open_max=False):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValidationError(f'{name} must be a finite real number')
    value = float(value)
    if minimum is not None and (value < minimum or (open_min and value == minimum)):
        raise ValidationError(f'{name} is below its permitted range')
    if maximum is not None and (value > maximum or (open_max and value == maximum)):
        raise ValidationError(f'{name} is above its permitted range')
    return value

def grid(values):
    try:
        values = tuple(real(x, 'concentration', minimum=0, open_min=True) for x in values)
    except TypeError as exc:
        raise ValidationError('concentrations must be a sequence') from exc
    if not values or any(a >= b for a,b in zip(values, values[1:])):
        raise ValidationError('concentrations must be nonempty, strictly increasing and distinct')
    return values

def _metadata_keys(value):
    if isinstance(value,dict):
        if any(not isinstance(k,str) for k in value):
            raise ValueError('metadata keys must be strings at every depth')
        for item in value.values(): _metadata_keys(item)
    elif isinstance(value,(list,tuple)):
        for item in value: _metadata_keys(item)

@dataclass(frozen=True)
class Interval:
    lower: float
    upper: float
    count: int
    lower_closed: bool = False

@dataclass(frozen=True, init=False)
class Reference:
    intervals: tuple
    thresholds: tuple
    n: int
    iid: bool
    label: str
    unit: str
    metadata: dict

    def __init__(self, thresholds, counts, *, iid=False, label='', unit='mg/L', metadata=None):
        thresholds = grid(thresholds)
        try:
            counts = tuple(integer(x, 'count') for x in counts)
        except TypeError as exc:
            raise ValidationError('counts must be a sequence') from exc
        if len(counts) != len(thresholds)+1:
            raise ValidationError('m thresholds require exactly m+1 counts, including both censored tails')
        edges = (0.,) + thresholds + (math.inf,)
        rows = tuple(Interval(a,b,n) for a,b,n in zip(edges,edges[1:],counts))
        self._set(rows, thresholds, iid, label, unit, metadata)

    def _set(self, rows, thresholds, iid, label, unit, metadata):
        if type(iid) is not bool:
            raise ValidationError('iid must be a Boolean declaration')
        if not isinstance(label,str) or not isinstance(unit,str) or not unit.strip():
            raise ValidationError('label and a nonempty unit must be strings')
        if metadata is not None and not isinstance(metadata,dict):
            raise ValidationError('metadata must be an object')
        n = integer(sum(x.count for x in rows), 'reference denominator', 1)
        try:
            _metadata_keys(metadata)
            meta = json.loads(json.dumps(metadata or {},allow_nan=False))
        except (ValueError,TypeError,OverflowError,RecursionError) as exc:
            raise ValidationError('metadata must contain only finite JSON-compatible values') from exc
        if 'original_n' in meta and integer(meta['original_n'],'original_n',1) != n:
            raise ValidationError('original_n conflicts with the measured-cohort denominator')
        for name,value in [('intervals',rows),('thresholds',thresholds),('n',n),('iid',iid),('label',label),('unit',unit),('metadata',meta)]:
            object.__setattr__(self,name,value)

    @classmethod
    def from_intervals(cls, rows, *, iid=False, label='', unit='mg/L', metadata=None):
        intervals=[]
        try:
            for row in rows:
                if isinstance(row,dict) and set(row)=={'point','count'}:
                    value=real(row['point'],'point',minimum=0,open_min=True)
                    intervals.append(Interval(value,value,integer(row['count'],'count'),True))
                    continue
                if not isinstance(row,dict) or not {'lower','upper','count'}<=set(row) or set(row)-{'lower','upper','count','lower_closed'}:
                    raise ValidationError('each observation requires point/count or lower/upper/count with optional lower_closed')
                closed=row.get('lower_closed',False)
                if type(closed) is not bool:raise ValidationError('lower_closed must be Boolean')
                lower = 0. if row['lower'] is None else real(row['lower'],'lower',minimum=0)
                upper = math.inf if row['upper'] is None else real(row['upper'],'upper',minimum=0,open_min=True)
                if lower>upper or (lower==upper and not closed) or (lower==0 and closed):
                    raise ValidationError('observation must be nonempty with positive outcome support')
                intervals.append(Interval(lower,upper,integer(row['count'],'count'),closed))
        except TypeError as exc:
            raise ValidationError('intervals must be a sequence') from exc
        obj=cls.__new__(cls)
        thresholds=tuple(sorted({v for r in intervals for v in (r.lower,r.upper) if 0<v<math.inf}))
        obj._set(tuple(intervals),thresholds,iid,label,unit,metadata)
        return obj

    def tail_counts(self, threshold):
        threshold=real(threshold,'threshold',minimum=0)
        return (sum(r.count for r in self.intervals if r.lower>threshold or (r.lower==threshold and not r.lower_closed)),
                sum(r.count for r in self.intervals if r.upper > threshold))
