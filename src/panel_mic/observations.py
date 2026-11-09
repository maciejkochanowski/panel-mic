"""Explicit measured-MIC parsing; unqualified rows remain in their denominator."""
import math
import re

from .models import integer,grid as concentration_grid
from .errors import ValidationError


def parse_mic(value,unit,*,numeric_is_complete=False,grid=None,count=1):
    """Normalize to mg/L while preserving raw text, operators and QC reasons.

    numeric_is_complete requires source evidence that a numeric readout is the
    fully specified standard assay outcome. Otherwise a containing interval is
    retained. This never recovers a noise-free biological MIC or a missing plate.
    """
    count=integer(count,'count')
    if type(numeric_is_complete) is not bool:
        raise ValidationError('numeric_is_complete must be Boolean')
    physical=None if grid is None else concentration_grid(grid)
    unknown={'lower':0.,'upper':None,'count':count}
    result={'original_value':None if value is None else str(value),
        'original_unit':None if unit is None else str(unit),'unit':'mg/L',
        'numeric_is_complete_declared':numeric_is_complete,'row':unknown,'QC_reason':None}
    units={'mg/l':1.,'ug/ml':1.,'mcg/ml':1.,'mg/ml':1000.,'g/l':1000.,
           'ug/l':.001,'ng/ml':.001,'mg/dl':10.}
    normalized=str(unit or '').strip().lower().replace('µ','u').replace('μ','u').replace(' ','')
    if normalized not in units:return {**result,'QC_reason':'UNKNOWN_MIC_UNIT'}
    text=str(value if value is not None else '').strip().replace('≤','<=').replace('≥','>=')
    if not text or text.lower() in {'na','nan','none','null','nd'}:
        return {**result,'QC_reason':'MISSING_MIC'}
    if '/' in text:return {**result,'QC_reason':'COMBINATION_COMPONENT_REQUIRED'}
    match=re.fullmatch(r'(<=|>=|<|>|=)?\s*(\d+(?:\.\d*)?|\.\d+)(?:([eE][+-]?\d+))?',text)
    if match is None:return {**result,'QC_reason':'UNPARSED_MIC'}
    op=match[1] or '=';v=float(match[2]+(match[3] or ''))*units[normalized]
    if not math.isfinite(v) or v<=0:return {**result,'QC_reason':'INVALID_MIC_VALUE'}
    if physical is not None and v not in physical:
        return {**result,'QC_reason':'MEASUREMENT_OFF_DECLARED_GRID'}
    if op=='=' and numeric_is_complete:return {**result,'row':{'point':v,'count':count}}
    if op in ('>','>='):
        row={'lower':v,'upper':None,'count':count}
        if op=='>=':row['lower_closed']=True
        return {**result,'row':row}
    previous=0.
    if op=='=' and physical is not None:
        index=physical.index(v);previous=physical[index-1] if index else 0.
    row={'lower':previous,'upper':v,'count':count}
    reason='NUMERIC_OUTCOME_RETAINED_AS_INTERVAL' if op=='=' else None
    return {**result,'row':row,'QC_reason':reason}
