"""Catch loss of known assay outcomes and incorrect endpoint semantics."""
import json
import math
import os
import subprocess
import sys

import pytest
import panel_mic
from panel_mic import Reference,ValidationError

def test_incomplete_declaration_and_overflow_cannot_create_a_finite_MIC():
    with pytest.raises(ValidationError):panel_mic.parse_mic(1,'mg/L',numeric_is_complete=1)
    r=panel_mic.parse_mic('1e309','mg/L',numeric_is_complete=True)
    assert r['QC_reason']=='INVALID_MIC_VALUE' and r['row']['upper'] is None


def test_points_and_closed_lower_endpoints_have_exact_tail_indicators():
    r=Reference.from_intervals([{'point':2.,'count':3},
        {'lower':2.,'upper':4.,'count':2},
        {'lower':2.,'upper':None,'lower_closed':True,'count':1},
        {'lower':0.,'upper':None,'count':1}])
    assert r.n==7
    assert r.tail_counts(2.)==(2,4)
    assert r.tail_counts(1.5)==(6,7)
    assert r.tail_counts(4.)==(0,2)
    assert r.tail_counts(0.)==(7,7)


@pytest.mark.parametrize('row',[
    {'point':0,'count':1},{'point':True,'count':1},
    {'point':math.inf,'count':1},{'point':1,'lower':0,'count':1},
    {'lower':2,'upper':2,'lower_closed':False,'count':1},
    {'lower':0,'upper':2,'lower_closed':True,'count':1},
    {'lower':1,'upper':2,'lower_closed':1,'count':1}])
def test_invalid_endpoint_contract_is_rejected(row):
    with pytest.raises(ValidationError):Reference.from_intervals([row])


def test_complete_numeric_readout_is_not_artificially_unresolved_off_grid():
    exact=panel_mic.parse_mic('2','ug/mL',numeric_is_complete=True,grid=[1,2,4])
    conservative=panel_mic.parse_mic('2','mg/L',grid=[1,2,4])
    r=Reference.from_intervals([exact['row']]);old=Reference.from_intervals([conservative['row']])
    assert r.tail_counts(1.5)==(1,1)
    assert old.tail_counts(1.5)==(0,1)
    assert exact['original_value']=='2'


def test_censoring_is_never_changed_into_a_point():
    for op,at,expected in [('<=',2.,(0,0)),('>',2.,(1,1)),('>=',2.,(0,1))]:
        r=Reference.from_intervals([panel_mic.parse_mic(op+'2','mg/L',numeric_is_complete=True)['row']])
        assert r.tail_counts(at)==expected
    r=Reference.from_intervals([panel_mic.parse_mic('>=2','mg/L')['row']])
    assert r.tail_counts(1.5)==(1,1)


@pytest.mark.parametrize('value,unit,reason',[
    ('R','mg/L','UNPARSED_MIC'),('0.5/9.5','mg/L','COMBINATION_COMPONENT_REQUIRED'),
    ('2','mm','UNKNOWN_MIC_UNIT'),(None,'mg/L','MISSING_MIC'),
    ('8','mg/L','MEASUREMENT_OFF_DECLARED_GRID')])
def test_unknown_measurements_remain_in_denominator_with_reason(value,unit,reason):
    r=panel_mic.parse_mic(value,unit,grid=[1,2,4],count=3)
    assert r['QC_reason']==reason
    cohort=Reference.from_intervals([r['row']])
    assert cohort.n==3 and cohort.tail_counts(2)==(0,3)


def test_unit_conversion_and_point_contract_survive_real_cli(tmp_path):
    row=panel_mic.parse_mic('.002','mg/mL',numeric_is_complete=True,count=1200)['row']
    assert Reference.from_intervals([row]).tail_counts(1.5)==(1200,1200)
    path=tmp_path/'known.json'
    path.write_text(json.dumps({'reference':{'intervals':[row],'iid':True},'upper_limit':1.5}))
    proc=subprocess.run([sys.executable,'-m','panel_mic','assess-upper',str(path)],
        capture_output=True,text=True,env={**os.environ,'PYTHONUTF8':'1'})
    assert proc.returncode==0,proc.stderr
    result=json.loads(proc.stdout)
    assert result['status']=='PANEL_TOO_NARROW'
    assert result['definite_exceedances']==result['possible_exceedances']==1200


def test_missing_spreadsheet_nan_keeps_a_serializable_unknown_record():
    record=panel_mic.parse_mic(float('nan'),'mg/L',count=4)
    json.dumps(record,allow_nan=False)
    assert Reference.from_intervals([record['row']]).n==4


def test_new_observation_forms_pass_the_installed_input_schema():
    import jsonschema
    from importlib.resources import files
    schema=json.loads(files('panel_mic').joinpath('schemas/input.json').read_text())
    payload={'reference':{'intervals':[{'point':2.,'count':5},
        {'lower':2.,'upper':None,'lower_closed':True,'count':1}],'iid':True},'upper_limit':1.5}
    jsonschema.validate(payload,schema)
