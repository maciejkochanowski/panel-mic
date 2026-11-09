import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]

def run(*args):
    env={**os.environ,'PYTHONPATH':str(ROOT/'src'),'PYTHONUTF8':'1'}
    return subprocess.run([sys.executable,'-m','panel_mic',*map(str,args)],capture_output=True,text=True,encoding='utf8',env=env)

def doc():
    return {'reference':{'thresholds':[1,2,4],'counts':[1290,8,5,0],'iid':True},'candidates':[1,2,4,8],'future_n':1000}

def test_cli_json_and_text(tmp_path):
    p=tmp_path/'in.json'; p.write_text(json.dumps(doc()))
    r=run('design-upper',p)
    assert r.returncode==0,r.stderr
    assert json.loads(r.stdout)['upper_limit']==4
    t=run('design-upper',p,'--text')
    assert t.returncode==0 and 'FOUND' in t.stdout and '1303' in t.stdout

def test_no_overwrite_or_input_collision(tmp_path):
    p=tmp_path/'in.json'; p.write_text(json.dumps(doc())); before=p.read_bytes()
    assert run('design-upper',p,'--output',p).returncode==2
    assert p.read_bytes()==before
    q=tmp_path/'out.json'; q.write_text('previous')
    assert run('design-upper',p,'--output',q).returncode==2
    assert q.read_text()=='previous'

def test_batch_retains_all_invalid_identified_cohorts(tmp_path):
    from copy import deepcopy
    valid={'id':'a','operation':'design-upper',**doc()}; bad=deepcopy(valid)
    bad['id']='b'; bad['reference']['counts'][0]=True
    p=tmp_path/'batch.json'; p.write_text(json.dumps({'cohorts':[valid,bad,{'id':'c','operation':'no-such-op'}]}))
    r=run('batch',p)
    assert r.returncode==2
    x=json.loads(r.stdout)
    assert x['cohorts_total']==3 and x['cohorts_invalid']==2
    assert [c['id'] for c in x['cohorts']]==['a','b','c']
    assert x['cohorts'][0]['status']=='FOUND' and x['cohorts'][1]['status']=='INVALID_INPUT'

def test_strict_json_and_unknown_fields(tmp_path):
    p=tmp_path/'bad.json'
    for text in ['{"future_n": NaN}','{"reference":{},"reference":{}}',json.dumps({**doc(),'future_n_typo':1000})]:
        p.write_text(text)
        r=run('design-upper',p)
        assert r.returncode==2 and 'Traceback' not in r.stderr

def test_json_schema_examples():
    import jsonschema
    schema=json.loads((ROOT/'src/panel_mic/schemas/input.json').read_text())
    for path in (ROOT/'examples').glob('*.json'):
        if path.name!='batch.json': jsonschema.validate(json.loads(path.read_text()),schema)

def test_inprocess_cli_and_duplicate_batch_ids(tmp_path,capsys):
    from panel_mic.cli import main,run_batch,explain,reference_from_dict
    from panel_mic import ValidationError
    import pytest
    p=tmp_path/'in.json';p.write_text(json.dumps(doc()))
    assert main(['design-upper',str(p)])==0
    assert json.loads(capsys.readouterr().out)['status']=='FOUND'
    assert main(['design-upper',str(p),'--text'])==0
    assert '1303' in capsys.readouterr().out
    output=tmp_path/'out.json'
    assert main(['design-upper',str(p),'--output',str(output)])==0
    assert json.loads(output.read_text())['upper_limit']==4
    for data in [None,{'cohorts':[]},{'cohorts':[{'id':'a'},{'id':'a'},None]}]:
        if data is None or data=={'cohorts':[]}:
            with pytest.raises(ValidationError):run_batch(data)
        else:
            r=run_batch(data);assert r['cohorts_invalid']==3
            assert 'Kolekcje: 3' in explain(r)
    for data in [None,{}, {'bogus':1},{'intervals':[],'counts':[]}]:
        with pytest.raises(ValidationError):reference_from_dict(data)

def test_direct_refusal_explanations_and_parameter_validation():
    from panel_mic import Reference,design_upper,frontier_upper,ValidationError
    from panel_mic.cli import explain
    import pytest
    r=Reference([1,2],[0,2000,0],iid=True)
    for result in [design_upper(r,[1,2],future_n=1000,unmodelled_share=.1),
                   design_upper(Reference([1],[200,0],iid=True),[1],future_n=1000),
                   design_upper(Reference([1],[2000,20],iid=True),[1,2,4],future_n=1000),
                   frontier_upper(r,[1,2],future_ns=[1000]),
                   frontier_upper(r,[1,2],future_ns=[1000],grid_prespecified=True)]:
        assert result['status'] in explain(result)
    for kw in [{'future_ns':[]},{'future_ns':None},{'future_ns':[1000],'drifts':[]},{'future_ns':[1000],'grid_prespecified':1}]:
        with pytest.raises(ValidationError):frontier_upper(r,[1,2],**kw)

def test_direct_cli_errors_and_iid_explanation(tmp_path,capsys):
    import pytest
    from panel_mic.cli import main,explain
    from panel_mic import Reference,design_upper
    with pytest.raises(SystemExit) as exc:main(['design-upper',str(tmp_path/'missing.json')])
    assert exc.value.code==2 and 'Traceback' not in capsys.readouterr().err
    result=design_upper(Reference([1],[2000,0]),[1],future_n=1000)
    assert 'nie zadeklarowano' in explain(result)
