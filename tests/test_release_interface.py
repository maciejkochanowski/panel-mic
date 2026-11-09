"""Release contracts not covered by the historical scientific tests."""
from copy import deepcopy
from importlib.resources import files
import inspect
import json
from pathlib import Path
import runpy

import jsonschema
import pytest

from panel_mic import Reference, ValidationError, __version__, design_upper
from panel_mic.cli import OPERATIONS, explain, main, run_operation

ROOT=Path(__file__).resolve().parents[1]
EXAMPLES={
    "assess-upper":"current-certification.json",
    "design-upper":"future-design.json",
    "diagnose-upper":"future-diagnostic.json",
    "protocol-assurance":"protocol-assurance.json",
}


def schema_for(operation):
    return json.loads(files("panel_mic").joinpath(
        "schemas/"+operation+".json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("operation",list(OPERATIONS))
def test_operation_schema_matches_signature_and_required_fields(operation):
    schema=schema_for(operation)
    signature=inspect.signature(OPERATIONS[operation])
    required=[name for name,p in signature.parameters.items()
              if p.default is inspect.Parameter.empty]
    assert set(schema["required"])==set(required)
    assert set(schema["properties"])==set(signature.parameters)
    jsonschema.Draft202012Validator.check_schema(schema)
    if operation=="frontier-upper":
        payload=json.loads((ROOT/"examples/frontier.json").read_text())
    else:
        payload=json.loads((ROOT/"examples/library"/EXAMPLES[operation]).read_text())
    jsonschema.validate(payload,schema)
    for key in required:
        incomplete=deepcopy(payload)
        incomplete.pop(key)
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(incomplete,schema)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate({**payload,"typo":1},schema)


def test_cli_text_default_and_polish_and_language_independent_json(capsys):
    path=ROOT/"examples/library/future-design.json"
    assert main(["design-upper",str(path),"--text"])==0
    english=capsys.readouterr().out
    assert "Reference observations: 392" in english
    assert "REFERENCE_N_TOO_SMALL" in english
    assert main(["design-upper",str(path),"--text","--language","pl"])==0
    polish=capsys.readouterr().out
    assert "Liczba obserwacji referencyjnych: 392" in polish
    outputs=[]
    for language in ("en","pl"):
        assert main(["design-upper",str(path),"--language",language])==0
        outputs.append(json.loads(capsys.readouterr().out))
    assert outputs[0]==outputs[1]
    assert outputs[0]["software_version"]==__version__


def test_protocol_provenance_is_in_machine_result():
    payload=json.loads((ROOT/"examples/library/protocol-assurance.json").read_text())
    result=run_operation("protocol-assurance",payload)
    assert result["software_version"]==__version__
    assert result["schema_version"].startswith("panel_mic.")
    assert "Stage success probabilities:" in explain(result,language="en")


def test_structurally_integral_float_remains_runtime_error():
    payload=json.loads((ROOT/"examples/library/future-design.json").read_text())
    payload["future_n"]=874.0
    jsonschema.validate(payload,schema_for("design-upper"))
    with pytest.raises(ValidationError,match="exact integer"):
        run_operation("design-upper",payload)


def test_legacy_presentation_default_and_unknown_language():
    result=design_upper(Reference([1],[392,0],iid=True),[1],future_n=874)
    assert "Wynik:" in explain(result)
    with pytest.raises(ValidationError):
        explain(result,language="de")


@pytest.mark.parametrize("filename",["tutorial.py","qc_profile.py"])
def test_english_python_examples_execute(filename,capsys):
    runpy.run_path(str(ROOT/"examples/library"/filename),run_name="__main__")
    result=json.loads(capsys.readouterr().out)
    assert result["software_version"]==__version__


def test_qc_example_rejects_mismatched_scope():
    namespace=runpy.run_path(str(ROOT/"examples/library/qc_profile.py"))
    checked=namespace["checked_caps"]
    target={"future_n":874,"threshold_mg_l":1.0,"process_id":"assay-v1"}
    profile={"future_n":874,"status":"QC_PROFILE_CONDITIONAL_ON_DECLARATIONS",
             "uniform_tail_caps":[[1,.01]]}
    record={"target":target,"profile":profile}
    assert checked(record,**target)==[[1,.01]]
    for changed in ({"future_n":875},{"threshold_mg_l":2.0},{"process_id":"assay-v2"}):
        with pytest.raises(ValueError):
            checked(record,**{**target,**changed})
