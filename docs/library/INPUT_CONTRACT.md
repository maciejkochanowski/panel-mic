# Input and output contract

The CLI accepts UTF-8 JSON objects. It rejects duplicate keys, NaN, infinity, unknown operation parameters, invalid counts, and input/output path collisions. It creates output files exclusively, preserving existing outputs. JSON output includes the software version and schema family for upper, protocol and batch operations.

## Operation schemas

Schemas are installed under panel_mic/schemas.

| CLI operation | Schema | Required fields |
| --- | --- | --- |
| assess-upper | assess-upper.json | reference, upper_limit |
| design-upper | design-upper.json | reference, candidates, future_n |
| frontier-upper | frontier-upper.json | reference, candidates, future_ns |
| diagnose-upper | diagnose-upper.json | reference, upper_limit |
| protocol-assurance | protocol-assurance.json | reference_n, control_n, future_n, possible_tail_reference, possible_tail_control |

input.json is a backward-compatible structural superset for the four upper operations; it does not establish that a particular operation has all required arguments. protocol.json is a compatibility copy of protocol-assurance.json. Each operation-specific schema is self-contained.

    import json
    from importlib.resources import files
    import jsonschema

    schema = json.loads(files("panel_mic").joinpath(
        "schemas/design-upper.json").read_text(encoding="utf-8"))
    payload = json.loads(open("examples/library/future-design.json").read())
    jsonschema.validate(payload, schema)

Schema validation is optional. Runtime validation remains authoritative and does not depend on jsonschema.

## Additional runtime constraints

JSON Schema's integer type describes numeric values, so some validators accept 2000.0 or 2e3 as integers. The CLI deliberately rejects these representations for discrete inputs: use 2000, not a floating-point JSON token. Python API counts must be integral objects, not floats or Booleans. No silent rounding or float-to-integer conversion occurs. The exact-integer upper limit is 2⁵³−1; some functions apply smaller computational caps.

Thresholds and candidates must be strictly increasing and distinct. Counts must have length len(thresholds)+1, with a positive total. Interval endpoints must define nonempty positive-support observations. Metadata must contain finite JSON-compatible values and string keys at every depth; original_n, if supplied, must equal the complete denominator.

The diagnostic requires target_power > alpha. Its scan cap is at most 10,000,000. QC and assurance future N also have explicit limits described in the API. Supplying a nonzero future_ambiguous_budget with QC tail caps is invalid. Schema validation alone cannot establish these cross-field, scientific or numerical conditions.

## Batch behavior and exit codes

Batch accepts exactly {"cohorts": [...]}. Every row has a unique, nonempty, trimmed string id, an operation, and that operation's input fields. Invalid rows remain with status INVALID_INPUT and their original row index. Duplicate IDs invalidate every affected row. Aggregate status is COMPLETED_WITH_ERRORS when any row is invalid.

Exit code 0 means evaluation completed, not that certification succeeded. Exit code 2 denotes invalid input, I/O failure, or a batch containing invalid rows.

## Scientific status interpretation

PANEL_REACHES_QUANTILE concerns current certification; FOUND concerns a planning rule's declared scope. UNDECIDED, REFERENCE_N_TOO_SMALL and REFERENCE_RANGE_UNRESOLVED retain lack of evidence. They must not be recoded as biological failure. IID_NOT_DECLARED prevents an unconditional inferential claim; conditional diagnostic calculations may still be reported.

GRID_NOT_PRESPECIFIED withholds the frontier's simultaneous guarantee. Fixed-scenario nested threshold selection in design_upper has a separate validity argument at its fixed transformed target. A blanket multiplicity penalty for every candidate there would misstate the method.

The --language options affect --text only. JSON keys, numeric results and status codes remain identical. For programmatic use, consume JSON or Python dictionaries rather than parsing explanatory text.
