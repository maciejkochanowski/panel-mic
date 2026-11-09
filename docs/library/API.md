# Python API and methodological scope

Import public functions from panel_mic. They return dictionaries except the reference classes, parser records and QC power scalar described below. Invalid input raises ValidationError, a ValueError subclass. Public exports are Reference, Interval, ValidationError, assess_upper, design_upper, frontier_upper, diagnose_upper, protocol_assurance, protocol_assurance_grid, qc_power_lower, qc_profile_from_counts and parse_mic. The package also exposes __version__.

## Observation model

**Reference(thresholds, counts, *, iid=False, label="", unit="mg/L", metadata=None)** represents counts in (0,t₁], (t₁,t₂], …, (tₘ,∞). Thresholds increase strictly; m thresholds require m+1 exact integer counts, including both tails. The total denominator must be positive. IID is a declaration, not a diagnostic result.

**Reference.from_intervals(rows, *, iid=False, label="", unit="mg/L", metadata=None)** accepts dictionaries with lower, upper and count, optionally lower_closed, or explicit point and count. A null lower endpoint means 0; a null upper endpoint means infinity. Intervals are open at the lower endpoint by default and closed at a finite upper endpoint. Closed equal endpoints represent a point; outcome support remains positive. Unknown (0,∞) observations stay in the denominator. No exact latent MIC is imputed.

**Interval(lower, upper, count, lower_closed=False)** is the immutable record stored in a Reference. Direct record construction is a low-level operation; use Reference.from_intervals to validate user inputs.

**reference.tail_counts(threshold)** returns (definite exceedances, possible exceedances) for the strict event X > threshold. Reference attributes include intervals, thresholds, n, iid, label, unit and copied metadata.

**parse_mic(value, unit, *, numeric_is_complete=False, grid=None, count=1)** converts supported units to mg/L and returns a row plus original_value, original_unit, unit, numeric_is_complete_declared and QC_reason. Unknown units, missing values, unsupported syntax and unresolved combination components retain unknown intervals with reasons. A physical grid is expressed in mg/L. Numeric complete points require explicit source justification; otherwise numeric readings remain intervals. Left censoring conservatively includes the finite upper endpoint. Preserve parser records alongside Reference inputs:

    from panel_mic import parse_mic, Reference
    records = [parse_mic(v, "mg/L", grid=[.5, 1, 2])
               for v in ["1", ">2", None]]
    reference = Reference.from_intervals(
        [r["row"] for r in records], iid=False,
        metadata={"original_n": 3, "parse_records": records},
    )
    assert reference.n == 3

## Current assessment

**assess_upper(reference, upper_limit, *, percentile=.99, alpha=.025)** evaluates a fixed upper threshold using definite and possible exceedances. Percentile is .99 or .999. Output includes status, tail_bounds, empirical_tail_interval and a quantile location interval with confidence 1−2α. A null upper endpoint is unbounded. PANEL_REACHES_QUANTILE, PANEL_TOO_NARROW and UNDECIDED distinguish positive evidence, negative evidence and insufficient information.

## Future planning

**design_upper(reference, candidates, *, future_n, percentile=.99, target_power=.8, alpha=.025, delta=.05, drift=0, unmodelled_share=0., future_ambiguous_budget=0)** selects the first permitted finite threshold meeting a fixed scenario. Drift is a nonnegative number of log₂ steps; the reference threshold is c/2ᴰ. Unmodelled population share η lies in [0,1). A hard ambiguity budget bounds unresolved outcomes in the future batch; it is not an estimated average.

Output includes q_star, q_eff, reference_cutoff, minimum_reference_zero_tail, candidate_diagnostics and certified_power_floor. The reference procedure controls the unconditional probability of selecting an inadequately powered design. It does not give a 1−δ posterior probability of adequacy conditional on issuing a recommendation. At a single fixed transformed target candidate tests are nested, with their own selection protection and no blanket candidate-count Bonferroni penalty. Changing the target or population-change scenario after examining reference data requires a separate argument or the frontier interface.

**frontier_upper(reference, candidates, *, future_ns, drifts=(0,), unmodelled_shares=(0.,), percentile=.99, target_power=.8, alpha=.025, delta=.05, grid_prespecified=False, future_ambiguous_budget=0)** constructs simultaneous tail bounds over unique transformed thresholds. A prespecified grid declaration is required. Every scenario is retained with its status, candidate bounds, selected threshold and power floor. No single monotone minimum N is assumed.

**diagnose_upper(reference, upper_limit, *, future_n=None, percentile=.99, target_power=.8, alpha=.025, delta=.05, drift=0, unmodelled_share=0., max_future_n=1000000, future_ambiguous_budget=0)** describes information obstructions at a fixed threshold. It returns conditional_analysis, complete future_n_intervals through the cap, and optional current_scenario. The target must exceed alpha; max_future_n cannot exceed 10,000,000. An optional future_n is evaluated separately from that scan cap. Empty intervals exclude only the scanned range. The IID gate controls the inferential guarantee. This is fixed-reference analysis, not optional stopping while collecting more reference data.

The article's practical primary scope is future N ≤1000. Larger computational limits do not extend the study's evidence.

## QC profiles

**qc_power_lower(future_n, q, ambiguity_tail_caps, *, percentile=.99, alpha=.025)** returns a conservative power floor under uniform bounds Pr(R≥r)≤γᵣ supplied as (r,γᵣ) pairs. It does not assume independence between unresolved and latent exceedance counts. Caps must describe the same threshold, future batch size, observation process and target population. The function treats caps as given and does not know their calibration confidence or transport validity. Future N is limited to 10,000,000.

**qc_profile_from_counts(unresolved_counts, *, future_n, thresholds, error=.05, iid_batches=False, same_target_process=False, thresholds_prespecified=False)** constructs simultaneous conservative CP bounds from independent fixed-size QC batches. Unknown batches cannot be discarded or replaced with zero. Counts are exact integers in [0,future_n]; prespecified count thresholds lie in [1,future_n]. Missing declarations yield QC_ASSUMPTIONS_NOT_DECLARED and no usable caps.

A successful profile retains future_n, batches, QC_error, declaration flags, scope, uniform_tail_caps, per_threshold_error and simultaneous_confidence_at_least. Save the whole profile and external process metadata. Never pass only its caps to a different N or assay and imply automatic transport.

The executable [QC example](../../examples/library/qc_profile.py) checks future N and explicit threshold/process metadata before use. It distinguishes hypothetical caps from a learned profile. QC error remains separately reported; no combined confidence is inferred automatically. Combining reference, QC and selection uncertainties needs an explicit allocation for the joint procedure.

## Fixed-protocol assurance

**protocol_assurance(reference_n, control_n, *, future_n, possible_tail_reference, possible_tail_control, independent_stages=False, percentile=.99, target_power=.8, alpha=.025, reference_error=.05, control_error=.025, threshold_family_size=1, future_ambiguous_budget=0, ambiguity_tail_caps=None, unmodelled_share=0.)** calculates the probability of both stage certificates passing for hypothetical observation-law probabilities. It does not fit these inputs or optimize sample allocation.

A nonzero hard ambiguity budget and QC caps cannot both be supplied. Stage error allocations protect the declared threshold family. Each stage's marginal probability is retained when the other stage cannot pass. With independent stages the joint probability is their product and a directed enclosure is supplied; otherwise joint_success_probability is null and joint_success_interval gives Fréchet bounds when both stages are feasible. A structurally failing stage implies joint probability zero. Schema and software versions support provenance. QC calibration confidence remains separate.

**protocol_assurance_grid(future_ns, *, minimum_joint_success=.8, **scenario)** evaluates distinct declared N using protocol_assurance. Pass fixed stage sizes and scenario inputs in scenario, omitting future_n. Output retains rows, evaluated_ns, admissible_future_n_intervals and complete_integer_scan. Only adjacent evaluated admissible integers are merged; unevaluated gaps are never filled. There is no cost optimizer or claim that every larger N improves assurance.

## Relation to the article

| Methodological component | Reusable interface |
| --- | --- |
| Strict-event ambiguity classification | Reference.tail_counts; parse_mic |
| Current quantile certificate | assess_upper |
| Fixed-scenario sharp-power planning | design_upper |
| Simultaneous scenario exploration | frontier_upper |
| Information obstructions and finite N windows | diagnose_upper |
| Uniform future ambiguity profiles | qc_power_lower; qc_profile_from_counts |
| Two-stage hypothetical assurance | protocol_assurance; protocol_assurance_grid |
| Restricted randomized comparator | Scientific supplement only |

Private numerical helpers are implementation details. Mathematical claims reside in the article and proofs; software supplies a reproducible operational implementation.
