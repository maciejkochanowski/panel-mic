# From observations to a planning decision

This tutorial shows why sufficient evidence for a current P99 statement may still be insufficient for an 80% future-power target. The inputs are a synthetic reconstruction of the article's sufficient-count example, not a redistributed patient-level cohort or independent validation.

## 1. Describe the observations

Here 392 observations lie in (0, 1] mg/L and none above 1 mg/L. The represented outcome is the prespecified assay measurement, not an unobserved noise-free MIC.

    from panel_mic import Reference, assess_upper, design_upper, diagnose_upper

    reference = Reference(
        [1.0], [392, 0], iid=True, unit="mg/L",
        metadata={"original_n": 392, "source": "synthetic count illustration"},
    )

The IID declaration must be justified by the sampling design. Repeated isolates, shared patients, selection and dependence require separate consideration. The software cannot establish these assumptions.

For raw textual MIC observations, use parse_mic and retain its QC records. The [API reference](API.md) describes interval and parser semantics. A missing measurement is an unknown interval, never an observation to drop.

## 2. Assess current certification

    current = assess_upper(reference, 1.0, percentile=.99, alpha=.025)
    assert current["status"] == "PANEL_REACHES_QUANTILE"
    print(current["tail_bounds"])

The upper tail-probability bound is approximately 0.00936627, below 0.01. Under the stated model the data therefore support the current P99 threshold statement.

    panel-mic assess-upper examples/library/current-certification.json --text

## 3. Ask the future-power question

    future = design_upper(
        reference, [1.0], future_n=874, target_power=.8,
        percentile=.99, alpha=.025, delta=.05,
    )
    assert future["status"] == "REFERENCE_N_TOO_SMALL"
    print(future["q_star"], future["minimum_reference_zero_tail"])

The maximum tail probability compatible with the future power requirement is approximately 0.00262896. At reference error δ=0.05, at least 1139 reference observations with zero possible exceedances are needed. The current 392 are insufficient. This is an information limitation, not evidence that the biological population necessarily fails the target.

Current assessment uses α=0.025, whereas planning uses a separate reference error δ=0.05. These answer different questions and should not be presented as interchangeable confidence statements. The article also reports the reference tail bound at 0.025 when comparing its two-stage protocol roles.

    panel-mic design-upper examples/library/future-design.json --text
    panel-mic diagnose-upper examples/library/future-diagnostic.json --text

The diagnostic checks future N only through its declared cap. Empty admissible intervals through N=1000 do not prove impossibility for every larger N. Feasible N need not form one monotone range.

## 4. Distinguish protocol assurance from certification power

    panel-mic protocol-assurance examples/library/protocol-assurance.json --text

With reference n=1139, control n=1402, future N=874, and hypothetical possible-exceedance probabilities 0.001 in both stages, the joint success probability is about 7.87% under stage independence. The 80% target is the required future-test power at a tail probability within the certified target. The reference procedure controls the unconditional probability of selecting an inadequately powered design; it does not give a 95% probability of adequacy conditional on an issued recommendation.

The possible-exceedance probabilities are hypothetical inputs to the observation-law calculation, not fitted estimates or biological guarantees. The two-stage protocol is one study design. Basic current certification does not require two stages.

Run the complete Python example:

    python examples/library/tutorial.py

It prints the software version and four results as JSON. Save the input files and output together. For complete simulations and data provenance use the separate scientific reproduction package.
