"""Scope-aware QC tutorial with synthetic batches, not biological validation."""
import json
from panel_mic import __version__, qc_profile_from_counts, qc_power_lower


def checked_caps(record, *, future_n, threshold_mg_l, process_id):
    """Example consumer check, not a new production profile type."""
    profile=record["profile"]
    expected={"future_n":future_n,"threshold_mg_l":threshold_mg_l,"process_id":process_id}
    if record["target"] != expected:
        raise ValueError("QC profile does not match the future N, threshold and process")
    if profile["future_n"] != future_n:
        raise ValueError("QC calibration batch size differs from future N")
    if profile["status"] != "QC_PROFILE_CONDITIONAL_ON_DECLARATIONS":
        raise ValueError("QC assumptions were not declared")
    return profile["uniform_tail_caps"]


def main():
    n=874
    profile=qc_profile_from_counts(
        [0]*200+[1]*3+[2], future_n=n, thresholds=[1,2,3],
        error=.05, iid_batches=True, same_target_process=True,
        thresholds_prespecified=True,
    )
    record={"target":{"future_n":n,"threshold_mg_l":1.0,"process_id":"synthetic-assay-v1"},
            "source":"constructed teaching batches","profile":profile}
    caps=checked_caps(record,future_n=n,threshold_mg_l=1.0,process_id="synthetic-assay-v1")
    learned=qc_power_lower(n,.001,caps)
    hypothetical=qc_power_lower(n,.001,[(1,.01)])
    print(json.dumps({
        "software_version":__version__,"calibration_record":record,
        "power_floor_conditional_on_learned_caps":learned,
        "hypothetical_given_cap_power_floor":hypothetical,
        "QC_error_separate_from_reference_error":profile["QC_error"],
        "combined_confidence_claim":None,
        "interpretation":"Conditional teaching analysis; no automatic confidence combination or QC transport."
    },indent=2,allow_nan=False))


if __name__=="__main__":
    main()
