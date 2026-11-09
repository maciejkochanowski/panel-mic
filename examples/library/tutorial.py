"""Synthetic count illustration; run with an installed panel-mic release."""
import json
from panel_mic import (
    Reference, __version__, assess_upper, design_upper, diagnose_upper,
    protocol_assurance,
)


def main():
    reference=Reference([1.0],[392,0],iid=True,unit="mg/L",
                        metadata={"original_n":392,"source":"synthetic count illustration"})
    current=assess_upper(reference,1.0,percentile=.99,alpha=.025)
    future=design_upper(reference,[1.0],future_n=874,percentile=.99,
                        target_power=.8,alpha=.025,delta=.05)
    diagnostic=diagnose_upper(reference,1.0,future_n=874,max_future_n=1000,
                              percentile=.99,target_power=.8,alpha=.025,delta=.05)
    assurance=protocol_assurance(1139,1402,future_n=874,
        possible_tail_reference=.001,possible_tail_control=.001,
        independent_stages=True)
    assert current["status"]=="PANEL_REACHES_QUANTILE"
    assert future["status"]=="REFERENCE_N_TOO_SMALL"
    assert future["minimum_reference_zero_tail"]==1139
    assert .00936 < current["tail_bounds"][1] < .00937
    assert .07868 < assurance["joint_success_probability"] < .07870
    print(json.dumps({"software_version":__version__,"input_scope":"synthetic illustration",
                      "current":current,"future":future,"diagnostic":diagnostic,
                      "protocol":assurance},indent=2,allow_nan=False))


if __name__=="__main__":
    main()
