# Changes

## 2.0.0rc1

Candidate accompanying the methodological article; publication and final verification are separate release gates.

- Centralize conservative Clopper–Pearson bounds and use directed numerical checks at decision boundaries.
- Verify zero-event minimum reference sample sizes against both adjacent integers.
- Preserve each protocol stage's marginal success probability when the other stage has insufficient information.
- Add software and schema provenance to protocol assurance results.
- Make CLI explanations English by default, with --language pl retaining Polish; JSON status codes and field names remain language-independent.
- Supply operation-specific input schemas and document runtime constraints beyond JSON Schema.
- Document every exported Python API, conditional QC use, and a current-certification versus future-power tutorial.
- Separate the immutable simulation-generation campaign from corrected analysis provenance in the full scientific package.

Existing Python operation signatures and CLI operation names are retained. New provenance fields are additive. The Python presentation helper panel_mic.cli.explain(result) retains its historical Polish default; pass language="en" for English.

Historical generation records in the scientific supplement describe the version that produced them. They are not test receipts for this candidate.
