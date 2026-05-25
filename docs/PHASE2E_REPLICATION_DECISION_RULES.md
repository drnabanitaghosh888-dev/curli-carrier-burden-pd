# PHASE2E Replication Decision Rules

## Ready Rules
A replication cohort is ready for processed-table replication analysis only if:
- sequencing type is shotgun metagenomics
- species abundance table is available
- sample metadata is available and linkable by sample ID
- PD/control labels are available
- access/reuse terms are acceptable

## Not Ready Rules
A replication cohort is not ready if:
- only 16S data are available
- processed tables are missing
- metadata cannot be linked to samples
- PD/control labels are missing
- access is restricted and unavailable

## Handling Missing Gene-Family Tables
- Species-level replication can proceed as the minimum viable replication tier.
- Missing gene-family/pathway tables should be recorded as limitations, not silent failures.

## Handling Species-Only Replication
- Evaluate directional consistency and effect-size uncertainty.
- Do not overstate significance when covariate-adjusted effects are unstable.

## Interpretation Logic
- Directionally consistent but non-significant replication can still support a weak hypothesis.
- Null replication weakens curli-specific claims and should be reported explicitly.
- Any processed-table replication remains proxy-level and requires later raw-read validation.
