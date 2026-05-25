# PHASE 2F Integrated-US Processed-Table Replication Plan

- Cohort: `INTEGRATED_US_MULTICOHORT_PD`
- Role: independent processed-table replication cohort
- Discovery cohort: Wallen
- Replication cohort: Integrated-US
- Raw-read validation: not performed in Phase 2F

## Scientific scope

This phase performs an ecological processed-table replication of the Curli Carrier Burden signal.
It does not claim gene-level or operon-level validation.
Curli burden is interpreted strictly as a taxonomic proxy signal.

## Exposure and outcomes

- Primary exposure: `CurliCarrierBurden`
- Primary predictor: `CurliCarrierBurden_log1p = log1p(CurliCarrierBurden)`
- Secondary exposure: `CurliCarrierPresence`
- Outcome: `PD_binary`

## Comparisons

- Primary: PD vs all controls (`PC + HC`)
- Sensitivity: PD vs PC only
- Sensitivity: PD vs HC only

## Covariates

- Primary adjusted covariates: `host_age`, `sex`, `host_body_mass_index`
- `donor_group` is not used as an adjustment covariate in the all-control model.
- `donor_group` is used only for sensitivity subsetting.

## Curli taxon matching policy

- Primary burden must use exact/binomial species-level matches.
- Genus fallback matches are excluded from the primary burden.
- Genus fallback may be tracked as exploratory transparency only.
