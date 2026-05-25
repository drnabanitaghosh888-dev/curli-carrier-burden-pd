# Phase 2D Statistical Model Plan

## Primary Comparison
- Compare `CurliCarrierBurden` between PD and Control using Mann-Whitney U.

## Regression Models
- Unadjusted logistic model:
  - `PD_binary ~ CurliCarrierBurden_log1p`
- Adjusted logistic model (if covariates available):
  - `PD_binary ~ CurliCarrierBurden_log1p + Age_at_collection + Sex + BMI`

## Secondary Presence/Absence Models
- Unadjusted logistic model:
  - `PD_binary ~ CurliCarrierPresence`
- Adjusted logistic model:
  - `PD_binary ~ CurliCarrierPresence + Age_at_collection + Sex + BMI`

## Exploratory Body-First Proxy Model
- Run only if constipation/body-first proxy variable is available and interpretable:
  - `constipation_binary ~ CurliCarrierBurden_log1p + Age_at_collection + Sex + BMI`

## Multiple Testing
- Apply Benjamini-Hochberg correction for model families and taxa-level exploratory tests.

## Reporting Rules
- Report effect direction, uncertainty, and corrected p-values.
- Use raw burden for non-parametric group comparison and log1p burden as the main regression predictor.
- Treat burden presence/absence as secondary evidence.
- Mark all outputs as processed-table proxy results.
- Do not present gene-level or operon-level conclusions from Phase 2D.
