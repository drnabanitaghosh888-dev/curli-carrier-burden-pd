# PHASE 2B Multi-Cohort Strategy

## Cohort Roles

- Discovery cohort: `WALLEN_PRJNA834801`
- Replication candidates:
  - `PALACIOS_PRODROMAL_PD`
  - `NISHIWAKI_MULTICOUNTRY_PD` (or verified source cohorts)
  - `MAO_CENTRAL_CHINA_PD`
  - `INTEGRATED_US_MULTICOHORT_PD`

All replication cohorts remain candidate status until verified.

Scientific rule: Only shotgun metagenomic datasets are eligible for the core
curli/csg analysis. 16S datasets lack gene-level resolution and are excluded
from curli gene detection, csg operon analysis, raw csg screening, strain-level
reconstruction, and primary association testing. They may be retained only as
optional background ecological context and should not be used in the main
analytical workflow.

## Planned Final Logic

1. Discovery in Wallen processed tables.
2. Replication across independent cohorts.
3. Selected raw-read validation on prioritized subset only.
4. Optional strain-focused reconstruction after validation.

## Curli Interpretation Guardrails

From Phase 1:
- Reliable for downstream screening interpretation: `csgD, csgE, csgF, csgG`.
- Cautious interpretation: `csgA/csgB` and `csgC/csgE` due to cross-reactivity context.
