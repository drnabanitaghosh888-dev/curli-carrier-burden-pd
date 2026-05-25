# Phase 2D Curli Carrier Burden Interpretation

## Definition
For sample `s`:

`CurliCarrierBurden_s = sum_i abundance_s,i * CurliCandidateConfidence_i`

Confidence weights:
- high = 1.0
- medium = 0.5
- low = 0.25
- unknown = 0.25

## Interpretation Boundaries
- This is a processed-table proxy signal.
- It is not proof of csg gene presence.
- It is not proof of curli expression.
- It is not proof of intact csg operon architecture.

## Biological Caution
Species-level abundance patterns provide ecological signal only. Gene-level and operon-level claims require later raw-read validation.

## Required Next Layer
Phase 3 remains blocked until run accession mapping is available. Raw-read subset screening is required before gene-level conclusions.
