# Making Hospital Context Visible to Support Learning Health Systems:

This repository is the public code-and-metadata companion for the paper:

**Making Hospital Context Visible to Support Learning Health Systems for Acute Care: A Transparent Observed-Data Workflow Using Linked National Survey Data**

## Archiving and citation

This repository is prepared for Zenodo archiving. Public releases contain code, metadata, mock templates, and non-disclosive manuscript outputs only. Proprietary AHA source data and restricted row-level derived datasets are not redistributed.

## What this repository contains

This release is designed for **GitHub + Zenodo DOI archiving** in a setting where the underlying source data are licensed and cannot be redistributed.

Included here:

- `scripts/` — Scripts 13–17 used for the paper
- `metadata/` — variable manifest, file-layout crosswalk, and selected-variable templates
- `mock/` — empty templates showing expected input structure
- `outputs/tables/` — non-disclosive aggregate outputs used in the manuscript
- `outputs/figures/` — manuscript figures
- `docs/` — data availability and reproducibility notes
- `environment.yml` — a portable environment specification
- `CITATION.cff` and `.zenodo.json` — repository citation and Zenodo metadata
- `LICENSE` — license for the code/documentation in this repository

## What this repository does **not** contain

The underlying AHA Annual Survey Database (ASDB) and AHA Information Technology data products are proprietary/licensed and are **not included** here.

Also excluded:

- raw AHA data files
- row-level merged analysis files
- hospital-level profile assignment outputs
- system-identifiable outputs
- local path configuration files used in the author’s computing environment

## Workflow order

The paper used the following script order:

1. `13_rebuild_weight_and_mask_missingness_v4.R`
2. `14_rebuild_masked_smc_network_v4.R`
3. `15_rebuild_weighted_masked_clustering_v5.R`
4. `16_rebuild_system_nesting_and_validation_masked_v7.R`
5. `17_update_manuscript_tables_and_figures_v3.R`

## Core design choices

- **Observed-data workflow**: no imputation of unreported capabilities
- **Inverse-probability weighting**: response weights estimated from consistently available structural variables
- **Masking**: pairwise similarity and score construction only from jointly observed values
- **Public release constraint**: only non-disclosive aggregate outputs are redistributed

## Weighting and masking summary

The weighting script defines an IT response indicator based on whether at least one selected IT variable is observed, then fits a parsimonious inverse-probability weighting model using consistently available structural variables such as state/region, bed-size-like variables, teaching-related variables, ownership/control variables, rural/urban indicators, and system/network-related variables. The resulting weights are trimmed at lower and upper quantiles before downstream use.

Masking means the workflow does not fill in missing capability values. Similarities and module scores are computed only from observed overlap. This preserves reported capabilities while avoiding synthetic values, at the cost of leaving some hospitals unassigned under the current completeness rule.

## How to use this repository

This package is suitable for:

- auditing the published workflow
- reviewing the variable curation logic
- reproducing manuscript figures and aggregate outputs
- adapting the code locally if you independently obtain authorized AHA access

It is **not** a turnkey rerun package unless you already have the licensed AHA data and can map them to the expected local layout.

## Recommended citation

Please cite both:
1. the associated manuscript, and
2. the Zenodo DOI for this repository release.

See `CITATION.cff` for machine-readable citation metadata.
