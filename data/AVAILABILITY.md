# Data availability

All raw data used in this study are public. Small reference inputs are included in this
repository (`data/raw/`, `data/processed/`). The following files are too large for the
repository; download them from the listed public sources and verify against the SHA-256
checksums below (computed 2026-11-17). Place them under `data/raw/` using the exact
file names.

## GDSC screening data (Sanger)

Source: `ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/`
(release 24 July 2022; the HTTPS mirror is the same directory).

| file | size | SHA-256 |
|---|---|---|
| `GDSC2_fitted_dose_response_24Jul22.csv` | 38 MB | 9d6dc82b9a072bb00b6139dcd6c4c11de1e047f210382fb726e92ab8dd0eae2a |
| `GDSC1_fitted_dose_response_24Jul22.csv` | 52 MB | 7101736651d0f0e3073e55f073f0eb970480706444c96635b7404cb993ee4987 |
| `GDSC2_public_raw_data_24Jul22.csv` | 2.0 GB | e915be2948b174982a9b64bea5fab00f6ec8c15f7baa76bc7e2f7a7df295500e |

The fitted tables are consumed by scripts 13–24; the 2 GB per-well raw file is needed
only by scripts 04 (extraction to `data/processed/coread_raw_wells.parquet`, which IS
included) and 31 (growth-rate-corrected rederivation).

## Cell Model Passports (Sanger)

Source: https://cellmodelpassports.sanger.ac.uk/downloads (bulk download).
`mutations_latest.csv.gz`, `model_list_latest.csv.gz`, `Cell_Lines_Details.xlsx` and
`growth_rate_20220907.csv` are small and are included in the repository directly.

## CCLE baseline transcriptomes

| file | size | SHA-256 |
|---|---|---|
| `CCLE_RNAseq_genes_rpkm_20180929.gct.gz` | 137 MB | 954c94233cac97a695567497dcd103bd01d342b507385dc7949777267cf2440a |

Source: Broad Institute CCLE portal (RPKM, release 2018-09-29), now redistributed through
the DepMap portal (https://depmap.org/portal/download/all/). Consumed by scripts 20–21.

## Third-party replication PharmacoSets (Orcestra)

Both are canonical PharmacoSets distributed by Orcestra (https://www.orcestra.ca/psets;
API `https://www.orcestra.ca/api/pset/canonical`). Direct download links are of the form
`https://storage.googleapis.com/orcestra-published-assets/<PSetName>/datasets/<file>.rds`.

| file | size | SHA-256 | Zenodo DOI |
|---|---|---|---|
| `PSet_CTRPv2_2015.rds` | 39 MB | 95c74304c08912f69d5e35e2863c1fe248cd0cbd21ade282a5e4e45200ee442e | 10.5281/zenodo.7826870 |
| `PSet_PRISM_2020.rds` | 82 MB | 023b2d20c1a029d18cd92558e09fae599a5952655383d3a3456ed15ba6c603a8 | 10.5281/zenodo.7826864 |

Consumed by scripts 32/32b. Parsing uses the Python package `rdata` 1.1.0
(`read_rds()` performs the conversion on read; do not call `.convert()`).

## cBioPortal cohorts

Clinical and genomic data are fetched live through the public cBioPortal REST API
(https://www.cbioportal.org/api): studies `coad_silu_2022` (AC-ICAM), `coadread_tcga`
(TCGA COADREAD) and `ccle_broad_2019` (BRAF V600 calls for the third-party replication).
The BRAF mutation table and sample lists fetched on 2026-11-16/17 are included in
`data/raw/` as JSON snapshots for provenance. No credentials are required.

## Licensing of source data

GDSC and Cell Model Passports data are released by the Wellcome Sanger Institute under
CC-BY 4.0 terms for academic use; CCLE/DepMap data are distributed under the DepMap
terms of use; CTRP and PRISM PharmacoSets are archived on Zenodo under CC-BY 4.0. Users
must comply with the respective source licences; this repository licenses only its own
code and derived tables (MIT, see LICENSE).
