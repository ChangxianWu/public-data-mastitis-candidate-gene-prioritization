# public-data-mastitis-candidate-gene-prioritization

Code for the manuscript:

> **A four-source framework for bovine mastitis candidate gene prioritization using exclusively
> public data: a cross-cohort integrative analysis**
> Changxian Wu, Min Wang, Fan Liu, Lingxi Ma (corresponding: lingxim26@163.com)
> College of Agriculture and Animal Husbandry Science and Technology,
> Ningxia Vocational and Technical College for Minorities, Wuzhong 751100, Ningxia, China

This repository contains **pipeline code and documentation only**. Input data, derived result
tables, figures and supplementary text are not included; where to find each of them is stated
below.

## What the pipeline does

Four evidence sources are integrated for bovine mastitis candidate genes:

1. **Cross-cohort differential expression** — five GPL2112 microarray cohorts enter a
   probe-level meta-analysis (24,128 probes; Stouffer + DerSimonian-Laird, Benjamini-Hochberg
   FDR), with two RNA-seq cohorts added at the gene level.
2. **Mammary cis-eQTL evidence** — overlap with CattleGTEx mammary significant cis-eQTLs and
   fine-mapping posterior inclusion probabilities (PIP).
3. **Network annotation** — STRING v11.5 bovine network over the symbol-resolvable subset.
4. **Granularity statement** — variant-level colocalization is infeasible with public data, so
   eQTL evidence is integrated at gene level with PIP weighting.

Priority score: `0.35 x normalized DE + 0.25 x eQTL presence + 0.30 x max PIP +
0.10 x normalized network degree`. eQTL presence is an inclusion criterion and is constant
within the integrated set; differential expression and PIP are the two rank-active components.

## Repository layout

```
code/                23 stage scripts (01-22, one sub-step) + postcheck_vm2.py
README.md            this file
LICENSE              MIT licence (code)
CITATION.cff         citation metadata
requirements.txt     Python dependencies
.gitignore
MANIFEST_SHA256.txt  checksum of every file in this repository
```

## Not included in this repository

| Item | Where to find it |
|---|---|
| GEO expression series | public GEO; accessions and download stage listed below |
| CattleGTEx mammary eQTL files | https://cattlegtex.farmgtex.org/ (see `CATTLE_GTEX_MAMMARY_DIR` below) |
| STRING v11.5 bovine network | https://stringdb-downloads.org/ (files and SHA-256 below) |
| Probe-level, gene-level and per-cohort result tables; the 286-gene tier; the Top-30 table | Supplementary Tables S1-S11 (Additional file 2 of the article) |
| Supplementary methods and results (S1-S11) | Additional file 1 of the article |
| Figures 1-7 | the article and its figure files |
| Competitor gene identifiers | the competitor's own supplementary material (Cai et al. 2024) |

The code regenerates every reported number from the public inputs; result tables are published
with the article rather than duplicated here.

## Data inputs

| Input | Access | Notes |
|---|---|---|
| GEO series | GSE24560, GSE25319, GSE24217, GSE93082, GSE15025, GSE131607, GSE75379, GSE9685, GSE34031 | downloaded by `code/01_geo_download.py` |
| CattleGTEx mammary eQTL | https://cattlegtex.farmgtex.org/ | set `CATTLE_GTEX_MAMMARY_DIR`; needs `Mammary_Gland_significant_eQTL.txt.gz` and `Mammary_Gland_finemapping_eQTL.txt.gz` (the significant file is an uncompressed TSV despite its name) |
| STRING v11.5 bovine network | https://stringdb-downloads.org/ | `9913.protein.links.v11.5.txt.gz`, SHA-256 `9dde88f89998853f9c5568c1e489cdbdee3dc1311b142d7e8a743412c872731f`; `9913.protein.info.v11.5.txt.gz`, SHA-256 `16d87247e6ff0c0e0f5ee4dc0c9dcd1caa4975626dd91e076102f6ccd02f6f66` |
| Animal QTLdb | https://www.animalgenome.org/QTLdb/ | interval layer deferred (portal under maintenance during the analysis) |
| Competitor supplementary tables | Cai et al. 2024, Genet Sel Evol 56:54 (PMID 39009986), MOESM1 | required only by `15_competitor_layer.py` / `21_competitor_reconcile.py` |

## Environment

Python 3.13 with `numpy`, `scipy`, `python-docx`, `openpyxl`: `pip install -r requirements.txt`.

## Stage order

```
01_geo_download -> 02_de_queue_select -> 03_parse_expr_matrices -> 04_de_meta -> 05_de_meta_final
-> 06_platform_annotation -> 07_gb2gene / 07b_gene_anchor -> 08_gene_meta_combined
-> 13_eqtl_integrate -> 09_string_network -> 10_priority_scoring -> 11_make_figures / 12_make_figures_ext
-> 14_loqo -> 15_competitor_layer -> 16_make_tables
-> 17_robustness / 18_enrichment -> 19_review_response / 20_exact_variance_sensitivity
-> 21_competitor_reconcile -> 22_string_full_network_check -> postcheck_vm2.py
```

Run each script from the project root; scripts read and write under `work/`. The eQTL stage
requires `CATTLE_GTEX_MAMMARY_DIR`; the network stage reads a frozen edge snapshot unless
`--refresh` is passed.

## Headline results (as reported in the article)

| Quantity | Value |
|---|---|
| Cohorts assembled / synthesized | 9 (339 samples) / 7 (257 samples) |
| Probe-level meta | 24,128 probes; 7,252 at BH q<0.05; 2,999 direction-consistent |
| Exact Welch-variance sensitivity | 2,663 at q<0.05; 66.5% of Tier A retained, direction unchanged |
| Gene-level meta | 24,597 genes; 764 tri-source; 563 at BH q<0.05 (21,442 with defined p) |
| DE x eQTL overlap | 6,203 genes; 3,395 integrated with fine-mapping |
| High-confidence tier | 286 genes (PIP >= 0.5, p < 0.05); 142 at PIP >= 0.9 |
| Priority scoring | 3,392 genes ranked; Top-30 all at gene-level q<0.05 |
| Competitor comparison | 150 published identifiers -> 93 Ensembl genes; 4 of 286 shared |

## Quality gate

`code/postcheck_vm2.py` verifies figure files, recomputes the headline numbers from the result
tables, checks reference consistency and the supplementary index, and compares the submission
package against its manifest. It expects the project directory layout described in the stage
order; it is not needed to run the analysis itself.

## Note on comments

Most scripts are documented in English. A few carry inline comments in Chinese, as written
during the original analysis; all identifiers, file names and outputs are in English.

## Licence

MIT for the code (see `LICENSE`). The article's result tables are released under CC BY 4.0.
Third-party data remain under their original terms (GEO, CattleGTEx/FarmGTEx, STRING).

## Citation

Please cite the manuscript above. Repository version 1.0 (2026-10-25).
