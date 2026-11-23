# Run edgeR QLF Differential Expression with Cellular Detection Rate (edgeRQLFDetRate)

Implements the benchmark-optimal `edgeRQLFDetRate` algorithm (Soneson &
Robinson, 2018; Duo et al., 2024 / Simpipe) for single-cell differential
expression analysis. Fits a quasi-likelihood negative binomial
generalized linear model including the cellular detection rate (DetRate:
fraction of expressed genes per cell) as a covariate.

## Usage

``` r
run_edger_qlf_detrate(
  counts,
  group,
  cdr = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0,
  min_cells = 3
)
```

## Arguments

- counts:

  Single-cell count matrix (genes x cells), or `SingleCellExperiment`
  object.

- group:

  Factor or vector specifying the comparison cell groups (must contain
  exactly 2 levels for pairwise DE).

- cdr:

  Optional numeric vector of cellular detection rates per cell. If NULL,
  computed automatically as fraction of non-zero genes.

- fdr_cutoff:

  FDR / adjusted p-value threshold for significance (default 0.05).

- logfc_cutoff:

  Absolute log2 fold-change threshold (default 0).

- min_cells:

  Minimum number of cells expressing a gene to test (default 3).

## Value

A data.frame with columns:

- gene:

  Gene identifier

- logFC:

  log2 fold change

- logCPM:

  log2 counts per million

- F_stat:

  Quasi-likelihood F-statistic (or t-statistic in fallback)

- PValue:

  Raw p-value

- FDR:

  Benjamini-Hochberg adjusted p-value

- is_deg:

  Logical flag indicating whether the gene meets the significance
  threshold

## Examples

``` r
data(example_scrna, package = "scSimEval")
top2 <- names(sort(table(example_scrna$cell_types), decreasing = TRUE))[1:2]
idx <- which(example_scrna$cell_types %in% top2)
res <- run_edger_qlf_detrate(example_scrna$ref[, idx], example_scrna$cell_types[idx])
head(res)
#>                      gene      logFC   logCPM     F_stat       PValue
#> Os05g0257100 Os05g0257100  1.2413647 15.72143  1.5688957 0.2105043876
#> Os06g0142400 Os06g0142400  0.4415879 15.32511  0.1843723 0.6676855771
#> Os02g0527300 Os02g0527300 -0.6532696 15.47727  0.3369859 0.5616353338
#> Os07g0429700 Os07g0429700 -3.3871019 16.41371 11.4313630 0.0007351286
#> Os01g0348900 Os01g0348900 -2.1212825 15.57154  4.3541055 0.0370380008
#> Os06g0142300 Os06g0142300 -1.9048782 15.43784  2.8054942 0.0940880639
#>                     FDR is_deg
#> Os05g0257100 0.45527693  FALSE
#> Os06g0142400 0.81703630  FALSE
#> Os02g0527300 0.77958337  FALSE
#> Os07g0429700 0.01709174   TRUE
#> Os01g0348900 0.17222670  FALSE
#> Os06g0142300 0.30019335  FALSE
```
