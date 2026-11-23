# Cell-Specific Mixing Score (CMS)

Tests the hypothesis that group-specific distance distributions of
k-nearest neighbor cells have the same underlying unspecified
distribution. Computes cell-specific p-values via Anderson-Darling tests
(if kSamples installed) or two-sample Kolmogorov-Smirnov /
Kruskal-Wallis tests, with optional delegation to CellMixS::cms.

## Usage

``` r
calc_cms(
  data = NULL,
  batch_info,
  k = NULL,
  n_pcs = 30,
  coords = NULL,
  cell_min = 4,
  unbalanced = FALSE,
  threads = 1,
  knn_res = NULL
)
```

## Arguments

- data:

  Count or log-normalized expression matrix (features x cells) or NULL
  if coords provided.

- batch_info:

  Factor or character vector of batch assignments for each cell.

- k:

  Neighborhood size k for k-NN graph. Default is
  min(table(batch_info))/2 or 30.

- n_pcs:

  Number of principal components to compute if coords is NULL. Default
  is 30.

- coords:

  Optional precomputed embedding matrix (cells x dimensions, e.g. PCA or
  UMAP).

- cell_min:

  Minimum number of cells per batch in neighborhood to perform test.
  Default is 4.

- unbalanced:

  Logical, whether to set pure neighborhoods to NA. Default is FALSE.

- threads:

  Number of CPU threads (default 1).

- knn_res:

  Optional precomputed k-NN list (\`idx\`, \`dists\`) from
  \`.knn_coords\`.

## Value

A named list:

- cms_scores:

  Numeric vector of CMS p-values for each cell (higher means better
  mixed)

- mean_cms:

  Mean CMS score across cells

- median_cms:

  Median CMS score across cells

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_cms(coords = coords, batch_info = batch)
#> $cms_scores
#>  [1] 0.1852750 0.3718450 0.9227450 0.0000000 0.3309650 0.4432350 0.4957950
#>  [8] 0.2241250 0.4432350 0.0955915 0.1502900 0.0861890 0.8772100 0.6694650
#> [15] 0.1240900 0.0000000 0.2751250 0.2155450 0.4957950 0.0459320 0.0233315
#> [22] 0.0784810 0.0459320 0.0000000 0.0233315 0.9160050 0.9885250 0.6077850
#> [29] 0.0967115 0.0000000 0.9777000 0.1027295 0.0000000 0.0000000 0.6453150
#> [36] 0.2721950 0.4631250 0.2174450 0.5535150 0.0000000 0.0967115 0.1654300
#> [43] 0.0000000 0.0000000 0.6694650 0.1240900 0.6071950 0.5308950 0.9907300
#> [50] 0.6079500
#> 
#> $mean_cms
#> [1] 0.3251409
#> 
#> $median_cms
#> [1] 0.216495
#> 
```
