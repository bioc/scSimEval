# Inverse Simpson Index (ISI) for Batch Mixing

Evaluates the diversity of batch labels within each cell's k-nearest
neighborhood using the Inverse Simpson Index (1 / sum(p_b^2)).
Implements a pure-R, lightweight formulation of LISI with optional
distance weighting.

## Usage

``` r
calc_isi(
  coords,
  batch_info,
  k = 30,
  weighted = TRUE,
  threads = 1,
  knn_res = NULL
)
```

## Arguments

- coords:

  Matrix of cell coordinates (cells x dimensions).

- batch_info:

  Factor or character vector of batch assignments.

- k:

  Number of nearest neighbors. Default is 30.

- weighted:

  Logical, whether to weight neighbor contributions by inverse distance.
  Default is TRUE.

- threads:

  Number of CPU threads (default 1).

- knn_res:

  Optional precomputed k-NN list (\`idx\`, \`dists\`) from
  \`.knn_coords\`.

## Value

A named list:

- isi_scores:

  Numeric vector of cell-level ISI values

- mean_isi:

  Mean ISI score (ranges from 1 to number of batches; higher indicates
  better mixing)

- median_isi:

  Median ISI score

## Examples

``` r
coords <- matrix(stats::rnorm(100), 50, 2)
batch <- factor(rep(c("B1", "B2"), length.out = 50))
calc_isi(coords, batch)
#> $isi_scores
#>  [1] 1.941166 1.999564 1.961253 1.958295 1.966258 1.998998 1.999933 1.980749
#>  [9] 1.989220 1.999493 1.982850 1.913597 1.990588 1.986243 1.951829 1.991336
#> [17] 2.000000 1.998703 1.999719 1.919800 1.999989 1.994835 1.997909 1.999591
#> [25] 1.996157 1.999995 1.999955 1.991767 1.935108 1.998598 1.952966 1.998121
#> [33] 1.915914 1.987390 1.930057 1.987932 1.928857 1.992956 1.981186 1.997833
#> [41] 1.997863 1.996592 1.966916 1.968121 1.997238 1.999933 1.995261 1.983099
#> [49] 1.984602 1.983342
#> 
#> $mean_isi
#> [1] 1.979794
#> 
#> $median_isi
#> [1] 1.990962
#> 
```
