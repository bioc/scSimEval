# Calculate Neighborhood Purity

Integrated from simPIC (Chugh et al., 2024; bluster::neighborPurity).
For each cell, evaluates the proportion of its k-nearest neighbors that
share the same cluster or cell-type identity.

## Usage

``` r
calc_neighborhood_purity(
  data,
  cluster_labels,
  k = NULL,
  is_distance = FALSE,
  threads = 1L
)
```

## Arguments

- data:

  Matrix of coordinates (cells x dimensions) or expression matrix
  (features x cells) or dist matrix.

- cluster_labels:

  Vector of cluster or cell-type labels.

- k:

  Number of nearest neighbors (default 5% of cells).

- is_distance:

  Logical, whether data is already a distance matrix. Default FALSE.

- threads:

  Integer number of CPU threads (default 1L).

## Value

A numeric vector of neighborhood purity scores (values from 0 to 1),
with an attribute "mean_purity" containing the overall average.

## Examples

``` r
data <- matrix(stats::rnorm(200), 20, 10)
cl <- factor(rep(c("A", "B"), each = 10))
dist_mat <- stats::dist(data)
calc_neighborhood_purity(data, cl)
#>  [1] 0.5714286 0.7142857 0.5000000 0.5000000 0.5714286 0.6666667 0.6666667
#>  [8] 0.5000000 0.5000000 0.6000000 0.7142857 0.6000000 0.3750000 0.6000000
#> [15] 0.6666667 0.8000000 0.6000000 0.3333333 0.5000000 0.7500000
#> attr(,"mean_purity")
#> [1] 0.5864881
```
