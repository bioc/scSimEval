# Calculate Posterior Excess Zero Weights (zingeR & ZINB-WaVE)

Calculates the posterior probability of zero counts being technical
dropouts (excess zeros) versus biological sampling zeros under a
Negative Binomial model.

## Usage

``` r
calc_excess_zero_weights(counts)
```

## Arguments

- counts:

  Count matrix (features x cells).

## Value

A named list summarizing posterior excess-zero probabilities under a
Negative Binomial model:

- `mean_excess_zero_weight`:

  Mean posterior excess-zero weight across features.

- `median_excess_zero_weight`:

  Median posterior excess-zero weight across features.

- `mean_zero_inflation`:

  Mean estimated excess-zero fraction (`pi`) beyond Negative Binomial
  expectations.

- `gene_excess_weights`:

  Numeric vector of per-feature posterior excess-zero weights.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_excess_zero_weights(example_scrna$ref)
#> $mean_excess_zero_weight
#> [1] 0.02805731
#> 
#> $median_excess_zero_weight
#> [1] 0
#> 
#> $mean_zero_inflation
#> [1] 0.0222688
#> 
#> $gene_excess_weights
#>   [1] 0.000000000 0.134244265 0.102636750 0.000000000 0.000000000 0.000000000
#>   [7] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [13] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [19] 0.000000000 0.315345876 0.000000000 0.000000000 0.000000000 0.000000000
#>  [25] 0.001378365 0.000000000 0.125474321 0.000000000 0.000000000 0.041770849
#>  [31] 0.000000000 0.002636009 0.093152687 0.000000000 0.000000000 0.000000000
#>  [37] 0.024948507 0.000000000 0.021477136 0.000000000 0.000000000 0.074772269
#>  [43] 0.000000000 0.000000000 0.000000000 0.000000000 0.101572086 0.093995210
#>  [49] 0.134377026 0.051553186 0.000000000 0.000000000 0.001689478 0.219156928
#>  [55] 0.000000000 0.069621654 0.000000000 0.000000000 0.000000000 0.037087678
#>  [61] 0.000000000 0.029125437 0.000000000 0.032197006 0.000000000 0.000000000
#>  [67] 0.000000000 0.000000000 0.000000000 0.024641412 0.000000000 0.000000000
#>  [73] 0.107666417 0.000000000 0.000000000 0.046730259 0.000000000 0.043869620
#>  [79] 0.000000000 0.132104547 0.000000000 0.000000000 0.114589296 0.000000000
#>  [85] 0.031596744 0.000000000 0.000000000 0.044165430 0.000000000 0.000000000
#>  [91] 0.133689928 0.000000000 0.000000000 0.000000000 0.169601921 0.000000000
#>  [97] 0.000000000 0.123365768 0.000000000 0.125496974
#> 
```
