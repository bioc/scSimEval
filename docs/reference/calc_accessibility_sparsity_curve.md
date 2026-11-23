# Fit Chromatin Accessibility-Sparsity Polynomial Curve (simATAC)

Fits a polynomial curve relating peak/bin mean accessibility to non-zero
cell proportion (NZP / detection frequency) as modeled in the simATAC
framework (Navidi et al., Genome Biology 2021): NZP = c0 + c1 \* mean +
c2 \* mean^2.

## Usage

``` r
calc_accessibility_sparsity_curve(data, poly_degree = 2)
```

## Arguments

- data:

  Count matrix (features/peaks x cells) or data.frame with 'mean' and
  'nzp'.

- poly_degree:

  Degree of polynomial (default: 2 for quadratic curve).

## Value

A named list summarizing the polynomial fit relating peak mean
accessibility to non-zero proportion (NZP):

- `c0`:

  Estimated intercept coefficient.

- `c1`:

  Estimated linear term coefficient.

- `c2`:

  Estimated quadratic term coefficient (or `0` if `poly_degree < 2`).

- `r_squared`:

  Coefficient of determination (\\R^2\\) of the polynomial regression.

- `spearman_cor`:

  Spearman correlation between mean accessibility and NZP.

- `pearson_cor`:

  Pearson correlation between mean accessibility and NZP.

- `poly_degree`:

  Polynomial degree used in the fit.

- `model`:

  The fitted `lm` object (omitted when too few valid peaks exist).

## Examples

``` r
data(example_scatac, package = "scSimEval")
calc_accessibility_sparsity_curve(example_scatac$ref)
#> $c0
#> [1] 0.02833341
#> 
#> $c1
#> [1] 0.3916258
#> 
#> $c2
#> [1] -0.08619656
#> 
#> $r_squared
#> [1] 0.8882643
#> 
#> $spearman_cor
#> [1] 0.9212479
#> 
#> $pearson_cor
#> [1] 0.9184381
#> 
#> $poly_degree
#> [1] 2
#> 
#> $model
#> 
#> Call:
#> stats::lm(formula = y ~ stats::poly(x, degree = eff_deg, raw = TRUE), 
#>     data = df)
#> 
#> Coefficients:
#>                                   (Intercept)  
#>                                       0.02833  
#> stats::poly(x, degree = eff_deg, raw = TRUE)1  
#>                                       0.39163  
#> stats::poly(x, degree = eff_deg, raw = TRUE)2  
#>                                      -0.08620  
#> 
#> 
```
