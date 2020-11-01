
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRCD

<!-- badges: start -->

<!-- badges: end -->

This package is for performing CD methods in the paper “Inferring causal
direction between two traits in the presence of horizontal pleiotropy
with GWAS summary
data”.

## Installation

<!-- You can install the released version of MRCD from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MRCD")
```
-->

Install development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xue-hr/MRCD")
```

## Example

We use LDL and CAD as an example.

``` r
library(MRCD)
## basic example code
```

Apply all 3 CD methods:

``` r
set.seed(1)
results = CD_3_methods(pruned)
```

Results of CD-Ratio:

``` r
results$CD_Ratio_result
#> $T1toT2
#>           K       se(K) 
#> 0.145074866 0.006803127 
#> 
#> $T2toT1
#>          K      se(K) 
#> 0.94092638 0.05439297 
#> 
#> $Q_T1toT2
#>         [,1]
#> [1,] 343.244
#> 
#> $Q_T2toT1
#>          [,1]
#> [1,] 498.7441
```

Results of CD-Egger:

``` r
results$CD_Egger_result
#> $T1toT2
#>          b0           K      se(b0)       se(K) 
#> 0.006055084 0.166677177 0.001990455 0.031729184 
#> 
#> $T2toT1
#>          b0           K      se(b0)       se(K) 
#> -0.02529651  3.29972149  0.00949699  0.63893717 
#> 
#> $Q_T1toT2
#>          [,1]
#> [1,] 22.49414
#> 
#> $Q_T2toT1
#>          [,1]
#> [1,] 22.00414
```

Results of CD-GLS:

``` r
results$CD_GLS_result
#> $T1toT2
#>          b0           K      se(b0)       se(K) 
#> 0.005959528 0.166556725 0.001965053 0.031333548 
#> 
#> $T2toT1
#>          b0           K      se(b0)       se(K) 
#> -0.02374490  3.11837442  0.00853066  0.58016321 
#> 
#> $Q_T1toT2
#>          [,1]
#> [1,] 22.72853
#> 
#> $Q_T2toT1
#>          [,1]
#> [1,] 23.96457
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:


```r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
