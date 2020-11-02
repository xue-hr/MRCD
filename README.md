
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

The list **pruned** contains two objects: **loci\_bed** is the reference
panel of 22 SNPs for 489 individuals, **sig\_part** is the summary
statistics for LDL and CAD:

``` r

pruned$loci_bed[1:6,1:6]
#>                 rs11206510 rs11591147 rs12740374 rs4245791 rs2351524 rs2282143
#> HG00096:HG00096          0          0          2         0         0         0
#> HG00097:HG00097          1          0          0         1         1         0
#> HG00099:HG00099          1          0          0         1         0         0
#> HG00101:HG00101          0          0          0         1         1         0
#> HG00102:HG00102          0          0          0         1         0         0
#> HG00103:HG00103          2          0          0         0         0         0
head(pruned$sig_part)
#>   chr       pos       rsid A1 A2   beta_LDL   se_LDL  N_LDL beta_CAD  se_CAD
#> 1   1  55496039 rs11206510  C  T -0.0695200 0.003555 294565 -0.06272 0.01138
#> 2   1  55505647 rs11591147  T  G -0.4752703 0.011494 265213 -0.22090 0.03500
#> 3   1 109817590 rs12740374  T  G -0.1618858 0.003197 294565 -0.10990 0.01022
#> 4   2  44074431  rs4245791  C  T  0.0721910 0.003093 270962  0.05444 0.00887
#> 5   2 203880992  rs2351524  T  C -0.0239770 0.004204 295826  0.10227 0.01228
#> 6   6 160557643  rs2282143  T  C  0.0582204 0.008825 247909  0.26186 0.03325
#>    N_CAD loci
#> 1 336860   34
#> 2 268736   34
#> 3 268733   67
#> 4 268741  160
#> 5 268745  253
#> 6 243575  734
```

Apply 3 CD methods to 22 SNPs from 12 loci:

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

Now apply 3 CD methods to 12 indpendent SNPs, each of them from a locus.
We do not need reference panel for independent SNPs.

``` r
set.seed(1)
sig_part = pruned$sig_part[c(2,3,4,5,6,9,11,13,14,15,18,21),]
independent_results = CD_3_methods_Independent(sig_part)
```

Results of CD-Ratio:

``` r
independent_results$CD_Ratio_result
#> $T1toT2
#>          K      se(K) 
#> 0.19413572 0.01050562 
#> 
#> $T2toT1
#>          K      se(K) 
#> 0.72765268 0.06337732 
#> 
#> $Q_T1toT2
#>          [,1]
#> [1,] 184.0483
#> 
#> $Q_T2toT1
#>          [,1]
#> [1,] 393.7104
```

Results of CD-Egger:

``` r
independent_results$CD_Egger_result
#> $T1toT2
#>          b0           K      se(b0)       se(K) 
#> 0.008331157 0.220582429 0.002173359 0.040930430 
#> 
#> $T2toT1
#>          b0           K      se(b0)       se(K) 
#> -0.02831589  3.21230872  0.00918193  0.59210624 
#> 
#> $Q_T1toT2
#>          [,1]
#> [1,] 11.93857
#> 
#> $Q_T2toT1
#>          [,1]
#> [1,] 11.99474
```

Results of CD-GLS:

``` r
independent_results$CD_GLS_result
#> $T1toT2
#>          b0           K      se(b0)       se(K) 
#> 0.008162980 0.229319510 0.002289542 0.041957061 
#> 
#> $T2toT1
#>           b0            K       se(b0)        se(K) 
#> -0.028377515  3.005738170  0.008756814  0.521305169 
#> 
#> $Q_T1toT2
#>          [,1]
#> [1,] 12.12634
#> 
#> $Q_T2toT1
#>          [,1]
#> [1,] 12.84648
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
