
# MiceLmer

<!-- badges: start -->
<!-- badges: end -->

The goal of MiceLmer is to get pooled estimates of random effects from lmer with multiple imputation based on mice.
You can apply these functions in both two-level and three-level multilevel models.

## Installation

You can install the development version of MiceLmer from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("toro-maguro/MiceLmer")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MiceLmer)
## ExtractRandomEffect(imp, "y ~ 1 + (1 | cluster_name)", "cluster_name")
## PlotRandomEffectDistribution(df_ExtractRandomEffect, "cluster_name")
## GetProportionOfRandomEffect(df_ExtractRandomEffect)
```

