John W Pickering
7 March 2023

<!-- README.md is generated from README.Rmd. Please edit that file -->

# rap

<!-- badges: start -->
<!-- badges: end -->

The rap package contains functions for generating statistical metrics
and visual means to assess the improvement in risk prediction of one
risk model over another. It includes the Risk Assessment Plot (hence
rap).

## Installation

You can install the development version or rap from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JohnPickering/rap")
```

## History and versions

rap began as Matlab code in 2012 after I wrote a paper
(<a href="https://cjasn.asnjournals.org/content/7/8/1355"
target="_blank">1</a>) for the Nephrology community on assessing the
added value of one biomarker to a clinical prediction model. I worked
with Professor Zoltan Endre on that paper. Dr David Cairns kindly
provided some R code for the Risk Assessment Plot. This formed the basis
of versions 0.1 to 0.4. Importantly, for those versions and the current
version all errors are mine (sorry) and not those of Professor Endre or
Dr Cairns. Since writing that paper I’ve come to consider some metrics
as not helpful. So, for the current version I have dropped some
statistical metrics that I believe are poor or wrongly applied. In
particularly, I dropped providing the total NRI (Net Reclassification
Improvement) and total IDI (Integrated Discrimination Improvement)
metrics. These should never be presented because they inappropriately
add together two fractions with differing denominators (NRI) or two
means (IDI). Instead, these the NRIs and IDIs for those with and without
the event of interest should be provided. Third, I have provided the
change in AUCs rather than a p-value because the change is much more
meaningful.

Version 1.03 were major changes:  
\* allowed as input logistic regression models from glm (stats) and lrm
(rms) as well as risk predictions calculated elsewhere. \* provided as
outputs in addition to the Risk Assessment Plot, a form of calibration
plot and decision curve.  
\* the output from the main functions CI.raplot, and CI.classNRI are now
lists that include the metrics for each bootstrap sample as well as the
summary metrics. CI.classNRI also produces confusion matrices for those
with and without the event of interest (separately). Bootstrapping is
used to determine confidence intervals.

Version 1.10 (current):  
\* addition of ROC plot. \* calibration plot now uses (best practice)
continuous curves (the old format is now “ggcalibrate_original()”).  
\* addition of precision recall curves. \* all plots can be for one or
two models.

## Example 1

This is a basic example for assessing the difference between two
logistic regression models:

``` r
library(rap)
## basic example code

baseline_risk <- data_risk$baseline    # or the baseline glm model itself
new_risk <- data_risk$new              # or the new glm model itself
outcome <- data_risk$outcome

assessment <- CI.raplot(x1 = baseline_risk, x2 = new_risk, y = outcome,
                        n.boot = 20, dp = 2) # Note the default is 1000 bootstraps (n.boot = 1000).  This can take quite some time to run, so when testing I use a smaller number of bootstraps.  

# View results  
## meta data  
(assessment$meta_data)
#>   Thresholds Confidence.interval Number.of.bootstraps Input.data.type
#> 1   baseline                  95                   20   User supplied
#>   X..decimal.places
#> 1                 2

## exact point estimates  
(assessment$Metrics)
#> $n
#> [1] 433
#> 
#> $n_event
#> [1] 86
#> 
#> $n_non_event
#> [1] 347
#> 
#> $Prevalence
#> [1] 0.1986143
#> 
#> $NRI_up_event
#> [1] 20
#> 
#> $NRI_up_nonevent
#> [1] 19
#> 
#> $NRI_down_event
#> [1] 13
#> 
#> $NRI_down_nonevent
#> [1] 79
#> 
#> $NRI_event
#> [1] 0.08139535
#> 
#> $NRI_nonevent
#> [1] 0.1729107
#> 
#> $IDI_event
#> [1] 0.1363479
#> 
#> $IDI_nonevent
#> [1] 0.03397117
#> 
#> $IP_baseline
#> [1] 0.1849132
#> 
#> $IS_baseline
#> [1] 0.2516693
#> 
#> $IP_new
#> [1] 0.1504058
#> 
#> $IS_new
#> [1] 0.388494
#> 
#> $Brier_baseline
#> [1] 0.1500246
#> 
#> $Brier_new
#> [1] 0.123228
#> 
#> $Brier_skill
#> [1] 17.86145
#> 
#> $AUC_baseline
#> [1] 0.6823772
#> 
#> $AUC_new
#> [1] 0.8227331
#> 
#> $AUC_difference
#> [1] 0.1403559

## bootstrap derived metrics with confidence intervals  
(assessment$Summary_metrics)
#> # A tibble: 22 × 2
#>    metric            statistics                  
#>    <chr>             <chr>                       
#>  1 n                 432 (CI: 427.9 to 440.52)   
#>  2 n_event           88.5 (CI: 71.22 to 97.57)   
#>  3 n_non_event       345.5 (CI: 333.42 to 361.15)
#>  4 Prevalence        0.2 (CI: 0.16 to 0.23)      
#>  5 NRI_up_event      19 (CI: 10.8 to 30.72)      
#>  6 NRI_up_nonevent   19.5 (CI: 11.95 to 27.57)   
#>  7 NRI_down_event    11 (CI: 6.48 to 18.1)       
#>  8 NRI_down_nonevent 69.5 (CI: 51.9 to 101.07)   
#>  9 NRI_event         0.09 (CI: 0 to 0.23)        
#> 10 NRI_nonevent      0.15 (CI: 0.08 to 0.22)     
#> # … with 12 more rows
```

## Graphical assessments

### The Risk Assessment Plot

``` r
ggrap(x1 = baseline_risk, x2 = new_risk, y = outcome)
```

<img src="man/figures/README-ggrap-1.png" width="100%" />

``` r

# for Single risks x2 = NULL
```

### The calibration curve

``` r
ggcalibrate(x1 = baseline_risk, x2 = new_risk, y = outcome)
```

<img src="man/figures/README-ggcalibrate-1.png" width="100%" />

### The decision curve

``` r
ggdecision(x1 = baseline_risk, x2 = new_risk, y = outcome)
```

<img src="man/figures/README-ggdecision-1.png" width="100%" />

### The precission-recall curve

``` r
ggprerec(x1 = baseline_risk, x2 = new_risk, y = outcome)
```

<img src="man/figures/README-ggrerec-1.png" width="100%" />

### The roc plot

``` r
ggroc(x1 = baseline_risk, x2 = new_risk, y = outcome)
```

<img src="man/figures/README-ggroc-1.png" width="100%" />

Note, there are additional options for the ROC plot including labelling
points and distinguishing areas of the plot that are diagnostic from
those that are not.

## Example 2

This is a basic example for assessing the difference in the results of
reclassification:

``` r
## basic example code

baseline_class <- data_class$base_class
new_class <- data_class$new_class
outcome_class <- data_class$outcome

class_assessment <- CI.classNRI(c1 = baseline_class, c2 = new_class, y = outcome_class, 
                                n.boot = 20, dp = 2) # Note the default is 2000 bootstraps (n.boot = 2000).  This can take quite some time to run, so when testing I use a smaller number of bootstraps.  

# View results  
## meta data  
(class_assessment$meta_data)
#>   Confidence.interval Number.of.bootstraps X..decimal.places
#> 1                  95                   20                 2

## exact point estimates and confusion matrices
(class_assessment$Metrics)
#> $n
#> [1] 444
#> 
#> $n_event
#> [1] 62
#> 
#> $n_non_event
#> [1] 382
#> 
#> $Prevalence
#> [1] 0.1396396
#> 
#> $NRI_up_event
#> [1] 21
#> 
#> $NRI_up_nonevent
#> [1] 94
#> 
#> $NRI_down_event
#> [1] 5
#> 
#> $NRI_down_nonevent
#> [1] 71
#> 
#> $NRI_event
#> [1] 0.2580645
#> 
#> $NRI_nonevent
#> [1] -0.06020942
#> 
#> $wNRI_event
#> NULL
#> 
#> $wNRI_nonevent
#> NULL
#> 
#> $confusion.matrix_event
#>          New
#> Baseline  class_1 class_2 class_3 class_4 class_5 class_6
#>   class_1       0       0       0       0       0       0
#>   class_2       0       2       3       0       0       0
#>   class_3       0       1       0       2       0       0
#>   class_4       0       0       2       9      11       1
#>   class_5       0       0       0       2      25       4
#>   class_6       0       0       0       0       0       0
#> 
#> $confusion.matrix_nonevent
#>          New
#> Baseline  class_1 class_2 class_3 class_4 class_5 class_6
#>   class_1       0       0       0       0       0       0
#>   class_2       9      52      17       3       1       0
#>   class_3       2      29      66      44       3       0
#>   class_4       0       2      21      52      25       0
#>   class_5       0       0       0       8      47       1
#>   class_6       0       0       0       0       0       0

## bootstrap derived metrics with confidence intervals  
(class_assessment$Summary_metrics)
#> # A tibble: 10 × 2
#>    metric            statistics                 
#>    <chr>             <chr>                      
#>  1 n                 444 (CI: 444 to 444)       
#>  2 n_event           62.5 (CI: 51.38 to 75.1)   
#>  3 n_non_event       381.5 (CI: 368.9 to 392.62)
#>  4 Prevalence        0.14 (CI: 0.12 to 0.17)    
#>  5 NRI_up_event      22 (CI: 14.48 to 31.57)    
#>  6 NRI_up_nonevent   94.5 (CI: 78.8 to 103.1)   
#>  7 NRI_down_event    5 (CI: 1.48 to 8)          
#>  8 NRI_down_nonevent 73 (CI: 59.42 to 83.62)    
#>  9 NRI_event         0.27 (CI: 0.16 to 0.42)    
#> 10 NRI_nonevent      -0.05 (CI: -0.1 to -0.01)
```
