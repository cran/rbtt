# rbtt (Robust bootstrapped t-test) 
[![Travis-CI Build Status](http://travis-ci.org/WannabeSmith/rbtt.svg?branch=master)](http://travis-ci.org/WannabeSmith/rbtt)

In datasets whose data-generating distributions are non-negative with excess zero observations, it can be difficult to find general-purpose statistical tests for comparing sample means while controlling type-I error rates. This R package allows users to perform a modified bootstrap-based t-test that aims to better control type-I error rates in these particular datasets.

To download and use this package, run the following:

```r
install.packages("devtools") # Unless you already have it
library(devtools)
devtools::install_github("WannabeSmith/rbtt")
```
To obtain details on how to use the package, run:

```r
help(rbtt)
```
