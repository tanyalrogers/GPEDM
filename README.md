
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GPEDM

<!-- badges: start -->
<!-- badges: end -->

**This is a very much a work in progress. Use at your own risk.**

This package contains functions for fiting hierarchical, separable
length scale GP models with automatic relevance determination (ARD) for
use in Empirical Dynamic Modeling (EDM) and other applications.

## Installation

To install the package:

``` r
devtools::install_github("tanyalrogers/GPEDM")
```

## Simple Example

These are truly trivial examples that I should replace with something
else. I was using them mainly for testing.

``` r
library(GPEDM)
set.seed(10)

#some trival data
l=20
x1=seq(-5,5,length.out = l)
x2=rnorm(l,0,1)
y=0.5*x1^2+rnorm(l)
xmat=cbind(x1,x2)
#create a missing value
ym=y
ym[5]=NA
#as a data frame
dftest=cbind.data.frame(response=y,responsem=ym,xmat,pop=rep(c(1,2),10),pop2=rep(c("site A","site B"),10))
dftest2=dftest[order(dftest$pop),]

#specify training data as vector and matrix
test=fitGP(yd=y,xd=xmat,scaling = "global",predictmethod = "loo")
#plot out of sample results
plot(test)
#> Plotting out of sample results.
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
#plot in sample results
plot(test, plotinsamp = T)
#> Plotting in sample results.
```

<img src="man/figures/README-example-2.png" width="100%" />

``` r
#plot conditional responses
con1=getconditionals(test, plot = T)
```

<img src="man/figures/README-example-3.png" width="100%" />

``` r
#specify training data with data frame, example with 2 populations and a missing data point
testhier=fitGP(data=dftest2,yd="responsem",xd=c("x1","x2"),pop="pop2",scaling="local", predictmethod="loo")
plot(testhier)
#> Plotting out of sample results.
```

<img src="man/figures/README-example-4.png" width="100%" />

``` r
cond2=getconditionals(testhier)
```

<img src="man/figures/README-example-5.png" width="100%" />

``` r
#specify training data as lags (E, tau)
yrand = rnorm(50,0,1)
testrand=fitGP(yd=yrand,E=2,tau=1)
#alternatively
yrandlags=makelags(yrand,E=2,tau=1)
testrand2=fitGP(yd=yrand,xd=yrandlags)
```

*Any advice on improving this package is appreciated.*
