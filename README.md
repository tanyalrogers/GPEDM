
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GPEDM

<!-- badges: start -->
<!-- badges: end -->

**Disclaimer: This is a very much a work in progress. Use at your own
risk.**

This package contains functions for fiting hierarchical, separable
length scale GP models with automatic relevance determination (ARD) for
use in Empirical Dynamic Modeling (EDM) and other applications. This is
an adaptation of code originally developed by Stephan Munch in MATLAB.

The main function is `fitGP` which is used to train the model and can
also produce predictions if desired. Use `predict.GP` to generate other
or additional predictions from a fitted model, `plot.GP` to plot
observed and predicted values and `getconditionals` to obtain
conditional reponses. Also available for use is the function `makelags`
which can be used to create delay vectors. See the (not yet created)
vignette for more detailed instructions.

## Installation

To install the package:

``` r
devtools::install_github("tanyalrogers/GPEDM")
```

If you are using Windows, you may need to install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).

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
#see a summary
summary(test)
#> Number of predictors: 2
#> Length scale parameters:
#>          phi1          phi2 
#>  4.625908e-01 1.727341e-303 
#> Observation variance (ve): 0.06311373
#> Function variance (tau): 3.235288
#> Number of populations: 1
#> In-sample R-squared: 0.9568493
#> Out-of-sample R-squared: 0.9026721
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
summary(testhier)
#> Number of predictors: 2
#> Length scale parameters:
#>      phi1      phi2 
#> 0.2656263 0.0000000 
#> Observation variance (ve): 0.07148576
#> Function variance (tau): 3.027286
#> Number of populations: 2
#> Dynamic correlation (rho): 0.901183
#> In-sample R-squared: 0.9578276
#> Out-of-sample R-squared: 0.8706187
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

## References

Munch, S. B., Poynor, V., and Arriaza, J. L. 2017. Circumventing
structural uncertainty: a Bayesian perspective on nonlinear forecasting
for ecology. Ecological Complexity, 32:134.

*Any advice on improving this package is appreciated.*
