
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
also produce predictions if desired. Use `summary.GP` to view a summary,
`predict.GP` to generate other or additional predictions from a fitted
model, `plot.GP` to plot observed and predicted values and
`getconditionals` to obtain conditional reponses. Also available for use
is the function `makelags` which can be used to create delay vectors.
See the (not yet created) vignette for more detailed instructions.

## Installation

To install the package:

``` r
install.packages("devtools") #if required
devtools::install_github("tanyalrogers/GPEDM")
```

If you are using Windows, you may need to install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), which allows
you to complile the C++ code.

## Example

Here are some simulated data from 2 populations with theta logistic
dynamics. The data contain some small lognormal process noise, and the
populations have different theta values.

``` r
library(GPEDM)
data("thetalog2pop")

pA=subset(thetalog2pop,Population=="PopA")
pB=subset(thetalog2pop,Population=="PopB")
N=nrow(pA)
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(Abundance~Time,data=pA,type="l",main="PopA")
plot(Abundance~Time,data=pB,type="l",main="PopB")
plot(pA$Abundance[1:(N-1)],pA$Abundance[2:N],
     xlab="Abundance t",ylab="Abundance t+1",main="PopA")
plot(pB$Abundance[1:(N-1)],pB$Abundance[2:N],
     xlab="Abundance t",ylab="Abundance t+1",main="PopB")
```

<img src="man/figures/README-data-1.png" width="100%" />

Here is how you might set up a hierarchical time-delay embedding model.
In this example, `Abundance` is the response variable (`yd`). We will
use an embedding dimension (`E`) of 3 and time delay (`tau`) of 1.
`Population` indicates which population the data are from, so it is
included under `pop`. Since the data are on somewhat different scales
and don’t necessarily represent the same ‘units’, we will use local
(within population) data scaling, as opposed to global. Just for fun, we
will also request leave-one-out predictions.

``` r
tlogtest=fitGP(data = thetalog2pop, yd = "Abundance", pop = "Population", E=3, tau=1, 
               scaling = "local", predictmethod = "loo")
summary(tlogtest)
#> Number of predictors: 3 
#> Abundance_1 Abundance_2 Abundance_3 
#> Length scale parameters:
#>          phi1          phi2          phi3 
#>  5.290040e-01 1.248172e-216 1.248172e-216 
#> Observation variance (ve): 0.01221358
#> Function variance (tau): 2.53986
#> Number of populations: 2
#> Dynamic correlation (rho): 0.3250464
#> In-sample R-squared: 0.9939452
#> Out-of-sample R-squared: 0.9919687
```

From the summary, we can see that ARD has (unsurprisingly) deemed lags 2
and 3 to be unimportant (length scales are 0), so E=1 is probably
sufficient. The dynamic correlation (rho) tells us the degree to which
the dynamics correlated (they are are rather dissimilar in this case).
Since the simulated data don’t contain that much noise, the R-squared
values are quite high.

We can examine the observed and predicted values using `plot`. Standard
error bands are included in the time series plots, although they’re a
little hard to see in this example.

``` r
plot(tlogtest)
#> Plotting out of sample results.
```

<img src="man/figures/README-plot1-1.png" width="100%" />

To get the in-sample predictions:

``` r
plot(tlogtest, plotinsamp = T)
#> Plotting in sample results.
```

<img src="man/figures/README-plot2-1.png" width="100%" />

The function `getconditionals` will compute and plot conditional
responses to each input variable (other input varibles set to their mean
value). From this we can also clearly see that lags 2 and 3 have no
impact, and we can see how the lag 1 dynamics of the 2 populations
differ.

``` r
con=getconditionals(tlogtest)
```

<img src="man/figures/README-conditionals-1.png" width="100%" />

We can use the `predict` function to get other types of predictions.

``` r
#sequential predictions (they should improve over time)
seqpred=predict(tlogtest,predictmethod = "sequential")
### This is broken ###
#plot(seqpred) #have to specify plot.GP for output of predict.GP, but not finding?
```

You could also supply new data for which to make predictions. In this
case, you would supply a new data frame (`datanew`), which should
contain columns `Abundance` and `Population`.

## Specifying training data

There are 4 ways that the training data for a model can be specified.
The above example is method 2. Method 1 allows for the most
customization of response and predictor variables including mixed
embeddings, inclusion of covariates, etc.

1.  data frame `data`, plus column names or indices for `yd` and `xd`  
2.  data frame `data`, plus column name or index for `yd`, and values
    for `E` and `tau`  
3.  vector `yd`, plus matrix or vector `xd`  
4.  vector `yd`, plus values for `E` and `tau`

The `pop` argument is optional in all of the above cases. If omitted, a
single population is assumed.

``` r
set.seed(10)
thetalog2pop$othervar=rnorm(nrow(thetalog2pop))
yvec=thetalog2pop$Abundance
popvec=thetalog2pop$Population
#function 'makelags' can be used to generate a lag matrix
#be sure to include 'pop' if data contain multiple pops to prevent crossover
xmat=makelags(yd=thetalog2pop[,c("Abundance","othervar")],pop=popvec,E=2,tau=1)
thetalog2pop2=cbind(thetalog2pop,xmat)

#Method 1
m1=fitGP(data=thetalog2pop2,yd="Abundance",xd=c("Abundance_1","othervar"),
         pop="Population",scaling="local")

#Method 3
#like 1, but your data aren't in a data frame
m3=fitGP(yd=yvec,xd=xmat,pop=popvec,scaling="local")

#Method 4
#like 2, but your data aren't in a data frame
m4=fitGP(yd=yvec,pop=popvec,E=2,tau=1,scaling="local")

summary(m1)
#> Number of predictors: 2 
#> Abundance_1 othervar 
#> Length scale parameters:
#>         phi1         phi2 
#> 0.5207983744 0.0003021264 
#> Observation variance (ve): 0.01031384
#> Function variance (tau): 2.539001
#> Number of populations: 2
#> Dynamic correlation (rho): 0.2929397
#> In-sample R-squared: 0.9947297
summary(m3)
#> Number of predictors: 4 
#> Abundance_1 Abundance_2 othervar_1 othervar_2 
#> Length scale parameters:
#>         phi1         phi2         phi3         phi4 
#> 5.090030e-01 2.128391e-19 2.752429e-04 1.830083e-44 
#> Observation variance (ve): 0.01026833
#> Function variance (tau): 2.671209
#> Number of populations: 2
#> Dynamic correlation (rho): 0.2873044
#> In-sample R-squared: 0.9948366
summary(m4)
#> Number of predictors: 2 
#> Length scale parameters:
#>          phi1          phi2 
#>  5.302715e-01 3.355234e-173 
#> Observation variance (ve): 0.01205433
#> Function variance (tau): 2.53792
#> Number of populations: 2
#> Dynamic correlation (rho): 0.3250085
#> In-sample R-squared: 0.9938883
```

## References

Munch, S. B., Poynor, V., and Arriaza, J. L. 2017. Circumventing
structural uncertainty: a Bayesian perspective on nonlinear forecasting
for ecology. Ecological Complexity, 32:134.

*Any advice on improving this package is appreciated.*
