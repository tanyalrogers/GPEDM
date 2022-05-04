#' Simulated time series from 2 theta-logistic populations
#'
#' Contains time series simulated from the discrete time theta-logisic model with noise:
#' \deqn{x_{t+1}=x_{t}e^{r(1-x_t^\theta)}e^{s\epsilon}}
#' \deqn{\epsilon \sim N(0,1)}
#' PopA has \eqn{\theta=3.5} and PopB has \eqn{\theta=2.5}. 
#' Both have \eqn{r=1.2} and \eqn{s=0.05} and time series length of 50, 
#' with random starting values and the first 200 iterates discarded.
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{Time}{Time index}
#'   \item{Population}{Population ID}
#'   \item{Abundance}{Abundance value}
#' }
#' @source The simulation code is in data-raw on the GPEDM Github page.
#' @keywords datasets
"thetalog2pop"

#' Simulated time series from 3 species Hastings-Powell model
#'
#' Contains time series simulated from the continuous time 3 species Hastings-Powell model:
#' \deqn{\frac{dX}{dt}=X(1-X)-aXY/(1+bX)}
#' \deqn{\frac{dY}{dt}=aXY/(1+bx)-cYZ/(1+dY)-mY}
#' \deqn{\frac{dZ}{dt}=cYZ/(1+dY)-\mu Z}
#' Parameters are \eqn{a=5, b=3, c=0.1, d=2, m=0.4, \mu=0.01}.  
#' Initial conditions are \eqn{X(0)=0.5, Y(0)=0.1, Z(0)=9}.  
#' Integrated for 500 time steps. Both have \eqn{r=1.2} and \eqn{s=0.05} and time series length of 50, 
#' with random starting values and the first 200 iterates discarded.
#'
#' @format A data frame with 501 rows and 4 variables:
#' \describe{
#'   \item{Time}{Time index}
#'   \item{X}{Abundance of species X}
#'   \item{Y}{Abundance of species Y}
#'   \item{Z}{Abundance of species Z}
#' }
#' @source The simulation code is in data-raw on the GPEDM Github page.
#' @keywords datasets
"HastPow3sp"