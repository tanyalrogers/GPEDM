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

#' NEFSC Fall Bottom Trawl Data (subset)
#' 
#' Catch-per-unit effort data for 3 species from the NOAA Northeast 
#' Fisheries Science Center Fall Bottom Trawl Survey for 4 regions
#' from 1973 to 2016 (44 years).
#' 
#' @format A data frame with 176 rows and 5 variables:
#' \describe{
#'   \item{YEAR}{Time index}
#'   \item{REGION}{Region: GME = Gulf of Maine, GEO = George's Bank,
#'      SNE = Southern New England, MAB = Mid-Atlantic Bight. This is 
#'      formatted as a factor (north to south).}
#'   \item{shortfin_squid}{CPUE of northern shortfin squid (*Illex illecebrosus*)}
#'   \item{longfin_squid}{CPUE of longfin squid (*Loligo pealeii*)}
#'   \item{silver_hake}{CPUE of silver hake (*Merluccius bilinearis*)}
#' }
#' @source https://www.fisheries.noaa.gov/inport/item/22560
#' @keywords datasets
"trawl"

#' Brown shrimp in the Gulf of Mexico
#' 
#' Annual catch-per-unit effort data for brown shrimp in the Gulf of Mexico from
#' the Southeast Area Monitoring and Assessment Program (SEAMAP) Summer
#' Groundfish Trawl Survey for 9 different statistical zones from 1987 to 2019
#' (33 years).
#' 
#' @format A data frame with 297 rows and 3 variables:
#' \describe{
#'   \item{year}{Time index}
#'   \item{zone}{Statistical zone (region)}
#'   \item{cpue}{Catch per unit effort}
#' }
#' 
#' @source SEAMAP Summer Groundfish Trawl Survey
#' @references 
#' Tsai CH, Munch SB, Masi MD, and Pollack AG. 2022. Predicting
#' nonlinear dynamics of short-lived penaeid shrimp species in the Gulf of
#' Mexico. Canadian Journal of Fisheries and Aquatic Sciences.
#' https://doi.org/10.1139/cjfas-2022-0029
#' @keywords datasets
"shrimp"

#' Simulated time series from a ricker model with fishing
#'
#' Contains a time series simulated from a Ricker model with harvesting:
#' \deqn{B_{t+1}=S_{t}e^{r-S_{t}/K+\epsilon_{b}}}
#' \deqn{X_{t+1}=bB_{t+1}}
#' \deqn{u_{t+1}=((1-a)u_{t}+au_{t}(B_{t+1}-B_{t})/B_{t})(1+\epsilon_{u})}
#' \deqn{C_{t+1}=B_{t+1}u_{t+1}}
#' \deqn{S_{t+1}=B_{t+1}-C_{t+1}}
#' \deqn{\epsilon_{b} \sim N(0,\sigma_{b})}
#' \deqn{\epsilon_{u} \sim Unif(-\sigma_{u},\sigma_{u})}
#' Where \eqn{B} is biomass, \eqn{X} is a survey index if biomass with catchability \eqn{b}, \eqn{u}
#' is the exploitation rate (which changes in response to changing biomass according to parameter
#' \eqn{a}), \eqn{C} is catch, and \eqn{S} is uncaught biomass (escapement).
#' Parameter values are \eqn{r=3}, \eqn{K=1000}, \eqn{b=0.01}, \eqn{a=0.1}, 
#' \eqn{\sigma_{b}=0.1}, and \eqn{\sigma_{u}=0.2}. 
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{Time}{Time index}
#'   \item{CPUE_index}{Index of abundance, catch per unit effort}
#'   \item{Catch}{Harvest}
#' }
#' @source The simulation code is in data-raw on the GPEDM Github page.
#' @keywords datasets
"RickerHarvest"
