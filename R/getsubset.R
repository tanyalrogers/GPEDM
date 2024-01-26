#' Subset training data using kmeans clustering
#'
#' Reduces the sizes of the training dataset for 'big data' GP applications. This
#' function is used internally in \code{fitGP} when \code{kmsubset=TRUE} (after
#' data scaling but before subsetting out missing data - rows not in the subset
#' are treated as 'missing'). However, this function can also be used on its own for
#' diagnostic purposes, or to pre-subset the training data if desired.
#' 
#' @details 
#' For a given training dataset, uses kmeans to divide the predictor rows into
#' clusters, and then subsamples (without replacement) from each of those clusters. 
#' The goal is to select points that are reasonably distributed across the attractor,
#' which appears to work better than points selected completely at random (nclust=1). 
#' 
#' The number of clusters must be specified, along with either the number of points per 
#' cluster or the total size of the subsample. 
#' 
#' If you want to use this function, note that xds (x) must be a matrix and the predictors
#' must be scaled. There can be missing values: the clustering and subsampling is only done
#' for complete rows (no missing x or y values). In other words, the indices that come out
#' will remove both incomplete rows and the rows not in the subset. 
#' 
#' @param xds Predictor variables as a matrix. These should be scaled! Can contain
#'   missing values.
#' @param yds Response variable as a vector. This is only used to exclude rows 
#'   with missing values for y, it is not used for clustering.
#' @param nclust The number of clusters.
#' @param clustsize The number of observations per cluster.
#' @param subtotal The total number of observations in the subset. The 
#'   number of observations per cluster will be subtotal/nclust. If subtotal
#'   is not an exact multiple of nclust, then the remainder is assigned to
#'   random clusters.
#'   
#' @return A vector of indices indicating the rows of xds and yds that should be included,
#'   which can be used for subsetting.
#' @export
#' @keywords functions
getsubset=function(xds, yds, nclust, clustsize=NULL, subtotal=NULL) {
  
  if(!is.null(clustsize) & !is.null(subtotal)) {
    stop("Provide either clustsize or subtotal, not both")
  }
  if(is.null(clustsize) & is.null(subtotal)) {
    stop("either clustsize or subtotal required")
  }
  
  completerows=complete.cases(cbind(yds,xds))
  X=as.matrix(xds[completerows,,drop=FALSE])

  
  clust=stats::kmeans(X,nclust)$cluster
  subrows=numeric()
  
  if(!is.null(clustsize)) { #take clustsize values from each cluster
    
    if(nrow(X)<nclust*clustsize) {
      stop("Number of non-missing data points is less than nclust*clustsize")
    }
    
    for(i in 1:nclust) {
      subrows=c(subrows,sample(which(completerows)[clust==i], size=clustsize, replace=FALSE))
    }
  }
  if(!is.null(subtotal)) { #take subtotal/nclust values from each cluster  
    
    if(nrow(X)<subtotal) {
      stop("Number of non-missing data points is less than subtotal")
    }

    clustsizevec=rep(floor(subtotal/nclust),times=nclust)
    rem=sample(1:nclust, subtotal %% nclust) #add remainder to random clusters
    clustsizevec[rem]=clustsizevec[rem]+1
    for(i in 1:nclust) {
      subrows=c(subrows,sample(which(completerows)[clust==i], size=clustsizevec[i], replace=FALSE))
    }
  }
  
  clustvec=yds*NA
  clustvec[completerows]=clust
  
  return(list(subrows=sort(subrows), clusters=clustvec))
}