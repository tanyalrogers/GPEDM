#' Get interated predictions from a GP model (or MSY from a "fisheries model")
#'
#' Obtains iterated predictions for given a GP model. Can also be used to obtain
#' MSY estimates from a "fisheries model".
#' 
#' The predicted value is inserted into the first lag of the next timepoint, and
#' the other lags are shifted right by one. The method thus assumes that all time lags
#' are present, evenly spaced, and in order; that all timepoints to be predicted are 
#' evenly spaced; and that the response variable and time lags are in the
#' same units.
#' 
#' Using this method requires the use of \code{data} with pre-generated lags
#' (option A1 in \code{\link{fitGP}}), that a \code{time} column is specified,
#' that \code{newdata} has exactly the same columns as \code{data}. The names of
#' the columns containing the lags (which will be iteratively updated) can be
#' indicated under \code{xlags}; if omitted, all columns originally supplied
#' under \code{x} are assumed to be lags. If your model contains covariates,
#' (values that are not time lags of the response variable) you will need to
#' specify \code{xlags}, which should not include the covariates. You don't have
#' to do this for "fisheries models" though, because lags (m) and covariates
#' (h,z) are have already been specified.
#'
#' Some care needs to be taken in the construction of \code{newdata}, which is
#' what is used for prediction. Predictions will be made for as many timepoints
#' as are in \code{newdata}, starting with the first timepoint.  Only the first
#' row of the \code{xlags} columns (for each population) needs to be complete.
#' The subsequent rows can be NA, or have values in them - those values will be
#' overwritten as the model iterates. If you have covariates, values for those
#' need to be supplied for those at all timepoints.
#' 
#' \strong{Obtaining MSY from a fisheries model:}
#'
#' If you are fitting a "fisheries model", the values of h should be set to some
#' constant value for all time points in \code{newdata}. At each iteration, the
#' values of m will be updated with the predictions, and the quantity m-bh will
#' be computed to obtain the next prediction. If covariates (z) are present,
#' values should be supplied for them. The value of \code{xlags} is assumed to be m.
#' 
#' 
#' @param object Output from \code{fitGP}.
#' @param newdata Data frame containing the same columns supplied in the
#'   original model.
#' @param xlags Which column names are the lags? Is not supplied, assumed to be all columns
#'   originally supplied under \code{x}.  
#' @return A list (class GPpred) with the same elements as \code{\link{predict.GP}}.
#' @examples 
#' xdata=makelags(data = thetalog2pop, y="Abundance", pop="Population", time = "Time", E=3, tau=1, append=T)
#' xdatatrain=subset(xdata, Time<=40)
#' xdatafore=subset(xdata, Time>40)
#' tlog=fitGP(data = xdatatrain, y = "Abundance", x=c("Abundance_1","Abundance_2","Abundance_3"), 
#'            pop = "Population", time = "Time", scaling = "local")
#' prediter=predict_iter(tlog,xdatafore)
#' plot(prediter)
#' @export
#' @keywords functions
predict_iter=function(object,newdata,xlags=NULL) { 
  
  if(is.null(object$inputs$time_names)) {
    stop("Model must include `data` and `time` to use this function.")
  }
  if(!is.null(object$inputs$E)) {
    stop("Lags must be pre-generated to use this function. See option A1 in help(fitGP).")
  }
  
  #will need modification for fisheries models
  
  if(is.null(xlags)) { #if not provided, assume all x columns are lags
    xlags=object$inputs$x_names
  }
  
  #get times
  timename=object$inputs$time_names
  newtimes=unique(newdata[,timename])
  
  if(is.null(object$inputs$pop_names)) {
    popname="pop"
    newdata$pop=1
  } else {
    popname=object$inputs$pop_names
  }
  up=unique(newdata[,popname])
  
  pred=list()
  
  for(i in seq_along(newtimes)) {
    newdatai=newdata[newdata[,timename]==newtimes[i],]
    #get prediction
    pred[[i]]=predict(object, newdata=newdatai)$outsampresults
    preds=pred[[i]]
    #create new data row for each population
    #assumes all lags are present and evenly spaced
    if(i+1 <= length(newtimes)) {
      for(j in 1:length(up)){
        predsj=preds[preds$pop==up[j],]$predmean
        newdata[newdata[,timename]==newtimes[i+1] & newdata[,popname]==up[j],xlags[1]]=predsj
        if(length(xlags>1)) {
          newdata[newdata[,timename]==newtimes[i+1] & newdata[,popname]==up[j],xlags[2:length(xlags)]]=
            newdata[newdata[,timename]==newtimes[i] & newdata[,popname]==up[j],xlags[1:(length(xlags)-1)]]
        }
      }
    }
  }
  #compile predictions
  outsampresults=do.call(rbind, pred)
  out=list(outsampresults=outsampresults)
  
  #return updated model and predictions
  if(object$inputs$y_names %in% colnames(newdata)) {
    if(!all(is.na(outsampresults$obs))) {
      out$outsampfitstats=c(R2=getR2(outsampresults$obs,outsampresults$predmean), 
                            rmse=sqrt(mean((outsampresults$obs-outsampresults$predmean)^2,na.rm=T)))
      if(length(unique(outsampresults$pop))>1) { #within site fit stats
        up=unique(outsampresults$pop)
        np=length(up)
        R2pop<-rmsepop<-numeric(np)
        names(R2pop)=up
        names(rmsepop)=up
        for(k in 1:np) {
          ind=which(outsampresults$pop==up[k])
          R2pop[k]=getR2(outsampresults$obs[ind],outsampresults$predmean[ind]) 
          rmsepop[k]=sqrt(mean((outsampresults$obs[ind]-outsampresults$predmean[ind])^2,na.rm=T))
        }
        out$outsampfitstatspop=list(R2pop=R2pop,rmsepop=rmsepop)
      }
    }
  }
  class(out)="GPpred"
  return(out)
}
