#' Get sequentially updated predictions from a GP model
#'
#' Obtains predictions for \code{newdata} given a GP model, but sequentially
#' adds each new observation to the training data and refits the model. This
#' simulates a real-time forecasting application.
#' 
#' In comparison to \code{\link{predict.GP}}, this prediction method is similar
#' to \code{predictmethod="sequential"} (leave-future-out), but in this case,
#' the future time points are not included in the training data. It is also
#' similar to the train/test split using \code{newdata}, but in this case, the
#' training data and model are updated with each timestep.
#' 
#' Using this method requires the use of \code{data} with pre-generated lags
#' (option A1 in \code{\link{fitGP}}), that a \code{time} column is specified,
#' that \code{newdata} has exactly the same columns as \code{data}, and that
#' \code{newdata} contains observed values. 
#' 
#' The original \code{data} data frame is obtained from the global environment,
#' so results will be not be accurate if it is changed. Also, if for whatever
#' reason you want to use \code{update} on the output of \code{predict_seq}, it
#' will probably not work unless you respecify \code{data} and
#' \code{initpars}. 
#'
#' @param object Output from \code{fitGP}.
#' @param newdata Data frame containing the same columns supplied in the
#'   original model.
#' @param restart.initpars Set initpars to pars of the previous model (F, default), or restart 
#'   with the default initspars each time (T). Starting at the last values can reduce the number 
#'   of iterations required, saving time, but might trap you in a local minimum.
#' @return A list (class GP and GPpred) with the same elements as \code{\link{fitGP}}. 
#'   The model information will for the final model including all of the original training
#'   data plus \code{newdata}. Out \code{outsampresults} and \code{outsampfitstats} will be
#'   for the sequentially updated predictions.
#' @export
#' @keywords functions
predict_seq=function(object,newdata,restart.initpars=F) { 
  
  if(is.null(object$inputs$time_names)) {
    stop("Model must include `data` and `time` to use this function.")
  }
  if(!is.null(object$inputs$E)) {
    stop("Lags must be pre-generated to use this function. See option A1 in help(fitGP).")
  }
  
  #get times
  timename=object$inputs$time_names
  newtimes=unique(newdata[,timename])
  
  #get original data from global env
  ogdata=eval(object$call$data)
  
  modelupdated=object
  dataupdated=ogdata
  pred=list()
  
  for(i in seq_along(newtimes)) {
    newdatai=newdata[newdata[,timename]==newtimes[i],]
    #get prediction
    pred[[i]]=predict(modelupdated, newdata=newdatai)$outsampresults
    #append new data
    dataupdated=rbind(dataupdated, newdatai)
    #set initpars
    if(restart.initpars) {
      inits=NULL
    } else {
      inits=modelupdated$pars
    }
    #refit model
    modelupdated=update(modelupdated, data=dataupdated, initpars=inits)
  }
  #compile predictions
  outsampresults=do.call(rbind, pred)
  modelupdated$outsampresults=outsampresults

  #return updated model and predictions
  modelupdated$outsampfitstats=c(R2=getR2(outsampresults$obs,outsampresults$predmean), 
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
    modelupdated$outsampfitstatspop=list(R2pop=R2pop,rmsepop=rmsepop)
  }

  return(modelupdated)
}