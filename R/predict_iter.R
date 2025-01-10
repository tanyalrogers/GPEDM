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
#' @param hrate For "fisheries models" a value from 0 to 1 indicating the harvest rate.
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
predict_iter=function(object,newdata,xlags=NULL,hrate=NULL) { 
  
  if(is.null(object$inputs$time_names)) {
    stop("Model must include `data` and `time` to use this function.")
  }
  if(!is.null(object$inputs$E)) {
    stop("Lags must be pre-generated to use this function. See option A1 in help(fitGP).")
  }
  
  if(!is.null(object$b)) {  
    if(!is.null(xlags)) {
      message("xlags in fisheries models set to `m` variables, input ignored")
    }
    xlags=object$inputs$m_names
    hlags=object$inputs$h_names
    b=object$b
  }
  
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
  
  if(is.null(names(b))) {
    names(b)=up
  }
  
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
        if(length(xlags)>1) {
          newdata[newdata[,timename]==newtimes[i+1] & newdata[,popname]==up[j],xlags[2:length(xlags)]]=
            newdata[newdata[,timename]==newtimes[i] & newdata[,popname]==up[j],xlags[1:(length(xlags)-1)]]
        }
        #compute next catch
        if(!is.null(hrate)) {
          newdata[newdata[,timename]==newtimes[i+1] & newdata[,popname]==up[j],hlags[1]]=predsj/b[as.character(up[j])]*hrate
          if(length(xlags)>1) {
            newdata[newdata[,timename]==newtimes[i+1] & newdata[,popname]==up[j],hlags[2:length(xlags)]]=
              newdata[newdata[,timename]==newtimes[i] & newdata[,popname]==up[j],hlags[1:(length(xlags)-1)]]
          }
        }
      }
    }
  }
  #compile predictions
  outsampresults=do.call(rbind, pred)
  out=list(outsampresults=outsampresults)
  if(!is.null(hrate)) {
    order=as.numeric(rownames(outsampresults))
    out$outsampresults=cbind(out$outsampresults,newdata[order,hlags,drop=F])
  }
  
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

#' Make lag matrix for iterated prediction
#'
#' Creates a lag matrix that can be used for making iterated predictions. 
#' Runs \code{makelags} with specified parameters and \code{forecast=TRUE},
#' then appends \code{nfore-1} empty rows for each pop with timesteps iterated by 1.
#' 
#' @inheritParams makelags
#' @param nfore Number of steps to forecast
#' @param ... Additional arguments passed to makelags
#' @return A data frame.
#' @export
#' @keywords functions
makelags_iter=function(data, nfore, y, pop=NULL, E, tau, time, ...) {
  skipcols=ifelse(is.null(pop), 1, 2)
  #use the forecast feature to create the first row
  foremat1=makelags(data, y=y, pop=pop, time=time, tau=tau, E=E, forecast=T, ...)
  #get remaining timepoints
  if(is.null(pop)) {
    foremat2=data.frame(time=(foremat1[1,time]+1):(foremat1[1,time]+nfore-1)) 
    colnames(foremat2)=time
  } else {
    foremat2=expand.grid(time=(foremat1[1,time]+1):(foremat1[1,time]+nfore-1), pop=unique(data[,pop]))
    colnames(foremat2)=c(time,pop)
  }
  #create empty matrix for future values
  forecols=colnames(foremat1)[-(1:skipcols)]
  foremat3=matrix(NA, nrow=nrow(foremat2), ncol=ncol(foremat1)-skipcols,
                  dimnames=list(NULL,forecols))
  #combine everything
  foremat=rbind(foremat1,cbind(foremat2,foremat3))
  return(foremat)
}

#' Wrapper for computing MSY from a "fisheries model"
#'
#' @param model A \code{fitGP_fish} model.
#' @param newdata Data frame for iterated prediction, such as made with
#'   \code{makelags_iter}.
#' @param hratevec Vector of harvest rates
#' @param tsave Number of timepoints to average over. Should not exceed nrow(newdata)
#' 
#' @export
#' @keywords functions
msy_wrapper=function(model, newdata, hratevec, tsave) {
  
  up=unique(model$inputs$pop)
  catch_1=model$inputs$h_names[1]
  
  msyfore=list()
  
  catchsavepop=expand.grid(pop=up,time=1:tsave,hrate=hratevec)
  catchsavepop$catch=NA
  catchsavepop$cpue=NA
  
  for(h in seq_along(hratevec)) {
    msyfore[[h]]=predict_iter(model, newdata = newdata, hrate = hratevec[h])$outsampresults
    msyfore[[h]]=cbind.data.frame(hrate=hratevec[h],msyfore[[h]])
    for(p in 1:length(up)){
      msyforei=subset(msyfore[[h]], pop==up[p])
      nfore=nrow(msyforei)
      catchsavepop$catch[catchsavepop$hrate==hratevec[h] & catchsavepop$pop==up[p]]=
        msyforei[(nfore-tsave+1):nfore,catch_1]
      catchsavepop$cpue[catchsavepop$hrate==hratevec[h] & catchsavepop$pop==up[p]]=
        msyforei$predmean[(nfore-tsave+1):nfore]
    }
  }
  catchforeall=do.call(rbind, msyfore)
  
  #sum over pops (if more than one)
  catchsavetotal=aggregate(cbind(catch,cpue)~hrate*time, data=catchsavepop, sum)
  
  #average over the last 10 timepoints (keeping split by pop)
  catchsavepopmean=aggregate(cbind(catch,cpue)~hrate*pop, data=catchsavepop, mean)
  
  #average over last 10 timepoints (sum of pops)
  catchsavetotalmean=aggregate(cbind(catch,cpue)~hrate, data=catchsavetotal, mean)
  
  fmsy=catchsavetotalmean$hrate[which.max(catchsavetotalmean$catch)]
  Bmsy=catchsavetotalmean$cpue[which.max(catchsavetotalmean$catch)]
  
  return(list(catchforeall=catchforeall, catchsavepop=catchsavepop, 
              catchsavetotal=catchsavetotal, catchsavepopmean=catchsavepopmean,
              catchsavetotalmean=catchsavetotalmean, fmsy=fmsy, Bmsy=Bmsy))
}

# Multi-population long form (for reference)
# up=unique(RHfore2pop$Region)
# 
# catchres=expand.grid(pop=up,time=1:tsave,hrate=hratevec)
# catchres$catch=NA
# catchres$cpue=NA
# for(h in seq_along(hratevec)) {
#   msyfore=predict_iter(fishfit2pop2, newdata = RHfore2pop, hrate = hratevec[h])$outsampresults
#   for(p in 1:length(up)){
#     msyforei=subset(msyfore, pop==up[p])
#     catchres$catch[catchres$hrate==hratevec[h] & catchres$pop==up[p]]=
#       msyforei$Catch_1[(nfore-tsave+1):nfore]
#     catchres$cpue[catchres$hrate==hratevec[h] & catchres$pop==up[p]]=
#       msyforei$predmean[(nfore-tsave+1):nfore]
#   }
# }
# 
# #average over the last 10 timepoints (keeping split by pop)
# catchresagg=aggregate(cbind(catch,cpue)~hrate*pop, data=catchres, mean)
# 
# #if you might want to add together the values the pops, assuming they are in the same units
# #sum over pops
# catchrestotal=aggregate(cbind(catch,cpue)~hrate*time, data=catchres, sum)
# #average over last 10 timepoints
# catchresaggtotal=aggregate(cbind(catch,cpue)~hrate, data=catchrestotal, mean)
# 
# (fmsy=catchresaggtotal$hrate[which.max(catchresaggtotal$catch)])
# (Bmsy=catchresaggtotal$cpue[which.max(catchresaggtotal$catch)])
# 
# #Total
# par(mfrow=c(1,2),mar=c(4,4,2,1))
# plot(catch~hrate, data=catchrestotal, col="gray")
# lines(catch~hrate, data=catchresaggtotal, lwd=2, col="red")
# abline(v=fmsy, col="red")
# plot(cpue~hrate, data=catchrestotal, col="gray")
# lines(cpue~hrate, data=catchresaggtotal, lwd=2, col="red")
# abline(v=fmsy, col="red")
# 
# #Split up by pop
# par(mfrow=c(1,2),mar=c(4,4,2,1))
# plot(catch~hrate, data=catchres[catchres$pop==up[1],], col="tomato", ylim=range(catchres$catch))
# points(catch~hrate, data=catchres[catchres$pop==up[2],], col="blue")
# lines(catch~hrate, data=catchresagg[catchresagg$pop==up[1],], lwd=2, col="tomato")
# lines(catch~hrate, data=catchresagg[catchresagg$pop==up[2],], lwd=2, col="blue")
# abline(v=fmsy, col="black")
# legend(x="topleft", legend = up, pch=1, col=c("tomato","blue"))
# plot(cpue~hrate, data=catchres[catchres$pop==up[1],], col="red")
# points(cpue~hrate, data=catchres[catchres$pop==up[2],], col="blue")
# lines(cpue~hrate, data=catchresagg[catchresagg$pop==up[1],], lwd=2, col="red")
# lines(cpue~hrate, data=catchresagg[catchresagg$pop==up[2],], lwd=2, col="blue")
# abline(v=fmsy, col="black")
# legend(x="topleft", legend = up, pch=1, col=c("red","blue"))
