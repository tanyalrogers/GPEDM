#' Fit a "fisheries" GP model
#'
#' Fits an alternate version of the GP model with one additonal parameter (catchability),
#' designed for use in fisheries applications.
#' 
#' This fits a GP model of the form
#' \deqn{y=f(m-bh,z)}
#' where \eqn{y} is catch per unit effort (CPUE), \eqn{m} are lags of CPUE, \eqn{h} are lags of
#' harvest (in numbers or biomass), \eqn{b} is catchability (scalar), and \eqn{z} are
#' optional covariates. CPUE is assumed proportional to total biomass (or numbers),
#' with proportionality constant \eqn{b}. The composite variable \eqn{m-bh} is the biomass or
#' number of individuals remaining after harvesting (escapement).
#' 
#' Parameter \eqn{b} is found using \code{optimize} applied to the posterior likelihood. 
#' Alternatively, a fixed value for \eqn{b} can be provided under \code{bfixed}.
#' 
#' Parameter \code{ytrans} applies a transformation to \eqn{y} before fitting the model. 
#' Inputs \code{y} and \code{m} should be in untransformed CPUE. \code{ytrans="none"} (the
#' default) will apply no tranformation, \code{ytrans="log"} with compute \eqn{log(y)},
#' \code{ytrans="gr1"} will compute \eqn{log(y_t/m_{t-1})}, and \code{ytrans="gr2"} will 
#' compute \eqn{log(y_t/(m_{t-1}-bh_{t-1})}.
#' 
#' Using this method requires the use of \code{data} with pre-generated lags
#' (option A1 in \code{\link{fitGP}}). For more details on fitting a fisheries
#' model and an example see the vignette. For more elaboration on the inputs,
#' (e.g. `scaling`, `augdata`) see \code{\link{fitGP}}.
#' 
#' @param data A data frame, required.
#' @param y The response variable (required). Typically CPUE.
#' @param m Lags of the response variable (required). Typically lags of CPUE.
#' @param h Lags of the variable to be multipled by b and subtracted from m (required).
#'   Typically harvest. The dimensions of m and h must match.
#' @param z Other predictor variables, e.g. covariates, that are unmodified (optional).
#' @param pop Identifies separate populations (optional, if not supplied, defaults to 1
#'   population). Population values can be either numeric, character, or factor. 
#' @param time The time index (recommended). If not supplied, defaults to a numeric index.
#' @param xname What the composite variable m-bh should be called. Defaults to "escapement".
#' @param ytrans Tranformation to apply to y before fitting (m remains untransformed). 
#'   Either "none" (default), "log", "gr1", or "gr2". R2 is calculated on y in its original
#'   units.
#' @param bfixed Fixes b and bypasses optimization.
#' 
#' @inheritParams fitGP
#' @inheritParams predict.GP
#' 
#' @return A list (class GP and GPpred) with the same elements as \code{\link{fitGP}} and with
#'   additonal element \code{b}, the names of m, h, z stored under \code{inputs}, and the 
#'   composite variable (escapement) included in the \code{insampresults} table. 
#' @export
#' @keywords functions

fitGP_fish=function(data,y,m,h,z=NULL,pop=NULL,time=NULL,
                scaling=c("global","local","none"),
                initpars=NULL,modeprior=1,fixedpars=NULL,rhofixed=NULL,
                rhomatrix=NULL,augdata=NULL,
                predictmethod=NULL,newdata=NULL,
                xname="escapement",ytrans=c("none","log","gr1","gr2"),bfixed=NULL) {
  
  #check that m and h have the same number of columns
  if(length(m)!=length(h)) {
    stop("m and h must have the same number of columns")
  }
  
  cl <- match.call()
  scaling <- match.arg(scaling)
  ytrans <- match.arg(ytrans)
  
  #extract variables
  md=data[,m,drop=F]
  hd=data[,h,drop=F] 
  
  if(!is.null(bfixed)) {
    b=bfixed
  } else {
    #optimize b
    bmax=min(md/hd, na.rm=T)
    boptim=optimize(GPlike_fish, interval = c(0, bmax), maximum = TRUE,
                    data=data,y=y,m=m,h=h,z=z,pop=pop,time=time,scaling=scaling,
                    initpars=initpars,modeprior=modeprior,fixedpars=fixedpars,rhofixed=rhofixed,
                    rhomatrix=rhomatrix,augdata=augdata,xname=xname,ytrans=ytrans)
    #get value of b
    b=boptim$maximum
  }
  
  #get final model fit
  
  #compute composite variable and append z to create new x
  x2=md-b*hd
  x=paste(xname,1:ncol(x2),sep="_")
  colnames(x2)=x
  data=cbind(data,x2)
  #get transformed y
  yt=paste(y,"trans",sep="_")
  yd=ytransfun(y=data[,y], m1=md[,1], e1=x2[,1], ytrans=ytrans)
  yd=data.frame(yd); colnames(yd)=yt
  data=cbind(data,yd)
  
  #do same for augdata if present
  if(!is.null(augdata)) {
    #extract variables
    maug=augdata[,m,drop=F]
    haug=augdata[,h,drop=F]
    #compute composite variable and append z to create new x
    x2aug=maug-b*haug
    colnames(x2aug)=x
    augdata=cbind(augdata,x2aug)
    #get transformed y
    ydaug=ytransfun(y=augdata[,y], m1=maug[,1], e1=x2aug[,1], ytrans=ytrans)
    ydaug=data.frame(ydaug); colnames(ydaug)=yt
    augdata=cbind(augdata,ydaug)
  }
  
  if(!is.null(z)) { x=c(x,z) }

  output=fitGP(data=data,y=yt,x=x,pop=pop,time=time,scaling=scaling,initpars=initpars,
             modeprior=modeprior,fixedpars=fixedpars,rhofixed=rhofixed,rhomatrix=rhomatrix,augdata=augdata)
  output$b=b
  
  #backtransform predictions
  ypred_trans=output$insampresults$predmean
  ypred=ytransfuninv(y=ypred_trans, m1=md[,1], e1=x2[,1], ytrans=ytrans)
  yobs=data[,y]
  
  #change insampresults 
  colnames(output$insampresults)=c("timestep","pop","predmean_trans","predfsd_trans",
                                   "predsd_trans","obs_trans")
  res2=cbind.data.frame(predmean=ypred, obs=yobs, data[,x,drop=F])
  output$insampresults=cbind.data.frame(output$insampresults,res2)
  
  #recalculate R2 values
  output$insampfitstats["R2"]=getR2(yobs,ypred)
  output$insampfitstats["rmse"]=sqrt(mean((yobs-ypred)^2,na.rm=T))
  pop=output$insampresults$pop
  if(length(unique(pop))>1) { #within site fit stats
    up=unique(pop)
    np=length(up)
    R2pop<-rmsepop<-numeric(np)
    names(R2pop)=up
    names(rmsepop)=up
    for(k in 1:np) {
      ind=which(pop==up[k])
      R2pop[k]=getR2(yobs[ind],ypred[ind])
      rmsepop[k]=sqrt(mean((yobs[ind]-ypred[ind])^2,na.rm=T))
    }
    output$insampfitstatspop=list(R2pop=R2pop,rmsepop=rmsepop)
  }
  
  output$inputs$m=md
  output$inputs$m_names=m
  output$inputs$h_names=h
  output$inputs$z_names=z
  output$inputs$y0_names=y
  output$inputs$ytrans=ytrans
  
  if(!is.null(predictmethod)|!is.null(newdata)) { #generate predictions if requested
    predictresults=predict.GP(output,predictmethod,newdata) 
    output=c(output,predictresults)
  }
  
  output$call=cl
  class(output)=c("GP","GPpred")
  return(output)
}


GPlike_fish=function(b,data,y,m,h,z,pop,time,scaling,initpars,modeprior,fixedpars,rhofixed,rhomatrix,augdata,xname="escapement",ytrans) {
  
  #extract variables
  md=data[,m,drop=F]
  hd=data[,h,drop=F]
  #compute composite variable and append z to create new x
  x2=md-b*hd
  x=paste(xname,1:ncol(x2),sep="_")
  colnames(x2)=x
  data=cbind(data,x2)
  #get transformed y
  yt=paste(y,"trans",sep="_")
  yd=ytransfun(y=data[,y], m1=md[,1], e1=x2[,1], ytrans=ytrans)
  yd=data.frame(yd); colnames(yd)=yt
  data=cbind(data,yd)
  
  #do same for augdata if present
  if(!is.null(augdata)) {
    #extract variables
    maug=augdata[,m,drop=F]
    haug=augdata[,h,drop=F]
    #compute composite variable and append z to create new x
    x2aug=maug-b*haug
    colnames(x2aug)=x
    augdata=cbind(augdata,x2aug)
    #get transformed y
    ydaug=ytransfun(y=augdata[,y], m1=maug[,1], e1=x2aug[,1], ytrans=ytrans)
    ydaug=data.frame(ydaug); colnames(ydaug)=yt
    augdata=cbind(augdata,ydaug)
  }
  
  if(!is.null(z)) { x=c(x,z) } 
  
  #fit model
  mfit=fitGP(data=data,y=yt,x=x,pop=pop,time=time,scaling=scaling,initpars=initpars,
             modeprior=modeprior,fixedpars=fixedpars,rhofixed=rhofixed,rhomatrix=rhomatrix,augdata=augdata)
  
  #return posterior likelihood
  return(mfit$insampfitstats["ln_post"])
}


#should predict_seq update b? This needs modification
#need to write MSY function
#what to do about z

#likelihood profile
# optl=10
# btest=seq(0,bmax,length.out=optl)
# lltest=numeric(optl)
# for(i in 1:optl) {
#   lltest[i]=GPlike_fish(btest[i],data=data,y=y,m=m,h=h,z=z,pop=pop,time=time,scaling=scaling,
#                         initpars=initpars,modeprior=modeprior,rhofixed=rhofixed,
#                         rhomatrix=rhomatrix,augdata=augdata,xname=xname)
# }
# plot(btest, lltest, type="l")

#transform
ytransfun=function(y, m1, e1, ytrans) {
  if(ytrans=="none") return(y)
  if(ytrans=="log") return(log(y))
  if(ytrans=="gr1") return(log(y/m1))
  if(ytrans=="gr2") return(log(y/e1))
}

#inverse transform
ytransfuninv=function(y, m1, e1, ytrans) {
  if(ytrans=="none") return(y)
  if(ytrans=="log") return(exp(y))
  if(ytrans=="gr1") return(exp(y)*m1)
  if(ytrans=="gr2") return(exp(y)*e1)
}