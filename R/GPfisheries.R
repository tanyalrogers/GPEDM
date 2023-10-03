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
#' Using this method requires the use of \code{data} with pre-generated lags
#' (option A1 in \code{\link{fitGP}}).
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
#' @param bfixed Fixes b and bypasses optimization.
#' 
#' @inheritParams fitGP
#' @inheritParams predict.GP
#' 
#' @return A list (class GP and GPpred) with the same elements as \code{\link{fitGP}} with
#'   additonal element \code{b}, the names of m, h, z stored under \code{inputs}, and the 
#'   composite variable (escapement) included in the \code{insampresults} table. 
#' @export
#' @keywords functions

fitGP_fish=function(data,y,m,h,z=NULL,pop=NULL,time=NULL,
                scaling=c("global","local","none"),
                initpars=NULL,modeprior=1,rhofixed=NULL,
                rhomatrix=NULL,augdata=NULL,
                predictmethod=NULL,newdata=NULL,
                xname="escapement",bfixed=NULL) {
  
  #check that m and h have the same number of columns
  if(length(m)!=length(h)) {
    stop("m and h must have the same number of columns")
  }
  
  cl <- match.call()
  scaling <- match.arg(scaling)
  
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
                    initpars=initpars,modeprior=modeprior,rhofixed=rhofixed,
                    rhomatrix=rhomatrix,augdata=augdata,xname=xname)
    #get value of b
    b=boptim$maximum
  }
  
  #get final model fit
  
  #compute composite variable and append z to create new x
  x2=md-b*hd
  x=paste(xname,1:ncol(x2),sep="_")
  colnames(x2)=x
  data=cbind(data,x2)
  
  #do same for augdata if present
  if(!is.null(augdata)) {
    #extract variables
    maug=augdata[,m,drop=F]
    haug=augdata[,h,drop=F]
    #compute composite variable and append z to create new x
    x2aug=maug-b*haug
    colnames(x2aug)=x
    augdata=cbind(augdata,x2aug)
  }
  
  if(!is.null(z)) { x=c(x,z) }

  output=fitGP(data=data,y=y,x=x,pop=pop,time=time,scaling=scaling,initpars=initpars,
             modeprior=modeprior,rhofixed=rhofixed,rhomatrix=rhomatrix,augdata=augdata)
  output$b=b
  output$insampresults=cbind(output$insampresults,data[,x,drop=F])
    output$inputs$m_names=m
  output$inputs$h_names=h
  output$inputs$z_names=z
  
  if(!is.null(predictmethod)|!is.null(newdata)) { #generate predictions if requested
    predictresults=predict.GP(output,predictmethod,newdata) 
    output=c(output,predictresults)
  }
  
  output$call=cl
  return(output)
}


GPlike_fish=function(b,data,y,m,h,z,pop,time,scaling,initpars,modeprior,rhofixed,rhomatrix,augdata,xname="escapement") {
  
  #extract variables
  md=data[,m,drop=F]
  hd=data[,h,drop=F]
  #compute composite variable and append z to create new x
  x2=md-b*hd
  x=paste(xname,1:ncol(x2),sep="_")
  colnames(x2)=x
  data=cbind(data,x2)
  
  #do same for augdata if present
  if(!is.null(augdata)) {
    #extract variables
    maug=augdata[,m,drop=F]
    haug=augdata[,h,drop=F]
    #compute composite variable and append z to create new x
    x2aug=maug-b*haug
    colnames(x2aug)=x
    augdata=cbind(augdata,x2aug)
  }
  
  if(!is.null(z)) { x=c(x,z) } 
  
  #fit model
  mfit=fitGP(data=data,y=y,x=x,pop=pop,time=time,scaling=scaling,initpars=initpars,
             modeprior=modeprior,rhofixed=rhofixed,rhomatrix=rhomatrix,augdata=augdata)
  
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
