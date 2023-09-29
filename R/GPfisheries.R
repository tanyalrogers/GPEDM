#' Fit a "fisheries" GP model
#'
#' Fits an alternate version of the GP model with one additonal parameter (catchability),
#' designed for use in fisheries applications.
#' 
#' This fits a GP model of the form
#' y=f(m-bh,z)
#' 
#' Where y is CPUE, m are lags of CPUE, h are lags of harvest, b is catchability (scalar),
#' and z are optional covariates. The dimensions of m and h must match.
#' 
#' Parameter b is found using \code{optimize} applied to the posterior likelihood.
#' 
#' Using this method requires the use of \code{data} with pre-generated lags
#' (option A1 in \code{\link{fitGP}}).
#' 
#' @param data A data frame, required.
#' @param y The response variable (required).
#' @param m Lags of the response variable (required).
#' @param h Lags of the variable to be multipled by b and subtracted from m (required).
#' @param z Other predictor variables, e.g. covariates, that are unmodified (optional).
#' @param pop Identifies separate populations (optional, if not supplied, defaults to 1
#'   population). Population values can be either numeric, character, or factor. 
#' @param time The time index (recommended). If not supplied, defaults to a numeric index.
#' @param xname What the composite variable m-bh should be called. Defaults to "escapement".
#' 
#' @inheritParams fitGP
#' @inheritParams predict.GP
#' 
#' @return A list (class GP and GPpred) with the same elements as fitGP with additonal
#'   element \code{b}, data_x, and names of m, h, z stored under inputs.
#' @export
#' @keywords functions

fitGP_fish=function(data,y,m,h,z=NULL,pop=NULL,time=NULL,
                scaling=c("global","local","none"),
                initpars=NULL,modeprior=1,rhofixed=NULL,
                rhomatrix=NULL,#augdata=NULL,
                predictmethod=NULL,newdata=NULL,
                xname="escapement") {
  
  #check that m and h have the same number of columns
  if(length(m)!=length(h)) {
    stop("m and h must have the same number of columns")
  }
  
  #cl <- match.call()
  
  
  GPlike_fish=function(b,data,y,m,h,z,pop,time,scaling,initpars,modeprior,rhofixed,rhomatrix,xname="escapement") {
    
    #extract variables
    md=data[,m,drop=F]
    hd=data[,h,drop=F]

    #compute composite variable and append z to create new x
    x2=md-b*hd
    x=paste(xname,1:ncol(x2),sep="_")
    colnames(x2)=x
    data=cbind(data,x2)
    
    if(!is.null(z)) { x=c(x,z) } 
    
    #fit model
    mfit=fitGP(data=data,y=y,x=x,pop=pop,time=time,scaling=scaling,initpars=initpars,
               modeprior=modeprior,rhofixed=rhofixed,rhomatrix=rhomatrix)
    
    #return posterior likelihood
    return(mfit$insampfitstats["ln_post"])
  }
  
  #optimize b
  bmax=min(m/h) #?
  boptim=optimize(GPlike_fish, interval = c(0, bmax), maximum = TRUE)
  
  #get value of b
  b=boptim$value
  
  #get final model fit
  
  #extract variables
  md=data[,m,drop=F]
  hd=data[,h,drop=F]
  
  #compute composite variable and append z to create new x
  x2=md-b*hd
  x=paste(xname,1:ncol(x2),sep="_")
  colnames(x2)=x
  data=cbind(data,x2)
  
  if(!is.null(z)) { x=c(x,z) }

  output=fitGP(data=data,y=y,x=x,pop=pop,time=time,scaling=scaling,initpars=initpars,
             modeprior=modeprior,rhofixed=rhofixed,rhomatrix=rhomatrix)
  output$b=b
  output$data_x=data[,x]
  object$inputs$m_names=m
  object$inputs$h_names=h
  object$inputs$z_names=z
  
  if(!is.null(predictmethod)|!is.null(newdata)) { #generate predictions if requested
    predictresults=predict.GP(output,predictmethod,newdata) 
    output=c(output,predictresults)
  }
  
  #output$call=cl
  return(output)
}


#modify predict functions so that if model list contains b, compute composite predictors in newdata
#should work for loo and newdata
#getconditionals should work fine

#should predict_seq update b? This needs modification
#need to write MSY function
#what to do about z

