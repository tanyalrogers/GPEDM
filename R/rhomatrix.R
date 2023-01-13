#' Get pairwise dynamic correlation matrix
#'
#' Fits a GP model to each pair of populations and constructs a matrix of pairwise
#' dynamic correlation (\code{rho}) values. The matrix can be passed to \code{fitGP}
#' under argument \code{rhomatrix}.
#' 
#' @details 
#' The \code{pop} argument is *required* and there must be more than 1 population.  
#' Function scales data prior to subsetting to each pair of populations.
#' 
#' @inheritParams fitGP
#'
#' @return A square matrix of pairwise dynamic correlations between all populations. Rows 
#' and columns are named with population names.
#' @seealso \code{\link{fitGP}}
#' @references Rogers, T. L. and Munch, S. B. 2020. Hidden similarities in the dynamics 
#'   of weakly synchronous marine metapopulation. Proceedings of the National Academy of 
#'   Sciences 117(1):479-485
#' @export
#' @keywords functions

getrhomatrix=function(data=NULL,y,x=NULL,pop,time=NULL,E=NULL,tau=NULL,
               scaling=c("global","local","none"),
               initpars=NULL,modeprior=1,augdata=NULL) {
  
  #this is duplicated from fitGP
  scaling <- match.arg(scaling)
  
  inputs=list()
  
  #if x is matrix, store names of predictors, if available
  if(!is.null(colnames(x))) {
    inputs$x_names=colnames(x)
  }
  
  #if data frame is supplied, take columns from it and store names
  if(!is.null(data)) {
    inputs$y_names=y
    y=data[,y]
    if(!is.null(x)) { 
      inputs$x_names=x
      x=data[,x] 
    }
    if(!is.null(pop)) { 
      inputs$pop_names=pop
      pop=data[,pop]
    }
    if(!is.null(time)) {
      inputs$time_names=time
      time=data[,time] 
    }
  }
  
  #if pop is not supplied, create vector of 1s (all same population)
  if(is.null(pop)) {
    pop=rep(1,length(y))
  }   
  #if time is not supplied, make a numeric index (only used in output)
  if(is.null(time)) { 
    up=unique(pop)
    time=y*0
    for(i in 1:length(up)) {
      time[pop==up[i]]=1:length(time[pop==up[i]])
    }
  }
  
  #if E and tau are supplied, generate lags of each x
  if(!is.null(E)) {
    if(is.null(x)) { #if x is not supplied, use y as x
      x=y
      if(!is.null(inputs$y_names)) inputs$x_names=inputs$y_names
    }
    x=makelags(y=x,pop=pop,E=E,tau=tau,yname=inputs$x_names)
    if(!is.null(inputs$x_names)) {
      inputs$x_names2=colnames(x)
      if(any(grepl("_1",inputs$x_names) | grepl("_2",inputs$x_names))) {
        message("It looks like x might already contain lags, are you sure you want to be using E and tau?")
      }
    }
    inputs$E=E
    inputs$tau=tau
  }
  
  #make sure x is a matrix, not vector or data frame
  x=as.matrix(x)
  d=ncol(x) #embedding dimension, or number of predictors
  
  #augmentation data
  if(!is.null(augdata)) {
    xaug=as.matrix(augdata[,inputs$x_names])
    yaug=augdata[,inputs$y_names]
    timeaug=augdata[,inputs$time_names]
    if(is.null(inputs$pop_names)) {
      popaug=rep(1,nrow(augdata))
    } else {
      popaug=augdata[,inputs$pop_names]
    }
    yd=c(y,yaug)
    xd=rbind(x,xaug)
    timed=c(time,timeaug)
    popd=c(pop,popaug)
    primary=c(rep(T,length(y)),rep(F,length(yaug)))
  } else { 
    #no augmentation data
    yd=y
    xd=x
    timed=time
    popd=pop
    primary=rep(T,length(y))
  }
  
  #rescale data
  if(scaling=="none") {
    ymeans=mean(yd,na.rm=T)
    ysds=sd(yd,na.rm=T)
    yds=yd
    xmeans=apply(xd,2,mean,na.rm=T)
    xsds=apply(xd,2,sd,na.rm=T)
    xds=xd
    #issue warnings if data are not scaled properly
    if(ymeans>0.1|ymeans<(-0.1)|ysds>1.1|ysds<(-0.9)) {
      warning("y is not scaled properly, model may be unreliable ", 
              "mean=",round(ymeans,3)," sd=",round(ysds,3),call. = F,immediate. = T)
    }
    if(any(xmeans>0.1)|any(xmeans<(-0.1))|any(xsds>1.1)|any(xsds<(-0.9))) {
      warning("one or more x is not scaled properly, model may be unreliable. ", 
              "means=",paste(round(xmeans,3),collapse=" "),
              " sds=",paste(round(xsds,3),collapse=" "),call. = F,immediate. = T)
    }
  }
  if(scaling=="global") {
    ymeans=mean(yd,na.rm=T)
    ysds=sd(yd,na.rm=T)
    yds=scale(yd)
    xmeans=apply(xd,2,mean,na.rm=T)
    xsds=apply(xd,2,sd,na.rm=T)
    xds=apply(xd,2,scale)
  }
  if(scaling=="local") {
    ymeans=tapply(yd,pop,mean,na.rm=T)
    ysds=tapply(yd,pop,sd,na.rm=T)
    xlist=split(as.data.frame(xd),pop)
    xmeans=lapply(xlist,function(x) apply(x,2,mean,na.rm=T))
    xsds=lapply(xlist,function(x) apply(x,2,sd,na.rm=T))
    up=unique(pop)
    xds=xd
    yds=yd
    for(i in 1:length(up)) {
      locmean=ymeans[as.character(up[i])==names(ymeans)]
      locsd=ysds[as.character(up[i])==names(ysds)]
      yds[pop==up[i]]=(yds[pop==up[i]]-locmean)/locsd
      for(j in 1:d) {
        locmean=xmeans[[which(as.character(up[i])==names(xmeans))]][j]
        locsd=xsds[[which(as.character(up[i])==names(xsds))]][j]
        xds[pop==up[i],j]=(xd[pop==up[i],j]-locmean)/locsd          
      }
    }
  }
  
  #end duplicated section
  
  up=unique(pop)
  np=length(up)
  if(np==1) {
    stop("There must be more than 1 population to get a pairwise matrix.")
  }
  rhomat=matrix(0, nrow = np, ncol = np, dimnames = list(up,up)) 
  
  for(i in 1:(np-1)) {
    for(j in (i+1):np) {
      popind=which(popd==up[i] | popd==up[j])
      ydsij=yds[popind]
      xdsij=xds[popind,]
      popdij=popd[popind]
      timedij=timed[popind]
      fittemp=suppressWarnings( #to not get scaling warnings
        fitGP(y=ydsij,x=xdsij,pop=popdij,time=timedij,
              scaling="none",initpars=initpars,modeprior=modeprior)
      )
      rhomat[i,j]=fittemp$pars["rho"]
    }
  }
  
  rhomatrix=rhomat+t(rhomat)+diag(np)
  return(rhomatrix)
}

#' Make a matrix postive definite
#'
#' If the output matrix from \code{\link{getrhomatrix}} is not positive
#' definite, this will convert it to a matrix that is.
#' 
#' @param mat A symmetrical square matrix.
#'
#' @return A version of `mat` that is positive definite.
#' @seealso \code{\link{getrhomatrix}}
#' @export
#' @keywords functions

posdef=function(mat) {
  eigd=eigen(t(mat)%*%mat)
  L=eigd$values
  V=eigd$vectors
  matpd=V%*%diag(sqrt(L))%*%t(V)
  colnames(matpd)=colnames(mat)
  rownames(matpd)=rownames(mat)
  return(matpd)
}
