#' Print a summary of a GP model
#'
#' Prints out some key information including hyperparameters and fit stats.
#'
#' @param object Output from \code{fitGP}.
#' @method summary GP
#' @export
#' 
summary.GP=function(object) {
  d=ncol(object$inputs$X)
  cat("Number of predictors:",d,"\n")
  cat("Length scale parameters:\n")
  if(!is.null(object$inputs$x_names2)) {
    print(data.frame(predictor=object$inputs$x_names2,posteriormode=round(object$pars[1:d],5)))
  } else if(!is.null(object$inputs$x_names)) {
    print(data.frame(predictor=object$inputs$x_names,posteriormode=round(object$pars[1:d],5)))
  } else {
    print(data.frame(posteriormode=round(object$pars[1:d],5)))
  }
  cat("Process variance (ve):",object$pars["ve"])
  cat("\nPointwise prior variance (sigma2):",object$pars["sigma2"])
  np=length(unique(object$inputs$pop))
  cat("\nNumber of populations:",np)
  if(np>1) {
    cat("\nDynamic correlation (rho):",object$pars["rho"])
  }
  cat("\nIn-sample R-squared:",object$insampfitstats["R2"],"\n")
  if(np>1) {
    cat("In-sample R-squared by population:\n")
    print(data.frame(R2=object$insampfitstatspop$R2))
  }
  if(!is.null(object$outsampfitstats)) {
    cat("Out-of-sample R-squared:",object$outsampfitstats["R2"])
    if(np>1) {
      cat("\nOut-of-sample R-squared by population:\n")
      print(data.frame(R2=object$outsampfitstatspop$R2))
    }
  }
}

#' Plot predictions from a GP model
#'
#' Plots observed and predicted values from a GP model. If the model includes
#' multiple populations, separate plots are produced for each population. The
#' standard deviations plotted are predfsd.
#'
#' @param x Output from \code{fitGP} or \code{predict.GP}.
#' @param plotinsamp Plot the in-sample results. Defaults to out-of-sample
#'   results if available.
#' @export

plot.GPpred=function(x, plotinsamp=F) {
  
  if(is.null(x$outsampresults) | plotinsamp) {
    dplot=x$insampresults
    message("Plotting in sample results.")
  } else {
    dplot=x$outsampresults
    message("Plotting out of sample results.")
  }
  
  up=unique(dplot$pop)
  np=length(up)  
  
  if(!is.null(x$inputs$y_names)) {
    yl=x$inputs$y_names
  } else {
    yl="y"
  }
  
  # par(mfrow=c(np,ifelse((ncol(x$inputs$x)==1 | is.null(x$inputs)),3,2)),mar=c(5,4,2,2))
  old.par <- par(mfrow=c(min(4,np),2),mar=c(5,4,2,2))
  on.exit(par(old.par),add = T,after = F)

  for(i in 1:np) {
    dploti=subset(dplot,pop==up[i])
    plot(dploti$timestep,dploti$predmean,type="o",ylab=yl,xlab="Time",
         ylim=range(dploti$predmean+dploti$predfsd,dploti$predmean-dploti$predfsd,dploti$obs,na.rm=T),main=up[i])
    lines(dploti$timestep,dploti$predmean+dploti$predfsd, lty=2)
    lines(dploti$timestep,dploti$predmean-dploti$predfsd, lty=2)
    points(dploti$timestep,dploti$obs,col="blue")
    legend(x = "bottomright",legend = c("obs","pred"),col=c("blue","black"),pch=1)
    
    if(!is.null(dploti$obs)) {    
      plot(dploti$obs,dploti$predmean,ylab=paste(yl,"pred"),xlab=paste(yl,"obs"),main=up[i])
      abline(a=0,b=1)
    }    
  }
  
}

#' Get conditional responses from a GP model
#'
#' This will obtain responses to each predictor variable across its range of
#' observed values, conditional on all other predictor variables being set to 0
#' (their mean value, if data are scaled properly). If automatic scaling is
#' used, predictors and responses are backtransformed to their original scale.
#' If the model includes multiple populations, responses are obtained for each
#' population.
#'
#' @param fit Output from \code{fitGP}.
#' @param xrange The range of *scaled* predictor values over which responses
#'   should be evaluated. Options are \code{"local"} (the range observed within
#'   each population), or \code{"global"} (the range across all populations).
#'   Irrelevant if there is only one population. The default (\code{"default"})
#'   will be \code{"local"} if \code{scaling="local"} was used, and
#'   \code{"global"} if \code{scaling="global"} or \code{scaling="none"} was
#'   used. This option can useful for either preventing or allowing for
#'   extrapolation beyond the range of the data in a given population.
#' @param extrap Percentage to extrapolate beyond the predictor range for each
#'   population (i.e. predictor values will range from \code{(1-extrap)*min(X)}
#'   to \code{(1+extrap)*max(X)} Defaults to 0.01.
#' @param nvals The number of values to evaluate across the predictor range for
#'   each population. Defaults to 25, which should be fine for most
#'   applications.
#' @param plot Produce a plot, or not (logical, defaults to TRUE).
#'
#' @return Returns (invisibly) a data frame containing the predictor values and conditional responses 
#'   (means and standard deviations).
#' @export
getconditionals=function(fit,xrange="default", extrap=0.01, nvals=25, plot=T) {

  #need to add 2d option*****
  
  #extract relevant stuff from model
  iKVs=fit$covm$iKVs
  phi=fit$pars[grepl("phi",names(fit$pars))]
  sigma2=fit$pars[names(fit$pars)=="sigma2"]
  rho=fit$pars[names(fit$pars)=="rho"]
  X=fit$inputs$X
  Y=fit$inputs$Y
  Pop=fit$inputs$Pop
  scaling=fit$scaling$scaling
  ymeans=fit$scaling$ymeans
  ysds=fit$scaling$ysds
  xmeans=fit$scaling$xmeans
  xsds=fit$scaling$xsds
  
  if(xrange=="default") {
    if(scaling=="local") { xrange2="local" } 
    else { xrange2="global" }
  } else { xrange2=xrange }
  
  up=unique(Pop)
  np=length(up)
  d=ncol(X)
  Tslp=nvals
  
  outlist=NULL
  for(k in 1:np) { #populations
    indi=which(Pop==up[k])
    xval<-predmean<-predsd<-xg<-matrix(0,nrow=Tslp,ncol=d)
    poppred=rep(up[k],Tslp)
    for(i in 1:d) { #predictors
      if(xrange2=="local") {
        xgi=seq((1-extrap)*min(X[indi,i]),(1+extrap)*max(X[indi,i]),length.out = Tslp)
      } else {
        xgi=seq((1-extrap)*min(X[,i]),(1+extrap)*max(X[,i]),length.out = Tslp)
      }
      xp=xg
      xp[,i]=xgi
      covmnew=getcov(phi,sigma2,rho,X,xp,Pop,poppred)
      Cs=covmnew$Cd #covariance matrix
      predmean[,i]=Cs%*%(iKVs%*%Y)
      predvar=numeric(length = Tslp)
      for(j in 1:Tslp) {
        predvar[j]=sigma2-Cs[j,]%*%iKVs%*%Cs[j,]
      }
      predsd[,i]=sqrt(predvar)
      xval[,i]=xgi
    }
    outlist[[k]]=list(poppred=poppred,xval=xval,predmean=predmean,predsd=predsd)
  }
  
  #unscale predictions and x
  if(scaling=="global") {
    for(i in 1:np) {
      for(j in 1:d) {
        outlist[[i]]$predmean[,j]=outlist[[i]]$predmean[,j]*ysds+ymeans
        outlist[[i]]$predsd[,j]=sqrt(outlist[[i]]$predsd[,j]^2*ysds^2)
        outlist[[i]]$xval[,j]=outlist[[i]]$xval[,j]*xsds[j]+xmeans[j]
      }
    }
  }
  if(scaling=="local") {
    for(i in 1:np) {
      for(j in 1:d) {
        locmean=ymeans[as.character(up[i])==names(ymeans)]
        locsd=ysds[as.character(up[i])==names(ysds)]
        outlist[[i]]$predmean[,j]=outlist[[i]]$predmean[,j]*locsd+locmean
        outlist[[i]]$predsd[,j]=sqrt(outlist[[i]]$predsd[,j]^2*locsd^2)
        locmean=xmeans[[which(as.character(up[i])==names(xmeans))]][j]
        locsd=xsds[[which(as.character(up[i])==names(xsds))]][j]
        outlist[[i]]$xval[,j]=outlist[[i]]$xval[,j]*locsd+locmean
      }
    }
  }

  if(!is.null(fit$inputs$x_names2)) { xlabels=fit$inputs$x_names2 } 
  else if(!is.null(fit$inputs$x_names)) { xlabels=fit$inputs$x_names }
  else { xlabels= paste0("x",1:d) }
  
  #combine into dataframe  #this needs work, put in longer format*****
  out=lapply(outlist,function(x) {
    cc2=cbind.data.frame(x$poppred,x$xval,x$predmean,x$predsd)
    colnames(cc2)=c("pop",xlabels,paste0(xlabels,"_yMean"),paste0(xlabels,"_ySD"))
    return(cc2)
  })
  out=do.call("rbind",out)
  
  if(plot) {
    if(!is.null(fit$inputs$y_names)) {
      yl=fit$inputs$y_names
    } else {
      yl="y"
    }
    old.par <- par(mfrow=c(min(4,np),min(4,d)),mar=c(5,4,2,2))
    on.exit(par(old.par),add = T,after = F)

    for(i in 1:np) {
      ylims=range(out[out$pop==up[i],grep("_yMean",colnames(out))]+out[out$pop==up[i],grep("_ySD",colnames(out))],
                  out[out$pop==up[i],grep("_yMean",colnames(out))]-out[out$pop==up[i],grep("_ySD",colnames(out))])
      pdata=outlist[[i]]
      for(j in 1:d) {
        plot(pdata$xval[,j],pdata$predmean[,j], type="l",xlab=xlabels[j],ylab=yl,main=up[i],ylim=ylims)
             #ylim=range(pdata$predmean[,j]+pdata$predsd[,j],pdata$predmean[,j]-pdata$predsd[,j]))
        lines(pdata$xval[,j],pdata$predmean[,j]+pdata$predsd[,j],lty=2)
        lines(pdata$xval[,j],pdata$predmean[,j]-pdata$predsd[,j],lty=2)
      }
    }
  }
  return(invisible(out))
}
