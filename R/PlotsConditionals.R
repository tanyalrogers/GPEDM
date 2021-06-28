#' Plot predictions from a GP model.
#'
#' Plots observed and predicted values from a GP model. If the model includes
#' multiple populations, separate plots are produced for each population.
#'
#' @param x Output from \code{fitGP} or \code{predict.GP}.
#' @param plotinsamp Plot the in-sample results. Defaults to out-of-sample
#'   results if available.
#' @export

plot.GP=function(x, plotinsamp=F) {
  old.par <- par(no.readonly = TRUE)
  
  if(is.null(x$outsampresults) | plotinsamp) {
    dplot=x$insampresults
    message("Plotting in sample results.")
  } else {
    dplot=x$outsampresults
    message("Plotting out of sample results.")
  }
  
  up=unique(dplot$pop)
  np=length(up)  
  
  # par(mfrow=c(np,ifelse((ncol(x$inputs$xd)==1 | is.null(x$inputs)),3,2)),mar=c(5,4,2,2))
  par(mfrow=c(np,2),mar=c(5,4,2,2))
  
  for(i in 1:np) {
    dploti=subset(dplot,pop==up[i])
    plot(dploti$timestep,dploti$predmean,type="o",ylab="y",xlab="time",
         ylim=range(dploti$predmean+dploti$predfsd,dploti$predmean-dploti$predfsd,dploti$obs,na.rm=T),main=up[i])
    lines(dploti$timestep,dploti$predmean+dploti$predfsd, lty=2)
    lines(dploti$timestep,dploti$predmean-dploti$predfsd, lty=2)
    points(dploti$timestep,dploti$obs,col="blue")
    legend(x = "bottomright",legend = c("obs","pred"),col=c("blue","black"),pch=1)
    
    plot(dploti$obs,dploti$predmean,ylab="y pred",xlab="y obs",main=up[i])
    abline(a=0,b=1)
    
  }
  
  on.exit(par(old.par))
  
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
#' @param plot Produce a plot, or not (logical).
#'
#' @return A data frame.
#' @export
getconditionals=function(fit, plot=T) {
  #produce sections across each input (all others inputs fixed to 0), for each population
  
  #need to add 2d option*****
  
  #extract relevant stuff from model
  iKVs=fit$covm$iKVs
  phi=fit$pars[grepl("phi",names(fit$pars))]
  tau=fit$pars[names(fit$pars)=="tau"]
  rho=fit$pars[names(fit$pars)=="rho"]
  X=fit$inputs$X
  Y=fit$inputs$Y
  Pop=fit$inputs$Pop
  scaling=fit$scaling$scaling
  ymeans=fit$scaling$ymeans
  ysds=fit$scaling$ysds
  xmeans=fit$scaling$xmeans
  xsds=fit$scaling$xsds
  
  up=unique(Pop)
  np=length(up)
  d=ncol(X)
  Tslp=25
  
  outlist=NULL
  for(k in 1:np) { #populations
    indi=which(Pop==up[k])
    xval<-predmean<-predsd<-xg<-matrix(0,nrow=Tslp,ncol=d)
    poppred=rep(up[k],Tslp)
    for(i in 1:d) { #predictors
      xgi=seq(0.95*min(X[,i]),1.05*max(X[,i]),length.out = Tslp)
      xdp=xg
      xdp[,i]=xgi
      covmnew=getcov(phi,tau,rho,X,xdp,Pop,poppred)
      Cs=covmnew$Cd #covariance matrix
      predmean[,i]=Cs%*%(iKVs%*%Y)
      predvar=numeric(length = Tslp)
      for(j in 1:Tslp) {
        predvar[j]=tau-Cs[j,]%*%iKVs%*%Cs[j,]
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
  
  if(!is.null(fit$inputs$xd_names)) { xlabels=fit$inputs$xd_names }
  else { xlabels= paste("x",1:d) }
  
  #combine into dataframe  #this needs work, put in longer format*****
  out=lapply(outlist,function(x) {
    cc2=cbind.data.frame(x$poppred,x$xval,x$predmean,x$predsd)
    colnames(cc2)=c("pop",xlabels,paste0(xlabels,"_yMean"),paste0(xlabels,"_ySD"))
    return(cc2)
  })
  out=do.call("rbind",out)
  
  if(plot) {
    old.par <- par(no.readonly = TRUE)
    par(mfrow=c(np,d),mar=c(5,4,2,2))
    for(i in 1:np) {
      pdata=outlist[[i]]
      for(j in 1:d) {
        plot(pdata$xval[,j],pdata$predmean[,j], type="l",xlab=xlabels[j],ylab="y",main=up[i],
             ylim=range(pdata$predmean[,j]+pdata$predsd[,j],pdata$predmean[,j]-pdata$predsd[,j]))
        lines(pdata$xval[,j],pdata$predmean[,j]+pdata$predsd[,j],lty=2)
        lines(pdata$xval[,j],pdata$predmean[,j]-pdata$predsd[,j],lty=2)
      }
    }
    on.exit(par(old.par))
  }
  return(out)
}
