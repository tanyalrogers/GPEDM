#' Fit an S-map model
#'
#' Fits a S-map (local linear) model with local weighting parameter theta found using
#' gradient descent. Includes implementation of the nonstationary S-map. 
#' 
#' There are several ways to specify the training data:
#' \itemize{
#'   \item A. supply data frame \code{data}, plus column names or indices for \code{y}, \code{x}, \code{pop}, and \code{time}.
#'   \item B. supply vectors for \code{y}, \code{pop}, and \code{time}, and a matrix or vector for \code{x}.
#' }
#' For each of the above 2 options, there are 3 options for specifying the predictors.
#' \enumerate{
#'   \item supplying \code{y} and \code{x} (omitting \code{E} and \code{tau}) will use the columns of \code{x} as predictors. This allows 
#'   for the most customization.
#'   \item supplying \code{y}, \code{E}, and \code{tau} (omitting \code{x}) will use \code{E} lags of \code{y} (with spacing \code{tau}) as predictors.
#'   This is equivalent to option 3 with \code{x}=\code{y}.
#'   \item supplying \code{y}, \code{x}, \code{E}, and \code{tau} will use \code{E} lags of *each column* of \code{x} (with spacing \code{tau}) as predictors.
#'   Do not use this option if \code{x} already contains lags, in that case use option 1.
#' }
#' Arguments \code{pop} and \code{time} are optional in all of the above cases. 
#' If \code{pop} and \code{time} are supplied, they should NOT contain missing values.
#' This function will also, optionally, produce predictions.
#' See Details for more information.
#' 
#' @details 
#' 
#' The data are assumed to be in time order. This is particularly important if \code{E} and \code{tau}
#' are supplied or \code{predictmethod="sequential"} is used.
#' 
#' The S-map model does not include a hierarchical stucture like the GP model. If multiple
#' populations are supplied, their dynamics are assumed to be identical. The \code{pop} 
#' argument is only relevant if you are supplying \code{E} and \code{tau} to generate
#' lags internally (to ensure there isn't crossover among populations), or if you want to
#' keep track (in the results tables) of which predictions are from which population. 
#' The model fitting code itself does not use the population designations.
#' 
#' \strong{Missing values:} 
#' 
#' Missing values for \code{y} and \code{x} are allowed. Any rows containing
#' missing values will be excluded from the model fitting procedure (reinserted
#' as NAs in the output). If \code{pop} and \code{time} are supplied, they
#' should NOT contain missing values. 
#' 
#' \strong{Hyperparameters:} 
#' 
#' Hyperparameter \code{theta} is the s-map local weighting parameter. A value of 0
#' corresponds to a globally linear model. The value can be fixed using \code{thetafixed},
#' otherwise it is estimated using gradient descent. 
#' 
#' The nonstationary implementation of s-map (ns-map) includes an additional
#' hyperparameter \code{delta}, which is the local weighting in time (as opposed
#' to state space, which is controlled by \code{theta}). In the standard s-map
#' (default \code{ns=FALSE}), \code{delta} is fixed to 0. If fitting a
#' nonstationary model (\code{ns=TRUE}), the value of \code{delta} can be fixed
#' using \code{deltafixed}, otherwise it is estimated using gradient descent,
#' along with \code{theta}. A dynamic linear model (DLM) can be approximated by
#' fixing \code{theta} to 0 and estimating \code{delta}.
#'  
#' \strong{Scaling:} 
#' 
#' Scaling in s-map is not as essential as in the GP model, since \code{theta}
#' is scaled interally to the average distance among points. If you are doing
#' univariate delay embedding, the results will be unaffected by the scaling of
#' the data. However, if you are fitting a multivariate model, you will need to
#' think about data scaling, since the weights are applied equally in all
#' dimensions, and the predictors are no longer in the same units. It is common
#' to standardize all of the predictors in a multivariate model.
#' 
#' For convenience, this function will scale the input data automatically by
#' default, and backtransform output to the original scale. Automatic scaling
#' can be \code{"global"} (default), or \code{"local"}. The latter will scale
#' variables within (as opposed to across) populations, so is obviously only
#' applicable if there is more than 1 population. You can also set
#' \code{scaling="none"}. 
#' 
#' Note that if you use automatic scaling, the s-map coefficients and sigma2 value
#' will be for the scaled data, not the original data.
#' 
#' If you are fitting a nonstationary model, the time variable is NOT automatically
#' scaled when using \code{scaling}, since it has a separate weighting parameter 
#' \code{delta}. The interpretability of \code{delta} (i.e. whether it expresses 
#' relative or absolute weighting in time), depends on how time has been scaled, 
#' so depending on your goals, you might choose to leave time in its original units, 
#' or scale it from 0 to 1.
#' 
#' \strong{Predictions:}
#' 
#' The s-map routine automatically does leave-one-out prediction, which is included in the
#' output, so there is no separate \code{predictmethod="loo"}. The in-sample predictions
#' are not automatically returned, since they are of limited utility, however you can
#' obtain them by setting \code{returninsamp=TRUE}.
#' 
#' There are several addtional options for producing predictions:
#' \enumerate{
#'   \item Use \code{predictmethod="lto"} for leave-timepoint-out prediction using the training data. This will leave
#'   out values with the same time index across multiple populations, rather than each individual datapoint. 
#'   If there is only one population, \code{"lto"} will be equivalent to \code{"loo"}.
#'   \item Use \code{predictmethod="sequential"} for sequential (leave-future-out) prediction using the training data.
#'   \item If data frame \code{data} was supplied, supply data frame \code{newdata} containing same column names. 
#'   Column for \code{y} is optional, unless \code{E} and \code{tau} were supplied in lieu of \code{x}.
#'   \item If vectors/matrices were supplied for \code{y}, \code{x}, etc, equivalent vector/matrices \code{xnew}, 
#'   \code{popnew} (if \code{pop} was supplied), and \code{timenew} (optional).
#'   \code{ynew} is optional, unless \code{E} and \code{tau} were supplied in lieu of \code{x}.
#' }
#' 
#' After fitting, predictions can also be made using \code{\link{predict.Smap}}.
#' Predictions can plotted using \code{\link{plot.Smappred}}.
#' Get conditional responses using \code{\link{getconditionals}}.
#'
#' @param data A data frame, or matrix with named columns.
#' @param y The response variable (required). If \code{data} is supplied, a column name
#'   (character) of index (numeric). If \code{data} is not supplied, a numeric vector.
#' @param x Predictor variables. If \code{data} is supplied, column names
#'   (character vector) of indices (numeric vector). If \code{data} is not supplied, a numeric matrix 
#'   (or vector, if there is only one predictor variable). If \code{x} is not supplied,
#'   values for \code{E} and \code{tau} must be provided to construct it internally.
#' @param pop Identifies separate populations (optional, if not supplied, defaults to 1
#'   population). Population values can be either numeric, character, or factor. 
#'   If \code{data} is supplied, a column name (character) of index (numeric). 
#'   If \code{data} is not supplied, a vector (numeric, character, or factor).
#' @param time A time index (optional, if not supplied, defaults to a numeric index). 
#'   Important: The time index is not used for model fitting (timesteps are 
#'   assumed to be evenly spaced) but supplying \code{time} will be add these values to the output table, 
#'   which may be useful for later plotting purposes. If \code{data} is supplied, a column name
#'   (character) of index (numeric). If \code{data} is not supplied, a numeric vector.
#' @param E Embedding dimension. If supplied, will be used to constuct lags of \code{x} (or
#'   lags of \code{y} if \code{x} is not supplied).
#' @param tau Time delay. If supplied, will be used to constuct lags of \code{x} (or
#'   lags of \code{y} if \code{x} is not supplied).
#' @param thetafixed Fixes the theta parameter, if desired (see Details).
#' @param exclradius Exclusion radius. How many adjacent time steps on either side of the focal
#'   point should be excluded when fitting the local linear model? Defaults to 0 (none). 
#'   This setting is retained when making other types of predictions. Must be a positive integer.
#' @param ns Logical. Use nonstationary s-map (or not). Defaults to FALSE. 
#' @param nsvar If \code{ns=TRUE}, the variable supplied under \code{time} is used by default.
#'   If for some reason you want to use a different variable, supply it here (column name of 
#'   \code{data}).
#' @param deltafixed If \code{ns=TRUE}, fixes the delta parameter, if desired (see Details).
#' @param scaling How the variables should be scaled (see Details). Scaling can be \code{"global"}, 
#'   \code{"local"} (only applicable if more than 1 pop), or \code{"none"} (default). Scaling is applied 
#'   to \code{y} and each \code{x}. All outputs are backtransformed to original scale.
#' @param initpars Starting values for hyperparameters (see Details) in the order
#'   \code{c(theta,delta)}. Should be a numeric vector with length 2. Defaults used 
#'   if not supplied: \code{c(0,0)}.
#' @param augdata A data frame with augmentation data (see Details).
#' @param returninsamp Return the in-sample results? Default to FALSE. The in-sample
#'   results from s-map are generally not that useful, but you can get them if you want.
#' @inheritParams predict.Smap
#'
#' @return A list (class Smap and Smappred) with the following elements:
#' \item{theta}{Maximum likelihood value of theta.}
#' \item{delta}{For nonstationary models, maximum likelihood value of delta. Will be fixed to 0 for standard s-map models.}
#' \item{sigma2}{Variance estimate for the model (mean squared deviations), from leave-one-out.}
#' \item{grad}{Likelihood gradients of hyperparameters. Can be used to assess convergence.}
#' \item{iter}{Number of iterations until convergence was reached. Currently fixed at a max of 200.}
#' \item{ns}{Whether the model is nonstationary.}
#' \item{inputs}{Inputs and scaled inputs. Note that if you use \code{E} and \code{tau},
#'   the names of the predictors in the input data frame will be stored under
#' \code{x_names}, and the names of the lagged predictors corresponding to the
#'   inverse length scales will be stored under \code{x_names2}, provided these names exist.}
#' \item{scaling}{Scaling information.}
#' \item{looresults}{Data frame with leave-one-out out-of-sample predictions.}
#' \item{loocoefs}{Data frame of s-map coefficients from leave-one-out. The last column is the intercept.} 
#' \item{loofitstats}{Fit statistics for leave-one-out out-of-sample predictions.}
#' \item{loofitstatspop}{If >1 population, fit statistics for leave-one-out out-of-sample predictions by population.}
#' \item{call}{Function call. Allows use of \code{\link{update}}.}
#' @seealso \code{\link{predict.Smap}}, \code{\link{plot.Smappred}}, \code{\link{getconditionals}}
#' @references 
#' @examples 
#' 
#' @export
#' @importFrom MASS ginv
#' @keywords functions

fitSmap=function(data=NULL,y,x=NULL,pop=NULL,time=NULL,E=NULL,tau=NULL,thetafixed=NULL,
                 exclradius=0,ns=FALSE,nsvar=time,deltafixed=NULL,
                 scaling=c("global","local","none"),
                 initpars=NULL, augdata=NULL, returninsamp=FALSE, parallel=FALSE,
                 predictmethod=NULL,newdata=NULL,xnew=NULL,popnew=NULL,timenew=NULL,ynew=NULL) {
  
  cl <- match.call()
  
  #input checks
  if((!is.null(E) & is.null(tau)) | (is.null(E) & !is.null(tau))) {
    stop("Both E and tau must be supplied if generating lags internally.")
  }
  if(is.null(x) & (is.null(E) | is.null(tau))) {
    stop("x and/or E and tau must be supplied")
  }

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
    if(is.factor(pop)) { #fixes bug that happens if pop is a factor, augdata is used, and local scaling is used
      popd=factor(c(as.character(pop),as.character(popaug)), levels=levels(pop))
    } else {
      popd=c(pop,popaug)
    }
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
    # #issue warnings if data are not scaled properly
    # if(ymeans>0.1|ymeans<(-0.1)|ysds>1.1|ysds<(-0.9)) {
    #   warning("y is not scaled properly, model may be unreliable ", 
    #           "mean=",round(ymeans,3)," sd=",round(ysds,3),call. = F,immediate. = T)
    # }
    # if(any(xmeans>0.1)|any(xmeans<(-0.1))|any(xsds>1.1)|any(xsds<(-0.9))) {
    #   warning("one or more x is not scaled properly, model may be unreliable. ", 
    #           "means=",paste(round(xmeans,3),collapse=" "),
    #           " sds=",paste(round(xsds,3),collapse=" "),call. = F,immediate. = T)
    # }
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
    ymeans=tapply(yd,popd,mean,na.rm=T)
    ysds=tapply(yd,popd,sd,na.rm=T)
    xlist=split(as.data.frame(xd),popd)
    xmeans=lapply(xlist,function(x) apply(x,2,mean,na.rm=T))
    xsds=lapply(xlist,function(x) apply(x,2,sd,na.rm=T))
    up=unique(popd)
    xds=xd
    yds=yd
    for(i in 1:length(up)) {
      locmean=ymeans[as.character(up[i])==names(ymeans)]
      locsd=ysds[as.character(up[i])==names(ysds)]
      yds[popd==up[i]]=(yds[popd==up[i]]-locmean)/locsd
      for(j in 1:d) {
        locmean=xmeans[[which(as.character(up[i])==names(xmeans))]][j]
        locsd=xsds[[which(as.character(up[i])==names(xsds))]][j]
        xds[popd==up[i],j]=(xd[popd==up[i],j]-locmean)/locsd          
      }
    }
  }
  
  #standardize time
  #Time=(Time-min(Time))/max(Time-min(Time))
  
  #if initpars is not supplied, use default starting parameters
  if(is.null(initpars)) {
    initpars=c(0,0)
  }
  #fixed parameters
  if(!ns) { #if model is stationary, fix delta to 0
    deltafixed=0
  }
  
  #remove missing values
  completerows=complete.cases(cbind(yds,xds))
  Y=yds[completerows]
  X=as.matrix(xds[completerows,,drop=FALSE])
  Pop=popd[completerows]
  Time=timed[completerows]
  Primary=primary[completerows]
  
  if(parallel) {
    doParallel::registerDoParallel(cores=5)
    on.exit(doParallel::stopImplicitCluster())
  }
  
  #optimize model
  output=fmingrad_Rprop_smap(initpars,Y,X,Pop,Time,thetafixed,deltafixed,exclradius,parallel)
  # contains: list(theta,delta,sigma2,df,nllpost,grad,ypred,ysd,ypred_loo,ysd_loo,coefs,coefs_loo,iter)
  
  #backfill missing values
  ypred<-ypred_loo<-ysd<-ysd_loo<-yd*NA
  ypred[completerows]=output$ypred
  ysd[completerows]=output$ysd
  ypred_loo[completerows]=output$ypred_loo
  ysd_loo[completerows]=output$ysd_loo
  
  if(is.null(inputs$x_names2)) {labs<-inputs$x_names} else {labs<-inputs$x_names2}
  coefs<-coefs_loo<-matrix(NA,nrow=nrow(xd),ncol=ncol(xd)+1,dimnames=list(NULL,c(labs,"int")))
  coefs[completerows,]=output$coefs
  coefs_loo[completerows,]=output$coefs_loo

  #remove predictions for augmentation data
  ypred=ypred[primary]
  ysd=ysd[primary]
  ypred_loo=ypred_loo[primary]
  ysd_loo=ysd_loo[primary]
  
  coefs=coefs[primary,]
  coefs_loo=coefs_loo[primary,]
  
  #unscale predictions
  if(scaling=="global") {
    ypred=ypred*ysds+ymeans
    ypred_loo=ypred_loo*ysds+ymeans
    ysd=sqrt(ysd^2*ysds^2)
    ysd_loo=sqrt(ysd_loo^2*ysds^2)
  }
  if(scaling=="local") {
    for(i in 1:length(up)) {
      locmean=ymeans[as.character(up[i])==names(ymeans)]
      locsd=ysds[as.character(up[i])==names(ysds)]
      ypred[pop==up[i]]=ypred[pop==up[i]]*locsd+locmean
      ypred_loo[pop==up[i]]=ypred_loo[pop==up[i]]*locsd+locmean
      ysd[pop==up[i]]=sqrt(ysd[pop==up[i]]^2*locsd^2)
      ysd_loo[pop==up[i]]=sqrt(ysd_loo[pop==up[i]]^2*locsd^2)
    }
  }
  
  #create additional outputs
  output$ns=ns
  output$exclradius=exclradius
  output$inputs=c(inputs,list(y=y,x=x,yd=yd,xd=xd,pop=pop,time=time,
                              completerows=completerows,
                              Y=Y,X=X,Pop=Pop,Time=Time,Primary=Primary))
  output$scaling=list(scaling=scaling,ymeans=ymeans,ysds=ysds,xmeans=xmeans,xsds=xsds)
  
  if(returninsamp) {
    output$insampresults=data.frame(timestep=time,pop=pop,predmean=ypred,predsd=ysd,obs=y)
    output$insampcoefs=as.data.frame(coefs)
    output$insampfitstats=c(R2=getR2(y,ypred),
                            rmse=sqrt(mean((y-ypred)^2,na.rm=T)))
  }
  
  output$looresults=data.frame(timestep=time,pop=pop,predmean=ypred_loo,predsd=ysd_loo,obs=y)
  output$loocoefs=as.data.frame(coefs_loo)
  output$loofitstats=c(R2=getR2(y,ypred_loo),
                       rmse=sqrt(mean((y-ypred_loo)^2,na.rm=T)),
                       ln_L=-output$nll, df=output$df)
  
  if(length(unique(pop))>1) { #within site fit stats
    up=unique(pop)
    np=length(up)
    R2pop<-rmsepop<-numeric(np)
    names(R2pop)=up
    names(rmsepop)=up
    for(k in 1:np) {
      ind=which(pop==up[k])
      R2pop[k]=getR2(y[ind],ypred_loo[ind]) #this should this be y instead of yd in GP, although doesn't matter since pop is the length of y
      rmsepop[k]=sqrt(mean((y[ind]-ypred_loo[ind])^2,na.rm=T))
    }
    output$loofitstatspop=list(R2pop=R2pop,rmsepop=rmsepop)
  }
  
  output[c("df","nllpost","ypred","ysd","ypred_loo","ysd_loo","coefs","coefs_loo")]=NULL 
  #take these out of main list

  if(!is.null(predictmethod)|!is.null(xnew)|!is.null(newdata)) { #generate predictions if requested
    predictresults=predict.Smap(output,predictmethod,newdata,xnew,popnew,timenew,ynew) 
    output=c(output,predictresults)
  }
  
  output$call=cl
  class(output)=c("Smap","Smappred")
  return(output)
}

#' S-map test for nonstationarity
#'
#' loop over values of E, fitSmap model, return table and delta agg value
#' 
#' 
#' @inheritParams fitSmap
#' @param covar Any additional predictors (these will not be lagged)
#' @param Emin Minimum E, defaults to 1.
#' @param Emax Maximum E (required). Should generally not exceed sqrt(n).
#' @param tau Time delay. Defaults to 1.
#' @param summaryonly Return just the summary tables, not the lists of model output for
#'   each value of E.
#' 
#' @return A list (class NStest) with the following elements:
#' \item{delta_agg}{The nonstationarity test statistic.}
#' \item{bestR2}{Best loo R2 value of all the models fit.}
#' \item{results}{Data frame of hyperparameters, logL, df, R2, and weights for stationary and 
#'   nonstationary model at each value of E.}
#' \item{series}{Data frame of time and y (for plotting purposes).}
#' \item{mods}{If \code{summaryonly=FALSE}, list of all stationary models.}
#' \item{mods_NS}{If \code{summaryonly=FALSE}, list of all nonstationary models.}
#' 
#' 
Smap_NStest=function(data,y,covar=NULL,pop=NULL,time, Emin=1, Emax, tau=1, thetafixed=NULL,
                     exclradius=0,nsvar=time,
                     scaling=c("global","local","none"),
                     initpars=NULL, augdata=NULL, summaryonly=TRUE, parallel=FALSE) {
  
  df=as.data.frame(data) #plot won't work with class 'tbl', must be data.frame
  dflag=GPEDM::makelags(df, y=y, E=Emax, tau=tau, append = T)
  
  if(parallel) {
    
    doParallel::registerDoParallel(cores=5)
    on.exit(doParallel::stopImplicitCluster())
    `%dopar%` <- foreach::`%dopar%`
    
    results=foreach::foreach(i=Emin:Emax, .combine = rbind) %dopar% {
      #Stationary
      mods=fitSmap(dflag, y=y, x=c(paste0(y,"_",Emin:i),covar), pop=pop, time=time, thetafixed=thetafixed,
                   exclradius=exclradius, ns = FALSE, initpars=initpars, augdata=augdata)
      theta=mods$theta
      logL=mods$loofitstats["ln_L"]
      dfs=mods$loofitstats["df"]
      R2=mods$loofitstats["R2"]
      #Nonstationary
      mods_NS=fitSmap(dflag, y=y, x=c(paste0(y,"_",Emin:i),covar), pop=pop, time=time, thetafixed=thetafixed,
                      exclradius=exclradius, ns = TRUE, nsvar=nsvar, initpars=initpars, augdata=augdata)
      theta_NS=mods_NS$theta
      logL_NS=mods_NS$loofitstats["ln_L"]
      dfs_NS=mods_NS$loofitstats["df"]
      R2_NS=mods_NS$loofitstats["R2"]
      delta=mods_NS$delta
      results=data.frame(E=i, theta=theta, theta_NS=theta_NS, delta=delta, logL=logL, logL_NS=logL_NS, dfs=dfs, dfs_NS=dfs_NS, R2=R2, R2_NS=R2_NS)
      
    }
  } else { #not parallel
    results=data.frame(E=Emin:Emax, theta=NA, theta_NS=NA, delta=NA, logL=NA, logL_NS=NA, dfs=NA, dfs_NS=NA, R2=NA, R2_NS=NA)
    mods <- mods_NS <- list() 
    for(i in Emin:Emax) {
      #Stationary
      mods[[i]]=fitSmap(dflag, y=y, x=c(paste0(y,"_",Emin:i),covar), pop=pop, time=time, thetafixed=thetafixed,
                        exclradius=exclradius, ns = FALSE, initpars=initpars, augdata=augdata)
      results$theta[i]=mods[[i]]$theta
      results$logL[i]=mods[[i]]$loofitstats["ln_L"]
      results$dfs[i]=mods[[i]]$loofitstats["df"]
      results$R2[i]=mods[[i]]$loofitstats["R2"]
      #Nonstationary
      mods_NS[[i]]=fitSmap(dflag, y=y, x=c(paste0(y,"_",Emin:i),covar), pop=pop, time=time, thetafixed=thetafixed,
                           exclradius=exclradius, ns = TRUE, nsvar=nsvar, initpars=initpars, augdata=augdata)
      results$theta_NS[i]=mods_NS[[i]]$theta
      results$logL_NS[i]=mods_NS[[i]]$loofitstats["ln_L"]
      results$dfs_NS[i]=mods_NS[[i]]$loofitstats["df"]
      results$R2_NS[i]=mods_NS[[i]]$loofitstats["R2"]
      results$delta[i]=mods_NS[[i]]$delta
    }
  }
  rownames(results)=NULL
  
  #This part is wrong and needs updating
  W_E=pmax(0,exp(results$logL_NS)-exp(results$logL))
  if(sum(W_E)!=0) W_E=round(W_E/sum(W_E),8)
  results$W_E=W_E
  delta_agg=round(sum(results$delta*W_E),8)
  
  bestR2=max(results$R2, results$R2_NS)
  
  series=df[,c(time,y)]
  
  if(summaryonly | parallel) {
    out=list(delta_agg=delta_agg, bestR2=bestR2, results=results, series=series)
    class(out)=c("NStest")
    return(out)
  } else {
    out=list(delta_agg=delta_agg, results=results, series=series, mods=mods, mods_NS=mods_NS)
    class(out)=c("NStest")
    return(out)
  }
}

fmingrad_Rprop_smap = function(initpars,Y,X,Pop,Time,thetafixed,deltafixed,exclradius,parallel) {
  
  p=length(initpars)
  
  #optimization parameters for Rprop
  Delta0 = 0.1*rep(1,p)
  Deltamin = 1e-6
  Deltamax = 50
  eta_minus = 0.5;eta_minus=eta_minus-1
  eta_plus = 1.2;eta_plus=eta_plus-1   
  maxiter=200
  
  #initialize 
  output=getlikegrad_smap(initpars,Y,X,Pop,Time,thetafixed,deltafixed,exclradius,returnpreds=F,parallel=parallel)
  nll=output$nll
  grad=output$grad
  s=drop(sqrt(grad%*%grad))
  
  #loop starts here
  pa=initpars
  iter=0
  del=Delta0
  df=10
  
  while (s>.0001 & iter<maxiter & df>.0000001) {
    
    #step 1-move
    panew=pa-sign(grad)*del
    panew=pmax(0,panew) #pars cannot be less than 0
    output_new=getlikegrad_smap(panew,Y,X,Pop,Time,thetafixed,deltafixed,exclradius,returnpreds=F,reusegrad=grad,parallel=parallel)
    nll_new=output_new$nll
    grad_new=output_new$grad
    #panew=output_new$parst #prevent parst from changing if using fixedpars (probably not necessary?)
    s=drop(sqrt(grad_new%*%grad_new))
    df=abs(nll_new/nll-1)
    
    #step 2 - update step size
    gc=grad*grad_new
    tdel=del*(1+eta_plus*(gc>0)*1+eta_minus*(gc<0)*1)
    del=ifelse(tdel<Deltamin,Deltamin,ifelse(tdel>Deltamax,Deltamax,tdel))
    
    pa=panew
    grad=grad_new
    nll=nll_new
    iter=iter+1
  }
  
  #print(iter)
  paopt=pa
  output_new=getlikegrad_smap(paopt,Y,X,Pop,Time,thetafixed,deltafixed,exclradius,returnpreds=T,reusegrad=grad,parallel=parallel)
  output_new$iter=iter
  
  return(output_new)
}

#Function to get likelihood, gradient, df, predictions
getlikegrad_smap = function(pars, Y, X, Pop, Time, thetafixed, deltafixed, exclradius, returnpreds, reusegrad=NULL, parallel) {

  theta=pars[1]
  delta=pars[2]
  
  thetagrad <- deltagrad <- TRUE
  
  #fix parameters
  if(!is.null(thetafixed)) {
    theta=thetafixed
    if(!is.null(reusegrad)) { 
      dl_dtheta=reusegrad[1] 
      thetagrad=FALSE}
  }
  if(!is.null(deltafixed)) {
    delta=deltafixed
    if(!is.null(reusegrad)) { 
      dl_ddelta=reusegrad[2] 
      deltagrad=FALSE}
  }
  
  n = nrow(X)
  dimnames(X) = NULL
  
  if(!parallel) {
    
    dSSE_dtheta = 0
    dDOF_dtheta = 0
    dSSE_ddelta = 0
    dDOF_ddelta = 0
    
    SSE = 0
    dof = 0
    
    if(returnpreds) {
      ypred_loo <- ypred <- numeric(length = n)
      varp_loo <- varp <- numeric(length = n)
      coefs <- coefs_loo <- matrix(nrow = n, ncol = ncol(X)+1)
    }

    for(i in 1:n) {
      
      #exclude adjacent points
      exclude=(i-exclradius):(i+exclradius)
      #exclude=exclude[exclude>0 & exclude<=n]
      exclude=exclude[exclude %in% which(Pop==Pop[i])]
      
      Xi = X[i,,drop=F]
      Yi = Y[i]
      Timei = Time[i]
      
      Xnoi = X[-exclude,,drop=F]
      Ynoi = Y[-exclude]
      Timenoi = Time[-exclude]
      
      #loo
      loo_out = NSMap(Xnoi, Ynoi, Timenoi, Xi, Timei, theta, delta, returnpreds, returngrad=TRUE, thetagrad=thetagrad, deltagrad=deltagrad, Yi=Yi)
      if(returnpreds) {
        ypred_loo[i] = loo_out$prediction
        varp_loo[i] = loo_out$varp
        coefs_loo[i,] = loo_out$coefs        
      }

      #in sample
      insamp_out = NSMap(X, Y, Time, Xi, Timei, theta, delta, returnpreds, returngrad=TRUE, thetagrad=thetagrad, deltagrad=deltagrad, insampind=i)
      if(returnpreds) {
        ypred[i] = insamp_out$prediction
        varp[i] = insamp_out$varp
        coefs[i,] = insamp_out$coefs
      }
      
      SSE = SSE + loo_out$SSE
      dof = dof + insamp_out$dof
      
      if(thetagrad) {
        dSSE_dtheta = dSSE_dtheta + loo_out$dSSE_dtheta
        dDOF_dtheta = dDOF_dtheta + insamp_out$dDOF_dtheta
      }
      if(deltagrad) {
        dSSE_ddelta = dSSE_ddelta + loo_out$dSSE_ddelta
        dDOF_ddelta = dDOF_ddelta + insamp_out$dDOF_ddelta
      }
    }
    
  } else { #parallel
    `%dopar%` <- foreach::`%dopar%`
    results=foreach::foreach(i=1:n) %dopar% {
      #exclude adjacent points
      exclude=(i-exclradius):(i+exclradius)
      #exclude=exclude[exclude>0 & exclude<=n]
      exclude=exclude[exclude %in% which(Pop==Pop[i])]
      
      Xi = X[i,,drop=F]
      Yi = Y[i]
      Timei = Time[i]
      
      Xnoi = X[-exclude,,drop=F]
      Ynoi = Y[-exclude]
      Timenoi = Time[-exclude]
      
      #loo
      loo_out = NSMap(Xnoi, Ynoi, Timenoi, Xi, Timei, theta, delta, returnpreds, returngrad=TRUE, Yi=Yi)
      # loo_dhdtheta = loo_out$dhdtheta
      # loo_dhddelta =loo_out$dhddelta
      if(returnpreds) {
        ypred_loo = loo_out$prediction
        varp_loo = loo_out$varp
        coefs_loo = loo_out$coefs
      }
      
      #in sample
      insamp_out = NSMap(X, Y, Time, Xi, Timei, theta, delta, returnpreds, returngrad=TRUE, insampind=i)
      # hat_vec = insamp_out$H
      # insamp_dhdtheta = insamp_out$dhdtheta
      # insamp_dhddelta = insamp_out$dhddelta
      if(returnpreds) {
        ypred = insamp_out$prediction
        varp = insamp_out$varp
        coefs = insamp_out$coefs
      }
      
      # residual = Yi - loo_out$prediction
      # 
      # SSEi = (residual) ^ 2
      
      SSEi = loo_out$SSE
      # dofi = hat_vec[i] #trace of hat matrix
      dofi = insamp_out$dof
      
      # dSSE_dthetai = -2 * residual * (loo_dhdtheta %*% Ynoi)
      # dSSE_ddeltai = -2 * residual * (loo_dhddelta %*% Ynoi)
      # dDOF_dthetai = insamp_dhdtheta[i]
      # dDOF_ddeltai = insamp_dhddelta[i]
      dSSE_dthetai = loo_out$dSSE_dtheta
      dSSE_ddeltai = loo_out$dSSE_ddelta
      dDOF_dthetai = insamp_out$dDOF_dtheta
      dDOF_ddeltai = insamp_out$dDOF_ddelta
      
      results=list(SSEi=SSEi, dofi=dofi,
                   dSSE_dthetai=dSSE_dthetai, dSSE_ddeltai=dSSE_ddeltai,
                   dDOF_dthetai=dDOF_dthetai, dDOF_ddeltai=dDOF_ddeltai)
      if(returnpreds) {
        results=c(list(ypred_loo=ypred_loo, varp_loo=varp_loo, coefs_loo=coefs_loo,
                       ypred=ypred, varp=varp, coefs=coefs),results)
      }
      results
    }
    
    if(returnpreds) {
      ypred_loo=sapply(results, function(x) x$ypred_loo)
      varp_loo=sapply(results, function(x) x$varp_loo)
      coefs_loo=t(sapply(results, function(x) x$coefs_loo))
      ypred=sapply(results, function(x) x$ypred)
      varp=sapply(results, function(x) x$varp)
      coefs=t(sapply(results, function(x) x$coefs))
    }
    
    SSE=sum(sapply(results, function(x) x$SSEi))
    dof=sum(sapply(results, function(x) x$dofi))
    
    dSSE_dtheta=sum(sapply(results, function(x) x$dSSE_dthetai))
    dSSE_ddelta=sum(sapply(results, function(x) x$dSSE_ddeltai))
    dDOF_dtheta=sum(sapply(results, function(x) x$dDOF_dthetai))
    dDOF_ddelta=sum(sapply(results, function(x) x$dDOF_ddeltai))
  }
  
  #likelihood gradient for theta and delta (max prevents divide by 0 errors)
  if(thetagrad) {
    dl_dtheta = -1 * (-n/2) * ( dSSE_dtheta / max(SSE, 10e-6) + dDOF_dtheta / max(n-dof, 10e-6))
  }
  if(deltagrad) {
    dl_ddelta = -1 * (-n/2) * ( dSSE_ddelta / max(SSE, 10e-6) + dDOF_ddelta / max(n-dof, 10e-6))
  }
  
  #log likelihood
  E = ((-n/2) * (log(max(SSE, 10e-6) / max(n-dof, 10e-6)) + log(2*pi) + 1))
  
  if(!returnpreds) { #exit here, return only likelihood and gradient
    return(list(nll = -E, grad = c(dl_dtheta, dl_ddelta)))
  }

  #compute sd and other stuff, return predictions
  sigma2_loo=SSE/n
  ysd_loo=sqrt(sigma2_loo*varp_loo)
  
  sigma2=mean((Y-ypred)^2)
  ysd=sqrt(sigma2*varp)
  
  out=list(theta=theta, delta=delta, sigma2=sigma2_loo, 
           df=dof, nll= -E, grad = c(-dl_dtheta, -dl_ddelta),
           ypred=ypred, ysd=ysd, ypred_loo=ypred_loo, ysd_loo=ysd_loo,
           coefs=coefs, coefs_loo=coefs_loo)
  return(out)
}

# Implementation of NSMap
#       X - (matrix) training data
#       Y - (vector) targets
#       T - (vector) time for each row in X
#       x - (matrix, 1 row) current state to predict
#       t - (scalar) current time to predict
#       theta - (scalar) hyperparameter
#       delta - (scalar) hyperparameter
NSMap=function(X, Y, Time, x, t, theta, delta, 
               returnpreds=TRUE, returngrad=FALSE, thetagrad=TRUE, deltagrad=TRUE, Yi=NULL, insampind=NULL) {
  
  n = nrow(X)
  norms = drop(sqrt(laGP::distance(x,X))) #la.norm(X - x, axis=1) #euclidean distance of each point from x
  d = mean(norms)
  
  W = exp(-(theta*norms)/d - delta*(Time-t)^2) #vector of weights for each point
  #should this be (Time-t)/n ?
  
  # augment the training data with a column of 1s to allow for intercepts
  M = cbind(X, rep(1, times=n))
  xaug = cbind(x, 1)
  
  # weighted penrose inverse
  pinv = MASS::ginv(W*M) #pinv = la.pinv(W*M)
  #pinv2 = solve(t(M) %*% diag(W)^2 %*% M) %*% t(M) %*% diag(W) #same
  
  Wmat = matrix(W,nrow = nrow(pinv), ncol=ncol(pinv), byrow = T)
  PinvW = multMat(pinv, Wmat) #same as t(t(pinv) * W)
  H = xaug %*% PinvW #hat matrix row (needed for dof calculation)
  prediction = innerProd(H,Y) # drop(H %*% Y) #prediction for point x
  # prediction = drop(xaug %*% pinv %*% diag(W) %*% Y) #same
  
  if(returnpreds) {
    coefs = PinvW %*% Y #smap coefs row
    varp=diag(H %*% t(H))+1 #to get sd=sqrt(sigma2*varp)
    if(!returngrad) {
      return(list(prediction=prediction, varp=varp, coefs=t(coefs))) #output for predict()
    }
  }
  
  #gradient calculation
  
  # compute the derivatives of the hat matrix wrt theta and delta
  if(thetagrad) {
    dWdtheta = -1 * W * norms / d 
    dWdthetaM = matrix(dWdtheta, nrow = nrow(pinv), ncol=ncol(pinv), byrow = T)
    dthetapinv = multMat(pinv, dWdthetaM) # t(dWdtheta * t(pinv))
    dhdtheta = 2 * xaug %*% (dthetapinv - dthetapinv %*% M %*% PinvW)
  }
  if(deltagrad) {
    dWddelta = -1 * W * ((Time-t)^2)
    dWddeltaM = matrix(dWddelta, nrow = nrow(pinv), ncol=ncol(pinv), byrow = T)
    ddeltapinv = multMat(pinv, dWddeltaM) # t(dWddelta * t(pinv))
    dhddelta = 2 * xaug %*% (ddeltapinv - ddeltapinv %*% M %*% PinvW)
  }
  #same (divided by 2) ** why are lines above multipled by 2? Is it because of W being squared?
  # dhdtheta2 = xaug %*% pinv %*% (diag(dWdtheta) - diag(dWdtheta) %*% M %*% pinv %*% diag(W))
  # dhddelta2 = xaug %*% pinv %*% (diag(dWddelta) - diag(dWddelta) %*% M %*% pinv %*% diag(W))
  
  if(!is.null(insampind)) { #in samp
    out=list(dof=H[insampind], dDOF_dtheta=ifelse(thetagrad,dhdtheta[insampind],0), dDOF_ddelta=ifelse(deltagrad,dhddelta[insampind],0))
  }
  if(!is.null(Yi)) { #loo
    residual = Yi - prediction
    SSE = (residual) ^ 2
    if(thetagrad) {
      dSSE_dtheta = -2 * residual * innerProd(dhdtheta, Y)
    }
    if(deltagrad) {
      dSSE_ddelta = -2 * residual * innerProd(dhddelta, Y)
    }
    out=list(SSE=SSE, dSSE_dtheta=ifelse(thetagrad,dSSE_dtheta,0), dSSE_ddelta=ifelse(deltagrad,dSSE_ddelta,0))
  }
  
  if(returnpreds) { 
    return(c(list(prediction=prediction, varp=varp, coefs=t(coefs)), out))
  } else {
    return(out)
  }
  
  # else{
  #     # use least squares to solve if hat matrix derivates aren't needed,
  #     # as this is much faster than computing the penrose inverse
  #     # and gives the same output
  #     prediction = xaug %*% la.lstsq( W * M, W * Y, rcond=None)[1]
  #     return(prediction=prediction)
  #   }
}

#' Get predictions from an S-map model
#'
#' Obtain predictions from a fitted S-map model.
#' 
#' The s-map routine automatically does leave-one-out prediction, which is included in the
#' output, so there is no separate \code{predictmethod="loo"}.
#' 
#' There are several additional options for producing predictions:
#' \enumerate{
#'   \item Use \code{predictmethod="lto"} for leave-timepoint-out prediction using the training data. This will leave
#'   out values with the same time index across multiple populations, rather than each individual datapoint. 
#'   If there is only one population, \code{"lto"} will be equivalent to \code{"loo"}.
#'   \item Use \code{predictmethod="sequential"} for sequential (leave-future-out) prediction using the training data.
#'   \item If data frame \code{data} was supplied, supply data frame \code{newdata} containing same column names. 
#'   Column for \code{y} is optional, unless \code{E} and \code{tau} were supplied in lieu of \code{x}.
#'   \item If vectors/matrices were supplied for \code{y}, \code{x}, etc, equivalent vector/matrices \code{xnew}, 
#'   \code{popnew} (if \code{pop} was supplied), and \code{timenew} (optional).
#'   \code{ynew} is optional, unless \code{E} and \code{tau} were supplied in lieu of \code{x}.
#' }
#'
#' @param object Output from \code{fitGP}.
#' @param predictmethod Using the training data, \code{lto} does 
#'   leave-timepoint-out prediction, and \code{sequential} does sequential (leave-future-out)
#'   prediction.
#' @param newdata Data frame containing the same columns supplied in the
#'   original model.
#' @param xnew New predictor matrix or vector. Not required if \code{newdata} is supplied.
#' @param popnew New population vector. Not required if \code{newdata} is supplied.
#' @param timenew New time vector. Not required if \code{newdata} is supplied.
#' @param ynew New response vector. Optional, unless \code{E} and \code{tau} were supplied in 
#'   lieu of \code{x}. Not required if \code{newdata} is supplied.
#' @param ... Other args (not used).
#' 
#' @return A list (class Smappred) with the following elements:
#' \item{outsampresults}{Data frame with out-of-sample predictions.}
#' \item{outsampcoefs}{Data frame with S-map coefficients for the out-of-sample predictions. The last column is the intercept.}
#' \item{outsampfitstats}{Fit statistics for out-of-sample predictions.
#'   Only computed if using \code{"lto"} or \code{"sequential"}, if \code{y} is found in \code{newdata},
#'   or if \code{ynew} supplied (i.e. if the observed values are known).}
#' \item{outsampfitstatspop}{If >1 population, fit statistics for out-of-sample predictions by population.}
#' @export
#' @keywords functions
predict.Smap=function(object,predictmethod=NULL,newdata=NULL,xnew=NULL,popnew=NULL,timenew=NULL,ynew=NULL, ...) {
  
  theta=object$theta
  delta=object$delta
  sigma2=object$sigma2
  exclradius=object$exclradius
  X=object$inputs$X
  Y=object$inputs$Y
  y=object$inputs$y
  Pop=object$inputs$Pop
  Time=object$inputs$Time
  scaling=object$scaling$scaling
  ymeans=object$scaling$ymeans
  ysds=object$scaling$ysds
  xmeans=object$scaling$xmeans
  xsds=object$scaling$xsds
  
  if(!is.null(newdata)|!is.null(xnew)) {
    
    #if data frame is supplied, take columns from it
    if(!is.null(newdata)) {
      if(is.character(object$inputs$y_names)) {
        if(object$inputs$y_names %in% colnames(newdata)) { ynew=newdata[,object$inputs$y_names] }
      }
      if(is.numeric(object$inputs$y_names)) {
        if(object$inputs$y_names <= ncol(newdata)) { ynew=newdata[,object$inputs$y_names] }
      }
      if(!is.null(object$inputs$x_names)) { xnew=newdata[,object$inputs$x_names] }
      if(!is.null(object$inputs$pop_names)) { popnew=newdata[,object$inputs$pop_names] }
      if(!is.null(object$inputs$time_names)) { timenew=newdata[,object$inputs$time_names] }
    }
    
    #get number of predictions
    if(!is.null(ynew)) {
      Tslp=length(ynew)
    } else {
      Tslp=nrow(as.matrix(xnew))
    }
    
    #if pop is not supplied, create vector of 1s (all same population)
    if(is.null(popnew)) {
      popnew=rep(1,Tslp)
    }
    #if time is not supplied, make a numeric index (only used in output)
    if(is.null(timenew)) { 
      up=unique(popnew)
      timenew=numeric(Tslp)
      for(i in 1:length(up)) {
        timenew[popnew==up[i]]=1:length(timenew[popnew==up[i]])
      }
    }
    
    #if E tau supplied, generate lags of xnew (or ynew if xnew not supplied).
    if(!is.null(object$inputs$E)) {
      if(is.null(xnew)) {
        xnew=ynew
      } 
      xnew=makelags(y=xnew,pop=popnew,E=object$inputs$E,tau=object$inputs$tau,yname=object$inputs$x_names)
    }
    
    #make sure x is a matrix, not vector or data frame
    xnew=as.matrix(xnew)
    
    #rescale inputs
    if(scaling=="none") {
      xnews=xnew
    }
    if(scaling=="global") {
      xnews=xnew
      for(j in 1:ncol(xnew)) {
        xnews[,j]=(xnew[,j]-xmeans[j])/xsds[j]
      }
    }
    if(scaling=="local") {
      up=unique(popnew)
      xnews=xnew
      for(i in 1:length(up)) {
        for(j in 1:ncol(xnew)) {
          locmean=xmeans[[which(as.character(up[i])==names(xmeans))]][j]
          locsd=xsds[[which(as.character(up[i])==names(xsds))]][j]
          xnews[popnew==up[i],j]=(xnew[popnew==up[i],j]-locmean)/locsd          
        }
      }
    }
    
    #remove missing values
    completerowsnew=complete.cases(xnews)
    #exit if no valid rows
    if(all(completerowsnew==FALSE)) {
      out=list(outsampresults=data.frame(timestep=timenew,pop=popnew,predmean=NA,predfsd=NA,predsd=NA))
      if(!is.null(ynew)) {
        out$outsampresults$obs=ynew
      }
      class(out)="Smappred"
      return(out)
    }
    Xnew=as.matrix(xnews[completerowsnew,,drop=F])
    Popnew=popnew[completerowsnew]
    Timenew=timenew[completerowsnew]
    
    ymean <- varp <- numeric(length = nrow(Xnew))
    coefs <- matrix(nrow = nrow(Xnew), ncol = ncol(X)+1)
    
    for(i in 1:nrow(Xnew)) {
      Xi = Xnew[i,,drop=F]
      Timei = Timenew[i]
      
      out = NSMap(X, Y, Time, Xi, Timei, theta, delta)
      ymean[i] = out$prediction
      varp[i] = out$varp
      coefs[i,] = out$coefs
    }
    
    #backfill missing values
    ypred<-ysd<-ynew*NA
    ypred[completerowsnew]=ymean
    ysd[completerowsnew]=sqrt(sigma2*varp)
    
    coefsnew<-matrix(NA,nrow=nrow(xnew),ncol=ncol(xnew)+1,dimnames=list(NULL,c(colnames(xnew),"int")))
    coefsnew[completerowsnew,]=coefs

  } else {
    
    # predictmethod=match.arg(predictmethod)
    
    if(predictmethod=="loo") {
      message("loo predictions are computed automatically in fitSmap, no need to include it under predictmethod")
      return(NULL)
    }
    
    if(predictmethod=="lto") {
      
      Time=object$inputs$Time #time index
      timevals=unique(Time) #assuming timepoints are in order
      nt=length(timevals) #number of timepoints
      Primary=object$inputs$Primary
      nd=length(which(Primary))
      ymean=numeric(length=nd)*NA
      varp=numeric(length=nd)*NA
      coefs <- matrix(nrow = nd, ncol = ncol(X)+1)
      
      for(i in 1:nt) {
        ind=which(Time[Primary]==timevals[i])
        train=which(Time!=timevals[i]) #this should also exclude any dups in the aug data
        # Xtrain=X[train,]
        # Ytrain=Y[train]
        # Timetrain=Time[train]
          
        for(j in ind) {
          #exclude adjacent points
          exclude=(j-exclradius):(j+exclradius)
          exclude=exclude[exclude %in% which(Pop==Pop[j])] 
          
          trainj=setdiff(train, exclude)
          
          Xtrain=X[trainj,]
          Ytrain=Y[trainj]
          Timetrain=Time[trainj]
          
          Xi = X[j,,drop=F]
          Timei = Time[j]
          
          out = NSMap(Xtrain, Ytrain, Timetrain, Xi, Timei, theta, delta)
          ymean[j] = out$prediction
          varp[j] = out$varp
          coefs[j,] = out$coefs
        }
      }
      
      #backfill missing values
      ypred<-ysd<-y*NA
      ypred[object$inputs$completerows[1:length(y)]]=ymean
      ysd[object$inputs$completerows[1:length(y)]]=sqrt(sigma2*varp)
      
      coefsnew<-matrix(NA,nrow=length(y),ncol=ncol(X)+1,dimnames=list(NULL,c(colnames(X),"int")))
      coefsnew[object$inputs$completerows[1:length(y)],]=coefs

      popnew=object$inputs$pop
      timenew=object$inputs$time
      ynew=y
    }
    
    if(predictmethod=="sequential") {
      
      Time=object$inputs$Time #time index
      timevals=unique(Time) #assuming timepoints are in order
      nt=length(timevals) #number of timepoints
      Primary=object$inputs$Primary
      nd=length(which(Primary))
      ymean=numeric(length=nd)*NA
      varp=numeric(length=nd)*NA
      coefs <- matrix(nrow = nd, ncol = ncol(X)+1)
      
      for(i in 3:nt) { #start with at least 3 time steps
        ind=which(Time[Primary]==timevals[i])
        train=which(Time<timevals[i])
        # Xtrain=X[train,]
        # Ytrain=Y[train]
        # Timetrain=Time[train]
        
        for(j in ind) {
          #exclude adjacent points
          exclude=(j-exclradius):(j+exclradius)
          exclude=exclude[exclude %in% which(Pop==Pop[j])] 
          
          trainj=setdiff(train, exclude)
          
          Xtrain=X[trainj,]
          Ytrain=Y[trainj]
          Timetrain=Time[trainj]
          
          Xi = X[j,,drop=F]
          Timei = Time[j]
          
          out = NSMap(Xtrain, Ytrain, Timetrain, Xi, Timei, theta, delta)
          ymean[j] = out$prediction
          varp[j] = out$varp
          coefs[j,] = out$coefs
        }
      }
      
      #backfill missing values
      ypred<-ysd<-y*NA
      ypred[object$inputs$completerows[1:length(y)]]=ymean
      ysd[object$inputs$completerows[1:length(y)]]=sqrt(sigma2*varp)
      
      coefsnew<-matrix(NA,nrow=length(y),ncol=ncol(X)+1,dimnames=list(NULL,c(colnames(X),"int")))
      coefsnew[object$inputs$completerows[1:length(y)],]=coefs
      
      popnew=object$inputs$pop
      timenew=object$inputs$time
      ynew=y
    }
  }
  
  #unscale predictions
  if(scaling=="global") {
    ypred=ypred*ysds+ymeans
    ysd=sqrt(ysd^2*ysds^2)
  }
  if(scaling=="local") {
    up=unique(popnew)
    for(i in 1:length(up)) {
      locmean=ymeans[as.character(up[i])==names(ymeans)]
      locsd=ysds[as.character(up[i])==names(ysds)]
      ypred[popnew==up[i]]=ypred[popnew==up[i]]*locsd+locmean
      ysd[popnew==up[i]]=sqrt(ysd[popnew==up[i]]^2*locsd^2)
    }
  }
  
  #probably need to output xnew (combine with table?, only if 1 predictor?)
  out=list(outsampresults=data.frame(timestep=timenew,pop=popnew,predmean=ypred,predsd=ysd))
  out$outsampcoefs=as.data.frame(coefsnew)
  if(!is.null(ynew)) {
    out$outsampresults$obs=ynew
    out$outsampfitstats=c(R2=getR2(ynew,ypred), 
                          rmse=sqrt(mean((ynew-ypred)^2,na.rm=T)))
    if(length(unique(popnew))>1) { #within site fit stats
      up=unique(popnew)
      np=length(up)
      R2pop<-rmsepop<-numeric(np)
      names(R2pop)=up
      names(rmsepop)=up
      for(k in 1:np) {
        ind=which(popnew==up[k])
        R2pop[k]=getR2(ynew[ind],ypred[ind]) 
        rmsepop[k]=sqrt(mean((ynew[ind]-ypred[ind])^2,na.rm=T))
      }
      out$outsampfitstatspop=list(R2pop=R2pop,rmsepop=rmsepop)
    }
  }
  class(out)="Smappred"
  return(out)
}

# data=subset(GPEDM::thetalog2pop, Population=="PopA")
# dlags=GPEDM::makelags(data, y = "Abundance", E=1, tau=1, append = T)
# X=as.matrix(dlags[-(1:3),4])
# Y=dlags$Abundance[-(1:3)]
# Time=dlags$Time[-(1:3)]/47
# 
# par(mfrow=c(2,2))
# thetagrid=seq(0,10, by=0.5)
# thetagrad <- nllike <- numeric(length = length(thetagrid))
# for(i in seq_along(thetagrid)) {
#   test=getlikegrad_smap(pars = c(thetagrid[i],0),X = X, Y = Y, Time = Time, thetafixed = NULL, deltafixed = 0, returngradonly = T)
#   nllike[i]=test$nll
#   thetagrad[i]=test$grad[1]
# }
# plot(thetagrid,nllike)
# plot(thetagrid,thetagrad); abline(h=0)
# 
# deltagrid=seq(-4,5, by=0.5)
# deltagrad <- dnllike <- numeric(length = length(deltagrid))
# for(i in seq_along(deltagrid)) {
#   test=getlikegrad_smap(pars = c(8,deltagrid[i]),X = X, Y = Y, Time = Time, thetafixed = 8, deltafixed = NULL, returngradonly = T)
#   dnllike[i]=test$nll
#   deltagrad[i]=test$grad[2]
# }
# plot(deltagrid,dnllike,ylab="nllike")
# plot(deltagrid,deltagrad); abline(h=0)
# 
# test=getlikegrad_smap(pars = c(0,0),X = X, Y = Y, Time = Time, thetafixed = 8, deltafixed = NULL, returngradonly = F)
# 
# par(mfrow=c(1,2))
# 
# plot(Time, Y, type="o")
# lines(Time, test$ypred, col="red")
# lines(Time, test$ypred_loo, col="blue")
# legend(x="topright", legend=c("obs","insamp","loo"), lty=1, col=c("black","red","blue"), cex=0.7)
# 
# plot(X[,1],Y)
# points(X[,1], test$ypred, col="red")
# points(X[,1], test$ypred_loo, col="blue")
# legend(x="topright", legend=c("obs","insamp","loo"), pch=1, col=c("black","red","blue"), cex=0.7)
# 
# test2=fmingrad_Rprop_smap(c(0,0), Y, X, Time, thetafixed = NULL, deltafixed = 0)
# test3=fmingrad_Rprop_smap(c(0,0), Y, X, Time, thetafixed = NULL, deltafixed = NULL)
# test4=fmingrad_Rprop_smap(c(0,0), Y, X, Time, thetafixed = 8, deltafixed = NULL)


#unresolved:
# - squaring of weighting kernel, matching rEDM