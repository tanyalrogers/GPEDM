#' Fit an S-map model
#'
#' Fits a S-map (local linear) model with local weighting parameter(s) found using
#' gradient descent.
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
#' \strong{Missing values:} 
#' 
#' Missing values for \code{y} and \code{x} are allowed. Any rows containing
#' missing values will be excluded from the model fitting procedure (reinserted
#' as NAs in the output). If \code{pop} and \code{time} are supplied, they
#' should NOT contain missing values. 
#' 
#' \strong{Hyperparameters:} 
#' 
#' 
#'  
#' \strong{Scaling:} 
#' 
#' The model priors assume the response variable and all predictor variables
#' have approximately zero mean and unit variance, therefore it is important to
#' scale the data properly when fitting these models. For convenience, this
#' function will scale the input data automatically by default, and
#' backtransform output to the original scale. Automatic scaling can be
#' \code{"global"} (default), or \code{"local"}. The latter will scale variables
#' within (as opposed to across) populations, so is obviously only applicable if
#' there is more than 1 population. You can also scale the data yourself
#' beforehand in whatever manner you wish and set \code{scaling="none"}. In this
#' case, you will obviously have to do the backtransformation yourself.
#' 
#' \strong{Predictions:}
#' 
#' There are several options for producing predictions:
#' \enumerate{
#'   \item Use \code{predictmethod="loo"} for leave-one-out prediction using the training data.
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
#' It should be noted that \code{"loo"} is not a "true" leave-one-out, for although each prediction is
#' made by removing one of the training points, the hyperparameters are fit using all of the training data.
#' The same goes for \code{"sequential"}.
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
#' @param scaling How the variables should be scaled (see Details). Scaling can be \code{"global"} 
#'   (default), \code{"local"} (only applicable if more than 1 pop), or \code{"none"}. Scaling is applied 
#'   to \code{y} and each \code{x}. For a mix of scalings, do it manually beforehand and 
#'   set scaling to \code{"none"}. All outputs are backtransformed to original scale.
#' @param initpars Starting values for hyperparameters (see Details) in the order
#'   \code{c(phi,ve,sigma2,rho)}, in constrained (non-transformed) space. Should be a numeric 
#'   vector with length \code{ncol(x)+3}. Defaults used if not supplied: \code{c(rep(0.1,ncol(x)), 0.001, 1-0.001, 0.5)}
#' @param augdata A data frame with augmentation data (see Details).
#' @inheritParams predict.Smap
#'
#' @return A list (class GP and GPpred) with the following elements:
#' \item{pars}{Posterior mode of hyperparameters (named vector).}
#' \item{parst}{Posterior mode of hyperparameters (named vector), transformed to real line (used for optimization).}
#' \item{grad}{Likelihood gradients of hyperparameters at posterior mode. Can be used to assess convergence.}
#' \item{covm}{List containing various covariance matrices and the inverse covariance matrix.}
#' \item{iter}{Number of iterations until convergence was reached. Currently fixed at a max of 200.}
#' \item{inputs}{Inputs and scaled inputs. Note that if you use \code{E} and \code{tau},
#' the names of the predictors in the input data frame will be stored under
#' \code{x_names}, and the names of the lagged predictors corresponding to the
#' inverse length scales will be stored under \code{x_names2}, provided these names exist.}
#' \item{scaling}{Scaling information.}
#' \item{insampresults}{Data frame with in-sample predictions. \code{predfsd} is the standard
#' deviation of the GP function, \code{predsd} includes process error.}
#' \item{insampfitstats}{Fit statistics for in-sample predictions. Includes R2, rmse, ln_post 
#'   (log posterior likelihood evaluated at the mode), lnL_LOO (generalized cross-validation approximate
#'   leave-one-out negative log likelihood), and df (estimated degrees of freedom, equal to 
#'   the trace of the smoother matrix). lnL_LOO does not include the prior, so is not directly
#'   comparable to ln_post. For diagnostics, also includes likelihood components ln_prior, SS, logdet.}
#' \item{insampfitstatspop}{If >1 population, fit statistics (R2 and rmse) for in-sample predictions by population.}
#' \item{outsampresults}{Data frame with out-of-sample predictions (if requested). \code{predfsd} is the standard
#' deviation of the GP function, \code{predsd} includes process error.}
#' \item{outsampfitstats}{Fit statistics for out-of-sample predictions (if requested). 
#'   Only computed if using \code{"loo"} or \code{"sequential"}, if \code{y} is found in \code{newdata},
#'   or if \code{ynew} supplied (i.e. if the observed values are known).}
#' \item{outsampfitstatspop}{If >1 population, fit statistics for out-of-sample predictions (if requested) by population.}
#' \item{call}{Function call. Allows use of \code{\link{update}}.}
#' @seealso \code{\link{predict.GP}}, \code{\link{plot.GPpred}}, \code{\link{getconditionals}}, \code{\link{getrhomatrix}}
#' @references Munch, S. B., Poynor, V., and Arriaza, J. L. 2017. Circumventing structural uncertainty: 
#'   a Bayesian perspective on nonlinear forecasting for ecology. Ecological Complexity, 32: 134.
#' @examples 
#' yrand=rnorm(20)
#' testgp=fitGP(y=yrand,E=2,tau=1,predictmethod = "loo")
#' @export
#' @importFrom laGP distance
#' @import Matrix
#' @keywords functions

fitSmap=function(data=NULL,y,x=NULL,pop=NULL,time=NULL,E=NULL,tau=NULL,thetafixed=NULL,
                 ns=FALSE,nsvar=time,deltafixed=NULL,
                 scaling=c("none","global","local"),
                 initpars=NULL,augdata=NULL,
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
  if(is.null(fixedpars)) {
    fixedpars=rep(NA,times=2)
  }
  
  # if(length(fixedpars)!=2) {
  #   stop("Length of fixedpars must be number of predictors + 2 (phi, ve, sigma2)")
  # }
  # fixedpars_index=which(!is.na(fixedpars))
  # initpars[fixedpars_index]=fixedpars[fixedpars_index]
  
  #remove missing values
  completerows=complete.cases(cbind(yds,xds))
  Y=yds[completerows]
  X=as.matrix(xds[completerows,,drop=FALSE])
  Pop=popd[completerows]
  Time=timed[completerows]
  Primary=primary[completerows]
  
  #optimize model
  output=fmingrad_Rprop_smap(initpars,Y,X,Time,thetafixed,deltafixed)
  # contains: list(theta,delta,nllpost,grad,lnL_LOO,df,mean,cov,iter)
  
  #backfill missing values
  ypred<-ypred_loo<-ysd<-ysd_loo<-yd*NA
  ypred[completerows]=output$ypred
  ysd[completerows]=output$ysd
  ypred_loo[completerows]=output$ypred_loo
  ysd_loo[completerows]=output$ysd_loo

  #remove predictions for augmentation data
  ypred=ypred[primary]
  ysd=ysd[primary]
  ypred_loo=ypred_loo[primary]
  ysd_loo=ysd_loo[primary]
  
  #unscale predictions
  if(scaling=="global") {
    ypred=ypred*ysds+ymeans
    yfsd=sqrt(yfsd^2*ysds^2)
    ysd=sqrt(ysd^2*ysds^2)
  }
  if(scaling=="local") {
    for(i in 1:length(up)) {
      locmean=ymeans[as.character(up[i])==names(ymeans)]
      locsd=ysds[as.character(up[i])==names(ysds)]
      ypred[pop==up[i]]=ypred[pop==up[i]]*locsd+locmean
      yfsd[pop==up[i]]=sqrt(yfsd[pop==up[i]]^2*locsd^2)
      ysd[pop==up[i]]=sqrt(ysd[pop==up[i]]^2*locsd^2)
    }
  }
  
  #create additional outputs
  output$inputs=c(inputs,list(y=y,x=x,yd=yd,xd=xd,pop=pop,time=time,
                              completerows=completerows,
                              Y=Y,X=X,Pop=Pop,Time=Time,Primary=Primary))
  output$scaling=list(scaling=scaling,ymeans=ymeans,ysds=ysds,xmeans=xmeans,xsds=xsds)
  output$insampresults=data.frame(timestep=time,pop=pop,predmean=ypred,predsd=ysd,
                                  predmean_loo=ypred_loo, predsd_loo=ysd_loo, obs=y)
  output$insampfitstats=c(R2=getR2(y,ypred),
                          rmse=sqrt(mean((y-ypred)^2,na.rm=T)),
                          ln_post=-output$nllpost, ln_prior=output$lp,
                          lnL_LOO=output$lnL_LOO,df=output$df,SS=output$SS,logdet=output$logdet)
  if(length(unique(pop))>1) { #within site fit stats
    up=unique(pop)
    np=length(up)
    R2pop<-rmsepop<-numeric(np)
    names(R2pop)=up
    names(rmsepop)=up
    for(k in 1:np) {
      ind=which(pop==up[k])
      R2pop[k]=getR2(yd[ind],ypred[ind])
      rmsepop[k]=sqrt(mean((yd[ind]-ypred[ind])^2,na.rm=T))
    }
    output$insampfitstatspop=list(R2pop=R2pop,rmsepop=rmsepop)
  }
  
  output[c("lnL_LOO","df","nllpost","mean","cov","SS","logdet","lp")]=NULL #take these out of main list
  #may eventually want to exclude some more outputs
  
  if(!is.null(predictmethod)|!is.null(xnew)|!is.null(newdata)) { #generate predictions if requested
    predictresults=predict.GP(output,predictmethod,newdata,xnew,popnew,timenew,ynew) 
    output=c(output,predictresults)
  }
  
  output$call=cl
  class(output)=c("Smap","Smappred")
  return(output)
}

#loop over values of E, fitSmap model, return table and delta agg value
Smap_NStest=function(data, y, x, pop, time, nsvar=time, Emin=1, Emax, tau, scaling) {
 #not done 
}

#Rprop to find optimum delta and/or theta (Kenneth's version)
# optimizeG=function(X, Y, Time, trainingSteps=40, hp=c(0,0), fixed=c(FALSE, FALSE)) {
#   
#   err=0
#   gradPrev = rep(1,times=length(hp))
#   deltaPrev = rep(1,times=length(hp))
#   
#   for (count in range(trainingSteps)) {
#     errPrev = err
#     out = getlikegrad_smap(X, Y, Time, hp[1], hp[2])
#     grad = out$grad
#     err = out$E
#     
#     if (abs(err-errPrev) < 0.01 | count == trainingSteps-1)
#       break
#     
#     # c(dweights, deltaPrev, gradPrev) = calculateHPChange(grad, gradPrev, deltaPrev)
#     rhoplus = 1.2 # if the sign of the gradient doesn't change, must be > 1
#     rhominus = 0.5 # if the sign DO change, then use this val, must be < 1
#     
#     grad = grad / sqrt(sum(grad^2)) #la.norm(grad) # NORMALIZE, because rprop ignores magnitude
#     
#     s = grad * gradPrev # ratio between -1 and 1 for each param
#     spos = ceiling(s) # 0 for - vals, 1 for + vals
#     sneg = -1 * (spos - 1)
#     
#     delta = ((rhoplus * spos) + (rhominus * sneg)) * (deltaPrev)
#     dweights = (delta) * ((np.ceil(grad) - 0.5) * 2) # make sure signs reflect the orginal gradient
#     
#     # floor and ceiling the hyperparameters
#     for (i in 1:2) {
#       if (!fixed[i])
#         hp[i] = max(0, hp[i] + dweights[i])
#     }
#   }
#   return(list(hp, err))
#   #needs to return more stuff
# }

fmingrad_Rprop_smap = function(initpars,Y,X,Time,thetafixed,deltafixed) {
  
  #this uses the sign of the gradient to determine descent directions 
  #and an adaptive step size - supposedly smoother convergence than
  #conjugate gradients for GP optimization.
  
  p=length(initpars)
  
  #optimization parameters for Rprop
  Delta0 = 0.1*rep(1,p)
  Deltamin = 1e-6
  Deltamax = 50
  eta_minus = 0.5;eta_minus=eta_minus-1
  eta_plus = 1.2;eta_plus=eta_plus-1   
  maxiter=200
  
  #pre-compute distance matrices
  D=list()
  for(i in 1:ncol(X)) {
    D[[i]]=laGP::distance(X[,i])
  }  
  
  #initialize 
  output=getlikegrad_smap(initpars,Y,X,Time,thetafixed,deltafixed,returngradonly=T)
  nllpost=output$nllpost
  grad=output$grad
  s=drop(sqrt(grad%*%grad))
  
  #loop starts here
  pa=initparst
  iter=0
  del=Delta0
  df=10
  
  while (s>.0001 & iter<maxiter & df>.0000001) {
    
    #step 1-move
    panew=pa-sign(grad)*del
    output_new=getlikegrad(panew,Y,X,Time,thetafixed,deltafixed,returngradonly=T)
    nllpost_new=output_new$nllpost
    grad_new=output_new$grad
    #panew=output_new$parst #prevent parst from changing if using fixedpars (probably not necessary?)
    s=drop(sqrt(grad_new%*%grad_new))
    df=abs(nllpost_new/nllpost-1)
    
    #step 2 - update step size
    gc=grad*grad_new
    tdel=del*(1+eta_plus*(gc>0)*1+eta_minus*(gc<0)*1)
    del=ifelse(tdel<Deltamin,Deltamin,ifelse(tdel>Deltamax,Deltamax,tdel))
    
    pa=panew
    grad=grad_new
    nllpost=nllpost_new
    iter=iter+1
  }
  
  #print(iter)
  paopt=pa
  output_new=getlikegrad(paopt,Y,X,Time,thetafixed,deltafixed,returngradonly=F)
  output_new$iter=iter
  
  return(output_new)
}


#Function to get likelihood, gradient, df

# finds the gradient of the likelihood function with respect to our hyperparameters theta and delta
getlikegrad_smap = function(pars, X, Y, Time, thetafixed, deltafixed, returngradonly) {
  # we should be able to pull this off with two passes, once for leave one out and again leave all in.
  
  theta=pars[1]
  delta=pars[2]
  
  #fix parameters
  if(!is.null(thetafixed)) {
    theta=thetafixed
  }
  if(!is.null(deltafixed)) {
    delta=deltafixed
  }
  
  n = dim(X)[1]
  dimnames(X) = NULL
  
  dSSE_dtheta = 0
  dDOF_dtheta = 0
  dSSE_ddelta = 0
  dDOF_ddelta = 0
  
  SSE = 0
  dof = 0
  
  ypred_loo <- ypred <- varp_loo <- varp <- numeric(length = n)
  
  for(i in 1:n) {
    Xi = X[i,,drop=F]
    Yi = Y[i]
    Timei = Time[i]
    
    Xnoi = X[-i,,drop=F]
    Ynoi = Y[-i]
    Timenoi = Time[-i]
    
    #loo
    loo_out = NSMap(Xnoi, Ynoi, Timenoi, Xi, Timei, theta, delta)
    loo_dhdtheta = loo_out$dhdtheta
    loo_dhddelta =loo_out$dhddelta
    ypred_loo[i] = loo_out$prediction
    varp_loo[i] = loo_out$varp
    
    #in sample
    insamp_out = NSMap(X, Y, Time, Xi, Timei, theta, delta)
    hat_vec = insamp_out$H
    insamp_dhdtheta = insamp_out$dhdtheta 
    insamp_dhddelta = insamp_out$dhddelta
    ypred[i] = insamp_out$prediction
    varp[i] = insamp_out$varp
    
    residual = Yi - loo_out$prediction
    
    SSE = SSE + (residual) ^ 2
    dof = dof + hat_vec[i] #trace of hat matrix
    
    dSSE_dtheta = dSSE_dtheta + -2 * residual * (loo_dhdtheta %*% Ynoi)
    dSSE_ddelta = dSSE_ddelta + -2 * residual * (loo_dhddelta %*% Ynoi)
    dDOF_dtheta = dDOF_dtheta + insamp_dhdtheta[i]
    dDOF_ddelta = dDOF_ddelta + insamp_dhddelta[i]
  }

  #likelihood gradient for theta and delta
  # this is ugly, but we have to include the max stuff to prevent divide by 0 errors
  dl_dtheta = (-n/2) * ( dSSE_dtheta / max(SSE, 10e-6) + dDOF_dtheta / max(n-dof, 10e-6))
  dl_ddelta = (-n/2) * ( dSSE_ddelta / max(SSE, 10e-6) + dDOF_ddelta / max(n-dof, 10e-6))
  
  #log likelihood
  E = ((-n/2) * (log(max(SSE, 10e-6) / max(n-dof, 10e-6)) + log(2*pi) + 1))
  
  if(returngradonly) { #exit here, return only likelihood and gradient
    return(list(nll = -E, grad = c(dl_dtheta, dl_ddelta)))
  }
  
  #compute sd and other stuff, return predictions
  sigma2_loo=SSE/n
  ysd_loo=sqrt(sigma2*varp_loo)
  
  sigma2=mean((Y-ypred)^2)
  ysd=sqrt(sigma2*varp_loo)
  
  out=list(theta=theta, delta=delta, sigma2=sigma2, df=dof, 
           nllpost= -E,grad = c(dl_dtheta, dl_ddelta),
           ypred=ypred, ysd=ysd, ypred_loo=ypred_loo, ysd_loo=ysd_loo,
           SS=SS,logdet=logdet,lp=lp)
  return(out)

}

# Implementation of NSMap! Note that T and t(where) must be standardized 
# to be between 0 and 1.
#   Parameters
#       X - (ndarray) training data, (n,p) array of state space variables
#       Y - (ndarray) labels
#       T - (ndarray) time for each row in X
#       x - (ndarray) current state to predict from
#       t - (scalar) current time to predict from
#       theta - (scalar) hyperparameter
#       delta - (scalar) hyperparameter
#   Returns
#       scalar prediction
#           OR
#       (scalar prediction, hat matrix row) if return_hat is true
#           OR
#       (scalar prediction, hat matrix row, derivative of h wrt theta, derivative of h wrt delta)
#               if return_hat_derivatives is True
NSMap=function(X, Y, Time, x, t, theta, delta) {
  
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
  
  H = xaug %*% t(t(pinv) * W) #hat matrix row (needed for dof calculation)
  prediction = drop(H %*% Y) #prediction for point x
  # prediction = drop(xaug %*% pinv %*% diag(W) %*% Y) #same
  varp=diag(H %*% t(H))+1 #to get sd=sqrt(sigma2*varp)
  
  # compute the derivatives of the hat matrix wrt theta and delta
  dWdtheta = -1 * W * norms / d 
  dWddelta = -1 * W * ((Time-t)^2)
  
  dthetapinv = t(dWdtheta * t(pinv))
  ddeltapinv = t(dWddelta * t(pinv))
  
  dhdtheta = 2 * xaug %*% (dthetapinv - dthetapinv %*% M %*% t(t(pinv) * W))
  dhddelta = 2 * xaug %*% (ddeltapinv - ddeltapinv %*% M %*% t(t(pinv) * W))
  
  #same (divided by 2) ** why are lines above multipled by 2? Is it because of W being squared?
  # dhdtheta2 = xaug %*% pinv %*% (diag(dWdtheta) - diag(dWdtheta) %*% M %*% pinv %*% diag(W))
  # dhddelta2 = xaug %*% pinv %*% (diag(dWddelta) - diag(dWddelta) %*% M %*% pinv %*% diag(W))
  
  return (list(prediction=prediction, varp=varp, H=H, dhdtheta=dhdtheta, dhddelta=dhddelta))
  
  # else{
  #     # use least squares to solve if hat matrix derivates aren't needed,
  #     # as this is much faster than computing the penrose inverse
  #     # and gives the same output
  #     prediction = xaug %*% la.lstsq( W * M, W * Y, rcond=None)[1]
  #     return(prediction=prediction)
  #   }
}

#unresolved:
# - confirm gradient calc is right
# - squaring of weighting kernel, matching rEDM

#Predict function for Smap
predict.Smap=function(object,predictmethod=c("loo","lto","sequential"),newdata=NULL,xnew=NULL,popnew=NULL,timenew=NULL,ynew=NULL, ...) {
  #not done
}

par(mfrow=c(2,2))
thetagrid=seq(0,5, by=0.2)
thetagrad <- nllike <- numeric(length = length(thetagrid))
for(i in seq_along(thetagrid)) {
  test=getlikegrad_smap(pars = c(thetagrid[i],0),X = X, Y = Y, Time = Time, thetafixed = NULL, deltafixed = 0, returngradonly = T)
  nllike[i]=test$nll
  thetagrad[i]=test$grad[1]
}
plot(thetagrid,nllike)
plot(thetagrid,thetagrad); abline(h=0)

deltagrid=seq(-1,5, by=0.5)
deltagrad <- nllike <- numeric(length = length(deltagrid))
for(i in seq_along(deltagrid)) {
  test=getlikegrad_smap(pars = c(0,deltagrid[i]),X = X, Y = Y, Time = Time, thetafixed = 2.2, deltafixed = NULL, returngradonly = T)
  nllike[i]=test$nll
  deltagrad[i]=test$grad[1]
}
plot(deltagrid,nllike)
plot(deltagrid,deltagrad); abline(h=0)

