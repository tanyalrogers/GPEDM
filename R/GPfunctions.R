#' Fit a GP model
#'
#' Fits a GP model with separable length scales and automatic relevance determination (ARD).
#' Also fits a hierarchical version of the GP model if distinct populations are indicated
#' using \code{pop}.
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
#' The model uses a squared exponential covariance function. See the (not yet written)
#' vignette for mathematical details.
#' 
#' There is one inverse length scale \code{phi} estimated for each predictor
#' variable (input dimension). Values near 0 indicate that the predictor has
#' little influence on the response variable, and it was likely dropped by ARD.
#' Higher values of \code{phi} indicate greater nonlinearity.
#' 
#' The estimated variance terms are \code{ve} (process variance) and 
#' \code{sigma2} (pointwise prior variance). 
#' 
#' Hyperparameter \code{rho} is the dynamic correlation in the hierarchical GP
#' model, indicating similarity of dynamics across populations. If there is only
#' 1 population (e.g. if \code{pop} is not supplied), \code{rho} is irrelevant
#' (not used by the model) and will revert to the mode of its prior (~0.5). The
#' value of \code{rho} can be fixed to any value from 0.0001 to 0.9999 using
#' \code{rhofixed}, otherwise it is estimated from the data. Alternatively, a matrix
#' of fixed pairwise rho values can be supplied using \code{rhomatrix}. In this case,
#' the single rho parameter will also not be used and will revert to the mode of the
#' prior (~0.5). A pairwise rho matrix can be generated using \code{\link{getrhomatrix}},
#' or you can create a custom one (e.g. based on geographical distance).
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
#' After fitting, predictions can also be made using \code{\link{predict.GP}}.
#' Predictions can plotted using \code{\link{plot.GPpred}}.
#' Get conditional responses using \code{\link{getconditionals}}.
#' 
#' It should be noted that \code{"loo"} is not a "true" leave-one-out, for although each prediction is
#' made by removing one of the training points, the hyperparameters are fit using all of the training data.
#' The same goes for \code{"sequential"}.
#' 
#' For the out-of-sample prediction methods \code{"loo"}, \code{"lto"}, and
#' \code{newdata}, the partial derivatives of the fitted GP function at each
#' time point with respect to each input can be obtained by setting
#' \code{returnGPgrad = T}. If you want the in-sample gradient, pass the
#' original (training) data as back in as \code{newdata}. It automatic scaling
#' was used, the gradient will also be backtransformed to the original units.
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
#' @param modeprior This value is used by the phi prior and sets the expected number of modes over the unit interval. 
#'   Defaults to 1.
#' @param fixedpars Fixes values of the hyperparameters phi, ve, and sigma2 (if desired). Should be a numeric 
#'   vector with length \code{ncol(x)+2} in the order \code{c(phi,ve,sigma2)} in constrained 
#'   (non-transformed) space. Hyperparameters to estimate should be input as NA. The value of rho is
#'   fixed separately using argument \code{rhofixed}.
#' @param rhofixed Fixes the rho parameter, if desired (see Details).
#' @param rhomatrix Symmetrical square matrix of fixed pairwise rho values to use, with 1's on the diagonal. 
#'   The rows and columns should be named with the population identifier. 
#'   The output of \code{\link{getrhomatrix}} can be used here.
#'   All populations in the dataset must appear in \code{rhomatrix}. Populations in 
#'   \code{rhomatrix} that are not in the dataset are allowed and will not be used.
#' @param augdata A data frame with augmentation data (see Details).
#' @param kmsubset Use kmeans subsetting to reduce training data size;
#'   potentially useful for large datasets. Defaults to FALSE. Arguments
#'   \code{nclust}, \code{clustsize}, and \code{subtotal} apply only if
#'   \code{kmsubset=TRUE}. If TRUE, any prediction-related inputs are ignored,
#'   and the full training dataset will be used as \code{newdata} (will appear under
#'   outsamp results). See \code{\link{getsubset}} for more details.
#' @inheritParams predict.GP
#' @inheritParams getsubset
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
#' \item{GPgrad}{If \code{returnGPgrad=T}, a data frame with the partial derivatives of the 
#'   function with respect to each input.}
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

fitGP=function(data=NULL,y,x=NULL,pop=NULL,time=NULL,E=NULL,tau=NULL,
               scaling=c("global","local","none"),
               initpars=NULL,modeprior=1,fixedpars=NULL,rhofixed=NULL,
               rhomatrix=NULL,augdata=NULL,
               predictmethod=NULL,newdata=NULL,xnew=NULL,popnew=NULL,timenew=NULL,ynew=NULL,returnGPgrad=FALSE,
               kmsubset=FALSE,nclust=NULL,clustsize=NULL,subtotal=NULL) {

  cl <- match.call()
  
  #input checks
  if((!is.null(E) & is.null(tau)) | (is.null(E) & !is.null(tau))) {
    stop("Both E and tau must be supplied if generating lags internally.")
  }
  if(is.null(x) & (is.null(E) | is.null(tau))) {
    stop("x and/or E and tau must be supplied")
  }
  if(!is.null(rhofixed)) {
    if(!is.null(rhomatrix)) {
      stop("Supply either rhofixed or rhomatrix, not both")
    }
    if(rhofixed<0.0001) {
      rhofixed=0.0001
      message("rhofixed must between 0.0001 and 0.9999, setting to 0.0001")
    }
    if(rhofixed>0.9999) {
      rhofixed=0.9999
      message("rhofixed must between 0.0001 and 0.9999, setting to 0.9999")
    }
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
    } else if(kmsubset & is.null(data)) {
      xog=x #save original x for weird edge case of method B3 (x and E supplied) and kmclust=T
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
  
  #if initpars is not supplied, use default starting parameters
  if(is.null(initpars)) {
    initpars=c(rep(0.1,d),0.001, 1-0.001, 0.5)
  }
  #fixed parameters
  if(is.null(fixedpars)) {
    fixedpars=rep(NA,times=d+2)
  }
  if(length(fixedpars)!=d+2) {
    stop("Length of fixedpars must be number of predictors + 2 (phi, ve, sigma2)")
  }
  fixedpars_index=which(!is.na(fixedpars))
  initpars[fixedpars_index]=fixedpars[fixedpars_index]
  
  #transform parameters (also at line 506!)
  vemin=0.0001;sigma2min=0.0001
  vemax=4.9999;sigma2max=4.9999
  rhomin=0;rhomax=1
  initparst=c(log(initpars[1:d]),logit(initpars[d+1],vemin,vemax), 
              logit(initpars[d+2],sigma2min,sigma2max), logit(initpars[d+3],rhomin,rhomax))
  
  #if only one population, set transformed rho to 0
  if(length(unique(popd))==1) {initparst[d+3]=0}
  
  #store rhomatrix if supplied, check formatting, pop representation
  if(!is.null(rhomatrix)) {
    inputs$rhomatrix=rhomatrix
    up=unique(popd)
    if(is.null(colnames(rhomatrix)) | is.null(rownames(rhomatrix))) {
      colnames(rhomatrix)=up
      rownames(rhomatrix)=up
      message("rhomatrix entries assumed in order: ", paste(up, collapse = " "), "\nName rows and columns to eliminate ambiguity.")
    }
    if(!all(as.character(up) %in% colnames(rhomatrix))) {
      stop("Population names are not all present in rhomatrix.")
    }
  }
  
  #remove missing values
  completerows=complete.cases(cbind(yds,xds))
  #if subsetting, treat values not in the subset as missing values
  if(kmsubset) {
    completerows=getsubset(xds,yds,nclust,clustsize,subtotal)
    #change prediction data to full training dataset
    if(!is.null(predictmethod)|!is.null(xnew)|!is.null(newdata)) {
      message("Ignoring prediction-related inputs and setting full training dataset as newdata.\nUse predict() to get other predictions.")
    }
    predictmethod=NULL
    if(!is.null(data)) { 
      newdata=data 
    } else {
      #this case may pose some problems if using autogenerated lags. I think the following fixes it.
      ynew=y
      popnew=pop
      timenew=time
      if(is.null(E)) { #if E not provided (B1), use x as xnew
        xnew=x
      }
      #for B2, x will be regenerated from y and E (x omitted)
      if(exists("xog")) { 
        xnew=xog #for B3 (x and E supplied), need to supply original x
      }
    }
  }
  Y=yds[completerows]
  X=as.matrix(xds[completerows,,drop=FALSE])
  Pop=popd[completerows]
  Time=timed[completerows]
  Primary=primary[completerows]
  
  #optimize model
  output=fmingrad_Rprop(initparst,Y,X,Pop,modeprior,fixedpars,rhofixed,rhomatrix)
  # contains: list(parst,pars,nllpost,grad,lnL_LOO,df,mean,cov,covm,iter)
  
  #backfill missing values
  ypred<-yfsd<-ysd<-yd*NA
  ypred[completerows]=output$mean
  yfsd[completerows]=sqrt(diag(output$cov))
  ysd[completerows]=sqrt(diag(output$cov)+output$pars["ve"])
  
  #remove predictions for augmentation data
  ypred=ypred[primary]
  yfsd=yfsd[primary]
  ysd=ysd[primary]
  
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
  output$insampresults=data.frame(timestep=time,pop=pop,predmean=ypred,predfsd=yfsd,predsd=ysd,obs=y)
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
  
  if(!is.null(predictmethod)|!is.null(xnew)|!is.null(newdata)|kmsubset) { #generate predictions if requested
    predictresults=predict.GP(output,predictmethod,newdata,xnew,popnew,timenew,ynew,returnGPgrad) 
    output=c(output,predictresults)
  }
  
  output$call=cl
  class(output)=c("GP","GPpred")
  return(output)
}

fmingrad_Rprop = function(initparst,Y,X,pop,modeprior,fixedpars,rhofixed,rhomatrix) {
  
  #this uses the sign of the gradient to determine descent directions 
  #and an adaptive step size - supposedly smoother convergence than
  #conjugate gradients for GP optimization.
  
  p=length(initparst)
  
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
  output=getlikegrad(initparst,Y,X,pop,modeprior,fixedpars,rhofixed,rhomatrix,D,returngradonly=T)
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
    output_new=getlikegrad(panew,Y,X,pop,modeprior,fixedpars,rhofixed,rhomatrix,D,returngradonly=T)
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
  output_new=getlikegrad(paopt,Y,X,pop,modeprior,fixedpars,rhofixed,rhomatrix,D,returngradonly=F)
  output_new$iter=iter
  
  return(output_new)
}

getlikegrad=function(parst,Y,X,pop,modeprior,fixedpars,rhofixed,rhomatrix,D,returngradonly) {
  
  d=ncol(X) #embedding dimension
  
  #transform parameters from real line to constrained space
  vemin=0.0001;sigma2min=0.0001
  vemax=4.9999;sigma2max=4.9999
  rhomin=0.0001;rhomax=1
  phi=exp(parst[1:d])
  ve=(vemax-vemin)/(1+exp(-parst[d+1]))+vemin
  sigma2=(sigma2max-sigma2min)/(1+exp(-parst[d+2]))+sigma2min
  rho=(rhomax-rhomin)/(1+exp(-parst[d+3]))+rhomin
  
  #if only 1 pop or rhomatrix exists, fix rho to 0.5
  if(length(unique(pop))==1 | !is.null(rhomatrix)) {
    rhofixed=0.5
  }
  
  #fix parameters
  if(!is.null(rhofixed)) {
    rho=rhofixed
    parst[d+3]=logit(rho,rhomin,rhomax)
  }
  if(any(!is.na(fixedpars[1:d]))) {
    phi[which(!is.na(fixedpars[1:d]))]=fixedpars[which(!is.na(fixedpars[1:d]))]
    parst[1:d]=log(phi)
  }
  if(!is.na(fixedpars[d+1])) {
    ve=fixedpars[d+1]
    parst[d+1]=logit(ve,vemin,vemax)
  }
  if(!is.na(fixedpars[d+2])) {
    sigma2=fixedpars[d+2]
    parst[d+2]=logit(sigma2,sigma2min,sigma2max)
  }
  
  #derivative for rescaled parameters wrt inputs -- for gradient calculation
  dpars=c(phi,(ve-vemin)*(1-(ve-vemin)/(vemax-vemin)),
          (sigma2-sigma2min)*(1-(sigma2-sigma2min)/(sigma2max-sigma2min)),
          (rho-rhomin)*(1-(rho-rhomin)/(rhomax-rhomin)))
  
  #specify priors 
  priors=getpriors(phi,ve,sigma2,rho,modeprior,sigma2max,vemax)
  lp=priors$lp #log prior
  dlp=priors$dlp #derivative of log prior
  
  #get covariance matrix
  covm=getcov(phi,sigma2,rho,X,X,pop,pop,rhomatrix)
  #D=covm$D #distance matrices
  Cd=covm$Cd #covariance matrix
  Id=diag(ncol(Cd))
  Sigma=Cd+ve*Id #covariance matrix with process noise
  covm$Sigma=Sigma

  #get inverse covariance
  Sigma=(Sigma+t(Sigma))/2
  icov=getcovinv(Sigma)
  iKVs=icov$iKVs
  L=icov$L
  covm$iKVs=iKVs
  
  a=iKVs%*%Y
  like=-0.5*Y%*%a-sum(log(diag(L))) #log likelihood
  
  #calculate gradient
  dl=0*dlp
  vQ=0.5*matrix2vector(a%*%t(a)-iKVs)
  
  # dl[1:d]<-sapply(D,function(x) {
  #   dC=-multMat(x,Cd)
  #   vdC=matrix2vector(dC)
  #   0.5*innerProd(vQ,vdC) })
  
  for(i in 1:d) {
    if(is.na(fixedpars[i])) {
      dC=-multMat(D[[i]],Cd)
      vdC=matrix2vector(dC)
      dl[i]=0.5*innerProd(vQ,vdC)
    }
  }

  # dl[1:d]=sapply(D,function(x) {
  #   dC=-x*Cd
  #   vdC=as.vector(dC)
  #   0.5*vQ%*%vdC
  # })
  
  # dC=NULL
  # for(i in 1:d) {
  #   dC[[i]]=-D[[i]]*Cd
  #   vdC=as.vector(dC[[i]])
  #   dl[i]=0.5*vQ%*%vdC
  # }
  
  # dC=NULL
  # for(i in 1:d) {
  #   #D=laGP::distance(X[,i])
  #   #dC[[i]]=-D*Cd
  #   dC[[i]]=-D[[i]]*Cd
  #   dl[i]=0.5*vQ%*%as.vector(dC[[i]])
  # }
  
  if(is.na(fixedpars[d+1])) {
    vdC=matrix2vector(Id)
    # dl[d+1]=vQ%*%vdC
    dl[d+1]=innerProd(vQ,vdC)
  }

  if(is.na(fixedpars[d+2])) {
    vdC=matrix2vector(Cd/sigma2)
    # dl[d+2]=vQ%*%vdC
    dl[d+2]=innerProd(vQ,vdC)
  }
  
  #if rhofixed exists, rhomatrix exists, or only 1 population, don't compute grad for rho
  if(!is.null(rhofixed)) {
    dl[d+3]=0
  } else {
    popsame=covm$popsame #population simularity matrix
    dC=Cd*(1-popsame)/rho
    vdC=matrix2vector(dC)
    # dl[d+3]=vQ%*%vdC
    dl[d+3]=innerProd(vQ,vdC)
  }

  # dC[[d+1]]=Id
  # dl[d+1]=vQ%*%as.vector(dC[[d+1]])
  # dC[[d+2]]=Cd/sigma2
  # dl[d+2]=vQ%*%as.vector(dC[[d+2]])
  # dC[[d+3]]=Cd*(1-popsame)/rho  
  # dl[d+3]=vQ%*%as.vector(dC[[d+3]])
  
  #J is gradient in parameter space - need gradient in transformed parameters
  J=dl+dlp
  neglgrad=-J*dpars #gradient of posteror nll
  neglpost=-(like+lp) #posterior nll
  
  if(returngradonly) { #exit here, return only likelihood and gradient
    return(list(nllpost=neglpost,grad=neglgrad))
  } 
  
  #likelihood components
  SS=Y%*%a
  logdet=sum(log(diag(L)))
  
  #transformed parameters
  pars=c(phi,ve,sigma2,rho)
  
  mpt=Cd%*%a  #in sample mean
  Cdt=Cd-Cd%*%iKVs%*%Cd  #in sample covariance
  
  #name stuff
  names(pars)=c(paste0("phi",1:length(phi)),"ve","sigma2","rho")
  names(parst)=c(paste0("phi",1:length(phi)),"ve","sigma2","rho")
  names(neglgrad)=c(paste0("phi",1:length(phi)),"ve","sigma2","rho")
  
  lnL_LOO=0.5*sum(log(diag(iKVs)))-0.5*sum(a^2/diag(iKVs))
  df=sum(diag((Cd%*%iKVs))) #trace
  
  out=list(pars=pars,parst=parst,
           nllpost=neglpost,grad=neglgrad,
           lnL_LOO=lnL_LOO,df=df,mean=mpt,cov=Cdt,
           covm=covm,SS=SS,logdet=logdet,lp=lp)
  return(out)
}

getpriors=function(phi,ve,sigma2,rho,modeprior,sigma2max,vemax){
  #length scale
  lam_phi=(2^(modeprior-1))^2*pi/2 #variance for gaussian - pi/2 means E(phi)=1
  lp_phi=-0.5*sum(phi^2)/lam_phi
  dlp_phi=-(phi^1)/lam_phi
  #sigma2
  a_sigma2=2
  b_sigma2=2 #beta
  lp_sigma2=(a_sigma2-1)*log(sigma2/sigma2max)+(b_sigma2-1)*log(1-sigma2/sigma2max)
  dlp_sigma2=(a_sigma2-1)/(sigma2)-(b_sigma2-1)/(sigma2max-sigma2)
  #process noise
  a_ve=2
  b_ve=2 #beta
  lp_ve=(a_ve-1)*log(ve/vemax)+(b_ve-1)*log(1-ve/vemax)  
  dlp_ve=(a_ve-1)/(ve)-(b_ve-1)/(vemax-ve)
  #correlation
  a_rho=2
  b_rho=2 #beta
  lp_rho=(a_rho-1)*log(rho)+(b_rho-1)*log(1-rho)  
  dlp_rho=(a_rho-1)/rho-(b_rho-1)/(1-rho)
  
  lp=lp_phi+lp_ve+lp_sigma2+lp_rho #log prior
  dlp=c(dlp_phi,dlp_ve,dlp_sigma2,dlp_rho) #derivative of log prior
  
  return(list(lp=lp,dlp=dlp))
}

getcov=function(phi,sigma2,rho,X,X2,pop,pop2,rhomat=NULL) {
  #construct base covariance matrix
  Tsl=nrow(X) #time series length
  Tslp=nrow(X2)
  
  #dexp=2 #Covariance function is exponential (dexp=1) or squared-exponential (dexp=2)
  
  #returns squared eucl distance
  lC0=-laGP::distance(t(t(X2)*sqrt(phi)), t(t(X)*sqrt(phi)))
  
  # lC0=matrix(0,nrow=Tslp,ncol=Tsl)
  # D=NULL
  # for(i in 1:length(phi)) {
  #   D[[i]]=abs(outer(X2[,i],X[,i],"-"))^dexp #distance matrix for coord i
  #   lC0=lC0-phi[i]*D[[i]]
  # }
  
  popsame=(outer(pop2,pop,"=="))*1 #should work for numeric, chr, or factor
  
  if(!is.null(rhomat)) { #use fixed rhomatrix (assumes rows/cols are named)
    up=unique(pop)
    np=length(up)
    Rmat=matrix(1,nrow=Tslp,ncol=Tsl)
    for(i in 1:np) {
      for(j in 1:np) {
        indi=which(pop2==up[i])
        indj=which(pop==up[j])
        Rmat[indi,indj]=rhomat[as.character(up[i]),as.character(up[j])]
      }
    }
    Cd=sigma2*exp(lC0)*Rmat
  } else { #use single rho
    Cd=sigma2*exp(lC0)*(popsame+rho*(1-popsame)) #covariance matrix without obs noise
  }
  
  # return(list(D=D,popsame=popsame,Cd=Cd))
  return(list(popsame=popsame,Cd=Cd))
}

getcovinv=function(Sigma) {
  # iKVs=solve(Sigma)
  
  # #chol algorithm
  # L=chol(Sigma)
  # Linv=solve(L)
  # iKVs=Linv%*%t(Linv) #inverse covariance matrix
  
  mm <- as(Sigma, "dpoMatrix")
  L=chol(mm)
  iKVs=chol2inv(L) #inverse covariance matrix
  return(list(L=as.matrix(L),iKVs=as.matrix(iKVs)))
}

#' Get predictions from a GP model
#'
#' Obtain predictions from a fitted GP model. There are several options:
#' \enumerate{
#'   \item (Default) Use \code{predictmethod="loo"} for leave-one-out prediction using the training data.
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
#' It should be noted that \code{"loo"} is not a "true" leave-one-out, for although each prediction is
#' made by removing one of the training points, the hyperparameters are fit using all of the training data.
#' The same goes for \code{"sequential"} and \code{"lto"}.
#'
#' @param object Output from \code{fitGP}.
#' @param predictmethod Using the training data, \code{loo} does leave-one-out prediction, \code{lto} does 
#'   leave-timepoint-out prediction, and \code{sequential} does sequential (leave-future-out)
#'   prediction.
#' @param newdata Data frame containing the same columns supplied in the
#'   original model.
#' @param xnew New predictor matrix or vector. Not required if \code{newdata} is supplied.
#' @param popnew New population vector. Not required if \code{newdata} is supplied.
#' @param timenew New time vector. Not required if \code{newdata} is supplied.
#' @param ynew New response vector. Optional, unless \code{E} and \code{tau} were supplied in 
#'   lieu of \code{x}. Not required if \code{newdata} is supplied.
#' @param returnGPgrad Return the gradient (derivative) of the GP model at each time point with 
#'   respect to each input. This is only computed for out-of-sample predictions using \code{newdata},
#'   \code{loo}, or \code{lto}. If you want the in-sample gradient, pass the original dataset as
#'   \code{newdata}. Defaults to FALSE.
#' @param ... Other args (not used).
#' @return A list (class GPpred) with the following elements:
#' \item{outsampresults}{Data frame with out-of-sample predictions (if requested). \code{predfsd} is the standard
#' deviation of the GP function, \code{predsd} includes process error.}
#' \item{outsampfitstats}{Fit statistics for out-of-sample predictions.
#'   Only computed if using a \code{predictmethod}, if \code{y} is found in \code{newdata},
#'   or if \code{ynew} supplied (i.e. if the observed values are known).}
#' \item{outsampfitstatspop}{If >1 population, fit statistics for out-of-sample predictions by population.}
#' \item{GPgrad}{If \code{returnGPgrad=T}, a data frame with the partial derivatives of the 
#'  function with respect to each input.}
#' @export
#' @keywords functions
predict.GP=function(object,predictmethod=c("loo","lto","sequential"),newdata=NULL,
                    xnew=NULL,popnew=NULL,timenew=NULL,ynew=NULL,
                    returnGPgrad=FALSE, ...) { 
  
  iKVs=object$covm$iKVs
  phi=object$pars[grepl("phi",names(object$pars))]
  sigma2=object$pars[names(object$pars)=="sigma2"]
  rho=object$pars[names(object$pars)=="rho"]
  ve=object$pars[names(object$pars)=="ve"]
  X=object$inputs$X
  Y=object$inputs$Y
  y=object$inputs$y
  rhomatrix=object$inputs$rhomatrix
  Pop=object$inputs$Pop
  scaling=object$scaling$scaling
  ymeans=object$scaling$ymeans
  ysds=object$scaling$ysds
  xmeans=object$scaling$xmeans
  xsds=object$scaling$xsds
  
  if(!is.null(newdata)|!is.null(xnew)|(!is.null(ynew)&!is.null(object$inputs$E))) {
    
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
      class(out)="GPpred"
      return(out)
    }
    Xnew=as.matrix(xnews[completerowsnew,,drop=F])
    Popnew=popnew[completerowsnew]
    
    #get covariance matrix
    covmnew=getcov(phi,sigma2,rho,X,Xnew,Pop,Popnew,rhomatrix)
    Cs=covmnew$Cd #covariance matrix
    
    #get predictions
    ymean=Cs%*%(iKVs%*%Y)
    yvar=numeric(length = nrow(Xnew))
    for (i in 1:nrow(Xnew)) {
      yvar[i]=sigma2-Cs[i,]%*%iKVs%*%Cs[i,]
    }
    
    #get GP gradient
    if(returnGPgrad) {
      GPgrad=Xnew*NA
      for(i in 1:nrow(Xnew)) {
        Xnewi=Xnew[i,]
        Csi=Cs[i,,drop=FALSE]
        GPgrad[i,]=getGPgrad(phi = phi, Xnew = Xnewi, X = X, Cdi = Csi, iKVs = iKVs, Y = Y)
      }
      gpgrad=xnew*NA
      if(!is.null(object$inputs$x_names2))
        {xnames=object$inputs$x_names2} else {xnames=object$inputs$x_names}
      colnames(gpgrad)=paste0("d_", xnames)
      gpgrad[completerowsnew,]=GPgrad
    }
    
    #backfill missing values
    ypred<-ysd<-yfsd<-numeric(Tslp)*NA
    ypred[completerowsnew]=ymean
    yfsd[completerowsnew]=sqrt(yvar)
    ysd[completerowsnew]=sqrt(yvar+ve)
    
  } else {
    
    predictmethod=match.arg(predictmethod)
    
    if(predictmethod=="loo") {
      Cd=object$covm$Cd
      Sigma=object$covm$Sigma
      Time=object$inputs$Time
      Primary=object$inputs$Primary
      nd=length(which(Primary))
      ymean=numeric(length=nd)
      yvar=numeric(length=nd)
      if(returnGPgrad) {
        GPgrad=matrix(NA, nrow = nd, ncol=ncol(X))
      }
      for(i in 1:nd) {
        #check for duplicates in aug data
        dups=which(paste0(Pop[i],Time[i])==paste0(Pop[!Primary],Time[!Primary]))+nd
        Cdi=Cd[i,-c(i,dups),drop=F]
        S_noi=Sigma[-c(i,dups),-c(i,dups)]
        Y_noi=Y[-c(i,dups)]
        S_noi=(S_noi+t(S_noi))/2
        icov_noi=getcovinv(S_noi)
        iKVs_noi=icov_noi$iKVs
        ymean[i]=Cdi%*%iKVs_noi%*%Y_noi
        yvar[i]=sigma2-Cdi%*%iKVs_noi%*%t(Cdi)
        if(returnGPgrad) {
          Xi=X[i,]
          Xnoi=X[-c(i,dups),]
          GPgrad[i,]=getGPgrad(phi = phi, Xnew = Xi, X = Xnoi, Cdi = Cdi, iKVs = iKVs_noi, Y = Y_noi)
        }
      }
      
      #backfill missing values
      ypred<-ysd<-yfsd<-y*NA
      ypred[object$inputs$completerows[1:length(y)]]=ymean
      yfsd[object$inputs$completerows[1:length(y)]]=sqrt(yvar)
      ysd[object$inputs$completerows[1:length(y)]]=sqrt(yvar+ve)
      
      if(returnGPgrad) {
        if(!is.null(object$inputs$x_names2))
          {xnames=object$inputs$x_names2} else {xnames=object$inputs$x_names}
        gpgrad=matrix(NA,nrow=length(y),ncol=ncol(X),dimnames = list(NULL, paste0("d_", xnames)))
        gpgrad[object$inputs$completerows[1:length(y)],]=GPgrad
      }
      
      popnew=object$inputs$pop
      timenew=object$inputs$time
      ynew=y
    }
    
    #for lfo, need to exclude aug values with future time index
    # if(predictmethod=="lfo") {
    #   Time=object$inputs$Time #time index
    #   timevals=unique(Time) #assuming timepoints are in order
    #   nt=length(timevals) #number of timepoints
    #   Cd=object$covm$Cd
    #   Sigma=object$covm$Sigma
    #   Primary=object$inputs$Primary
    #   nd=length(which(Primary))
    #   ymean=numeric(length=nd)*NA
    #   yvar=numeric(length=nd)*NA
    #   for(i in 3:nt) {
    #     ind=which(Time[Primary]==timevals[i])
    #     train=which(Time<timevals[i])
    #     Cdi=Cd[ind,train,drop=F]
    #     S_noi=Sigma[train,train]
    #     Y_noi=Y[train]
    #     icov_noi=getcovinv(S_noi)
    #     iKVs_noi=icov_noi$iKVs
    #     ymean[ind]=Cdi%*%iKVs_noi%*%Y_noi
    #     yvar[ind]=diag(sigma2-Cdi%*%iKVs_noi%*%t(Cdi))
    #   }
    #   
    #   #backfill missing values
    #   ypred<-ysd<-yfsd<-y*NA
    #   ypred[object$inputs$completerows[1:length(y)]]=ymean
    #   yfsd[object$inputs$completerows[1:length(y)]]=sqrt(yvar)
    #   ysd[object$inputs$completerows[1:length(y)]]=sqrt(yvar+ve)
    #   
    #   popnew=object$inputs$pop
    #   timenew=object$inputs$time
    #   ynew=y
    # }
    
    if(predictmethod=="lto") {
      Time=object$inputs$Time #time index
      timevals=unique(Time) #assuming timepoints are in order
      nt=length(timevals) #number of timepoints
      Cd=object$covm$Cd
      Sigma=object$covm$Sigma
      Primary=object$inputs$Primary
      nd=length(which(Primary))
      ymean=numeric(length=nd)*NA
      yvar=numeric(length=nd)*NA
      if(returnGPgrad) {
        GPgrad=matrix(NA, nrow = nd, ncol=ncol(X))
      }
      for(i in 1:nt) {
        ind=which(Time[Primary]==timevals[i])
        train=which(Time!=timevals[i]) #this should also exclude any dups in the aug data
        Cdi=Cd[ind,train,drop=F]
        S_noi=Sigma[train,train]
        Y_noi=Y[train]
        S_noi=(S_noi+t(S_noi))/2
        icov_noi=getcovinv(S_noi)
        iKVs_noi=icov_noi$iKVs
        ymean[ind]=Cdi%*%iKVs_noi%*%Y_noi
        yvar[ind]=diag(sigma2-Cdi%*%iKVs_noi%*%t(Cdi))
        if(returnGPgrad) {
          for(j in ind) {
            Xi=X[j,]
            Xnoi=X[train,]
            Cdii=Cd[j,train,drop=F]
            GPgrad[j,]=getGPgrad(phi = phi, Xnew = Xi, X = Xnoi, Cdi = Cdii, iKVs = iKVs_noi, Y = Y_noi)
          }
        }
      }
      
      #backfill missing values
      ypred<-ysd<-yfsd<-y*NA
      ypred[object$inputs$completerows[1:length(y)]]=ymean
      yfsd[object$inputs$completerows[1:length(y)]]=sqrt(yvar)
      ysd[object$inputs$completerows[1:length(y)]]=sqrt(yvar+ve)
      
      if(returnGPgrad) {
        if(!is.null(object$inputs$x_names2))
          {xnames=object$inputs$x_names2} else {xnames=object$inputs$x_names}
        gpgrad=matrix(NA,nrow=length(y),ncol=ncol(X),dimnames = list(NULL, paste0("d_", xnames)))
        gpgrad[object$inputs$completerows[1:length(y)],]=GPgrad
      }
      
      popnew=object$inputs$pop
      timenew=object$inputs$time
      ynew=y
    }
    
    if(predictmethod=="sequential") {
      
      Time=object$inputs$Time #time index
      up=unique(Pop)
      np=length(up) #number of populations
      timevals=unique(Time) #assuming timepoints are in order
      nt=length(timevals) #number of timepoints
      S=object$covm$Cd #sigma2*exp(lC0).*(popsame+rho*(1-popsame))
      
      Primary=object$inputs$Primary
      nd=length(which(Primary))
      
      #mint=min(Time);maxt=max(Time)
      #timevals=mint:maxt #min to max time index
      #nt=maxt-mint+1
      
      # ymean1=matrix(0,nrow=nd,ncol=nt)
      # ysd1=(sigma2+ve)*matrix(1,nrow=nd,ncol=nt)
      
      ymean1=matrix(0,nrow=length(Y),ncol=nt)
      ysd1=(sigma2+ve)*matrix(1,nrow=length(Y),ncol=nt)
      
      if(returnGPgrad) {
        message("returnGPgrad not available for predictmethod='sequential'")
        returnGPgrad=FALSE
      }
      
      for(i in 1:(nt-1)) {
        ind=which(Time==timevals[i])
        I=diag(length(ind))
        SS=S[ind,ind]+ve*I;
        SS=(SS+t(SS))/2 #prevents Matrix from thinking it's not symmetric, numerical issue possibly
        icov_ss=getcovinv(SS)
        iK=icov_ss$iKVs
        k=iK%*%(Y[ind]-ymean1[ind,i])
        ymean1[,i+1]=ymean1[,i]+S[,ind]%*%k
        S=S-S[,ind]%*%iK%*%S[ind,]
        ysd1[,i+1]=sqrt(abs(diag(S))+ve)
      }
      ymean=rep(NA,length(Y))
      ysd2=rep((sigma2+ve),length(Y))
      for(i in 1:nt) {
        for(k in 1:np) {
          ind=which(Time==timevals[i] & Pop==up[k])
          if(length(ind)>0) {
            ymean[ind]=ymean1[ind,i]
            ysd2[ind]=ysd1[ind,i]
          }
        }
      }
      
      #backfill missing values
      ypred<-ysd<-yfsd<-y*NA
      ypred[object$inputs$completerows[1:length(y)]]=ymean[Primary]
      yfsd[object$inputs$completerows[1:length(y)]]=sqrt(ysd2^2-ve)[Primary]
      ysd[object$inputs$completerows[1:length(y)]]=ysd2[Primary]
      
      popnew=object$inputs$pop
      timenew=object$inputs$time
      ynew=y
    }
    
  }
  
  #unscale predictions
  if(scaling=="global") {
    ypred=ypred*ysds+ymeans
    yfsd=sqrt(yfsd^2*ysds^2)
    ysd=sqrt(ysd^2*ysds^2)
    if(returnGPgrad) {
      for(i in 1:ncol(gpgrad)) {
        gpgrad[,i]=gpgrad[,i]*ysds/object$scaling$xsds[i]
      }
    }
  }
  if(scaling=="local") {
    up=unique(popnew)
    for(i in 1:length(up)) {
      locmean=ymeans[as.character(up[i])==names(ymeans)]
      locsd=ysds[as.character(up[i])==names(ysds)]
      ypred[popnew==up[i]]=ypred[popnew==up[i]]*locsd+locmean
      yfsd[popnew==up[i]]=sqrt(yfsd[popnew==up[i]]^2*locsd^2)
      ysd[popnew==up[i]]=sqrt(ysd[popnew==up[i]]^2*locsd^2)
      if(returnGPgrad) {
        for(j in 1:ncol(gpgrad)) {
          xsds=object$scaling$xsds
          locsdx=xsds[[which(as.character(up[i])==names(xsds))]][j]
          gpgrad[popnew==up[i],j]=gpgrad[popnew==up[i],j]*locsd/locsdx
        }
      }
    }
  }
  
  #probably need to output xnew (combine with table?, only if 1 predictor?)
  out=list(outsampresults=data.frame(timestep=timenew,pop=popnew,predmean=ypred,predfsd=yfsd,predsd=ysd))
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
  if(returnGPgrad) {
    out$GPgrad=as.data.frame(gpgrad)
  }
  class(out)="GPpred"
  return(out)
}

logit=function(x, min=0, max=1) {
  log((x-min)/(max-x))
}

#' Calculate R-squared
#'
#' Calculates an R-squared value. Vectors for observed and predicted values 
#' must be the same length. Missing values are allowed. Any observed-predicted
#' pairs with missing values are removed before computations are done.
#' 
#' @details
#' Returned R-squared might be negative. This indicates that the prediction is 
#' worse than using the mean.
#'
#' @param obs Vector of observed values.
#' @param pred Vector of predicted values.
#' @return The R-squared value.
#' @export
#' @keywords functions
getR2=function(obs, pred) {
  d=na.omit(cbind(obs, pred))
  R2=1-sum((d[,1]-d[,2])^2)/sum((d[,1]-mean(d[,1]))^2)
  return(R2)
}

#' Generate delay vectors
#'
#' Create a lag matrix for a time series (or suite of time series) using given E
#' and tau values. Like in fitGP, either a data frame and column names or vectors/matrices
#' can be provided. Use \code{forecast=TRUE} to generate a forecast matrix. Use 
#' \code{vtimestep=TRUE} to use the variable timestep method.
#' 
#' @details 
#' 
#' When using the standard (fixed timestep) method, the response variable and
#' covariates can all be input under \code{y}, and lags will be generated for
#' all variables.
#'
#' When using the variable timestep method, it is necessary to differentiate
#' between the response variable and covariates, as they are handled
#' differently. The response variable should go under \code{y} and the
#' covariates should go under \code{x}. The output matrix will include lags of
#' Tdiff (time difference).
#'
#' An augmentation matrix for use with the variable timestep method can be
#' generated by setting \code{vtimestep=TRUE} and \code{augment=TRUE}. Under
#' default behavior, the augmentation matrix will include only the Tdiff
#' combinations observed in the original vtimestep matrix, up to \code{nreps}.
#' If you supply a vector \code{Tdiff_fore}, then the augmentation matrix will
#' include or all possible combinations of the Tdiff values supplied in
#' \code{Tdiff_fore}, even if they weren't in the original vtimestep matrix,
#' up to \code{nreps}. 
#' 
#' If generating forecasts, the output matrix will include a column for
#' population if there is more than one, and will include a time column if
#' \code{time} is supplied. The time increment is based on the minimum
#' difference between timepoints. If generating forecasts using the variable
#' timestep method, a vector of time units to forecast beyond the last timestep
#' can be provided under \code{Tdiff_fore}.
#' 
#'
#' @param data A data frame, or matrix with named columns.
#' @param y A vector containing a time series (usually the response variable).
#'   Alternatively, a matrix or data frame where each column is a time series
#'   (usually the response variable and covariates). In the case of multiple
#'   time series, lags will be generated for each variable. If \code{data} is
#'   supplied, the column names or indices.
#' @param pop A vector of populations (optional, if not supplied defaults to 1
#'   population). If \code{data} is supplied, the column name or index.
#' @param E Embedding dimension. Required.
#' @param tau Time delay. Required.
#' @param yname Optional, name of the variable if \code{y} is input as a vector, or variables
#'   if \code{y} is a matrix with unnamed columns. Otherwise not used.
#'   Defaults to "var" if NULL.
#' @param forecast Produce a forecast matrix instead of the standard training/test matrix.
#' @param vtimestep Use variable timestep method. 
#' @param x If using \code{vtimestep=TRUE}, put the response variable under use \code{y}, and
#'   covariates under \code{x} (see details). If \code{vtimestep=FALSE}, will be appended
#'   to \code{y}.
#' @param time Used to generate forecast time if \code{forecast=TRUE} and to calculate
#'   Tdiff if \code{vtimestep=TRUE}. If not supplied and \code{vtimestep=TRUE}, a numeric index is used.
#' @param augment If \code{vtimestep=TRUE}, produce augmentation lag matrix.
#' @param Tdiff_max If \code{vtimestep=TRUE}, the maximum Tdiff value to allow in production
#'   of the lag matrix or augmentation lag matrix. 
#' @param Tdiff_fore If \code{vtimestep=TRUE} and \code{forecast=TRUE}, vector of time units to 
#'   forecast beyond the last timestep. Defaults to the minimum time difference between consecutive
#'   timepoints, multipled by tau.
#' @param nreps If \code{augment=TRUE}, the max number of delay vectors for each
#'   Tdiff value.
#' @param append Return \code{data} with the lags appended (defaults to FALSE). Only
#'   relevant if \code{forecast=FALSE} and \code{augment=FALSE}.
#' @return A matrix with named columns, the appended number indicating the time lag. If
#'   \code{y} has named columns, named columns in the lag matrix will match. If
#'   generating forecasts, the output matrix will include a column for
#'   population if there is more than one, and will include a time column if
#'   \code{time} is supplied. If using the variable timestep method, the output matrix
#'   will include lags of Tdiff (time difference).
#' @export
#' @examples
#' set.seed(1)
#' yrand <- rnorm(20)
#' site <- rep(c("a","b"),each=10)
#' dfrand <- data.frame(firstvar=rnorm(20),secondvar=rnorm(20))
#' makelags(y=yrand,E=3,tau=1)
#' makelags(y=dfrand,E=2,tau=2)
#' makelags(y=dfrand,pop=site,E=2,tau=1)
#' makelags(y=yrand,pop=site,E=2,tau=3,forecast = TRUE, yname="SomeName",time=c(1:10,1:10))
#' dfrand2 <- cbind.data.frame(Time=c(1:10,1:10),Site=site,dfrand)
#' makelags(data=dfrand2, y=c("firstvar","secondvar"),pop="Site",E=2,tau=3)
#' makelags(data=dfrand2, y=c("firstvar","secondvar"),pop="Site",E=2,tau=3,
#' forecast = TRUE,time="Time")
#' @keywords functions
makelags=function(data=NULL,y,pop=NULL,E,tau,yname=NULL,
                  forecast=FALSE,vtimestep=FALSE,x=NULL,time=NULL,
                  augment=FALSE,Tdiff_max=NULL,Tdiff_fore=NULL,nreps=NULL,
                  append=FALSE) {
  
  #input checks
  if(!is.wholenumber(E) | !is.wholenumber(tau) | E<0 | tau<0) {
    stop("E and tau must be positive integers")
  }
  if((vtimestep==FALSE | forecast==TRUE) & augment==TRUE) {
    stop("augment==TRUE can only be used with vtimestep==TRUE and forecast==FALSE")
  }
  if(!is.null(Tdiff_fore) & !(vtimestep==TRUE & (forecast==TRUE | augment==TRUE))) {
    message("Tdiff_fore not used; only used with vtimestep==TRUE and either forecast==TRUE or augment==TRUE")
  }
  if(!is.null(Tdiff_max) & vtimestep==FALSE) {
    message("Tdiff_max not used; only used with vtimestep==TRUE")
  }
  if(!is.null(nreps) & (vtimestep==FALSE | augment==FALSE)) {
    message("nreps not used; only used with vtimestep==TRUE and augment==TRUE")
  }
  if(tau>1 & vtimestep==TRUE) {
    message("Using tau>1 with vtimestep method is not recommended")
  }
  
  #default names in case where data is not used and names are needed
  pname="pop"
  tname="time"
  xname="xvar"
  
  #if data frame is supplied, take columns from it and store names
  if(!is.null(data)) {
    yname=y
    y=data[,y]
    if(!is.null(x)) { 
      xname=x
      x=data[,x] 
    }
    if(!is.null(pop)) { 
      pname=pop
      pop=data[,pop]
    }
    if(!is.null(time)) {
      tname=time
      time=data[,time] 
    }
  }
  
  if(is.null(yname)) yname="var"
  
  y=as.matrix(y)
  num_vars=ncol(y)
  
  #in the case where y is a vector or unnamed matrix
  #assign names or reassign name from data
  if(is.null(colnames(y))) {
    if(num_vars==length(yname)) {
      colnames(y)=yname
    } else {
      colnames(y)=paste0(yname, seq_len(num_vars))
    }
  }
  #same for x (usually only relevant if there is only one x)
  if(!is.null(x)) {
    x=as.matrix(x)
    if(is.null(colnames(x))) {
      if(num_vars==length(xname)) {
        colnames(x)=xname
      } else {
        colnames(x)=paste0(xname, seq_len(ncol(x)))
      }
    }
  }
  
  #if pop is not supplied, create vector of 1s (all same population)
  if(is.null(pop)) {
    pop=rep(1,nrow(y))
  }
  up=unique(pop)
  
  #standard method
  if(vtimestep==FALSE) {
    
    #if covariates are specified
    if(!is.null(x)) { y=cbind(y,x) }
    
    output=makelags_subrt(y=y,pop=pop,E=E,tau=tau,forecast=forecast)
    
    if(forecast==FALSE) {
      if(append==TRUE & !is.null(data)) {
        output=cbind.data.frame(data, output)
      }
      return(output)
    }
    
    if(forecast==TRUE) {
      #append pop if more than 1
      if(length(up)>1) {
        popfore=rep(up,each=tau)
        output=cbind.data.frame(popfore,output)
        colnames(output)[1]=pname
      }
      #append time if supplied
      if(!is.null(time)) {
        timefore=numeric(length(up)*tau)
        popfore=rep(up,each=tau)
        for(k in 1:length(up)) { #populations
          ind=which(pop==up[k])
          indfore=which(popfore==up[k])
          timei=time[ind]
          Tdiff=diff(timei[(length(timei)-1):length(timei)])
          timefore[indfore]=timei[length(timei)]+(1:tau)*Tdiff
        }
        output=cbind.data.frame(timefore,output)
        colnames(output)[1]=tname
      }
      return(output)
    }
  }
  
  #variable timestep method
  if(vtimestep==TRUE) {
    
    #if time is not supplied, make a numeric index
    if(is.null(time)) { 
      time=y*0
      for(i in 1:length(up)) {
        time[pop==up[i]]=1:length(time[pop==up[i]])
      }
    }
    
    #if covariates are specified
    if(!is.null(x)) { xy=cbind(y,x) 
    } else { xy=y }
    
    outlist=NULL
    for(k in 1:length(up)) {
      pind=which(pop==up[k])
      xyi=xy[pind,,drop=F]
      timei=time[pind]
      nm=which(apply(xyi, 1, function(x) all(!is.na(x))))
      out=matrix(nrow = nrow(xyi), ncol = E*(ncol(xyi)+1))
      vnames=colnames(xyi)
      colnames(out)=paste(rep(c(vnames,"Tdiff"),each=E),rep(seq(tau,tau*E,tau), times=length(vnames)+1), sep = "_")
      
      if(forecast) {
        if(is.null(Tdiff_fore)) {
          #Tdiff_fore=diff(timei[(length(timei)-1):length(timei)])
          Tdiff_fore=min(diff(timei))*tau
        }
        foredf=matrix(NA,ncol = ncol(xyi), nrow = length(Tdiff_fore))
        colnames(foredf)=colnames(xyi)
        xy2=rbind(xyi, foredf)
        time2=c(timei, timei[length(timei)]+Tdiff_fore)
        out=matrix(nrow = nrow(xyi)+length(Tdiff_fore), ncol = E*(ncol(xyi)+1))
        colnames(out)=paste(rep(c(vnames,"Tdiff"),each=E),rep(seq(tau,tau*E,tau), times=length(vnames)+1), sep = "_")
        indfore=(nrow(xyi)+1):nrow(xy2)
      } else {
        xy2=xyi
        time2=timei
      }
      
      for(i in 1:nrow(xy2)) {
        if(i>E*tau) {
          index=numeric(length = E+1)
          index[1]=i
          indiff=i-nm
          for(Ei in 1:E) {
            newind=nm[which.min(indiff[indiff>=tau*Ei & nm<=index[Ei]-tau])]
            if(length(newind)!=0) {
              index[Ei+1]=newind
              out[i,0:(length(vnames)-1)*E+Ei]=as.numeric(xyi[index[Ei+1],])
              out[i,length(vnames)*E+Ei]=time2[index[Ei]]-time2[index[Ei+1]]
            }
          }
        }
      }
      
      if(forecast) {
        out=cbind.data.frame(time2[indfore],out[indfore,,drop=F])
        colnames(out)[1]=tname
        if(length(up)>1) {
          out=cbind.data.frame(rep(up[k],nrow(out)),out)
          colnames(out)[1]=pname
        }
        outlist[[k]]=out
      } else {
        outlist[[k]]=out
      }
    }
    output=do.call(rbind,outlist)
    
    #append to data
    if(augment==FALSE & append==TRUE & !is.null(data)) {
      output=cbind.data.frame(data, output)
    }
    
    if(!is.null(Tdiff_max)) {
      cols=grep("Tdiff_",colnames(output))
      output[, cols][output[, cols] > Tdiff_max] <- NA
    }
    
    #augmentation
    if(augment==TRUE) {
      
      if(is.null(nreps)) {
        message("defaulting to nreps=10")
        nreps=10
      }
      
      cols=grep("Tdiff_",colnames(output))
      auglist=NULL
      
      for(k in 1:length(up)) {
        pind=which(pop==up[k])
        xyi=xy[pind,,drop=F]
        yi=y[pind,]
        timei=time[pind]
        outputi=outlist[[k]]
        
        #get existing unique combinations of Tdiffs
        outputi2=cbind(outputi,yi)
        outbase=as.data.frame(na.omit(outputi2))[,cols]
        outbase$combo=apply(outbase,1,paste,collapse="")
        outtab=as.data.frame(table(combo=outbase$combo))
        outtab=unique(merge(outbase,outtab))
        outtab=outtab[,-1]
        
        #use Tdiff_fore to get missing combos
        if(!is.null(Tdiff_fore)) {
          Tgrid=expand.grid(c(rep(list(Tdiff_fore),E)))
          colnames(Tgrid)=paste0("Tdiff_",1:E)
          outtab=merge(outtab,Tgrid,all = T)
          outtab$Freq[is.na(outtab$Freq)]=0
        }
        
        if(!is.null(Tdiff_max)) {
          outtab=outtab[apply(outtab[,1:E],1,max)<=Tdiff_max,]
        }
        
        #iterate though all combos until nmax or out of xyi possibilities
        #retain valid lag vectors and associated time index
        augmat=outputi[0,]
        newrow0=outputi[1,]
        timeaug=numeric()
        yaug=numeric()
        outtab$Freq_new=NA
        rownames(outtab)=NULL
        for(i in 1:nrow(outtab)) {
          freq=outtab$Freq[i]
          remainingrows=nrow(xyi)
          randind=sample(1:nrow(xyi),nrow(xyi),replace = F)
          randindi=1
          while(freq<nreps & randindi<=nrow(xyi)) {
            if(randind[randindi]>E*tau & !is.na(yi[randind[randindi]])) {
              index=numeric(length = E+1)
              index[1]=randind[randindi]
              tdiffs=as.numeric(outtab[i,grep("Tdiff_",colnames(outtab))])
              newrow=newrow0
              for(Ei in 1:E) {
                index[Ei+1]=index[Ei]-tdiffs[Ei]
                if(index[Ei+1]>0) {
                  newrow[0:(length(vnames)-1)*E+Ei]=as.numeric(xyi[index[Ei+1],])
                  newrow[length(vnames)*E+Ei]=time[index[Ei]]-time[index[Ei+1]]
                }
              }
              if(all(!is.na(newrow)) & !duplicated(rbind(outputi,newrow))[nrow(outputi)+1]) {
                augmat=rbind(augmat,newrow)
                timeaug=c(timeaug,time[randind[randindi]])
                yaug=c(yaug,yi[randind[randindi]])
                freq=freq+1
              }
            }
            randindi=randindi+1
          }
          outtab$Freq_new[i]=freq
        }
        rownames(augmat)=NULL
        out=cbind.data.frame(timeaug,yaug,augmat)
        colnames(out)[1:2]=c(tname,yname)
        if(length(up)>1) {
          out=cbind.data.frame(rep(up[k],nrow(out)),out)
          colnames(out)[1]=pname
        }
        auglist[[k]]=out
        cat("Population ",up[k],"\n")
        if(nrow(out)==0) {
          message("Additional combinations not available. Check nreps or Tdiff_fore for more options.")
        }
        print(outtab)        
      }
      output=do.call(rbind,auglist)
    }
    return(output)
  }
}

makelags_subrt=function(y, #matrix
                        pop,E,tau,forecast) {
  up=unique(pop)
  num_vars=ncol(y)
  
  output=matrix(NA, nrow = nrow(y), ncol = num_vars*E)
  outputfore=matrix(NA, nrow = length(up)*tau, ncol = num_vars*E)
  popfore=rep(up,each=tau)
  
  col_names=character(num_vars*E)
  
  for(k in 1:length(up)) { #populations
    col_index=1
    for(j in 1:num_vars) { #variables
      ind=which(pop==up[k])
      indfore=which(popfore==up[k])
      ts=y[ind,j]
      if(length(ts)<=E*tau) stop("time series length <= E*tau for pop ",up[k])
      for (i in 1:E) {
        tslag=c(rep_len(NA, tau*i),ts[1:(length(ts) - tau*i)])
        output[ind,col_index]=tslag
        tsfore=ts[(length(ts) - tau*i+1):(length(ts)- tau*(i-1))]
        outputfore[indfore,col_index]=tsfore
        col_names[col_index]=paste0(colnames(y)[j],"_", i * tau)
        col_index=col_index + 1
      }
    }
  }
  colnames(output)=col_names
  colnames(outputfore)=col_names
  
  if(forecast==FALSE) return(output)
  else return(outputfore)
}

is.wholenumber=function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

getGPgrad=function(phi, Xnew, X, Cdi, iKVs, Y) {
  grad = t(-2 * t(t(Xnew-t(X)) %*% diag(phi)) %*% (t(Cdi) * iKVs %*% Y))
}