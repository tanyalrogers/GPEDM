#' Fit a GP model
#'
#' Fits a GP model with separable length scales and automatic relevance determination (ARD).
#' Also fits a hierarchical version of the GP model if distinct populations are indicated
#' using \code{pop}.
#' There are several ways to specify the training data:
#' \enumerate{
#'   \item data frame \code{data}, plus column names or indices for \code{yd} and \code{xd}
#'   \item data frame \code{data}, plus column name or index for \code{yd}, and values for \code{E} and \code{tau} 
#'   \item vector \code{yd}, plus matrix or vector \code{xd}
#'   \item vector \code{yd}, plus values for E and tau
#' }
#' Arguments \code{pop} and \code{time} are optional in all of the above cases. 
#' They should be column names or indices in cases 1 and 2, vectors in cases 3 and 4.
#' If \code{pop} and \code{time} are supplied, they should NOT contain missing values.
#' This function will also, optionally, produce predictions.
#' See Details for more information.
#' 
#' @details 
#' 
#' The data are assumed to be should be in time order. This is particularly important if the E/tau option
#' or \code{predictmethod="sequential"} is used.
#' 
#' \strong{Missing values:} 
#' 
#' Missing values for \code{yd} and \code{xd} are allowed. Any rows containing
#' missing values will be excluded from the model fitting procedure (reinserted
#' as NAs in the output). If \code{pop} and \code{time} are supplied, they
#' should NOT contain missing values. 
#' 
#' \strong{Hyperparameters:} 
#' 
#' The model uses the following covariance function:
#' 
#' (Insert equation here.)
#' 
#' There is one inverse length scale \code{phi} estimated for each predictor variable. 
#' Values near 0 indicate that the predictor has little influence on the response variable, 
#' and it was likely dropped by ARD. Higher values of \code{phi} indicate greater nonlinearity.
#' 
#' The estimated variance terms are \code{ve} (observation variance) and 
#' \code{tau} (pointwise function variance). 
#' 
#' Hyperparameter \code{rho} is the dynamic
#' correlation in the hierarchical GP model, indicating similarity of dynamics
#' across populations. If there is only 1 population (e.g. if \code{pop} is not
#' supplied), \code{rho} is irrelevant (not used by the model) and will revert
#' to the mode of the prior (~0.5). The value of \code{rho} can be fixed to any
#' value from 0.0001 to 1 using \code{rhofixed}, otherwise it is estimated from the data.
#' See Munch et al. (2017) for more details.
#'  
#' \strong{Scaling:} 
#' 
#' The model priors assume the response variable and all predictor
#' variables have approximately zero mean and unit variance, therefore it is
#' important to scale the data properly when fitting these models. For
#' convenience, this function will scale the input data automatically by
#' default, and backtransform output to the original scale. Automatic scaling
#' can be \code{"global"} (default), or \code{"local"}. The latter will scale
#' variables within (as opposed to across) populations, so is obviously only
#' applicable if there is more than 1 population. You can also scale the data
#' yourself beforehand in whatever manner you wish and set
#' \code{scaling="none"}. In this case, you will obviously have to do the
#' backtransformation yourself.
#' 
#' \strong{Predictions:}
#' 
#' There are several options for producing predictions:
#' \enumerate{
#'   \item Use \code{predictmethod="loo"} for leave-one-out prediction using the training data
#'   \item Use \code{predictmethod="sequential"} for sequential (leave-future-out) prediction using the training data
#'   \item If data frame \code{data} was supplied, supply data frame \code{datanew} containing same column names. 
#'   Column for \code{yd} is optional, unless \code{E} and \code{tau} were supplied in lieu of \code{xd}.
#'   \item If vectors/matrices were supplied for \code{yd}, \code{xd}, etc, equivalent vector/matrices \code{xnew}, 
#'   \code{popnew} (if \code{pop} was supplied), and \code{timenew} (optional).
#'   \code{ynew} is optional, unless \code{E} and \code{tau} were supplied in lieu of \code{xd}.
#' }
#' 
#' After fitting, predictions can also be made using \code{\link{predict.GP}}.
#' Predictions can plotted using \code{\link{plot.GP}}.
#' Get conditional responses using \code{\link{getconditionals}}.
#'
#' @param data A data frame, or matrix with named columns.
#' @param yd The response variable (required). If \code{data} is supplied, a column name
#'   (character) of index (numeric). If \code{data} is not supplied, a numeric vector.
#' @param xd Predictor variables. If \code{data} is supplied, column names
#'   (character vector) of indices (numeric vector). If \code{data} is not supplied, a numeric matrix 
#'   (or vector, if there is only one predictor variable). If \code{xd} is not supplied,
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
#' @param E Embedding dimension. Required if \code{xd} is not supplied.
#' @param tau Time delay. Required if \code{xd} is not supplied.
#' @param scaling How the variables should be scaled (see Details). Scaling can be \code{"global"} 
#'   (default), \code{"local"} (only applicable if more than 1 pop), or \code{"none"}. Scaling is applied 
#'   to \code{yd} and each \code{xd}. For a mix of scalings, do it manually beforehand and 
#'   set scaling to \code{"none"}. All outputs are backtransformed to original scale.
#' @param initpars Starting values for hyperparameters (see Details) in the order
#'   \code{c(phi,ve,tau,rho)}, in constrained (non-transformed) space. Should be a numeric 
#'   vector with length \code{ncol(xd)+3}. Defaults used if not supplied: \code{c(rep(0.1,ncol(xd)), 0.001, 1-0.001, 0.5)}
#' @param modeprior This value is used by the phi prior and sets the expected number of modes over the unit interval. 
#'   Defaults to 1.
#' @param rhofixed Fixes the rho parameter, if desired (see Details).
#' @inheritParams predict.GP
#'
#' @return A list (class GP) with the following elements:
#' \item{pars}{Posterior mode of hyperparameters.}
#' \item{parst}{Posterior mode of hyperparameters, transformed to real line (used for optimization).}
#' \item{grad}{Likelihood gradients of hyperparameters at posterior mode. Can be used to assess convergence.}
#' \item{covm}{List containing various covariance matrices and the inverse covariance matrix.}
#' \item{iter}{Number of iterations until convergence was reached. Currently fixed at a max of 200.}
#' \item{inputs}{Inputs and scaled inputs.}
#' \item{scaling}{Scaling information.}
#' \item{insampresults}{Data frame with in-sample predictions. \code{predfsd} is the standard
#' deviation of the GP function, \code{predsd} includes observation error.}
#' \item{insampfitstats}{Fit statistics for in-sample predictions.}
#' \item{outsampresults}{Data frame with out-of-sample predictions (if requested). \code{predfsd} is the standard
#' deviation of the GP function, \code{predsd} includes observation error.}
#' \item{outsampfitstats}{Fit statistics for out-of-sample predictions (if requested). 
#'   Only computed if using \code{"loo"} or \code{"sequential"}, if \code{yd} is found in \code{datanew},
#'   or if \code{ynew} supplied (i.e. if the observed values are known).}
#' @seealso \code{\link{predict.GP}}, \code{\link{plot.GP}}, \code{\link{getconditionals}}
#' @references Munch, S. B., Poynor, V., and Arriaza, J. L. 2017. Circumventing structural uncertainty: 
#'   a Bayesian perspective on nonlinear forecasting for ecology. Ecological Complexity, 32: 134.
#' @examples 
#' yrand=rnorm(20)
#' testgp=fitGP(yd=yrand,E=2,tau=1,predictmethod = "loo")
#' @export
#' @importFrom laGP distance
#' @import Matrix

fitGP=function(data=NULL,yd,xd=NULL,pop=NULL,time=NULL,E=NULL,tau=NULL,scaling="global",
               initpars=NULL,modeprior=1,rhofixed=NULL,
               predictmethod=NULL,datanew=NULL,xnew=NULL,popnew=NULL,timenew=NULL,ynew=NULL) {
  
  inputs=list()
  
  #store names of predictors, if available
  if(!is.null(colnames(xd))) {
    inputs$xd_names=colnames(xd)
  }
  
  #if data frame is supplied, take columns from it and store names
  if(!is.null(data)) {
    inputs$yd_names=yd
    yd=data[,yd]
    if(!is.null(xd)) { 
      inputs$xd_names=xd
      xd=data[,xd] 
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
    pop=rep(1,length(yd))
  }   
  #if time is not supplied, make a numeric index (only used in output)
  if(is.null(time)) { 
    up=unique(pop)
    time=yd*0
    for(i in 1:length(up)) {
      time[pop==up[i]]=1:length(time[pop==up[i]])
    }
  }
  #if xd is not supplied (only yd), generate xd from lags using supplied E and tau
  if(is.null(xd)) {
    xd=makelags(yd,pop,E,tau)
  }
  
  #make sure xd is a matrix, not vector or data frame
  xd=as.matrix(xd)
  d=ncol(xd) #embedding dimension, or number of predictors
  
  #rescale data
  if(scaling=="none") {
    ymeans=mean(yd,na.rm=T)
    ysds=sd(yd,na.rm=T)
    yds=yd
    xmeans=apply(xd,2,mean,na.rm=T)
    xsds=apply(xd,2,sd,na.rm=T)
    xds=xd
    #issue warnings if data are not scaled properly
    if(ymeans>0.01|ymeans<(-0.01)|ysds>1.01|ysds<(-0.99)) {
      warning("yd is not scaled properly, model may be unreliable ", 
              "mean=",round(ymeans,3)," sd=",round(ysds,3),call. = F,immediate. = T)
    }
    if(any(xmeans>0.01)|any(xmeans<(-0.01))|any(xsds>1.01)|any(xsds<(-0.99))) {
      warning("one or more xd is not scaled properly, model may be unreliable. ", 
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
  
  #if initpars is not supplied, use default starting parameters
  if(is.null(initpars)) {
    initpars=c(rep(0.1,d),0.001, 1-0.001, 0.5)
  }
  #transform parameters
  initparst=c(log(initpars[1:d]),logit(initpars[d+1]), logit(initpars[d+2]), logit(initpars[d+3]))
  
  #if only one population, set rho to 0
  if(length(unique(pop))==1) {initparst[d+3]=0}
  
  #remove missing values
  completerows=complete.cases(cbind(yds,xds))
  Y=yds[completerows]
  X=as.matrix(xds[completerows,])
  Pop=pop[completerows]
  
  #optimize model
  output=fmingrad_Rprop(initparst,Y,X,Pop,modeprior,rhofixed)
  # contains: list(parst,pars,nllpost,grad,lnL_LOO,df,mean,cov,covm,iter)
  
  #backfill missing values
  ypred<-yfsd<-ysd<-yd*NA
  ypred[completerows]=output$mean
  yfsd[completerows]=sqrt(diag(output$cov))
  ysd[completerows]=sqrt(diag(output$cov)+output$pars["ve"])
  
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
  output$inputs=c(inputs,list(yd=yd,xd=xd,pop=pop,time=time,completerows=completerows,Y=Y,X=X,Pop=Pop))
  output$scaling=list(scaling=scaling,ymeans=ymeans,ysds=ysds,xmeans=xmeans,xsds=xsds)
  output$insampresults=data.frame(timestep=time,pop=pop,predmean=ypred,predfsd=yfsd,predsd=ysd,obs=yd)
  output$insampfitstats=c(R2=1-sum((yd-ypred)^2,na.rm=T)/sum((mean(yd,na.rm=T)-yd)^2,na.rm=T),
                          rmse=sqrt(sum((yd-ypred)^2,na.rm=T)),
                          nllpost=output$nllpost,
                          lnL_LOO=output$lnL_LOO,df=output$df)
  output[c("lnL_LOO","df","nllpost","mean","cov")]=NULL #take these out of main list
  #may eventually want to exclude some more outputs
  
  if(!is.null(predictmethod)|!is.null(xnew)|!is.null(datanew)) { #generate predictions if requested
    predictresults=predict.GP(output,predictmethod,datanew,xnew,popnew,timenew,ynew) 
    output$outsampresults=predictresults$outsampresults
    if(!is.null(predictresults$outsampfitstats)) {
      output$outsampfitstats=predictresults$outsampfitstats
    }
  }
  
  class(output)="GP"
  return(output)
}

fmingrad_Rprop = function(initparst,Y,X,pop,modeprior,rhofixed) {
  
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
  output=getlikegrad(initparst,Y,X,pop,modeprior,rhofixed,D,returngradonly=T)
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
    output_new=getlikegrad(panew,Y,X,pop,modeprior,rhofixed,D,returngradonly=T)
    nllpost_new=output_new$nllpost
    grad_new=output_new$grad
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
  output_new=getlikegrad(paopt,Y,X,pop,modeprior,rhofixed,D,returngradonly=F)
  output_new$iter=iter
  
  return(output_new)
}

getlikegrad=function(parst,Y,X,pop,modeprior,rhofixed,D,returngradonly) {
  
  d=ncol(X) #embedding dimension
  
  #transform parameters from real line to constrained space
  vemin=0.01;taumin=0.01
  vemax=4.99;taumax=4.99
  rhomin=0.0001;rhomax=1
  phi=exp(parst[1:d])
  ve=(vemax-vemin)/(1+exp(-parst[d+1]))+vemin
  tau=(taumax-taumin)/(1+exp(-parst[d+2]))+taumin
  rho=(rhomax-rhomin)/(1+exp(-parst[d+3]))+rhomin
  
  if (!is.null(rhofixed)) {rho=rhofixed}
  
  #derivative for rescaled parameters wrt inputs -- for gradient calculation
  dpars=c(phi,(ve-vemin)*(1-(ve-vemin)/(vemax-vemin)),
          (tau-taumin)*(1-(tau-taumin)/(taumax-taumin)),
          (rho-rhomin)*(1-(rho-rhomin)/(rhomax-rhomin)))
  
  #specify priors 
  priors=getpriors(phi,ve,tau,rho,modeprior,taumax,vemax)
  lp=priors$lp #log prior
  dlp=priors$dlp #derivative of log prior
  
  #get covariance matrix
  covm=getcov(phi,tau,rho,X,X,pop,pop)
  #D=covm$D #distance matrices
  popsame=covm$popsame #population simularity matrix
  Cd=covm$Cd #covariance matrix
  Id=diag(ncol(Cd))
  Sigma=Cd+ve*Id #covariance matrix with observation noise
  covm$Sigma=Sigma

  #get inverse covariance
  
  icov=getcovinv(Sigma)
  iKVs=icov$iKVs
  L=icov$L
  covm$iKVs=iKVs
  
  a=iKVs%*%Y
  like=-0.5*Y%*%a-sum(log(diag(L))) #log likelihood
  
  #calculate gradient
  dl=0*dlp
  vQ=0.5*matrix2vector(a%*%t(a)-iKVs)
  
  dl[1:d]<-sapply(D,function(x) {
    dC=-multMat(x,Cd)
    vdC=matrix2vector(dC)
    0.5*innerProd(vQ,vdC) })
  
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
  
  vdC=matrix2vector(Id)
  # dl[d+1]=vQ%*%vdC
  dl[d+1]=innerProd(vQ,vdC)
  vdC=matrix2vector(Cd/tau)
  # dl[d+2]=vQ%*%vdC
  dl[d+2]=innerProd(vQ,vdC)
  dC=Cd*(1-popsame)/rho
  vdC=matrix2vector(dC)
  # dl[d+3]=vQ%*%vdC
  dl[d+3]=innerProd(vQ,vdC)
  
  # dC[[d+1]]=Id
  # dl[d+1]=vQ%*%as.vector(dC[[d+1]])
  # dC[[d+2]]=Cd/tau
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
  
  #transformed parameters
  pars=c(phi,ve,tau,rho)
  
  mpt=Cd%*%a  #in sample mean
  Cdt=Cd-Cd%*%iKVs%*%Cd  #in sample covariance
  
  #name stuff
  names(pars)=c(paste0("phi",1:length(phi)),"ve","tau","rho")
  names(parst)=c(paste0("phi",1:length(phi)),"ve","tau","rho")
  names(neglgrad)=c(paste0("phi",1:length(phi)),"ve","tau","rho")
  
  lnL_LOO=0.5*sum(log(diag(iKVs)))-0.5*sum(a^2/diag(iKVs))
  df=sum(diag((Cd%*%iKVs))) #trace
  
  out=list(pars=pars,parst=parst,
           nllpost=neglpost,grad=neglgrad,
           lnL_LOO=lnL_LOO,df=df,mean=mpt,cov=Cdt,
           covm=covm)
  return(out)
}

getpriors=function(phi,ve,tau,rho,modeprior,taumax,vemax){
  #length scale
  lam_phi=(2^(modeprior-1))^2*pi/2 #variance for gaussian - pi/2 means E(phi)=1
  lp_phi=-0.5*sum(phi^2)/lam_phi
  dlp_phi=-(phi^1)/lam_phi
  #tau
  a_tau=2
  b_tau=2 #beta
  lp_tau=(a_tau-1)*log(tau/taumax)+(b_tau-1)*log(1-tau/taumax)
  dlp_tau=(a_tau-1)/(tau)-(b_tau-1)/(taumax-tau)
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
  
  lp=lp_phi+lp_ve+lp_tau+lp_rho #log prior
  dlp=c(dlp_phi,dlp_ve,dlp_tau,dlp_rho) #derivative of log prior
  
  return(list(lp=lp,dlp=dlp))
}

getcov=function(phi,tau,rho,X,X2,pop,pop2) {
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
  Cd=tau*exp(lC0)*(popsame+rho*(1-popsame)) #covariance matrix without obs noise
  
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
#'   \item Use \code{predictmethod="loo"} for leave-one-out prediction using the training data
#'   \item Use \code{predictmethod="sequential"} for sequential (leave-future-out) prediction using the training data
#'   \item If data frame \code{data} was supplied, supply data frame \code{datanew} containing same column names. 
#'   Column for \code{yd} is optional, unless \code{E} and \code{tau} were supplied in lieu of \code{xd}.
#'   \item If vectors/matrices were supplied for \code{yd}, \code{xd}, etc, equivalent vector/matrices \code{xnew}, 
#'   \code{popnew} (if \code{pop} was supplied), and \code{timenew} (optional).
#'   \code{ynew} is optional, unless \code{E} and \code{tau} were supplied in lieu of \code{xd}.
#' }
#'
#' @param object Output from \code{fitGP}.
#' @param predictmethod \code{loo} (default) does leave-one-out prediction using
#'   the training data. \code{sequential} does sequential (leave-future-out)
#'   prediction using the training data.
#' @param datanew Data frame containing the same columns supplied in the
#'   original model.
#' @param xnew New preditor matrix or vector. Not required if \code{datanew} is supplied.
#' @param popnew New population vector. Not required if \code{datanew} is supplied.
#' @param timenew New time vector. Not required if \code{datanew} is supplied.
#' @param ynew New response vector. Optional, unless \code{E} and \code{tau} were supplied in lieu of \code{xd}. Not required if \code{datanew} is supplied.
#' @return A list with the following elements:
#' \item{outsampresults}{Data frame with out-of-sample predictions (if requested). \code{predfsd} is the standard
#' deviation of the GP function, \code{predsd} includes observation error.}
#' \item{outsampfitstats}{Fit statistics for out-of-sample predictions.
#'   Only computed if using \code{"loo"} or \code{"sequential"}, if \code{yd} is found in \code{datanew},
#'   or if \code{ynew} supplied (i.e. if the observed values are known).}
#' @export
predict.GP=function(object,predictmethod="loo",datanew=NULL,xnew=NULL,popnew=NULL,timenew=NULL,ynew=NULL) { 
  
  iKVs=object$covm$iKVs
  phi=object$pars[grepl("phi",names(object$pars))]
  tau=object$pars[names(object$pars)=="tau"]
  rho=object$pars[names(object$pars)=="rho"]
  ve=object$pars[names(object$pars)=="ve"]
  X=object$inputs$X
  Y=object$inputs$Y
  yd=object$inputs$yd
  Pop=object$inputs$Pop
  scaling=object$scaling$scaling
  ymeans=object$scaling$ymeans
  ysds=object$scaling$ysds
  xmeans=object$scaling$xmeans
  xsds=object$scaling$xsds
  
  if(!is.null(datanew)|!is.null(xnew)) {
    
    #if data frame is supplied, take columns from it
    if(!is.null(datanew)) {
      ynew=datanew[,object$inputs$yd_names]
      if(!is.null(object$inputs$xd_names)) { xnew=datanew[,object$inputs$xd_names] }
      if(!is.null(object$inputs$pop_names)) { popnew=datanew[,object$inputs$pop_names] }
      if(!is.null(object$inputs$time_names)) { timenew=datanew[,object$inputs$time_names] }
    }
    
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
    #if xd is not supplied (only yd), generate xd from lags using supplied E and tau
    if(is.null(xnew)) {
      xnew=makelags(ynew,popnew,object$inputs$E,object$inputs$tau)
    }
    
    #make sure xd is a matrix, not vector or data frame
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
    Xnew=as.matrix(xnews[completerowsnew,])
    Popnew=popnew[completerowsnew]
    
    #get covariance matrix
    covmnew=getcov(phi,tau,rho,X,Xnew,Pop,Popnew)
    Cs=covmnew$Cd #covariance matrix
    
    #get predictions
    ymean=Cs%*%(iKVs%*%Y)
    yvar=numeric(length = nrow(Xnew))
    for (i in 1:nrow(Xnew)) {
      yvar[i]=tau-Cs[i,]%*%iKVs%*%Cs[i,]
    }
    
    #backfill missing values
    ypred<-ysd<-yfsd<-ynew*NA
    ypred[completerowsnew]=ymean
    yfsd[completerowsnew]=sqrt(yvar)
    ysd[completerowsnew]=sqrt(yvar+ve)
    
  } else {
    
    if(predictmethod=="loo") {
      yd=object$inputs$yd
      Cd=object$covm$Cd
      Sigma=object$covm$Sigma
      nd=ncol(Cd)
      ymean=numeric(length=nd)
      yvar=numeric(length=nd)
      for(i in 1:nd) {
        Cdi=Cd[i,-i]
        S_noi=Sigma[-i,-i]
        Y_noi=Y[-i]
        icov_noi=getcovinv(S_noi)
        iKVs_noi=icov_noi$iKVs
        ymean[i]=Cdi%*%iKVs_noi%*%Y_noi
        yvar[i]=tau-Cdi%*%iKVs_noi%*%Cdi
      }
      
      #backfill missing values
      ypred<-ysd<-yfsd<-yd*NA
      ypred[object$inputs$completerows]=ymean
      yfsd[object$inputs$completerows]=sqrt(yvar)
      ysd[object$inputs$completerows]=sqrt(yvar+ve)
      
      popnew=object$inputs$pop
      timenew=object$inputs$time
      ynew=yd
    }
    
    if(predictmethod=="sequential") {
      
      td=object$inputs$time[object$inputs$completerows] #time index
      up=unique(Pop)
      np=length(up) #number of populations
      tgrid=unique(td) #assuming timepoints are in order
      nt=length(tgrid) #number of timepoints
      S=object$covm$Cd #tau*exp(lC0).*(popsame+rho*(1-popsame))
      
      #mint=min(td);maxt=max(td)
      #tgrid=mint:maxt #min to max time index
      #nt=maxt-mint+1
      
      ymean1=matrix(0,nrow=length(Y),ncol=nt)
      ysd1=(tau+ve)*matrix(1,nrow=length(Y),ncol=nt)
      
      for(i in 1:(nt-1)) {
        it=which(td==tgrid[i])
        I=diag(length(it))
        SS=S[it,it]+ve*I
        icov_ss=getcovinv(SS)
        iK=icov_ss$iKVs
        k=iK%*%(Y[it]-ymean1[it,i])
        ymean1[,i+1]=ymean1[,i]+S[,it]%*%k
        S=S-S[,it]%*%iK%*%S[it,]
        ysd1[,i+1]=sqrt(abs(diag(S))+ve)
      }
      ymean=rep(NA,length(Y))
      ysd2=rep((tau+ve),length(Y))
      for(i in 1:nt) {
        for(k in 1:np) {
          ind=which(td==tgrid[i] & Pop==up[k])
          if(length(ind)>0) {
            ymean[ind]=ymean1[ind,i]
            ysd2[ind]=ysd1[ind,i]
          }
        }
      }
      
      #backfill missing values
      ypred<-ysd<-yfsd<-yd*NA
      ypred[object$inputs$completerows]=ymean
      yfsd[object$inputs$completerows]=sqrt(ysd2^2-ve)
      ysd[object$inputs$completerows]=ysd2
      
      popnew=object$inputs$pop
      timenew=object$inputs$time
      ynew=yd
    }
    
  }
  
  #unscale predictions
  if(scaling=="global") {
    ypred=ypred*ysds+ymeans
    yfsd=sqrt(yfsd^2*ysds^2)
    ysd=sqrt(ysd^2*ysds^2)
  }
  if(scaling=="local") {
    up=unique(popnew)
    for(i in 1:length(up)) {
      locmean=ymeans[as.character(up[i])==names(ymeans)]
      locsd=ysds[as.character(up[i])==names(ysds)]
      ypred[popnew==up[i]]=ypred[popnew==up[i]]*locsd+locmean
      yfsd[popnew==up[i]]=sqrt(yfsd[popnew==up[i]]^2*locsd^2)
      ysd[popnew==up[i]]=sqrt(ysd[popnew==up[i]]^2*locsd^2)
    }
  }
  
  #probably need to output xnew (combine with table?, only if 1 predictor?)
  out=list(outsampresults=data.frame(timestep=timenew,pop=popnew,predmean=ypred,predfsd=yfsd,predsd=ysd))
  if(!is.null(ynew)) {
    out$outsampresults$obs=ynew
    out$outsampfitstats=c(R2=1-sum((ynew-ypred)^2,na.rm=T)/sum((mean(ynew,na.rm=T)-ynew)^2,na.rm=T),
                          rmse=sqrt(sum((ynew-ypred)^2,na.rm=T)))
    if(length(unique(popnew))>1) { #within site fit stats
      up=unique(popnew)
      np=length(up)
      R2pop<-rmsepop<-numeric(np)
      names(R2pop)=up
      names(rmsepop)=up
      for(k in 1:np) {
        ind=which(popnew==up[k])
        R2pop[k]=1-sum((ynew[ind]-ypred[ind])^2,na.rm=T)/sum((mean(ynew[ind],na.rm=T)-ynew[ind])^2,na.rm=T)
        rmsepop[k]=sqrt(sum((ynew[ind]-ypred[ind])^2,na.rm=T))
      }
      out$outsampfitstatspop=list(R2pop=R2pop,rmsepop=rmsepop)
    }
  }
  return(out)
}

logit=function(x) {
  log(x/(1-x))
}

#' Generate delay vectors
#'
#' Create a lag matrix for a time series (or suite of time series) using given E
#' and tau values.
#'
#' @param yd A vector containing a time series. Alternatively, a matrix or data
#'   frame where each column is a time series. In the case of multiple time
#'   series, lags will be generated for each variable.
#' @param pop A vector of populations (optional).
#' @param E Embedding dimension.
#' @param tau Time delay.
#' @return A matrix with named columns, the appended number indicating tau. If
#'   \code{yd} has named columns, named columns in the lag matrix will match.
#' @export
#' @examples
#' yrand <- rnorm(20)
#' site <- rep(c("a","b"),each=10)
#' dfrand <- data.frame(firstvar=rnorm(20),secondvar=rnorm(20))
#' lags1 <- makelags(yrand,E=3,tau=1)
#' lags2 <- makelags(dfrand,E=2,tau=2)
#' lags3 <- makelags(dfrand,pop=site,E=2,tau=1)
#' 
makelags=function(yd,pop=NULL,E,tau) {
  
  yd=as.matrix(yd)
  num_vars=ncol(yd)
  if(is.null(colnames(yd))) {
    colnames(yd)=paste0("var", seq_len(num_vars))
  }
  if(is.null(pop)) {
    pop=rep(1,nrow(yd))
  }
  up=unique(pop)
  
  output=matrix(NA, nrow = nrow(yd), ncol = num_vars*E)
  col_names=character(num_vars*E)
  
  for(k in 1:length(up)) { #populations
    col_index=1
    for(j in 1:num_vars) { #variables
      ind=which(pop==up[k])
      ts=yd[ind,j]
      for (i in 1:E) {
        ts=c(rep_len(NA, tau),ts[1:(length(ts) - tau)])
        output[ind,col_index]=ts
        col_names[col_index]=paste0(colnames(yd)[j],"_", i * tau)
        col_index=col_index + 1
      }
    }
  }
  colnames(output)=col_names
  return(output)
}

