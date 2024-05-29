## code to prepare simulated data

# 2 populations with slightly different theta-logistic dynamics

thetalogistic=function(N,Nb,r,theta,s) {
  y=numeric(length=N+Nb)
  y[1]=runif(1)
  for(i in 2:(N+Nb)) {
    y[i]=y[i-1]*exp(r*(1-y[i-1]^theta))*exp(rnorm(1)*s)
  }
  y=y[(Nb+1):(N+Nb)]
  return(y)
}

set.seed(20)

N=50
Nb=200
p3=thetalogistic(N,Nb,1.2,3.5,0.05)
p4=thetalogistic(N,Nb,1.2,2.5,0.05)
par(mfrow=c(1,2))
plot(p3[1:(N-1)],p3[2:N])
plot(p4[1:(N-1)],p4[2:N])

thetalog2pop=data.frame(Time=rep(1:N,2),
                       Population=rep(c("PopA","PopB"),each=50),
                       Abundance=c(p3,p4))
tlogtest=fitGP(data = thetalog2pop, yd = "Abundance", pop = "Population", E=3, tau=1, 
                 scaling = "local", predictmethod = "loo")
con=getconditionals(tlogtest)
summary(tlogtest)
plot(tlogtest)

usethis::use_data(thetalog2pop, overwrite = TRUE)

## Hastings-Powell 3 species model
library(deSolve)

parms=list(a=5.0, b=3.0, c=0.1, d=2.0, m=0.4, mu=0.01)
Y0=c(0.8, 0.1, 9)
times=seq(0,500,1)

solve.HP=function(t, y, parms) {
  a=parms$a; b=parms$b; c=parms$c; d=parms$d; m=parms$m; mu=parms$mu
  X=y[1]
  Y=y[2]
  Z=y[3]
  dX=X*(1-X)-a*X*Y/(1+b*X)
  dY=a*X*Y/(1+b*X)-c*Y*Z/(1+d*Y)-m*Y
  dZ=c*Y*Z/(1+d*Y)-mu*Z
  return(list(c(dX, dY, dZ)))
}
HPsim <- ode(y = Y0, times=times, func=solve.HP, parms, method="ode45")

plot(HPsim[,1],HPsim[,2], type="l")
plot(HPsim[,1],HPsim[,3], type="l")
plot(HPsim[,1],HPsim[,4], type="l")

HastPow3sp=as.data.frame(HPsim)
colnames(HastPow3sp)=c("Time","X","Y","Z")

usethis::use_data(HastPow3sp, overwrite = TRUE)

# Stationary and Nonstationary series
NSseries=read.csv("data-raw/round1_all.csv", stringsAsFactors = F)
usethis::use_data(NSseries, overwrite = TRUE)

# Other models (not used)

ricker=function(N,Nb,r,s) {
  y=numeric(length=N+Nb)
  y[1]=runif(1)
  for(i in 2:(N+Nb)) {
    y[i]=y[i-1]*exp(r*(1-y[i-1]))*exp(rnorm(1)*s)
  }
  y=y[(Nb+1):(N+Nb)]
  return(y)
}
