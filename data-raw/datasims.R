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
