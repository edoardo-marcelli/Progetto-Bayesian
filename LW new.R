############################
#
# STOCHASTIC VOLATILITY
#
# LW FILTER
#
#############################
#
#   y(t) = exp(h(t)/2)*v(t)
#   x(t) = alpha + beta*x(t-1) + tau*w(t)
#
#   v(t) ~ N(0,1)
#   w(t) ~ N(0,1)
#   x(0) ~ N(0,100)
#
#   PRIOR ON alpha,beta,tau2

# Simulating the data
# -------------------
set.seed(12345)
n     =  500
alpha = -0.0031
beta  =  0.9951
tau2  =  0.0074
tau   = sqrt(tau2)
y     = rep(0,n)
x     = rep(0,n)
x[1]  = alpha/(1-beta)
y[1]  = rnorm(1,0,exp(x[1]/2))
for (t in 2:n){
  x[t] = alpha + beta*x[t-1] + rnorm(1,0,tau)
  y[t] = rnorm(1,0,exp(x[t]/2))
}
alpha.true = alpha
beta.true  = beta
tau2.true  = tau2

# Data and prior hyperparameters 
# ------------------------------
m0      = 0.0
C0      = 0.1
sC0     = sqrt(C0)
ealpha  = alpha
valpha  = 0.01
ebeta    = beta
vbeta    = 0.01
nu      = 3
lambda  = tau2

# Liu and West filter
# -------------------
set.seed(8642)
N      = 50000

LWfun<-function(data,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda){

xs     = rnorm(N,m0,sqrt(C0))
pars   = cbind(rnorm(N,ealpha,sqrt(valpha)),rnorm(N,ebeta,sqrt(vbeta)),
               log(1/rgamma(N,nu/2,nu*lambda/2)))
delta  = 0.75
a      = (3*delta-1)/(2*delta)
h2     = 1-a^2
parss  = array(0,c(N,3,n))
xss    = NULL
ws     = NULL
ESS    = NULL
w      = rep(1/N,N)
for (t in 1:length(data)){
  mpar        = apply(pars,2,mean)
  vpar        = var(pars)
  ms          = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
  mus         = pars[,1]+pars[,2]*xs 
  weight      = w*dnorm(data[t],0,exp(mus/2))
  k           = sample(1:N,size=N,replace=T,prob=weight)
  ms1         = ms[k,] + matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
  xt          = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
  w           = dnorm(y[t],0,exp(xt/2))/dnorm(y[t],0,exp(mus/2))
  w           = w/sum(w)
  ESS         = 1/sum(w^2)
  
  
  if(ESS<(N/2)){
    index<-sample(N,size=N,replace=T,prob=w)
    xs<-xt[index]
    pars<-ms1[index,]
    w<-rep(1/N,N)
  }else{
    xs<-xt
    pars<-ms1
  }
  
  xss         = rbind(xss,xs)
  parss[,,t]  = pars 
  ws          = rbind(ws,w)
}
return(list(xss=xss,parss=parss,ws=ws))
}

l<-LWfun(y,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda)
