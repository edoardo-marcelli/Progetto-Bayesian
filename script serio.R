#@@@@@@@@@@@@@@@@@@@@##
#SEQUENTIAL MONTECARLO#
#@@@@@@@@@@@@@@@@@@@@##
library(ggplot2)

#Generate a random walk plus noise
#---------------------------------
#n = sample dim
#v = obs noise
#w = state noise
#x0 = init state

dlm.sim = function(n,sigma,tau,x0){
  x    = rep(0,n)
  y        = rep(0,n)
  x[1] = rnorm(1,x0,tau)
  y[1] = rnorm(1,x[1],sigma)
  for (t in 2:n){
    x[t] = rnorm(1,x[t-1],tau)
    y[t] = rnorm(1,x[t],sigma)
  }
  return(list(x=x,y=y))
}

#Filtering, using data
#---------------------
#y = observed data

DLM = function(y,sigma2,tau2,m0,C0){
  n = length(y)
  m = rep(0,n)
  C = rep(0,n)
  for (t in 1:n){
    if (t==1){
      a = m0
      R = C0 + tau2
    }else{
      a = m[t-1]
      R = C[t-1] + tau2
    }
    A = R/(R+sigma2)
    m[t] = (1-A)*a + A*y[t]
    C[t] = A*sigma2
  }
  return(list(m=m,C=C))
}


########################################################

#Generate a random walk plus noise
#---------------------------------
n=500
tau2=1
sigma2=1
sigma=sqrt(sigma2)
tau=sqrt(tau2)
x0=0

sim<-dlm.sim(n,sigma,tau,x0)
y<-sim$y
x<-sim$x

#Kalman Filter
#-------------
m0=0
C0=100
filtval<-DLM(y,sigma2,tau2,m0,C0)
m<-filtval$m
C<-filtval$C
Low = m + qnorm(0.025)*sqrt(C)
Up = m + qnorm(0.975)*sqrt(C)

#Put everything in a data frate (for ggplot)
timeframe<-c(1:length(y))
DLM.df<-data.frame(timeframe,y,x,m,C,Low,Up)

#Plots (observetions in black and filtered states in red)
ggplot(DLM.df,aes(x=timeframe))+
  geom_line(aes(y=y))+
  geom_line(aes(y=m),col="red")+
  geom_ribbon(aes(ymin = Low, ymax = Up),
              fill="red",alpha=0.16) +
  labs(x="Time",
       y="")+
  ggtitle("Kalman Filter")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

#Sequential Importance Sampling
#------------------------------
#N = number of samples
N   = 1000
#matrix of {(x,w)_i} where i=1 to N. For n times
x<-array(0,c(n,N))
w<-array(0,c(n,N))
wnorm<-array(0,c(n,N))
ESS<-c()

#first sample {(x0,w0)_i}
  x0   = rnorm(N,m0,sqrt(C0))
  w0   = rep(1/N,N)
#second sample{(x1,w1)_i}
  x[1,] = rnorm(N,x0,tau)         #sample from N(x_{0},tau)
  w[1,] = w0*dnorm(y[1],x[1,],sigma)  #update weight
  wnorm[1,]   = w[1,]/sum(w[1,])               #normalized weight
  ESS[1]  = 1/sum(wnorm[1,]^2)
#other n-1 samples {(xt,wt)_i}
for(t in c(2,n)){
  
  x[t,]    = rnorm(N,x[t-1,],tau)         #sample from N(x_{0},tau)
  w[t,]    = w[t-1,]*dnorm(y[t],x[t,],sigma)  #update weight
  wnorm[t,]   = w[t,]/sum(w[t,])               #normalized weight
  ESS[t]  = 1/sum(wnorm[t,]^2)            #effective sample size
  
  }
#Nota: il SIS così costruito (Chopin-Papastiliopoluos) non produce ESS
  

#Sequential Importance Sampling with Resampling
#----------------------------------------------
#N = number of samples
N   = 1000
#matrix of {(x,w)_i} where i=1 to N. For n times
x<-array(0,c(n,N))
w<-array(0,c(n,N))
wnorm<-array(0,c(n,N))
ESS<-c()
  
  #first sample {(x0,w0)_i}
  x0   = rnorm(N,m0,sqrt(C0))
  w0   = rep(1/N,N)
  #second sample{(x1,w1)_i}
  x[1,] = rnorm(N,x0,tau)         #sample from N(x_{0},tau)
  w[1,] = w0*dnorm(y[1],x[1,],sigma)  #update weight
  wnorm[1,]   = w[1,]/sum(w[1,])               #normalized weight
  ESS[1]  = 1/sum(wnorm[1,]^2)
  #other n-1 samples {(xt,wt)_i}
  for(t in c(2,n)){
    
    x[t,]    = rnorm(N,x[t-1,],tau)         #sample from N(x_{0},tau)
    w[t,]    = w[t-1,]*dnorm(y[t],x[t,],sigma)  #update weight
    wnorm[t,]   = w[t,]/sum(w[t,])               #normalized weight
    ESS[t]  = 1/sum(wnorm[t,]^2)            #effective sample size
    
  }

  
  
 


