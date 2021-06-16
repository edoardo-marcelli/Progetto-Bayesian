#@@@@@@@@@@@@@@@@@@@@##
#   SMC: FILTERING    #
#@@@@@@@@@@@@@@@@@@@@##
library(ggplot2)
q025=function(x){quantile(x,0.025)}
q975=function(x){quantile(x,0.975)}


#Generate a random walk plus noise
#---------------------------------
#n = sample dim
#sigma = obs noise variance
#tau = state noise variance
#x0 = init state

dlm.sim = function(n,sigma,tau,x0){
  x    = rep(0,n)  #state process
  y    = rep(0,n)  #obs process
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

DLM<-function(data,sig2,tau2,m0,C0){
  n  = length(data)
  m  = rep(0,n)
  C  = rep(0,n)
  s = rep(0,n)
  S = rep(0,n)
  B  = rep(0,n-1)
  for (t in 1:n){
    if (t==1){
      a = m0
      R = C0 + tau2
    }else{
      a = m[t-1]
      R = C[t-1] + tau2
    }
    A = R/(R+sig2)
    m[t] = (1-A)*a + A*y[t]
    C[t] = A*sig2
    B[t-1] = C[t]/(C[t]+tau2)
  }
  s[n] = m[n]
  S[n] = C[n]
  for (t in (n-1):1){
    s[t] = (1-B[t])*m[t]+B[t]*s[t+1]
    S[t] = (1-B[t])*C[t]+B[t]^2*S[t+1]
  }
  return(list(m=m,C=C,s=s,S=S))
}


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
#This function produce sequential importance samples
#for Space state model that have linear gaussian state process
#y=data
#N=n°sample produced
#m0,c0,tau and sigma specified as above
#It returns xs,ws and ess which are respectively the filtered states
#the relative weight and the effective sample size for any period.

SISfun<-function(data,N,m0,C0,tau,sigma,smooth){
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  for(t in 1:length(data)){
    x    = rnorm(N,x,tau)                   #sample from N(x_{0},tau)
    w    = w*dnorm(data[t],x,sigma)         #update weight
    wnorm= w/sum(w)                         #normalized weight
    ESS  = 1/sum(wnorm^2)                   #effective sample size
    draw = sample(x,size=N,replace=T,prob=wnorm)
    xs = rbind(xs,draw)
    ws = rbind(ws,wnorm)
    ess =rbind(ess,ESS)
  }

  if(smooth){
    xb = matrix(0,n,N)
    for (j in 1:N){
      xb[n,j] = xs[n,j]
      for (t in (n-1):1){
        w = dnorm(xb[t+1,j],xs[t,],tau)
        w = w/sum(w)
        xb[t,j] = sample(xs[t,],size=1,prob=w)
      }
    }
    return(list(xs=xs,ws=ws,ess=ess,xb=xb))
  }else{}
  
    return(list(xs=xs,ws=ws,ess=ess))
}
#Nota: il SIS così costruito (Chopin-Papastiliopoluos) 
#ha problemi con i pesi infatti wnorm dopo un po di iteraioni 
#è un vettore (0,0,...,0,1) dopodichè da errore

#Sequential Importance Sampling with Resampling
#----------------------------------------------
SISRfun<-function(data,N,m0,C0,tau,sigma,smooth){
  if(missing(smooth)){smooth=F}
  ws  = NULL
  xs  = NULL
  ess = NULL
  x   = rnorm(N,m0,sqrt(C0))
  w   = rep(1/N,N)
  n   = length(data)
  for (t in 1:n){
    x1  = rnorm(N,x,tau)
    w   = dnorm(data[t],x1,sigma)
    w   = w/sum(w)
    ESS = 1/sum(w^2)
    x   = sample(x1,size=N,replace=T,prob=w)
    xs  = rbind(xs,x)
    ws  = rbind(ws,w/sum(w))
    ess = c(ess,ESS)
  }
  if(smooth){
    xb = matrix(0,n,N)
    for (j in 1:N){
      xb[n,j] = xs[n,j]
      for (t in (n-1):1){
        w = dnorm(xb[t+1,j],xs[t,],tau)
        w = w/sum(w)
        xb[t,j] = sample(xs[t,],size=1,prob=w)
      }
    }
    return(list(xs=xs,ws=ws,ess=ess,xb=xb))
  }else{}
  
  return(list(xs=xs,ws=ws,ess=ess))
}

#ESS-based Adaptive Resampling
#-----------------------------
g<-ESSARfun(y,1000,0,100,1,1)
s<-SISRfun(y,1000,0,100,1,1)
ESSARfun<-function(data,N,m0,C0,tau,sigma){
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  wnorm= w/sum(w)
  ESS  = 1/sum(wnorm^2)  
  for(t in 1:length(data)){
    
    if(ESS<(N/2)){
      x  =  rnorm(N,draw,tau)
      w   = dnorm(data[t],x,sigma)
      wnorm= w/sum(w)
      draw   = sample(x,size=N,replace=T,prob=wnorm)
      xs = rbind(xs,draw)
    }else{
      x    = rnorm(N,x,tau) 
      w    = w*dnorm(data[t],x,sigma)
      wnorm= w/sum(w)
      draw = sample(x,size=N,replace=T,prob=wnorm)
      xs = rbind(xs,draw)}
    
    ESS  = 1/sum(wnorm^2)  
    ws = rbind(ws,wnorm)
    ess =rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

#ESS-based adaptive resampling specificato da libro Petris
#---------------------------------------------------------
BPFfun<-function(data,N,m0,C0,tau,sigma){
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  for(t in 1:length(data)){
    
    x    = rnorm(N,x,tau) 
    w    = w*dnorm(data[t],x,sigma)
    wnorm= w/sum(w)
    ESS  = 1/sum(wnorm^2)  
    
    if(ESS<(N/2)){
      index<-sample(N,size=N,replace=T,prob=wnorm)
      x<-x[index]
      wnorm<-rep(1/N,N)
      
    }else{}
    
    draw   = sample(x,size=N,replace=T,prob=wnorm)
    
    xs = rbind(xs,draw)
    ws = rbind(ws,wnorm)
    ess =rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

#Auxiliary Particle Filter (ESS-based resampling)
#------------------------------------------------
APFfun<-function(data,N,m0,C0,tau,sigma){
  ws  = NULL
  xs  = NULL
  ess = NULL
  x   = rnorm(N,m0,sqrt(C0))
  w   = rep(1/N,N)
  for (t in 1:n){
    w0  = w*dnorm(data[t],x,sigma)
    k   = sample(1:N,size=N,replace=TRUE,prob=w0)
    x1  = rnorm(N,x[k],tau)
    lw  = dnorm(data[t],x1,sigma,log=TRUE)-dnorm(data[t],x[k],sigma,log=TRUE)
    w   = exp(lw-max(lw))
    w   = w/sum(w)
    ESS = 1/sum(w^2)
    
    if(ESS<(N/2)){
      index<-sample(N,size=N,replace=T,prob=w)
      x1<-x1[index]
      w<-rep(1/N,N)
    }else{}
    
    x   = sample(x1,size=N,replace=T,prob=w)
    xs  = rbind(xs,x)
    ws  = rbind(ws,w)
    ess = c(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

#APF no ESS based resampling
#---------------------------
NAMEfun<-function(data,N,m0,C0,tau,sigma,smooth){
  if(missing(smooth)){smooth=F}
  ws  = NULL
  xs  = NULL
  ess = NULL
  x   = rnorm(N,m0,sqrt(C0))
  w   = rep(1/N,N)
  for (t in 1:n){
    w0  = w*dnorm(data[t],x,sigma)
    k   = sample(1:N,size=N,replace=TRUE,prob=w0)
    x1  = rnorm(N,x[k],tau)
    lw  = dnorm(data[t],x1,sigma,log=TRUE)-dnorm(data[t],x[k],sigma,log=TRUE)
    w   = exp(lw-max(lw))
    w   = w/sum(w)
    ESS = 1/sum(w^2)
    x   = sample(x1,size=N,replace=T,prob=w)
    xs  = rbind(xs,x)
    ws  = rbind(ws,w/sum(w))
    ess = c(ess,ESS)
  }
  
  if(smooth){
    xb = matrix(0,n,N)
    for (j in 1:N){
      xb[n,j] = xs[n,j]
      for (t in (n-1):1){
        w = dnorm(xb[t+1,j],xs[t,],tau)
        w = w/sum(w)
        xb[t,j] = sample(xs[t,],size=1,prob=w)
      }
    }
    return(list(xs=xs,ws=ws,ess=ess,xb=xb))
  }else{}
  
  return(list(xs=xs,ws=ws,ess=ess))
}

#.................................................
#Filtervalues is a function that computes the mean
#and the variance of a SISR already estimated (post estimation command)
Filtervalues<-function(fun){
  xhat<-sapply(1:n,function(i)
    weighted.mean(fun$xs[i,],fun$ws[i,]))
  sdhat<-sapply(1:n,function(i)
    sqrt(weighted.mean((fun$xs[i,]-xhat[i])^2,fun$ws[i,])))
  
  return(list(mean=xhat,sd=sdhat))
  
}
#Example of Comparison (do in GGPLOT if possible)
h<-SISRfun(y,N,m0,C0,tau,sigma)
Var<-Filtervalues(h)$var
plot(sqrt(C),type="l",ylim=c(0.4,1))
par(new=T)
plot(Var,type="l",col="red",ylim=c(0.4,1))

#Filterplot is a function similar to Filtervalues 
#returning plots of median and quantiles
Filterplot<-function(data,fun,title){
  require(ggplot2)
  mx = apply(fun$xs,1,median)
  lx = apply(fun$xs,1,q025)
  ux = apply(fun$xs,1,q975)
  
  timeframe<-c(1:length(data))
  df<-data.frame(timeframe,data,mx,lx,ux)
  
  ggplot(df,aes(x=timeframe))+
    geom_line(aes(y=data))+
    geom_line(aes(y=mx),col="red")+
    geom_ribbon(aes(ymin = lx, ymax = ux),
                fill="red",alpha=0.16) +
    labs(x="Time",
         y="")+
    ggtitle(title)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
}

Smoothplot<-function(data,fun,title){
    require(ggplot2)
    mx = apply(fun$xb,1,median)
    lx = apply(fun$xb,1,q025)
    ux = apply(fun$xb,1,q975)
    
    timeframe<-c(1:length(data))
    df<-data.frame(timeframe,data,mx,lx,ux)
    
    ggplot(df,aes(x=timeframe))+
      geom_line(aes(y=data))+
      geom_line(aes(y=mx),col="blue")+
      geom_ribbon(aes(ymin = lx, ymax = ux),
                  fill="red",alpha=0.16) +
      labs(x="Time",
           y="")+
      ggtitle(title)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
  
}


