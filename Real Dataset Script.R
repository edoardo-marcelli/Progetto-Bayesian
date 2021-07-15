#############################
#
#   y(t) = exp(x(t)/2)*v(t)
#   x(t) = alpha + beta*x(t-1) + tau*w(t)
#
#   v(t) ~ N(0,1)
#   w(t) ~ N(0,1)
#   x(0) ~ N(0,100)
#
#   For LW, we have a prior on alpha,beta,tau2
#
#   alpha ~ N(ealpha,valpha)
#   beta  ~ N(ebeta,vbeta)
#   tau2  ~ IG(nu/2,nu*lambda/2)
#   
#
#Sequential Importance Sampling
#---------------
SISfun<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    x<-rnorm(N,alpha+beta*x,tau)
    w1<-w*dnorm(data[t],0,exp(x/2))
    
    w = w1/sum(w1)
    ESS  = 1/sum(w^2)
    
    if(ESS<0){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x[index]
      w<-rep(1/N,N)
    }else{}
    
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

# Basis Particle Filter (Resampling at each step)
#---------------
SVBAPFfun<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    x<-rnorm(N,alpha+beta*x,tau)
    w1<-w*dnorm(data[t],0,exp(x/2))
    
    w = w1/sum(w1)
    ESS  = 1/sum(w^2)
    
    if(ESS>-1){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x[index]
      w<-rep(1/N,N)
    }else{}
    
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

# Bootstrap Particle Filter
#---------------
SVPFfun<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    x<-rnorm(N,alpha+beta*x,tau)
    w1<-w*dnorm(data[t],0,exp(x/2))
    
    w = w1/sum(w1)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x[index]
      w<-rep(1/N,N)
    }else{}
    
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

#
#-----------------------
#Guided Particle Filter Optimal Kernel
#------------------------------
SVPFoptfun<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    xprev<-x
    x<-rnorm(N,alpha+beta*x+(tau^2)/4*((data[t]^2)*exp(-(alpha+beta*x))-2),tau)
    w1<-w*dnorm(data[t],0,exp(x/2))*dnorm(x,alpha+beta*xprev,tau)/dnorm(x,mean=alpha+beta*xprev+(tau^2)/4*((data[t]^2)*exp(-(alpha+beta*xprev))-2),tau)
    
    w = w1/sum(w1)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x[index]
      w<-rep(1/N,N)
    }else{}
    
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

#Auxiliary Particle Filter
#-------------------------
SVAPFfun<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    weight = w*dnorm(data[t],0,exp(x/2))
    k   = sample(1:N,size=N,replace=TRUE,prob=weight)
    x1   = rnorm(N,alpha+beta*x[k],tau)
    lw  = dnorm(data[t],0,exp(x1/2),log=TRUE)-dnorm(data[t],0,exp(x[k]/2),log=TRUE)
    w   = exp(lw)
    w   = w/sum(w)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x1[index]
      w<-rep(1/N,N)
    }else{}
    
    x <- x1
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
    
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

#Auxiliary Particle Filter (PETRONE p 218)
#-------------------------
SVAPFfun_2<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    x_prev=x
    x_hat = alpha+beta*x #means of the transition from the previous particle 
    I_k = sample(1:N,prob=w*dnorm(data[t],0,sqrt(exp((x_hat)/2))))
    for(k in 1:N){
      x[k]=rnorm(size=1,mean=alpha+beta*x_prev[I_k],sd=tau)
      w[k]=dnorm(data[t],0,sqrt(exp((x[k])/2)))/dnorm(data[t],0,sqrt(exp((x_hat)/2)))
    }
    w = w/sum(w)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x[index]
      w<-rep(1/N,N)
    }else{}
    
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
    
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

#Auxiliary Particle Filter Optimal Kernel
#----------------------------------------
SVAPFoptfun<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    xprev<-x
    
    weight = w*dnorm(data[t],0,exp(x/2))
    k   = sample(1:N,size=N,replace=TRUE,prob=weight)
    x1 <-rnorm(N,alpha+beta*x[k]+(tau^2)/4*((data[t]^2)*exp(-(alpha+beta*x[k]))-2),tau)
    lw  = dnorm(data[t],0,exp(x1/2),log=TRUE)-dnorm(data[t],0,exp(x[k]/2),log=TRUE)
    w   = exp(lw)
    w   = w/sum(w)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x1[index]
      w<-rep(1/N,N)
    }else{}
    
    w1<-w*dnorm(data[t],0,exp(x/2))*dnorm(x,alpha+beta*xprev,tau)/dnorm(x,mean=alpha+beta*xprev+(tau^2)/4*((data[t]^2)*exp(-(alpha+beta*xprev))-2),tau)
    x <- x1
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
    
  }
  return(list(xs=xs,ws=ws,ess=ess))
}


#Liu West
#--------
SVLWfun<-function(data,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda){
  
  xs = rnorm(N,m0,sqrt(C0))
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
  vpars   = NULL
  for (t in 1:length(data)){
    
    #mpar        = apply(pars,2,mean) #alternatively use this 
    mpar<-c()
    for(i in 1:3){
      mpar[i]    = weighted.mean(pars[,i],w)
    }
    
    vpar     = var(pars)
    ms          = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    mus         = pars[,1]+pars[,2]*xs 
    weight      = w*dnorm(data[t],0,exp(mus/2))
    k           = sample(1:N,size=N,replace=T,prob=weight)
    ms1         = ms[k,] + matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    xt          = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    w           = dnorm(data[t],0,exp(xt/2))/dnorm(data[t],0,exp(mus[k]/2))
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
    
    vpars        = rbind(vpars,vpar)
    xss         = rbind(xss,xs)
    parss[,,t]  = pars 
    ws          = rbind(ws,w)
  }
  return(list(xs=xss,pars=parss,ws=ws,vpars=vpars))
}

#Set hyperparameters
#-------------------
#n     =  length(y)
alpha = -0.0031
beta  =  0.9951
tau2  =  0.0074
m0      = 0
C0      = 100
sC0     = sqrt(C0)
ealpha  = alpha
valpha  = 0.01
ebeta    = beta
vbeta    = 0.01
nu      = 3
lambda  = tau2
tau   = sqrt(tau2)

#Real dataset

library(openxlsx)
realdatasetSPX<-read.xlsx("Dataset 17-21.xlsx", sheet = 1)
y<-realdatasetSPX[,2]
n = length(y)
volatilitySPX<-read.xlsx("Dataset 17-21.xlsx", sheet = 1)
x<-volatilitySPX[,3]
alpha.true = alpha
beta.true  = beta
tau2.true  = tau2
timeframe<-c(1:length(y))
dfsv<-data.frame(timeframe,y,x)


#Filtering
#---------
set.seed(12345)
N=10000
svsis<-SISfun(y,N,m0,C0,alpha,beta,tau)
svbapf<-SVBAPFfun(y,N,m0,C0,alpha,beta,tau)
svpf<-SVPFfun(y,N,m0,C0,alpha,beta,tau)
svapf<-SVAPFfun(y,N,m0,C0,alpha,beta,tau)
svlw<-SVLWfun(y,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda)
svopt<-SVPFoptfun(y,N,m0,C0,alpha,beta,tau)
svapfopt<-SVAPFoptfun(y,N,m0,C0,alpha,beta,tau)

#Filtering Values and Plot
#-------------------------
Filtervalues<-function(fun){
  xhat<-sapply(1:n,function(i)
    weighted.mean(fun$xs[i,],fun$ws[i,]))
  sdhat<-sapply(1:n,function(i)
    sqrt(weighted.mean((fun$xs[i,]-xhat[i])^2,fun$ws[i,])))
  
  return(list(mean=xhat,sd=sdhat))
  
}

Filtplot<-function(dataframe,fun,title){
  
  Filt<-Filtervalues(fun)
  mean<-Filt$mean
  mean1=exp(mean/2)
  sd<-Filt$sd
  dataframe<-data.frame(dataframe,mean1,sd)
  
  ggplot(dataframe,aes(x=timeframe))+
    geom_line(aes(y=x, col="True States", linetype="True States"))+
    geom_line(aes(y=mean1, col="Filtered States", linetype="Filtered States"))+
    geom_ribbon(aes(ymin = mean1-1.96*sd, ymax = mean1+1.96*sd),
                fill="red",alpha=0.16) +
    scale_color_manual("",
                       values=c("True States" = "black",
                                "Filtered States"= "red"))+
    scale_linetype_manual("",
                          values=c("True States" = 3,
                                   "Filtered States"=1))+
    labs(x="Time",
         y="")+
    ggtitle(title)+
    theme_bw()+
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) +
    theme(plot.title = element_text(hjust = 0.5))
}

#Plots
#-----
library(ggplot2)
library(ggpubr)
plot1<-Filtplot(dfsv,svpf,"Particle Filter")
plot2<-Filtplot(dfsv,svapf, "Auxiliary Particle Filter")
plot3<-Filtplot(dfsv,svlw,"Liu and West Filter")
plot4<-Filtplot(dfsv,svopt,"Opt Ker Particle Filter")
plot5<-Filtplot(dfsv,svapfopt,"Opt Ker Auxiliary Particle Filter")
plot6<-Filtplot(dfsv,svsis,"No Resampling")
plot7<-Filtplot(dfsv,svbapf,"Always Resampling")
ggarrange(plot1,plot2)
ggarrange(plot1)
ggarrange(plot2)
ggarrange(plot3)
ggarrange(plot4)
ggarrange(plot5)
ggarrange(plot6,plot1)
ggarrange(plot7)

## RMSE MAE comparison

realisedx<-x
Errorvol<-matrix(NA,ncol=6,nrow=2)
colnames(Errorvol)<-c("N","BPF", "GPFOPT", "APF","LWF", "SIS", "PF")
rownames(Errorvol)<-c("RMSE","MAE")
Errorvol[,1]<-c(10000)
RMSE<-function(x,xhat){sqrt(mean((x-xhat)^2))}
MAE<-function(x,xhat){mean(abs(x-xhat))}
comparablevol<-function(fun){exp((Filtervalues(fun)$mean)/2)}
Errorvol[1,2]<-RMSE(realisedx,comparablevol(svpf))
Errorvol[1,3]<-RMSE(realisedx,comparablevol(svopt))
Errorvol[1,4]<-RMSE(realisedx,comparablevol(svapf))
Errorvol[1,5]<-RMSE(realisedx,comparablevol(svlw))
Errorvol[1,6]<-RMSE(realisedx,comparablevol(svsis))
Errorvol[1,7]<-RMSE(realisedx,comparablevol(svbapf))
Errorvol[2,2]<-MAE(realisedx,comparablevol(svpf))
Errorvol[2,3]<-MAE(realisedx,comparablevol(svopt))
Errorvol[2,4]<-MAE(realisedx,comparablevol(svapf))
Errorvol[2,5]<-MAE(realisedx,comparablevol(svlw))
Errorvol[2,6]<-MAE(realisedx,comparablevol(svsis))
Errorvol[2,7]<-MAE(realisedx,comparablevol(svbapf))

#Kalman filter
m00=0
c00=100
nu2=1
sig2=1
DLM<-function(data,sig2,nu2,m00,C00){
  n  = length(data)
  m  = rep(0,n)
  C  = rep(0,n)
  f  = rep(0,n)
  Q  = rep(0,n)
  for (t in 1:n){
    if (t==1){
      a = m00
      R = C00 + nu2
      j = a
      K = R + sig2
    }else{
      a = m[t-1]
      R = C[t-1] + nu2
      j = a
      K = R + sig2
    }
    A = R/(R+sig2)
    m[t] = (1-A)*a + A*y[t]
    C[t] = A*sig2
    f[t] = j
    Q[t] = K
  }
  return(list(m=m,C=C,f=f,Q=Q))
}
filtval<-DLM(y,sig2,nu2,m00,c00)
f<-filtval$f
Q<-filtval$Q
Low = f + qnorm(0.025)*sqrt(Q)
Up = f + qnorm(0.975)*sqrt(Q)

#Put everything in a data frate (for ggplot)
timeframeKF<-c(1:length(y))
DLM.df<-data.frame(timeframeKF,y,x,f,Q,Low,Up)

library(ggplot2)
library(ggpubr)
ggplot(DLM.df,aes(x=timeframeKF))+
  geom_line(aes(y=y, col="Observations", linetype="Observations"))+
  geom_line(aes(y=x, col="True States", linetype="True States"))+
  geom_line(aes(y=f, col="Filtered States", linetype="Filtered States"))+
  geom_ribbon(aes(ymin = Low, ymax= Up),
              fill="red",alpha=0.16) +
  scale_color_manual("",
                     values=c("Observations" = "black",
                              "True States" = "black",
                              "Filtered States"= "red"))+
  scale_linetype_manual("",
                        values=c("Observations" = 1,
                                 "True States" = 3,
                                 "Filtered States"=1))+
  labs(x="Time",
       y="")+
  ggtitle("Kalman Filter")+
  theme_bw()+
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
        legend.background = element_rect(fill = "white", colour = "gray30")) +
  theme(plot.title = element_text(hjust = 0.5))



