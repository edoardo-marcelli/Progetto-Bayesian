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
BAPFfun<-function(data,N,m0,C0,alpha,beta,tau,r){
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
BPFfun<-function(data,N,m0,C0,alpha,beta,tau,r){
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
GPFoptfun<-function(data,N,m0,C0,alpha,beta,tau,r){
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
APFfun<-function(data,N,m0,C0,alpha,beta,tau,r){
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
APFoptfun<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      epsilon=alpha+beta*x+((tau^2)/4)*((data[t]^2)*exp(-alpha-beta*x)-2)
      weight = w*exp((1/(2*tau^2))*((epsilon)^2-(alpha+beta*x)^2)-((data[t]^2)/2)*exp(-alpha-beta*x)*(1+alpha+beta*x))
      k   = sample(1:N,size=N,replace=TRUE,prob=weight)
      x1   = rnorm(N,alpha+beta*x[k]+(tau^2)/4*((data[t]^2)*exp(-(alpha+beta*x[k]))-2),tau)
      w   = rep(1/N,N)
    }else{
      weight = rep(1/N,N)
      k   = sample(1:N,size=N,replace=FALSE,prob=weight)
      x1   = rnorm(N,alpha+beta*x[k]+(tau^2)/4*((data[t]^2)*exp(-(alpha+beta*x[k]))-2),tau)
      w   <-  w*dnorm(data[t],0,exp(x1/2))*dnorm(x1,alpha+beta*x[k],tau)/
        dnorm(x1,mean=alpha+beta*x[k]+(tau^2)/4*((data[t]^2)*exp(-(alpha+beta*x[k]))-2),tau)
    }
    
    w   = w/sum(w)
    x <- x1
    
    
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
    
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

#Liu West
#--------
LWfun<-function(data,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda){
  
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
alpha =  0
beta  =  0.99
tau2  =  0.05
m0      = 0
C0      = 100 #N(0,100) is an informative prior about the mean of the volatility. It is used by Jacquier et al. (2004), and it is calibrated on percentage returns
sC0     = sqrt(C0)
ealpha  = alpha
valpha  = 0.01
ebeta    = beta
vbeta    = 0.01
nu      = 5
lambda  = tau2*3/5
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
set.seed(12345)
N=10000
svbapf<-BAPFfun(y,N,m0,C0,alpha,beta,tau)
set.seed(12345)
N=10000
svbpf<-BPFfun(y,N,m0,C0,alpha,beta,tau)
set.seed(12345)
N=10000
svapf<-APFfun(y,N,m0,C0,alpha,beta,tau)
set.seed(12345)
N=10000
svgpfopt<-GPFoptfun(y,N,m0,C0,alpha,beta,tau)
set.seed(12345)
N=10000
svlw<-LWfun(y,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda)
set.seed(12345)
N=10000
svapfopt<-APFoptfun(y,N,m0,C0,alpha,beta,tau)

#Filtering Values
#-------------------------
Filtervalues<-function(fun){
  xhat<-sapply(1:n,function(i)
    weighted.mean(fun$xs[i,],fun$ws[i,]))
  sdhat<-sapply(1:n,function(i)
    sqrt(weighted.mean((fun$xs[i,]-xhat[i])^2,fun$ws[i,])))
  
  return(list(mean=xhat,sd=sdhat))
  
}

#Graphical comparison of filtered volatility

Filtplot<-function(dataframe,fun,title){
  
  Filt<-Filtervalues(fun)
  mean<-Filt$mean
  mean1=exp(mean/2)
  sd<-Filt$sd ##qua c'è il problema che la sd si riferisca allo stato non trasformato...
  ## un'idea potrebbe essere di trasformare prima tutte le particles, calcolando solo poi media e standard deviation...ma questo porta a una diversa media
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
plot1<-Filtplot(dfsv,svbpf,"Bootstrap Particle Filter")
plot2<-Filtplot(dfsv,svapf, "Auxiliary Particle Filter")
plot3<-Filtplot(dfsv,svlw,"Liu and West Filter")
plot4<-Filtplot(dfsv,svgpfopt,"Opt Ker Guided Particle Filter")
plot5<-Filtplot(dfsv,svsis,"No Resampling")
plot6<-Filtplot(dfsv,svbapf,"Always Resampling")
plot7<-Filtplot(dfsv,svapfopt,"Opt Ker Auxiliary Particle Filter")
ggarrange(plot1,plot2)
ggarrange(plot3,plot4)
ggarrange(plot5,plot6)
plot1
plot2
plot3
plot4
plot5
plot6
plot7

## RMSE MAE comparison

realisedx<-x
Errorvol<-matrix(NA,ncol=7,nrow=2)
colnames(Errorvol)<-c("N","BPF", "GPFOPT", "APF","LWF", "SIS", "BAPF")
rownames(Errorvol)<-c("RMSE","MAE")
Errorvol[,1]<-c(5000)
RMSE<-function(x,xhat){sqrt(mean((x-xhat)^2))}
MAE<-function(x,xhat){mean(abs(x-xhat))}
comparablevol<-function(fun){exp((Filtervalues(fun)$mean)/2)}
Errorvol[1,2]<-RMSE(realisedx,comparablevol(svbpf))
Errorvol[1,3]<-RMSE(realisedx,comparablevol(svgpfopt))
Errorvol[1,4]<-RMSE(realisedx,comparablevol(svapf))
Errorvol[1,5]<-RMSE(realisedx,comparablevol(svlw))
Errorvol[1,6]<-RMSE(realisedx,comparablevol(svsis))
Errorvol[1,7]<-RMSE(realisedx,comparablevol(svbapf))
Errorvol[2,2]<-MAE(realisedx,comparablevol(svbpf))
Errorvol[2,3]<-MAE(realisedx,comparablevol(svgpfopt))
Errorvol[2,4]<-MAE(realisedx,comparablevol(svapf))
Errorvol[2,5]<-MAE(realisedx,comparablevol(svlw))
Errorvol[2,6]<-MAE(realisedx,comparablevol(svsis))
Errorvol[2,7]<-MAE(realisedx,comparablevol(svbapf))

# Save mean and sd series of PF with 50000 particles
library(openxlsx)
Filt<-Filtervalues(svbpf)
mean<-Filt$mean
sd<-Filt$sd
list1<-list('mean'=mean,'sd'=sd)
write.xlsx(list1,file="svbpf50000p.xlsx")

#RMSE and MAE comparison of particle filters (but not LW) vs Bootstrap PF with 50000 particles

library(openxlsx)
RMSE<-function(x,xhat){sqrt(mean((x-xhat)^2))}
MAE<-function(x,xhat){mean(abs(x-xhat))}
bpf50000mean<-read.xlsx("svbpf50000p.xlsx", sheet = 1)
bpf50000sd<-read.xlsx("svbpf50000p.xlsx", sheet = 2)
bpf100000mean<-read.xlsx("svbpf100000p.xlsx", sheet = 1)
bpf100000sd<-read.xlsx("svbpf100000p.xlsx", sheet = 2)
apf50000mean<-read.xlsx("svapf50000p.xlsx", sheet = 1)
apf50000sd<-read.xlsx("svapf50000p.xlsx", sheet = 2)
apf100000mean<-read.xlsx("svapf100000p.xlsx", sheet = 1)
apf100000sd<-read.xlsx("svapf100000p.xlsx", sheet = 2)
z50mean<-bpf50000mean[,1]
z50sd<-bpf50000sd[,1]
z100mean<-bpf100000mean[,1]
z100sd<-bpf100000sd[,1]
w50mean<-apf50000mean[,1]
w100mean<-apf100000mean[,1]
w50sd<-apf50000sd[,1]
w100sd<-apf100000sd[,1]
dfsv1<-data.frame(timeframe,y,z50mean,z50sd,w50mean,z50sd)
Errorcomp<-matrix(NA,ncol=6,nrow=2)
colnames(Errorcomp)<-c("N","BPF", "GPFOPT", "APF", "SIS", "BAPF")
rownames(Errorcomp)<-c("RMSE","MAE")
Errorcomp[,1]<-c(10000)
RMSE<-function(x,xhat){sqrt(mean((x-xhat)^2))}
MAE<-function(x,xhat){mean(abs(x-xhat))}
comparablemean<-function(fun){Filtervalues(fun)$mean}
Errorcomp[1,2]<-RMSE(z50mean,comparablemean(svbpf))
Errorcomp[1,3]<-RMSE(z50mean,comparablemean(svgpfopt))
Errorcomp[1,4]<-RMSE(z50mean,comparablemean(svapf))
Errorcomp[1,5]<-RMSE(z50mean,comparablemean(svsis))
Errorcomp[1,6]<-RMSE(z50mean,comparablemean(svbapf))
Errorcomp[2,2]<-MAE(z50mean,comparablemean(svbpf))
Errorcomp[2,3]<-MAE(z50mean,comparablemean(svgpfopt))
Errorcomp[2,4]<-MAE(z50mean,comparablemean(svapf))
Errorcomp[2,5]<-MAE(z50mean,comparablemean(svsis))
Errorcomp[2,6]<-MAE(z50mean,comparablemean(svbapf))
RMSEbench50000p<-RMSE(z50mean,w50mean)
MAEDbench50000p<-MAE(z50mean,w50mean)
difference50000p<-c(RMSEbench50000p,MAEDbench50000p)
RMSEbench100000p<-RMSE(z100mean,w100mean)
MAEDbench100000p<-MAE(z100mean,w100mean)
difference100000p<-c(RMSEbench100000p,MAEDbench100000p)

# Graphical comparison of filtered states

Filtplot1<-function(dataframe,fun,title){
  
  Filt<-Filtervalues(fun)
  mean<-Filt$mean
  sd<-Filt$sd 
  dataframe<-data.frame(dataframe,mean,sd)
  
  ggplot(dataframe,aes(x=timeframe))+
    geom_line(aes(y=z50mean, col="50000p BT Filtered States", linetype="50000p BT Filtered States"))+
    geom_line(aes(y=mean, col="Filtered States", linetype="Filtered States"))+
    geom_ribbon(aes(ymin = z50mean -1.96*z50sd, ymax = z50mean +1.96*z50sd),
                fill="black",alpha=0.16) +
    geom_ribbon(aes(ymin = mean-1.96*sd, ymax = mean+1.96*sd),
                fill="red",alpha=0.16) +
    scale_color_manual("",
                       values=c("50000p BT Filtered States" = "black",
                                "Filtered States"= "red"))+
    scale_linetype_manual("",
                          values=c("50000p BT Filtered States" = 1,
                                   "Filtered States"=1))+
    labs(x="Time",
         y="")+
    ggtitle(title)+
    theme_bw()+
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) +
    theme(plot.title = element_text(hjust = 0.5))
}

library(ggplot2)
library(ggpubr)
plot8<-Filtplot1(dfsv1,svbpf,"Bootstrap Particle Filter")
plot9<-Filtplot1(dfsv1,svapf,"Auxiliary Particle Filter")
plot10<-Filtplot1(dfsv1,svgpfopt,"Opt Ker Guided Particle Filter")
plot11<-Filtplot1(dfsv1,svsis,"No Resampling")
plot12<-Filtplot1(dfsv1,svbapf,"Always Resampling")
plot8
plot9
plot10
plot11
plot12

# Plot the two 50000 particles benchmarks
Filtplot2<-function(dataframe,title){
  
  ggplot(dataframe,aes(x=timeframe))+
    geom_line(aes(y=z50mean, col="50000p BT Filtered States", linetype="50000p BT Filtered States"))+
    geom_ribbon(aes(ymin = z50mean-1.96*z50sd, ymax = z50mean+1.96*z50sd),
                fill="black",alpha=0.16) +
    geom_line(aes(y=w50mean, col="50000p AP Filtered States", linetype="50000p AP Filtered States"))+
    geom_ribbon(aes(ymin = w50mean-1.96*w50sd, ymax = w50mean+1.96*w50sd),
                fill="red",alpha=0.16) +
    scale_color_manual("",
                       values=c("50000p BT Filtered States" = "black",
                                "50000p AP Filtered States"= "red"))+
    scale_linetype_manual("",
                          values=c("50000p BT Filtered States" = 1,
                                   "50000p AP Filtered States"=1))+
    labs(x="Time",
         y="")+
    ggtitle(title)+
    theme_bw()+
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot13<-Filtplot2(dfsv1,"Benchmark comparison")
plot13


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




