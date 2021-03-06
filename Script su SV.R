#############################
#
#   y(t) = exp(x(t)/2)*v(t)
#   x(t) = alpha + beta*x(t-1) + tau*w(t)
#
#   v(t) ~ N(0,1)
#   w(t) ~ N(0,1)
#   x(0) ~ N(0,100)
#
#   PRIOR ON alpha,beta,tau2

#Particle Filter
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
#Particle Filter Optimal Kernel
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
n     =  50
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

#Simulate sv process
#-------------------
set.seed(12345)
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

timeframe<-c(1:length(y))
dfsv<-data.frame(timeframe,y,x)

#Filtering
#---------
set.seed(12345)
N=5000
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
  sd<-Filt$sd
  dataframe<-data.frame(dataframe,mean,sd)
  
  ggplot(dataframe,aes(x=timeframe))+
    geom_line(aes(y=y, col="Observations", linetype="Observations"))+
    geom_line(aes(y=x, col="True States", linetype="True States"))+
    geom_line(aes(y=mean, col="Filtered States", linetype="Filtered States"))+
    geom_ribbon(aes(ymin = mean-1.96*sd, ymax = mean+1.96*sd),
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
    ggtitle(title)+
    theme_bw()+
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) +
    theme(plot.title = element_text(hjust = 0.5))
}

#Plots
#-----
plot1<-Filtplot(dfsv,svpf,"Particle Filter")
plot2<-Filtplot(dfsv,svapf,"Auxiliary Particle Filter")
#set.seed(123)
#svlw<-SVLWfun(y,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda)
plot3<-Filtplot(dfsv,svlw,"Liu e West Filter")
plot4<-Filtplot(dfsv,svopt,"Opt Ker")
plot6<-Filtplot(dfsv,svapfopt,"Opt Ker")

ggarrange(plot1,plot2,plot3)

ggarrange(plot1,plot4)







