################################
#
# FUNZIONI CON FORECAST
#
###############################



DLM<-function(data,sig2,tau2,m0,C0){
  n  = length(data)
  m  = rep(0,n)
  C  = rep(0,n)
  f  = rep(0,n)
  Q  = rep(0,n)
  for (t in 1:n){
    if (t==1){
      a = m0
      R = C0 + tau2
      j = a
      K = R + sig2
    }else{
      a = m[t-1]
      R = C[t-1] + tau2
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

PFfun<-function(data,N,m0,C0,tau,sigma,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ys<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    x<-rnorm(N,x,tau)
    y<-rnorm(N,x,sigma)
    w1<-w*dnorm(data[t],x,sigma)
    
    w = w1/sum(w1)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x[index]
      w<-rep(1/N,N)
    }else{}
    
    ys = rbind(ys,y)
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess,ys=ys))
}

APFfun<-function(data,N,m0,C0,tau,sigma,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    weight = w*dnorm(data[t],x,sigma)
    k   = sample(1:N,size=N,replace=TRUE,prob=weight)
    x1   = rnorm(N,x[k],tau)
    lw  = dnorm(data[t],x1,sigma,log=TRUE)-dnorm(data[t],x[k],sigma,log=TRUE)
    w   = exp(lw)
    w   = w/sum(w)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x1<-x1[index]
      w<-rep(1/N,N)
    }else{}
    
    x <- x1
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
    
  }
  return(list(xs=xs,ws=ws,ess=ess))
}

LWfun<-function(data,N,m0,C0,alphav,betav,alphaw,betaw,delta,unif,r){
  if(missing(r)){r=2}else{}
  xs     = rnorm(N,m0,sqrt(C0))
  if(unif==T){
  pars   = cbind(runif(N,0,10),runif(N,0,10))}else{}
  pars   = cbind(rgamma(N,shape=alphav,scale=betav),rgamma(N,shape=alphaw,scale=betaw))
  a      = (3*delta-1)/(2*delta)
  h2     = 1-a^2
  parss  = array(0,c(N,2,n))
  xss    = NULL
  ws     = NULL
  ess    = NULL
  w      = rep(1/N,N)
  for (t in 1:length(data)){
    meanV = weighted.mean(pars[,1],w)
    varV  = weighted.mean((pars[,1]-meanV)^2,w)
    meanW = weighted.mean(pars[,2],w)
    varW  = weighted.mean((pars[,2]-meanW)^2,w)
    
    muV = a*pars[,1]+(1-a)*meanV
    sigma2V = (1-a^2)*varV
    alphaV = muV^2/sigma2V
    betaV = muV/sigma2V
    
    muW = a*pars[,1]+(1-a)*meanW
    sigma2W = (1-a^2)*varW
    alphaW = muW^2/sigma2W
    betaW = muW/sigma2W
    
    weight      = w*dnorm(data[t],xs,sqrt(muV))
    k           = sample(1:N,size=N,replace=T,prob=weight)
    
    pars[,1]<-rgamma(N,shape=alphaV[k],rate=betaV[k])
    pars[,2]<-rgamma(N,shape=alphaW[k],rate=betaW[k])
    
    xsprevious<-xs[k]
    xs = rnorm(N,xs[k],sqrt(pars[,2]))
    
    w           = exp(dnorm( data[t],xs,sqrt(pars[,1]),log=T)-
                        dnorm( data[t],xsprevious,sqrt(muV[k]),log=T))
    w           = w/sum(w)
    ESS         = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      xs<-xs[index]
      pars<-pars[index,]
      w<-rep(1/N,N)
    }else{
      xs<-xs
      pars<-pars
    }
    
    
    xss         = rbind(xss,xs)
    parss[,,t]  = pars 
    ws          = rbind(ws,w)
    ess         = rbind(ess,ESS)
  }
  return(list(xss=xss,parss=parss,ws=ws,ess=ess))
}


# STOCHASTIC VOLATILITY


SVPFfun<-function(data,N,m0,C0,alpha,beta,tau,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ys<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
   
  for(t in 1:length(data)){
    
    x<-rnorm(N,alpha+beta*x,tau)
    y<-rnorm(N,0,exp(x/2))
    w1<-w*dnorm(data[t],0,exp(x/2))
    
    w = w1/sum(w1)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x[index]
      w<-rep(1/N,N)
    }else{}
    
    ys = rbind(ys,y)
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess,ys=ys))
}

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
    lw  = dnorm(data[t],0,exp(x/2),log=TRUE)-dnorm(data[t],0,exp(x/2),log=TRUE)
    w   = exp(lw)
    w   = w/sum(w)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x1<-x1[index]
      w<-rep(1/N,N)
    }else{}
    
    x <- x1
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
    
  }
  return(list(xs=xs,ws=ws,ess=ess))
}


SVLWfun<-function(data,N,m0,C0,ealpha,valpha,ebeta,vbeta,nu,lambda){

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
return(list(xs=xss,pars=parss,ws=ws))
}



















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

n=50
tau2=1
sigma2=1
sigma=sqrt(sigma2)
tau=sqrt(tau2)
x0=0

sim<-dlm.sim(n,sigma,tau,x0)
y<-sim$y
x<-sim$x

timeframe<-c(1:length(y))
df<-data.frame(timeframe,y,x)

#Kalman Filter
#-------------
m0=0
C0=100
filtval<-DLM(y,sigma2,tau2,m0,C0)
f<-filtval$f
Q<-filtval$Q
Low = f + qnorm(0.025)*sqrt(Q)
Up = f + qnorm(0.975)*sqrt(Q)

#Put everything in a data frate (for ggplot)
timeframe<-c(1:length(y))
DLM.df<-data.frame(timeframe,y,x,f,Q,Low,Up)


ggplot(DLM.df,aes(x=timeframe))+
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


Forecastvalues<-function(fun){
  yhat<-sapply(1:n,function(i)
    weighted.mean(fun$ys[i,],fun$ws[i,]))
  sdhat<-sapply(1:n,function(i)
    sqrt(weighted.mean((fun$ys[i,]-yhat[i])^2,fun$ws[i,])))
  
  return(list(mean=yhat,sd=sdhat))
  
}
svfpf<-SVPFfun(y,N,0,100,0,1,1)


Forecastplot<-function(dataframe,fun,title){
  
  For<-Forecastvalues(fun)
  mean<-For$mean
  sd<-For$sd
  dataframe<-data.frame(dataframe,mean,sd)
  
  ggplot(dataframe,aes(x=timeframe))+
    geom_line(aes(y=y, col="Observations", linetype="Observations"))+
    geom_line(aes(y=x, col="True States", linetype="True States"))+
    geom_line(aes(y=mean, col="Filtered States", linetype="Filtered States"))+
    geom_ribbon(aes(ymin = mean-1.96*sd, ymax = mean+1.96*sd),
                fill="blue",alpha=0.16) +
    scale_color_manual("",
                       values=c("Observations" = "black",
                                "True States" = "black",
                                "Filtered States"= "blue"))+
    scale_linetype_manual("",
                          values=c("Observations" = 1,
                                   "True States" = 3,
                                   "Filtered States"=1))+
    coord_cartesian(ylim=c(-4,4))+
    labs(x="Time",
         y="")+
    ggtitle(title)+
    theme_bw()+
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) +
    theme(plot.title = element_text(hjust = 0.5))
}

Forecastplot(dfsv,svfpf,title="Sequential Importance Sampling Filter")
