## model
#   SV
#   y(t) = a+exp(x(t)/2)*v(t)
#   x(t) = a.v + b.v*x(t-1) + s.v*w(t)
#
#   v(t) ~ N(0,1)
#   w(t) ~ N(0,1)
#   x(0) ~ N(0,100)
#
#   PRIOR ON alpha,beta,tau2
## packages
library(tidyverse)
## functions
SVPFfun<-function(data,N,m0,C0,a,a.v,b.v,s.v,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  xp.1<-rnorm(N,a.v+b.v*x,s.v)
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    x<-rnorm(N,a.v+b.v*x,s.v)
    w1<-w*dnorm(data[t],a,exp(x/2))
    
    w = w1/sum(w1)
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
      index<-sample(N,size=N,replace=T,prob=w)
      x<-x[index]
      w<-rep(1/N,N)
    }else{}
    
    xs  = rbind(xs,x)
    ws  = rbind(ws,w)
    ess = rbind(ess,ESS)
  }
  return(list(xs=xs,ws=ws,ess=ess,xp.1=xp.1))
}
# SVPFfun.f<-function(data,N,m0,C0,a,a.v,b.v,s.v,r){
#   if(missing(r)){r=2}else{}
#   xs<-NULL
#   xf<-NULL
#   ws<-NULL
#   ess<-NULL
#   x  = rnorm(N,m0,sqrt(C0))
#   w  = rep(1/N,N)
#   
#   for(t in 1:length(data)){
#     
#     x<-rnorm(N,a.v+b.v*x,s.v)
#     x.f1<-x
#     w1<-w*dnorm(data[t],a,exp(x/2))
#     x.f<-rnorm(N,a.v+b.v*x,s.v)
#     w = w1/sum(w1)
#     ESS  = 1/sum(w^2)
#     
#     if(ESS<(N/r)){
#       index<-sample(N,size=N,replace=T,prob=w)
#       x<-x[index]
#       w<-rep(1/N,N)
#     }else{}
#     
#     xs = rbind(xs,x)
#     xs.f = rbind(xs.f,x.f) 
#     ws = rbind(ws,w)
#     ess =rbind(ess,ESS)
#   }
#   return(list(xs=xs,xs.f=xs.f,ws=ws,ess=ess,x.f1=x.f1))
# }
## data
data<-read.csv("sp500_17-21.csv", header=T)
data$t<-as.Date(data$t)
y<-data$r
#y<-rnorm(996,mean=1,sd=1) #test
n<-length(y)

# find parameters for constant and stochastic vol models
mean=sum(y)/n
res<-(y-mean)^2
lres<-log(res)
l.lres<-c(NA,lres[1:n-1])
resdf=data.frame(lres,l.lres)
lm_res<-lm(lres~1+l.lres,data=resdf)
c_res<-lm_res[[1]]
sigma_res<-sqrt(mean(summary(lm_res)$residuals^2))
meansq_res<-sqrt(sum(res)/n)

a<-mean
s<-meansq_res

a.v<-c_res[1]
b.v<-c_res[2]
s.v<-sigma_res

m0=0
C0=100
# calculate cumulative marginal likelihood
# CV model
l_cv<-dnorm(y,a,s)
ll_cv<-log(l_cv)
cll_cv<-ll_cv
for (i in 2:n){
  cll_cv[i]=ll_cv[i]+cll_cv[i-1]
}
# SV model
N=1000
pf<-SVPFfun(data=y,N=N,m0=m0,C0=C0,a=a,a.v=a.v,b.v=b.v,s.v=s.v,r=4)
xs<-pf$xs
ws<-pf$ws
# make particles for predictive state distribution
xp<-xs
xp[1,1:N]<-pf$xp.1
for (i in 2:n){
  xp[i,1:N]=rnorm(N,a.v+b.v*xs[i-1,1:N],s.v)
}
wp<-ws
wp[1,1:N]=1/N
wp[2:n,1:N]=ws[1:n-1,1:N]
l_sv<-rep(0,n)
for (i in 1:n){
  for (j in 1:N){
  l_sv[i]=l_sv[i]+dnorm(y[i],a,exp(xp[i,j]/2))*wp[i,j]
  }
}
ll_sv<-log(l_sv)
cll_sv<-ll_sv
for (i in 2:n){
  cll_sv[i]=ll_sv[i]+cll_sv[i-1]
}
# Difference
cll_diff=cll_cv-cll_sv
ll_diff=ll_cv-ll_sv
cll_diff.df=data.frame(data$t,cll_diff)
colnames(cll_diff.df)<-c('t','cll_diff')

cll_plot<-ggplot(data = cll_diff.df, aes(x = t, y = cll_diff))+
geom_line()

save(cll_plot, file="cll_plot.Rda")
