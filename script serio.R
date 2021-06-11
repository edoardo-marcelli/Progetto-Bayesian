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

dlm.sim = function(n,v,w,x0){
  x    = rep(0,n)
  y        = rep(0,n)
  x[1] = rnorm(1,x0,w)
  y[1] = rnorm(1,x[1],v)
  for (t in 2:n){
    x[t] = rnorm(1,x[t-1],w)
    y[t] = rnorm(1,x[t],v)
  }
  return(list(x=x,y=y))
}

#Filtering, using data
#---------------------
#y = observed data

DLM = function(y,v2,w2,m0,C0){
  n = length(y)
  m = rep(0,n)
  C = rep(0,n)
  for (t in 1:n){
    if (t==1){
      a = m0
      R = C0 + v2
    }else{
      a = m[t-1]
      R = C[t-1] + v2
    }
    A = R/(R+w2)
    m[t] = (1-A)*a + A*y[t]
    C[t] = A*w2
  }
  return(list(m=m,C=C))
}


########################################################

#Generate a random walk plus noise
#---------------------------------
n=500
v2=1
w2=1
v=sqrt(v2)
w=sqrt(w2)
x0=0

sim<-dlm.sim(n,v,w,x0)
y<-sim$y
x<-sim$x

#Kalman Filter
#-------------
m0=0
C0=100
filtval<-DLM(y,v2,w2,m0,C0)
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
#matrix of (x,w) 
wss = array(0,c(1,n,N))
xss = array(0,c(1,n,N))
Ess = matrix(0,1,n)

#
  x   = rnorm(N,m0,sqrt(C0))
  w   =  rep(1/N,N)
  ws  = NULL
  xs  = NULL
  ess = NULL
  for (t in 1:n){
    x    = rnorm(N,x,w)
    w    = w*dnorm(y[1],x,v)
    w1   = w/sum(w)
    ESS  = 1/sum(w1^2)
    draw = sample(x,size=N,replace=T,prob=w1)
    xs   = rbind(xs,draw)
    ws   = rbind(ws,w1)
    ess  = c(ess,ESS)
  }
  wss[i,,] = ws
  xss[i,,] = xs
  Ess[i,] = ess
}
mx = matrix(0,n,nsig)
Lx = matrix(0,n,nsig)
Ux = matrix(0,n,nsig)
for (i in 1:nsig){
  mx[,i] = apply(xss[i,,],1,median)
  Lx[,i] = apply(xss[i,,],1,q025)
  Ux[,i] = apply(xss[i,,],1,q975)
}


