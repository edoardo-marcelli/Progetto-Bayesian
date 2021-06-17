library(usethis)

use_git_config(user.name = "edoardo-marcelli", 
               user.email = "edoardomarcelli2@gmail.com")
use_git_config(user.name = "michelebolla", 
               user.email = "michele.bolla@studbocconi.it")



#@@@@@@@@@@@@@@@@@@@@##
#SEQUENTIAL MONTECARLO#
#@@@@@@@@@@@@@@@@@@@@##
library(ggplot2)

#Generate a random walk plus noise
#---------------------------------
#n = sample dim
#v = obs noise
#w = state noise
#theta0 = init state

dlm.sim = function(n,v,w,theta0){
  theta    = rep(0,n)
  y        = rep(0,n)
  theta[1] = rnorm(1,theta0,w)
  y[1] = rnorm(1,theta[1],v)
  for (t in 2:n){
    theta[t] = rnorm(1,theta[t-1],w)
    y[t] = rnorm(1,theta[t],v)
  }
  return(list(theta=theta,y=y))
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
theta0=0

sim<-dlm.sim(n,v,w,theta0)
y<-sim$y
theta<-sim$theta

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
DLM.df<-data.frame(timeframe,y,theta,m,C,Low,Up)

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

#vediamo se funziona!!

#ciao sono Michele
# test (Alexei)
# alexei teast from local
