---
header-includes: 
- \newcommand{\R}{\textnormal{\sffamily\bfseries R}}
- \def\code#1{\texttt{#1}}
- \newtheorem{algorithm}{Algorithm}[section]

title: "R Notebook"

output:
  pdf_document: 
   latex_engine: xelatex
   fig_width: 5
   fig_height: 3.5
   fig_caption: true
   extra_dependencies: ["float"]
   number_sections: true 
---






\textbf{Implementation}

In this section we present a practical illustration of the algorithms discussed. For the sake of simplicity we use a random walk plus noise model, i.e. the most basic form of a linear Gaussian state-space model. 
\begin{align}
y_{t}|x_{t} & \sim N(x_{t},\sigma^{2}) \\
x_{t}|x_{t-1} & \sim N(x_{t-1},\tau^{2}) \\
x_{0} & \sim N(m_{0},C_{0})
\end{align}
As already mentioned before, in this case the filtering distribution can be computed in closed form solutions using the Kalman filter. However, this toy example will be used also to illustrate more involved filtering strategies described in this work. We believe indeed that it represents a useful starting point to understand the logic of the algorithms which may be eventually replicated when dealing with more complex models. \
The filtering strategy is applied to 50 simulated data. Figure XX shows the simulated true states sequence assuming as data generating process the Equation (2) with $\tau^{2}=1$ and simulated observed sequence process form Equation (1) with $\sigma^{2}=1$.
\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-3-1} 

}

\caption{Simulated random walk plus noise model}\label{fig:unnamed-chunk-3}
\end{figure}
The Kalman Filter for this model is implemented as described below.

\begin{algorithm} Kalman Filter for Random Walk plus Noise Model
\begin{itemize}
\item Initialize $\theta_{0} \sim N(m_{0},C_{0})$
\item For $t=1,...,n$:
\begin{enumerate}
\item Compute the one-step-ahead state predictive distribution at time $t-1$, notice that if $t=1$ no data have been observed yet and therefore $x_{1}|x_{0}\sim N(a_{1},R_{1})$ otherwise
\begin{align*}
x_{t}|y_{1:t-1} & \sim N(a_{t},R_{t})\\
a_{t} & = m_{t-1}\\
R_{t} & = C_{t-1}+\tau^2
\end{align*}
\item Compute the filtering distribution at time $t$ as $p(x_{t}|y_{1:t}) \propto p(x_{t}|y_{1:t-1})p(y_{t}|x_{t})$, i.e. the product of the one-step-ahead state predictive distribution and the likelihood
\begin{align*}
x_{t}|y_{1:t} & \sim N(m_{t},C_{t}) \\
m_{t} & = \big(1-\frac{R_{t}}{R_{t}+\sigma^2}\big)a_{t}+\frac{R_{t}}{R_{t}+\sigma^2}y_{t} \\
C_{t} & = \frac{R_{t}}{R_{t}+\sigma^2}\sigma^2
\end{align*}
\end{enumerate}
\end{itemize}
\end{algorithm}
Our \code{DLM} function implement in \R replicate this steps.
\
\hrule
\hrule
\
\code{DLM(data,sig2,tau2,m0,C0)}
\
\hrule
\textbf{Arguments}

\code{data} \ \ the observed process. It has to be a vector or a univariate time series.\
\code{sig2} \ \ the variance $\sigma^{2}$ in Equation (1)
\
\code{tau2} \ \ the variance $\tau^{2}$ in Equation (2)
\
\code{m0} \ \ central value of the normal prior state distribution
\
\code{C0} \ \ variance of the normal prior state distribution
\hrule
\hrule

```r
DLM<-function(data,sig2,tau2,m0,C0){
  n  = length(data)
  m  = rep(0,n)
  C  = rep(0,n)
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
  }
  return(list(m=m,C=C))
}
```
In Figure XX below filtered states estimated using Kalman Filter with $x_{0} \sim N(0,100)$ and $\sigma^{2}=\tau^{2}=1$ are compared to the true states values.
Notice how closely the filtered states follow the observations and the goodness of the approximation of the true states. 95 percent credible intervals are computed as
\[[E(\theta_{t}|y_{1:t})-z_{1-\alpha/2}\sqrt{V(\theta_{t}|y_{1:t})},E(\theta_{t}|y_{1:t})+z_{1-\alpha/2}\sqrt{V(\theta_{t}|y_{1:t})}]\]
\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-5-1} 

}

\caption{Kalman Filtered States with credible interval (in red)}\label{fig:unnamed-chunk-5}
\end{figure}
The discussion of Kalman Filter will continue in Section XX where we compare it to the Particle Filter.

\
\section{Implementation}
Consider again the random walk plus noise model of section XX,the challenge here is to estimate the filtered states using a Sequential Importance Sampling algorithm. The main idea of the SIS applied to this univariate linear gaussian model is described in the following algorithm. We indicate with $n$ the sample size of time observations and with $N$ the generated sample size for each step of the Sequenial Monte Carlo. 

\begin{algorithm} SIS filter for Random Walk plus Noise Model
\begin{itemize}
\item Let $\{(x_{0},w_{0})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{0}|y_{0})$ such that, for example, $E(g(x_{0})|y_{0}) \approx \sum_{i=1}^{N}w_{0}^{(i)}g(x_{0}^{(i)})$. In particular, initialize $(x_{0}^{(1)},...,x_{0}^{(N)})$ form $N(m_{0},C_{0})$ and set $w_{0}^{(i)}=N^{-1} \ \forall \ i=1,...,N$.
\item For $t=1,...,n$:
\begin{enumerate}
\item Draw $x_{t}^{(i)} \sim N(x_{t-1}^{(i)},\tau^2) \ \ i=1,...,N$ such that $\{(x_{t},w_{t-1})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{t}|y_{t-1})$
\item Set $w_{t}^{(i)} = w_{t-1}^{(i)}f_{N}(y_{t};x_{t}^{(i)},\sigma^2) \ \ i=1,...,N$ such that $\{(x_{t},w_{t})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{t}|y_{t})$
\item Set $p(x_{t}|y_{t})=\sum_{i=1}^{N}w_{t}^{(i)}\delta_{x_{t}^{(i)}}$
\end{enumerate}
\end{itemize}
\end{algorithm}
Our \code{SISfun} function implemented in \R replicate this steps.
\
\hrule
\hrule
\
\code{SISfun(data,N,m0,C0,tau,sigma)}
\
\hrule
\textbf{Arguments}

\code{data} \ \ the observed process. It has to be a vector or a univariate time series.\
\code{N} \ \ number of particles generated at each step
\
\code{m0} \ \ central value of the normal prior state distribution
\
\code{C0} \ \ variance of the normal prior state distribution
\
\code{tau} \ \ the standard deviation $\tau$ in Equation (2)
\
\code{sigma} \ \ the standard deviation $\sigma$ in Equation (1)
\hrule
\hrule

```r
SISfun<-function(data,N,m0,C0,tau,sigma){
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
  for(t in 1:length(data)){
    x    = rnorm(N,x,tau)                   #sample from N(x_{t-1},tau)
    w    = w*dnorm(data[t],x,sigma)         #update weight
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    
    wnorm= w/sum(w)                         #normalized weight
    ESS  = 1/sum(wnorm^2)                   #effective sample size
    
    ess =rbind(ess,ESS)
  }
  
  return(list(xs=xs,ws=ws,ess=ess))
}
```
We have already discussed the reasons why the SIS algorithm does not provide a good strategy in the filtering problem. We provide a graphical intuition of what happens when we use such filtering strategy on a simulated dataset. We decide to set $N=1000$,$m_{0}=0$,$C_{0}=100$ and $\tau=\sigma=1$. The results shown in the following two plots shows a clear degeneration of the effective sample size and bad fit of filtered states with respect to the true values.


\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-8-1} 

}

\caption{a) SIS Filtered States with credible interval (in red). b) Effective sample size.}\label{fig:unnamed-chunk-8}
\end{figure}

\section{Implementation}
In this section, Bootstrap Particle Filter will be used to estimate filtered states of the Random Walk plus Noise introduced in section XX. The overall strategy replicate the SIS filter with the addition of a ESS-based resampling step that applies when the effective sample size is smaller than a predetermined threshold opportunely chosen (in our example we decide to follow a common rule of thumb consisting in setting the threshold at $N/2$). The steps are presented in the following algorithm.

\begin{algorithm} BPF for Random Walk plus Noise Model
\begin{itemize}
\item Let $\{(x_{0},w_{0})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{0}|y_{0})$ such that, for example, $E(g(x_{0})|y_{0}) \approx \sum_{i=1}^{N}w_{0}^{(i)}g(x_{0}^{(i)})$. In particular, initialize $(x_{0}^{(1)},...,x_{0}^{(N)})$ form $N(m_{0},C_{0})$ and set $w_{0}^{(i)}=N^{-1} \ \forall \ i=1,...,N$.
\item For $t=1,...,n$:
\begin{enumerate}
\item Draw $x_{t}^{(i)} \sim N(x_{t-1}^{(i)},\tau^2) \ \ i=1,...,N$ such that $\{(x_{t},w_{t-1})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{t}|y_{t-1})$
\item Set $w_{t}^{(i)} = w_{t-1}^{(i)}f_{N}(y_{t};x_{t}^{(i)},\sigma^2) \ \ i=1,...,N$ such that $\{(x_{t},w_{t})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{t}|y_{t})$
\item if $ESS<N/2$ then
\begin{enumerate}
\item Draw a sample of size N, $(x_{t}^{(1)},...,x_{t}^{(N)})$, from the discrete distribution $P(x_{t}=x_{t}^{(i)})=w_{t}^{(i)},\ \ i=1,...,N$
\item Reset the weights: $w_{t}^{(i)}=N^{-1}$, $i=1,...,N$.
\end{enumerate}
\item Set $p(x_{t}|y_{t})=\sum_{i=1}^{N}w_{t}^{(i)}\delta_{x_{t}^{(i)}}$
\end{enumerate}
\end{itemize}
\end{algorithm}
These steps are resumed in our \code{PFfun} function.
\
\hrule
\hrule
\
\code{PFfun(data,N,m0,C0,tau,sigma,r)}
\
\hrule
\textbf{Arguments}

\code{data} \ \ the observed process. It has to be a vector or a univariate time series.\
\code{N} \ \ number of particles generated at each step
\
\code{m0} \ \ central value of the normal prior state distribution
\
\code{C0} \ \ variance of the normal prior state distribution
\
\code{tau} \ \ the standard deviation $\tau$ in Equation (2)
\
\code{sigma} \ \ the standard deviation $\sigma$ in Equation (1)
\
\code{r} \ \  if present the threshold is set equal to $N/r$ otherwise, if missing, the threshold is set equal to $N/2$
\hrule
\hrule

```r
PFfun<-function(data,N,m0,C0,tau,sigma,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  w  = rep(1/N,N)
   
  for(t in 1:length(data)){
    
    x<-rnorm(N,x,tau)
    w1<-w*dnorm(data[t],x,sigma)
    
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
```
The estimated states of the Boostrap Particle Filter togheter with the effective sample size are shown in Figure XX. We decide to set $N=1000$,$m_{0}=0$,$C_{0}=100$ and $\tau=\sigma=1$. Notice how the resampling step allows the effective sample size not to drop, improving results.


\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-11-1} 

}

\caption{a) PF Filtered States with credible interval (in red). b) Effective sample size (in black) with threshold (in yellow).}\label{fig:unnamed-chunk-11}
\end{figure}

Since the Random Walk plus Noise model allows for closed form solutions, we want to compare the results of Boostrap Particle Filter(BPF) with Kalman Filter(KF). As we can see from figure XX, when the number of particles generated at any interaction increases the results for both the estimated mean and the variance tend to converge [[SPIEGARE PERCHE']]. Moreover, comparing the Root Mean Square Errors of the filtered states with respect to the true state values, when the sample size increases the BPF decreases and, for large N, it reaches the accuracy of the Kalman Filter.




\begin{longtable}[t]{cccc}
\caption{\label{tab:unnamed-chunk-13}Root Mean Square Errors}\\
\toprule
N & Threshold & KF & BPF\\
\midrule
100 & 0.5 & 0.879 & 0.888\\
1000 & 0.5 & 0.879 & 0.886\\
10000 & 0.5 & 0.879 & 0.878\\
\bottomrule
\end{longtable}



\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-14-1} 

}

\caption{Comparison Bootstrap Particle Filter(BPF) and Kalman Filter (KF) for increasing number of generated particles (N)}\label{fig:unnamed-chunk-14}
\end{figure}
\
\section{Implementation}
\
The Bootstrap Particle Filter Approach for the Random Walk plus Noise Model described in Section XX can be improved accounting for the observations in the importance transition density and it consists of generating $x_{t}$ from its conditional distribution given $x_{t-1}$ and $y_{t}$. In the Normal model we are considering, the optimal proposal will be a Normal density as well with mean and variance given by 
\begin{align*}
\mu_{opt}=E(x_{t}|x_{t-1},y_{t})&=x_{t-1}+\frac{\tau^{2}}{\tau^{2}+\sigma^{2}}(y_{t}-x_{t-1})\\
\sigma_{opt}^{2}=V(x_{t}|x_{t-1},y_{t})&=\frac{\tau^{2}\sigma^{2}}{\tau^{2}+\sigma^{2}}
\end{align*}
On the other hand, the incremental weights, using this importance transition density, are proportional to the conditional density of $y_{t}$ given $x_{t-1}=x_{t-1}^{(i)}$ , i.e $N(x_{t-1}^{(i)},\tau^{2}+\sigma^{2})$, evaluated at $y_{t}$. In other words the algorithm implemented in \R is:
\begin{algorithm} GPF for Random Walk plus Noise Model
\begin{itemize}
\item Let $\{(x_{0},w_{0})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{0}|y_{0})$ such that, for example, $E(g(x_{0})|y_{0}) \approx \sum_{i=1}^{N}w_{0}^{(i)}g(x_{0}^{(i)})$. In particular, initialize $(x_{0}^{(1)},...,x_{0}^{(N)})$ from $N(m_{0},C_{0})$ and set $w_{0}^{(i)}=N^{-1} \ \forall \ i=1,...,N$.
\item Compute $\sigma_{opt}^{2}$ 
\item For $t=1,...,n$:
\begin{enumerate}
\item Compute $\mu_{opt}$
\item Draw $x_{t}^{(i)} \sim N(\mu_{opt},\sigma_{opt}^{2}) \ \ i=1,...,N$ such that $\{(x_{t},w_{t-1})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{t}|y_{t-1})$
\item Set $w_{t}^{(i)} = w_{t-1}^{(i)}f_{N}(y_{t};x_{t}^{(i)},\sigma^2+\tau^{2}) \ \ i=1,...,N$ such that $\{(x_{t},w_{t})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{t}|y_{t})$
\item if $ESS<N/2$ then
\begin{enumerate}
\item Draw a sample of size N, $(x_{t}^{(1)},...,x_{t}^{(N)})$, from the discrete distribution $P(x_{t}=x_{t}^{(i)})=w_{t}^{(i)},\ \ i=1,...,N$
\item Reset the weights: $w_{t}^{(i)}=N^{-1}$, $i=1,...,N$.
\end{enumerate}
\item Set $p(x_{t}|y_{t})=\sum_{i=1}^{N}w_{t}^{(i)}\delta_{x_{t}^{(i)}}$
\end{enumerate}
\end{itemize}
\end{algorithm}
The \code{GPFfun} function resume this passages.
\
\hrule
\hrule
\
\code{GPFfun(data,N,m0,C0,tau,sigma,r)}
\
\hrule
\textbf{Arguments}

\code{data} \ \ the observed process. It has to be a vector or a univariate time series.\
\code{N} \ \ number of particles generated at each step
\
\code{m0} \ \ central value of the normal prior state distribution
\
\code{C0} \ \ variance of the normal prior state distribution
\
\code{tau} \ \ the standard deviation $\tau$ in Equation (2)
\
\code{sigma} \ \ the standard deviation $\sigma$ in Equation (1)
\
\code{r} \ \  if present the threshold is set equal to $N/r$ otherwise, if missing, the threshold is set equal to $N/2$
\hrule
\hrule

```r
GPFfun<-function(data,N,m0,C0,tau,sigma,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  importancesd<-sqrt(tau - tau^2 /(tau + sigma))
  predsd <- sqrt(sigma+tau)
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    
    means<-x+(tau/(tau+sigma))*(data[t]-x)
    x<-rnorm(N,means,importancesd)
    w1<-w*dnorm(data[t],x,predsd)
    
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
```



\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-17-1} 

}

\caption{a) GPF Filtered States with credible interval (in red). b) Effective sample size (in black) with threshold (in yellow).}\label{fig:unnamed-chunk-17}
\end{figure}



Let's provide directly a brief comparison between the Bootstrap Particle Filter (BPF) and the Guided Particle Filter (GPF). As we can see from the Figure XX, the Guided Particle Filter provides point estimates that are slightly better with respect to the ones of the BPF, and this is confirmed by the Table showing the RMSE. Moreover, also the variance is better, suggesting for higher precision.




\begin{longtable}[t]{cccc}
\caption{\label{tab:unnamed-chunk-19}Root Mean Square Errors}\\
\toprule
N & Threshold & BPF & GPF\\
\midrule
1000 & 0.50 & 0.916 & 0.880\\
1000 & 0.25 & 0.895 & 0.862\\
1000 & 0.10 & 0.869 & 0.881\\
\bottomrule
\end{longtable}


\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-20-1} 

}

\caption{Comparison Bootstrap Particle Filter(BPF) and Guided Particle Filter (GPF), number of generated particles N=1000}\label{fig:unnamed-chunk-20}
\end{figure}




\textbf{Implementation}

For illustration purposes, we are going implement an auxiliary particle filter with auxiliary function $g(x_{t-1})=E(x_{t}|x_{t-1})=x_{t-1}$.  

\begin{algorithm} APF for Random Walk plus Noise Model
\begin{itemize}
\item Let $\{(x_{0},w_{0})^{(i)}\}_{i=1}^{N}$ summarizes $p(x_{0}|y_{0})$ such that, for example, $E(g(x_{0})|y_{0}) \approx \sum_{i=1}^{N}w_{0}^{(i)}g(x_{0}^{(i)})$. In particular, initialize $(x_{0}^{(1)},...,x_{0}^{(N)})$ from $N(m_{0},C_{0})$ and set $w_{0}^{(i)}=N^{-1} \ \forall \ i=1,...,N$.
\item For $t=1,...,n$:
\begin{enumerate}
\item For $k=1,...,N$:
\begin{enumerate}
\item Draw $I_{k}$ with $P(I_{k}) \propto w_{t-1}^{(i)}f(y_{t}|g(x_{t-1}^{(i)})) $
\item Draw $x_{t}^{(k)} \sim N(x_{t-1}^{(I_{k})},\tau^2)$
\item Set  $\tilde{w}_{t}^{(k)} = \frac{f_{N}(y_{t}|x_{t}^{(k)})}{f_{N}(y_{t}|g(x_{t-1}^{(I_{k})}))}$
\end{enumerate}
\item Normalize the weights: $w_{t}^{(i)}=\frac{\tilde{w}_{t}^{(i)}}{\sum_{j=1}^{N}(\tilde{w}_{t}^{(j)})}$
\item Compute $ESS=\Bigg(\sum_{i=1}^{N}(w_{t}^{(i)})^{2}\Bigg)^{-1}$
\item if $ESS<N/2$ then
\begin{enumerate}
\item Draw a sample of size N, $(x_{t}^{(1)},...,x_{t}^{(N)})$, from the discrete distribution $P(x_{t}=x_{t}^{(i)})=w_{t}^{(i)},\ \ i=1,...,N$
\item Reset the weights: $w_{t}^{(i)}=N^{-1}$, $i=1,...,N$.
\end{enumerate}
\item Set $p(x_{t}|y_{1:t})=\sum_{i=1}^{N}w_{t}^{(i)}\delta_{x_{t}^{(i)}}$
\end{enumerate}
\end{itemize}
\end{algorithm}
The \code{APFfun} function resume this passages.
\
\hrule
\hrule
\
\code{APFfun(data,N,m0,C0,tau,sigma,r)}
\
\hrule
\textbf{Arguments}

\code{data} \ \ the observed process. It has to be a vector or a univariate time series.\
\code{N} \ \ number of particles generated at each step
\
\code{m0} \ \ central value of the normal prior state distribution
\
\code{C0} \ \ variance of the normal prior state distribution
\
\code{tau} \ \ the standard deviation $\tau$ in Equation (2)
\
\code{sigma} \ \ the standard deviation $\sigma$ in Equation (1)
\
\code{r} \ \  if present the threshold is set equal to $N/r$ otherwise, if missing, the threshold is set equal to $N/2$
\hrule
\hrule


```r
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
```



\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-23-1} 

}

\caption{a) APF Filtered States with credible interval (in red). b) Effective sample size (in black) with threshold (in yellow).}\label{fig:unnamed-chunk-23}
\end{figure}

Let's compare the Auxiliary Particle Filter (APF) and the Bootstrap Particle Filter (BPF).

\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-24-1} 

}

\caption{Comparison Bootstrap Particle Filter(BPF) and Guided Particle Filter (GPF), number of generated particles N=1000}\label{fig:unnamed-chunk-24}
\end{figure}




\begin{longtable}[t]{cccc}
\caption{\label{tab:unnamed-chunk-26}Root Mean Square Errors}\\
\toprule
N & Threshold & BPF & APF\\
\midrule
1000 & 0.50 & 0.885 & 0.875\\
1000 & 0.25 & 0.893 & 0.892\\
1000 & 0.10 & 0.855 & 0.870\\
\bottomrule
\end{longtable}



```r
APFoptfun<-function(data,N,m0,C0,tau,sigma,r){
  if(missing(r)){r=2}else{}
  xs<-NULL
  ws<-NULL
  ess<-NULL
  x  = rnorm(N,m0,sqrt(C0))
  importancesd<-sqrt(tau - tau^2 /(tau + sigma))
  predsd <- sqrt(sigma+tau)
  w  = rep(1/N,N)
  
  for(t in 1:length(data)){
    ESS  = 1/sum(w^2)
    
    if(ESS<(N/r)){
    weight = w*dnorm(data[t],x,predsd)
    k   = sample(1:N,size=N,replace=TRUE,prob=weight)
    }else{
    weight = rep(1/N,N)
    k   = sample(1:N,size=N,replace=TRUE,prob=weight)
    }
    
    means<-x[k]+(tau/(tau+sigma))*(data[t]-x[k])
    x1   = rnorm(N,means,importancesd)
    lw  = dnorm(data[t],x1,predsd,log=TRUE)-dnorm(data[t],x[k],predsd,log=TRUE)
    w   = exp(lw)
    w   = w/sum(w)
    x <- x1
   
    
    xs = rbind(xs,x)
    ws = rbind(ws,w)
    ess =rbind(ess,ESS)
    
  }
  return(list(xs=xs,ws=ws,ess=ess))
}
```



\begin{figure}[ht]

{\centering \includegraphics{draft-implementation-2_files/figure-latex/unnamed-chunk-29-1} 

}

\caption{a) APF Filtered States with credible interval (in red). b) Effective sample size (in black) with threshold (in yellow).}\label{fig:unnamed-chunk-29}
\end{figure}



\
\section{Implementation}
\
Consider the linear Gaussian example of section XX, but this time with unknown variances $\tau^{2}$ and $\sigma^{2}$. Thus, let $\psi=(\sigma^{2},\tau^{2})$ be the unknown parameter vector and assign a gamma prior for its components,
\begin{align*}
\sigma^{2}  & \sim G(\alpha_{v},\beta_{v}) \\
\tau^{2}  & \sim G(\alpha_{w},\beta_{w})
\end{align*}
Alternatively, assign them a uniform prior if we have no knowledge of the hyperparameters. The algorithm follows the following steps.
\begin{algorithm} LWF for Random Walk plus Noise Model
\begin{itemize}
\item Initialize $(x_{0}^{(1)},...,x_{0}^{(N)})$ from $N(m_{0},C_{0})$, $({\sigma^{2}}^{(1)},...,{\sigma^{2}}^{(N)})$ from $G(\alpha_{v},\beta_{v})$ and $({\tau^{2}}^{(1)},...,{\tau^{2}}^{(N)})$ from $G(\alpha_{w},\beta_{w})$. Set $w_{0}^{(i)}=N^{-1} \ \forall \ i=1,...,N$. Therefore $\psi^{(i)}=({\sigma^{2}}^{(i)},{\tau^{2}}^{(i)})$, and
$\hat{\pi}_{0}=p(x_{0}|y_{0})=\sum_{i=1}^{N}w_{0}^{(i)}\delta_{(x_{0}^{(i)},\psi^{(i)})}$
\item For $t=1,...,n$:.
\begin{enumerate}
\item Compute $\overline{\psi}=E_{\hat{\pi}_{t-1}}(\psi)$ and $\Sigma=V_{\hat{\pi}_{t-1}}(\psi)$. For $i=1,...,N$, set
\begin{align*}
m^{(i)} & = a\psi^{(i)}+(1-a)\hat{\psi} \\
\hat{x}_{t}^{(i)} & = E(x_{t}|x_{t-1}=x_{t-1}^{(i)},\psi=\psi^{(k)})
\end{align*}
\item For $k=1,...,N$:
\begin{itemize}
\item Draw $I_{k}$, with $P(I_{k}=i) \propto w_{t-1}^{(i)}f_{N}(y_{t}|g(x_{t-1}^{(i)}),\psi=m^{(i)}) $
\item Draw ${\sigma^{2}}^{(k)}$ from $G(\alpha_{v}^{(I_k)},\beta_{v}^{(I_k)})$ 
\item Draw ${\tau^{2}}^{(k)}$ from $G(\alpha_{w}^{(I_k)},\beta_{w}^{(I_k)})$
\item Draw $x_{t}^{(k)}$ from $N(x_{t-1}^{(I_k)},\psi=\psi^{(k)})$
\item Set  $\tilde{w}_{t}^{(k)} = \frac{f_{N}(y_{t}|x_{t}^{(k)},\psi=\psi^{(k)})}{f_{N}(y_{t}|g(x_{t-1}^{(I_{k})},\psi=m^{(I_k)}))}$
\end{itemize}
\item Compute $ESS=\Bigg(\sum_{i=1}^{N}(w_{t}^{(i)})^{2}\Bigg)^{-1}$
\item if $ESS<N/2$ then
\begin{enumerate}
\item Draw a sample of size N, $(x_{t}^{(1)},...,x_{t}^{(N)})$, from the discrete distribution $P(x_{t}=x_{t}^{(i)})=w_{t}^{(i)},\ \ i=1,...,N$
\item Reset the weights: $w_{t}^{(i)}=N^{-1}$, $i=1,...,N$.
\end{enumerate}
\item Set $\hat{\pi}_{t}=p(x_{t}|y_{1:t})=\sum_{i=1}^{N}w_{t}^{(i)}\delta_{(x_{t}^{(i)},\psi^{(i)})}$
\end{enumerate}
\end{itemize}
\end{algorithm}
Our \code{LWfun} function goes through the illustrated steps.
\
\hrule
\hrule
\
\code{LWfun(data,N,m0,C0,alphav,betav,alphaw,betaw,delta,unif,r)}
\
\hrule
\textbf{Arguments}

\code{data} \ \ the observed process. It has to be a vector or a univariate time series.\
\code{N} \ \ number of particles generated at each step
\
\code{m0} \ \ central value of the normal prior state distribution
\
\code{C0} \ \ variance of the normal prior state distribution
\
\code{alphav, betav} \ \  Gamma prior hyperparameters on $\sigma^{2}$
\
\code{alphaw,betaw} \ \  Gamma prior hyperparameters on $\tau^{2}$
\
\code{delta} \ \ hyperparameter delta value
\
\code{unif} \ \ if True then it sets a Uniform $(0,10)$ prior on $\sigma^{2}$ and $\tau^{2}$
\
\code{r} \ \  if present the threshold is set equal to $N/r$ otherwise, if missing, the threshold is set equal to $N/2$
\hrule
\hrule

```r
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
    
    muW = a*pars[,2]+(1-a)*meanW
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
```

