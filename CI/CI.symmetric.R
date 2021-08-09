F_order<-function(xs,ws){
  t=dim(xs)[1]
  n=dim(xs)[2]
  F_order=matrix(NA,t,n)
  for (i_t in 1:t){
    order=order(xs[i_t,])
    ws.o=ws[i_t,][order]
    F_order[i_t,1]=ws.o[1]
    for (i_n in 1:(n-1)){
      F_order[i_t,i_n+1]=F_order[i_t,i_n]+ws.o[i_n+1]
    }
  }
  return(F_order)
}

CI.symmetric<-function(xs,ws,p=0.05){
  t=dim(xs)[1]
  n=dim(xs)[2]
  CI=matrix(NA,t,2)
  F_order=F_order(xs,ws)
  lower.mat=(F_order>p/2)
  higher.mat=(F_order>1-p/2)
  for (i_t in 1:t){
    xs.o=xs[i_t,][order(xs[i_t,])]
    if (lower.mat[i_t,1]==T){
      lower=xs.o[1]
    }else{
      for (i_n in 1:(n-1)){
        if (lower.mat[i_t,i_n+1]==T & lower.mat[i_t,i_n]==F){
          lower=(xs.o[i_n]+xs.o[i_n+1])/2
        }
      }
    }
    
    for (i_n in (n-1):1){
      if (higher.mat[i_t,i_n+1]==T & higher.mat[i_t,i_n]==F){
        higher=(xs.o[i_n+1]+xs.o[i_n])/2
      }
    }

    CI[i_t,]=c(lower,higher)
  }
  return(CI)
}