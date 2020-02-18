
According to previous experiments, simple average (sa) has been the best performing method.In order to verify the correctness of our calculation process, we try to make only one-step prediction, which is consistent with the method in the optimal prediction pool.
* Data: M4 quarterly data (1000 randomly)
* One-step prediction
# Code
```
library(tsfeatures)
library(M4metaresults)
library(M4comp2018)
library(M4metalearning)
library(forecast)
library(tseries)
load("BMA_PP.RData")
load("BMA_FF.RData")
load("M4_q1.RData")
```

```
for (a in 1:1000) {
y<-M4_q1[[a]]$x
p<-PP[[a]]
features_y<-FF[[a]]

#feature based
log_score1<-log(0.5*p[1,1]+0.5*p[1,2])+log(0.5*p[2,1]+0.5*p[2,2])+log(0.5*p[3,1]+0.5*p[3,2])+log(0.5*p[4,1]+0.5*p[4,2])+log(0.5*p[5,1]+0.5*p[5,2])+log(0.5*p[6,1]+0.5*p[6,2])+log(0.5*p[7,1]+0.5*p[7,2])+log(0.5*p[8,1]+0.5*p[8,2])+log(0.5*p[9,1]+0.5*p[9,2])
log_score<-function(beta){
  for (i in 10:(length(y))) {
    w<-1/(1+exp(-features_y[i-9,feature_select[[11]]]%*%beta))
    log_score1<-log_score1+log(w*p[i,1]+(1-w)*p[i,2])
  }
  return(-log_score1)
}

set.seed(2019-02-14)
#w_max<-optim(fn=log_score,par=runif((length(feature_select[[f]])+1), 0, 0),method="SANN")
w_max<-optim(fn=log_score,par=runif(6, 0, 0),method="SANN")
if(w_max$convergence!=0){
  print(a)
  }
beta_optim<-w_max$par

w<-1/(1+exp(-features_y[length(features_y[,1]),feature_select[[11]]]%*%beta_optim))
y_hat_feature<-w*(M4_q1[[a]]$ets_fore[1])+(1-w)*(M4_q1[[a]]$ari_fore[1])
M4_q1[[a]]$y_hat_feature<-y_hat_feature

#optimal pool
y1<-0
y_score<-function(x){
  for (i in 1:(length(y))) {
    y1<-y1+log(x*p[i,1]+(1-x)*p[i,2])
  }
  return(-y1)
}

w_optim<-optimize(y_score,c(0,1),tol = 0.0001)
w_optim<-w_optim$minimum
M4_q1[[a]]$y_hat_w<-w_optim*M4_q1[[a]]$ets_fore[1]+(1-w_optim)*M4_q1[[a]]$ari_fore[1]

#SA
M4_q1[[a]]$average<-0.5*(M4_q1[[a]]$ets_fore[1]+M4_q1[[a]]$ari_fore[1])

#error
ff<-rbind(M4_q1[[a]]$y_hat_feature,M4_q1[[a]]$y_hat_w,M4_q1[[a]]$average,M4_q1[[a]]$ari_fore[1],M4_q1[[a]]$ets_fore[1])
row.names(ff)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
M4_q1[[a]]$ff<-ff

lentry <- M4_q1[[a]]
insample <- lentry$x
ff<- lentry$ff
frq<-stats::frequency(insample)
insample <- as.numeric(insample)
outsample <- as.numeric(lentry$xx)
masep <- mean(abs(utils::head(insample, -frq) - utils::tail(insample, -frq)))
repoutsample <- matrix(rep(outsample, each = nrow(ff)), nrow = nrow(ff))[,1]
smape_err <- 200 * abs(ff - repoutsample)/(abs(ff) + abs(repoutsample))
mase_err <- abs(ff - repoutsample)/masep
M4_q1[[a]]$mase_err <- mase_err
M4_q1[[a]]$smape_err <- smape_err

ari_score<-dnorm(outsample[1],mean =ff[4,1],sd=M4_q1[[a]]$sigma_ari )
ets_score<-dnorm(outsample[1],mean =ff[5,1],sd=M4_q1[[a]]$sigma_ets )
ave_score<-0.5*dnorm(outsample[1],mean =ff[4,1],sd=M4_q1[[a]]$sigma_ari)+0.5*dnorm(outsample[1],mean =ff[5,1],sd=M4_q1[[a]]$sigma_ets)
hat_w_score<-(1-w_optim)*dnorm(outsample[1],mean =ff[4,1],sd=M4_q1[[a]]$sigma_ari)+w_optim*dnorm(outsample[1],mean =ff[5,1],sd=M4_q1[[a]]$sigma_ets)
hat_feature_score<-(1-w)*dnorm(outsample[1],mean =ff[4,1],sd=M4_q1[[a]]$sigma_ari)+w*dnorm(outsample[1],mean =ff[5,1],sd=M4_q1[[a]]$sigma_ets)
score<-c(hat_feature_score,hat_w_score,ave_score,ari_score,ets_score)
names(score)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
M4_q1[[a]]$score<-score
}
```
# Result
