# Test 1 One-step prediction
According to previous experiments, simple average (SA) has been the best performing method. In order to verify the correctness of our calculation process, we try to make only one-step prediction, which is consistent with the method in the optimal prediction pool.
* Data: M4 quarterly data (1000 randomly)
## Code
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
## Result
```
> summery
          y_hat_feature    y_hat_w    average   ari_fore   ets_fore
mase_err      0.5563622  0.5575904  0.5523294  0.5598119  0.5724977
smape_err     5.4888827  5.4729528  5.3708209  5.3339046  5.6922742
log_score    -6.6685471 -6.6634350 -6.6399594 -6.6942272 -6.6775508
```
* The method of optimal prediction pool is still inferior to the SA.
* For one-step prediction, there is no obvious advantage in taking features into consideration.
## Discuss
Plot the histogram of the weights in the first two methods.

![hist_weight](/plot/hist_weight.png)

Select the data with the highest frequency in the histogram of the two methods, that is, the data with weight less than or equal to 0.1, and calculate the performance of the different methods in these data. However, the results imply that SA is still the best.
* the optimal pool
```
optim_select<-c()
for (i in 1:1000){
  if(M4_q1[[i]]$w_optim<=0.1){
    optim_select<-c(optim_select,i)
  }
}

> length(optim_select)
[1] 397
```
```
             y_hat_w    average   ari_fore  ets_fore
mase_err   0.5262258  0.5197022  0.5264386  0.541758
smape_err  5.2769434  5.1338279  5.2966418  5.505232
log_score -6.5272805 -6.4907724 -6.5288729 -6.530911
```
* feature based
```
feature_select<-c()
for (i in 1:1000){
  if(M4_q1[[i]]$w_feature<=0.1){
    feature_select<-c(feature_select,i)
  }
}

> length(feature_select)
[1] 439
```
```
          y_hat_feature    y_hat_w    average   ari_fore   ets_fore
mase_err       0.478759  0.4770845  0.4673669  0.4791699  0.4848627
smape_err      4.833397  4.8080849  4.6521682  4.8417760  4.8846863
log_score     -6.465986 -6.4522726 -6.4190756 -6.4841300 -6.4491572
```
## Questions
* The calculation of probability density 

# Test 2  Weights are not updated
```
for (a in 1:1000) {
  w<-as.numeric(M4_q1[[a]]$w_feature)
  y_hat_feature<-w*(as.matrix(M4_q1[[a]]$ets_fore))+(1-w)*(as.matrix(M4_q1[[a]]$ari_fore))
  M4_q1[[a]]$y_hat_feature<-t(y_hat_feature)
  ...
}
```
```
> print(summery)
[1]   1.103422   9.370611 -85.723277
1.106628   9.417214 -86.822581
```
The results are roughly consistant with the method for updating weights (1.1036	9.351	-85.73).
# Test 3 Yearly data
1000 of M4 yearly data were randomly selected, and the feature-based method only considered the selected six features.
*The predictions for individual data were poor, with log score -Inf, so these data were removed for comparison.
```
> summery
          y_hat_feature     y_hat_w     average    ari_fore    ets_fore
mase_err       3.328958    3.325678    3.270496    3.349635    3.359202
smape_err     14.122818   14.153858   13.963836   14.379974   14.291773
log_score   -119.456855 -118.539386 -114.111131 -136.114439 -126.652622
```
```
score_output<-c()
for (j in 1:5){
  score0<-c()
  for (i in 1:1000){
    if(i %in% c(8,112,159,303,399,477,478,640,848)){
      score0<-score0
    }else{
      score0<-c(score0,M4_y1[[i]]$score[j])
    }
  }
  score0<-mean(score0)
  score_output<-c(score_output,score0)
}

summery<-rbind(mase_err,smape_err,score_output)
colnames(summery)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
rownames(summery)<-c("mase_err","smape_err","log_score")
summery
```
# Test 4 Function log score
* optimal pool
```
y1<-0
y_score<-function(x){
  for (i in 1:(length(y))) {
    y1<-y1+log(x*p[i,1]+(1-x)*p[i,2])
  }
  return(-y1)
}
w_optim<-optimize(y_score,c(0,1),tol = 0.0001)
> w_optim
$minimum
[1] 0.4066993

$objective
[1] 731.5824
```
![logscore-w](/plot/logscore-w.png)
```
y<-M4_q1[[1]]$x
p<-PP[[1]]
features_y<-FF[[1]]
log_score1<-log(0.5*p[1,1]+0.5*p[1,2])+log(0.5*p[2,1]+0.5*p[2,2])+log(0.5*p[3,1]+0.5*p[3,2])+log(0.5*p[4,1]+0.5*p[4,2])+log(0.5*p[5,1]+0.5*p[5,2])+log(0.5*p[6,1]+0.5*p[6,2])+log(0.5*p[7,1]+0.5*p[7,2])+log(0.5*p[8,1]+0.5*p[8,2])+log(0.5*p[9,1]+0.5*p[9,2])
log_score<-function(beta){
  for (i in 10:(length(y))) {
    w<-1/(1+exp(-features_y[i-9,feature_select[[4]]]%*%beta))
    log_score1<-log_score1+log(w*p[i,1]+(1-w)*p[i,2])
  }
  return(-log_score1)
}

par(mfrow=c(1,1))
plot(log_score,xlim=c(-10,10),xlab="beta")
```
![logscore-beta](/plot/logscore-beta.png)

```
set.seed(2019-02-14)
w_max<-optim(fn=log_score,par=runif(1, 0, 0),method="SANN")
> w_max
$par
[1] -8.692382
$value
[1] 731.5224
$counts
function gradient 
   10000       NA 
$convergence
[1] 0
```
* 2 features
```
$par
[1] -18.232871  -6.426983
$value
[1] 731.5516
$counts
function gradient 
   10000       NA 
$convergence
[1] 0
```
* 3 features
```
$par
[1]  0.9803708 -3.7436639 -6.5726971
$value
[1] 731.5546
$counts
function gradient 
   10000       NA 
$convergence
[1] 0
```
* 4 features
```
$par
[1]   0.1111206   0.8026376 -19.9875921  -9.0262070
$value
[1] 731.5535
$counts
function gradient 
   10000       NA 
$convergence
[1] 0
```
* 5 features
```
$par
[1]   0.7042163  -8.1092401  -2.1315472 -14.9969504  -2.0160211
$value
[1] 731.5459
$counts
function gradient 
   10000       NA 
```
* 6 features
```
$par
[1]   0.3171425 -13.8511234  -4.7591751  -4.0357003 -14.6981076   4.5653223
$value
[1] 731.4784
$counts
function gradient 
   10000       NA 
$convergence
[1] 0

#
$par
[1] -0.9933166 -6.5409451 -2.4041121 -7.6890374 -2.7938049 -6.3608746  7.3838571
$value
[1] 731.4674
$counts
function gradient 
   10000       NA 
$convergence
[1] 0
```
