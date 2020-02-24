```
library(gratis)
library(tsfeatures)
library(M4metaresults)
library(M4comp2018)
library(M4metalearning)
library(forecast)
library(tseries)
install.packages("sgd")
library(sgd)
```
Randomly select 1000 m4 quarterly data
```
##data(m4 quarterly data)
M4_q<-M4[23001:47000]
set.seed(2019-12-31)
indices <- sample(length(M4_q))
M4_q1 <- M4_q[indices[1:1000]]
```
Compute the results of the five methods
* Our method (with features)
* John and Gianni's method
* Simple model averaging
* ARIMA
* ETS
```
for (a in 1:1000) {
## determine ets
ets1<-ets(M4_q1[[a]]$x,model="ZZZ")
ets_model<-ets1$components
ets_fore<-forecast(ets1,8)
M4_q1[[a]]$ets_fore<-ets_fore$mean

## determine arima
arima1<-auto.arima(M4_q1[[a]]$x)
d<-arimaorder(arima1)

## add predicted values of arima into M4_y1[[1]]
ari_fore<-forecast(arima1,8)
M4_q1[[a]]$ari_fore<-ari_fore$mean

y<-M4_q1[[a]]$x
p<-matrix(nrow = length(y),ncol = 2)

# p[t,1] of model ets
# 这个判断是否合理？
if("M" %in% ets_model){
  sigma_ets<-(sd(ets1$residuals))*(mean(ets1$fitted))
}else{
  sigma_ets<-(sd(ets1$residuals))
}

# PDF 是否合理？
for (t in 1:(length(y))) {
  p[t,1]<-dnorm(y[t],mean =ets1$fitted[t],sd=sigma_ets )
}

# p[t,2] of model arima
sigma_ari<-sd(arima1$residuals)
for (t in 1:(length(y))) {
  p[t,2]<-dnorm(y[t],mean =arima1$fitted[t],sd=sigma_ari )
}

# optimal solution of w
features_y<-c()
for (i in 2:length(y)) {
  myfeatures <- tsfeatures(y[1:i])
  myfeatures<-data.matrix(myfeatures)
  features_y<-rbind(features_y,myfeatures)
}
#令缺失值等于0
for (i in 1:length(features_y[,1])) {
  for (j in 1:length(features_y[1,])) {
    if(is.na(features_y[i,j])){
      features_y[i,j]<-0
    }
    else{
      features_y[i,j]<-features_y[i,j]
    }
  }
}

#权重为【t=1到上一期】特征的函数
log_score1<-0
log_score<-function(beta){
  for (i in 3:(length(y))) {
    w<-1/(1+exp(-features_y[i-2,]%*%beta))
    log_score1<-log_score1+log(w*p[i,1]+(1-w)*p[i,2])
  }
  return(-log_score1)
}
#optim默认求最小值
#SANN是一种模拟退火的方法，与通常的数学函数的算法不同，该算法是一种概率算法，对于各种复杂的情况，尤其是很多不可微的函数，该算法可以找到最优解，但效率上比不上其他数学算法。
set.seed(2019-01-09)
w_max<-optim(fn=log_score,par=runif(16, min = 0, max = 0),method="SANN")
beta_optim<-w_max$par

# 不加特征
y1<-0
y_score<-function(x){
  for (i in 1:(length(y))) {
    y1<-y1+log(x*p[i,1]+(1-x)*p[i,2])
  }
  return(-y1)
}
#plot(y_score,type = "l",xlim=c(0,1),xlab="w",ylab="log_score")
w_optim<-optimize(y_score,c(0,1),tol = 0.0001)
w_optim<-w_optim$minimum


#y_hat_output 8次预测结果
y_hat_feature<-c()
for (i in 1:8) {
  w<-1/(1+exp(-features_y[length(features_y[,1]),]%*%beta_optim))
  y_hat<-w*(M4_q1[[a]]$ets_fore[i])+(1-w)*(M4_q1[[a]]$ari_fore[i])
  y<-c(y,y_hat)
  features <- tsfeatures(y[1:length(y)])
  features<-data.matrix(features)
  features_y<-rbind(features_y,features)
  y_hat_feature<-c(y_hat_feature,y_hat)
   }

M4_q1[[a]]$y_hat_feature<-y_hat_feature
M4_q1[[a]]$y_hat_w<-w_optim*M4_q1[[a]]$ets_fore+(1-w_optim)*M4_q1[[a]]$ari_fore
M4_q1[[a]]$average<-0.5*(M4_q1[[a]]$ets_fore+M4_q1[[a]]$ari_fore)

#errors
ff<-rbind(M4_q1[[a]]$y_hat_feature,M4_q1[[a]]$y_hat_w,M4_q1[[a]]$average,M4_q1[[a]]$ari_fore,M4_q1[[a]]$ets_fore)
row.names(ff)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
M4_q1[[a]]$ff<-ff

lentry <- M4_q1[[a]]
insample <- lentry$x
ff<- lentry$ff
frq<-stats::frequency(insample)
insample <- as.numeric(insample)
outsample <- as.numeric(lentry$xx)
masep <- mean(abs(utils::head(insample, -frq) - utils::tail(insample, -frq)))
repoutsample <- matrix(rep(outsample, each = nrow(ff)), nrow = nrow(ff))
smape_err <- 200 * abs(ff - repoutsample)/(abs(ff) + abs(repoutsample))
mase_err <- abs(ff - repoutsample)/masep
M4_q1[[a]]$mase_err <- rowMeans(mase_err)
M4_q1[[a]]$smape_err <- rowMeans(smape_err)

#log predictive score
ari_score_sum<-0
for(j in 1:8){
  ari_score<-dnorm(outsample[j],mean =ff[4,j],sd=sigma_ari )
  browser()
  ari_score_sum<-ari_score_sum+log(ari_score)
} 
ets_score_sum<-0
for(j in 1:8){
  ets_score<-dnorm(outsample[j],mean =ff[5,j],sd=sigma_ets )
  ets_score_sum<-ets_score_sum+log(ets_score)
}  
ave_score_sum<-0
for(j in 1:8){
  ave_score<-0.5*dnorm(outsample[j],mean =ff[4,j],sd=sigma_ari)+0.5*dnorm(outsample[j],mean =ff[5,j],sd=sigma_ets)
  ave_score_sum<-ave_score_sum+log(ave_score)
}  
hat_w_score_sum<-0
for(j in 1:8){
  hat_w_score<-(1-w_optim)*dnorm(outsample[j],mean =ff[4,j],sd=sigma_ari)+w_optim*dnorm(outsample[j],mean =ff[5,j],sd=sigma_ets)
  hat_w_score_sum<-hat_w_score_sum+log(hat_w_score)
} 
hat_feature_score_sum<-0
for(j in 1:8){
  w<-1/(1+exp(-features_y[(length(features_y[,1])-9+j),]%*%beta_optim))
  hat_feature_score<-(1-w)*dnorm(outsample[j],mean =ff[4,j],sd=sigma_ari)+w*dnorm(outsample[j],mean =ff[5,j],sd=sigma_ets)
  hat_feature_score_sum<-hat_feature_score_sum+log(hat_feature_score)
} 
score<-c(hat_feature_score_sum,hat_w_score_sum,ave_score_sum,ari_score_sum,ets_score_sum)
names(score)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
M4_q1[[a]]$score<-score
}

```
Prediction performance of the data set 
* mase error
* smape error
* log score
```
mase_err<-c()
for (j in 1:5){
  mase_err0<-c()
  for (i in 1:1000){
  mase_err0<-c(mase_err0,M4_q1[[i]]$mase_err[j])
  }
  mase_err0<-mean(mase_err0)
  mase_err<-c(mase_err,mase_err0)
}


smape_err<-c()
for (j in 1:5){
  smape_err0<-c()
  for (i in 1:1000){
    smape_err0<-c(smape_err0,M4_q1[[i]]$smape_err[j])
  }
  smape_err0<-mean(smape_err0)
  smape_err<-c(smape_err,smape_err0)
}


score_output<-c()
for (j in 1:5){
  score0<-c()
  for (i in 1:1000){
    score0<-c(score0,M4_q1[[i]]$score[j])
  }
  score0<-mean(score0)
  score_output<-c(score_output,score0)
}

summery<-rbind(mase_err,smape_err,score_output)
colnames(summery)<-c("y_hat_feature","y_hat_w","average","ari_fore","ets_fore")
rownames(summery)<-c("mase_err","smape_err","log_score")
summery

```
# Results
```
> summery
          y_hat_feature    y_hat_w    average  ari_fore   ets_fore
mase_err       1.108739   1.109655   1.096166  1.116552   1.139252
smape_err      9.496338   9.441946   9.327331  9.569106   9.785150
log_score    -91.760041 -86.160908 -83.786251      -Inf -96.348629
```
# Discuss
In the article 173 data, the log score of arima is -Inf.
```
a1 <- arima1 %>% forecast(h = 8) %>% accuracy(M4_q1[[173]]$xx)
a1[,c("RMSE","MAE","MAPE","MASE")]
#                   RMSE       MAE      MAPE       MASE
# Training set  21.32889  17.74928 0.5567828 0.09182668
# Test set     555.05254 465.49995 7.1022832 2.40828394
a2 <- ets1 %>% forecast(h = 8) %>% accuracy(M4_q1[[173]]$xx)
a2[,c("RMSE","MAE","MAPE","MASE")]
#                 RMSE       MAE     MAPE      MASE
# Training set  35.17056  28.39818 0.888371 0.1469192
# Test set     219.44364 185.01128 2.825963 0.9571638
```
If delete M4_q1[[173]]
```
> summery
          y_hat_feature    y_hat_w    average   ari_fore   ets_fore
mase_err       1.107438   1.108356   1.096537   1.115259   1.139435
smape_err      9.498367   9.443921   9.334495   9.571208   9.792170
log_score    -91.394933 -85.989661 -83.677604 -97.209176 -96.258060
```

