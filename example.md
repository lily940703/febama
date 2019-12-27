# A simple example of feature-based bayesian model averaging
* model 1  AR(1)  
* model 2  random walk 

### a ts
```
x <- generate_ts(n.ts = 1, freq = 1, nComp = 2, n = 20)
y<-x$N1$x
#加初值
y<-c(0,y)
y<-c(0.00000,21.72881,26.80634,28.72464,33.37096,34.12180,40.02512,
43.12906,42.79923,43.28925,47.78483,50.43229,52.01138,56.37153,59.94831,
64.44778,64.94955,69.85086,72.23654,78.41327,80.10752)
```
### matrix p(20*2)
```
#p为n*2矩阵，n为ts长度
p<-matrix(nrow = 20,ncol = 2)
#t=1的p11，p12
p[1,1]<-dnorm(y[2], 0, 10)
p[1,2]<-dnorm(y[2], 0, 10)
#postrior_seta
likelihood<-1
postrior_seta<-function(seta,t){
  prior<-dnorm(seta, 0, 1)
  for (i in 1:(t-1)) {
    likelihood<-likelihood*dnorm(y[i+1], seta*y[i], 10)
    i<-i+1
  }
  postrior<-prior*likelihood
  return(postrior)
}

#proposal
#产生1个均值是x，标准差为1的随机数
proposal_seta <- function(x) rnorm(1, x, 1)

step_mcmc <- function(x, postrior, proposal,t) {
  ## Pick new point
  xp <- proposal(x)
  ## Acceptance probability:
  alpha <- min(1, postrior(xp,t) / postrior(x,t))
  # test<-postrior(xp,t) / postrior(x,t)
  # cat()
  ## Accept new point with probability alpha:
  #runif(1)默认生成1个[0,1]上的均匀分布随机数
  if (runif(1) < alpha)
    x <- xp
  ## Returning the point:
  x
}

run_mcmc <- function(x, postrior, proposal, nsteps,t) {
  res <- matrix(NA, nsteps, length(x))
  for (i in seq_len(nsteps))
    res[i,] <- x <- step_mcmc(x, postrior, proposal,t)
  return(res)
}
# p[t,1] of model A1
for (t in 2:(length(y)-1)) {
  res <- run_mcmc(1, postrior_seta, proposal_seta, 4000,t=t)
  res<-unique(res[-(1:100),])
  if(length(res)<20){
    warning("samples of seta < 20")
  }
  res<-sample(res, 20)
  #acceptance <- 1-mean(duplicated(res))
  p_sum<-0
  for (i in 1:length(res)) {
    p_sum<-p_sum+dnorm(y[t+1], res[i]*y[t], 4)
  }
  p[t,1]<-p_sum/length(res)
}

# p[t,2] of model A2
for (t in 2:(length(y)-1)) {
  p[t,2]<-dnorm(y[t+1], y[t], 10)
}
```
### optimal solution of w
```
# log_score<-function(w){
#   for (i in 1:(length(y)-1)) {
#     log_score1<-log_score1+log(w*p[i,1]+(1-w)*p[i,2])
#   }
#   return(log_score1)
# }

features_y<-c()
for (i in 4:21) {
  myfeatures <- tsfeatures(y[2:i])
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
#前两期features都为0，这部分没有进行优化，权重都是0.5
x1<-runif(16,min=0,max=0)
features_y<-rbind(x1,features_y)
features_y<-rbind(x1,features_y)
#权重为【上一期】特征的函数
log_score1<-log(0.5*p[1,1]+(1-0.5)*p[1,2])
log_score<-function(beta){
  for (i in 2:(length(y)-1)) {
    w<-1/(1+exp(-features_y[i-1,]%*%beta))
    log_score1<-log_score1+log(w*p[i,1]+(1-w)*p[i,2])
  }
  return(-log_score1)
}
#optim默认求最小值
#SANN是一种模拟退火的方法，与通常的数学函数的算法不同，该算法是一种概率算法，对于各种复杂的情况，尤其是很多不可微的函数，该算法可以找到最优解，但效率上比不上其他数学算法。
set.seed(2019-12-19)
w_max<-optim(fn=log_score,par=runif(16, min = 0, max = 0),method="SANN")

```
```
# > w_max
# $par
# [1]   4.824908   8.850448  19.520589  -3.942189  17.707244   3.119714   2.685472   3.057554  -7.766437   8.456901
# [11]  -6.119628  19.859726 -26.134054  12.350602   1.344076  18.172853
# 
# $value
# [1] 293.1942
# 
# $counts
# function gradient 
# 10000       NA 
# 
# $convergence
# [1] 0
# 
# $message
# NULL
```
### forecast
```
proposal_y <- function(x) rnorm(1, x, 10)
beta_optim<-w_max$par
#每步预测权重需要随特征变化更新
density_y<-function(y1,seta,y_old){
  w<-1/(1+exp(-features_y[length(features_y[,1]),]%*%beta_optim))
  density_y<-w*dnorm(y1, seta*y_old, 10)+(1-w)*dnorm(y1, y_old, 10)
  return(density_y)
}
step_mcmc_y <- function(x, postrior, proposal,seta,y_old) {
  ## Pick new point
  xp <- proposal(x)
  ## Acceptance probability:
  alpha <- min(1, postrior(xp,seta,y_old) / postrior(x,seta,y_old)) 
  ## Accept new point with probability alpha:
  #runif(1)默认生成1个[0,1]上的均匀分布随机数
  if (runif(1) < alpha)
    x <- xp
  ## Returning the point:
  x
}
run_mcmc_y <- function(x, postrior, proposal, nsteps,seta,y_old) {
  res <- matrix(NA, nsteps, length(x))
  for (i in seq_len(nsteps))
    res[i,] <- x <- step_mcmc_y(x, postrior, proposal,seta,y_old)
  return(res)
}


#y_hat_res存6次预测的抽样结果矩阵
y_hat_res<-list()
#y_hat_output 6次预测结果
y_hat_output<-c()
y<-c(0.00000,21.72881,26.80634,28.72464,33.37096,34.12180,40.02512,43.12906,42.79923,43.28925,47.78483,50.43229,52.01138,56.37153,59.94831,64.44778,64.94955,69.85086,72.23654,78.41327,80.10752)

for (i in 1:6) {
  y_hat_all<-c()
  t<-20+i
  seta_f <- run_mcmc(1, postrior_seta, proposal_seta, 5000,t=t)
  #plot(seta_f,type = 'l')
  seta_f<-unique(seta_f[-(1:100),])
  if(length(seta_f)<20){
    warning("samples of seta < 20")
  }
  seta_f<-sample(seta_f, 20)
  for (j in 1:length(seta_f)) {
    seta<-seta_f[j]
    y_old<-y[length(y)]
    y_hat <- run_mcmc_y(y[length(y)], density_y, proposal_y, 5000, seta = seta, y_old = y_old)
    y_hat<-unique(y_hat[-(1:100),])
    if(length(y_hat)<20){
      warning("samples of y_hat < 20")
    }
    y_hat<-sample(y_hat, 20)
    y_hat_all<-rbind(y_hat_all,y_hat)
  }  
  y_hat_all_seta<-cbind(seta_f,y_hat_all)
  rownames(y_hat_all_seta)<-NULL
  y_hat_res[[i]]<-y_hat_all_seta
  y<-c(y,mean(y_hat_all_seta[,2:21]))
  features <- tsfeatures(y[2:length(y)])
  features<-data.matrix(features)
  features_y<-rbind(features_y,features)
  y_hat_output<-c(y_hat_output,mean(y_hat_all_seta[,2:21]))
}
y_hat_output
# [1]  82.82455  87.38487  91.34578  97.56553 101.35953 106.73193

```
