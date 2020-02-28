### Update code of gen_train_data.R

### Some queations
* 1 Standardization of time series and features  
```
 y <- data[[i_ts]]$x

  y01 = scale(y, center = TRUE, scale = TRUE)
  y_mean = attr(y01, "scaled:center")
  y_sd = attr(y01, "scaled:scale")
  y01 = as.numeric(y01)

  y_true = data[[i_ts]]$xx
  y01_true = as.numeric(scale(y_true, center = y_mean, scale = y_sd))

```
```
if(!is.null(features_y))
    {
      myts <- list(list(x=ts(y_new, frequency = frequency)))
      myfeatures <- THA_features(myts)[[1]]$features
      myfeatures <- data.matrix(myfeatures)
      myfeatures_scaled = scale(myfeatures, center = features_y_mean, scale = features_y_sd)
      # null in myfeatures_scaled
      myfeatures_scaled[is.na( myfeatures_scaled)] <- 0
      ## features_y_hat[t, ] <- myfeatures_scaled
    } else
    {
      myfeatures_scaled = NULL
    }
```
* 2 The formula for calculating weights
```
log_score<-function(beta, features, prob, intercept){

  if(intercept) features = cbind(rep(1, nrow(prob)), features)

  exp_lin = exp(features%*%beta)

  w <- exp_lin/(1+rowSums(exp_lin)) # T-by-(n-1)

  ## Full (T-by-n) matrix. To keep identification, only first (n-1) are connected with
  ## features. TODO: Common features?
  w_full = cbind(w, 1 - rowSums(w)) # T-by-n

  out = sum(log(rowSums(w_full * prob)))
  return(out)
}
```

### Some tests  

* 1 one step prediction (optimal pool, SA, ETS, ARIMA) 
```
> print(data.frame(lpds = lpds, lpds_simple = lpds_simple,
+                  lpds_ets = lpds_ets , lpds_arima = lpds_arima ))

       lpds lpds_simple  lpds_ets lpds_arima
1 0.9384428   0.2625277 -1.154993  0.9978105

> optimal_w
 [1] 9.999999e-01 5.266312e-01 3.438070e-01 6.713627e-01 1.278169e-01 3.026081e-01 6.436414e-09 9.999996e-01
 [9] 5.147197e-02 3.322197e-01
```
* 2 one step prediction (feature-based) 
```
w_max<-try(optim(fn=log_score,
                  par=runif(43, min = 0, max = 0),
                  features = features_y,
                  prob = exp(y_lpd),
                  intercept = TRUE,
                  method="SANN",
                  control = list(fnscale = -1))              )
  if(is(w_max, "try-error")) browser()
  if(w_max$convergence!=0){
    cat("The optimization does not converge in data", a)
  }
  beta_optim <- w_max$par
```
```
> print(data.frame(lpds_feature = lpds))
  lpds_feature
1     1.135728

> feature_w
 [1]  1.806752e-51  1.000000e+00  1.000000e+00  1.000000e+00  7.775492e-18  6.913349e-31  1.000000e+00
 [8] 7.137164e-114  1.676356e-43  8.341157e-51
 
 > print(data.frame(lpds_feature = lpds))
  lpds_feature
1   -0.6459095

> feature_w
 [1]  1.377138e-33  9.513313e-03  3.215037e-51  1.758064e-63  1.555745e-50  6.886911e-14  3.401801e-13
 [8] 5.249988e-134  1.000000e+00  1.000000e+00
```
The optimization is convergent, but the optimal solution obtained each time is absolutely different, resulting in great fluctuation of the calculated weights.  
So, when adding so much features, the current method is obviously ineffective.

* 3 Eight-step prediction (optimal pool, SA, ETS, ARIMA) 
```
> print(data.frame(lpds = lpds, lpds_simple = lpds_simple,
+                  lpds_ets = lpds_ets , lpds_arima = lpds_arima ))

       lpds lpds_simple  lpds_ets lpds_arima
1 -25.16407   -42.45068 -56.53136  -37.76633
```
### Plan for the following week
*  Expanding the number of data to verify the advantage of the method of optimal pool (including mase and smape error).  
**Then parallel is needed.** How to simplify the fit and forecast procedure?   
*  Add features gradually.  
**Then change to a better optimization.**

