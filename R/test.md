### Update code of gen_train_data.R

### Some queations

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

