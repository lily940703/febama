# Default model config template 

```{r}
num_models = 6 # The number of forecasting models used in our method

model_conf_default = list(
  frequency = 12
  , ets_model = "ANN" # Simple exponential smoothing with additive errors
  , forecast_h = 18 # Forecasting horizon
  , train_h = 1 # Out of sample training, 1 means rolling forecast
  , history_burn = 36 # Let the model to start with at least this length of historical data.
  , PI_level = 90 # Predictive Interval level, used to extract out-of-sample variance from forecasting models.
  , roll = NULL # The length of historical data if use rolling samples, NULL means with full samples. 
  , feature_window = NULL # The length of moving window when computing features, NULL means using full historical data.
  , features_used = rep(list(c("entropy", "arch_acf", "alpha", "beta", "unitroot_kpss")), num_models - 1) # Features taken into consideration.
  , fore_model = c("ets_fore",  "naive_fore", "rw_drift_fore","auto.arima_fore",
                   "snaive_fore", "stlm_ar_fore") # Forecasting models used in our method.
  , lpd_features_parl = list(par = F, ncores = 1) # Parallel setup when computing log predictive density and features, par = F means no parallel.
  ## Variable selection settings. By default, every model shares the same
  ## settings. Otherwise, write the full list, same applies to priArgs, this would allow
  ## for some models with only intercept.  Variable selection candidates, NULL: no
  ## variable selection use the full covariates provided by $init. ("all-in", "all-out",
  ## "random", or user-input)
  , varSelArgs = rep(list(list(cand = "2:end", init = "all-in")), num_models - 1)
  
  , priArgs = rep(list(list("beta" = list(type = "cond-mvnorm",
                                          mean = 0, covariance = "identity", shrinkage = 10),
                            "betaIdx" = list(type = "beta", alpha0 = 1, beta0 = 1))), num_models - 1)
  
  , algArgs = list(initOptim = TRUE, # Use LBFGS to optimize initial values
                   algName = "sgld", 
                   nIter = 1, # Number of iterations, if nIter = 1, only use the initial values (obtained by MAP) of beta.
                   "sgld" = list(max_batchSize = 108,
                                 nEpoch = 10,
                                 burninProp = 0.4, # burnin proportion within each SGLD.
                                 stepsize = 0.1,
                                 gama = 0.55,
                                 a = 0.4,
                                 b = 10) Default parameter setting in SGLD algorithm.
  )
  
)
```
