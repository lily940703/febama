#' Defualt parameter settings of FEBAMA framework
#' 
#' @export

model_conf_default <- function(){
  num_models = 4
  parameter_list = list(
    frequency = 12 # Monthly data
    , ets_model = "AAN"
    , ts_scale = T # Whether the time series needs standardization.
    , forecast_h = 18 # Forecasting horizon
    , train_h = 1  # Out of sample training, 1 means rolling forecast
    , history_burn = 25 # Let the model to start with at least this length of historical data.
    , PI_level = 90 # Predictive Interval level, used to extract out-of-sample variance from forecasting models.
    , roll = NULL # The length of rolling samples, larger than history_burn.
    , feature_window = NULL # The length of moving window when computing features
    , features_used = rep(list(c("x_acf1","diff1_acf1", "entropy", 
                                 "alpha", "beta", "unitroot_kpss")), num_models - 1)
    , fore_model = c("ets_fore",  "naive_fore", "rw_drift_fore","auto.arima_fore")
    , lpd_features_parl = list(par = F, ncores = 1) # Whether parallel when computing predictive densities and features.
    
    ## Variable selection settings. By default, every model shares the same
    ## settings. Otherwise, write the full list, same applies to priArgs, this would allow
    ## for some models with only intercept.  Variable selection candidates, NULL: no
    ## variable selection use the full covariates provided by $init. ("all-in", "all-out",
    ## "random", or user-input)
    , varSelArgs = rep(list(list(cand = "2:end", init = "all-in")), num_models - 1)
    
    , priArgs = rep(list(list("beta" = list(type = "cond-mvnorm",
                                            mean = 0, covariance = "identity", shrinkage = 10),
                              "betaIdx" = list(type = "beta", alpha0 = 1, beta0 = 1))), num_models - 1)
    
    , algArgs = list(initOptim = TRUE, 
                     algName = "MAP", 
                     nIter = 1, 
                     "sgld" = list(max_batchSize = 108,
                                   nEpoch = 10,
                                   burninProp = 0.4, 
                                   stepsize = 0.1,
                                   gama = 0.55,
                                   a = 0.4,
                                   b = 10)
    )
    
  )
  parameter_list
}

