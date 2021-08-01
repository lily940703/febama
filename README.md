# `febama`: Feature-based Bayesian Forecasting Model Averaging

**Short Introduction**: In this work, we propose a novel framework for density forecast combination by constructing time-varying weights based on time series features, which is called FEature-based BAyesian forecasting Model Averaging (FEBAMA). Our framework estimates weights in the forecast combination via Bayesian log predictive scores, in which the optimal forecasting combination is determined by time-series features from historical information. In particular, we use an automatic Bayesian variable selection method to weight the importance of different features. To this end, our approach has better interpretability compared to other black-box forecasting combination schemes.  

## Installation
You can install the package `febama` from GitHub Repository with:
```
devtools::install_github("lily940703/febama")
```
## Usage
This part explains how to generate forecasts based on FEBAMA framework.

### Packages we need
```
library(febama)
library("tsfeatures")
library("M4metalearning")
library("forecast")
library("mvtnorm")
library(parallel)
library(doParallel)
library(foreach)
```

### Example data
```
# M3 monthly
library(Mcomp)
data_m3_mon <- Filter(function(l) l$period == "MONTHLY", M3)

# Pick a time series randomly
set.seed(2021-7-29)
id = sample(1:length(data_m3_mon), 1)
data_example = data_m3_mon[[id]]
```

### Training period
```
# 1. Training period
# 1.1 Compute predictive densities and features
model_conf_default = model_conf_default()
lpd_features = lpd_features_multi(data = data_example, model_conf = model_conf_default)

# If parallel
model_conf_curr = model_conf_default
model_conf_curr$lpd_features_parl = list(par = T, ncores = 2)
lpd_features = lpd_features_multi(data = data_example, model_conf = model_conf_curr)

# The features used
lpd_features <- feature_clean(list(lpd_features))[[1]]
fe <- lpd_features$feat
fm <- lpd_features$feat_mean
fs <- lpd_features$feat_sd
lpd_features$feat<- fe[, unique(unlist(model_conf_curr$features))]
lpd_features$feat_mean <- fm[unique(unlist(model_conf_curr$features))]
lpd_features$feat_sd <- fs[unique(unlist(model_conf_curr$features))]

# Up to now, we obtain log predictive densities and features for training.
head(lpd_features$lpd)
head(lpd_features$feat)

# 1.2 Inference procedure
# (1) For FEBAMA
parameters = febama_mcmc(data = lpd_features, model_conf = model_conf_curr)
# (2) For FEBAMA+VS
model_conf_curr$algArgs$nIter = 50
parameters_vs = febama_mcmc(data = lpd_features, model_conf = model_conf_curr)
```

### Forecasting period
```
# 2. Forecasting period
# (1) For FEBAMA
model_conf_curr$algArgs$nIter = 1
data_example_fore = forecast_feature_results_multi(ts = data_example, model_conf = model_conf_curr,
                                                   data = lpd_features, beta_out = parameters)

cat("The average LS and MASE of the example data based on FEBAMA method are \n", data_example_fore$err_feature[1:2])

# (2) For FEBAMA+VS
model_conf_curr$algArgs$nIter = 50
data_example_fore_vs = forecast_feature_results_multi(ts = data_example, model_conf = model_conf_curr,
                               data = lpd_features, beta_out = parameters_vs)
cat("The average LS and MASE of the example data based on FEBAMA+VS method are \n", data_example_fore_vs$err_feature[1:2])
```
