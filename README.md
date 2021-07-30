# `febama`: Feature-based Bayesian Forecasting Model Averaging

**Short Introduction**: In this work, we propose a novel framework for density forecast combination by constructing time-varying weights based on time series features, which is called FEature-based BAyesian forecasting Model Averaging (FEBAMA). Our framework estimates weights in the forecast combination via Bayesian log predictive scores, in which the optimal forecasting combination is determined by time-series features from historical information. In particular, we use an automatic Bayesian variable selection method to weight the importance of different features. To this end, our approach has better interpretability compared to other black-box forecasting combination schemes.  

## Installation
You can install the package fuma from GitHub Repository with:
```
devtools::install_github("lily940703/febama")

```
