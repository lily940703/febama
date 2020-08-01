# `febama`: Feature-based Bayesian Forecasting Model Averaging

**Short Introduction**: In this work, we propose a feature-based Bayesian forecasting model averaging framework (_febama_). Our Bayesian framework estimates weights of the feature-based forecasting combination via a Bayesian log predictive score, in which the optimal forecasting combination is connected and determined by time-series features from historical information. In particular, we utilize the prior knowledge of the coefficients of time-series features. We use an efficient Bayesian variable selection method to weight important features that may affect the forecasting combinations. To this end, our approach has better interpretability compared to other black-box forecasting combination schemes. Our framework is more computational efficient because the log predictive score and time-series features are calculated in the offline phase. We apply our framework to stock market data and M4 competition data. Based on our structure, a simple maximum-a-posteriori scheme outperforms the optimal prediction pools (Geweke and Amisano, 2011) or simple averaging, and Bayesian variable selection further enhanced the forecasting performance. 

## TODO

- [ ] Do we need to have full cumulative time series features? We could just calculate the
      feature with a moving window.
- [ ] Write a gradient function for the log score
- [ ] Replace optimization function with SGD
- [ ] Parallel the code
- [ ] Package the code
