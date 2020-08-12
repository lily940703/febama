#' Forecasting with our method.
#'
#' @param ts: A list (the time series to be predicted})
#' \describe{
#'   \item{x}{historical data}
#'   \item{xx}{real data in the forecasting horizon}
#' }
#' @param model_conf
#' @param data:A list 
#' \describe{
#'   \item{lpd}{Log probability densities }
#'   \item{feat}{Features}
#'   \item{feat_mean}{the mean of features from historical data}
#'   \item{feat_sd}{the sd of features from historical data}
#' }
#' @param beta_out: the output of function \code{febama_mcmc}
#' @return A list  
#' \describe{
#'   \item{err_feature}{log score, mase, smape measure of forecasts}
#'   \item{ff_feature}{forecasts in the forecasting horizon}
#'   \item{w_time_varying}{time-varying weights in the forecasting horizon}
#' }
#' export

forecast_feature_results_multi <-function(ts, model_conf = model_conf_curr, intercept = T,
                                          data, beta_out)
{
    feature_window = model_conf$feature_window
    roll = model_conf$roll
    frequency = model_conf$frequency
    history_burn = model_conf$history_burn
    ets_model = model_conf$ets_model
    forecast_h = model_conf$forecast_h
    train_h = model_conf$train_h
    PI_level = model_conf$PI_level
    fore_model = model_conf$fore_model


    ## forecasting
    y_hat_matrix <- matrix(ncol = forecast_h, nrow = 1)

    y <- ts$x
    y01 = scale(y, center = TRUE, scale = TRUE)
    y_mean = attr(y01, "scaled:center")
    y_sd = attr(y01, "scaled:scale")
    y01 = as.numeric(y01)
    y_true = ts$xx
    y01_true = as.numeric(scale(y_true, center = y_mean, scale = y_sd))
    y_new = y01
    y_new_nonsd = as.numeric(y)
    
    if(length(model_conf$features_used[[1]]) == 0){
        features_y = NULL
    }else{
        features_y = data$feat
    }
    
    features_y_mean = data$feat_mean
    features_y_sd = data$feat_sd

    lpds = 0
    pred_densities = matrix(NA, forecast_h, length(fore_model))
    colnames(pred_densities) <- unlist(fore_model)
    w_full_mean_h <- c()

    w_full_all <- list()
    for (t in 1:forecast_h)
    {
        ## Update features
        if(!is.null(feature_window)){
            y_new_nonsd1 <- tail(y_new_nonsd, feature_window)
        }else{
            y_new_nonsd1 <- y_new_nonsd
        }

        if (!is.null(features_y))
        {
            myts <- list(list(x = ts(y_new_nonsd1, frequency = frequency)))
            myfeatures <- M4metalearning::THA_features(myts)[[1]]$features
            myfeatures <- data.matrix(myfeatures)
            myfeatures <- myfeatures[, colnames(myfeatures) %in% colnames(features_y)]
            myfeatures_scaled = scale(t(myfeatures),
                                      center = features_y_mean, scale = features_y_sd)
        } else
        {
            myfeatures_scaled = NULL
        }

        if(model_conf$algArgs$nIter > 1){
            w_full <- w_get(beta_out, myfeatures_scaled, model_conf = model_conf)
        }else if(model_conf$algArgs$nIter == 1){
            w_full <- w_get_ini(beta_out, myfeatures_scaled, model_conf = model_conf)
        }
        
        w_full_all[[t]] <- w_full
        w_full_mean <- colMeans(w_full)
        w_full_mean_h <- cbind(w_full_mean_h, w_full_mean)

        ## forecast
        if(is.null(roll)){
            y_new1 <- y_new
        }else{
            y_new1 <- tail(y_new, roll)
        }

        multi_fore <- lapply(fore_model, function(method){
            method_fun <- get(method)
            mean_sd <- method_fun (y_new1, train_h, PI_level)
            return(mean_sd)
        })

        y_pred_multi <- sum (t(w_full_mean) * (sapply(multi_fore, function(mean_sd){
            return(mean_sd[[1]])
        })))

        y_new = c(y_new, y_pred_multi)
        y_new_nonsd = c(y_new_nonsd, (y_pred_multi * y_sd + y_mean))
        y_hat_matrix[1, t] <- y_pred_multi * y_sd + y_mean

        ## The predictive log score
        pd_multi <- sapply(multi_fore, function(mean_sd){
            dnorm(y01_true[t], mean = mean_sd[[1]], sd = mean_sd[[2]], log = F)
        })
        lpds_multi = log(sum(pd_multi * w_full_mean))

        if(lpds_multi < log(1e-323)){
            lpds_multi = log(1e-323)
        }else{
            lpds_multi = lpds_multi
        }
        lpds = lpds + lpds_multi
    }

    ts$ff_feature <- y_hat_matrix
    colnames(w_full_mean_h) <- seq(1, forecast_h, 1)
    rownames(w_full_mean_h) <- unlist(fore_model)
    ts$w_time_varying <- w_full_mean_h
    ts$w_detail <- w_full_all

    ## mase smape
    ff <- ts$ff_feature
    insample <- as.numeric(ts$x)
    frq <- stats::frequency(insample)
    outsample <- as.numeric(ts$xx)
    masep <- mean(abs(utils::head(insample,-frq) - utils::tail(insample,-frq)))
    smape_err <- 200 * abs(ff - outsample) / (abs(ff) + abs(outsample))
    mase_err <- abs(ff - outsample) / masep
    mase_err_h <- rowMeans(mase_err)
    smape_err_h <- rowMeans(smape_err)

    ts$err_feature <- cbind(lpds, mase_err_h, smape_err_h)
    return(ts)
}

#' Compute the weights of models based on beta(after optimization by MCMC) and features
#' export
w_get <- function(beta_out, myfeatures_scaled, model_conf){
    
    num_models_updated = length(model_conf$fore_model)-1
    exp_lin <- matrix(nrow = model_conf$algArgs$nIter-1, ncol = num_models_updated+1)
    for(iComp in 1:num_models_updated)
    {
        betaCurr = beta_out$beta[[iComp]]
        betaIdxCurr = beta_out$betaIdx[[iComp]]
        features_used_curr = model_conf$features_used[[iComp]]
        features0 = cbind(1, myfeatures_scaled[, features_used_curr, drop = FALSE])
        
        me = features0 %*% t(betaCurr[-1,])
        me[me>709] <- 709 # avoid overflow
        exp_lin[, iComp] = exp(me)
    }
    
    exp_lin[, num_models_updated + 1] = 1 # assume last model is 1
    weights = exp_lin/ rowSums(exp_lin) # T-by-n, assuming last is deterministic

    return(weights)
}

#' Compute the weights of models based on beta and features
#' If only use the initial value (MAP) in \code{febama_mcmc}
#' export
w_get_ini <- function(beta_out, myfeatures_scaled, model_conf){
    
    num_models_updated = length(model_conf$fore_model)-1
    exp_lin <- matrix(nrow = 1, ncol = num_models_updated+1)
    for(iComp in 1:num_models_updated)
    {
        betaCurr = beta_out$beta[[iComp]]
        betaIdxCurr = beta_out$betaIdx[[iComp]]
        features_used_curr = model_conf$features_used[[iComp]]
        features0 = cbind(1, myfeatures_scaled[, features_used_curr, drop = FALSE])
        
        me = features0 %*% t(betaCurr)
        me[me>709] <- 709 # avoid overflow
        exp_lin[, iComp] = exp(me)
    }
    
    exp_lin[, num_models_updated + 1] = 1 # assume last model is 1
    weights = exp_lin/ rowSums(exp_lin) # T-by-n, assuming last is deterministic
    
    return(weights)
}

#' Output comprehensive forecast performance (average)
#' @param data: A list of multiple ts with \item{err_feature}
#' @return A vector of "Log score","Mase" and "Smape".
#' export
forecast_feature_performance<-function(data)
{
    ## log score
    performance<-c()
    for (i_ts in 1:length(data)) {
        performance<-rbind(performance,data[[i_ts]]$err_feature)
    }

    performance_out<- t(colMeans(performance))

    colnames(performance_out)<-c("Log score","Mase","Smape")
    return(performance_out)
}
