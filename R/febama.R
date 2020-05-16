forecast_feature_results_multi <-function(data, model_conf, intercept = T,
                                          lpd_features, beta_pre)
{
    ##attach(model_conf)
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

    y <- data$x
    y01 = scale(y, center = TRUE, scale = TRUE)
    y_mean = attr(y01, "scaled:center")
    y_sd = attr(y01, "scaled:scale")
    y01 = as.numeric(y01)
    y_true = data$xx
    y01_true = as.numeric(scale(y_true, center = y_mean, scale = y_sd))
    y_new = y01
    y_new_nonsd = as.numeric(y)

    features_y = lpd_features$feat
    features_y_mean = lpd_features$feat_mean
    features_y_sd = lpd_features$feat_sd

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
            ##用非标准化的数据计算特征
            myts <- list(list(x = ts(y_new_nonsd1, frequency = frequency)))
            myfeatures <- THA_features(myts)[[1]]$features
            myfeatures <- data.matrix(myfeatures)
            myfeatures <- myfeatures[, colnames(myfeatures) %in% colnames(features_y)]
            myfeatures_scaled = scale(t(myfeatures),
                                      center = features_y_mean, scale = features_y_sd)
        } else
        {
            myfeatures_scaled = NULL
        }

        ## Update predictive weights
        w_get <- function(beta_pre, myfeatures_scaled){
            myfeatures_scaled <- myfeatures_scaled[, beta_pre$features_select]
            if (intercept){
                myfeatures_scaled = cbind(1, t(myfeatures_scaled))
            }
            exp_lin = exp(matrix(myfeatures_scaled, nrow = 1) %*% beta_pre$beta)
            ##  avoid Inf
            exp_lin[exp_lin > exp(709)] <- exp(709)
            w <- exp_lin / (1 + rowSums(exp_lin))
            w_full = cbind(w, 1 - rowSums(w))
            return(w_full)
        }

        w_full <- sapply(beta_pre, w_get, myfeatures_scaled = myfeatures_scaled)
        w_full_all[[t]] <- w_full
        w_full_mean <- rowMeans(w_full)
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

        y_pred_multi <- sum (w_full_mean * (sapply(multi_fore, function(mean_sd){
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

    data$ff_feature <- y_hat_matrix
    colnames(w_full_mean_h) <- seq(1, forecast_h, 1)
    rownames(w_full_mean_h) <- unlist(fore_model)
    data$w_time_varying <- w_full_mean_h
    data$w_detail <- w_full_all

    ## mase smape
    ff <- data$ff_feature
    insample <- as.numeric(data$x)
    frq <- stats::frequency(insample)
    outsample <- as.numeric(data$xx)
    masep <- mean(abs(utils::head(insample,-frq) - utils::tail(insample,-frq)))
    smape_err <- 200 * abs(ff - outsample) / (abs(ff) + abs(outsample))
    mase_err <- abs(ff - outsample) / masep
    mase_err_h <- rowMeans(mase_err)
    smape_err_h <- rowMeans(smape_err)

    data$err_feature <- cbind(lpds, mase_err_h, smape_err_h)
    return(data)
}

forecast_feature_performance<-function(data)
{
    ## log score
    performance<-c()
    for (i_ts in 1:length(data)) {
        performance<-rbind(performance,data[[i_ts]]$err_feature)
    }
    performance_out1<- sum(performance[,1])

    performance_out2<- colMeans(performance[,2:3])

    performance_out<-cbind(performance_out1,t(performance_out2))
    colnames(performance_out)<-c("Log score","Mase","Smape")
    return(performance_out)
}
