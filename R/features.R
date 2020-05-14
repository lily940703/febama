#' Compute log predictive density and features
#' @param data a ts data with $x
#' @param model_conf model conf
#' @return a list including two matrices of log probability densities and features.
#' @author Feng Li
#' @export
lpd_feature_multi <- function(data, model_conf) {

    feature_window = model_conf$feature_window
    roll = model_conf$roll
    frequency = model_conf$frequency
    history_burn = model_conf$history_burn
    ets_model = model_conf$ets_model
    forecast_h = model_conf$forecast_h
    train_h = model_conf$train_h
    PI_level = model_conf$PI_level
    fore_model = model_conf$fore_model

    ## A single historical data
    y <- data$x
    y1 = scale(y, center = TRUE, scale = TRUE)
    y_mean = attr(y1, "scaled:center")
    y_sd = attr(y1, "scaled:scale")

    y1 = as.numeric(y1)

    ## Calculate historical log predictive density
    num_model <- length(fore_model)
    log_pred_densities <-
        matrix(nrow = length(y) - history_burn - train_h + 1,
               ncol = num_model)
    colnames(log_pred_densities) <- unlist(fore_model)

    for (t in (history_burn):(length(y) - train_h))
    {
        ## TODO: We may simplify this to assume the fit and forecast procedure is
        ## invariant with certain time period to save time.
        if(is.null(roll)){
            y01 <- y1[1:t]
        }else if (t < roll){
            y01 <- y1[1:t]
        }else{
            y01 <- y1[(t-roll+1):t]
        }

        ## To keep numeric stability, we calculate log P(y_pred)
        use_model <- lapply(fore_model, function(method){
            method_fun <- get(method)
            mean_sd <- method_fun (y01, train_h, PI_level)
            return(mean_sd)
        })
        log_pred_den <- lapply(use_model, function(mean_sd){
            lpd <- sum(dnorm(y1[(t + 1):(t + train_h)],
                             mean = mean_sd[[1]],sd = mean_sd[[2]],log = TRUE ))
            return(lpd)
        })
        log_pred_den <- as.numeric(log_pred_den)
        log_pred_den[log_pred_den < log(1e-323)] <- log(1e-323)
        log_pred_den[log_pred_den > log(1e+308)] <- log(1e+308)
        log_pred_densities[(t - history_burn + 1), ] <- log_pred_den
    }

    ## Calculate historical features
    features_y <- matrix(nrow = length(y) - train_h - history_burn + 1,
                         ncol = 42)
    myts <- list(list(x = ts(y[1:history_burn], frequency = frequency)))
    myfeatures <- THA_features(myts)[[1]]$features
    names <- colnames(myfeatures)
    colnames(features_y) <- names

    if(is.null(feature_window)){
        for (t in (history_burn):(length(y) - train_h))
        {
            myts <- list(list(x = ts(y[1:t], frequency = frequency)))
            myfeatures <- THA_features(myts)[[1]]$features
            myfeatures <- data.matrix(myfeatures)
            features_y[(t - history_burn + 1),] <- myfeatures
        }
    }
    else{
        features_y <- features_window(y, window = feature_window,
                                      history_burn, train_h)
    }

    features_y_scaled = scale(features_y, center = TRUE, scale = TRUE)

    lpd_feature <- list(lpd = log_pred_densities, feat = features_y_scaled)
    return(lpd_feature)
}

features_window <- function(y, window, history_burn, train_h){
    if(window < history_burn){
        stop("The length of window needs to be larger than history_burn")
    }
    features_y<-c()
    for (t in (history_burn):(length(y) - train_h))
    {
        if(t <= window){
            myts <-list(list(x=ts(y[1:t], frequency = 1)))
            myfeatures <- THA_features(myts)[[1]]$features
            myfeatures <- data.matrix(myfeatures)
        }else{
            myts <-list(list(x=ts(y[(t-window+1):t], frequency = 1)))
            myfeatures <- THA_features(myts)[[1]]$features
            myfeatures <- data.matrix(myfeatures)
        }
        features_y<-rbind(features_y,myfeatures)
    }
    return(features_y)
}



#' Delete the features with NaN and add attributes
#'
#' @param lpd_feature
#' @return final features
#' @author Feng Li
#' @export
feature_clean <- function(lpd_feature){
    for (i in 1:length(lpd_feature)) {
        ind <- which(is.nan(lpd_feature[[i]]$feat), arr.ind = T)[,2]
        ind0 <- ind[!duplicated(ind)]
        lpd_feature[[i]]$feat_mean <- attr(lpd_feature[[i]]$feat, "scaled:center")[-ind0]
        lpd_feature[[i]]$feat_sd <- attr(lpd_feature[[i]]$feat, "scaled:scale")[-ind0]
        lpd_feature[[i]]$feat <- lpd_feature[[i]]$feat[, -ind0]
    }
    return(lpd_feature)
}
