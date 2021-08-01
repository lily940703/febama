#' Compute log predictive densities and features
#' 
#' Compute log predictive densities and features for training in FEBAMA framework.
#' 
#' @param data A list with related information of a time series. Historical data \code{x} is neccessary.
#' @param model_conf Parameter settings of FEBAMA framework. Defualt \code{model_conf_default()}.
#' 
#' @return A list with log predictive densities and features.
#' @export
lpd_features_multi <- function(data, model_conf) {

    y <- data$x
    if (model_conf$ts_scale == T){
        y1 = scale(y, center = TRUE, scale = TRUE)
        y_mean = attr(y1, "scaled:center")
        y_sd = attr(y1, "scaled:scale")
        y1 = as.numeric(y1)
    }else{
        y1 = as.numeric(y)
    }
    
    ## if some subseries are constant, reset history_burn to
    ## ensure the subseries can be scaled when computing features
    burn = 0
    for (t in 1:length(y1)) {
        if(y1[t] == y1[1]){ 
            burn = burn+1
        }else{
                break;
            }
    }
    if(burn >= model_conf$history_burn -1){
        history_burn = burn+2
    }else{
        history_burn = model_conf$history_burn
    }
    
    train_h = model_conf$train_h
    
    if(model_conf$lpd_features_parl$par == F){
        lpd_features = lpd_feat(t_seq = c((history_burn):(length(y) - train_h)), 
                                ts_sd = y1, ts_nosd = y, model_conf = model_conf, 
                                history_burn = history_burn)
        lpd_features$feat = scale(lpd_features$feat, center = TRUE, scale = TRUE)
        
    }else{
        if(is.na(model_conf$lpd_features_parl$ncores)){
            ncores = parallel::detectCores()
        }else{
            ncores = model_conf$lpd_features_parl$ncores
        }
        
        t_seq = c((history_burn):(length(y) - train_h))
        num_block = ceiling(length(t_seq)/ncores)
        for (ncores0 in 1:ncores) {
            if (ncores0*num_block >= length(t_seq) & 
                (ncores0-1)*num_block < length(t_seq)){
               break
            }
        }
        ncores = ncores0
        t_seqs = list()
        for (i in 1:ncores) {
            if(i != ncores){
                t_seqs[[i]] = t_seq[(num_block*(i-1) + 1): (num_block*i)]
            }else{
                t_seqs[[i]] = t_seq[(num_block*(i-1) + 1): length(t_seq)]
            }
        }
        # lpd_features0 = mclapply(t_block, lpd_feat, ts_sd = y1, ts_nosd = y, 
        #                         model_conf = model_conf, mc.cores = ncores)
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
       
        lpd_features0 = foreach::foreach(i = 1:ncores, .packages = c("rugarch","M4metalearning"), 
                                         .export = c("lpd_feat", model_conf$fore_model))%dopar%
            ## "R.utils"
            # withTimeout({
            #     lpd_feat(t_seq = t_seqs[[i]], ts_sd = y1, ts_nosd = y,
            #          model_conf = model_conf, history_burn = history_burn)
            # }, timeout = 30.00, substitute = FALSE, onTimeout = "error",
            # TimeoutException = function(ex) {
            #     options(warn = 1)
            #     warning("Timeout in part: ", i)
            # })
            lpd_feat(t_seq = t_seqs[[i]], ts_sd = y1, ts_nosd = y,
                     model_conf = model_conf, history_burn = history_burn)
        stopCluster(cl)
        # lpd_features0 = lapply(t_seqs, lpd_feat, ts_sd = y1, ts_nosd = y, 
        #                        model_conf = model_conf)
        lpds = lapply(lpd_features0, function(x) x$lpd)
        feats = lapply(lpd_features0, function(x) x$feat)
        lpd_features = list(lpd = do.call(rbind, lpds), feat = do.call(rbind, feats))
        lpd_features$feat = scale(lpd_features$feat, center = TRUE, scale = TRUE)
        rm(lpds)
        rm(feats)
    }
    
    return(lpd_features)
}

                       
lpd_feat = function(t_seq, ts_sd, ts_nosd, model_conf, history_burn ){
    feature_window = model_conf$feature_window
    roll = model_conf$roll
    frequency = model_conf$frequency
    ets_model = model_conf$ets_model
    forecast_h = model_conf$forecast_h
    train_h = model_conf$train_h
    PI_level = model_conf$PI_level
    fore_model = model_conf$fore_model
    
    y1 = ts_sd
    log_pred_densities <- matrix(nrow = length(t_seq), ncol = length(fore_model))
    colnames(log_pred_densities) <- unlist(fore_model)
    for (t in t_seq)
    {
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
        rm(use_model)
        log_pred_den <- as.numeric(log_pred_den)
        log_pred_den[log_pred_den < log(1e-323)] <- log(1e-323)
        log_pred_den[log_pred_den > log(1e+308)] <- log(1e+308)
        log_pred_densities[(t - t_seq[1] + 1), ] <- log_pred_den
        # options(warn = 1)
        # warning("The forecasting models of time ", t, " finished!")
    }
    
    ## Calculate historical features
    y = ts_nosd
    features_y <- matrix(nrow = length(t_seq), ncol = 42)
    myts <- list(list(x = ts(y[1:history_burn], frequency = frequency)))
    colnames(features_y) <- colnames(THA_features(myts)[[1]]$features)
    
    if(is.null(feature_window)){
        for (t in t_seq)
        {
            myts <- list(list(x = ts(y[1:t], frequency = frequency)))
            myfeatures <- THA_features(myts)[[1]]$features
            myfeatures <- data.matrix(myfeatures)
            features_y[(t - t_seq[1] + 1),] <- myfeatures
        }
    }else{
        for (t in t_seq)
        {
            if(t <= feature_window){
                myts <-list(list(x=ts(y[1:t], frequency = 1)))
                myfeatures <- THA_features(myts)[[1]]$features
                myfeatures <- data.matrix(myfeatures)
            }else{
                myts <-list(list(x=ts(y[(t-feature_window+1):t], frequency = 1)))
                myfeatures <- THA_features(myts)[[1]]$features
                myfeatures <- data.matrix(myfeatures)
            }
            features_y[(t - t_seq[1] + 1),] <- myfeatures
        }
    }
    lpd_features <- list(lpd = log_pred_densities, feat = features_y)
    rm(log_pred_densities)
    rm(features_y)
    return(lpd_features)
}


#' Delete the features with NaN and add attributes
#'
#' @param lpd_features A list of several outputs of \code{lpd_features_multi}.
#' 
#' @return A list with log predictive densities and features. 
#' @export
feature_clean <- function(lpd_features){
    for (i_ts in 1:length(lpd_features)) {
        NA_ind <- apply(lpd_features[[i_ts]]$feat, 2, anyNA)
        lpd_features[[i_ts]]$feat_mean <- attr(lpd_features[[i_ts]]$feat, "scaled:center")[!NA_ind]
        lpd_features[[i_ts]]$feat_sd <- attr(lpd_features[[i_ts]]$feat, "scaled:scale")[!NA_ind]
        lpd_features[[i_ts]]$feat <- lpd_features[[i_ts]]$feat[, !NA_ind]
    }
    return(lpd_features)
}

