 # beta <- matrix(rep(5,20), ncol = 4)
 # I <- c(0,1,0,1,0,1,1)
 # prior (beta, I, sig = 10)
 # gradient_prior (beta, I, sig = 10)
 # a1 <- jacobian(func = dmvnorm , x = as.vector(beta[,1]), method="Richardson",
 #          mean=matrix(0, 5, 1),sigma = 10*diag(5))
 # a1 * prior (beta[,2:4], I, sig = 10)

log_posterior <- function(data, beta, I, prior = prior,
                          logLik = logscore, sig = 10){
    pri <- prior(beta, I, sig = sig)
    features_select <- which(I==1)
    prob <- exp(data$lpd)
    prob[prob == 0] <- 1e-323
    features <- data$feat
    LS <- logLik(beta = beta, features = features, features_select = features_select,
                 prob = prob, intercept = TRUE)
    log_post <- log(pri) + LS
    return(log_post)
}
