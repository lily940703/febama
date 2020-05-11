# p(beta | I) * p(I)
# beta = [beta1, beta2,..., beta(n-1)]
# beta1,...,beta(n-1), iid,  ~ N(0,sigma), sigma = diag{10,...,10}
prior <- function(beta, I, sig = 10){
  pri <- apply(beta, 2, dmvnorm, mean=matrix(0, dim(beta)[1], 1), sigma = sig*diag(dim(beta)[1]))
  pri <- exp(sum(log(pri)))
  # avoid prior=0
  if (pri < 1e-323) {pri <- 1e-323}
  return(pri)
}

# the gradient of the prior can be computed in this way only when
# beta1,...,beta(n-1), iid,  ~ N(0,sigma)
gradient_prior <- function(beta, I, sig){
  coef = -prior(beta, I, sig)
  sigma_inver = solve(sig*diag(dim(beta)[1]))
  mu = matrix(0, dim(beta)[1], 1)
  beta_gra = apply(beta, 2, function(beta){
    sigma_inver %*% (beta-mu)
  })
  if(is.vector(beta_gra)) beta_gra <- t(as.matrix(beta_gra))
  return(coef * beta_gra)
}

# beta <- matrix(rep(5,20), ncol = 4)
# I <- c(0,1,0,1,0,1,1)
# prior (beta, I, sig = 10)
# gradient_prior (beta, I, sig = 10)
# a1 <- jacobian(func = dmvnorm , x = as.vector(beta[,1]), method="Richardson",
#          mean=matrix(0, 5, 1),sigma = 10*diag(5))
# a1 * prior (beta[,2:4], I, sig = 10)


log_posterior <- function(data, beta, I, prior = prior,
                          logLik = log_score, sig = 10){
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
