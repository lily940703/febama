log_score <- function(beta, features, features_select = NULL, prob, intercept = T){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)
  num_model <- dim(prob)[2]
  # 避免Inf/Inf
  me <- features0 %*% beta
  me[me>709] <- 709
  exp_lin = exp(me)
  deno = matrix (rep((1+rowSums(exp_lin)), num_model-1), ncol = num_model-1)
  w <- exp_lin/ deno # T-by-(n-1)
  w_full = cbind(w, 1 - rowSums(w)) # T-by-n
  out = sum(log(rowSums(w_full * prob)))
  return(out)
}

gradient_logscore <- function(beta, features, features_select = NULL, prob, intercept){
  if(is.null(features_select)){
    features0<-features
  }else{
    features0<-data.matrix(features[,features_select])
  }
  if(intercept) features0 = cbind(rep(1, nrow(prob)), features0)

  ex = exp(features0 %*% beta)
  ex[ex > exp(700)] = exp(700)

  ex_sum = rowSums(ex)
  n = dim(prob)[2]
  gradient0<-function(i, t, p){
    m = p[t,i] *(1 + ex_sum[t] - ex[t,i])
    if(m > 1e+308) m = 1e+308
    s = sum(p[t,-c(i,n)] * ex[t,-i]) - p[t,n]
    if(s > 1e+308) s = 1e+308
    q = m - s
    if(q > 1e+308) q = 1e+308
    if(q < -1e+308) q = -1e+308
    a = ex[t,i] * (q) * features0[t,]
    a[a>1e+308] = 1e+308;a[a < (-1e+308)] = -1e+308
    b = (sum(ex[t,] * p[t, 1:(n-1)]) + p[t,n] ) * (1 + ex_sum[t])
    if(b > 1e+308) b = 1e+308
    if(b < 1e-323) b = 1e-323
    out0 = a/b
    return(out0)
  }

  out_t<-c()
  out_i<-c()
  for (i in 1:(n-1)) {
    for (t in 1:length(prob[,1])) {
      out_t<-cbind(out_t, gradient0(i, t, prob))
    }
    out_i<-cbind(out_i, rowSums(out_t))
  }
  return(out_i)
}
