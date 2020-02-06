inverse_log_logistic = function(u, t0, alpha, beta, lambda = 1){


  alpha * (
    (1 + (t0/alpha) ^ beta) /
      ((1 - u)^ (1 / lambda)) - 1) ^ (1 / beta) - t0

}
