#' @export
emplik_tcf2 <- function(out_predict_tcf2, out_boot_tcf2, n, true_cov_tcf2,
                        true_adj_tcf2){
  ll_tcf2_est <- mapply(function(u, v) {
    ll_emp(n = n, theta_true = u, theta_est = v)
  }, u = true_cov_tcf2, v = out_predict_tcf2$out_theta_2_est$theta_2_est)
  res <- list()
  res$ll_tcf2_est <- ll_tcf2_est
  res$ll_tcf2_adj <- ll_emp(n = n, theta_true = true_adj_tcf2,
                            theta_est = out_predict_tcf2$theta_2_est_adj)
  res_lr_cov_tcf2 <- out_boot_tcf2$w_cov_tcf2 * res$ll_tcf2_est
  res_lr_adj_tcf2 <- out_boot_tcf2$w_adj_tcf2 * res$ll_tcf2_adj
  res$res_lr_cov_tcf2 <- res_lr_cov_tcf2
  res$res_lr_adj_tcf2 <- res_lr_adj_tcf2
  return(res)
}

#' @export
ci_ll_emp_fun <- function(theta, theta_est, n, w_adj, qc) {
  ll_est <- ll_emp(n = n, theta_true = theta, theta_est = theta_est)
  if(is.na(ll_est)) ll_est <- Inf
  ll_est_adj <- ll_est
  if(!is.infinite(ll_est)){
    ll_est_adj <- w_adj * ll_est_adj
  }
  return(ll_est_adj - qc)
}

#' @export
ci_tcf2 <- function(out_pred_tcf2, out_tcf2_boot, ci_level, n) {
  cov_tcf2_est <- out_pred_tcf2$out_theta_2_est$theta_2_est
  qc <- qchisq(ci_level, 1)
  m <- length(cov_tcf2_est)
  out_ci <- sapply(1:m, function(i){
    lwr <- uniroot(f = ci_ll_emp_fun, interval = c(0, cov_tcf2_est[i]),
                   theta_est = cov_tcf2_est[i], n = n, qc = qc,
                   w_adj = out_tcf2_boot$w_cov_tcf2[i], tol = 1e-8)$root
    upr <- uniroot(f = ci_ll_emp_fun, interval = c(cov_tcf2_est[i], 1),
                   theta_est = cov_tcf2_est[i], n = n, qc = qc,
                   w_adj = out_tcf2_boot$w_cov_tcf2[i], tol = 1e-8)$root
    return(c(lwr, upr))
  })
  res <- out_pred_tcf2
  res$out_theta_2_est$lwr <- out_ci[1,]
  res$out_theta_2_est$upr <- out_ci[2,]
  ###
  adj_tcf2_est <- out_pred_tcf2$theta_2_est_adj
  lwr_adj_tcf2 <- uniroot(f = ci_ll_emp_fun,
                          interval = c(0, adj_tcf2_est),
                          theta_est = adj_tcf2_est, n = n, qc = qc,
                          w_adj = out_tcf2_boot$w_adj_tcf2, tol = 1e-8)$root
  upr_adj_tcf2 <- uniroot(f = ci_ll_emp_fun,
                          interval = c(adj_tcf2_est, 1),
                          theta_est = adj_tcf2_est, n = n, qc = qc,
                          w_adj = out_tcf2_boot$w_adj_tcf2, tol = 1e-8)$root
  res$theta_2_est_adj <- c("estimate" = res$theta_2_est_adj,
                           "lwr" = lwr_adj_tcf2, "upr" = upr_adj_tcf2)
  return(res)
}
