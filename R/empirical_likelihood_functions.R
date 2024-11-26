#' @export
ll_emp <- function(n, theta_true, theta_est){
  temp1 <- (theta_est * log(theta_est/theta_true) +
              ((1 - theta_est) * log((1 - theta_est)/(1 - theta_true))))
  temp1 <- 2 * n * temp1
  if (is.nan(temp1)) temp1 <- Inf
  return(temp1)
}

#' @export
ll_emp_vec <- Vectorize(FUN = ll_emp, vectorize.args = "theta_est")

