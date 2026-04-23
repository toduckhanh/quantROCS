#' @export
fit_tcf2 <- function(formula_q1, formula_q3, formula_tcf2, data_1, data_2,
                     data_3, theta_10, theta_30, family, ...) {
  res <- list()
  res$formula_tcf2_org <- formula_tcf2
  res_rq_1 <- rq(formula_q1, tau = theta_10, data = data_1, ...)
  res_rq_3 <- rq(formula_q3, tau = 1 - theta_30, data = data_3, ...)
  t1_est <- as.numeric(
    model.matrix(formula_q1, data = data_2) %*% res_rq_1$coefficients
  )
  t2_est <- as.numeric(
    model.matrix(formula_q3, data = data_2) %*% res_rq_3$coefficients
  )
  frame_data_2 <- model.frame(formula_tcf2, data = data_2)
  response <- model.response(frame_data_2)
  working_response <- as.numeric((response > t1_est) * (response <= t2_est))
  formula_tcf2 <- update.formula(formula_tcf2, new = working_response ~ .)
  data_2$working_response <- working_response
  fit_logis <- glm(formula_tcf2, family = family, data = data_2, ...)
  res$t1_est <- t1_est
  res$t2_est <- t2_est
  res$res_rq_1 <- res_rq_1
  res$res_rq_3 <- res_rq_3
  res$fit_logis <- fit_logis
  res$working_response <- working_response
  res$n <- nrow(data_2)
  return(res)
}

#' @export
predict_tcf2 <- function(out_fit_tcf2, newdata) {
  theta_2_est <- predict(out_fit_tcf2$fit_logis, newdata = newdata,
                         type = "response")
  theta_2_est_adj <- mean(out_fit_tcf2$working_response)
  res <- list()
  res$out_theta_2_est <- newdata
  res$out_theta_2_est$theta_2_est <- theta_2_est
  res$theta_2_est_adj <- theta_2_est_adj
  return(res)
}

#' @export
fit_rq13_boot <- function(out_fit_tcf2, theta_10, theta_30, B, bsmethod, ...) {
  res_rq_1_boot <- try(boot.rq(x = out_fit_tcf2$res_rq_1$x,
                               y = out_fit_tcf2$res_rq_1$y,
                               tau = theta_10, R = B,
                               bsmethod = bsmethod, ...),
                       silent = TRUE)
  res_rq_3_boot <- try(boot.rq(x = out_fit_tcf2$res_rq_3$x,
                               y = out_fit_tcf2$res_rq_3$y,
                               tau = 1 - theta_30,
                               R = B, bsmethod = bsmethod, ...),
                       silent = TRUE)
  if(inherits(res_rq_1_boot, "try-error") ||
     inherits(res_rq_3_boot, "try-error")) res <- NULL
  else {
    res <- list()
    res$res_rq_1_boot <- res_rq_1_boot$B
    res$res_rq_3_boot <- res_rq_3_boot$B
  }
  return(res)
}

#' @export
fit_tcf2_boot <- function(out_fit_tcf2, out_predict_tcf2, theta_10, theta_30, B,
                          method = c("mb", "np"),
                          type = c("mean", "median", "variance"),
                          newdata, bsmethod, ...) {
  # mb = model-based; np = non-parametric
  method <- match.arg(method)
  data_2_mb <- out_fit_tcf2$fit_logis$data
  n <- out_fit_tcf2$n
  out_boot <- list()
  out_boot$method <- method
  if (method == "mb") {
    prob_mb <- predict(out_fit_tcf2$fit_logis, type = "response")
    working_response_mb <- replicate(B, rbinom(n, size = 1, prob = prob_mb))
    theta_2_est_adj_mb <- colMeans(working_response_mb)
    theta_2_est_adj_mb[theta_2_est_adj_mb == 0] <- 1 - n/(n + 0.5)
    theta_2_est_adj_mb[theta_2_est_adj_mb == 1] <- n/(n + 0.5)
    formula_tcf2_mb <- update.formula(out_fit_tcf2$fit_logis$formula,
                                      new = working_response_mb ~ .)
    theta_2_est_mb <- apply(working_response_mb, 2, function(x) {
      data_2_mb$working_response_mb <- x
      fit_logis_mb <- glm(formula_tcf2_mb,
                          family = out_fit_tcf2$fit_logis$family,
                          data = data_2_mb)
      res <- predict(fit_logis_mb, newdata = newdata, type = "response")
      res[res == 0] <- 1 - n/(n + 0.5)
      res[res == 1] <- n/(n + 0.5)
      return(res)
    })
    out_boot$theta_2_est_boot <- theta_2_est_mb
    out_boot$theta_2_est_adj_boot <- theta_2_est_adj_mb
  } else {
    out_rq13_bts <- fit_rq13_boot(out_fit_tcf2 = out_fit_tcf2,
                                  theta_10 = theta_10, theta_30 = theta_30,
                                  B = B, bsmethod = bsmethod)
    theta_2_est_np <- matrix(0, nrow = nrow(newdata), ncol = B)
    theta_2_est_adj_np <- numeric(B)
    for (i in 1:B) {
      id_2 <- sample(1:n, size = n, replace = TRUE)
      data_2_bst <- data_2_mb[id_2, ]
      t1_est_bts <- as.numeric(
        model.matrix(out_fit_tcf2$res_rq_1$formula, data = data_2_bst) %*%
          out_rq13_bts$res_rq_1_boot[i, ]
      )
      t2_est_bts <- as.numeric(
        model.matrix(out_fit_tcf2$res_rq_3$formula, data = data_2_bst) %*%
          out_rq13_bts$res_rq_3_boot[i, ]
      )
      frame_data_2_bts <- model.frame(out_fit_tcf2$formula_tcf2_org,
                                      data = data_2_bst)
      response_bts <- model.response(frame_data_2_bts)
      working_response_bst <- as.numeric(
        (response_bts > t1_est_bts) * (response_bts <= t2_est_bts)
      )
      data_2_bst$working_response <- working_response_bst
      fit_logis_bst <- glm(out_fit_tcf2$fit_logis$formula,
                           family = out_fit_tcf2$fit_logis$family,
                           data = data_2_bst)
      temp_pred <- predict(fit_logis_bst, newdata = newdata,
                           type = "response")
      temp_pred[temp_pred == 0] <- 1 - n/(n + 0.5)
      temp_pred[temp_pred == 1] <- n/(n + 0.5)
      theta_2_est_np[, i] <- temp_pred
      temp_adj <- mean(data_2_bst$working_response)
      if(temp_adj == 0) temp_adj <- 1 - n/(n + 0.5)
      if(temp_adj == 1) temp_adj <- n/(n + 0.5)
      theta_2_est_adj_np[i] <- temp_adj
    }
    out_boot$theta_2_est_boot <- theta_2_est_np
    out_boot$theta_2_est_adj_boot <- theta_2_est_adj_np
  }
  out_boot$sd_theta_2_est <- apply(out_boot$theta_2_est_boot, 1, sd)
  out_boot$sd_theta_2_est_adj <- sd(out_boot$theta_2_est_adj_boot)
  theta_2_est_ls <- split(out_boot$theta_2_est_boot,
                          row(out_boot$theta_2_est_boot))
  theta_2_est <- out_predict_tcf2$out_theta_2_est$theta_2_est
  type <- match.arg(type)
  w_cov_tcf2 <- switch (type,
    mean = 1/colMeans(mapply(function(u, v) {
      ll_emp_vec(n = n, theta_true = u, theta_est = v)
    }, u = as.list(theta_2_est), v = theta_2_est_ls)),
    median = qchisq(0.5, 1)/apply(mapply(function(u, v) {
      ll_emp_vec(n = n, theta_true = u, theta_est = v)
    }, u = as.list(theta_2_est), v = theta_2_est_ls), 2, median),
    variance = theta_2_est*(1 - theta_2_est)/(n*out_boot$sd_theta_2_est^2)
  )
  out_boot$w_cov_tcf2 <- w_cov_tcf2
  theta_2_adj_est <- out_predict_tcf2$theta_2_est_adj
  w_adj_tcf2 <- switch (type,
    mean = 1/mean(sapply(out_boot$theta_2_est_adj_boot, function (x){
      ll_emp(n = n, theta_true = theta_2_adj_est,
             theta_est = x)
    })),
    median = qchisq(0.5, 1)/median(sapply(out_boot$theta_2_est_adj_boot, function (x){
      ll_emp(n = n, theta_true = theta_2_adj_est,
             theta_est = x)
    })),
    variance = theta_2_adj_est*(1 - theta_2_adj_est)/(n*out_boot$sd_theta_2_est_adj^2)
  )
  out_boot$w_adj_tcf2 <- w_adj_tcf2
  return(out_boot)
}
