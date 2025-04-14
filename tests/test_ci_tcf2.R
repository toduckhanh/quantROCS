# 1. load required package ----
require(quantROCS)

# 2. setup parameters ----
theta_10 <- 0.3
theta_30 <- 0.3
x_eval <- seq(0.05, 0.95, length.out = 25)

# 3. simulate data ----
n <- 150
z1 <- runif(n, 0, 1)
rho_1n <- runif(n, 0, 1)
y1 <- 1 + qnorm(rho_1n) + (- qnorm(rho_1n) - 0.5) * z1
data_1 <- data.frame(y1 = y1, z = z1)

z3 <- runif(n, 0, 1)
rho_3n <- runif(n, 0, 1)
y3 <- 2 + qnorm(1 - rho_3n) + (3 - qnorm(1 - rho_3n)) * z3
data_3 <- data.frame(y3 = y3, z = z3)

z2 <- runif(n, 0, 1)
y2 <- 1.5 + 0.5*z2 + rlogis(n, 0, 1)
data_2 <- data.frame(y2 = y2, z = z2)

# 4. fitting models ----
out_test <- fit_tcf2(formula_q1 = y1 ~ z, formula_q3 = y3 ~ z,
                     formula_tcf2 = y2 ~ z, data_1 = data_1, data_2 = data_2,
                     data_3 = data_3, theta_10 = 0.3, theta_30 = 0.3,
                     family = binomial(link = "logit"))
out_test

out_test_predict <- predict_tcf2(out_fit_tcf2 = out_test,
                                 newdata = data.frame(z = x_eval))
out_test_predict

out_test_boot_mcmb <- fit_tcf2_boot(
  out_fit_tcf2 = out_test, out_predict_tcf2 = out_test_predict,
  theta_10 = 0.3, theta_30 = 0.3, B = 200, method = "np",
  newdata = data.frame(z = x_eval), bsmethod = "mcmb"
)
out_test_boot_mcmb

theta_seq <- seq(0, 1, by = 0.05)

ll_6 <- sapply(theta_seq, function(x){
  ci_ll_emp_fun(theta = x,
                theta_est = out_test_predict$out_theta_2_est$theta_2_est[6],
                n = n, w_adj = out_test_boot_mcmb$w_cov_tcf2[6],
                qc = qchisq(0.95, 1))
})

my_fun1 <- function(theta, theta_est, n, w_adj){
  ll_est <- ll_emp(n = n, theta_true = theta, theta_est = theta_est)
  if(is.na(ll_est)) ll_est <- Inf
  ll_est_adj <- ll_est
  if(!is.infinite(ll_est)){
    ll_est_adj <- w_adj * ll_est_adj
  }
  return(list("-2LLR" = ll_est_adj))
}

emplik::findUL(step = 0.01, fun = my_fun1,
               MLE = out_ci_tcf2$out_theta_2_est$theta_2_est[3],
               theta_est = out_test_predict$out_theta_2_est$theta_2_est[3],
               n = n, w_adj = out_test_boot_mcmb$w_cov_tcf2[3],
               level = qchisq(0.95, 1))

plot(theta_seq, exp((ll_6 + qchisq(0.95, 1))), type = "l")
abline(h = 0, col = "blue", lty = 2)

out_ci_tcf2 <- ci_tcf2(out_pred_tcf2 = out_test_predict,
                       out_tcf2_boot = out_test_boot_mcmb,
                       ci_level = 0.95, n = n)

plot(x_eval, out_ci_tcf2$out_theta_2_est$theta_2_est, type = "l",
     xlim = c(0, 1), ylim = c(0, 1), ylab = "Cov TCF2")
lines(x_eval, out_ci_tcf2$out_theta_2_est$lwr, col = "blue")
lines(x_eval, out_ci_tcf2$out_theta_2_est$upr, col = "blue")


out_ci_tcf2_2 <- ci_tcf2_optim(out_pred_tcf2 = out_test_predict,
                               out_tcf2_boot = out_test_boot_mcmb,
                               ci_level = 0.95, n = n)
out_ci_tcf2_2

plot(x_eval, out_ci_tcf2_2$out_theta_2_est$theta_2_est, type = "l",
     xlim = c(0, 1), ylim = c(0, 1), ylab = "Cov TCF2")
lines(x_eval, out_ci_tcf2_2$out_theta_2_est$lwr, col = "blue")
lines(x_eval, out_ci_tcf2_2$out_theta_2_est$upr, col = "blue")

