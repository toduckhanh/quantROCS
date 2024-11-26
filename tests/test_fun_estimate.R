require(quantROCS)

theta_10 <- 0.3
theta_30 <- 0.3
x_eval <- seq(0.05, 0.95, length.out = 25)

fun_simu_model_tcf2 <- function(n, theta_10, theta_30) {
  z2 <- runif(n, 0, 1)
  y2 <- 1.5 + 0.5*z2 + rlogis(n, 0, 1)
  t1_0 <- 1 + qnorm(theta_10) + (- qnorm(theta_10) - 0.5) * z2
  t2_0 <- 2 + qnorm(1 - theta_30) + (3 - qnorm(1 - theta_30)) * z2
  y2_wk <- as.numeric((y2 >= t1_0) * (y2 <= t2_0))
  res_logis <- glm(y2_wk ~ z2, family = binomial(link = "logit"))
  theta_2_est <- predict(
    res_logis,
    newdata = data.frame(z2 = seq(0.05, 0.95, length.out = 25)),
    type = "response"
  )
  theta_2_est_adj <- mean(y2_wk)
  return(list(res_logis$coefficients, theta_2_est, theta_2_est_adj))
}

system.time({
  out_md_logis_03_100 <- sapply(1:5000, FUN = function(i) {
    fun_simu_model_tcf2(5000, theta_10 = 0.3, theta_30 = 0.3)
  })
})

true_coef_03 <- apply(
  apply(out_md_logis_03_100, 2, function(x) x[[1]]), 1, mean
)

true_cov_tcf2_03 <- apply(
  apply(out_md_logis_03_100, 2, function(x) x[[2]]), 1, mean
)

true_adj_tcf2_03 <- mean(apply(out_md_logis_03_100, 2, function(x) x[[3]]))

n <- 50
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

out_test <- fit_tcf2(formula_q1 = y1 ~ z, formula_q3 = y3 ~ z,
                     formula_tcf2 = y2 ~ z, data_1 = data_1, data_2 = data_2,
                     data_3 = data_3, theta_10 = 0.3, theta_30 = 0.3,
                     family = binomial(link = "logit"))

frame_data_2_bts <- model.frame(out_test$fit_logis$formula,
                                data = out_test$fit_logis$data)
response_bts <- model.response(frame_data_2_bts)

out_test_predict <- predict_tcf2(out_fit_tcf2 = out_test,
                                 newdata = data.frame(z = x_eval))
out_test_predict


out_test_boot_xy <- fit_tcf2_boot(
  out_fit_tcf2 = out_test, out_predict_tcf2 = out_test_predict,
  theta_10 = 0.3, theta_30 = 0.3, B = 200, method = "np",
  newdata = data.frame(z = x_eval), bsmethod = "xy"
)
out_test_boot_xy

out_test_boot_mcmb <- fit_tcf2_boot(
  out_fit_tcf2 = out_test, out_predict_tcf2 = out_test_predict,
  theta_10 = 0.3, theta_30 = 0.3, B = 200, method = "np",
  newdata = data.frame(z = x_eval), bsmethod = "mcmb"
)
out_test_boot_mcmb


lr_test_tcf2_mcmb <- emplik_tcf2(out_predict_tcf2 = out_test_predict,
                                 out_boot_tcf2 = out_test_boot_mcmb, n = 50,
                                 true_cov_tcf2 = true_cov_tcf2_03,
                                 true_adj_tcf2 = true_adj_tcf2_03)
lr_test_tcf2_mcmb$res_lr_cov_tcf2 <= qchisq(0.95, 1)

lr_test_tcf2_xy <- emplik_tcf2(out_predict_tcf2 = out_test_predict,
                               out_boot_tcf2 = out_test_boot_xy, n = 50,
                               true_cov_tcf2 = true_cov_tcf2_03,
                               true_adj_tcf2 = true_adj_tcf2_03)
lr_test_tcf2_xy$res_lr_cov_tcf2 <= qchisq(0.95, 1)




theta_2_est_bp_ls <- split(out_test_boot_mcmb$theta_2_est_boot,
                           row(out_test_boot_mcmb$theta_2_est_boot))
w_cov <- 1/colMeans(mapply(function(u, v) {
  ll_emp_vec(n = n, theta_true = u, theta_est = v)
}, u = as.list(out_test$theta_2_est), v = theta_2_est_bp_ls))


B <- 200
data_2_bst <- sapply(1:B, function(i) {
  id_2 <- sample(1:n, size = n, replace = TRUE)
  return(data_2[id_2, ])
}, simplify = FALSE)

out_rq13_bts_xy <- fit_rq13_boot(out_fit_tcf2 = out_test, theta_10 = 0.3,
                                 theta_30 = 0.3, B = B, bsmethod = "xy")
out_rq13_bts_xy$res_rq_1_boot

t1_est_bts <- as.numeric(
  model.matrix(out_test$res_rq_1$formula, data = data_2_bst[[1]]) %*%
    out_rq13_bts_xy$res_rq_1_boot[1,]
)

data_2_bst <- data_2_mb[id_2, ]

t2_est_bts <- as.numeric(
  model.matrix(out_fit_tcf2$res_rq_3$formula, data = data_2_bst) %*%
    out_rq13_bts$res_rq_3_boot[i,]
)
frame_data_2_bts <- model.frame(out_fit_tcf2$fit_logis$formula,
                                data = data_2_bst)
response_bts <- model.response(frame_data_2_bts)
working_response_bst <- as.numeric(
  (response_bts > t1_est_bts) * (response_bts <= t2_est_bts)
)



sd_theta_2_est_bp <- apply(theta_2_est_bp, 1, sd)
sd_theta_2_est_adj_bp <- sd(theta_2_est_adj_bp)
sd_theta_2_est_mcmb <- apply(theta_2_est_mcmb, 1, sd)
sd_theta_2_est_adj_mcmb <- sd(theta_2_est_adj_mcmb)
## empirical-likelihood for covariate est with bootstrap bp
theta_2_est_bp_ls <- split(theta_2_est_bp, row(theta_2_est_bp))
w_cov_bp <- 1/colMeans(mapply(function(u, v) {
  ll_emp_vec(n = n, theta_true = u, theta_est = v)
}, u = as.list(res_est[[4]]), v = theta_2_est_bp_ls))
## empirical-likelihood for covariate est with bootstrap mcmb
theta_2_est_mcmb_ls <- split(theta_2_est_mcmb, row(theta_2_est_mcmb))
w_cov_mcmb <- 1/colMeans(mapply(function(u, v) {
  ll_emp_vec(n = n, theta_true = u, theta_est = v)
}, u = as.list(res_est[[4]]), v = theta_2_est_mcmb_ls))
## empirical-likelihood for adjusted est with bootstrap bp
ll_theta2_est_adj_bp <- sapply(theta_2_est_adj_bp, function (x){
  ll_emp(n = n, theta_true = res_est[[5]], theta_est = x)
})
w_bp <- 1/mean(ll_theta2_est_adj_bp)
## empirical-likelihood for adjusted est with bootstrap mcmb
ll_theta2_est_adj_mcmb <- sapply(theta_2_est_adj_mcmb, function (x){
  ll_emp(n = n, theta_true = res_est[[5]], theta_est = x)
})
w_mcmb <- 1/mean(ll_theta2_est_adj_mcmb)


hist(out_test_boot_xy$theta_2_est_adj_boot)
hist(out_test_boot_mcmb$theta_2_est_adj_boot)

hist(out_test_boot_xy$theta_2_est_boot[3,])
hist(out_test_boot_mcmb$theta_2_est_boot[3,])

u1 <- apply(out_test_boot_xy$theta_2_est_boot, 1,
            function(x) quantile(x, probs = c(0.025, 0.975)))
u2 <- apply(out_test_boot_mcmb$theta_2_est_boot, 1,
            function(x) quantile(x, probs = c(0.025, 0.975)))

plot(x_eval, out_test$theta_2_est, type = "l", ylim = c(0.2, 1))
lines(x_eval, u1[1,], col = "blue", lty = 2)
lines(x_eval, u1[2,], col = "blue", lty = 2)
lines(x_eval, u2[1,], col = "red", lty = 2)
lines(x_eval, u2[2,], col = "red", lty = 2)

###
predict(out_test$fit_logis, type = "response")

out_rq13_bts_xy <- fit_rq13_boot(out_fit_tcf2 = out_test, theta_10 = 0.3,
                                 theta_30 = 0.3, B = 200, bsmethod = "xy")
out_rq13_bts_xy$res_rq_1_boot

# first bootstrap
theta_2_est_bst <- matrix(0, nrow = length(x_eval), ncol = 200)
theta_2_est_adj_bst <- numeric(200)

for (i in 1:200) {
  id_2 <- sample(1:n, size = n, replace = TRUE)
  data_bst_2 <- out_test$fit_logis$data[id_2, ]
  t1_est_bts <- as.numeric(
    model.matrix(out_test$res_rq_1$formula, data = data_bst_2) %*%
      out_rq13_bts_xy$res_rq_1_boot[i,]
  )
  t2_est_bts <- as.numeric(
    model.matrix(out_test$res_rq_3$formula, data = data_bst_2) %*%
      out_rq13_bts_xy$res_rq_3_boot[i,]
  )
  frame_data_2_bts <- model.frame(out_test$fit_logis$formula, data = data_bst_2)
  response_bts <- model.response(frame_data_2_bts)
  working_response_bst <- as.numeric(
    (response_bts > t1_est_bts) * (response_bts <= t2_est_bts)
  )
  data_bst_2$working_response <- working_response_bst
  fit_logis_bst <- glm(out_test$fit_logis$formula,
                       family = binomial(link = "logit"),
                       data = data_bst_2)
  theta_2_est_bst[,i] <- predict(fit_logis_bst,
                                 newdata = data.frame(z = x_eval),
                                 type = "response")
  theta_2_est_adj_bst[i] <- mean(data_bst_2$working_response)
}


hist(theta_2_est_adj_bst)
hist(theta_2_est_bst[3,])

u1 <- apply(theta_2_est_bst, 1,
            function(x) quantile(x, probs = c(0.025, 0.975)))

plot(x_eval, out_test$theta_2_est, type = "l", ylim = c(0.2, 1))
lines(x_eval, u1[1,], col = "blue", lty = 2)
lines(x_eval, u1[2,], col = "blue", lty = 2)

#####

formula_q1 <- y1 ~ z
formula_q3 <- y3 ~ z
formula_tcf2 <- y2 ~ z

model.frame(formula_q1, data = data_2)
model.matrix(formula_q1, data = data_2)

t1_est <- predict(res_rq_1, newdata = data.frame(z = data_2$z))

cbind(t1_est, model.matrix(formula_q1, data = data_2) %*% res_rq_1$coefficients)

fit_rq13_boot(out_fit_tcf2 = out_test, theta_10 = 0.3, theta_30 = 0.3, B = 200,
              bsmethod = "mcmb")

fit_rq13_boot(out_fit_tcf2 = out_test, theta_10 = 0.3, theta_30 = 0.3, B = 200,
              bsmethod = "wxy")

out_test_boot <- fit_tcf2_boot(out_test, theta_10 = 0.3, theta_30 = 0.3,
                               B = 200, method = "mb",
                               newdata = data.frame(z = x_eval))
rowMeans(out_test_boot[[1]])
mean(out_test_boot[[2]])


plot(x_eval, out_test$theta_2_est, ylim = c(0, 1), type = "l")

out_test_tcf2_est <- sapply(1:1000, FUN = function(i){
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
  out_test <- fit_tcf2(formula_q1 = y1 ~ z, formula_q3 = y3 ~ z,
                       formula_tcf2 = y2 ~ z, data_1 = data_1, data_2 = data_2,
                       data_3 = data_3, theta_10 = 0.3, theta_30 = 0.3,
                       newdata = data.frame(z = x_eval))
  return(out_test$theta_2_est)
})

out_test_tcf2_est

plot(x_eval, rowMeans(out_test_tcf2_est), ylim = c(0, 1), type = "l")

plot(x_eval, out_test_tcf2_est[,1], ylim = c(0, 1), type = "l", col = "gray60")
for (i in 2:1000) {
  lines(x_eval, out_test_tcf2_est[,i], col = "gray60")
}
lines(x_eval, rowMeans(out_test_tcf2_est), col = "blue")

#### NOT RUN #############
formula_q1 <- y1 ~ z
formula_q3 <- y3 ~ z
formula_tcf2 <- y2 ~ z

model.frame(formula_q1, data = data_2)
model.matrix(formula_q1, data = data_2)

res_rq_1 <- rq(formula_q1, tau = theta_10, data = data_1)
boot.rq(x = res_rq_1$x, y = res_rq_1$y, tau = theta_10, R = 200,
        bsmethod = "xy")

t1_est <- predict(res_rq_1, newdata = data.frame(z = data_2$z))

cbind(t1_est, model.matrix(formula_q1, data = data_2) %*% res_rq_1$coefficients)

res_rq_3 <- rq(formula_q3, tau = 1 - theta_30, data = data_3)
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
fit_logis <- glm(formula_tcf2, family = binomial(link = "logit"),
                 data = data_2)

x_eval <- seq(0.05, 0.95, length.out = 25)
newdata <- data.frame(z = x_eval)

theta_2_est <- predict(
  fit_logis,
  newdata = newdata,
  type = "response"
)

theta_2_est
