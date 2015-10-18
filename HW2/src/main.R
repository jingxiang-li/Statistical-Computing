rm(list = ls())

require(MASS)

set.seed(5701)
# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent = FALSE) {
    vars <- codetools::findGlobals(f)
    found <-
        !vapply(vars, exists, logical(1), envir = as.environment(2))
    if (!silent && any(found)) {
        warning("global variables used: ", paste(names(found)[found], collapse =
                                                     ', '))
        return(invisible(FALSE))
    }
    
    ! any(found)
}

# 1 -----------------------------------------------------------------------


#' compute the confidence interval of sigma^2 given i.i.d. observations
#'
#' @param x vector of i.i.d. observations
#' @param alph this is 1 - confidence level
#'
#' @return a vector (left, right), (1 - alpha) confidence interval of sigma^2
#' @export
confint_sigma_sq = function(x, alpha = 0.05) {
    n = length(x)
    var_sample = var(x)
    left = (n - 1) * var_sample / qchisq(1 - alpha / 2, df = n - 1)
    right = (n - 1) * var_sample / qchisq(alpha / 2, df = n - 1)
    c(left, right)
}

#' simulate the confidence interval of sigma^2 for i.i.d. normal random
#' observations
#'
#' @param reps number of replications for the simulation study
#' @param n sample size of normal realization
#' @param mu true mean of the normal distribution
#' @param sigma true standard deviation of the normal distribution
#' @param alpha alpha level used for constructing confidence interval of sigma^2
#' @param score_conf confidence level for constructing the confidence interval
#'   of the coverage probability
#'
#' @return a vector (left, right), score_conf level confidence interval for the
#'   coverage probability
#' @export
sim_confint_sigma_sq_norm = function(reps, n, mu, sigma, alpha, score_conf) {
    x_sim = rnorm(n = n * reps, mean = mu, sd = sigma)
    x_sim = matrix(x_sim, nrow = reps, ncol = n)
    int_sim = apply(x_sim, 1, confint_sigma_sq, alpha = alpha)
    indicator_sim = apply(int_sim, 2, function(x) {
        x[1] < sigma ^ 2 && x[2] > sigma ^ 2
    })
    test_result = prop.test(
        sum(indicator_sim), n = reps,
        p = 1 - alpha, conf.level = score_conf
    )
    test_result$conf.int
}

#' simulate the confidence interval of sigma^2 for i.i.d. Exp(1) random
#' observations
#'
#' @param number of replications for the simulation study
#' @param n sample size of normal realization
#' @param alpha alpha level used for constructing confidence interval of sigma^2
#' @param score_conf confidence level for constructing the confidence interval
#'   of the coverage probability
#'
#' @return a vector (left, right), score_conf level confidence interval for the
#'   coverage probability
#' @export
sim_confint_sigma_sq_exp = function(reps, n, alpha, score_conf) {
    x_sim = rexp(n = n * reps)
    x_sim = matrix(x_sim, nrow = reps, ncol = n)
    int_sim = apply(x_sim, 1, confint_sigma_sq, alpha = alpha)
    indicator_sim = apply(int_sim, 2, function(x) {
        x[1] < 1 && x[2] > 1
    })
    test_result = prop.test(
        sum(indicator_sim), n = reps,
        p = 1 - alpha, conf.level = score_conf
    )
    test_result$conf.int
}

#' simulate the confidence interval of sigma^2 for multivariate random
#' observations, with mean = rep(mu, n), covariance matrix var_mat, where
#' var_mat[i, j] = sigma ^ 2 * rho ^ abs(i - j)
#'
#' @param reps number of replications for the simulation study
#' @param n sample size of normal realization
#' @param mu true mean of the normal distribution
#' @param sigma true marginal standard deviation of the normal distribution
#' @param rho decay parameter used for constructing the covariace matrix
#' @param alpha alpha level used for constructing confidence interval of sigma^2
#' @param score_conf confidence level for constructing the confidence interval
#'   of the coverage probability
#'
#' @return a vector (left, right), score_conf level confidence interval for the
#'   coverage probability
#' @export
sim_confint_sigma_sq_multinorm = function(reps, n, mu, sigma, rho, alpha, score_conf) {
    mean_vec = rep(mu, n)
    var_mat = matrix(0, n, n)
    for (i in 1:n)
        for (j in 1:n) {
            var_mat[i, j] = sigma ^ 2 * rho ^ abs(i - j)
        }
    x_sim = mvrnorm(reps, mean_vec, var_mat)
    int_sim = apply(x_sim, 1, confint_sigma_sq, alpha = alpha)
    indicator_sim = apply(int_sim, 2, function(x) {
        x[1] < sigma ^ 2 && x[2] > sigma ^ 2
    })
    test_result = prop.test(
        sum(indicator_sim), n = reps,
        p = 1 - alpha, conf.level = score_conf
    )
    test_result$conf.int
}

## 1.a
reps = 10000
n = 10
mu = 68
sigma = 3
alpha = 0.05
score_conf = 0.99
sim_confint_sigma_sq_norm(reps, n, mu, sigma, alpha, score_conf)

## 1.b
reps = 10000
n_vec = c(10, 50, 500)
alpha_vec = c(0.01, 0.05)
score_conf = 0.99
result = matrix(0, length(n_vec) * length(alpha_vec), 4)
colnames(result) = c("n", "alpha", "left", "right")
k = 1
for (n in n_vec)
    for (alpha in alpha_vec) {
        result[k, 1] = n
        result[k, 2] = alpha
        result[k, c(3, 4)] = sim_confint_sigma_sq_exp(reps, n, alpha, score_conf)
        k = k + 1
    }
result

#1.c
reps = 10000
n_vec = c(10, 50, 500)
alpha_vec = c(0.01, 0.05)
score_conf = 0.99
mu = 68
sigma = 3
rho = 0.7
result = matrix(0, length(n_vec) * length(alpha_vec), 4)
colnames(result) = c("n", "alpha", "left", "right")
k = 1
for (n in n_vec)
    for (alpha in alpha_vec) {
        result[k, 1] = n
        result[k, 2] = alpha
        result[k, c(3, 4)] =
            sim_confint_sigma_sq_multinorm(reps, n, mu, sigma, rho, alpha, score_conf)
        k = k + 1
    }
result




# 2 -----------------------------------------------------------------------


#' simulate p-values from a t-test scenario where two samples are drawn from 
#' normal distribution
#'
#' @param reps number of replications for this simulation
#' @param n1 sample size of sample #1
#' @param n2 sample size of sample #2
#' @param mu1 mean value of sample #1
#' @param mu2 mean value of sample #2
#' @param sigma1 standard deviation of sample #1
#' @param sigma2 standard deviation of sample #2
#'
#' @return a vector of p-values calculated from a series of t-tests
sim_pval_ttest = function(reps, n1, n2, mu1, mu2, sigma1, sigma2) {
    x1 = rnorm(reps * n1, mu1, sigma1)
    x1 = matrix(x1, reps, n1)
    x2 = rnorm(reps * n2, mu2, sigma2)
    x2 = matrix(x2, reps, n2)
    p_values = numeric(reps)
    for (r in 1:reps) {
        test_result = t.test(x1[r,], x2[r,], var.equal = TRUE)
        p_values[r] = test_result$p.value
    }
    p_values
}

#' calculate the power of a test given a simulated p-value vector
#'
#' @param vector of simulated p-value based on the test 
#' @param significance level to determine the rejection region
#'
#' @return the power of the test
calc_power = function(pvals, alpha) {
    sum(pvals > (1 - alpha / 2) | pvals < (alpha / 2)) / length(pvals)
}

## 2.a
n1 = n2 = 50
mu1 = mu2 = 68
sigma1 = sigma2 = 3
reps = 1000

p_values = sim_pval_ttest(reps, n1, n2, mu1, mu2, sigma1, sigma2)
hist(p_values)


## 2.b
n1 = n2 = 50
mu1 = 68
sigma1 = sigma2 = 3
alpha = 0.05
reps = 1000

length_mu2 = 5
mu2_vec = seq(from = mu1, by = 1, length.out = length_mu2)
power_result = numeric(length_mu2)
for (i in 1:length_mu2) {
    p_values = sim_pval_ttest(reps, n1, n2, mu1, mu2_vec[i], sigma1, sigma2)
    power_result[i] = calc_power(p_values, alpha)
}
power_result[power_result < 1]
mu2_vec[power_result < 1]

## 2.c
mu1 = 68
mu2 = 68.5
sigma1 = sigma2 = 3
alpha = 0.05
reps = 1000

length_n = 10
n_vec = seq(from = 1050, by = 20, length.out = length_n)
power_result = numeric(length_n)
for (i in 1:length_n) {
    p_values = sim_pval_ttest(reps, n_vec[i], n_vec[i], mu1, mu2, sigma1, sigma2)
    power_result[i] = calc_power(p_values, alpha)
}
n_vec[which.min(abs(power_result - 0.95))]

## 2.d

#' Draw qqplot (modified from Adam's code)
#'
#' @param x.list target sample
#' @param quant.func quantile function for the target distribution
#' @param ... parameters to be passed to the quant.func()
#'
#' @return no return value
adam.qqplot = function (x.list, quant.func, ...) {
    n <- length(x.list)
    probs <- ppoints(n)
    plot(
        quant.func(probs, ...), quantile(x.list, probs),
        xlab = "Expected percentile", ylab = "Data percentile",
        main = "Q-Q Plot"
    )
    abline(0, 1)
}

#2.d.i

n1 = n2 = 100
mu1 = mu2 = 68
sigma1 = 3
sigma2 = 6
reps = 1000

p_values = sim_pval_ttest(reps, n1, n2, mu1, mu2, sigma1, sigma2)
par(mfrow = c(1, 2))
adam.qqplot(p_values, qunif)
hist(p_values)
quantile(p_values, c(0.01, 0.05))

#2.d.ii

n1 = 20
n2 = 100
mu1 = mu2 = 68
sigma1 = 3
sigma2 = 6
reps = 1000

p_values = sim_pval_ttest(reps, n1, n2, mu1, mu2, sigma1, sigma2)
par(mfrow = c(1, 2))
adam.qqplot(p_values, qunif)
hist(p_values)
quantile(p_values, c(0.01, 0.05))

#2.d.iii

n1 = 100
n2 = 20
mu1 = mu2 = 68
sigma1 = 3
sigma2 = 6
reps = 1000

p_values = sim_pval_ttest(reps, n1, n2, mu1, mu2, sigma1, sigma2)
par(mfrow = c(1, 2))
adam.qqplot(p_values, qunif)
hist(p_values)
quantile(p_values, c(0.01, 0.05))

##2.e

#' simulate p-values in a t-test scenario where data are generated from a 
#' multivariate normal distribution
#'
#' @param reps number of replications in the simulation study
#' @param n sample size for each group
#' @param mu1 marginal mean for group 1
#' @param mu2 marginal mean for group 2
#' @param sigma marginal standard deviation for both groups
#' @param rho decay parameter for the covariance between two grouops
#'
#' @return a vector of simulated p-values
sim_pval_ttest_multinorm = function(reps, n, mu1, mu2, sigma, rho) {
    mu_vec = c(mu1, mu2)
    var_mat = matrix(0, 2, 2)
    var_mat[1, 1] = var_mat[2, 2] = sigma ^ 2
    var_mat[1, 2] = var_mat[2, 1] = sigma ^ 2 * rho
    p_values = numeric(reps)
    for (r in 1:reps) {
        data_xy = mvrnorm(n, mu_vec, var_mat)
        test_result = t.test(data_xy[, 1], data_xy[, 2], var.equal = TRUE)
        p_values[r] = test_result$p.value
    }
    p_values
}

#2.e.1

n = 20
mu1 = mu2 = 68
sigma = 3
rho = 0.01
reps = 1000

p_values = sim_pval_ttest_multinorm(reps, n, mu1, mu2, sigma, rho)
par(mfrow = c(1, 2))
adam.qqplot(p_values, qunif)
hist(p_values)
quantile(p_values, c(0.01, 0.05))

#2.e.2

n = 20
mu1 = mu2 = 68
sigma = 3
rho = 0.1
reps = 1000

p_values = sim_pval_ttest_multinorm(reps, n, mu1, mu2, sigma, rho)
par(mfrow = c(1, 2))
adam.qqplot(p_values, qunif)
hist(p_values)
quantile(p_values, c(0.01, 0.05))


# 3 -----------------------------------------------------------------------


#' simulate the mean square errors for likelihood estimator and shrinkage 
#' estimator for the mean of a Exponential distribution
#'
#' @param reps number of replications for this simulation
#' @param n sample size
#' @param mu mean of the Exponential distribution
#'
#' @return a vector containing simulated mean square errors for the two estimators
#' (likelihood, shrinkage)
sim_diff_lkh_shk = function(reps, n, mu) {
    x = rexp(reps * n, 1 / mu)
    x = matrix(x, reps, n)
    a = n / (n + 1)
    erros = apply(x, 1, function(x, mu) {
        c((mean(x) - mu) ^ 2, (a * mean(x) - mu) ^ 2)
    }, mu)
    apply(erros, 1, mean)
}

n_vec = c(5, 10, 50)
mu_vec = c(0.5, 1, 10)
reps = 1000
result = matrix(0, length(n_vec) * length(mu_vec), 6)
colnames(result) = c("n", "mu", "sim_lkh_err", "sim_shk_est",
                     "lkh_err", "shk_est")
k = 1
for (n in n_vec)
    for (mu in mu_vec) {
        result[k, 1] = n
        result[k, 2] = mu
        result[k, c(3, 4)] = sim_diff_lkh_shk(reps, n, mu)
        result[k, 5] = mu ^ 2 / n
        result[k, 6] = mu ^ 2 / (n + 1)
        k = k + 1
    }
round(result, 3)
