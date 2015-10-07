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