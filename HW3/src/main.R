rm(list = ls())
require(MASS)



# p2 ----------------------------------------------------------------------

# parameter Settings
n = 20
p = 5
beta_star = c(1, 1, 2, 1,-2)
sigma_star = 1

# make covariance matrix for generating the design matrix
rho = 0.5
Sigma = matrix(rho, nrow = p - 1, ncol = p - 1)
for (i in 1:(p - 1))
  Sigma[i, i] = 1

# determine the confidence level of the resulting confidence interval
alpha = 0.05
Z = qnorm(1 - alpha / 2)

# new observation, we will calculate confidence intervel for x_0 * beta_{*}
x_0 = c(1, 0, 1, 0, 1)
y_0 = crossprod(x_0, beta_star)

# randomly generate the Design Matrix, intercept inclusive
X = cbind(1, mvrnorm(
  n = n, mu = rep(0, p - 1), Sigma = Sigma
))
std_xbeta = sigma_star * sqrt(t(x_0) %*% solve(crossprod(X)) %*% x_0)
qr_X = qr(X)

# Simulate the confidence interval for x_0 * beta_{*}
# In each iteration, randomly generate y, estimate beta_hat, derive the
# confidence interval for x_{0} * beta_{*} and record if the it
# covers the true value y_{0}
n_reps = 10000
result = rep(FALSE, n_reps)
for (i in 1:n_reps) {
  y = X %*% beta_star + rnorm(n) * sigma_star
  beta_hat = qr.coef(qr_X, y)
  y_hat = crossprod(x_0, beta_hat)
  lower_bound = y_hat - Z * std_xbeta
  upper_bound = y_hat + Z * std_xbeta
  if (lower_bound < y_0 && y_0 < upper_bound)
    result[i] = TRUE
}

# report the estimated coverage probability and
# a 99% score confidence interval for the coverage probability
mean(result)
prop.test(sum(result), n_reps, conf.level = 0.99, correct = FALSE)



# p3 -----------------------------------------------------------------------

# wrapper function for simulating the confidence interval of the coverage probabiliy
sim_exp_res = function(n, beta_star, x_0, alpha, n_reps) {
  # parameter Settings
  p = length(beta_star)
  sigma_star = 1
  
  # make covariance matrix for generating the design matrix
  rho = 0.5
  Sigma = matrix(rho, nrow = p - 1, ncol = p - 1)
  for (i in 1:(p - 1))
    Sigma[i, i] = 1
  
  # determine the confidence level of the resulting confidence interval
  Z = qnorm(1 - alpha / 2)
  
  # new observation, we will calculate confidence intervel for x_0 * beta_{*}
  y_0 = crossprod(x_0, beta_star)
  
  # randomly generate the Design Matrix, intercept inclusive
  X = cbind(1, mvrnorm(
    n = n, mu = rep(0, p - 1), Sigma = Sigma
  ))
  std_xbeta = sigma_star * sqrt(t(x_0) %*% solve(crossprod(X)) %*% x_0)
  qr_X = qr(X)
  
  # Simulate the confidence interval for x_0 * beta_{*}
  # In each iteration, randomly generate y, estimate beta_hat, derive the
  # confidence interval for x_{0} * beta_{*} and record if the it
  # covers the true value y_{0}
  result = rep(FALSE, n_reps)
  for (i in 1:n_reps) {
    y = X %*% beta_star + rexp(n) - 1
    beta_hat = qr.coef(qr_X, y)
    y_hat = crossprod(x_0, beta_hat)
    lower_bound = y_hat - Z * std_xbeta
    upper_bound = y_hat + Z * std_xbeta
    if (lower_bound < y_0 && y_0 < upper_bound)
      result[i] = TRUE
  }
  
  interval = prop.test(sum(result), n_reps, conf.level = 0.99, correct = FALSE)$conf.int
  covers = interval[1] < (1 - alpha) &&  (1 - alpha) < interval[2]
  return (c(interval, covers))
}

beta_star = c(1, 1, 2, 1,-2)
x_0 = c(1, 0, 1, 0, 1)
n_reps = 10000

n_vec = c(10, 50, 200)
alpha_vec = c(0.01, 0.05)

i = 1
result = matrix(
  data = 0, nrow = length(n_vec) * length(alpha_vec), ncol = 5
)
colnames(result) = c("n", "alpha", "lower", "upper", "covers")
for (n in n_vec)
  for (alpha in alpha_vec) {
    result[i, 1] = n
    result[i, 2] = alpha
    result[i, 3:5] = sim_exp_res(n, beta_star, x_0, alpha, n_reps)
    i = i + 1
  }

result


# p4 ----------------------------------------------------------------------


# using cross validation to estimate the mean squared prediction error
# for a linear regression model
# X: Design Matrix
# y: Response Vector
# K: Number of folds
olscv <- function(X, y, K) {
  n <- length(y)
  
  # determine the cross validation index for each observation
  index = sample(1:n) %% K
  
  # sum of squared prediction erros
  SE = 0
  
  for (i in 0:(K - 1)) {
    # determine the train and test set
    X_train = X[index != i,]
    y_train = y[index != i]
    X_test = X[index == i,]
    y_test = y[index == i]
    
    # train the model
    m.lm <- lm.fit(x = X_train, y = y_train)
    beta = m.lm$coefficients
    
    # make prediction
    y_predict = X_test %*% beta
    
    # update the squared prediction errors
    SE = SE + sum((y_test - y_predict) ^ 2)
  }
  
  # return the mean-squared prediction error
  return(SE / n)
}



# p5ai ----------------------------------------------------------------------

rm(list = ls())

# function for simulating the power of the F test
# H0: beta_4 = 0, H1: otherwise
sim_power = function(n, beta_4, rho, alpha, reps) {
  # parameter settings
  beta_star = c(1, 1, 1, beta_4)
  sigma_star = 1
  p = length(beta_star)
  d = 1
  
  # Sigma matrix for generating Design Matrix X
  Sigma = matrix(rho, nrow = p - 1, ncol = p - 1)
  for (i in 1:(p - 1))
    Sigma[i, i] = 1
  
  # Randomly generate the design matrix X
  X = cbind(1, mvrnorm(
    n = n, mu = rep(0, p - 1), Sigma = Sigma
  ))
  
  # cache the qr decomposition result for X and X_null
  # where X_null is the one that without the last column (variable for beta_4)
  qr_X_null = qr(X[, 1:3])
  qr_X = qr(X)
  
  # store the p_value for the F-test in the pval.list array
  pval.list = rep(0, reps)
  for (r in 1:reps) {
    # randomly generate y
    y = X %*% beta_star + sigma_star * rnorm(n)
    
    # for null model, calc rss0
    beta_null = qr.coef(qr_X_null, y)
    rss0 = sum((y - X[, 1:3] %*% beta_null) ^ 2)
    
    # for full model, calc rss1
    beta = qr.coef(qr_X, y)
    rss1 = sum((y - X %*% beta) ^ 2)
    
    # calc F statistic and store current p-value
    f = ((rss0 - rss1) / d) / (rss1 / (n - p))
    pval.list[r] = 1 - pf(f, d, n - p)
  }
  
  # return the simulated power
  return(mean(pval.list < alpha))
}


# parameter settings
alpha = 0.01
reps = 3000

# possible values for beta-4, rho and n
beta_4_vec = c(0.5, 1)
rho_vec = c(0.3, 0.9)
n_vec = c(seq(10, 30, 2),
          seq(40, 70, 3),
          seq(300, 400, 10))

# result matrix
result = matrix(0, nrow = length(beta_4_vec) * length(rho_vec), ncol = 4)
colnames(result) = c("beta_4", "rho", "n", "power")

# grid search
i = 1
for (beta_4 in beta_4_vec)
  for (rho in rho_vec) {
    # find the n that gives the minimum difference abs(power - 0.8)
    # use diff to store the current minimum difference
    # use selected_n and selected_power to store the n that gives the minimum
    # difference and the corresponding simulated power
    diff = 3
    selected_n = -1
    selected_power = -1
    for (n in n_vec) {
      power = sim_power(n, beta_4, rho, alpha, reps)
      if (abs(power - 0.8) < diff) {
        selected_n = n
        selected_power = power
        diff = abs(power - 0.8)
      }
    }
    # store result into the result matrix
    result[i,] = c(beta_4, rho, selected_n, selected_power)
    i = i + 1
  }

result


# p5aii ---------------------------------------------------------------------

## For a model with design matrix X and
## measured response vector y,
## the following function returns a vector, where the first
## entry is the realized AIC and the second entry is the
## realized BIC
get.ic = function(X, y) {
  m <- lm(y ~ X - 1)
  c(AIC(m), BIC(m))
}

## parameter set derived from part i
para_set = matrix(0, nrow = 4, ncol = 3)
colnames(para_set) = c("beta_4", "rho", "n")
para_set[1,] = c(0.5, 0.3, 49)
para_set[2,] = c(0.5, 0.9, 380)
para_set[3,] = c(1, 0.3, 22)
para_set[4,] = c(1, 0.9, 70)

## result table to be reported
result_table = cbind(para_set, matrix(0, nrow = 4, ncol = 3))
colnames(result_table)[4:6] = c("AIC", "BIC", "CV")

## iterate over each parameter set
for (i in 1:4) {
  # current parameter setting
  beta_4 = para_set[i, 1]
  rho = para_set[i, 2]
  n = para_set[i, 3]
  beta_star = c(1, 1, 1, beta_4)
  p = length(beta_star)
  sigma_star = 1
  K = 5
  
  # Sigma matrix for generating Design Matrix X
  Sigma = matrix(rho, nrow = p - 1, ncol = p - 1)
  for (ii in 1:(p - 1))
    Sigma[ii, ii] = 1
  
  # Randomly generate the design matrix X
  X = cbind(1, mvrnorm(
    n = n, mu = rep(0, p - 1), Sigma = Sigma
  ))
  
  # result matrix for each paramter setting
  reps = 10000
  result = matrix(FALSE, nrow = reps, ncol = 3)
  colnames(result) = c("AIC", "BIC", "CV")
  
  for (r in 1:reps) {
    # randomly generate y
    y = X %*% beta_star + sigma_star * rnorm(n)
    
    # calculate criterions for both full model and null model
    criterion_full = c(get.ic(X, y), olscv(X, y, 5))
    criterion_null = c(get.ic(X[, 1:(p - 1)], y), olscv(X[, 1:(p - 1)], y, 5))
    
    # determine if the criterion is able to select the correct model
    result[r,] = criterion_full < criterion_null
  }
  
  # update the result table
  result_table[i, 4:6] = apply(result, 2, mean)
}

result_table


# p5aiii ------------------------------------------------------------------

# current parameter setting
n = 100
beta_4 = 0
rho = 0.3
alpha = 0.01
beta_star = c(1, 1, 1, beta_4)
p = length(beta_star)
sigma_star = 1
K = 5

# Sigma matrix for generating Design Matrix X
Sigma = matrix(rho, nrow = p - 1, ncol = p - 1)
for (ii in 1:(p - 1))
  Sigma[ii, ii] = 1

# Randomly generate the design matrix X
X = cbind(1, mvrnorm(
  n = n, mu = rep(0, p - 1), Sigma = Sigma
))

# result matrix for each paramter setting
reps = 10000
result = matrix(FALSE, nrow = reps, ncol = 3)
colnames(result) = c("AIC", "BIC", "CV")

for (r in 1:reps) {
  # randomly generate y
  y = X %*% beta_star + sigma_star * rnorm(n)
  
  # calculate criterions for both full model and null model
  criterion_full = c(get.ic(X, y), olscv(X, y, 5))
  criterion_null = c(get.ic(X[, 1:(p - 1)], y), olscv(X[, 1:(p - 1)], y, 5))
  
  # determine if the criterion is able to select the correct model
  result[r,] = criterion_full < criterion_null
}

# update the result table
result_table = matrix(0, nrow = 1, ncol = 4)
colnames(result_table) = c("AIC", "BIC", "CV", "F-test")
result_table[1:3] = apply(result, 2, mean)
result_table[4] = sim_power(n, beta_4, rho, alpha, reps)
result_table



# p5b ---------------------------------------------------------------------

rm(list = ls())
sim_power = function(n, beta_4, rho, reps) {
  # current parameter setting
  alpha = 0.01
  beta_star = c(1, 1, 1, beta_4)
  p = length(beta_star)
  d = 2
  sigma_star = 1
  
  # Sigma matrix for generating Design Matrix X
  Sigma = matrix(rho, nrow = p - 1, ncol = p - 1)
  for (ii in 1:(p - 1))
    Sigma[ii, ii] = 1
  
  # Randomly generate the design matrix X
  X = cbind(1, mvrnorm(
    n = n, mu = rep(0, p - 1), Sigma = Sigma
  ))
  qr_X = qr(X)
  
  # Generate X_null from X
  X_null = cbind(1, apply(X[,2:4], 1, sum))
  qr_X_null = qr(X_null)
  
  pval.list = rep(0, reps)
  for (r in 1:reps) {
    # randomly generate y
    y = X %*% beta_star + sigma_star * rnorm(n)
    
    # for null model, calc rss0
    beta_null = qr.coef(qr_X_null, y)
    rss0 = sum((y - X_null %*% beta_null) ^ 2)
    
    # for full model, calc rss1
    beta = qr.coef(qr_X, y)
    rss1 = sum((y - X %*% beta) ^ 2)
    
    # calc F statistic and store current p-value
    f = ((rss0 - rss1) / d) / (rss1 / (n - p))
    pval.list[r] = 1 - pf(f, d, n - p)
  }
  
  # return the simulated power
  return(mean(pval.list < alpha))
}

beta_4_vec = c(1.5, 2.5)
rho_vec = c(0.3, 0.9)
n_vec = c(seq(17, 50, 3),
          seq(60, 120, 10),
          seq(800, 1000, 20))
reps = 3000

# result matrix
result = matrix(0, nrow = length(beta_4_vec) * length(rho_vec), ncol = 4)
colnames(result) = c("beta_4", "rho", "n", "power")

# grid search
i = 1
for (beta_4 in beta_4_vec)
  for (rho in rho_vec) {
    # find the n that gives the minimum difference abs(power - 0.8)
    # use diff to store the current minimum difference
    # use selected_n and selected_power to store the n that gives the minimum
    # difference and the corresponding simulated power
    diff = 3
    selected_n = -1
    selected_power = -1
    for (n in n_vec) {
      power = sim_power(n, beta_4, rho, reps)
      if (abs(power - 0.85) < diff) {
        selected_n = n
        selected_power = power
        diff = abs(power - 0.85)
      }
    }
    # store result into the result matrix
    result[i,] = c(beta_4, rho, selected_n, selected_power)
    i = i + 1
  }

result