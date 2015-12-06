require(glmnet)
require(ncvreg)
require(MASS)

lasso.feature <- function(X, y) {
  n = nrow(X)
  p = ncol(X)
  m.cv = cv.glmnet(X, y, nfolds = 5)
  lambda.best = m.cv$lambda.min
  m = glmnet(X, y, lambda = lambda.best)
  m.coef = as.numeric(coef(m))[2:(p + 1)]
  coef.nonzero = 1 * (m.coef > 1e-07 | m.coef < -1e-07)
  return(coef.nonzero)
}

ncv.feature <- function(X, y, penalty) {
  n = nrow(X)
  p = ncol(X)
  m.cv = cv.ncvreg(X, y, penalty = penalty, nfolds = 5)
  lambda.best = m.cv$lambda.min
  m = ncvreg(X, y, penalty = penalty, lambda = lambda.best)
  m.coef = as.numeric(coef(m))[2:(p + 1)]
  coef.nonzero = 1 * (m.coef > 1e-07 | m.coef < -1e-07)
  return(coef.nonzero)
}

mcp.feature <- function(X, y) {
  return(ncv.feature(X, y, "MCP"))
}

scad.feature <- function(X, y) {
  return(ncv.feature(X, y, "SCAD"))
}

stepwise.feature <- function(X, y, k) {
  n = nrow(X)
  p = ncol(X)
  data.set = data.frame(y, X)
  m.lm <- lm(y ~ ., data = data.set)
  m <- step(m.lm, trace = 0, k = k)
  n.features = length(coef(m))
  feature.index = as.numeric(gsub("X", "", names(coef(m))[2:n.features]))
  coef.nonzero = numeric(p)
  coef.nonzero[feature.index] = 1
  return(coef.nonzero)
}

stepwise.aic.feature <- function(X, y) {
  return(stepwise.feature(X, y, 2))
}

stepwise.bic.feature <- function(X, y) {
  return(stepwise.feature(X, y, log(nrow(X))))
}

plot.feature.freq <- function(beta, X, func.feature, nreps, ...) {
  n = nrow(X)
  p = ncol(X)
  coef.freq = numeric(p)
  for (r in 1:nreps) {
    y = X %*% beta + rnorm(n)
    coef.nonzero = func.feature(X, y)
    coef.freq = coef.freq + coef.nonzero
  }
  barplot(coef.freq, col = (beta != 0), 
          ylab = "selection frequency", xlab = "features", ...)
  return(coef.freq)
}

