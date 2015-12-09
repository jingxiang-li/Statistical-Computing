
# helper functions --------------------------------------------------------


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
          ylab = "Selection Frequency", xlab = "Features", 
          ylim = c(0, nreps), space = 0.2, ...)
  axis(1, labels = seq(5, 50, 5), at = seq(5, 50, 5) * 1.2)
  return(coef.freq)
}


# main function -----------------------------------------------------------


rm(list = ls())
source("./helper.R")
require(corrplot)

# scenario 1 --------------------------------------------------------------

n = 50
p = 50
nreps = 500
mu = rep(0, p)
Sigma = matrix(0, nrow = p, ncol = p)
diag(Sigma) = 1
X = mvrnorm(n, mu, Sigma)
beta = numeric(p)
beta[(1:p) %% 5 != 0] = 0
beta[(1:p) %% 5 == 0] = 1

cairo_pdf(filename = "./tex/fig/s1.cor.pdf", width = 6, height = 6)
corrplot(cor(X), method = "square", is.corr = TRUE, tl.pos = FALSE)
dev.off()

cairo_pdf(filename = "./tex/fig/s1.lasso.pdf", width = 8, height = 4)
freq.lasso = plot.feature.freq(beta, X, lasso.feature, nreps, main = "Scenario1, Feature Selection Frequency for LASSO Regression")
dev.off()
cairo_pdf(filename = "./tex/fig/s1.mcp.pdf", width = 8, height = 4)
freq.mcp = plot.feature.freq(beta, X, mcp.feature, nreps, main = "Scenario1, Feature Selection Frequency for MCP Regression")
dev.off()
cairo_pdf(filename = "./tex/fig/s1.scad.pdf", width = 8, height = 4)
freq.scad = plot.feature.freq(beta, X, scad.feature, nreps, main = "Scenario1, Feature Selection Frequency for SCAD Regression")
dev.off()

p1 = function(x) {
  sum(x[-seq(5, 50, 5)]) / (40 * nreps)
}
p2 = function(x) {
  sum(x[seq(5, 50, 5)]) / (10 * nreps)
}
c(p1(freq.lasso), p2(freq.lasso))
c(p1(freq.scad), p2(freq.scad))
c(p1(freq.mcp), p2(freq.mcp))



# scenario 2 --------------------------------------------------------------


n = 50
p = 50
nreps = 500
mu = rep(0, p)
Sigma = matrix(0, nrow = p, ncol = p)
rho = 0.9
for (i in 1 : p)
  for (j in 1 : p) {
    Sigma[i, j] = rho^abs(i - j)
  }

X = mvrnorm(n, mu, Sigma)
beta = numeric(p)
beta[(1:p) %% 5 != 0] = 0
beta[(1:p) %% 5 == 0] = 1

cairo_pdf(filename = "./tex/fig/s2.cor.pdf", width = 6, height = 6)
corrplot(cor(X), method = "square", is.corr = TRUE, tl.pos = FALSE)
dev.off()
cairo_pdf(filename = "./tex/fig/s2.lasso.pdf", width = 8, height = 4)
freq.lasso = plot.feature.freq(beta, X, lasso.feature, nreps, main = "Scenario2, Feature Selection Frequency for LASSO Regression")
dev.off()
cairo_pdf(filename = "./tex/fig/s2.mcp.pdf", width = 8, height = 4)
freq.mcp = plot.feature.freq(beta, X, mcp.feature, nreps, main = "Scenario2, Feature Selection Frequency for MCP Regression")
dev.off()
cairo_pdf(filename = "./tex/fig/s2.scad.pdf", width = 8, height = 4)
freq.scad = plot.feature.freq(beta, X, scad.feature, nreps, main = "Scenario2, Feature Selection Frequency for SCAD Regression")
dev.off()

p1 = function(x) {
  sum(x[-seq(5, 50, 5)]) / (40 * nreps)
}
p2 = function(x) {
  sum(x[seq(5, 50, 5)]) / (10 * nreps)
}
c(p1(freq.lasso), p2(freq.lasso))
c(p1(freq.scad), p2(freq.scad))
c(p1(freq.mcp), p2(freq.mcp))



# scenario 3 --------------------------------------------------------------


n = 50
p = 50
nreps = 500
mu = rep(0, p)
Sigma = matrix(1, nrow = p, ncol = p)
rho = 0.9
for (i in 1 : p)
  for (j in 1 : p) {
    if (i != j)
      Sigma[i, j] = rho
  }

X = mvrnorm(n, mu, Sigma)
beta = numeric(p)
beta[(1:p) %% 5 != 0] = 0
beta[(1:p) %% 5 == 0] = 1

cairo_pdf(filename = "./tex/fig/s3.cor.pdf", width = 6, height = 6)
corrplot(cor(X), method = "square", is.corr = TRUE, tl.pos = FALSE)
dev.off()

cairo_pdf(filename = "./tex/fig/s3.lasso.pdf", width = 8, height = 4)
freq.lasso = plot.feature.freq(beta, X, lasso.feature, nreps, main = "Scenario3, Feature Selection Frequency for LASSO Regression")
dev.off()
cairo_pdf(filename = "./tex/fig/s3.mcp.pdf", width = 8, height = 4)
freq.mcp = plot.feature.freq(beta, X, mcp.feature, nreps, main = "Scenario3, Feature Selection Frequency for MCP Regression")
dev.off()
cairo_pdf(filename = "./tex/fig/s3.scad.pdf", width = 8, height = 4)
freq.scad = plot.feature.freq(beta, X, scad.feature, nreps, main = "Scenario3, Feature Selection Frequency for SCAD Regression")
dev.off()

p1 = function(x) {
  sum(x[-seq(5, 50, 5)]) / (40 * nreps)
}
p2 = function(x) {
  sum(x[seq(5, 50, 5)]) / (10 * nreps)
}
c(p1(freq.lasso), p2(freq.lasso))
c(p1(freq.scad), p2(freq.scad))
c(p1(freq.mcp), p2(freq.mcp))