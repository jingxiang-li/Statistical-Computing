setwd("./src")
rm(list = ls())
set.seed(5701)

require(MASS)


# 1.a ---------------------------------------------------------------------

# ridge regression for a single rugularization paramter lam (lambda)
# X, design matrix, first column to be all 1
# y, response vector
# lam, regularization parameter lambda
lm.ridge = function(X, y, lam) {
  n = nrow(X)
  p = ncol(X)
  xbar = apply(X[,-1], 2, mean)
  ybar = mean(y)
  Xtilde = scale(X[,-1], center = xbar, scale = FALSE)
  ytilde = y - ybar
  q = min(c(n - 1, p - 1))
  out = svd(Xtilde, nu = q, nv = q)
  H = diag(out$d / (out$d ^ 2 + lam))
  bhatm1 = out$v %*% H %*% t(out$u) %*% ytilde
  bhat1 = ybar - sum(xbar * bhatm1)
  beta.hat = c(bhat1, bhatm1)
  return(beta.hat)
}


# 1.b ---------------------------------------------------------------------

# given a beta vector, calculate the gradient of a ridge regression set
# X, design matrix, first column all 1
# y, response vector
# lam, lambda
# beta, coefficient vector
lm.ridge.grad = function(X, y, lam, beta) {
  grad.loss = -2 * crossprod(X, y - X %*% beta)
  grad.reg = 2 * lam * beta
  grad.reg[1] = 0
  return(grad.loss + grad.reg)
}


n = 50
p = 5
sigma.star = 1

beta.star = rnorm(p, mean = 0, sd = sqrt(0.05))
X = cbind(1, matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1))
y = X %*% beta.star + rnorm(n) * sigma.star
lambda = 1

beta.hat = lm.ridge(X, y, lambda)
grad = lm.ridge.grad(X, y, lambda, beta.hat)
cat("lambda :", lambda,
    "\nbeta.hat :", beta.hat,
    "\ngradient :", grad)


# 1.c ---------------------------------------------------------------------

# given a vector of lambda's, return a coefficient matrix
# for each lambda value as a matrix
# X, design matrix, first column all 1
# y, response vector
# lam.vec, a vector of lambda value to be tuned
lm.ridge1 = function(X, y, lam.vec) {
  n = nrow(X)
  p = ncol(X)
  xbar = apply(X[,-1], 2, mean)
  ybar = mean(y)
  Xtilde = scale(X[,-1], center = xbar, scale = FALSE)
  ytilde = y - ybar
  q = min(c(n - 1, p - 1))
  out = svd(Xtilde, nu = q, nv = q)
  beta.hat.matrix = matrix(0, nrow = p, ncol = length(lam.vec))
  for (j in 1:length(lam.vec))
  {
    H = diag(out$d / (out$d ^ 2 + lam.vec[j]))
    bhatm1 = out$v %*% H %*% t(out$u) %*% ytilde
    bhat1 = ybar - sum(xbar * bhatm1)
    beta.hat.matrix[, j] = c(bhat1, bhatm1)
  }
  return(beta.hat.matrix)
}

# run cv to find the best lambda
# X, design matrix, first column all 1
# y, response vector
# lam.vec, a vector of lambda value to be tuned
# K, number of folds
lm.ridge.cv = function(X, y, lam.vec, K) {
  n = nrow(X)
  p = ncol(X)
  cv.err = rep(0, length(lam.vec))
  
  order_index = sample.int(n)
  for (k in 0:(K - 1)) {
    test.index = order_index %% K == k
    X.test = X[test.index,]
    y.test = y[test.index]
    X.train = X[!test.index,]
    y.train = y[!test.index]
    
    beta.mat = lm.ridge1(X.train, y.train, lam.vec)
    err.mat = y.test %*% t(rep(1, length(lam.vec))) -
      X.test %*% beta.mat
    
    cv.err = cv.err + apply(err.mat, 2, function(x) {
      sum(x ^ 2)
    })
  }
  
  cv.err = cv.err / n
  best.lam = lam.vec[which.min(cv.err)]
  
  result = list()
  result$b = lm.ridge(X, y, best.lam)
  result$best.lam = best.lam
  result$cv.error = cv.err
  return(result)
}


# 1.d ----------------------------------------------------------------------

n = 100
p = 50
theta = 0.7
lambda.vec = 10 ^ seq(-8, 8)
K.vec = c(0, 5, 10, n)

Sigma = matrix(0, nrow = p - 1, ncol = p - 1)
for (i in 1:(p - 1))
  for (j in 1:(p - 1)) {
    Sigma[i, j] = theta ^ abs(i - j)
  }

X = cbind(1, mvrnorm(n, mu = rep(0, p - 1), Sigma = Sigma))
beta.star = rnorm(p, mean = 0, sd = 1 / sqrt(p))
y.expected = X %*% beta.star

nreps = 200
loss1.mat = matrix(0, nrow = nreps, ncol = length(K.vec))
loss2.mat = matrix(0, nrow = nreps, ncol = length(K.vec))
for (r in 1:nreps) {
  # print(r)
  y = y.expected + rnorm(n)
  for (i in 1:length(K.vec)) {
    K = K.vec[i]
    if (0 == K) {
      beta.hat = lm.fit(X, y)$coefficients
    } else {
      m.ridge = lm.ridge.cv(X, y, lambda.vec, K)
      beta.hat = m.ridge$b
    }
    loss1.mat[r, i] = sum((beta.star - beta.hat) ^ 2)
    loss2.mat[r, i] = sum((y.expected - X %*% beta.hat) ^ 2)
  }
}

cat("mean :", c(apply(loss1.mat, 2, mean),
                apply(loss2.mat, 2, mean)))
cat("SE :", c(sqrt(apply(loss1.mat, 2, var) / nreps),
              sqrt(apply(loss2.mat, 2, var) / nreps)))


# 2.a ---------------------------------------------------------------------

# calculate the gradient of logistic regression
logis.grad = function(beta, X, y, n.list, lambda, t.list, ...) {
  prob.list = 1 / (1 + exp(-X %*% beta))
  grad.loss = -crossprod(X, n.list * (y - prob.list))
  
  beta1 = beta[2:length(beta)]
  grad.reg = c(0, lambda * (beta1 - t.list))
  return(grad.loss + grad.reg)
}

# train the logistic regression
glm.logis = function(X, y, n.list, lam, t.list, tol, maxit, ...) {
  n = nrow(X)
  p = ncol(X)
  beta = runif(p, -1, 1)

  for (iter in 1:maxit) {
    grad = logis.grad(beta, X, y, n.list, lam, t.list)
    eta = 1 / sqrt((lambda + 1) * iter)
    beta.new = beta - eta * grad
    beta = beta.new
    
    if (sum(abs(grad)) < tol)
      break
  }
  return(list(b = beta, total.iterations = iter))
}

n = 100
p = 10

n.list = rpois(n, 9) + 1
X = cbind(1, 
          matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1))
beta.star = rnorm(p, sd = sqrt(1 / p))
prob.list = 1 / (1 + exp(-X %*% beta.star))
y = numeric(n)
for (i in 1:n) {
  y[i] = rbinom(1, n.list[i], prob.list[i]) / n.list[i]
}
t.list = rep(0, p - 1)
lambda = 4
tol = 1e-07
maxit = 1e+04

m = glm.logis(X, y, n.list, lambda, t.list, tol, maxit)
logis.grad(m$b, X, y, n.list, lambda, t.list)
data.frame(beta.star, beta.hat = m$b)



