rm(list = ls())
source("./helper.R")

n = 500
p = 50
mu = rep(0, p)
Sigma = matrix(0, nrow = p, ncol = p)
diag(Sigma) = 1
X = mvrnorm(n, mu, Sigma)
beta = numeric(p)
beta[(1:p) %% 5 != 0] = 0
beta[(1:p) %% 5 == 0] = 1
nreps = 1000

freq.lasso = plot.feature.freq(beta, X, lasso.feature, nreps, main = "Simulated Feature Selection Frequency for LASSO Regression")
freq.mcp = plot.feature.freq(beta, X, mcp.feature, nreps, main = "Simulated Feature Selection Frequency for MCP Regression")
freq.scad = plot.feature.freq(beta, X, scad.feature, nreps, main = "Simulated Feature Selection Frequency for SCAD Regression")
freq.aic = plot.feature.freq(beta, X, stepwise.aic.feature, nreps, main = "Simulated Feature Selection Frequency for Stepwise AIC Regression")
freq.bic = plot.feature.freq(beta, X, stepwise.bic.feature, nreps, main = "Simulated Feature Selection Frequency for Stepwise BIC Regression")

save(freq.lasso, freq.mcp, freq.scad, freq.aic, freq.bic, file = "./sim_freq.RData")
