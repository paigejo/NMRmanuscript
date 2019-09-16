# get PC prior for the BYM2 and SPDE total variance
library(INLA)
set.seed(123)
N = 10000
U = 1
prec = inla.pc.rprec(alpha = 0.01, u = U, n = N)
sig = 1/sqrt(prec)
z = rnorm(N, mean = rep(0, N), sd = sig)
effect = exp(z)
quantile(effect, c(0.025, 0.975))
#      2.5%     97.5% 
# 0.5016154 1.8766063 

set.seed(123)
N = 10000
U = .15
prec = inla.pc.rprec(alpha = 0.01, u = U, n = N)
sig = 1/sqrt(prec)
z = rnorm(N, mean = rep(0, N), sd = sig)
effect = exp(z)
quantile(effect, c(0.025, 0.975))
#        2.5%     97.5% 
#   0.9016866 1.0990210 

# get PC prior for the cluster effect variance
# Mortality = 0.04 # Kenya's first year mortality rate is approximately 4% for 2005-2009
# mu = log(Mortality/(1-Mortality))
# N = 1000000
# U = 3
# alpha = 0.01
# prec = inla.pc.rprec(alpha = alpha, u = U, n = N)
# sig = 1/sqrt(prec)
# z = rnorm(N, mean = rep(mu, N), sd = sig)
# p = 1/(1+exp(-z))
# y = matrix(rbinom(n = 2*N, size = 1, prob = rep(p,each = 2)), ncol = 2, byrow = T)
# print(cor(y)[1,2])
# # 0.1130728

# Kenya's secondary education completion rate is approximately 20% and 40% for 
# urban and rural areas respectively 2010-2014, women aged 20-29
edRate = 0.4
mu = log(edRate/(1-edRate))
N = 1000000
U = 3
alpha = 0.01
prec = inla.pc.rprec(alpha = alpha, u = U, n = N)
sig = 1/sqrt(prec)
z = rnorm(N, mean = rep(mu, N), sd = sig)
p = 1/(1+exp(-z))
y = matrix(rbinom(n = 2*N, size = 1, prob = rep(p,each = 2)), ncol = 2, byrow = T)
print(cor(y)[1,2])
# [1] 0.1021545
# nearly the same for edRate = 0.4

edRates = seq(from=0.1, to=0.5, l=9)
clustVars = seq(from= 0.01, to=3, l=15)
iccs = matrix(nrow = length(edRates), ncol = length(clustVars))
for(i in 1:length(edRates)) {
  for(j in 1:length(clustVars)) {
    edRate = edRates[i]
    mu = log(edRate/(1-edRate))
    
    N = 1000000
    clustVar = clustVars[j]
    z = rnorm(N, mean = rep(mu, N), sd = sqrt(clustVar))
    p = expit(z)
    y = matrix(rbinom(n = 2*N, size = 1, prob = rep(p,each = 2)), ncol = 2, byrow = T)
    iccs[i,j] = cor(y)[1,2]
  }
}
gridVals = expand.grid(list(edRate=edRates, clustVar=clustVars))
breaksGrid = list(x=seq(0.05, 0.55, l=10), y= seq(0.005, 3.005, l=16))
quilt.plot(gridVals, c(iccs), grid=breaksGrid, nx=length(edRates), ny=length(clustVars), main="Intraclass correlation", xlab="SCR", ylab="Cluster Variance", ylim=c(0.01, 3))


# exp(-x*12) = 1.04625