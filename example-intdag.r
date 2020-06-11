
source("intdag.r")


## hub graph
library(mvtnorm)
n <- 300
p <- 20
q <- 20

U <- matrix(0,p,p)
U[1,-1] <- 1
W <- matrix(0,nrow = q, ncol = p)
diag(W)[1:p] <- 1

X <- matrix(NA,nrow = n, ncol = q)
for(i in 1:n) {
    z <- rnorm(1)
    X[i,] <- rbinom(q, 1, exp(z)/(1+exp(z)))
    X[i,] <- ifelse(X[i,]==1,1,-1)
}
Y <- (X%*%W + rmvnorm(n,sigma=diag(seq(from=0.5,to=1,length.out=p),p,p)))%*%solve(diag(p) - U)

m <- intdag.mle(Y=Y,X=X,tau=0.4,gamma=1)
m$U - U
m$W - W




## random graph
library(mvtnorm)
set.seed(2020)
n <- 300
p <- 20
q <- 20
sparsity <- 1/p
U <- matrix(rbinom(p*p,1,sparsity)*sign(runif(p*p,min=-1,max=1))*1,p,p)
U[lower.tri(U, diag = T)] <- 0
W <- matrix(0,nrow = q, ncol = p)
diag(W)[1:p] <- 1

X <- matrix(NA,nrow = n, ncol = q)
for(i in 1:n) {
    z <- rnorm(1)
    X[i,] <- rbinom(q, 1, exp(z)/(1+exp(z)))
    X[i,] <- ifelse(X[i,]==1,1,-1)
}
Y <- (X%*%W + rmvnorm(n,sigma=diag(seq(from=0.5,to=1,length.out=p),p,p)))%*%solve(diag(p) - U)

m <- intdag.mle(Y=Y,X=X,tau=0.4,gamma=1)
m$U - U
m$W - W
