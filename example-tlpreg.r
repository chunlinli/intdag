source("tlpreg.r")

library(ncvreg)
library(picasso)

# example
n <- 1000
p <- 10000
X = matrix(rnorm(n*p),n,p)
y = X[,1] + X[,2] + rnorm(n)
for(j in 1:p) {
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}

t0 <- Sys.time()
b1 <- picasso(X=X, Y=y, lambda = 0.01, method = "mcp", standardize = F, intercept = F)$beta
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
b2 <- tlpreg0(y, X, tau=0.05, gamma=0.2, tol=1e-4)
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
b3 <- as.numeric(ncvreg(X,y,lambda=0.01)$beta[-1])
t1 <- Sys.time()
t1-t0
b1[1:20]
b2[1:20]
b3[1:20]

## example 
X = matrix(rnorm(200*50),200,50)
y = X[,1] + X[,4]*0.7 + rnorm(200)
for(j in 1:20) {
    X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
y <- y - mean(y)

library(glmnet) 
as.numeric(glmnet(X,y,lambda=0.5, intercept = FALSE, standardize = FALSE)$beta)
lasso(y,X,lambda=0.5,tol=1e-7)
cv.lasso(y,X)
cv.glmnet(X,y)$lambda.min
tlpreg(y,X)


library(ncvreg)
n <- 500
p <- 10
X = matrix(rnorm(n*p),n,p)
y = X[,1] + X[,4]*0.7 + rnorm(n)
for(j in 1:20) {
    X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
y <- y - mean(y)
tlpreg(y,X)
as.numeric(ncvreg(X,y,lambda = 0.2)$beta)[-1]












library(ncvreg)
library(glmnet)

n <- 400
p <- 15000
X = matrix(rnorm(n*p),n,p)
y = X[,1:10]%*%rep(1,10) + rnorm(n)
for(j in 1:p) {
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
y <- y - mean(y)

t0 <- Sys.time()
b1 <- as.numeric(glmnet(X,y,lambda=0.001, intercept = FALSE, standardize = FALSE)$beta)
t1 <- Sys.time()
t1-t0

b1[1:20]
sum(b1!=0)

t0 <- Sys.time()
b2 <- lasso(y,X,lambda=0.005,tol=5e-5)
t1 <- Sys.time()
t1-t0

b2[1:20]
sum(b2!=0)

t0 <- Sys.time()
b3 <- as.numeric(ncvreg(X,y,lambda=0.0005, penalty = "lasso")$beta[-1])
t1 <- Sys.time()
t1-t0

b3[1:20]
