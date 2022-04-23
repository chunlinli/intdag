\name{glmtlp-package}
\alias{glmtlp}
\docType{package}
\title{
\packageTitle{glmtlp}
}
\description{
\packageDescription{glmtlp}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{glmtlp}
\packageIndices{glmtlp}
}
\author{
\packageAuthor{glmtlp}
}
\keyword{ package, 
          lasso }
\examples{
set.seed(2021)
n <- 200
p <- 50
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
y <- 1 + 0.5*(X[,1] - X[,2] + X[,10] - X[,50]) + rnorm(n) 
m3 <- lasso(X=X, y=y, delta=2.0, tol = min(7e-5, sum((y - mean(y))^2)/2*(1e-7)))
}
