#dyn.load("./src/glmtlp")

# remove the initial guess of b.

tlpreg <- function(y, X, lambda=NULL, tau=0.3 * sqrt(log(p)/n), 
                   wobs = rep(1.0,n), pen.fac=rep(1.0,ncol(X)), 
                   delta = 2.0, tol=1e-4, dc.maxit=50, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    
    r <- y - mean(y)
    if (is.null(lambda)) {
        lambda.max <- max(abs(crossprod(X, r))) / n
        lambda <- exp(seq(
            from = log(lambda.max),
            to = log(0.05 * lambda.max),  
            length.out = 100
        ))
    }
    nlambda <- as.integer(length(lambda))

    b0 <- rep(mean(y), nlambda)
    b <- matrix(0, p, nlambda)

    # call regularized version
    .Call("rtlp", b0, b, r, X, wobs, pen.fac, lambda, nlambda, tau, 
          n, p, delta, tol, as.integer(dc.maxit), as.integer(cd.maxit))

    list(b = b, b0 = b0, lambda = lambda, tau = tau)
}


l0reg <- function(y, X, s=NULL, lambda=NULL, tau=0.3 * sqrt(log(p)/n), 
                   wobs = rep(1.0,n), pen.fac=rep(1.0,ncol(X)), 
                   delta = 2.0, tol=1e-4, dc.maxit=20, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))

    # initialize K sequence
    if(is.null(s))
        s <- 1:min(p, as.integer(n/log(p)))
    ns <- as.integer(length(s))
    
    r <- y - mean(y)
    if (is.null(lambda)) {
        lambda.max <- max(abs(crossprod(X, r))) / n
        lambda <- exp(seq(
            from = log(lambda.max),
            to = log(ifelse(n>p,.001,.05) * lambda.max),  
            length.out = 100
        ))
    }
    nlambda <- as.integer(length(lambda))

    b0 <- rep(mean(y), ns)
    b <- matrix(0, p, ns)

    # call regularized version
    .Call("rl0", b0, b, r, X, wobs, pen.fac, s, ns, lambda, nlambda, tau, 
          n, p, delta, tol, as.integer(dc.maxit), as.integer(cd.maxit))

    list(b = b, b0 = b0, s = s, lambda = lambda, tau = tau)
}
