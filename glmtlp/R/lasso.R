
# remove the initial guess of b.

lasso <- function(y, X, lambda = NULL,
                  wobs = rep(1.0, n), pen.fac = rep(1.0, p),
                  delta = 2.0, tol = 1e-4, cd.maxit = 1e+4) {
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
    
    # call C interface: rlasso.cc
    .Call(
        "rlasso", b0, b, r, X, wobs, pen.fac, lambda, nlambda,
        n, p, delta, tol, as.integer(cd.maxit)
    )

    list(b = b, b0 = b0, lambda = lambda)
}

# b0 should be a vec: length(b0) == length(lambda)
