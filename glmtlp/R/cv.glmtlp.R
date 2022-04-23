# regularized version 
cv.tlpreg <- function(y, X, lambda = NULL, tau=0.3*sqrt(log(p)/n), 
                      wobs = rep(1.0,n), pen.fac = rep(1.0,ncol(X)), 
                       nfold=10, delta=2.0,
                      tol=1e-4, dc.maxit=50, cd.maxit=1e+4, ncores=1) {
    require(doParallel)
    registerDoParallel(cores = ncores)

    n <- nrow(X)
    p <- ncol(X)
    if (n < nfold) stop("nfold cannot be larger than n.")

    r <- y - mean(y)
    if (is.null(lambda)) {
        lambda.max <- max(abs(crossprod(X, r))) / n
        lambda <- exp(seq(
            from = log(lambda.max),
            to = log(0.05 * lambda.max),  
            length.out = 100
        ))
    }
    
    # nfold-cross-validation
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    cv <- foreach(fold = 1:nfold, .combine = "cbind", .packages=c("glmtlp")) %dopar% {
        b <- tlpreg(y=y[obs.fold!=fold], X=X[obs.fold!=fold,], lambda=lambda,
                    tau=tau, wobs = wobs, pen.fac = pen.fac, delta = 2.0,
                    tol = tol, dc.maxit = dc.maxit, cd.maxit = cd.maxit)$b
        
        apply(b, 2, function(beta) {
            res <- (y[obs.fold==fold] - X[obs.fold==fold,] %*% beta)
            mean((res - mean(res))^2)
        })
    }

    # output CV graph
    cvm <- rowMeans(cv)

    m <- tlpreg(y=y, X=X, lambda=lambda,
                    tau=tau, wobs = wobs, pen.fac = pen.fac, delta = 2.0,
                    tol = tol, dc.maxit = dc.maxit, cd.maxit = cd.maxit)

    list(m=m, lambda.min=lambda[which.min(cvm)], cvm = cvm, lambda=lambda, tau=tau)
}




cv.l0reg <- function(y, X, s=NULL, lambda = NULL, tau=0.3*sqrt(log(p)/n), 
                      wobs = rep(1.0,n), pen.fac = rep(1.0,ncol(X)), 
                       nfold=10, delta=2.0,
                      tol=1e-4, dc.maxit=20, cd.maxit=1e+4, ncores=1) {
    require(doParallel)
    registerDoParallel(cores = ncores)

    n <- nrow(X)
    p <- ncol(X)
    if (n < nfold) stop("nfold cannot be larger than n.")

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
            length.out = 50
        ))
    }
    
    # nfold-cross-validation
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    cv <- foreach(fold = 1:nfold, .combine = "cbind", .packages=c("glmtlp")) %dopar% {
        b <- l0reg(y=y[obs.fold!=fold], X=X[obs.fold!=fold,], s=s, lambda=lambda,
                    tau=tau, wobs = wobs, pen.fac = pen.fac, delta = delta,
                    tol = tol, dc.maxit = dc.maxit, cd.maxit = cd.maxit)$b
        
        apply(b, 2, function(beta) {
            res <- (y[obs.fold==fold] - X[obs.fold==fold,] %*% beta)
            mean((res - mean(res))^2)
        })
    }

    # output CV graph
    cvm <- rowMeans(cv)

    m <- l0reg(y=y, X=X, s=s, lambda=lambda,
                    tau=tau, wobs = wobs, pen.fac = pen.fac, delta = delta,
                    tol = tol, dc.maxit = dc.maxit, cd.maxit = cd.maxit)

    list(m=m, s.min=s[which.min(cvm)], cvm = cvm, s=s, tau=tau)
}