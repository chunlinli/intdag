## reference: Li, Shen, Pan "Simultaneous inference of directed relations with interventions"
## author: Chunlin Li (li000007@umn.edu)

dyn.load("ctlpreg.so")

lasso0 <- function(y, X, b.init=rep(0,ncol(X)), lambda, pen.fac=rep(1,ncol(X)), tol=1e-5, cd.maxit=1e+4) {

    n <- nrow(X)
    p <- ncol(X)
    .C('classo', y = as.double(y),
                 X = as.double(X),
                 b0 = as.double(numeric(1)),
                 b = as.double(b.init),
                 r = as.double(y),
                 n = as.integer(n),
                 p = as.integer(p),
                 lambda = as.double(lambda),
                 pen_fac = as.integer(pen.fac),
                 tol = as.double(tol),
                 cd_maxit = as.integer(cd.maxit))$b
}

tlpreg0 <- function(y, X, b.init=rep(0,ncol(X)), tau=0.01, gamma=0.5, pen.fac=rep(1,ncol(X)), tol=1e-5, dc.maxit=20, cd.maxit=1e+4) {

    n <- nrow(X)
    p <- ncol(X)
    .C('ctlpreg0', y = as.double(y),
                   X = as.double(X),
                   b0 = as.double(numeric(1)),
                   b = as.double(b.init),
                   r = as.double(y),
                   n = as.integer(n),
                   p = as.integer(p),
                   tau = as.double(tau),
                   gamma = as.double(gamma),
                   pen_fac = as.integer(pen.fac),
                   tol = as.double(tol),
                   dc_maxit = as.integer(dc.maxit),
                   cd_maxit = as.integer(cd.maxit))$b
}

lasso <- function(y, X, b.init=rep(0,ncol(X)), lambda, pen.fac=rep(1,ncol(X)), tol=1e-4, cd.maxit=1e+4) {

    # C code.
    p <- ncol(X)
    n <- nrow(X)
    b.curr <- b.init
    b.next <- b.curr
    for(it in 1:cd.maxit) {
        b0 <- mean(y - X%*%b.next)
        for(j in 1:p) {
            # use working residuals
            b.next[j] <- sum(X[,j] * (y - b0 - X[,-j]%*%b.next[-j]))
            b.next[j] <- sign(b.next[j])*max(0,(abs(b.next[j])-pen.fac[j]*lambda*n)/sum(X[,j]*X[,j]))
        }
        if(sum((b.curr - b.next)^2) <= tol^2) break
        b.curr <- b.next
    }
    if (it == cd.maxit) warning("The coordinate descent algorithm does not converge.")
    b.next
}

tlpreg <- function(y, X, b.init=rep(0,ncol(X)), 
                   K=NULL, tau=0.05, gamma=1, pen.fac=rep(1,ncol(X)), 
                   tol=1e-4, dc.maxit=20, cd.maxit=1e+4) {
    # C code.
    b.curr <- b.init
    b.next <- b.curr
    for(it in 1:dc.maxit) {
        b.next <- lasso(y=y, X=X, b.init=b.next, lambda=gamma*tau, 
                        pen.fac=pen.fac*as.numeric(abs(b.next) < tau), 
                        tol=tol, cd.maxit=cd.maxit)
        if(sum((b.curr - b.next)^2) <= tol^2) break
        b.curr <- b.next
    }
    if (it == dc.maxit) warning("The difference of convex algorithm does not converge.")
    if (is.null(K)) list(b0=mean(y - X%*%b.next), b=b.next)
    
    # truncate / relaxation

    list(b0=mean(y - X%*%b.next), b=b.next)
}


cv.lasso <- function(y, X, b.init = rep(0,ncol(X)), pen.fac = rep(1,ncol(X)), lambda = NULL, 
                     nfold=10, tol=1e-4, cd.maxit=1e+4) {  
    # tuning lambda
    p <- ncol(X)
    n <- nrow(X)
    if (n < nfold) stop("nfold cannot be larger than n.")
    if(is.null(lambda)) {
        lambda.max <- max(abs(t(X)%*%(y-mean(y))))/n
        lambda <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p<n, 1e-4, .01)*lambda.max), length.out=100))
    } 

    # nfold-cross-validation
    obs.fold <- rep(1:nfold,n)[sample(1:n,n)]
    cvm <- matrix(0,length(lambda),nfold)
    for(k in 1:nfold) {
        X.tr <- X[obs.fold!=k,]
        y.tr <- y[obs.fold!=k]
        X.te <- X[obs.fold==k,]
        y.te <- y[obs.fold==k]
        
        # warm starts CV: C code 
        for(j in 1:length(lambda)) {
            b <- lasso(y.tr, X.tr, b.init=b.init, lambda=lambda[j], tol=tol, cd.maxit=cd.maxit)
            b.init <- b
            cvm[j,k] <- mean((y.te - X.te%*%b)^2)
        }
    }
    
    # output CV graph
    cvm <- apply(cvm, 1, mean)
    list(lambda.min=lambda[which.min(cvm)], cvm = cvm, lambda=lambda)
}


cv.tlpreg <- function(y, X, b.init = rep(0,ncol(X)), pen.fac = rep(1,ncol(X)), 
                      tau=0.05, K=NULL, gamma=NULL, constrained=FALSE, 
                      nfold=10, tol=1e-4, dc.maxit=10, cd.maxit=1e+4) {
    
    # tune K and gamma 
    n <- nrow(X)
    p <- ncol(X)
    if (n < nfold) stop("nfold cannot be larger than n.")
    if(is.null(lambda)) {
        lambda.max <- max(abs(t(X)%*%(y-mean(y))))/n
        lambda <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p<n, 1e-4, .01)*lambda.max), length.out=100))
    } 
    
    # nfold-cross-validation
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    cvm <- matrix(0,length(lambda),nfold)
    for(k in 1:nfold) {
        X.tr <- X[obs.fold!=k,]
        y.tr <- y[obs.fold!=k]
        X.te <- X[obs.fold==k,]
        y.te <- y[obs.fold==k]
        
        # warm starts CV: C code 
        for(j in 1:length(lambda)) {
            b <- lasso(y.tr, X.tr, b.init=b.init, pen.fac = pen.fac, lambda=lambda[j], tol=tol, cd.maxit=cd.maxit)
            b.init <- b
            b <- tlpreg(y.tr,X.tr, b.init=b, tau=tau, pen.fac = pen.fac, lambda=lambda[j], 
                        tol=tol, dc.maxit=dc.maxit, cd.maxit = cd.maxit)$b
            cvm[j,k] <- mean((y.te - X.te%*%b)^2)
        }
    }
    
    # output CV graph
    cvm <- apply(cvm, 1, mean)
    list(lambda.min=lambda[which.min(cvm)], cvm = cvm, lambda=lambda)
}






