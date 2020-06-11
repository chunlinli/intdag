## reference: Li, Shen, Pan "Simultaneous inference of directed relations with interventions"
## author: Chunlin Li (li000007@umn.edu)

## tlpnet: estimate a gaussian directed acyclic graph with specified constraints 
## X: data matrix
## A, Lambda: initial estimate, A must be a DAG !!!
## D: index of hypothesized edges, no penalty is imposed
## tau: threshold par in TLP
## mu: sparsity par
## rho: ADMM par
## tol, rel_tol: absolute/relative tolerance
## dc_max_iter, admm_max_iter: DC/ADMM max iteration
## trace_obj: boolean value, indicate whether to print the obj function after each DC iter

#intdag.pmle <- function(Y, X, V.init = matrix(0,ncol(X),ncol(Y)), pen.fac = rep(1,ncol(X)), 
#                        K, tau = 1e-1, gamma=0.3, tol = tau*1e-3, cd.maxit = 1e+4) {
    # parallel r package: foreach
    
    
#}


# cross-validation 
#cv.intdag.pmle <- function() {
#    # DC+CD+warm start
#    
#}

source("tlpreg.r")

intdag.pmle <- function(X, Y, tau, gamma) {
    p <- ncol(Y)
    q <- ncol(X)
    V <- apply(Y, 2, function(y){
        m <- tlpreg0(y=y, X=X, tau = tau, gamma=gamma) 
        })
    V
}


intdag.ans <- function(V, D=matrix(0,ncol(V),ncol(V)), CorrMat = NULL) {
    
    # replace with C code
    p <- ncol(V)
    q <- nrow(V)
    if (q < p) stop("No sufficient interventions: q < p.")
    
    x.idx <- 1:q
    y.idx <- 1:p
    Pi <- matrix(0,p,p) 
    Phi <- matrix(0,q,p) 
    V.abs <- abs(V)
    V.nz <- ifelse(V!=0,1,0)
    
    L0.col <- apply(as.matrix(V.nz[x.idx, y.idx]), 2, sum)
    B <- y.idx[L0.col == 0]
    if(length(B)>0 & is.null(CorrMat)) 
        stop(paste("The following node(s) have no intervention", B)) # Can output the index with no interventions.
    for (j in B) {
        l <- x.idx[which.max(CorrMat[x.idx,j])]
        V.nz[l,j] <- 1
    }
    
    L0.row <- apply(as.matrix(V.nz[x.idx, y.idx]), 1, sum)
    for(k in 1:p) {
        
        # A contains all indexes of X which intervenes a leaf.
        A <- x.idx[L0.row == min(L0.row[L0.row>0])]
        j0 <- y.idx[which.max(V.abs[A[1], y.idx])]
        A0 <- c() # A0 contains all indexes of X which intervenes leaf Y[j0].
        for(l in A) {
            j <- y.idx[which.max(V.abs[l, y.idx])]
            if (j == j0) {
                A0 <- c(A0,l)
            }
        }
        
        # The upper bound of length(A0)
        A0.max.length <- length(x.idx)-length(y.idx)+1
        if (length(A0) > A0.max.length) {
            # Can use p-values to sort.
            # Truncate.
            A0 <- A0[1:A0.max.length]
        }
        
        # X[A0] intervenes Y[j0]
        Phi[A0,j0] <- 1
        
        # If Y[i] is removed and V[l,i] != 0, then Y[j0] is ancestor of Y[i].
        for (i in 1:p) {
            for (l in A0) {
                if (! i%in%y.idx & V.nz[l,i]==1) {
                    Pi[j0, i] <- 1
                }
            }
        }
        
        # Remove (Y[j0], X[A0]) from graph.
        x.idx <- setdiff(x.idx, A0)
        y.idx <- setdiff(y.idx, j0)
        if (length(y.idx)==0 | length(x.idx)==0) break
        
        L0.col <- apply(as.matrix(V.nz[x.idx, y.idx]), 2, sum)
        B <- y.idx[L0.col == 0]
        if(length(B)>0 & is.null(CorrMat)) 
            stop(paste("The following node(s) have no intervention", B)) 
        for (j in B) {
            l <- x.idx[which.max(CorrMat[x.idx,j])]
            V.nz[l,j] <- 1
        }
        
        L0.row <- apply(as.matrix(V.nz[x.idx, y.idx]), 1, sum)
    }
    
    # Adjust ordering 
    
    list(Pi=Pi, Phi=Phi)
}



intdag.mle <- function(X, Y, tau, gamma) {

    n <- nrow(X)
    q <- ncol(X)
    p <- ncol(Y)

    X.m <- apply(X, 2, function(x) (x - mean(x))/sd(x))
    Y.m <- apply(Y, 2, function(y) (y - mean(y))/sd(y))
    CorrMat <- crossprod(X.m,Y.m)/n

    ## estimate V
    V <- intdag.pmle(Y=Y, X=X, tau=tau, gamma=gamma) 

    ## estimate Pi (ancestral relations) and Phi (interventional relations) 
    m <- intdag.ans(V=V, CorrMat=CorrMat)
    Pi <- m$Pi
    Phi <- m$Phi
    
    U <- matrix(0,p,p)
    W <- matrix(0,q,p)

    for(j in 1:p) {
        idx.y <- which(Pi[,j] != 0)
        idx.x <- which(Phi[,j]!= 0)
        if (length(idx.y)>0) {
            Z <- cbind(rep(1,n),Y[,idx.y],X[,idx.x])
            b <- as.numeric(solve(crossprod(Z), crossprod(Z,Y[,j])))[-1]
            U[idx.y,j] <- b[1:length(idx.y)]
            W[idx.x,j] <- b[-(1:length(idx.y))]
        } else {
            idx.x <- which(Phi[,j]!= 0)
            Z <- cbind(rep(1,n),X[,idx.x])
            W[idx.x,j] <- as.numeric(solve(crossprod(Z), crossprod(Z,Y[,j])))[-1]
        }
    }

    list(U = U, W = W)
}




# MLE inference 

#intdag.mle <- function(X,Y,V = NULL) {

#    if (is.null(V)) V <- pmle(X,Y)




    # out put V, D, phi, df, Lr, pval
    # warning if phi is not injective: pval may not be correct 
#    list(U=U, V=V, W=W)
#}




# Data perturbation inference
#intdag.dp <- function() {
    
#}


