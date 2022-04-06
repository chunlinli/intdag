#library(glmtlp)
library(ncvreg) 
# currently glmtlp has some numerical issues, use ncvreg for regression instead

v_estimation <- function(y, x, tau = 0.2 * sqrt(log(q) / n), model_selection = c("cv", "bic")) {
    p <- ncol(y)
    q <- ncol(x)
    if (q < p) stop("No sufficient interventions: q < p.")

    if (nrow(y) == nrow(x)) {
        n <- nrow(y)
    } else {
        stop("Dimensions of y and x do not match.")
    }

    v <- matrix(0, nrow = q, ncol = p)
    err_var <- rep(0, p)
    lambda <- rep(0, p)

    if (model_selection == "cv") {
        for (j in 1:p) {
            m <- cv.ncvreg(X = x, y = y[, j], nfolds = 5)
            v[, j] <- as.numeric(m$fit$beta[-1, m$min])
            err_var[j] <- m$fit$loss[m$min] / n
            lambda[j] <- m$lambda.min
        }
    } else {
        for (j in 1:p) {
            m <- ncvreg(X = x, y = y[, j])
            hbic <- m$loss + (log(n) + 2 * log(p)) * colSums(m$beta[-1, ] != 0)
            min <- which.min(hbic)
            v[, j] <- as.numeric(m$beta[-1, min])
            err_var[j] <- m$loss[min] / n
            lambda[j] <- m$lambda[min]
        }
    }

    list(v = v, err_var = err_var, lambda = lambda)
}

topological_order <- function(v, thresh=0.15) {
    p <- ncol(v)
    q <- nrow(v)
    if (q < p) stop("No sufficient interventions: q < p.")

    v[abs(v) < thresh] <- 0

    an_mat <- matrix(0, p, p)
    in_mat <- matrix(0, q, p)
    iv_mat <- matrix(0, q, p)

    removed_x <- rep(FALSE, q)
    removed_y <- rep(FALSE, p)

    v_abs <- abs(v)
    v_nz <- (v_abs != 0)

    # check if there is a primary variable without intervention.
    while (!all(removed_y)) {

        # leaf-instrument pairs
        iv_targets <- rowSums(as.matrix(v_nz[, !removed_y]))
        one <- min(iv_targets[iv_targets > 0 & !removed_x])
        leaf_iv <- which(!removed_x & iv_targets == one)
        if (length(leaf_iv) == 0) break
        leaf <- rep(NA, length(leaf_iv))
        leaf_iter <- 1
        for (l in leaf_iv) {
            j <- which(v_abs[l, ] == max(v_abs[l, !removed_y]))[1]
            iv_mat[l, j] <- 1
            in_mat[l, j] <- 1
            leaf[leaf_iter] <- j
            leaf_iter <- leaf_iter + 1
        }
        leaf <- unique(leaf)

        # leaf-noninstrument pairs
        for (j in leaf) {
            leaf_interventions <- which(!removed_x & v_abs[, j] != 0)
            leaf_noniv <- setdiff(leaf_interventions, leaf_iv)
            in_mat[leaf_noniv, j] <- 1
        }

        # ancestral relations
        for (j in leaf) {
            j_instrument <- which(iv_mat[, j] != 0)
            j_descendant <- vector("list", length = length(j_instrument))
            for (l in 1:length(j_instrument)) {
                l2 <- j_instrument[l]
                j_descendant[[l]] <- which(removed_y & v_abs[l2, ] != 0)
            }
            # j_descendant <- Reduce(intersect, j_descendant)
            j_descendant <- Reduce(union, j_descendant)
            an_mat[j, j_descendant] <- 1
        }

        # peeling-off
        removed_y[leaf] <- TRUE
        removed_x[leaf_iv] <- TRUE
    }

    # reconstruction of topological order
    an_mat <- (solve(diag(p) - an_mat) != 0) - diag(p)
    in_mat <- 1 * (in_mat %*% (diag(p) + an_mat) > 0)

    list(an_mat = an_mat, in_mat = in_mat, iv_mat = iv_mat)
}

causal_discovery <- function(y, x, an_mat, in_mat, tau = 0.2 * sqrt(log(q) / n), model_selection = c("cv", "bic")) {
    p <- ncol(y)
    q <- ncol(x)
    if (q < p) stop("No sufficient interventions: q < p.")

    if (nrow(y) == nrow(x)) {
        n <- nrow(y)
    } else {
        stop("Dimensions of y and x do not match.")
    }

    u <- matrix(0, nrow = p, ncol = p)
    w <- matrix(0, nrow = q, ncol = p)


    for (j in 1:p) {
        ancestor <- which(an_mat[, j] != 0)
        intervention <- which(in_mat[, j] != 0)

        if (length(ancestor) > 0) {
            z <- cbind(y[, ancestor], x[, intervention])

            if (model_selection == "cv") {
                m <- cv.ncvreg(X = z, y = y[, j], nfolds = 5)
                beta <- as.numeric(m$fit$beta[-1, m$min])
            } else {
                m <- ncvreg(X = z, y = y[, j])
                hbic <- m$loss + (log(n) + 2 * log(p)) * colSums(m$beta[-1, ] != 0)
                min <- which.min(hbic)
                beta <- as.numeric(m$beta[-1, min])
            }

            u[ancestor, j] <- beta[1:length(ancestor)]
            w[intervention, j] <- beta[(length(ancestor) + 1:length(intervention))]
        } else {
            if (model_selection == "cv") {
                m <- cv.ncvreg(X = x[, intervention], y = y[, j], nfolds = 5)
                beta <- as.numeric(m$fit$beta[-1, m$min])
            } else {
                m <- ncvreg(X = x[, intervention], y = y[, j])
                hbic <- m$loss + (log(n) + 2 * log(p)) * colSums(m$beta[-1, ] != 0)
                min <- which.min(hbic)
                beta <- as.numeric(m$beta[-1, min])
            }

            w[intervention, j] <- beta
        }
    }

    list(u = u, w = w)
}


