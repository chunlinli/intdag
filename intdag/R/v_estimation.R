
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
        for (j in seq_len(p)) {
            m <- cv.ncvreg(X = x, y = y[, j], nfolds = 5, lambda.min = .05)
            v[, j] <- as.numeric(m$fit$beta[-1, m$min])
            err_var[j] <- m$fit$loss[m$min] / n
            lambda[j] <- m$lambda.min
        }
    } else {
        for (j in seq_len(p)) {
            m <- ncvreg(X = x, y = y[, j], lambda.min = .05)
            hbic <- m$loss + (log(n) + 2 * log(p)) * colSums(m$beta[-1, ] != 0)
            min <- which.min(hbic)
            v[, j] <- as.numeric(m$beta[-1, min])
            err_var[j] <- m$loss[min] / n
            lambda[j] <- m$lambda[min]
        }
    }

    list(v = v, err_var = err_var, lambda = lambda)
}

v_estimation_internal <- function(y, x, v_init, lambda) {
    p <- ncol(y)
    q <- ncol(x)
    v <- matrix(0, nrow = q, ncol = p)

    for (j in seq_len(p)) {
        m <- ncvfit(X = x, y = y[, j], init = v_init[, j], lambda = lambda[j])
        v[, j] <- as.numeric(m$beta[-1])
    }
    list(v = v)
}
