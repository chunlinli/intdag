
causal_discovery <- function(y, x, an_mat, in_mat,
                             model_selection = c("cv", "bic"),
                             penalty = c("lasso", "MCP")) {
    model_selection <- match.arg(model_selection)
    penalty <- match.arg(penalty)

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
                m <- cv.ncvreg(
                    X = z, y = y[, j], penalty = penalty,
                    nfolds = 5, lambda.min = .05
                )
                beta <- as.numeric(m$fit$beta[-1, m$min])
            } else {
                m <- ncvreg(
                    X = z, y = y[, j],
                    penalty = penalty, lambda.min = .05
                )
                hbic <- m$loss +
                    (log(n) + 2 * log(p)) * colSums(m$beta[-1, ] != 0)
                min <- which.min(hbic)
                beta <- as.numeric(m$beta[-1, min])
            }

            u[ancestor, j] <- beta[seq_len(ancestor)]
            w[intervention, j] <- beta[
                (length(ancestor) + seq_len(intervention))
            ]
        } else {
            if (model_selection == "cv") {
                m <- cv.ncvreg(
                    X = x[, intervention], y = y[, j], penalty = penalty,
                    nfolds = 5, lambda.min = .05
                )
                beta <- as.numeric(m$fit$beta[-1, m$min])
            } else {
                m <- ncvreg(
                    X = x[, intervention], y = y[, j],
                    penalty = penalty, lambda.min = .05
                )
                hbic <- m$loss +
                    (log(n) + 2 * log(p)) * colSums(m$beta[-1, ] != 0)
                min <- which.min(hbic)
                beta <- as.numeric(m$beta[-1, min])
            }

            w[intervention, j] <- beta
        }
    }

    list(u = u, w = w)
}
