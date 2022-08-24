
causal_discovery_internal <- function(y, x, an_mat, in_mat, nlambda = 5,
                                      penalty = c("MCP", "lasso")) {
    penalty <- match.arg(penalty)

    p <- ncol(y)
    q <- ncol(x)
    if (q < p) stop("No sufficient interventions: q < p.")

    if (nrow(y) != nrow(x)) {
        stop("Dimensions of y and x do not match.")
    }

    u_list <- vector(mode = "list", length = nlambda)
    w_list <- vector(mode = "list", length = nlambda)
    lambda_list <- rep(0, nlambda)

    for (k in seq_len(nlambda)) {
        u_list[[k]] <- matrix(0, nrow = p, ncol = p)
        w_list[[k]] <- matrix(0, nrow = q, ncol = p)
    }

    for (j in seq_len(p)) {
        ancestor <- which(an_mat[, j] != 0)
        intervention <- which(in_mat[, j] != 0)

        if (length(ancestor) > 0) {
            z <- cbind(y[, ancestor], x[, intervention])

            m <- ncvreg(
                X = z, y = y[, j],
                penalty = penalty, lambda.min = .001, nlambda = nlambda,
            )

            lambda_list <- lambda_list + m$lambda

            for (k in seq_len(nlambda)) {
                beta <- as.numeric(m$beta[-1, k])
                u_list[[k]][ancestor, j] <- beta[seq_len(length(ancestor))]
                w_list[[k]][intervention, j] <- beta[
                    (length(ancestor) + seq_len(length(intervention)))
                ]
            }
        } else {
            m <- ncvreg(
                X = x[, intervention], y = y[, j],
                penalty = penalty, lambda.min = .001, nlambda = nlambda,
            )

            lambda_list <- lambda_list + m$lambda

            for (k in seq_len(nlambda)) {
                beta <- as.numeric(m$beta[-1, k])
                w_list[[k]][intervention, j] <- beta
            }
        }
    }

    lambda_list <- lambda_list / p
    list(u_list = u_list, w_list = w_list, lambda_list = lambda_list)
}
