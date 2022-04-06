
causal_inference <- function(y, x, an_mat, in_mat, f_mat, method = c("asymptotic", "dp"), test_type = c("edge", "test")) {
    # based on an_mat

    p <- ncol(y)
    q <- ncol(x)
    if (q < p) stop("No sufficient interventions: q < p.")

    if (nrow(y) == nrow(x)) {
        n <- nrow(y)
    } else {
        stop("Dimensions of y and x do not match.")
    }

    if (test_type == "edge") {
        if (method == "asymptotic") {
            # d_mat should be acyclic
            d_mat <- f_mat
            likelihood_ratio <- 0

            for (j in 1:p) {
                test <- which(d_mat[, j] != 0)
                ancestor <- setdiff(which(an_mat[, j] != 0), test)
                intervention <- which(in_mat[, j] != 0)

                if (length(test) > 0 && n > length(test) + length(ancestor) + length(intervention) + 1) {
                    if (length(ancestor) > 0) {

                        # add sparsity here


                        z_h1 <- cbind(y[, test], y[, ancestor], x[, intervention])
                        beta_h1 <- solve(crossprod(z_h1), crossprod(z_h1, y[, j]))
                        res_h1 <- y[, j] - z_h1 %*% beta_h1




                        z_h0 <- cbind(y[, ancestor], x[, intervention])
                        beta_h0 <- solve(crossprod(z_h0), crossprod(z_h0, y[, j]))
                        res_h0 <- y[, j] - z_h0 %*% beta_h0
                    } else {
                        z_h1 <- cbind(y[, test], x[, intervention])
                        beta_h1 <- solve(crossprod(z_h1), crossprod(z_h1, y[, j]))
                        res_h1 <- y[, j] - z_h1 %*% beta_h1

                        z_h0 <- x[, intervention]
                        beta_h0 <- solve(crossprod(z_h0), crossprod(z_h0, y[, j]))
                        res_h0 <- y[, j] - z_h0 %*% beta_h0
                    }
                    sig2 <- sum(res_h1^2) / (n - ncol(z_h1))
                    likelihood_ratio <- likelihood_ratio + (sum(res_h0^2) - sum(res_h1^2)) / sig2
                }
            }

            p_value <- pchisq(likelihood_ratio, df = sum(d_mat != 0), lower.tail = FALSE)

            list(likelihood_ratio = likelihood_ratio, df = sum(d_mat != 0), p_value = p_value)
        } else {
            # noise
            # refit






            # compute DP likelihood ratio
        }
    } else {


        # test type is path
        if (method == "asymptotic") {
            # d_mat should be acyclic
            d_mat <- f_mat
            path_length <- sum(d_mat)
            likelihood_ratios <- rep(0, path_length)
            k <- 1

            for (j in 1:p) {
                test <- which(d_mat[, j] != 0)
                ancestor <- setdiff(which(an_mat[, j] != 0), test)
                intervention <- which(in_mat[, j] != 0)

                if (length(test) > 0 && n > length(test) + length(ancestor) + length(intervention) + 1) {
                    if (length(ancestor) > 0) {

                        # add sparsity here
                        z_h1 <- cbind(y[, test], y[, ancestor], x[, intervention])
                        beta_h1 <- solve(crossprod(z_h1), crossprod(z_h1, y[, j]))
                        res_h1 <- y[, j] - z_h1 %*% beta_h1

                        z_h0 <- cbind(y[, ancestor], x[, intervention])
                        beta_h0 <- solve(crossprod(z_h0), crossprod(z_h0, y[, j]))
                        res_h0 <- y[, j] - z_h0 %*% beta_h0
                    } else {
                        z_h1 <- cbind(y[, test], x[, intervention])
                        beta_h1 <- solve(crossprod(z_h1), crossprod(z_h1, y[, j]))
                        res_h1 <- y[, j] - z_h1 %*% beta_h1

                        z_h0 <- x[, intervention]
                        beta_h0 <- solve(crossprod(z_h0), crossprod(z_h0, y[, j]))
                        res_h0 <- y[, j] - z_h0 %*% beta_h0
                    }
                    sig2 <- sum(res_h1^2) / (n - ncol(z_h1))
                    likelihood_ratios[k] <- (sum(res_h0^2) - sum(res_h1^2)) / sig2
                    k <- k + 1
                }
            }

            p_value <- max(pchisq(likelihood_ratios, df = 1, lower.tail = FALSE))

            list(likelihood_ratios = likelihood_ratios, p_value = p_value)
        } else {
            # noise
            # refit
            # compute DP likelihood ratio
        }
    }
}
