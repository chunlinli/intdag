
v_estimation_internal <- function(y, x, v_init, lambda) {
    p <- ncol(y)
    q <- ncol(x)
    v <- matrix(0, nrow = q, ncol = p)
    
    for (j in 1:p) {
        m <- ncvfit(X = x, y = y[, j], init = v_init[, j], lambda = lambda[j])
        v[, j] <- as.numeric(m$beta[-1])
    }
    list(v = v)
}

causal_inference_dp <- function(y, x, an_mat, in_mat, f_mat, v_out, mc_size = 500) {
    p <- ncol(y)
    q <- ncol(x)
    if (q < p) stop("No sufficient interventions: q < p.")

    if (nrow(y) == nrow(x)) {
        n <- nrow(y)
    } else {
        stop("Dimensions of y and x do not match.")
    }

    d_mat <- f_mat
    hypothesized_nodes <- which(colSums(d_mat) != 0)
    lr_seq <- rep(0, length(hypothesized_nodes))

    for (j in hypothesized_nodes) {
        test <- which(d_mat[, j] != 0)
        ancestor <- setdiff(which(an_mat[, j] != 0), test)
        intervention <- which(in_mat[, j] != 0)

        if (n > length(test) + length(ancestor) + length(intervention) + 1) {
            if (length(ancestor) > 0) {
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
            lr_seq[j] <- (sum(res_h0^2) - sum(res_h1^2)) / sig2
        }
    }

    likelihood_ratio <- sum(lr_seq)
    if (likelihood_ratio == 0) {
        list(likelihood_ratio = 0, p_value = 1, mc_size = mc_size)
    }

    lr_dp <- rep(0, mc_size)

    # need parallel dp
    for (mc in 1:mc_size) {

        # generate dp sample
        err_dp <- matrix(rnorm(n*p, sd = sqrt(v_out$err_var)), nrow=n, byrow = TRUE) 
        y_dp <- y + err_dp

        v_out_dp <- v_estimation_internal(y = y_dp, x = x, v_init = v_out$v, lambda = v_out$lambda)
        top_out_dp <- topological_order(v_out_dp$v)
        an_mat_dp <- top_out_dp$an_mat
        in_mat_dp <- top_out$in_mat

        for (j in hypothesized_nodes) {
            test <- which(d_mat[, j] != 0)
            A1 <- union(test, which(an_mat[, j] != 0))
            A2 <- union(test, which(an_mat_dp[, j] != 0))
            B1 <- which(in_mat[, j] != 0)
            B2 <- which(in_mat_dp[, j] != 0)
            under_select_nodes <- setdiff(A1, A2)
            over_select_nodes <- setdiff(A2, A1)
            under_select_intervention <- setdiff(B1, B2)

            ancestor <- setdiff(which(an_mat_dp[, j] != 0), test)
            intervention <- which(in_mat_dp[, j] != 0)

            if (length(under_select_nodes) == 0 && length(under_select_intervention) == 0 && sum(an_mat[j, over_select_nodes]) == 0 && n > length(test) + length(ancestor) + length(intervention) + 1) {
                if (length(ancestor) > 0) {
                    z_h1 <- cbind(y[, test], y[, ancestor], x[, intervention])
                    z_h0 <- cbind(y[, ancestor], x[, intervention])
                } else {
                    z_h1 <- cbind(y[, test], x[, intervention])
                    z_h0 <- x[, intervention]
                }
                P1 <- z_h1 %*% tcrossprod(solve(crossprod(z_h1)), z_h1)
                P0 <- z_h0 %*% tcrossprod(solve(crossprod(z_h0)), z_h0)
                eP1e <- crossprod(err_dp[, j], P1 %*% err_dp[, j])
                eP0e <- crossprod(err_dp[, j], P0 %*% err_dp[, j])
                lr_dp[mc] <- lr_dp[mc] + (eP1e - eP0e) / ((sum(err_dp[, j]^2) - eP1e) / (n - length(z_h1)))
            } else {
                lr_dp[mc] <- lr_dp[mc] + lr_seq[j] #rf(1, df1=length(test), df2=n-length(z_h1))
            }
        }
    }
    p_value <- sum(lr_dp >= likelihood_ratio) / sum(lr_dp > 0)
    list(likelihood_ratio = likelihood_ratio, p_value = p_value, mc_size = mc_size)
}