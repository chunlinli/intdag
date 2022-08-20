library(mvtnorm)

graph_generation <- function(p = 100, q = 500,
                             graph_type = c("random", "hub"),
                             iv_sufficient = FALSE) {
    if (graph_type == "random") {
        sparsity <- 1 / p
        u <- matrix(rbinom(p * p, 1, sparsity), p, p)
        u[lower.tri(u, diag = TRUE)] <- 0

        w <- matrix(0, nrow = q, ncol = p)
        for (k in 1:(p - 1)) {
            w[k, k] <- 1
            if (!iv_sufficient) w[k, k + 1] <- 1

            w[k + p, k] <- 1
            w[k + p, k + 1] <- 1
        }
        w[p, p] <- 1
        w[2 * p, p] <- 1
    } else {
        # generate graph
        num_of_hub <- 2
        idx <- rep(1:num_of_hub, length.out = p)
        u <- matrix(0, p, p)
        for (k in 1:num_of_hub) {
            u[k, idx == k] <- 1
        }
        diag(u) <- 0

        w <- matrix(0, nrow = q, ncol = p)
        for (k in 1:(p - 1)) {
            w[k, k] <- 1
            if (!iv_sufficient) w[k, k + 1] <- 1

            w[k + p, k] <- 1
            w[k + p, k + 1] <- 1
        }
        w[p, p] <- 1
        w[2 * p, p] <- 1
    }

    list(u = u, w = w)
}


