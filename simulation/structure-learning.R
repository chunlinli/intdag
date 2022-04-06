

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

structure_learning_simulation <- function(p_seq = 100, q_seq = 500,
                                          n_seq = c(200, 250, 300, 350, 400),
                                          n_sim = 100,
                                          graph_type = c("random", "hub"),
                                          iv_sufficient = FALSE, rho = 0.5, seed = 1234) {
    set.seed(seed)

    library(mvtnorm)
    library(progress)
    source("intdag/intdag.R")

    result_file <- paste0(
        "./", graph_type,
        "_iv_sufficient_", iv_sufficient, ".csv"
    )
    if (file.exists(result_file)) {
        file.remove(result_file)
    }

    result <- c()

    for (p in p_seq) {
        for (q in q_seq) {
            graph <- graph_generation(
                p = p, q = q,
                graph_type = graph_type, iv_sufficient = iv_sufficient
            )

            for (n in n_seq) {
                cat(paste0(
                    "p = ", p, ", q = ",
                    q, ", n = ", n, ", type = ", graph_type, "\n"
                ))

                pb <- progress_bar$new(
                    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                    total = n_sim,
                    complete = "=", # Completion bar character
                    incomplete = "-", # Incomplete bar character
                    current = ">", # Current bar character
                    clear = FALSE, # If TRUE, clears the bar when finish
                    width = 100
                ) # Width of the progress bar

                for (sim in 1:n_sim) {
                    pb$tick()

                    x <- matrix(rnorm(n * q), nrow = n, ncol = q)
                    if (rho != 0) {
                        for (j in 2:q) {
                            x[, j] <- sqrt(1 - rho^2) * x[, j] + rho * x[, j - 1]
                        }
                    }
                    y <- (x %*% graph$w + rmvnorm(n, sigma = diag(seq(from = 0.5, to = 1, length.out = p), p, p))) %*% solve(diag(p) - graph$u)

                    v_out <- v_estimation(y = y, x = x, model_selection = "bic")
                    top_out <- topological_order(v_out$v)
                    an_mat <- top_out$an_mat
                    in_mat <- top_out$in_mat
                    discovery_out <- causal_discovery(y = y, x = x, an_mat = an_mat, in_mat = in_mat)

                    result <- rbind(
                        result,
                        c(
                            p, q, n, sim,
                            sum(abs((abs(discovery_out$u) > 0.15) - (graph$u != 0)))
                        )
                    )
                    colnames(result) <- c("p", "q", "n", "sim", "shd")
                    write.csv(result, result_file, row.names = FALSE)
                }
                cat("\n")
            }
        }
    }

    cat("Finished! \n")
}

structure_learning_simulation(graph_type = "random", iv_sufficient = FALSE)
structure_learning_simulation(graph_type = "random", iv_sufficient = TRUE)
structure_learning_simulation(graph_type = "hub", iv_sufficient = TRUE)
structure_learning_simulation(graph_type = "hub", iv_sufficient = FALSE)