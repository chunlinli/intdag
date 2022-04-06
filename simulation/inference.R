

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

inference_simulation <- function(p_seq = 100, q_seq = 500,
                                 n_seq = 200,
                                 n_sim = c(500, 50),
                                 signal_seq = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                                 graph_type = c("random", "hub"),
                                 iv_sufficient = FALSE, codimension_seq = c(1, 15),
                                 test_type = c("edge", "path"), rho = 0.5, seed = 1234) {
    set.seed(seed)

    library(mvtnorm)
    library(progress)
    source("intdag/intdag.R")

    res_dir <- "simulation/inference_result"
    if (!dir.exists(res_dir)) dir.create(res_dir)

    result_file <- paste0(
        test_type,
        "_inference_", graph_type,
        "_iv_sufficient_", iv_sufficient, ".csv"
    )
    if (file.exists(result_file)) {
        file.remove(result_file)
    }

    result <- c()

    if (test_type == "edge") {
        for (p in p_seq) {
            for (q in q_seq) {
                graph <- graph_generation(
                    p = p, q = q,
                    graph_type = graph_type, iv_sufficient = iv_sufficient
                )
                graph$u[, 20] <- 0

                for (signal in signal_seq) {
                    graph$u[1, 20] <- signal

                    for (n in n_seq) {
                        cat(paste0(
                            "p = ", p, ", q = ", q,
                            ", n = ", n,
                            ", type = ", graph_type,
                            ", iv = ", iv_sufficient,
                            ", test = ", test_type,
                            "\n"
                        ))

                        n_rep <- ifelse(signal == 0, n_sim[1], n_sim[2])

                        pb <- progress_bar$new(
                            format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                            total = n_rep,
                            complete = "=", # Completion bar character
                            incomplete = "-", # Incomplete bar character
                            current = ">", # Current bar character
                            clear = FALSE, # If TRUE, clears the bar when finish
                            width = 100
                        ) # Width of the progress bar



                        for (sim in 1:n_rep) {
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

                            for (codim in codimension_seq) {
                                f_mat <- matrix(0, p, p)
                                f_mat[1:codim, 20] <- 1

                                infer_out <- causal_inference(y, x, an_mat, in_mat, f_mat, method = "asymptotic", test_type = "edge")
                                result <- rbind(
                                    result,
                                    c(
                                        p, q, signal, n, sim, codim,
                                        infer_out$p_value, "DP-MLR"
                                    )
                                )
                                colnames(result) <- c("p", "q", "signal", "n", "sim", "codim", "p_value", "method")
                                write.csv(result, sprintf("%s/%s", res_dir, result_file), row.names = FALSE)

                                infer_out <- causal_inference(y, x, an_mat, in_mat, f_mat, method = "asymptotic", test_type = "edge")
                                result <- rbind(
                                    result,
                                    c(
                                        p, q, signal, n, sim, codim,
                                        infer_out$p_value, "LR"
                                    )
                                )
                                colnames(result) <- c("p", "q", "signal", "n", "sim", "codim", "p_value", "method")
                                write.csv(result, sprintf("%s/%s", res_dir, result_file), row.names = FALSE)

                                infer_out <- causal_inference(y, x, graph$u, graph$w, f_mat, method = "asymptotic", test_type = "edge")
                                result <- rbind(
                                    result,
                                    c(
                                        p, q, signal, n, sim, codim,
                                        infer_out$p_value, "OLR"
                                    )
                                )
                                colnames(result) <- c("p", "q", "signal", "n", "sim", "codim", "p_value", "method")
                                write.csv(result, sprintf("%s/%s", res_dir, result_file), row.names = FALSE)
                            }
                        }
                        cat("\n")
                    }
                }
            }
        }
    } else {
        # path test
        codim <- NA
        for (p in p_seq) {
            for (q in q_seq) {
                graph <- graph_generation(
                    p = p, q = q,
                    graph_type = graph_type, iv_sufficient = iv_sufficient
                )

                f_mat <- matrix(0, p, p)
                f_mat[1, 5] <- 1
                f_mat[5, 10] <- 1
                f_mat[10, 15] <- 1
                f_mat[15, 20] <- 1

                for (signal in signal_seq) {
                    graph$u[f_mat != 0] <- signal

                    for (n in n_seq) {
                        cat(paste0(
                            "p = ", p, ", q = ", q,
                            ", n = ", n,
                            ", type = ", graph_type,
                            ", iv = ", iv_sufficient,
                            ", test = ", test_type,
                            "\n"
                        ))

                        n_rep <- ifelse(signal == 0, n_sim[1], n_sim[2])

                        pb <- progress_bar$new(
                            format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                            total = n_rep,
                            complete = "=", # Completion bar character
                            incomplete = "-", # Incomplete bar character
                            current = ">", # Current bar character
                            clear = FALSE, # If TRUE, clears the bar when finish
                            width = 100
                        ) # Width of the progress bar

                        for (sim in 1:n_rep) {
                            pb$tick()

                            x <- matrix(rnorm(n * q), nrow = n, ncol = q)
                            y <- (x %*% graph$w + rmvnorm(n, sigma = diag(seq(from = 0.5, to = 1, length.out = p), p, p))) %*% solve(diag(p) - graph$u)

                            v_out <- v_estimation(y = y, x = x, model_selection = "bic")
                            top_out <- topological_order(v_out$v)
                            an_mat <- top_out$an_mat
                            in_mat <- top_out$in_mat

                            infer_out <- causal_inference(y, x, an_mat, in_mat, f_mat, method = "asymptotic", test_type = "path")
                            result <- rbind(
                                result,
                                c(
                                    p, q, signal, n, sim, codim, test_type,
                                    infer_out$p_value, "DP-MLR"
                                )
                            )
                            colnames(result) <- c("p", "q", "signal", "n", "sim", "codim", "test", "p_value", "method")
                            write.csv(result, sprintf("%s/%s", res_dir, result_file), row.names = FALSE)

                            infer_out <- causal_inference(y, x, an_mat, in_mat, f_mat, method = "asymptotic", test_type = "path")
                            result <- rbind(
                                result,
                                c(
                                    p, q, signal, n, sim, codim, test_type,
                                    infer_out$p_value, "LR"
                                )
                            )
                            colnames(result) <- c("p", "q", "signal", "n", "sim", "codim", "test", "p_value", "method")
                            write.csv(result, sprintf("%s/%s", res_dir, result_file), row.names = FALSE)

                            infer_out <- causal_inference(y, x, graph$u, graph$w, f_mat, method = "asymptotic", test_type = "path")
                            result <- rbind(
                                result,
                                c(
                                    p, q, signal, n, sim, codim, test_type,
                                    infer_out$p_value, "OLR"
                                )
                            )
                            colnames(result) <- c("p", "q", "signal", "n", "sim", "codim", "test", "p_value", "method")
                            write.csv(result, sprintf("%s/%s", res_dir, result_file), row.names = FALSE)
                        }
                        cat("\n")
                    }
                }
            }
        }
    }

    cat("Finished! \n")
}