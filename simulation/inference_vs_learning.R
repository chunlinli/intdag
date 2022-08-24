
# This simulation compares the performance of inference and learning algorithms
# repeat 100 times
# p = 30, q = 100, n = 100, 200, 400, 800
# focused coefs: 1/sqrt(n)
# small nuisance coefs: 1/sqrt(n)
# large nuisance coefs: 0.4, 0.7, 1
# compare: type I error (false disovery), power
# intervention setup: A
# graph types: random, hub
comp_simulation <- function(p_seq = c(50), q_seq = c(100),
                            n_seq = c(100, 200, 300, 400),
                            n_sim = 300,
                            magnitude_seq = 1,
                            hypothesis = c("Null", "Alternative"),
                            graph_type = c("random", "hub"),
                            iv_sufficient = TRUE, rho = .5,
                            seed = 1234) {
    set.seed(seed)

    library(mvtnorm)
    library(progress)
    library(ncvreg)

    source("intdag/R/v_estimation.R")
    source("intdag/R/topological_order.R")
    source("intdag/R/causal_discovery_internal.R")
    source("intdag/R/causal_inference.R")

    res_dir <- "simulation/additional_result"
    if (!dir.exists(res_dir)) dir.create(res_dir)

    result_file <- paste0(
        graph_type,
        "_iv_sufficient.csv"
    )
    if (file.exists(result_file)) {
        file.remove(result_file)
    }

    result <- c()

    for (p in p_seq) {
        for (q in q_seq) {
            for (magnitude in magnitude_seq) {
                graph <- graph_generation(
                    p = p, q = q,
                    magnitude_large = magnitude,
                    # magnitude_small = 0,
                    graph_type = graph_type, iv_sufficient = iv_sufficient
                )

                graph$u[1, 20] <- NA

                for (n in n_seq) {
                    for (h in hypothesis) {
                        graph$u[1:15, 20] <- 0
                        graph$u[1, 20] <- 2 * sqrt(1 / n) * (h == "Alternative")

                        cat(paste0(
                            "p = ", p, ", q = ", q,
                            ", n = ", n, ", hypothesis = ", h,
                            ", magnitude = ", magnitude,
                            ", type = ", graph_type, "\n"
                        ))

                        pb <- progress_bar$new(
                            format = paste0(
                                "(:spin) [:bar] :percent ",
                                "[Elapsed time: :elapsedfull || ",
                                "Estimated time remaining: :eta]"
                            ),
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
                                    x[, j] <- sqrt(1 - rho^2) * x[, j] +
                                        rho * x[, j - 1]
                                }
                            }
                            y <- (x %*% graph$w +
                                rmvnorm(n,
                                    sigma = diag(
                                        seq(from = 0.5, to = 1, length.out = p),
                                        p, p
                                    )
                                )) %*% solve(diag(p) - graph$u)

                            # save(y, x, file = "data.RData")

                            v_out <- v_estimation(
                                y = y, x = x,
                                model_selection = "bic"
                            )

                            top_out <- topological_order(v_out$v, thresh = .15)
                            an_mat <- top_out$an_mat
                            in_mat <- top_out$in_mat
                            discovery_out <- causal_discovery_internal(
                                y = y, x = x,
                                an_mat = an_mat,
                                in_mat = in_mat,
                                penalty = "MCP"
                            )

                            f_mat <- matrix(0, nrow = p, ncol = p)
                            f_mat[1, 20] <- 1
                            infer_out1 <- causal_inference(y, x,
                                an_mat, in_mat, f_mat,
                                method = "asymptotic", test_type = "edge"
                            )

                            f_mat[1:15, 20] <- 1
                            infer_out15 <- causal_inference(y, x,
                                an_mat, in_mat, f_mat,
                                method = "asymptotic", test_type = "edge"
                            )
                            result <- rbind(
                                result,
                                c(
                                    p, q, n, h, sim, magnitude, graph_type,
                                    (discovery_out$u_list[[1]][1, 20] != 0) * 1,
                                    (discovery_out$u_list[[2]][1, 20] != 0) * 1,
                                    (discovery_out$u_list[[3]][1, 20] != 0) * 1,
                                    (discovery_out$u_list[[4]][1, 20] != 0) * 1,
                                    (discovery_out$u_list[[5]][1, 20] != 0) * 1,
                                    any(discovery_out$u_list[[1]][1:15, 20] != 0) * 1,
                                    any(discovery_out$u_list[[2]][1:15, 20] != 0) * 1,
                                    any(discovery_out$u_list[[3]][1:15, 20] != 0) * 1,
                                    any(discovery_out$u_list[[4]][1:15, 20] != 0) * 1,
                                    any(discovery_out$u_list[[5]][1:15, 20] != 0) * 1,
                                    discovery_out$lambda_list[1],
                                    discovery_out$lambda_list[2],
                                    discovery_out$lambda_list[3],
                                    discovery_out$lambda_list[4],
                                    discovery_out$lambda_list[5],
                                    infer_out1$p_value,
                                    infer_out15$p_value
                                )
                            )

                            colnames(result) <- c(
                                "p", "q", "n", "hypothesis", "sim", "magnitude", "graph",
                                "select1_1",
                                "select1_2",
                                "select1_3",
                                "select1_4",
                                "select1_5",
                                "select15_1",
                                "select15_2",
                                "select15_3",
                                "select15_4",
                                "select15_5",
                                "lambda1",
                                "lambda2",
                                "lambda3",
                                "lambda4",
                                "lambda5",
                                "pvalue1",
                                "pvalue15"
                            )
                            write.csv(result,
                                sprintf("%s/%s", res_dir, result_file),
                                row.names = FALSE
                            )
                        }
                        cat("\n")
                    }
                }
            }
        }
    }

    cat("Finished! \n")
}


graph_generation <- function(p = 30, q = 100,
                             magnitude_large = .7,
                             graph_type = c("random", "hub"),
                             iv_sufficient = TRUE) {
    if (graph_type == "random") {
        sparsity <- 2 / p
        u <- matrix(rbinom(p * p, 1, sparsity) * magnitude_large, p, p)

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
    }

    list(u = u, w = w)
}


comp_simulation(graph_type = "random", seed = 1110)
comp_simulation(graph_type = "hub", seed = 1110)
