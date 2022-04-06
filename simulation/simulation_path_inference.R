source("./simulation/inference.R")

inference_simulation(graph_type = "random", iv_sufficient = TRUE, test_type = "path", n_sim = c(200,50), n_seq = c(200,300))
inference_simulation(graph_type = "random", iv_sufficient = FALSE, test_type = "path", n_sim = c(200,50), n_seq = c(200,300))
inference_simulation(graph_type = "hub", iv_sufficient = TRUE, test_type = "path", n_sim = c(200,50), n_seq = c(200,300))
inference_simulation(graph_type = "hub", iv_sufficient = FALSE, test_type = "path", n_sim = c(200,50), n_seq = c(200,300))