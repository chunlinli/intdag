source("./simulation/inference.R")

inference_simulation(graph_type = "random", iv_sufficient = TRUE, test_type = "edge", n_seq = c(200,300), seed = 1110)
inference_simulation(graph_type = "random", iv_sufficient = FALSE, test_type = "edge", n_seq = c(200,300), seed = 1110)
inference_simulation(graph_type = "hub", iv_sufficient = TRUE, test_type = "edge", n_seq = c(200,300), seed = 1110)
inference_simulation(graph_type = "hub", iv_sufficient = FALSE, test_type = "edge", n_seq = c(200,300), seed = 1110)