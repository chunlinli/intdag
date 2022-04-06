
library(ggplot2)
library(dplyr)

df1 <- read.csv("simulation/inference_result/path_inference_hub_iv_sufficient_TRUE.csv")
df1$graph <- "Hub"
df1$setup <- "Setup A"

df2 <- read.csv("simulation/inference_result/path_inference_hub_iv_sufficient_FALSE.csv")
df2$graph <- "Hub"
df2$setup <- "Setup B"

df3 <- read.csv("simulation/inference_result/path_inference_random_iv_sufficient_TRUE.csv")
df3$graph <- "Random"
df3$setup <- "Setup A"

df4 <- read.csv("simulation/inference_result/path_inference_random_iv_sufficient_FALSE.csv")
df4$graph <- "Random"
df4$setup <- "Setup B"

data_all <- rbind(df1, df2, df3, df4)

fig_dir <- "simulation/fig"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
fig_name <- "path_inference"

df <- data_all %>%
    group_by(n, graph, setup, signal, method) %>%
    summarise(reject_rate = mean((p_value < 0.05), na.rm = TRUE))

df$graph_setup <- interaction(df$graph, df$setup)
levels(df$graph_setup) <- c("Hub, Setup A", "Random, Setup A", "Hub, Setup B", "Random, Setup B")

idx0 <- which(df$method == "DP-MLR")
df$method[idx0] <- "DP-LR (Ours)"
idx0 <- which(df$method == "LR")
df$method[idx0] <- "LR (Ours)"
idx0 <- which(df$method == "OLR")
df$method[idx0] <- "OLR (Ideal)"

library(scales)
trans_p_value <- trans_new(
  name = "logp",
  transform = function(x) log(x*20),
  inverse = function(x) exp(x)/20,
  breaks = log_breaks()
)

output <- ggplot(df, aes(signal, reject_rate, colour = method)) +
    facet_grid(~ graph_setup) +
    geom_line(aes(linetype = as.factor(n)), size = 0.5) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", size = 0.3) +
    ylim(c(0,1)) + 
    labs(
        y = "rejection rate",
        x = expression(paste("coefficients ", U[1 ~ 5], ", ", U[5 ~ 10], ", ", U[10 ~ 15], ", ", U[15 ~ 20])),
        linetype = "sample size n"
    ) +
    theme_bw() +
    ggtitle("p = 100, q = 500") +
    theme(legend.position = "bottom", plot.title = element_text(size=11)) 
    
ggsave(sprintf("%s/%s.pdf", fig_dir, fig_name), output, width = 8, height = 3.22)