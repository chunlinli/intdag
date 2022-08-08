
library(ggplot2)
library(dplyr)
set.seed(1110)

df1 <- read.csv("simulation/learning_result/hub_iv_sufficient_TRUE_valid.csv")
df1$setup <- "Setup C"
df1$graph <- "Hub"
df1$method <- "Ours"

df2 <- read.csv("simulation/learning_result/hub_iv_sufficient_TRUE.csv")
df2$setup <- "Setup A"
df2$graph <- "Hub"
df2$method <- "Ours"

df3 <- read.csv("simulation/learning_result/hub_iv_sufficient_FALSE.csv")
df3$setup <- "Setup B"
df3$graph <- "Hub"
df3$method <- "Ours"

df4 <- read.csv("simulation/learning_result/2SPLS_hub_iv_sufficient_TRUE_valid.csv")
df4$setup <- "Setup C"
df4$graph <- "Hub"
df4$method <- "2SPLS"

df5 <- read.csv("simulation/learning_result/2SPLS_hub_iv_sufficient_TRUE.csv")
df5$setup <- "Setup A"
df5$graph <- "Hub"
df5$method <- "2SPLS"

df6 <- read.csv("simulation/learning_result/2SPLS_hub_iv_sufficient_FALSE.csv")
df6$setup <- "Setup B"
df6$graph <- "Hub"
df6$method <- "2SPLS"

data_hub <- rbind(df1, df2, df3, df4, df5, df6)

df1 <- read.csv("simulation/learning_result/random_iv_sufficient_TRUE_valid.csv")
df1$setup <- "Setup C"
df1$graph <- "Random"
df1$method <- "Ours"

df2 <- read.csv("simulation/learning_result/random_iv_sufficient_TRUE.csv")
df2$setup <- "Setup A"
df2$graph <- "Random"
df2$method <- "Ours"

df3 <- read.csv("simulation/learning_result/random_iv_sufficient_FALSE.csv")
df3$setup <- "Setup B"
df3$graph <- "Random"
df3$method <- "Ours"

df4 <- read.csv("simulation/learning_result/2SPLS_random_iv_sufficient_TRUE_valid.csv")
df4$setup <- "Setup C"
df4$graph <- "Random"
df4$method <- "2SPLS"

df5 <- read.csv("simulation/learning_result/2SPLS_random_iv_sufficient_TRUE.csv")
df5$setup <- "Setup A"
df5$graph <- "Random"
df5$method <- "2SPLS"

df6 <- read.csv("simulation/learning_result/2SPLS_random_iv_sufficient_FALSE.csv")
df6$setup <- "Setup B"
df6$graph <- "Random"
df6$method <- "2SPLS"

data_random <- rbind(df1, df2, df3, df4, df5, df6)

data_all <- rbind(data_hub, data_random)
data_all$shd <- data_all$shd + 0.05 * rchisq(nrow(data_all), df = 1) # jitter (positive for log transform)

fig_dir <- "simulation/fig"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
fig_name <- "structure_learning"

library(scales)
df <- data_all %>%
    group_by(n, graph, setup, method) %>%
    summarise(shd_avg = mean(shd, na.rm = TRUE), shd_se = sd(shd, na.rm = TRUE))

output <- ggplot(df, aes(n, shd_avg, group = interaction(setup, method), colour = method)) +
    facet_wrap(~graph) +
    geom_line(aes(linetype = setup)) +
    scale_y_continuous(trans = log2_trans(), breaks = c(1, 2, 4, 32, 512)) +
    labs(y = "structural Hamming distance", x = "sample size n", linetype = paste0("intervention\n", "setup")) +
    theme_bw() +
    ggtitle("p = 100, q = 500") + theme(plot.title = element_text(size=11))  

ggsave(sprintf("%s/%s.pdf", fig_dir, fig_name), output, width = 8, height = 4)
