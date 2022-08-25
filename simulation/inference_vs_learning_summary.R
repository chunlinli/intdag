library(ggplot2)
library(dplyr)

set.seed(1110)

data1 <- read.csv("simulation/additional_result/hub_iv_sufficient.csv")
data2 <- read.csv("simulation/additional_result/random_iv_sufficient.csv")
data <- rbind(data1, data2)

df_lambda <- data %>% summarise(
  lambda1 = mean(lambda1),
  lambda2 = mean(lambda2),
  lambda3 = mean(lambda3),
  lambda4 = mean(lambda4),
  lambda5 = mean(lambda5)
)

df1_1 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select1_2, na.rm = TRUE))
df1_1$method <- paste0(
  "selection: ",
  formatC(df_lambda[[1]], format = "f", digits = 3)
)

df1_2 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select1_2, na.rm = TRUE))
df1_2$method <- paste0(
  "selection: ",
  formatC(df_lambda[[2]], format = "f", digits = 3)
)

df1_3 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select1_3, na.rm = TRUE))
df1_3$method <- paste0(
  "selection: ",
  formatC(df_lambda[[3]], format = "f", digits = 3)
)

df1_4 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select1_4, na.rm = TRUE))
df1_4$method <- paste0(
  "selection: ",
  formatC(df_lambda[[4]], format = "f", digits = 3)
)

df1_5 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select1_5, na.rm = TRUE))
df1_5$method <- paste0(
  "selection: ",
  formatC(df_lambda[[5]], format = "f", digits = 3)
)

df1_select <- rbind(df1_1, df1_2, df1_3, df1_4, df1_5)

df1_test <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(pvalue1 < 0.05, na.rm = TRUE))
df1_test$method <- "testing"

df1 <- rbind(df1_select, df1_test)

fig_dir <- "simulation/fig"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
fig_name <- "additional"

print(df1, n = 160)

output <- ggplot(df1, aes(n, res, group = interaction(graph, method), colour = method)) +
  facet_wrap(~hypothesis) +
  geom_point(aes(shape = graph), size = 2.5) +
  geom_line() +
  labs(y = "power", x = "sample size n", shape = paste0("graph type")) +
  theme_bw() +
  ggtitle("p = 50, q = 100") +
  theme(plot.title = element_text(size = 11))

output

ggsave(sprintf("%s/%s_1.pdf", fig_dir, fig_name), output, width = 8, height = 5)





df15_1 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select15_2, na.rm = TRUE))
df15_1$method <- paste0(
  "kappa = ",
  formatC(1, format = "d")
)

df15_2 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select15_2, na.rm = TRUE))
df15_2$method <- paste0(
  "kappa = ",
  formatC(2, format = "d")
)

df15_3 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select15_3, na.rm = TRUE))
df15_3$method <- paste0(
  "kappa = ",
  formatC(3, format = "d")
)

df15_4 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select15_4, na.rm = TRUE))
df15_4$method <- paste0(
  "kappa = ",
  formatC(4, format = "d")
)

df15_5 <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(select15_5, na.rm = TRUE))
df15_5$method <- paste0(
  "kappa = ",
  formatC(5, format = "d")
)

df15_select <- rbind(df15_1, df15_2, df15_3, df15_4, df15_5)

df15_test <- data %>%
  group_by(n, graph, hypothesis) %>%
  summarise(res = mean(pvalue15 < 0.05, na.rm = TRUE))
df15_test$method <- "alpha = 0.05"

df15 <- rbind(df15_select, df15_test)

fig_dir <- "simulation/fig"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
fig_name <- "additional"

print(df15, n = 160)

output <- ggplot(df15, aes(n, res, group = interaction(graph, method), colour = method)) +
  facet_wrap(~hypothesis) +
  geom_point(aes(shape = graph), size = 2.5) +
  geom_line() +
  labs(y = "rejection rate", x = "sample size n", shape = paste0("graph type")) +
  theme_bw() +
  ggtitle("p = 50, q = 100") +
  theme(plot.title = element_text(size = 11))

output

ggsave(sprintf("%s/%s_15.pdf", fig_dir, fig_name), output, width = 10, height = 5)

