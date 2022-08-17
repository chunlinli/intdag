
library(ggplot2)
library(dplyr)

set.seed(1110)

df1 <- read.csv("simulation/additional_result/hub_iv_sufficient.csv")
df2 <- read.csv("simulation/additional_result/random_iv_sufficient.csv")
data1 <- rbind(df1, df2)
data2 <- data1
data1$method <- "structure\n learning"
data1$res1 <- data1$select1
data1$res15 <- data1$select15
data2$method <- "inference"
data2$res1 <- ifelse(data2$pvalue1 < 0.05, 1, 0)
data2$res15 <- ifelse(data2$pvalue15 < 0.05, 1, 0)
data <- rbind(data1, data2)

#data$shd <- data_all$shd + 0.05 * rchisq(nrow(data_all), df = 1) # jitter (positive for log transform)

fig_dir <- "simulation/fig"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
fig_name <- "additional"

library(scales)
df <- data %>%
  group_by(n, graph,method,hypothesis) %>%
  summarise(res = mean(res1, na.rm = TRUE))



output <- ggplot(df, aes(n, res, group=interaction(graph,method), colour = method)) +
  facet_wrap(~hypothesis) +
  geom_point(aes(shape = graph), size = 2.5) +
  geom_line() +
  #scale_y_continuous(trans = log2_trans(), breaks = c(1, 2, 4, 32, 512)) +
  labs(y = "power", x = "sample size n", shape = paste0("graph type")) +
  theme_bw() +
  ggtitle("p = 50, q = 100") + theme(plot.title = element_text(size=11))

output

ggsave(sprintf("%s/%s.pdf", fig_dir, fig_name), output, width = 8, height = 5)



df3 <- data %>%
  group_by(n, graph, magnitude, method, hypotheis) %>%
  summarise(shd_mean = mean(shd_nuisance, na.rm = TRUE))



print(df3, n = 130)


# fdr
df1 <- data1 %>%
  group_by(n, graph,magnitude) %>%
  summarise(fdr_mean = mean(fdr, na.rm = TRUE))
print(df1, n=80)



