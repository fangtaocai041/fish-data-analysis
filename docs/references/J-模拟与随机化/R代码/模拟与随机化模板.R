# ============================================================
# J — 模拟与随机化 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(vegan)
library(boot)         # 自举
library(ggplot2)

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(42)
n <- 50

# 一个有空间结构的数据集
x <- runif(n, 0, 100)
y <- runif(n, 0, 100)
group <- factor(rep(c("A", "B"), each = n/2))

# 响应变量（处理组效应 + 随机误差）
response <- 10 + ifelse(group == "A", 0, 3) + rnorm(n, 0, 4)

df <- data.frame(response, group, x, y)

# ============================================================
# J1 — 零模型
# ============================================================

# --- 1a. 物种共现矩阵的零模型 ---
set.seed(123)
spe_mat <- matrix(
  rbinom(200, 1, 0.3),
  nrow = 10,
  dimnames = list(paste0("Site", 1:10), paste0("sp", 1:20))
)

# 使用 quasiswap 保留行列和
null_mods <- nullmodel(spe_mat, "quasiswap")
null_sims <- simulate(null_mods, nsim = 99)

# 计算观测的 C-score (共现指数)
obs_cscore <- apply(spe_mat, 2, function(x) {
  apply(spe_mat, 2, function(y) {
    (colSums(matrix(x)) - 1) * (colSums(matrix(y)) - 1) / (nrow(spe_mat) * (nrow(spe_mat) - 1) / 2)
  })
})
# 简化：用 vegan 的 oecosimu
oecosimu(spe_mat, nestedchecker, "quasiswap", nsimul = 999)

# --- 1b. C-score 检验 ---
cs_out <- oecosimu(spe_mat, function(x) {
  com <- t(x)
  n <- ncol(com)
  cs <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cs <- cs + sum(com[, i] == 1 & com[, j] == 1) * 
                   sum(com[, i] == 1 | com[, j] == 1)
    }
  }
  return(cs)
}, "quasiswap", nsimul = 999)
cs_out  # p-value 判断是否显著共现/排斥

# ============================================================
# J2 — 置换检验
# ============================================================

# --- 2a. 手动置换检验 ---
# 检验两组间均值差异
obs_diff <- mean(df$response[df$group == "A"]) -
            mean(df$response[df$group == "B"])

n_perm <- 9999
perm_diffs <- numeric(n_perm)
for (i in 1:n_perm) {
  perm_group <- sample(df$group)
  perm_diffs[i] <- mean(df$response[perm_group == "A"]) -
                   mean(df$response[perm_group == "B"])
}

# 双侧 p 值
p_value <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
p_value

# 可视化置换分布
ggplot(data.frame(diff = perm_diffs), aes(x = diff)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = obs_diff, color = "red", size = 1.5) +
  geom_vline(xintercept = -obs_diff, color = "red", size = 1.5, linetype = 2) +
  labs(title = paste0("置换检验 (p = ", round(p_value, 4), ")"),
       x = "组间均值差（置换分布）", y = "频数") +
  theme_bw()

# ============================================================
# J3 — 自举 (Bootstrap)
# ============================================================

# --- 3a. 单变量自举 ---
boot_mean <- function(data, indices) {
  mean(data[indices])
}

boot_obj <- boot(df$response, statistic = boot_mean, R = 2000)
boot.ci(boot_obj, conf = 0.95, type = c("basic", "perc", "bca"))

# --- 3b. 回归模型自举 ---
boot_reg <- function(data, indices) {
  d <- data[indices, ]
  fit <- lm(response ~ group, data = d)
  return(coef(fit))
}

boot_reg_obj <- boot(df, statistic = boot_reg, R = 2000)
boot_reg_obj
boot.ci(boot_reg_obj, index = 2, conf = 0.95)  # groupB 系数的 CI

# --- 3c. 可视化 bootstrap 分布 ---
boot_df <- data.frame(
  intercept = boot_reg_obj$t[, 1],
  groupB_effect = boot_reg_obj$t[, 2]
)
ggplot(boot_df, aes(x = groupB_effect)) +
  geom_histogram(bins = 30, fill = "coral", alpha = 0.7) +
  geom_vline(xintercept = boot_reg_obj$t0[2], color = "darkred", size = 1.5) +
  labs(title = "Bootstrap 分布：处理效应", x = "系数估计", y = "频数") +
  theme_bw()

# ============================================================
# J4 — 贝叶斯 (使用 brms)
# ============================================================

# 安装: install.packages("brms")
# brms 依赖于 rstan，需要较长时间安装

# --- 4a. 简单贝叶斯线性模型 ---
library(brms)

# brm_mod <- brm(
#   response ~ group,
#   data = df,
#   prior = set_prior("normal(0, 10)", class = "b"),
#   chains = 4, iter = 2000, warmup = 1000,
#   seed = 42
# )
# summary(brm_mod)
# plot(brm_mod)
# posterior_interval(brm_mod)

# 手动贝叶斯（简化版）
# 用 MCMC 估计均值差
n_iter <- 10000
mu_a <- 0; sigma_a <- 10  # 先验
mcmc_means <- matrix(NA, n_iter, 2)
for (i in 1:n_iter) {
  # 从后验采样（共轭先验）
  y_a <- df$response[df$group == "A"]
  y_b <- df$response[df$group == "B"]
  n_a <- length(y_a); n_b <- length(y_b)
  
  # 假设方差已知（用样本方差）
  s2 <- var(df$response)
  
  # 均值的后验（正态-正态共轭）
  post_var_a <- 1 / (1/sigma_a^2 + n_a/s2)
  post_mean_a <- post_var_a * (mu_a/sigma_a^2 + sum(y_a)/s2)
  mcmc_means[i, 1] <- rnorm(1, post_mean_a, sqrt(post_var_a))
  
  post_var_b <- 1 / (1/sigma_a^2 + n_b/s2)
  post_mean_b <- post_var_b * (mu_a/sigma_a^2 + sum(y_b)/s2)
  mcmc_means[i, 2] <- rnorm(1, post_mean_b, sqrt(post_var_b))
}

# 后验分布
diff_post <- mcmc_means[, 2] - mcmc_means[, 1]
quantile(diff_post, c(0.025, 0.5, 0.975))
mean(diff_post > 0)  # 后验概率

# ============================================================
# J5 — 简单 IBM（基于个体的模拟）
# ============================================================

# 模拟一个简单的种群动态
simulate_population <- function(n0 = 50, r = 0.1, K = 100,
                                 years = 50, env_noise = 0.05) {
  n <- numeric(years)
  n[1] <- n0
  for (t in 2:years) {
    # 逻辑斯蒂增长 + 环境噪声
    expected <- n[t-1] * exp(r * (1 - n[t-1]/K))
    n[t] <- rlnorm(1, meanlog = log(expected), sdlog = env_noise)
  }
  return(n)
}

# 多个模拟
sims <- replicate(100, simulate_population())
sim_df <- as.data.frame(sims)
sim_df$year <- 1:nrow(sim_df)
sim_long <- pivot_longer(sim_df, -year, names_to = "sim", values_to = "size")

ggplot(sim_long, aes(x = year, y = size, group = sim)) +
  geom_line(alpha = 0.1, color = "steelblue") +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               color = "red", size = 2) +
  labs(title = "100 次逻辑斯蒂种群模拟",
       subtitle = "r = 0.1, K = 100, σ_env = 0.05",
       y = "种群大小") +
  theme_bw()
