# ============================================================
# A — 数据基础：数据探索与预处理 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

# ---- 加载包 ----
library(tidyverse)   # 数据操作 + 可视化
library(vegan)       # 生态学数据分析（含标准化函数）
library(corrplot)    # 相关矩阵可视化
library(car)         # Box-Cox 转换
library(missForest)  # 随机森林插补（可选）
library(VIM)         # 缺失值可视化
library(psych)       # 描述统计

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(42)
n <- 100
fake_data <- data.frame(
  site = paste0("S", 1:n),
  pH = rnorm(n, mean=7.0, sd=0.8),
  DO = rnorm(n, mean=8.5, sd=2.0),
  temp = rnorm(n, mean=22, sd=4),
  NO3 = rlnorm(n, meanlog=1, sdlog=1.2),   # 偏态分布
  species_richness = rpois(n, lambda=15)    # 计数数据
)
# 手动加一些缺失值
fake_data$pH[sample(1:n, 5)] <- NA
fake_data$NO3[sample(1:n, 8)] <- NA
# 手动加异常值
fake_data$DO[3] <- 25
fake_data$temp[50] <- -5

# ============================================================
# A1 — 异常值检测
# ============================================================

# --- 1a. Boxplot 规则 ---
detect_outliers_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower <- q[1] - 1.5 * iqr
  upper <- q[2] + 1.5 * iqr
  return(x < lower | x > upper)
}

outlier_flags <- apply(fake_data[, c("DO", "temp", "NO3")], 2,
                       detect_outliers_iqr)
apply(outlier_flags, 2, which)

# --- 1b. 可视化异常值 ---
fake_data %>%
  pivot_longer(cols = c(DO, temp, NO3), names_to = "var") %>%
  ggplot(aes(x = var, y = value)) +
  geom_boxplot(outlier.color = "red", outlier.size = 3) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_bw() +
  labs(title = "异常值检测 — Boxplot 法")

# ============================================================
# A2 — 正态性检验与转换
# ============================================================

# --- 2a. Shapiro-Wilk 检验 ---
shapiro.test(fake_data$temp)      # p > 0.05 → 正态
shapiro.test(fake_data$NO3)       # p < 0.05 → 非正态

# --- 2b. 多变量同时检验 ---
apply(fake_data[, c("pH", "DO", "temp", "NO3")], 2,
      function(x) shapiro.test(x)$p.value)

# --- 2c. Box-Cox 转换 ---
bc <- car::powerTransform(fake_data$NO3 + 0.01)  # +0.01 避免 log(0)
fake_data$NO3_bc <- car::bcPower(fake_data$NO3 + 0.01, bc$lambda)

# 转前 vs 转后
par(mfrow = c(1, 2))
hist(fake_data$NO3, main = "转换前", col = "steelblue")
hist(fake_data$NO3_bc, main = paste0("Box-Cox (λ = ", round(bc$lambda, 3), ")"),
     col = "coral")
par(mfrow = c(1, 1))

# ============================================================
# A3 — 数据标准化（群落数据常用）
# ============================================================

# 模拟物种丰度矩阵（10样方 × 15物种）
set.seed(123)
spe_mat <- matrix(rpois(150, lambda = c(rep(5, 5), rep(20, 5), rep(2, 5))),
                  nrow = 10,
                  dimnames = list(paste0("Site", 1:10),
                                  paste0("sp", 1:15)))

# --- 3a. Hellinger 转换（最常用） ---
spe_hel <- vegan::decostand(spe_mat, method = "hellinger")

# --- 3b. Chord 转换 ---
spe_chord <- vegan::decostand(spe_mat, method = "normalize")  # chord = normalize per row

# --- 3c. 标准化方法比较 ---
methods <- c("total", "max", "pa", "hellinger", "wisconsin")
lapply(methods, function(m) {
  decostand(spe_mat, method = m)
}) |> setNames(methods)

# ============================================================
# A4 — 缺失值处理
# ============================================================

# --- 4a. 缺失模式可视化 ---
VIM::aggr(fake_data[, c("pH", "DO", "temp", "NO3")],
          numbers = TRUE, sortVars = TRUE)

# --- 4b. 简单均值填充 ---
impute_mean <- function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}
fake_data_mean <- fake_data
fake_data_mean$pH <- impute_mean(fake_data_mean$pH)
fake_data_mean$NO3 <- impute_mean(fake_data_mean$NO3)

# --- 4c. MICE 多重插补（推荐） ---
library(mice)
imp <- mice(fake_data[, c("pH", "DO", "temp", "NO3")],
            m = 5, method = "pmm", seed = 42)
fake_data_mice <- complete(imp, 1)  # 选取第1套插补

# ============================================================
# A5 — 探索性数据分析 (EDA)
# ============================================================

# --- 5a. 描述统计 ---
psych::describe(fake_data[, c("pH", "DO", "temp", "NO3", "species_richness")])

# --- 5b. 相关性矩阵 ---
cor_mat <- cor(fake_data[, c("pH", "DO", "temp", "NO3")],
               use = "pairwise.complete.obs")
corrplot::corrplot(cor_mat, method = "number", type = "upper")

# --- 5c. 散点图矩阵 ---
GGally::ggpairs(fake_data[, c("pH", "DO", "temp", "NO3")])

# ============================================================
# 导出清洗后的数据
# ============================================================
# write.csv(fake_data_mice, "data_cleaned.csv", row.names = FALSE)
