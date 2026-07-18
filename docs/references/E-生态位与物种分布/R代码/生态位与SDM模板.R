# ============================================================
# E — 生态位与物种分布模型 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(ecospat)      # 生态位比较
library(dismo)        # SDM 框架
library(maxnet)       # MaxEnt 实现
library(ggplot2)

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(42)
n <- 200  # 背景点数量

# 两个环境变量
env <- data.frame(
  temp = rnorm(n, mean = 20, sd = 5),
  precip = rnorm(n, mean = 1000, sd = 300)
)
# 假设的物种存在概率（响应曲面）
prob <- with(env, {
  z <- -1 + 0.05*temp - 0.002*temp^2 + 0.002*precip - 0.000001*precip^2
  plogis(z)
})
presence <- rbinom(n, 1, prob)

# 分离存在点和背景点
occ <- env[presence == 1, ]
bg <- env[presence == 0, ]
n_occ <- nrow(occ)
cat("存在点:", n_occ, "  背景点:", nrow(bg))

# ============================================================
# E1 — 生态位度量与比较
# ============================================================

# --- 1a. 生态位重叠 ---
# 模拟两个物种
sp1_env <- data.frame(
  temp = rnorm(50, 18, 3),
  precip = rnorm(50, 1000, 200)
)
sp2_env <- data.frame(
  temp = rnorm(50, 22, 4),
  precip = rnorm(50, 800, 250)
)

# 计算 PCA 环境空间
# combine sp1 + sp2 背景环境
all_env <- rbind(sp1_env, sp2_env)
pca_env <- prcomp(all_env, scale = TRUE, center = TRUE)

# 投影到 PCA 空间
sp1_proj <- predict(pca_env, sp1_env)
sp2_proj <- predict(pca_env, sp2_env)

# 网格化 PCA 空间
ecospat::ecospat.grid.clim.dyn(
  glob = all_env,
  glob1 = all_env,
  sp1 = sp1_env,
  sp2 = sp2_env,
  R = 100
)

# 计算 Schoener's D 生态位重叠
niche_overlap <- ecospat::ecospat.niche.overlap(
  ecospat::ecospat.grid.clim.dyn(sp1_env, all_env, R = 100),
  ecospat::ecospat.grid.clim.dyn(sp2_env, all_env, R = 100),
  cor = TRUE
)
niche_overlap$D  # Schoener's D
niche_overlap$I  # Hellinger-based I

# ============================================================
# E2 — MaxEnt 物种分布模型
# ============================================================

# --- 2a. 准备数据 ---
# 存在 = 1, 背景 = 0
occ_data <- occ
bg_data <- bg[sample(1:nrow(bg), min(1000, nrow(bg))), ]
pa_data <- rbind(
  cbind(occ_data, pa = 1),
  cbind(bg_data, pa = 0)
)

# --- 2b. 拟合 MaxEnt 模型 ---
max_mod <- maxnet(
  p = pa_data$pa,
  data = pa_data[, c("temp", "precip")],
  f = maxnet.formula(pa_data$pa, pa_data[, c("temp", "precip")])
)

# --- 2c. 变量响应曲线 ---
plot(max_mod, type = "cloglog")

# --- 2d. 变量重要性 ---
# 计算置换重要性
imp <- maxnet::maxnet.genericModel(max_mod)
# 查看系数
coef(max_mod)

# --- 2e. 预测到空间 ---
# 模拟预测栅格
pred_grid <- expand.grid(
  temp = seq(5, 35, length = 50),
  precip = seq(200, 2000, length = 50)
)
pred_grid$suitability <- predict(max_mod, pred_grid, type = "cloglog")

ggplot(pred_grid, aes(x = temp, y = precip, fill = suitability)) +
  geom_raster() +
  scale_fill_viridis_c(name = "适宜度") +
  labs(title = "MaxEnt 预测的物种分布适宜度") +
  theme_bw()

# ============================================================
# E3 — 模型评估
# ============================================================

# --- Bootstrap 评估 ---
evaluate_model <- function(model, test_pa, test_env) {
  pred <- predict(model, test_env, type = "cloglog")
  # 使用 dismo::evaluate
  eval <- dismo::evaluate(
    p = pred[test_pa == 1],
    a = pred[test_pa == 0]
  )
  return(data.frame(
    AUC = eval@auc,
    TSS = max(eval@TPR + eval@TNR - 1),
    threshold = threshold(eval)
  ))
}

# 分割训练/测试
set.seed(123)
train_idx <- sample(1:nrow(pa_data), size = round(0.7 * nrow(pa_data)))
train <- pa_data[train_idx, ]
test <- pa_data[-train_idx, ]

max_mod_cv <- maxnet(train$pa, train[, c("temp", "precip")],
                     f = maxnet.formula(train$pa, train[, c("temp", "precip")]))
evaluate_model(max_mod_cv, test$pa, test[, c("temp", "precip")])
