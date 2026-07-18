# ============================================================
# I — 统计建模与机器学习 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(lme4)         # GLMM
library(mgcv)         # GAM
library(glmmTMB)      # 灵活 GLMM（零膨胀等）
library(randomForest) # 随机森林
library(DHARMa)       # 残差诊断
library(ggplot2)

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(42)
n <- 200
n_sites <- 20

# 嵌套结构：观测 (n=200) 嵌套在样方 (n_sites=20)
df <- data.frame(
  site = factor(rep(paste0("Site", 1:n_sites), each = n / n_sites)),
  pH = rnorm(n, 7, 0.8),
  temp = rnorm(n, 22, 4),
  depth = rnorm(n, 5, 2)
)

# 模拟响应：物种丰富度（计数）
# 固定效应：pH + temp + depth + temp:depth
# 随机效应：site 的随机截距
site_effects <- rnorm(n_sites, 0, 2)
df$richness <- with(df, {
  mu <- exp(2.0 + 0.3 * pH - 0.05 * temp +
              0.1 * depth + 0.01 * temp * depth +
              site_effects[site])
  rpois(n, mu)
})

# ============================================================
# I1 — GLM
# ============================================================

# --- 1a. Poisson GLM ---
glm_pois <- glm(richness ~ pH * temp + depth, data = df, family = poisson)
summary(glm_pois)

# --- 1b. 过离散检查 ---
deviance(glm_pois) / df.residual(glm_pois)  # > 1.5 → 过离散

# --- 1c. 负二项 GLM ---
library(MASS)
glm_nb <- glm.nb(richness ~ pH * temp + depth, data = df)
summary(glm_nb)

# ============================================================
# I2 — GLMM (混合效应模型)
# ============================================================

# --- 2a. Poisson GLMM ---
glmm_pois <- glmer(richness ~ pH + temp + depth + (1 | site),
                    data = df, family = poisson)
summary(glmm_pois)

# --- 2b. 负二项 GLMM ---
glmm_nb <- glmer.nb(richness ~ pH + temp + depth + (1 | site),
                     data = df)
summary(glmm_nb)

# --- 2c. glmmTMB（更灵活：零膨胀 + NB）---
library(glmmTMB)
glmm_zinb <- glmmTMB(richness ~ pH + temp + depth + (1 | site),
                     data = df,
                     ziformula = ~1,        # 零膨胀部分
                     family = nbinom2)
summary(glmm_zinb)

# --- 2d. 模型比较 ---
AIC(glm_pois, glm_nb, glmm_pois, glmm_nb, glmm_zinb)

# --- 2e. DHARMa 残差诊断 ---
sim_res <- simulateResiduals(glmm_nb)
plot(sim_res)
testDispersion(sim_res)

# ============================================================
# I3 — GAM (广义加性模型)
# ============================================================

# --- 3a. GAM 用 Poisson 族 ---
gam_mod <- gam(richness ~ s(pH, k = 5) + s(temp, k = 5) +
                 s(depth, k = 5) +
                 te(temp, depth, k = 4),   # 张量积平滑（交互）
               data = df, family = poisson, method = "REML")
summary(gam_mod)

# --- 3b. GAM 可视化 ---
plot(gam_mod, pages = 1, shade = TRUE, scale = 0)

# --- 3c. GAM 空间平滑示例 ---
# 如果有空间坐标
# gam_spatial <- gam(richness ~ s(x, y, k = 20) + s(pH) + s(temp),
#                    data = df, family = poisson)

# ============================================================
# I5 — 随机森林
# ============================================================

# --- 5a. 随机森林回归 ---
rf_mod <- randomForest(richness ~ pH + temp + depth,
                       data = df, ntree = 500,
                       importance = TRUE)
rf_mod
# 解释方差
plot(rf_mod)  # 树的数量 vs 误差

# --- 5b. 变量重要性 ---
importance(rf_mod)
varImpPlot(rf_mod)

# --- 5c. 偏依赖图 ---
partialPlot(rf_mod, df, "pH", main = "pH 对丰富度的偏依赖")
partialPlot(rf_mod, df, "temp", main = "温度对丰富度的偏依赖")

# ============================================================
# 模型预测和比较
# ============================================================

# 交叉验证预测
set.seed(123)
train_idx <- sample(1:n, size = round(0.7 * n))
train <- df[train_idx, ]
test <- df[-train_idx, ]

# GLM
glm_fit <- glm(richness ~ pH + temp + depth, data = train, family = poisson)
pred_glm <- predict(glm_fit, test, type = "response")

# GAM
gam_fit <- gam(richness ~ s(pH) + s(temp) + s(depth),
               data = train, family = poisson)
pred_gam <- predict(gam_fit, test, type = "response")

# RF
rf_fit <- randomForest(richness ~ pH + temp + depth, data = train)
pred_rf <- predict(rf_fit, test)

# 比较 RMSE
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
data.frame(
  GLM = rmse(test$richness, pred_glm),
  GAM = rmse(test$richness, pred_gam),
  RF = rmse(test$richness, pred_rf)
)
