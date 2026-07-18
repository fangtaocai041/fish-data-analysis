# ============================================================
# D — 空间生态学 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(vegan)
library(spdep)       # 空间自相关
library(gstat)       # 半方差图
library(adespatial)  # MEM/PCNM
library(ape)         # Moran's I
library(ggplot2)

# ============================================================
# 0. 模拟空间数据
# ============================================================
set.seed(42)
n <- 40
coords <- data.frame(
  x = runif(n, 0, 100),
  y = runif(n, 0, 100)
)
rownames(coords) <- paste0("Site", 1:n)

# 模拟物种数据（含空间梯度）
easting_effect <- scale(coords$x)[, 1]
northing_effect <- scale(coords$y)[, 1]

spe <- matrix(0, nrow = n, ncol = 15,
              dimnames = list(paste0("Site", 1:n), paste0("sp", 1:15)))
for (i in 1:15) {
  # 有些物种有空间梯度，有些没有
  if (i <= 5) {
    mu <- 10 + 3 * easting_effect - 2 * northing_effect
  } else if (i <= 10) {
    mu <- 5 - 4 * northing_effect
  } else {
    mu <- 8
  }
  mu <- pmax(mu, 0.5)
  spe[, i] <- rpois(n, lambda = mu)
}

# ============================================================
# D1 — 空间自相关
# ============================================================

# --- 1a. 构建空间权重矩阵 (KNN, k=4) ---
nb <- knn2nb(knearneigh(as.matrix(coords), k = 4))
listw <- nb2listw(nb, style = "W")

# --- 1b. 全局 Moran's I ---
# 对 α多样性
alpha <- specnumber(spe)
moran.test(alpha, listw)

# 对 Shannon 多样性
shannon <- diversity(spe)
moran.test(shannon, listw)

# --- 1c. 局部 Moran's I (LISA) ---
local_m <- localmoran(alpha, listw)
head(local_m)

# --- 1d. 半方差图 ---
# 任意一个环境变量（如模拟环境梯度）
env_gradient <- 10 + 2 * easting_effect + rnorm(n, 0, 1)

g <- gstat(formula = env_gradient ~ 1, locations = ~ x + y)
v <- variogram(g)
plot(v, main = "半方差图")

# 拟合理论半方差模型
v_fit <- fit.variogram(v, vgm("Sph"))
plot(v, v_fit, main = "拟合球形半方差模型")

# --- 1e. Mantel 检验（空间 vs 群落距离）---
geo_dist <- dist(coords)
comm_dist <- vegdist(spe, method = "bray")
mantel(comm_dist, geo_dist, method = "pearson", permutations = 9999)

# ============================================================
# D2 — MEM (空间特征变量)
# ============================================================

# --- 2a. db-MEM（默认方法）---
# 使用距离阈值（保持连通的最小距离）
thresh <- give.thresh(dist(coords))
mem <- dbmem(coords, thresh = thresh[1])
summary(mem)

# --- 2b. 选择显著的MEM ---
spe_hel <- decostand(spe, method = "hellinger")
rda_mem <- rda(spe_hel, as.data.frame(mem))
anova(rda_mem)                  # 整体显著性
sel_mem <- ordistep(rda(spe_hel ~ 1, data = as.data.frame(mem)),
                    scope = formula(rda_mem),
                    direction = "forward", pstep = 999)

# --- 2c. 空间 vs 环境 方差分解 ---
# 假设有一些环境变量
env <- data.frame(
  pH = rnorm(n, 7, 0.5),
  Temp = rnorm(n, 22, 3)
)

# VPA: 环境 | 空间 | 共享
vpa_spatial <- varpart(spe_hel, env,
                       as.data.frame(mem[, attr(sel_mem$terms, "term.labels")]))
plot(vpa_spatial)

# ============================================================
# D3 — 空间回归
# ============================================================

# --- 3a. 回归残差的空间自相关检验 ---
env_data <- data.frame(
  pH = rnorm(n, 7, 0.5),
  Temp = rnorm(n, 22, 3)
)
lm_mod <- lm(alpha ~ pH + Temp, data = env_data)
lm_resid <- residuals(lm_mod)
moran.test(lm_resid, listw)  # 显著 = 残差仍有空间结构

# --- 3b. 空间自回归 (SAR) ---
library(spatialreg)
sar_mod <- lagsarlm(alpha ~ pH + Temp, data = env_data, listw)
summary(sar_mod)

# --- 3c. 空间 GAM ---
library(mgcv)
gam_spatial <- gam(alpha ~ pH + Temp + s(x, y, k = 10),
                   data = cbind(env_data, coords))
summary(gam_spatial)
moran.test(residuals(gam_spatial), listw)  # 检查残差

# ============================================================
# D4 — 点格局分析 (简要)
# ============================================================
library(spatstat)
# 模拟随机点 + 聚集点
pp_random <- rpoispp(50, win = owin(c(0, 100), c(0, 100)))
pp_clust <- rMatClust(20, 0.1, 5)

par(mfrow = c(1, 2))
plot(pp_random, main = "完全随机")
plot(pp_clust, main = "聚集")
par(mfrow = c(1, 1))

# Ripley's K
K_random <- Kest(pp_random)
K_clust <- Kest(pp_clust)
plot(K_random, main = "K 函数 — 随机")
plot(K_clust, main = "K 函数 — 聚集")
