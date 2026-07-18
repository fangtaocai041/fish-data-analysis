# ============================================================
# W3 — 功能多样性对干扰的响应：完整分析 Pipeline
# 方陶文库 · 数量生态学 · 综合工作流
#
# 串联：A (数据清洗) → B4 (功能多样性) → I (GLM/GAM)
# ============================================================
library(tidyverse)
library(vegan)
library(FD)
library(fundiversity)
library(mgcv)
library(ggplot2)
library(patchwork)
library(corrplot)

# ============================================================
# 0. 模拟数据：沿干扰梯度的鱼类群落 + 功能性状
# ============================================================
set.seed(789)
n_sites <- 40
n_spp <- 18

spp_names <- paste0("sp", 1:n_spp)

# --- 干扰梯度 ---
disturbance <- c(
  seq(0.1, 0.4, length.out = 12),   # 低干扰
  seq(0.5, 1.0, length.out = 16),   # 中高干扰
  seq(1.0, 2.0, length.out = 12)    # 高干扰
)[1:n_sites]

env <- data.frame(
  site = paste0("S", 1:n_sites),
  disturbance = disturbance,
  habitat_het = 1 / (disturbance + 0.3) + rnorm(n_sites, 0, 0.1),  # 干扰降低生境异质性
  productivity = 100 - 40 * disturbance + rnorm(n_sites, 0, 8)       # 干扰降低生产力
)
rownames(env) <- env$site

# --- 功能性状 ---
traits <- data.frame(
  row.names = spp_names,
  body_size    = rnorm(n_spp, mean = 20, sd = 10),          # 体型 (cm)
  trophic_level = round(runif(n_spp, 1.5, 4.5), 1),        # 营养级
  fecundity    = rlnorm(n_spp, meanlog = 8, sdlog = 1.5),  # 繁殖力
  tolerance    = rbeta(n_spp, 2, 3) * 10,                   # 耐受力 (0-10)
  habitat_use  = sample(c("底栖","中上层","广栖"), n_spp, 
                        replace = TRUE, prob = c(0.4, 0.4, 0.2))
)

# --- 模拟群落矩阵 ---
# 高干扰 → 耐受种占优 → 功能多样性低
# 低干扰 → 各种生态位共存 → 功能多样性高
spe <- matrix(0, nrow = n_sites, ncol = n_spp,
              dimnames = list(env$site, spp_names))

for (i in 1:n_sites) {
  d <- disturbance[i]
  for (j in 1:n_spp) {
    # 耐受种在高干扰下仍有较高丰度
    tol_effect <- traits$tolerance[j] / 5
    # 大体型在低干扰下更丰盛
    size_effect <- if(d > 1) -0.3 * traits$body_size[j] / 20 else 0.1 * traits$body_size[j] / 20
    lambda <- exp(2 + tol_effect * d - 2 * d + size_effect)
    spe[i, j] <- rpois(1, lambda = pmax(lambda, 0.05))
  }
}

# ============================================================
# A — 数据准备
# ============================================================

cat("\n========== A: 数据概览 ==========\n")
cat("样方数:", n_sites, "  物种数:", n_spp, "\n")
cat("干扰梯度范围:", round(min(disturbance), 2), "-", round(max(disturbance), 2), "\n")

# 稀有物种过滤
spp_total <- colSums(spe)
spe <- spe[, spp_total >= 5]
cat("过滤后物种数:", ncol(spe), "\n")

# 性状数据对齐
traits <- traits[colnames(spe), ]

# 处理分类性状 → 转为哑变量
traits_num <- traits
if ("habitat_use" %in% names(traits_num)) {
  hab_dummies <- model.matrix(~ habitat_use - 1, data = traits_num)
  colnames(hab_dummies) <- c("hab_benthic", "hab_pelagic", "hab_general")
  traits_num <- cbind(traits_num[, setdiff(names(traits_num), "habitat_use")],
                       hab_dummies)
}

# 性状标准化
traits_scaled <- decostand(traits_num, method = "standardize")

# ============================================================
# B4 — 功能多样性计算
# ============================================================

cat("\n========== B4: 功能多样性 ==========\n")

# 检查是否能算 FRic
if (ncol(spe) > ncol(traits_scaled) + 1) {
  fd_out <- dbFD(traits_scaled, spe, w.abun = TRUE, calc.FRic = TRUE)
  cat("✅ 物种数 > 性状数+1，FRic 可计算\n")
} else {
  fd_out <- dbFD(traits_scaled, spe, w.abun = TRUE, calc.FRic = FALSE)
  cat("⚠️ S ≤ T+1，跳过 FRic，仅用 FDis\n")
}

fd_df <- data.frame(
  site = env$site,
  disturbance = env$disturbance,
  FRic = if ("FRic" %in% names(fd_out)) fd_out$FRic else NA,
  FEve = fd_out$FEve,
  FDiv = fd_out$FDiv,
  FDis = fd_out$FDis,
  RaoQ = fd_out$RaoQ
)

# 同时算物种多样性作为对比
fd_df$richness <- specnumber(spe)
fd_df$Shannon <- diversity(spe, "shannon")

cat("功能多样性指数范围:\n")
cat(sprintf("  FDis:  %.2f - %.2f\n", min(fd_df$FDis, na.rm=TRUE), max(fd_df$FDis, na.rm=TRUE)))
if (!all(is.na(fd_df$FRic))) {
  cat(sprintf("  FRic:  %.2f - %.2f\n", min(fd_df$FRic, na.rm=TRUE), max(fd_df$FRic, na.rm=TRUE)))
}
cat(sprintf("  RaoQ:  %.2f - %.2f\n", min(fd_df$RaoQ, na.rm=TRUE), max(fd_df$RaoQ, na.rm=TRUE)))

# ============================================================
# I — 统计建模
# ============================================================

cat("\n========== I: 功能多样性 ~ 干扰梯度模型 ==========\n")

# I1: GAM — 功能离散度与干扰的关系
gam_fdis <- gam(FDis ~ s(disturbance, k = 5), data = fd_df, method = "REML")
cat(sprintf("GAM (FDis ~ disturbance): EDF = %.2f, R² = %.3f, p = %.4f\n",
            summary(gam_fdis)$edf, summary(gam_fdis)$r.sq,
            summary(gam_fdis)$s.pv))

# I1: GAM — RaoQ与干扰的关系
gam_raoq <- gam(RaoQ ~ s(disturbance, k = 5), data = fd_df, method = "REML")
cat(sprintf("GAM (RaoQ ~ disturbance): EDF = %.2f, R² = %.3f, p = %.4f\n",
            summary(gam_raoq)$edf, summary(gam_raoq)$r.sq,
            summary(gam_raoq)$s.pv))

# 对比：物种多样性与干扰
gam_shan <- gam(Shannon ~ s(disturbance, k = 5), data = fd_df, method = "REML")
cat(sprintf("GAM (Shannon ~ disturbance): EDF = %.2f, R² = %.3f, p = %.4f\n",
            summary(gam_shan)$edf, summary(gam_shan)$r.sq,
            summary(gam_shan)$s.pv))

# 功能和物种多样性是否"解耦"？
fd_df$FD_resid <- residuals(gam_fdis)
fd_df$S_resid  <- residuals(gam_shan)
decoupling_r <- cor(fd_df$FD_resid, fd_df$S_resid, use = "complete.obs")
cat(sprintf("\n功能 vs 物种多样性残差相关: r = %.3f\n", decoupling_r))
if (abs(decoupling_r) < 0.3) {
  cat("→ ⚡ 功能与物种多样性'解耦'！功能多样性提供了物种数之外的独特信息\n")
}

# ============================================================
# 可视化
# ============================================================

# P1: 功能多样性 vs 干扰梯度
p1 <- ggplot(fd_df, aes(x = disturbance, y = FDis)) +
  geom_point(size = 3, alpha = 0.7, color = "steelblue") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), 
              color = "darkred", fill = "grey80") +
  labs(title = "功能离散度 (FDis) 沿干扰梯度",
       x = "干扰强度", y = "FDis") +
  theme_bw()

# P2: 功能 vs 物种多样性
p2 <- ggplot(fd_df, aes(x = Shannon, y = FDis)) +
  geom_point(size = 3, alpha = 0.7, aes(color = disturbance)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  scale_color_viridis_c() +
  labs(title = "功能多样性 vs 物种多样性",
       subtitle = paste("r =", round(cor(fd_df$Shannon, fd_df$FDis, use="c"), 3)),
       x = "Shannon 指数", y = "FDis") +
  theme_bw()

# P3: RaoQ vs 干扰
p3 <- ggplot(fd_df, aes(x = disturbance, y = RaoQ)) +
  geom_point(size = 3, alpha = 0.7, color = "coral") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              color = "darkblue", fill = "grey80") +
  labs(title = "Rao's Q 沿干扰梯度",
       x = "干扰强度", y = "Rao's Q") +
  theme_bw()

# P4: 功能多样性指数间相关性
fd_cor <- cor(fd_df[, c("FRic","FEve","FDiv","FDis","RaoQ")], use = "complete.obs")
corrplot(fd_cor, method = "number", type = "upper",
         title = "功能多样性指数间的相关性", mar = c(0,0,2,0))

# 合并图
final_fig <- (p1 | p3) / (p2) +
  plot_annotation(
    title = "W3: 功能多样性对干扰梯度的响应",
    caption = "方陶文库 · 数量生态学"
  )

ggsave("Fig_W3_功能多样性.pdf", final_fig, width = 12, height = 10)
cat("\n✅ W3 分析完成！\n")
