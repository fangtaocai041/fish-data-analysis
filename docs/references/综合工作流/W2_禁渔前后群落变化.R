# ============================================================
# W2 — 禁渔前后群落结构变化：完整分析 Pipeline
# 方陶文库 · 数量生态学 · 综合工作流
#
# 串联：A → B2 (β分解) → C3 (PERMANOVA) → C4 (SIMPER/IndVal) 
#       → G3 (EWS) → I3 (GAM)
# ============================================================
library(tidyverse)
library(vegan)
library(betapart)
library(indicspecies)
library(earlywarnings)
library(mgcv)
library(ggplot2)
library(patchwork)

# ============================================================
# 0. 模拟数据：模拟禁渔前(2010-2019)和禁渔后(2020-2025)的鱼类监测
# ============================================================
set.seed(123)
n_pre  <- 40   # 禁渔前样方数
n_post <- 30   # 禁渔后样方数
n_spp  <- 20

# 物种名录（长江常见鱼类）
spp_names <- c("鲢","鳙","草鱼","青鱼","鲤","鲫","鳊","鲌",
               "刀鲚","鳜","黄颡鱼","鲶","鳤","铜鱼","赤眼鳟",
               "银鲴","黄尾鲴","似鳊","花䱻","麦穗鱼")

# --- 模拟禁渔前群落 ---
# 禁渔前：少数优势种(鲢/鳙/鲫/鲤)占绝对优势，多样性低
lambda_pre <- c(50, 40, 10, 8, 30, 35, 5, 3,  # 8个常见种，前4种超优势
                2, 1, rep(0.5, 10))           # 12个稀有种
spe_pre <- matrix(0, nrow = n_pre, ncol = n_spp,
                  dimnames = list(paste0("Pre_S", 1:n_pre), spp_names))
for (i in 1:n_pre) {
  spe_pre[i, ] <- rpois(n_spp, lambda = lambda_pre * runif(n_spp, 0.5, 1.5))
}

# --- 模拟禁渔后群落 ---
# 禁渔后：优势种衰退，多样性和均匀度上升
lambda_post <- c(20, 18, 15, 12, 25, 28, 10, 8,     # 优势种下降
                 10, 6, 5, 4, 3, 3, 3, 2, 2, 2, 2, 2)  # 稀有种恢复
spe_post <- matrix(0, nrow = n_post, ncol = n_spp,
                   dimnames = list(paste0("Post_S", 1:n_post), spp_names))
for (i in 1:n_post) {
  spe_post[i, ] <- rpois(n_spp, lambda = lambda_post * runif(n_spp, 0.5, 1.5))
}

# 合并
spe_all <- rbind(spe_pre, spe_post)
period <- factor(c(rep("禁渔前", n_pre), rep("禁渔后", n_post)))

# ============================================================
# A — 数据清洗 & EDA
# ============================================================

cat("\n========== A: 数据概览 ==========\n")
cat("总样方数:", nrow(spe_all), "\n")
cat("总物种数:", ncol(spe_all), "\n")
cat("禁渔前样方:", n_pre, "  禁渔后样方:", n_post, "\n")

# 稀有物种过滤（总出现次数 < 3 的物种）
spp_total <- colSums(spe_all)
rare_spp <- names(spp_total[spp_total < 3])
if (length(rare_spp) > 0) {
  cat("移除稀有物种 (总丰度<3):", paste(rare_spp, collapse=", "), "\n")
  spe_all <- spe_all[, spp_total >= 3]
}

# Hellinger 转换（用于后续 RDA/GAM）
spe_hel <- decostand(spe_all, method = "hellinger")

# ============================================================
# B2 — β多样性分解 (Baselga)
# ============================================================

cat("\n========== B2: β多样性分解 ==========\n")

# 总β多樣性
beta_total <- beta.multi(spe_all, index.family = "sorensen")
cat(sprintf("总β多样性 (Sørensen):    %.3f\n", beta_total$beta.SOR))
cat(sprintf("  周转组分 (β_SIM):       %.3f (%.1f%%)\n",
            beta_total$beta.SIM,
            beta_total$beta.SIM / beta_total$beta.SOR * 100))
cat(sprintf("  嵌套组分 (β_SNE):       %.3f (%.1f%%)\n",
            beta_total$beta.SNE,
            beta_total$beta.SNE / beta_total$beta.SOR * 100))

# 分时段计算β多样性
beta_pre  <- beta.multi(spe_pre, index.family = "sorensen")
beta_post <- beta.multi(spe_post, index.family = "sorensen")

cat(sprintf("\n禁渔前 β多样性: %.3f (周转: %.1f%%)\n",
            beta_pre$beta.SOR,
            beta_pre$beta.SIM / beta_pre$beta.SOR * 100))
cat(sprintf("禁渔后 β多样性: %.3f (周转: %.1f%%)\n",
            beta_post$beta.SOR,
            beta_post$beta.SIM / beta_post$beta.SOR * 100))

# ============================================================
# C3 — PERMANOVA + betadisper
# ============================================================

cat("\n========== C3: 组间差异检验 ==========\n")

bc_dist <- vegdist(spe_all, method = "bray")

# Step 1: 同质性检验
bd <- betadisper(bc_dist, period)
bd_anova <- anova(bd)
cat(sprintf("betadisper: F = %.3f, p = %.4f\n",
            bd_anova$`F value`[1], bd_anova$`Pr(>F)`[1]))

if (bd_anova$`Pr(>F)`[1] < 0.05) {
  cat("⚠️ 组内离散度不相等！PERMANOVA显著需谨慎解释\n")
} else {
  cat("✅ 组内离散度齐性\n")
}

# Step 2: PERMANOVA
perm <- adonis2(bc_dist ~ period, permutations = 9999)
cat(sprintf("PERMANOVA: R² = %.3f, F = %.3f, p = %.4f\n",
            perm$R2[1], perm$F[1], perm$`Pr(>F)`[1]))

# ============================================================
# C4 — SIMPER + IndVal
# ============================================================

cat("\n========== C4: 关键变化物种 ==========\n")

# SIMPER
simp <- simper(spe_all, period, permutations = 999)
top_spp <- head(summary(simp$禁渔前_禁渔后, ordered = TRUE), 8)

cat("SIMPER — 对组间差异贡献最大的物种:\n")
for (i in 1:nrow(top_spp)) {
  cat(sprintf("  %-10s 贡献率=%.1f%%  禁渔前均丰=%.1f  禁渔后均丰=%.1f  p=%.3f\n",
              rownames(top_spp)[i],
              top_spp[i, "average"] / sum(top_spp$average) * 100,
              top_spp[i, "ava"], top_spp[i, "avb"],
              top_spp[i, "p"]))
}

# IndVal
indval <- multipatt(spe_all, period, func = "r.g", 
                    control = how(nperm = 9999))
cat("\nIndVal — 各时期的指示种:\n")
print(summary(indval, indvalcomp = TRUE))

# ============================================================
# G3 — 早期预警信号 (EWS)
# ============================================================

cat("\n========== G3: 早期预警信号检查 ==========\n")

# 模拟时间序列（年际总丰度）
set.seed(456)
years <- 2010:2025
n_years <- length(years)

# 构造总丰度序列：禁渔前稳定偏低 → 2020禁渔 → 逐渐恢复
total_abundance <- c(
  # 2010-2013: 下降趋势
  seq(800, 550, length.out = 4) + rnorm(4, 0, 30),
  # 2014-2016: 低谷波动
  rnorm(3, 450, 25),
  # 2017-2019: 小幅回升（禁渔前已有恢复信号？）
  seq(480, 600, length.out = 3) + rnorm(3, 0, 20),
  # 2020-2025: 禁渔后快速恢复
  seq(650, 1200, length.out = 6) + rnorm(6, 0, 40)
)

ts_data <- data.frame(year = years, abundance = total_abundance)

# EWS 分析
ews_out <- generic_ews(timeseries = ts_data$abundance, 
                        winsize = 0.5, detrending = "gaussian")

cat("早期预警信号指标:\n")
cat(sprintf("  AR(1):   %.3f\n", ews_out$ar1))
cat(sprintf("  SD:      %.3f\n", ews_out$sd))
cat(sprintf("  Skewness: %.3f\n", ews_out$sk))
cat(sprintf("  Kurtosis: %.3f\n", ews_out$kurt))
cat("→ 如果 AR(1) 和 SD 在禁渔前持续增大 → 可能有临界减慢信号\n")

# ============================================================
# I3 — GAM 时间趋势（优势种年际变化）
# ============================================================

cat("\n========== I3: GAM 时间趋势 ==========\n")

# 选前3个优势种
dominant_spp <- names(sort(colMeans(spe_all), decreasing = TRUE)[1:3])
cat("分析优势种:", paste(dominant_spp, collapse = ", "), "\n")

# 为每个优势种拟合 GAM
gam_results <- list()
for (sp in dominant_spp) {
  # 构造该物种在各年的平均丰度（模拟）
  set.seed(which(spp_names == sp))
  spp_ts <- c(
    rnorm(5, mean = sort(runif(5, 30, 80), decreasing = TRUE), sd = 10),  # 禁渔前：下降
    rnorm(6, mean = sort(runif(6, 15, 60), decreasing = FALSE), sd = 10), # 禁渔后：恢复
    rnorm(5, mean = sort(runif(5, 40, 90), decreasing = FALSE), sd = 10)  # 持续恢复
  )
  
  gam_mod <- gam(spp_ts ~ s(1:16, k = 5), method = "REML")
  gam_results[[sp]] <- list(
    model = gam_mod,
    fitted = fitted(gam_mod),
    edf = summary(gam_mod)$edf
  )
  
  cat(sprintf("  %s: EDF = %.1f (1=线性, >1=非线性)\n", 
              sp, gam_results[[sp]]$edf))
}

# ============================================================
# 可视化
# ============================================================

# P1: NMDS 禁渔前后比较
nmds <- metaMDS(bc_dist, k = 2, try = 100)
nmds_scores <- as.data.frame(scores(nmds)$sites)
nmds_scores$period <- period

p1 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = period)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 1.2) +
  scale_color_manual(values = c("禁渔前" = "#D73027", "禁渔后" = "#4575B4")) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
           label = paste("Stress =", round(nmds$stress, 3))) +
  labs(title = "NMDS: 禁渔前后群落结构", 
       subtitle = paste("PERMANOVA p =", round(perm$`Pr(>F)`[1], 4))) +
  theme_bw()

# P2: 时间序列 + GAM 趋势
p2 <- ggplot(ts_data, aes(x = year, y = abundance)) +
  geom_point(size = 3, color = "steelblue") +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), 
              se = TRUE, fill = "grey70") +
  geom_vline(xintercept = 2019.5, linetype = "dashed", 
             color = "red", linewidth = 1) +
  annotate("text", x = 2019.5, y = max(ts_data$abundance) * 0.95,
           label = "禁渔开始", hjust = -0.1, color = "red") +
  labs(title = "总丰度年际变化 (GAM平滑)", y = "总丰度") +
  theme_bw()

# P3: α 多样性变化
alpha_df <- data.frame(
  period = period,
  Shannon = diversity(spe_all, "shannon"),
  Simpson = diversity(spe_all, "simpson")
)

p3 <- ggplot(alpha_df, aes(x = period, y = Shannon, fill = period)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
  scale_fill_manual(values = c("禁渔前" = "#D73027", "禁渔后" = "#4575B4")) +
  labs(title = "Shannon 多样性: 禁渔前 vs 后", y = "Shannon 指数") +
  theme_bw() + theme(legend.position = "none")

# P4: β多样性分解饼图
beta_viz <- data.frame(
  组分 = c("周转", "嵌套"),
  比例 = c(beta_total$beta.SIM, beta_total$beta.SNE) / beta_total$beta.SOR
)

p4 <- ggplot(beta_viz, aes(x = "", y = 比例, fill = 组分)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("周转" = "#2E86AB", "嵌套" = "#F18F01")) +
  labs(title = "β多样性分解", subtitle = paste0("总β = ", round(beta_total$beta.SOR, 3))) +
  theme_void()

# 合并图
final_fig <- (p1 | p3) / (p2 | p4) +
  plot_annotation(
    title = "W2: 禁渔对鱼类群落结构的影响",
    caption = "方陶文库 · 数量生态学"
  )

ggsave("Fig_W2_禁渔效应.pdf", final_fig, width = 14, height = 10)
cat("\n✅ W2 分析完成！结果图已保存\n")
