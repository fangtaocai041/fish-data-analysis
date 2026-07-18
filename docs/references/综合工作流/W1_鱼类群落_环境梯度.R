# ============================================================
# W1 — 鱼类群落对环境梯度的响应：完整分析 Pipeline
# 方陶文库 · 数量生态学 · 综合工作流
#
# 串联：A → B → C → D (VPA) → 可视化
# ============================================================
# 加载包
library(tidyverse)
library(vegan)
library(adespatial)   # MEM
library(corrplot)
library(ggplot2)
library(patchwork)    # 拼图

# ============================================================
# 0. 模拟数据（模拟长江某支流的鱼类调查）
# ============================================================
set.seed(42)
n_sites <- 30

# --- 环境梯度 ---
env <- data.frame(
  site = paste0("S", 1:n_sites),
  pH = c(rnorm(10, 6.8, 0.3),   # 上游
         rnorm(10, 7.2, 0.4),   # 中游
         rnorm(10, 7.8, 0.5)),  # 下游
  temp = c(rnorm(10, 18, 2),    # 上游冷
           rnorm(10, 22, 2),    # 中游
           rnorm(10, 26, 2)),   # 下游暖
  DO = c(rnorm(10, 9, 1),       # 上游高氧
         rnorm(10, 7, 1),       # 中游
         rnorm(10, 5.5, 1)),    # 下游低氧
  depth = c(rnorm(10, 1.5, 0.5),# 上游浅
            rnorm(10, 3, 1),    # 中游
            rnorm(10, 5, 1.5)), # 下游深
  zone = factor(rep(c("上游", "中游", "下游"), each = 10))
)
rownames(env) <- env$site

# --- 模拟鱼类群落（15个物种，3个生态类群）---
# 上游偏好种：喜冷、高氧、浅水
# 下游偏好种：喜暖、耐低氧、深水
# 广布种：到处都有
spp_names <- c("宽鳍鱲", "马口鱼", "拉氏鱥",           # 上游 1-3
               "中华花鳅", "棒花鱼", "麦穗鱼",         # 上游-中游 4-6
               "鲫", "鲤", "鲢", "鳙",                   # 广布 7-10
               "黄颡鱼", "鳊",                            # 中-下游 11-12
               "刀鲚", "鲌", "鳜")                       # 下游 13-15

spe <- matrix(0, nrow = n_sites, ncol = 15,
              dimnames = list(env$site, spp_names))

for (i in 1:n_sites) {
  z <- env[i, ]
  lambda <- numeric(15)
  # 上游种（1-3）：pH低、温度低、DO高、深度浅
  lambda[1:3] <- exp(1.0 - 0.3*(z$pH - 7) - 0.2*(z$temp - 22) + 0.3*(z$DO - 7) - 0.5*(z$depth - 3))
  # 上游-中游种（4-6）：温和偏好
  lambda[4:6] <- exp(0.5 - 0.1*(z$pH - 7) + 0.1*(z$temp - 22) + 0.1*(z$DO - 7))
  # 广布种（7-10）：稳定
  lambda[7:10] <- exp(2.0)
  # 中-下游种（11-12）
  lambda[11:12] <- exp(0.5 + 0.2*(z$temp - 22) - 0.2*(z$DO - 7) + 0.3*(z$depth - 3))
  # 下游种（13-15）：pH高、温度高、DO低、深度深
  lambda[13:15] <- exp(0.5 + 0.3*(z$pH - 7) + 0.3*(z$temp - 22) - 0.3*(z$DO - 7) + 0.4*(z$depth - 3))

  spe[i, ] <- rpois(15, lambda = pmax(lambda, 0.1))
}

# ============================================================
# A — 数据探索与预处理
# ============================================================

# A1 描述统计
psych::describe(env[, c("pH", "temp", "DO", "depth")])

# A2 异常值检查
spe_long <- as.data.frame(spe) %>%
  rownames_to_column("site") %>%
  pivot_longer(-site, names_to = "species", values_to = "abundance")

ggplot(spe_long, aes(x = species, y = abundance)) +
  geom_boxplot(outlier.color = "red") +
  theme_bw() + labs(title = "各物种丰度分布") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# A3 Hellinger 转换（群落数据）
spe_hel <- decostand(spe, method = "hellinger")

# ============================================================
# B — α多样性
# ============================================================

alpha <- data.frame(
  site = rownames(spe),
  S = specnumber(spe),
  Shannon = diversity(spe, index = "shannon"),
  Simpson = diversity(spe, index = "simpson"),
  zone = env$zone
)

# α多样性沿环境梯度的变化
p_alpha <- ggplot(alpha, aes(x = zone, y = Shannon, fill = zone)) +
  geom_boxplot() + geom_jitter(width = 0.2, size = 3) +
  scale_fill_manual(values = c("上游" = "#2E86AB", "中游" = "#A23B72",
                                "下游" = "#F18F01")) +
  labs(title = "α多样性沿河流梯度变化", y = "Shannon 指数") +
  theme_bw() + theme(legend.position = "none")

# α多样性与环境因子的关系
alpha_env <- merge(alpha, env, by = "site")
ggplot(alpha_env, aes(x = temp, y = Shannon, color = zone)) +
  geom_point(size = 3) + geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("上游" = "#2E86AB", "中游" = "#A23B72",
                                 "下游" = "#F18F01")) +
  labs(title = "Shannon 指数 vs 水温") + theme_bw()

# ============================================================
# C — 多元统计分析
# ============================================================

# C1a 梯度长度检查 → 决定 RDA vs CCA
dcaxis <- decorana(spe_hel)
cat("DCA 第一轴长度: ", round(dcaxis$eval[1], 2), " SD\n")
# < 4 SD → RDA

# C1b RDA
env_z <- decostand(env[, c("pH", "temp", "DO", "depth")],
                   method = "standardize")
rda_mod <- rda(spe_hel ~ ., data = env_z)
summary(rda_mod)

# RDA 显著性检验
set.seed(123)
rda_global <- anova(rda_mod, permutations = 9999)
cat("RDA 全局检验 p = ", rda_global$`Pr(>F)`[1], "\n")

# RDA 逐项检验
rda_terms <- anova(rda_mod, by = "terms", permutations = 9999)
rda_terms

# RDA 可视化
plot(rda_mod, scaling = 2, main = "RDA：鱼类群落 ~ 环境因子")

# C1c 变量选择
rda_sel <- ordistep(rda(spe_hel ~ 1, data = env_z),
                    scope = formula(rda_mod),
                    direction = "forward", pstep = 999)
rda_sel$anova

# C2a PERMANOVA
dist_bc <- vegdist(spe, method = "bray")
adonis2(dist_bc ~ zone, data = env, permutations = 9999)

# C2b betadisper（同质性检验）
bd <- betadisper(dist_bc, env$zone)
anova(bd)  # p > 0.05 接受同质性

# C2c SIMPER
simp <- with(env, simper(spe, zone))
summary(simp$上游_vs_下游, ordered = TRUE)

# C3a NMDS 可视化
nmds <- metaMDS(dist_bc, k = 2, try = 100)
nmds_scores <- as.data.frame(scores(nmds)$sites)
nmds_scores$zone <- env$zone

p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = zone)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = c("上游" = "#2E86AB", "中游" = "#A23B72",
                                 "下游" = "#F18F01")) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
           label = paste("Stress =", round(nmds$stress, 3))) +
  labs(title = "NMDS 排序 (Bray-Curtis)") +
  theme_bw()

# C3b 环境因子投影到 NMDS
ef <- envfit(nmds, env_z, permutations = 999)
ef_df <- as.data.frame(ef$vectors$arrows * sqrt(ef$vectors$r))
ef_df$var <- rownames(ef_df)

p_nmds_env <- p_nmds +
  geom_segment(data = ef_df,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "darkred", linewidth = 1) +
  geom_text(data = ef_df,
            aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = var),
            color = "darkred", size = 4)

# ============================================================
# D — 空间结构分解 (MEM + VPA)
# ============================================================

# D1 模拟空间坐标
set.seed(1)
coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100),
  row.names = env$site
)
# 沿河流方向排序（模拟上下游）
coords$x <- as.numeric(env$zone) * 15 + rnorm(n_sites, 0, 3)
coords$y <- rnorm(n_sites, 50, 10)

# D2 MEM 构建
thresh <- give.thresh(dist(coords))
mem <- dbmem(coords, thresh = thresh[1], MEM.autocor = "positive")

# 选择显著 MEM
rda_mem <- rda(spe_hel ~ ., data = as.data.frame(mem))
sel_mem <- ordistep(rda(spe_hel ~ 1, data = as.data.frame(mem)),
                    scope = formula(rda_mem),
                    direction = "forward", pstep = 999)
sel_vars <- attr(sel_mem$terms, "term.labels")

# D3 方差分解 (VPA): 环境 vs 空间
if (length(sel_vars) >= 2) {
  vpa <- varpart(spe_hel, env_z, mem[, sel_vars])
} else {
  vpa <- varpart(spe_hel, env_z, mem)
}
plot(vpa, Xnames = c("环境", "空间"), bg = c("steelblue", "coral"))
title("VPA：环境 vs 空间对群落变异的解释")

# VPA 显著性检验
env_part <- rda(spe_hel, env_z)
space_part <- rda(spe_hel, mem[, sel_vars])
env_cond_space <- rda(spe_hel, env_z, Condition(as.matrix(mem[, sel_vars])))
space_cond_env <- rda(spe_hel, mem[, sel_vars], Condition(as.matrix(env_z)))

cat("\n=== VPA 显著性检验 ===\n")
cat("环境部分（纯）:", anova(env_cond_space, permutations = 999)$`Pr(>F)`[1], "\n")
cat("空间部分（纯）:", anova(space_cond_env, permutations = 999)$`Pr(>F)`[1], "\n")

# ============================================================
# C4 — 指示种分析
# ============================================================
library(indicspecies)
indval <- multipatt(spe, env$zone, func = "r.g",
                    control = how(nperm = 9999))
summary(indval, indvalcomp = TRUE)

# ============================================================
# 发表级可视化拼图
# ============================================================

# 合并图
final_fig <- (p_alpha | p_nmds_env) /
  wrap_elements(~plot(vpa, Xnames = c("环境", "空间"),
                      bg = c("#2E86AB", "#F18F01")))

ggsave("Fig1_鱼类群落分析.pdf", final_fig, width = 12, height = 10)
cat("✅ 分析完成！结果图已保存\n")

# ============================================================
# 结果汇总表
# ============================================================
results <- data.frame(
  Component = c("RDA 总解释率", "环境纯效应", "空间纯效应",
                "环境×空间", "未解释"),
  Explained = round(c(
    sum(rda_mod$CCA$eig) / rda_mod$tot.chi,
    vpa$part$indfract[1, "Adj.R.square"],
    vpa$part$indfract[2, "Adj.R.square"],
    vpa$part$indfract[3, "Adj.R.square"],
    vpa$part$indfract[4, "Adj.R.square"]
  ), 3)
)
cat("\n=== 方差分解结果 ===\n")
print(results)

cat("\n=== 主要指示种 ===\n")
print(indval$sign)
