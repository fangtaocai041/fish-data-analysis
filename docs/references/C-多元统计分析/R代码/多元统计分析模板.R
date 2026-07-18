# ============================================================
# C — 多元统计分析 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(vegan)
library(ggplot2)
library(indicspecies)  # 指示种分析

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(42)
n_sites <- 30
n_species <- 20
n_env <- 5

# 物种数据（群落矩阵）
spe <- matrix(0, nrow = n_sites, ncol = n_species,
              dimnames = list(paste0("Site", 1:n_sites),
                              paste0("sp", 1:n_species)))
# 模拟3个群落类型
comm_type <- factor(rep(c("TypeA", "TypeB", "TypeC"), each = 10))
for (i in 1:n_sites) {
  base <- switch(as.character(comm_type[i]),
                 TypeA = rpois(n_species, lambda = c(rep(10, 5), rep(1, 15))),
                 TypeB = rpois(n_species, lambda = c(rep(1, 10), rep(8, 10))),
                 TypeC = rpois(n_species, lambda = rep(3, 20)))
  spe[i, ] <- base
}

# 环境变量
env <- data.frame(
  pH = rnorm(n_sites, 7, 0.5),
  Temp = rnorm(n_sites, 22, 3),
  DO = rnorm(n_sites, 8, 1.5),
  Depth = rnorm(n_sites, 5, 2),
  NO3 = rlnorm(n_sites, 1, 0.8)
)
rownames(env) <- rownames(spe)

# ============================================================
# C1 — 排序分析
# ============================================================

# --- 1a. 计算梯度长度 → 决定线性 or 单峰 ---
spe_hel <- decostand(spe, method = "hellinger")
dcaxis <- decorana(spe_hel)
dcaxis$eval  # DCA 轴特征值
# 第一轴长度 (SD) < 4 → RDA 优于 CCA

# --- 1b. PCA (非约束线性) ---
pca <- rda(spe_hel)  # rda() 无约束时 = PCA
summary(pca)         # 查看特征值、累积解释率
screeplot(pca)       # 碎石图
biplot(pca, scaling = 2, main = "PCA 双标图")

# --- 1c. RDA (约束线性排序) ---
env_z <- decostand(env, method = "standardize")
rda_mod <- rda(spe_hel ~ ., data = env_z)
summary(rda_mod)
anova(rda_mod)       # 整体显著性
vif.cca(rda_mod)     # 共线性诊断

# --- 1d. 变量选择 (前向选择) ---
rda_sel <- ordistep(rda(spe_hel ~ 1, data = env_z),
                    scope = formula(rda(spe_hel ~ ., data = env_z)),
                    direction = "forward", pstep = 999)
rda_sel

# --- 1e. 方差分解 (VPA) ---
# 假设两组环境因子
env1 <- env_z[, c("pH", "Temp")]
env2 <- env_z[, c("DO", "Depth", "NO3")]
vpa <- varpart(spe_hel, env1, env2)
plot(vpa)
vpa$part$indfract  # 各部分的独立解释率

# --- 1f. 方差分解显著性检验 ---
anova(rda(spe_hel, env1))
anova(rda(spe_hel, env2))
anova(rda(spe_hel, env1, Condition(as.matrix(env2))))  # 偏RDA

# --- 1g. CCA (单峰约束排序) ---
cca_mod <- cca(spe ~ ., data = env_z)
summary(cca_mod)
anova(cca_mod, by = "terms")  # 逐变量检验

# --- 1h. NMDS ---
bc_dist <- vegdist(spe, method = "bray")
nmds <- metaMDS(bc_dist, k = 2, try = 100)
nmds$stress  # <0.15 为好

nmds_scores <- as.data.frame(scores(nmds)$sites)
nmds_scores$comm_type <- comm_type

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = comm_type)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
           label = paste("Stress =", round(nmds$stress, 3))) +
  theme_bw() + labs(title = "NMDS 排序")

# 同时用 envfit 将环境因子投影到 NMDS 上
ef <- envfit(nmds, env_z, permutations = 999)
ef
plot(ef)  # 需要基础R图形

# ============================================================
# C2 — 聚类分析
# ============================================================

# --- 2a. 层次聚类 (UPGMA) ---
hc <- hclust(bc_dist, method = "average")  # UPGMA
plot(hc, main = "UPGMA 聚类树")
rect.hclust(hc, k = 3, border = 2:4)

# --- 2b. 聚类+排序组合图 ---
plot(nmds$points, type = "n", main = "聚类 + NMDS")
text(nmds$points, labels = rownames(spe), col = cutree(hc, k = 3))

# ============================================================
# C3 — 组间差异检验
# ============================================================

# --- 3a. PERMANOVA ---
adonis2(spe ~ comm_type, method = "bray")  # 整体检验

# --- 3b. 配对 PERMANOVA ---
# 需要 pairwiseAdonis 包或手动
pairwise_perm <- function(dist_mat, group) {
  groups <- unique(group)
  result <- list()
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      sel <- group %in% c(groups[i], groups[j])
      sub_dist <- as.dist(as.matrix(dist_mat)[sel, sel])
      sub_group <- factor(group[sel])
      mod <- adonis2(sub_dist ~ sub_group)
      result[[paste(groups[i], "vs", groups[j])]] <- mod$`Pr(>F)`[1]
    }
  }
  return(data.frame(comparison = names(result), p.value = unlist(result)))
}
pairwise_perm(bc_dist, comm_type)

# --- 3c. SIMPER ---
simp <- with(spe, simper(spe, comm_type))
summary(simp$TypeA_TypeB)  # TypeA vs TypeB 的关键差异物种

# --- 3d. ANOSIM ---
anosim(bc_dist, comm_type, permutations = 999)

# ============================================================
# C4 — 指示种分析
# ============================================================

# 多水平模式分析
indval <- multipatt(spe, comm_type, func = "r.g",
                    control = how(nperm = 999))
summary(indval)  # 每个群落的指示种及其统计量
