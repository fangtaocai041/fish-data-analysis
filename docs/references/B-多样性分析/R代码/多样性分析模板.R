# ============================================================
# B — 多样性分析 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(vegan)
library(betapart)    # β多样性分解
library(FD)          # 功能多样性
library(picante)     # 系统发育多样性
library(ggplot2)

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(42)
# 10个样方 × 20个物种
spe_mat <- matrix(
  rpois(200, lambda = c(rep(3, 10), rep(15, 5), rep(1, 5))),
  nrow = 10,
  dimnames = list(paste0("Site", 1:10), paste0("sp", 1:20))
)

# 分组信息
group <- factor(c(rep("Control", 4), rep("Treatment", 4), rep("Reference", 2)))

# ============================================================
# B1 — α 多样性
# ============================================================

# --- 1a. 常用指数 ---
alpha <- data.frame(
  site = rownames(spe_mat),
  S = specnumber(spe_mat),
  Shannon = diversity(spe_mat, index = "shannon"),
  Simpson = diversity(spe_mat, index = "simpson"),
  invSimpson = diversity(spe_mat, index = "invsimpson"),
  Pielou = diversity(spe_mat, index = "shannon") / log(specnumber(spe_mat)),
  Fisher = fisher.alpha(spe_mat)
)
alpha

# --- 1b. 稀疏曲线 ---
rarecurve(spe_mat, step = 1, col = group,
          xlab = "样本量", ylab = "物种数",
          main = "稀疏曲线")

# --- 1c. Chao1 估计 ---
estimateR(spe_mat)  # 返回每个样方的 Chao1, ACE 等

# --- 1d. α多样性箱线图 ---
alpha$group <- group
ggplot(alpha, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2) +
  labs(title = "α多样性比较", y = "Shannon index") +
  theme_bw() + theme(legend.position = "none")

# --- 1e. α多样性差异检验 ---
# 两组: t.test()  多组: aov() / kruskal.test()
summary(aov(Shannon ~ group, data = alpha))

# ============================================================
# B2 — β 多样性
# ============================================================

# --- 2a. Bray-Curtis 相异矩阵 ---
bc_dist <- vegdist(spe_mat, method = "bray")

# --- 2b. β多样性分解 (Baselga 2012) ---
# 经典 Jaccard 分解: β_jtu (周转) + β_jne (嵌套)
beta_multi <- beta.multi(spe_mat, index.family = "jaccard")
beta_multi

# 样方对分解
beta_pair <- beta.pair(spe_mat, index.family = "sorensen")
names(beta_pair)

# --- 2c. PERMANOVA（组间 β多样性差异检验）---
adonis2(spe_mat ~ group, method = "bray")

# --- 2d. 同质性检验（PERMANOVA 前提假设）---
mod <- betadisper(bc_dist, group)
boxplot(mod)  # 组内离散度可视化
anova(mod)    # p > 0.05 接受同质性

# --- 2e. NMDS 可视化 ---
nmds <- metaMDS(bc_dist, k = 2, try = 100)
nmds_scores <- as.data.frame(scores(nmds)$sites)
nmds_scores$group <- group

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = group)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.95) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
           label = paste0("Stress: ", round(nmds$stress, 3))) +
  theme_bw()

# ============================================================
# B4 — 功能多样性 (示例)
# ============================================================

# 模拟性状数据（3个性状 × 20物种）
set.seed(123)
traits <- data.frame(
  row.names = paste0("sp", 1:20),
  size = rnorm(20, mean = 10, sd = 3),
  leaf_area = rlnorm(20, meanlog = 2, sdlog = 0.5),
  SLA = rnorm(20, mean = 20, sd = 5)
)
# 性状标准化
traits_scaled <- decostand(traits, method = "standardize")

# 计算功能多样性
fd_out <- dbFD(traits_scaled, spe_mat, w.abun = TRUE)

# 提取指数
fd_indices <- data.frame(
  site = names(fd_out$FRic),
  FRic = fd_out$FRic,
  FEve = fd_out$FEve,
  FDiv = fd_out$FDiv,
  FDis = fd_out$FDis,
  RaoQ = fd_out$RaoQ
)
fd_indices

# ============================================================
# B5 — 系统发育多样性 (示例)
# ============================================================

# 模拟系统发育树
library(ape)
set.seed(456)
phy <- rcoal(20)
phy$tip.label <- paste0("sp", 1:20)

# Faith's PD
pd_out <- pd(spe_mat, phy, include.root = TRUE)
pd_out

# ses.MPD (标准化平均谱系距离)
mpd_out <- ses.mpd(spe_mat, cophenetic(phy),
                   null.model = "taxa.labels",
                   abundance.weighted = TRUE, runs = 999)
mpd_out[, c("ntaxa", "mpd.obs", "mpd.obs.z", "mpd.rand.mean", "runs")]
