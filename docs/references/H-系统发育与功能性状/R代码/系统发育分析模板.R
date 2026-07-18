# ============================================================
# H — 系统发育与功能性状 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(ape)          # 系统发育树操作
library(phytools)     # 系统发育比较方法
library(caper)        # PGLS
library(picante)      # 系统发育多样性
library(geiger)       # 系统发育信号 + 模型拟合

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(42)
n_spp <- 30

# 随机生成系统发育树
phy <- rcoal(n_spp)
phy$tip.label <- paste0("sp", 1:n_spp)

# 模拟性状进化（布朗运动）
traits <- data.frame(
  row.names = phy$tip.label,
  body_size = fastBM(phy, sig2 = 1, mu = 10),
  metabolic_rate = fastBM(phy, sig2 = 0.5, mu = 5),
  trophic_level = fastBM(phy, sig2 = 0.3, mu = 3)
)
# 加一些噪声（非系统发育部分）
traits$body_size <- traits$body_size + rnorm(n_spp, 0, 2)
traits$metabolic_rate <- traits$metabolic_rate + rnorm(n_spp, 0, 1)

# 系统发育树可视化
plot(phy, cex = 0.7)
axisPhylo()

# ============================================================
# H1 — 系统发育信号
# ============================================================

# --- 1a. Pagel's λ ---
# 用 phylolm 或 geiger 计算
lambda_body <- phylosig(phy, traits$body_size, method = "lambda")
lambda_body  # λ ≈ 1 = 强信号，≈ 0 = 无信号

lambda_meta <- phylosig(phy, traits$metabolic_rate, method = "lambda")
lambda_meta

# --- 1b. Blomberg's K ---
K_body <- phylosig(phy, traits$body_size, method = "K")
K_body       # K > 1 = 比布朗运动更保守，K < 1 = 趋同

K_meta <- phylosig(phy, traits$metabolic_rate, method = "K")
K_meta

# ============================================================
# H2 — 系统发育比较方法
# ============================================================

# --- 2a. PIC (系统发育独立对比) ---
pic_bodysize <- pic(traits$body_size, phy)
pic_metabolic <- pic(traits$metabolic_rate, phy)

# PIC 回归（强制过原点）
pic_model <- lm(pic_metabolic ~ pic_bodysize - 1)
summary(pic_model)

# --- 2b. PGLS (系统发育广义最小二乘) ---
# 使用 caper 包
comp_data <- comparative.data(phy, traits,
                               names.col = "row.names",
                               vcv = TRUE, na.omit = FALSE)

pgls_mod <- pgls(metabolic_rate ~ body_size, data = comp_data,
                 lambda = "ML")  # λ 的最大似然估计
summary(pgls_mod)

# --- 2c. 使用 phylolm（更灵活） ---
library(phylolm)
pglm_mod <- phylolm(metabolic_rate ~ body_size,
                    data = traits, phy = phy,
                    model = "lambda")  # 布朗运动下的 λ 模型
summary(pglm_mod)

# --- 2d. 比较 OLS vs PGLS ---
ols_mod <- lm(metabolic_rate ~ body_size, data = traits)
summary(ols_mod)

# 比较系数变化
coef(ols_mod)
coef(pglm_mod)

# ============================================================
# H3 — 系统发育多样性
# ============================================================

# 模拟群落矩阵
set.seed(123)
comm <- matrix(
  rbinom(10 * n_spp, 1, 0.3),
  nrow = 10,
  dimnames = list(paste0("Site", 1:10), phy$tip.label)
)
# 转成丰度
comm <- comm * matrix(rpois(10 * n_spp, 5), nrow = 10)

# --- 3a. Faith's PD ---
pd(comm, phy, include.root = TRUE)

# --- 3b. MPD / MNTD ---
# 系统发育距离矩阵
phylo_dist <- cophenetic(phy)

mpd_values <- mpd(comm, phylo_dist, abundance.weighted = TRUE)
mntd_values <- mntd(comm, phylo_dist, abundance.weighted = TRUE)

data.frame(site = rownames(comm), MPD = mpd_values, MNTD = mntd_values)

# --- 3c. ses.MPD (标准化效应大小) ---
ses_mpd <- ses.mpd(comm, phylo_dist,
                   null.model = "taxa.labels",
                   abundance.weighted = TRUE,
                   runs = 999)
ses_mpd[, c("ntaxa", "mpd.obs", "mpd.obs.z", "mpd.rand.mean",
            "mpd.obs.rank", "runs", "mpd.obs.p")]
# mpd.obs.z < -1.96 或 > 1.96 = 显著聚集或分散
# mpd.obs.p < 0.05 = 显著偏离零模型

# --- 3d. NTI / NRI ---
ses_mntd <- ses.mntd(comm, phylo_dist,
                     null.model = "taxa.labels",
                     abundance.weighted = TRUE,
                     runs = 999)
NRI <- -1 * ses_mpd$mpd.obs.z
NTI <- -1 * ses_mntd$mntd.obs.z
data.frame(site = rownames(comm), NRI, NTI)
# NTI > 2 = 系统发育聚集, NTI < -2 = 系统发育分散
