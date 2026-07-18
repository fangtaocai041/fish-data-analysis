# ============================================================
# 🔧 数量生态学 · 综合诊断工具箱
# 方陶文库 · 覆盖 A-J 十类的核心诊断函数
#
# 使用：source("诊断工具箱.R") 后直接调用各函数
# ============================================================
library(vegan)
library(tidyverse)

# ============================================================
# A — 数据质量一键诊断
# ============================================================

#' 综合数据质量报告
#' @param df data.frame 待检查的数据
#' @param spe_threshold 群落数据稀有物种过滤阈值
diagnose_data_quality <- function(df, spe_threshold = 3) {
  cat("\n========== 数据质量诊断 ==========\n")
  
  # 1. 结构
  cat(sprintf("观测数: %d  变量数: %d\n", nrow(df), ncol(df)))
  
  # 2. 缺失
  miss_pct <- colMeans(is.na(df)) * 100
  high_miss <- miss_pct[miss_pct > 5]
  if (length(high_miss) > 0) {
    cat("⚠️ 高缺失变量 (>5%):\n")
    for (i in seq_along(high_miss)) {
      cat(sprintf("  %s: %.1f%% → %s\n", names(high_miss)[i], high_miss[i],
                  if(high_miss[i] > 50) "建议删除" else "建议MICE插补"))
    }
  } else cat("✅ 缺失率正常 (<5%)\n")
  
  # 3. 数值变量异常值
  num_cols <- names(df)[sapply(df, is.numeric)]
  cat(sprintf("\n数值变量 (%d个):\n", length(num_cols)))
  for (v in num_cols) {
    x <- df[[v]]
    x <- x[is.finite(x)]
    if (length(x) < 5) next
    sk <- abs(moments::skewness(x))
    q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
    iqr <- q[2] - q[1]
    n_out <- sum(x < q[1] - 1.5*iqr | x > q[2] + 1.5*iqr, na.rm = TRUE)
    flag <- if(sk > 2) "🚫 严重偏态" else if(sk > 1) "⚠️ 偏态" else "✅"
    cat(sprintf("  %-15s skew=%.2f %s outliers(IQR)=%d\n", v, sk, flag, n_out))
  }
  
  # 4. 高相关性
  if (length(num_cols) >= 2) {
    cor_mat <- cor(df[num_cols], use = "pairwise.complete.obs")
    high_cor <- which(abs(cor_mat) > 0.7 & upper.tri(cor_mat), arr.ind = TRUE)
    if (nrow(high_cor) > 0) {
      cat(sprintf("\n⚠️ 高相关对 (|r|>0.7): %d对\n", nrow(high_cor)))
    }
  }
  cat("========================================\n")
}

# ============================================================
# B — 采样充分性 + 多样性诊断
# ============================================================

#' 采样充分性诊断
diagnose_sampling <- function(spe) {
  cat("\n========== 采样充分性诊断 ==========\n")
  
  # 稀疏曲线是否接近渐近线
  rare <- rarecurve(spe, plot = FALSE)
  # 检查最后一步的斜率
  for (i in seq_along(rare)) {
    last_step <- diff(tail(rare[[i]], 2))
    flag <- if(abs(last_step) < 0.5) "✅ 接近渐近线" else "⚠️ 仍在上升"
    cat(sprintf("  %s: 最后一步增量=%.1f %s\n", 
                names(rare)[i], last_step, flag))
  }
  
  # Chao1 可靠性
  est <- estimateR(spe)
  cat("\nChao1 可靠性:\n")
  for (i in 1:nrow(spe)) {
    chao <- est["S.chao1", i]
    se <- est["se.chao1", i]
    if (!is.na(se) && chao > 0) {
      cv <- se / chao
      flag <- if(cv > 0.5) "🚫 CV>50%" else if(cv > 0.3) "⚠️" else "✅"
      cat(sprintf("  %s: Chao1=%.0f SE=%.0f CV=%.0f%% %s\n",
                  rownames(spe)[i], chao, se, cv*100, flag))
    }
  }
  cat("========================================\n")
}

# ============================================================
# C — 排序方法选择 + PERMANOVA 完整链路
# ============================================================

#' 排序方法自动决策
decide_ordination <- function(spe_hel, env = NULL) {
  dca <- decorana(spe_hel)
  axis1 <- dca$evals[1]
  
  cat(sprintf("\nDCA 第1轴: %.2f SD → ", axis1))
  if (axis1 < 3) cat("推荐 RDA/PCA (线性)\n")
  else if (axis1 > 4) cat("推荐 CCA/CA (单峰)\n")
  else cat("过渡区 (3-4)，两种都跑\n")
  
  if (!is.null(env)) {
    if (axis1 < 3) {
      mod <- rda(spe_hel ~ ., data = env)
    } else {
      mod <- cca(spe_hel ~ ., data = env)
    }
    vif_vals <- vif.cca(mod)
    cat("VIF检查:\n")
    for (i in seq_along(vif_vals)) {
      cat(sprintf("  %s: %.1f %s\n", names(vif_vals)[i], vif_vals[i],
                  if(vif_vals[i] > 10) "⚠️" else "✅"))
    }
  }
}

#' PERMANOVA 完整链路（含 betadisper + 事后比较）
permanova_full <- function(spe, group, dist_method = "bray", nperm = 9999) {
  dist_mat <- vegdist(spe, method = dist_method)
  
  # Step 1
  bd <- betadisper(dist_mat, group)
  bd_p <- anova(bd)$`Pr(>F)`[1]
  cat(sprintf("\n[1] betadisper: p = %.4f %s\n", bd_p,
              if(bd_p < 0.05) "⚠️ 离散度不同！" else "✅ 齐性"))
  
  # Step 2
  perm <- adonis2(dist_mat ~ group, permutations = nperm)
  cat(sprintf("[2] PERMANOVA: R²=%.3f F=%.2f p=%.4f\n",
              perm$R2[1], perm$F[1], perm$`Pr(>F)`[1]))
  
  # Step 3 (仅多组时)
  groups <- unique(group)
  if (length(groups) > 2) {
    cat("[3] 配对比较 (Bonferroni):\n")
    for (i in 1:(length(groups)-1)) {
      for (j in (i+1):length(groups)) {
        sel <- group %in% c(groups[i], groups[j])
        sub_d <- as.dist(as.matrix(dist_mat)[sel, sel])
        sub_p <- adonis2(sub_d ~ factor(group[sel]), permutations = nperm)$`Pr(>F)`[1]
        adj_p <- min(sub_p * choose(length(groups), 2), 1)
        cat(sprintf("  %s vs %s: p=%.4f (adj: %.4f)\n",
                    groups[i], groups[j], sub_p, adj_p))
      }
    }
  }
}

# ============================================================
# D — 空间自相关诊断
# ============================================================

#' 残差空间自相关检查
check_spatial_autocorr <- function(residuals, coords, k = 8) {
  library(spdep)
  nb <- knn2nb(knearneigh(as.matrix(coords), k = k))
  listw <- nb2listw(nb, style = "W")
  mi <- moran.test(residuals, listw)
  cat(sprintf("\n残差 Moran's I=%.3f p=%.4f → %s\n",
              mi$estimate[1], mi$p.value,
              if(mi$p.value < 0.05) "⚠️ 有显著空间自相关，需校正！" else "✅ 残差空间独立"))
  invisible(mi)
}

# ============================================================
# E — SDM 采样偏差 + AUC 陷阱检查
# ============================================================

#' SDM 快速评估
diagnose_sdm <- function(pred_presence, pred_bg, n_presence) {
  cat("\n========== SDM 诊断 ==========\n")
  cat(sprintf("存在点数: %d → %s\n", n_presence,
              if(n_presence < 20) "🚫 太少" else if(n_presence < 50) "⚠️ 偏少" else "✅"))
  
  # 简单 AUC 近似
  if (requireNamespace("pROC", quietly = TRUE)) {
    library(pROC)
    all_pred <- c(pred_presence, pred_bg)
    all_true <- c(rep(1, length(pred_presence)), rep(0, length(pred_bg)))
    auc_val <- roc(all_true, all_pred)$auc
    cat(sprintf("AUC = %.3f %s\n", auc_val,
                if(auc_val > 0.95) "⚠️ 可能虚高（背景范围过大）"
                else if(auc_val > 0.7) "✅" else "⚠️ 偏低"))
  }
  cat("========================================\n")
}

# ============================================================
# I — GLM 过离散 + 零膨胀诊断
# ============================================================

#' GLM 过离散诊断
diagnose_glm_dispersion <- function(model) {
  if (!inherits(model, "glm")) stop("需要 glm 对象")
  
  disp <- sum(resid(model, type = "pearson")^2) / df.residual(model)
  cat(sprintf("\nDispersion = %.2f → ", disp))
  
  if (disp < 1.2) {
    cat("✅ 无过离散 (Poisson OK)\n")
  } else if (disp < 1.5) {
    cat("⚠️ 轻微过离散 (考虑 quasi-Poisson)\n")
  } else if (disp < 3) {
    cat("🚫 严重过离散！换 NB: glm.nb() 或 glmmTMB(family=nbinom2)\n")
  } else {
    cat("🚫 极度过离散！必须用 NB + 考虑零膨胀\n")
  }
  invisible(disp)
}

#' 零膨胀检查
diagnose_zero_inflation <- function(y) {
  n <- length(y)
  n_zero <- sum(y == 0)
  pct <- n_zero / n * 100
  
  cat(sprintf("零值比例: %.1f%% → ", pct))
  if (pct < 20) cat("✅ 正常\n")
  else if (pct < 40) cat("⚠️ 偏高，检查是否超预期\n")
  else cat("🚫 零膨胀！用 zero-inflated 或 hurdle 模型\n")
  invisible(pct)
}

#' GLMM 随机效应组数检查
diagnose_glmm_groups <- function(n_groups, group_name = "random effect") {
  cat(sprintf("%s 组数: %d → ", group_name, n_groups))
  if (n_groups < 5) cat("🚫 太少！改用固定效应\n")
  else if (n_groups < 10) cat("⚠️ 勉强OK，方差估计可能不稳\n")
  else cat("✅ 足够\n")
}

# ============================================================
# J — 贝叶斯收敛诊断速查
# ============================================================

#' brms 模型快速诊断
diagnose_brms <- function(brms_fit) {
  cat("\n========== 贝叶斯诊断 ==========\n")
  
  # R-hat
  rhats <- brms::rhat(brms_fit)
  max_rhat <- max(rhats, na.rm = TRUE)
  cat(sprintf("Max R-hat: %.4f → %s\n", max_rhat,
              if(max_rhat < 1.01) "✅ 收敛" else "🚫 未收敛！"))
  
  # ESS
  ess <- brms::neff_ratio(brms_fit)
  min_ess <- min(ess, na.rm = TRUE)
  cat(sprintf("Min ESS ratio: %.2f → %s\n", min_ess,
              if(min_ess > 0.1) "✅" else "⚠️"))
  
  cat("========================================\n")
}

# ============================================================
# 一键全诊断（适合快速检查完整分析流程）
# ============================================================

#' 群落数据完整分析前诊断
diagnose_community_full <- function(spe, env = NULL, group = NULL, coords = NULL) {
  cat("\n╔══════════════════════════════════╗")
  cat("\n║  数量生态学 · 一键全诊断       ║")
  cat("\n╚══════════════════════════════════╝\n")
  
  # A: 数据质量
  diagnose_data_quality(as.data.frame(spe))
  
  # B: 采样
  diagnose_sampling(spe)
  
  # C: 排序决策
  spe_hel <- decostand(spe, "hellinger")
  decide_ordination(spe_hel, env)
  
  # C: PERMANOVA (如果有分组)
  if (!is.null(group)) {
    permanova_full(spe, group)
  }
  
  # D: 空间 (如果有坐标)
  if (!is.null(coords) && !is.null(env)) {
    # 检查残差
    rda_mod <- rda(spe_hel ~ ., data = env)
    check_spatial_autocorr(residuals(rda_mod), coords)
  }
  
  cat("\n✅ 全诊断完成\n")
}

cat("✅ 诊断工具箱加载完毕！\n")
cat("可用函数:\n")
cat("  diagnose_data_quality()  — A类数据质量\n")
cat("  diagnose_sampling()      — B类采样充分性\n")
cat("  decide_ordination()      — C类排序决策\n")
cat("  permanova_full()         — C类PERMANOVA链路\n")
cat("  check_spatial_autocorr() — D类空间自相关\n")
cat("  diagnose_sdm()           — E类SDM评估\n")
cat("  diagnose_glm_dispersion()— I类过离散\n")
cat("  diagnose_zero_inflation()— I类零膨胀\n")
cat("  diagnose_glmm_groups()   — I类GLMM组数\n")
cat("  diagnose_brms()          — J类贝叶斯收敛\n")
cat("  diagnose_community_full()— 一键全诊断\n")
