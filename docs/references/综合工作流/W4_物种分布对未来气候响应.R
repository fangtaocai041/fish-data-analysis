# ============================================================
# W4 — 物种分布对未来气候的响应：完整分析 Pipeline
# 方陶文库 · 数量生态学 · 综合工作流
#
# 串联：E (SDM) → I5 (集成预测) → J (不确定性)
# ============================================================
library(dismo)
library(maxnet)
library(raster)
library(ENMeval)
library(randomForest)
library(gbm)
library(ggplot2)
library(patchwork)
library(viridis)

# ============================================================
# 0. 模拟数据
# ============================================================
set.seed(2024)
n_occurrence <- 80   # 存在点
n_background <- 500  # 背景点

# 模拟中国中部某区域环境栅格（简化为数据框）
n_grid <- 2000
env_grid <- data.frame(
  bio1 = rnorm(n_grid, mean = 18, sd = 5),      # 年均温
  bio12 = rlnorm(n_grid, meanlog = 7, sdlog = 0.5), # 年降水
  bio5 = rnorm(n_grid, mean = 32, sd = 4),       # 最热月最高温
  elev = rnorm(n_grid, mean = 500, sd = 300)     # 海拔
)

# 模拟当前气候下的存在点（偏好中等温度+高降水）
env_occurrence <- data.frame(
  bio1 = rnorm(n_occurrence, mean = 18, sd = 2),
  bio12 = rlnorm(n_occurrence, meanlog = 7.2, sdlog = 0.3),
  bio5 = rnorm(n_occurrence, mean = 30, sd = 3),
  elev = rnorm(n_occurrence, mean = 400, sd = 150)
)

# 模拟未来气候（温度上升2°C，降水减少10%）
env_future <- data.frame(
  bio1 = env_grid$bio1 + 2 + rnorm(n_grid, 0, 0.5),
  bio12 = env_grid$bio12 * 0.9 * rlnorm(n_grid, 0, 0.1),
  bio5 = env_grid$bio5 + 3 + rnorm(n_grid, 0, 0.5),
  elev = env_grid$elev  # 海拔不变
)

# ============================================================
# E — 物种分布模型
# ============================================================

cat("\n========== E: 物种分布模型 ==========\n")

# 检查变量相关性
env_cor <- cor(env_grid)
cat("环境变量相关性:\n")
print(round(env_cor, 2))
high_cor <- which(abs(env_cor) > 0.8 & upper.tri(env_cor), arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("⚠️ |r| > 0.8 的变量对，建议保留一个\n")
}

# E2a: MaxEnt 模型
cat("\n--- MaxEnt ---\n")
maxent_mod <- maxnet(
  p = c(rep(1, n_occurrence), rep(0, n_background)),
  data = rbind(env_occurrence, env_grid[sample(1:n_grid, n_background), ]),
  maxnet.formula(p ~ bio1 + bio12 + bio5 + elev,
                 classes = "lqph"),
  regmult = 1.0
)
cat("MaxEnt 拟合完成\n")

# 预测当前分布
pred_current <- predict(maxent_mod, env_grid, type = "cloglog")
# 预测未来分布
pred_future  <- predict(maxent_mod, env_future, type = "cloglog")

cat(sprintf("当前适宜生境比例 (>0.5): %.1f%%\n",
            mean(pred_current > 0.5) * 100))
cat(sprintf("未来适宜生境比例 (>0.5): %.1f%%\n",
            mean(pred_future > 0.5) * 100))

# E2b: 随机森林（集成对比）
cat("\n--- 随机森林 ---\n")
train_data <- rbind(
  data.frame(presence = 1, env_occurrence),
  data.frame(presence = 0, env_grid[sample(1:n_grid, n_occurrence * 2), ])
)
rf_mod <- randomForest(
  as.factor(presence) ~ bio1 + bio12 + bio5 + elev,
  data = train_data, ntree = 500
)
cat(sprintf("RF OOB error: %.1f%%\n", rf_mod$err.rate[500, 1] * 100))

rf_current <- predict(rf_mod, env_grid, type = "prob")[, 2]
rf_future  <- predict(rf_mod, env_future, type = "prob")[, 2]

# ============================================================
# E3 — 模型评估
# ============================================================

cat("\n========== E3: 模型评估 ==========\n")

# 用交叉验证评估（简化版）
# AUC (简单模拟)
set.seed(999)
cv_pred <- predict(maxent_mod, env_occurrence, type = "cloglog")
bg_pred  <- predict(maxent_mod, env_grid[sample(1:n_grid, n_occurrence), ], type = "cloglog")

# 简单的AUC近似
all_pred <- c(cv_pred, bg_pred)
all_true <- c(rep(1, n_occurrence), rep(0, n_occurrence))

if (requireNamespace("pROC", quietly = TRUE)) {
  library(pROC)
  auc_val <- roc(all_true, all_pred)$auc
  cat(sprintf("AUC (approx): %.3f\n", auc_val))
} else {
  cat("AUC: 需要 pROC 包\n")
}

# ============================================================
# I5/J — 集成预测 + 不确定性
# ============================================================

cat("\n========== I5/J: 集成与不确定性 ==========\n")

# 多模型集合均值
ensemble_current <- (pred_current + rf_current) / 2
ensemble_future  <- (pred_future + rf_future) / 2

# 模型间的一致性（标准差越小 = 模型越一致）
uncertainty_current <- apply(cbind(pred_current, rf_current), 1, sd)
uncertainty_future  <- apply(cbind(pred_future, rf_future), 1, sd)

cat(sprintf("当前预测不确定性 (SD 均值): %.3f\n", mean(uncertainty_current)))
cat(sprintf("未来预测不确定性 (SD 均值): %.3f\n", mean(uncertainty_future)))

# 分布变化分析
# 类别：稳定、扩张、收缩、新适宜
threshold <- 0.5
current_binary <- ensemble_current > threshold
future_binary  <- ensemble_future > threshold

change <- ifelse(current_binary & future_binary, "稳定",
          ifelse(!current_binary & future_binary, "扩张",
          ifelse(current_binary & !future_binary, "收缩", "不适宜")))

change_summary <- table(change)
cat("\n未来分布变化:\n")
for (nm in names(change_summary)) {
  cat(sprintf("  %s: %d (%.1f%%)\n", nm, change_summary[nm],
              change_summary[nm] / n_grid * 100))
}

# ============================================================
# 可视化
# ============================================================

# P1: 当前 vs 未来适宜度分布（密度图）
pred_df <- rbind(
  data.frame(period = "当前", suitability = ensemble_current),
  data.frame(period = "未来", suitability = ensemble_future)
)
p1 <- ggplot(pred_df, aes(x = suitability, fill = period)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("当前" = "#2E86AB", "未来" = "#F18F01")) +
  labs(title = "生境适宜度分布：当前 vs 未来",
       x = "适宜度", y = "密度") +
  theme_bw()

# P2: 变化类别
change_df <- data.frame(
  category = names(change_summary),
  count = as.numeric(change_summary),
  pct = as.numeric(change_summary) / n_grid * 100
)
p2 <- ggplot(change_df, aes(x = reorder(category, -pct), y = pct, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("稳定" = "#4575B4", "扩张" = "#1A9641",
                                "收缩" = "#D73027", "不适宜" = "grey80")) +
  labs(title = "分布变化类别", x = "", y = "面积比例 (%)") +
  theme_bw() + theme(legend.position = "none")

# P3: 不确定性地图（简化散点图）
uncertainty_df <- data.frame(
  bio1 = env_grid$bio1,
  bio12 = env_grid$bio12,
  uncertainty = uncertainty_future
)
p3 <- ggplot(uncertainty_df, aes(x = bio1, y = bio12, color = uncertainty)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_viridis(option = "inferno") +
  labs(title = "未来预测不确定性",
       subtitle = paste("均值 SD =", round(mean(uncertainty_future), 3)),
       x = "年均温 (°C)", y = "年降水 (log)") +
  theme_bw()

# P4: 变量重要性（MaxEnt）
var_imp <- data.frame(
  variable = c("bio1", "bio12", "bio5", "elev"),
  importance = abs(coef(maxent_mod)[-1])  # 排除截距
)
var_imp <- var_imp[order(var_imp$importance, decreasing = TRUE), ]
p4 <- ggplot(var_imp, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "MaxEnt 变量重要性", x = "", y = "|系数|") +
  theme_bw()

# 合并图
final_fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "W4: 物种分布对未来气候的响应",
    caption = "方陶文库 · 数量生态学"
  )

ggsave("Fig_W4_物种分布气候响应.pdf", final_fig, width = 14, height = 10)
cat("\n✅ W4 分析完成！\n")
