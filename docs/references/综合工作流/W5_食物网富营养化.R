# ============================================================
# W5 — 食物网结构对富营养化的响应：完整分析 Pipeline
# 方陶文库 · 数量生态学 · 综合工作流
#
# 串联：F (网络分析) → J1 (零模型检验)
# ============================================================
library(igraph)
library(bipartite)
library(vegan)
library(ggplot2)
library(patchwork)

# ============================================================
# 0. 模拟数据：沿营养梯度的湖泊食物网
# ============================================================
set.seed(555)
n_lakes <- 15
n_spp_per_lake <- 20

# 营养状态梯度
trophic_level <- c(
  seq(10, 30, length.out = 5),   # 贫营养
  seq(35, 65, length.out = 5),   # 中营养
  seq(70, 100, length.out = 5)   # 富营养
)
lakes <- paste0("Lake", 1:n_lakes)

# ============================================================
# F1 — 构建食物网 & 计算拓扑指标
# ============================================================

cat("\n========== F1: 食物网拓扑分析 ==========\n")

# 为每个湖泊生成食物网矩阵并计算指标
net_results <- data.frame(
  lake = lakes,
  trophic = trophic_level,
  S = NA,           # 物种数
  L = NA,           # 连接数
  C = NA,           # connectance
  mean_trophic = NA, # 平均营养级
  max_trophic = NA,  # 最大食物链长度
  modularity = NA,   # 模块性
  nestedness = NA    # WNODF (近似)
)

for (i in 1:n_lakes) {
  tl <- trophic_level[i]
  
  # 富营养化 → 物种数减少、食物链缩短、简化
  n_spp <- max(8, round(20 - tl * 0.12 + rnorm(1, 0, 2)))
  
  # 构建随机食物网矩阵（捕食者×猎物）
  n_pred <- max(3, round(n_spp * 0.3))
  n_prey <- n_spp - n_pred
  
  # 物种间互作概率：富营养化 → 泛化种多 → 连接度高但简单
  connect_prob <- 0.15 + tl * 0.003 + rnorm(1, 0, 0.02)
  connect_prob <- pmin(pmax(connect_prob, 0.1), 0.5)
  
  # 捕食矩阵
  mat <- matrix(
    rbinom(n_pred * n_prey, 1, connect_prob),
    nrow = n_pred, ncol = n_prey
  )
  rownames(mat) <- paste0("P", 1:n_pred)
  colnames(mat) <- paste0("p", 1:n_prey)
  
  # 转为 igraph
  # 边列表：捕食(P) → 猎物(p)
  edge_list <- which(mat == 1, arr.ind = TRUE)
  if (nrow(edge_list) > 0) {
    edges <- data.frame(
      from = rownames(mat)[edge_list[, 1]],
      to = colnames(mat)[edge_list[, 2]]
    )
    # 避免自环 + 确保方向
    g <- graph_from_data_frame(edges, directed = TRUE)
    
    # 计算指标
    net_results$S[i] <- vcount(g)
    net_results$L[i] <- ecount(g)
    net_results$C[i] <- ecount(g) / (vcount(g) * (vcount(g) - 1))  # 有向图 L / (S*(S-1))
    
    # 营养级（简化：使用 PageRank 作为近似）
    net_results$mean_trophic[i] <- 1 + mean(page_rank(g)$vector) * 3  # 缩放
    net_results$max_trophic[i] <- diameter(g, directed = TRUE, unconnected = TRUE) + 1
    
    # 模块性
    comm <- cluster_fast_greedy(as.undirected(g))
    net_results$modularity[i] <- modularity(comm)
    
  } else {
    net_results$S[i] <- n_spp
    net_results$modularity[i] <- 0
  }
}

# ============================================================
# F — 网络指标与营养梯度的关系
# ============================================================

cat("\n========== F: 网络结构沿营养梯度 ==========\n")

# 相关性分析
cat("网络指标 vs 营养梯度的 Spearman 相关:\n")
metrics <- c("S", "C", "mean_trophic", "max_trophic", "modularity")
for (m in metrics) {
  r <- cor(net_results[[m]], net_results$trophic, 
           method = "spearman", use = "complete.obs")
  cat(sprintf("  %-20s: ρ = %+.3f\n", m, r))
}

# ============================================================
# J1 — 零模型检验
# ============================================================

cat("\n========== J1: 零模型检验 ==========\n")

# 选贫营养湖泊 vs 富营养湖泊的比较
oligo_lakes <- which(trophic_level <= 30)
eutro_lakes <- which(trophic_level >= 70)

if (length(oligo_lakes) >= 2 && length(eutro_lakes) >= 2) {
  # 比较 connectance
  C_oligo <- net_results$C[oligo_lakes]
  C_eutro <- net_results$C[eutro_lakes]
  
  obs_diff <- mean(C_eutro, na.rm = TRUE) - mean(C_oligo, na.rm = TRUE)
  
  # 零模型：随机重排营养状态标签
  n_perm <- 9999
  null_diffs <- numeric(n_perm)
  set.seed(42)
  for (p in 1:n_perm) {
    shuffled <- sample(n_lakes)
    null_oligo <- net_results$C[shuffled[1:length(oligo_lakes)]]
    null_eutro <- net_results$C[shuffled[(length(oligo_lakes)+1):(length(oligo_lakes)+length(eutro_lakes))]]
    null_diffs[p] <- mean(null_eutro, na.rm = TRUE) - mean(null_oligo, na.rm = TRUE)
  }
  
  p_value <- (sum(abs(null_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
  
  cat(sprintf("贫营养湖泊 mean C: %.4f\n", mean(C_oligo, na.rm = TRUE)))
  cat(sprintf("富营养湖泊 mean C: %.4f\n", mean(C_eutro, na.rm = TRUE)))
  cat(sprintf("观测差异: %+.4f\n", obs_diff))
  cat(sprintf("置换检验 p = %.4f (n=%d)\n", p_value, n_perm))
  
  if (p_value < 0.05) {
    cat("→ ✅ 富营养化显著改变了食物网连接度\n")
  } else {
    cat("→ 未检测到显著差异（可能样本量太小）\n")
  }
}

# ============================================================
# 可视化
# ============================================================

# P1: 各指标沿营养梯度的变化
p1 <- ggplot(net_results, aes(x = trophic, y = C)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(title = "Connectance 沿营养梯度",
       x = "营养状态指数", y = "Connectance") +
  theme_bw()

p2 <- ggplot(net_results, aes(x = trophic, y = S)) +
  geom_point(size = 3, color = "coral") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(title = "物种数 沿营养梯度",
       x = "营养状态指数", y = "物种数 (S)") +
  theme_bw()

p3 <- ggplot(net_results, aes(x = trophic, y = max_trophic)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(title = "食物链长度 沿营养梯度",
       x = "营养状态指数", y = "最大营养级") +
  theme_bw()

p4 <- ggplot(net_results, aes(x = trophic, y = modularity)) +
  geom_point(size = 3, color = "purple") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(title = "模块性 沿营养梯度",
       x = "营养状态指数", y = "Modularity (Q)") +
  theme_bw()

# 合并图
final_fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "W5: 食物网结构对富营养化的响应",
    caption = "方陶文库 · 数量生态学"
  )

ggsave("Fig_W5_食物网富营养化.pdf", final_fig, width = 12, height = 10)
cat("\n✅ W5 分析完成！\n")
