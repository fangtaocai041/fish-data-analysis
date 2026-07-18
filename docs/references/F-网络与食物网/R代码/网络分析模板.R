# ============================================================
# F — 网络与食物网 R 代码模板
# 方陶文库 · 数量生态学
# ============================================================

library(tidyverse)
library(igraph)       # 网络分析
library(bipartite)    # 二分网络（互惠网络）

# ============================================================
# 0. 模拟数据
# ============================================================

# --- 食物网（邻接矩阵）---
set.seed(42)
n_species <- 15
# 模拟食物网：越高营养级物种索引越大
adj <- matrix(0, nrow = n_species, ncol = n_species)
for (i in 2:n_species) {
  # 随机吃一些低营养级物种
  lower_sp <- 1:(i-1)
  if (length(lower_sp) > 0) {
    prey <- sample(lower_sp, size = sample(1:min(5, length(lower_sp)), 1))
    adj[i, prey] <- 1
  }
}
rownames(adj) <- paste0("sp", 1:n_species)
colnames(adj) <- paste0("sp", 1:n_species)

# ============================================================
# F1 — 食物网分析
# ============================================================

# --- 1a. 转换为 igraph 对象（有向）---
fw_graph <- graph_from_adjacency_matrix(adj, mode = "directed")

# --- 1b. 基本拓扑指标 ---
vcount(fw_graph)       # 节点数 = 物种数
ecount(fw_graph)       # 边数 = 互作数
graph.density(fw_graph)  # 连接度 (Connectance)

# --- 1c. 营养级 ---
# 计算最短路径（从生产者/最底层开始）
dist_matrix <- distances(fw_graph, mode = "in")
# 每个节点的最长入链长度 ≈ 营养级
tl <- apply(dist_matrix, 2, max, na.rm = TRUE)
data.frame(species = V(fw_graph)$name, trophic_level = tl)

# --- 1d. 网络模块性 ---
fw_undirected <- as.undirected(fw_graph, mode = "collapse")
modules <- cluster_louvain(fw_undirected)
modularity(modules)
plot(modules, fw_undirected,
     main = paste("食物网模块，Q =", round(modularity(modules), 3)))

# --- 1e. 关键节点识别 ---
data.frame(
  species = V(fw_graph)$name,
  degree_centrality = degree(fw_graph, mode = "all"),
  betweenness = betweenness(fw_graph, directed = TRUE),
  closeness = closeness(fw_graph, mode = "all")
)

# ============================================================
# F2 — 共现网络（从物种丰度矩阵）
# ============================================================

# --- 2a. 模拟物种丰度 ---
n_sites <- 50
n_spp <- 20
abund <- matrix(rpois(n_sites * n_spp, lambda = 5),
                nrow = n_sites,
                dimnames = list(paste0("Site", 1:n_sites),
                                paste0("sp", 1:n_spp)))

# --- 2b. Spearman 相关矩阵 ---
cor_mat <- cor(abund, method = "spearman")
# 只保留显著相关的对
cor_pval <- apply(cor_mat, 1, function(x) {
  apply(abund, 2, function(y) {
    cor.test(abund[, 1], y, method = "spearman")$p.value
  })
})

# 构建邻接矩阵：|相关系数| > 0.5 且 p < 0.05
adj_cooc <- (abs(cor_mat) > 0.5) * (cor_pval < 0.05)
diag(adj_cooc) <- 0  # 去掉自环

# --- 2c. 构建网络 ---
cooc_graph <- graph_from_adjacency_matrix(
  adj_cooc * cor_mat,  # 保留相关系数作为权重
  mode = "undirected", weighted = TRUE
)

# --- 2d. 网络可视化 ---
V(cooc_graph)$size <- degree(cooc_graph) * 3
plot(cooc_graph,
     main = "物种共现网络",
     vertex.label.cex = 0.7,
     edge.width = abs(E(cooc_graph)$weight) * 2,
     edge.color = ifelse(E(cooc_graph)$weight > 0, "steelblue", "coral"))

# --- 2e. 网络拓扑指标 ---
data.frame(
  n_species = vcount(cooc_graph),
  n_links = ecount(cooc_graph),
  connectance = graph.density(cooc_graph),
  mean_degree = mean(degree(cooc_graph)),
  modularity = modularity(cluster_louvain(cooc_graph)),
  transitivity = transitivity(cooc_graph)  # 聚类系数
)

# ============================================================
# F3 — 二分网络（如传粉者-植物互作）
# ============================================================

# --- 3a. 模拟二分网络 ---
# 行 = 植物，列 = 传粉者
n_plants <- 10
n_pollinators <- 15
interact <- matrix(0, nrow = n_plants, ncol = n_pollinators,
                   dimnames = list(paste0("Plant", 1:n_plants),
                                   paste0("Pol", 1:n_pollinators)))

for (i in 1:n_plants) {
  # 每个植物与部分传粉者互作
  n_partners <- sample(2:5, 1)
  interact[i, sample(1:n_pollinators, n_partners)] <- sample(1:10, n_partners)
}

# --- 3b. 网络级指标 ---
networklevel(interact, index = c("connectance", "nestedness",
                                  "modularity", "H2"))

# --- 3c. 物种级指标 ---
specieslevel(interact)

# --- 3d. 零模型检验（网络嵌套性是否非随机）---
null_m <- nullmodel(interact, N = 100, method = "r00")
obs_nest <- nested(interact, method = "NODF")
null_nest <- sapply(null_m, nested, method = "NODF")
p_value <- sum(null_nest >= obs_nest) / length(null_nest)
cat("观测嵌套性:", obs_nest,
    "\n零模型均值:", mean(null_nest),
    "\nP值:", p_value)

# --- 3e. 网络可视化 ---
plotweb(interact, method = "cca",
        col.high = "darkgreen", col.low = "gold",
        bor.col.high = "darkgreen", bor.col.low = "gold",
        main = "植物-传粉者互作网络")
