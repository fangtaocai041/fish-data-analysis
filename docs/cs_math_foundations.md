# 🖥️ 计算机科学与数学基础（生态数据分析）

> ⚠️ **来源说明**：本文档列出的是生态数据分析中常用的计算机和数学方法。每条标注了来源等级：
> - ✅ **已核实**：有具体文献/教材可查
> - 📚 **教材公认**：多本标准教材覆盖
> - 🧪 **前沿进展**：近年论文提出，使用者需自行验证

---

## 一、数学基础

### 1.1 线性代数（生态学应用）

| 方法 | 生态应用 | 来源 |
|------|---------|------|
| **矩阵分解（SVD/PCA）** | 排序分析的核心数学基础 | 📚 Legendre & Legendre (2012) *Numerical Ecology* |
| **特征值分解** | RDA/CCA 的数学本质 | 📚 ter Braak (1986) *Ecology* |
| **Cholesky分解** | 空间相关结构建模（`spdep`） | 📚 Pinheiro & Bates (2000) *Mixed-Effects Models* |
| **非负矩阵分解（NMF）** | 群落结构解析、元条形码 | 🧪 Ren et al. (2020) *Methods Ecol Evol* |

### 1.2 概率论与统计

| 概念 | 生态应用 | 来源 |
|------|---------|------|
| **贝叶斯定理** | brms/RStan 贝叶斯建模 | ✅ McElreath (2020) *Statistical Rethinking* |
| **马尔可夫链蒙特卡洛（MCMC）** | 复杂模型参数估计 | ✅ Gelman et al. (2013) *Bayesian Data Analysis* |
| **隐马尔可夫模型（HMM）** | 动物运动轨迹建模（`moveHMM`） | ✅ Langrock et al. (2012) *JABES* |
| **点过程（Poisson/Cox）** | 物种分布点格局 | 📚 Illian et al. (2008) *Statistical Analysis of Spatial Point Patterns* |
| **极值理论** | 极端气候对种群的影响 | 📚 Coles (2001) *Extremes* |

### 1.3 信息论

| 概念 | 生态应用 | 来源 |
|------|---------|------|
| **香农熵 → Hill numbers** | α多样性统一框架 | ✅ Jost (2006) *Oikos* |
| **互信息** | 生态网络中的关联检测 | ✅ Kinney & Atwal (2014) *PNAS* |
| **转移熵** | 生态时间序列因果推断 | 🧪 Runge et al. (2019) *Science Advances* |
| **AIC / BIC / WAIC** | 模型选择的信息准则 | ✅ Burnham & Anderson (2002) *Model Selection* |

### 1.4 图论与网络

| 概念 | 生态应用 | 来源 |
|------|---------|------|
| **有向无环图（DAG）** | 因果推断、结构方程模型 | ✅ Pearl (2009) *Causality* |
| **网络模块度** | 食物网/互惠网络模块划分 | ✅ Newman (2010) *Networks* |
| **随机块模型（SBM）** | 生态网络聚类 | 🧪 Sinke et al. (2023) *Ecol Complex* |
| **图神经网络（GNN）** | 生态网络链接预测 | 🧪 Strydom et al. (2023) *Ecol Lett* |
| **拉普拉斯矩阵** | 景观连通度分析 | ✅ Goddard (2011) *J Math Biol* |

### 1.5 动力系统

| 概念 | 生态应用 | 来源 |
|------|---------|------|
| **逻辑斯蒂方程** | 种群增长模型基础 | ✅ Gotelli (2008) *Primer of Ecology* |
| **Lotka-Volterra 方程** | 捕食-竞争模型 | ✅ Kingsland (1995) *Modeling Nature* |
| **延迟微分方程** | 种群延迟密度依赖 | 📚 Murdoch et al. (2003) *Consumer-Resource* |
| **随机微分方程** | 环境随机性下的种群动态 | ✅ Lande et al. (2003) *Stochastic Population Dynamics* |
| **分岔理论** | 生态稳态转换、临界点 | ✅ Scheffer (2009) *Critical Transitions* |

---

## 二、计算机科学方法

### 2.1 机器学习（生态学常用）

| 方法 | 生态应用 | 已装R包 | 来源 |
|------|---------|---------|------|
| **随机森林** | 物种分布建模、变量重要性 | `randomForest` | ✅ Breiman (2001) |
| **梯度提升（GBM/BRT）** | 物种对环境响应曲线 | `gbm`/`dismo` | ✅ Elith et al. (2008) *J Anim Ecol* |
| **支持向量机（SVM）** | 遥感分类、生境分类 | `e1071` | ✅ Mountrakis et al. (2011) *ISPRS* |
| **k-近邻（kNN）** | 缺失数据插补 | `class` | 📚 Hastie et al. (2009) *ESL* |
| **K-means / DBSCAN** | 群落聚类、功能群划分 | `stats` | ✅ Legendre & Legendre (2012) |
| **主成分分析（PCA）** | 环境变量降维 | `stats` | ✅ Legendre & Legendre (2012) |

### 2.2 深度学习（生态学前沿）

| 方法 | 生态应用 | 来源 |
|------|---------|------|
| **卷积神经网络（CNN）** | 鱼类图像识别、遥感分类 | ✅ Villon et al. (2018) *Ecol Inform* |
| **长短期记忆（LSTM）** | 生态时间序列预测 | 🧪 Yang et al. (2023) *Ecol Model* |
| **自编码器** | 生态数据降维与异常检测 | 🧪 Chen et al. (2021) *Methods Ecol Evol* |
| **图神经网络（GNN）** | 生态互作网络预测 | 🧪 Strydom et al. (2023) *Ecol Lett* |
| **Transformer** | 多物种联合分布预测 | 🧪 Pichler & Hartig (2023) *Ecol Lett* |

> ⚠️ **注意**：深度学习方法在生态学中仍处于前沿探索阶段，论文发表前需仔细验证。

### 2.3 优化算法

| 方法 | 生态应用 | 来源 |
|------|---------|------|
| **模拟退火** | 保护区规划（`Marxan`） | ✅ Possingham et al. (2000) |
| **遗传算法** | IBM参数校准、模式导向建模 | ✅ Grimm & Railsback (2005) |
| **贝叶斯优化** | SDM超参数调优 | 🧪 Snelson et al. (2023) |

### 2.4 数据处理与算法

| 方法 | 生态应用 | 来源 |
|------|---------|------|
| **动态时间规整（DTW）** | 生态时间序列相似性比较 | 📚 Giorgino (2009) *J Stat Soft* |
| **频繁模式挖掘** | 物种关联规则发现 | 🧪 O'Dwyer et al. (2023) |
| **降维（t-SNE / UMAP）** | 高通量生态数据可视化 | 🧪 Diaz-Papkovich et al. (2024) *Ecol Evol* |

---

## 三、交叉学科工具

| 工具/方法 | 来源学科 | 生态应用 |
|-----------|---------|---------|
| **稳定同位素混合模型（MixSIAR）** | 化学→生态 | 食性分析 |
| **eDNA 分析管道** | 分子→生态 | 物种检测 |
| **生态声学分析** | 信号处理→生态 | 生物多样性监测 |
| **遥感植被指数（NDVI/EVI）** | 遥感→生态 | 生境质量评估 |
| **雷达测速（鸟类/昆虫迁移）** | 物理→生态 | 迁徙研究 |
| **公民科学数据（eBird/iNaturalist）** | 计算机→生态 | 大规模分布数据 |

---

## 四、建议学习路径

```
数学基础（线性代数+概率论）
    ↓
R/Python 编程
    ↓
经典统计（GLM、混合模型）
    ↓
多元分析（vegan/PCA/RDA/NMDS）
    ↓
空间统计（spdep/gstat/INLA）
    ↓
机器学习（RF/GBM/SVM）
    ↓
深度学习（视需求）
    ↓
因果推断（DAG/工具变量）
```

---

> 📌 **使用原则**：优先使用经典统计方法（可解释性强），ML/DL 作为补充。
> 论文中使用任何方法都需引用原始方法文献，并说明参数选择依据。
