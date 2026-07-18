# 📊 鱼类生态数据分析项目

> 基于数量生态学 10 大模块 + 综合工作流 + workspace 分析管线  
> 整合 R/Python 生态分析工具，支持可复现研究

---

## 快速开始

```bash
# 入口命令（在 Reasonix 会话中）
/fish-data-analysis
```

## 项目结构

```
├── data/                          # 数据
│   ├── raw/                       #   原始数据（只读）
│   ├── processed/                 #   清洗后数据
│   └── external/                  #   外部来源数据
│
├── src/                           # 源代码
│   ├── R/                         #   R 分析函数
│   └── python/                    #   Python 分析脚本
│
├── analyses/                      # 可复现分析报告（.Rmd）
│
├── config/                        # 配置文件
│   └── analysis_config.yaml       #   分析参数配置
│
├── docs/                          # 文档
│   └── references/                #   参考资料
│       ├── A-数据基础/            #     10个模块（R代码+深度详解）
│       ├── B-多样性分析/
│       ├── C-多元统计分析/
│       ├── ...
│       ├── 综合工作流/            #     5个完整分析管道
│       ├── 诊断工具箱/            #     综合诊断R函数
│       ├── 方法适用边界/          #     方法选择与约束
│       ├── 快速速查表.md          #     核心阈值一页纸
│       └── workspace脚本/         #     workspace分析管线
│
├── output/                        # 输出产物
│   ├── figures/                   #   图表
│   ├── tables/                    #   表格
│   └── reports/                   #   报告
│
├── tests/                         # 测试
├── scripts/                       # 工具脚本
├── .github/                       # CI/CD 配置
│
├── CHANGELOG.md                   # 变更日志
├── LICENSE                        # MIT 许可证
└── README.md                      # 本文件
```

## 分析流程

### 标准工作流

```mermaid
graph LR
    A[数据准备] --> B[探索性分析]
    B --> C[方法选择]
    C --> D[统计分析]
    D --> E[可视化]
    E --> F[报告输出]
```

### 方法选择路径

```
数据特征 → 决策树（docs/references/Z-决策树/）
  → 梯度长度 < 3 SD → RDA（docs/references/C-多元统计分析/）
  → 梯度长度 3-4 SD → CCA
  → 梯度长度 > 4 SD → CA/DCA
  → 组间比较 → PERMANOVA（含 betadisper 前置检验）
  → 计数数据 → 过离散检测（>1.5 换负二项）
  → 环境变量 → VIF 共线性诊断（>10 剔除）
  → 空间数据 → Moran's I 空间自相关检验
```

## 分析约束速查

| 检查项 | 阈值 | 参考 |
|--------|------|------|
| VIF 共线性 | >10 剔除 | `docs/references/快速速查表.md` |
| 过离散 (dispersion) | >1.5 换负二项 | `docs/references/I-统计建模与机器学习/` |
| 相关系数 \|r\| | >0.7 删一个变量 | `docs/references/快速速查表.md` |
| PERMANOVA 前置 | 必须做 betadisper | `docs/references/C-多元统计分析/` |
| 残差空间自相关 | 必须检验 Moran's I | `docs/references/D-空间生态学/` |
| 样本量 | 每参数≥10观测 | `docs/references/I-统计建模与机器学习/` |
| RDA 样点数 | >变量数×3 | `docs/references/C-多元统计分析/` |

## 软件依赖

- **R 4.6.0**: vegan, ade4, MASS, mgcv, nlme, lme4, sp, raster, terra, spdep, gstat, ggplot2, dplyr, corrplot 等
- **Python 3.13**: pandas, numpy, scipy, matplotlib, seaborn, sklearn, statsmodels
- **MCP 工具**: GBIF, FishBase, PubMed, CNKI, Europe PMC

## 发布规范

- 🗺️ 物种分布图须含九段线
- 📝 出版前申请审图号
- 🐟 已有 shp 数据（D:\ArcGIS\）边界合规
