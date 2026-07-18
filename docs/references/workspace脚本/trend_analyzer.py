#!/usr/bin/env python3
"""trend_analyzer.py — 研究趋势分析引擎 (v1.0)

独立可执行脚本 + 可导入模块。
对物种知识库中的文献进行趋势分析，输出:
  §7.1 年份×方向分布矩阵
  §7.2 方法学跃迁路径
  §7.3 核心作者/研究组识别
  §7.4 主题词演变轨迹

对应: LIT_SEARCH_PROTOCOL.md §7 研究趋势分析

用法:
  python scripts/trend_analyzer.py "珠星三块鱼"                      # 交互预览
  python scripts/trend_analyzer.py "珠星三块鱼" --json               # JSON 输出
  python scripts/trend_analyzer.py "珠星三块鱼" --compare "滩头三块鱼" # 双物种对比

  from scripts.trend_analyzer import TrendAnalyzer
  analyzer = TrendAnalyzer()
  result = analyzer.analyze("珠星三块鱼")
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

_WORKSPACE = Path(__file__).resolve().parent.parent
KB_PATH = _WORKSPACE / "config" / "root_config" / "species_kb.yaml"
GRAPH_PATH = _WORKSPACE / "config" / "root_config" / "species_graph.yaml"


# ═══════════════════════════════════════════════════
# 方法学分类（按关键词推断）
# ═══════════════════════════════════════════════════

METHODOLOGY_KEYWORDS = {
    "形态测量": ["morpholog", "morphometry", "鳍条", "鳞式", "可量性状", "pharyngeal teeth", "齿式", "retina", "ganglion"],
    "寄生虫学": ["trematode", "digenea", "parasite", "helminth", "吸虫", "绦虫", "寄生虫", "Metagonimus", "Zoogonidae"],
    "生态观测": ["habitat", "spawning", "trophic", "cascade", "predation", "alarm", "distribution", "feeding", "diet", "产卵", "栖息地", "食性"],
    "生理生化": ["physiology", "spectral", "vitelline", "enzyme", "sperm", "osmoregulation", "激素", "酶活"],
    "分子标记": ["COI", "DNA barcoding", "mitochondrial", "cyt b", "RAPD", "微卫星", "microsatellite", "molecular marker", "PCR"],
    "系统发育": ["phylogeny", "phylogeograph", "speciation", "taxonomic", "systematic", "系统发育", "分类"],
    "基因组学": ["genome", "transcriptome", "genomic", "RNA-seq", "gene expression", "mitogenome", "全基因组", "转录组"],
    "毒理环境": ["heavy metal", "Cs-137", "radiocesium", "Fukushima", "pollution", "重金属", "福岛", "放射性"],
}

# 年代分期
ERAS = [
    (1970, 1989, "奠基期"),
    (1990, 2004, "拓展期"),
    (2005, 2014, "转型期"),
    (2015, 2026, "爆发期"),
]


class TrendAnalyzer:
    """研究趋势分析引擎"""

    def __init__(self, kb_path: Path = KB_PATH, graph_path: Path = GRAPH_PATH):
        self.kb_path = kb_path
        self.graph_path = graph_path

    def load_papers(self, species_name: str) -> List[dict]:
        """从统一数据源加载论文列表"""
        from scripts.kb_loader import get_papers
        return get_papers(species_name)

    def classify_methodology(self, title: str, journal: str = "",
                             category: str = "") -> str:
        """根据标题推断方法学分类"""
        text = (title + " " + journal + " " + category).lower()
        for method, keywords in METHODOLOGY_KEYWORDS.items():
            for kw in keywords:
                if kw.lower() in text:
                    return method
        # fallback to research direction
        dir_map = {
            "genomics": "基因组学",
            "parasitology": "寄生虫学",
            "genetics": "分子标记/系统发育",
            "ecology": "生态观测",
            "morphology": "形态测量",
            "physiology": "生理生化",
            "toxicology": "毒理环境",
        }
        return dir_map.get(category, "未分类")

    def compute_era_distribution(self, papers: List[dict]) -> Dict[str, Any]:
        """计算年份×方向分布矩阵"""
        direction_counter = defaultdict(lambda: Counter())
        era_counter = defaultdict(lambda: defaultdict(int))

        for p in papers:
            year = p.get("year", 0)
            if not isinstance(year, (int, float)) or year == 0:
                continue
            year = int(year)
            cat = p.get("category", "others")
            method = self.classify_methodology(
                p.get("title", ""), p.get("journal", ""), cat
            )

            direction_counter[cat][method] += 1
            for start, end, label in ERAS:
                if start <= year <= end:
                    era_counter[label][cat] += 1
                    break

        return {
            "direction_counter": dict(direction_counter),
            "era_counter": {k: dict(v) for k, v in era_counter.items()},
        }

    def identify_methodology_shift(self, papers: List[dict]) -> List[dict]:
        """识别方法学跃迁路径"""
        timeline = defaultdict(lambda: Counter())

        for p in papers:
            year = p.get("year", 0)
            if not isinstance(year, (int, float)) or year == 0:
                continue
            era = ""
            for start, end, label in ERAS:
                if start <= year <= end:
                    era = label
                    break
            if not era:
                continue
            method = self.classify_methodology(
                p.get("title", ""), p.get("journal", ""), p.get("category", "")
            )
            timeline[era][method] += 1

        # Convert to list format (保持年代顺序)
        era_order = {label: start for start, end, label in ERAS}
        shifts = []
        for era in sorted(timeline.keys(), key=lambda e: era_order.get(e, 9999)):
            methods = timeline[era]
            total = sum(methods.values())
            top = methods.most_common(5)
            shifts.append({
                "era": era,
                "total_papers": total,
                "top_methods": [{"name": m, "count": c, "pct": round(c / total * 100)} for m, c in top],
            })
        return shifts

    def identify_key_authors(self, papers: List[dict]) -> List[dict]:
        """识别核心作者"""
        author_counter = Counter()
        for p in papers:
            authors = p.get("authors", [])
            if authors and isinstance(authors, list):
                for a in authors:
                    if a and isinstance(a, str) and len(a) > 3:
                        author_counter[a] += 1
        return [{"name": n, "papers": c}
                for n, c in author_counter.most_common(15)]

    def analyze(self, species_name: str, verbose: bool = True) -> dict:
        """全量趋势分析"""
        papers = self.load_papers(species_name)
        if not papers:
            return {"status": "not_found", "species": species_name, "papers": 0}

        # 基础统计
        valid_papers = [p for p in papers if isinstance(p.get("year"), (int, float)) and p.get("year", 0) > 0]
        years = sorted(set(int(p["year"]) for p in valid_papers if p.get("year")))
        year_range = f"{years[0]}-{years[-1]}" if years else "N/A"

        # 方向分布
        cat_counter = Counter()
        for p in papers:
            cat_counter[p.get("category", "others")] += 1

        # 分析
        era_dist = self.compute_era_distribution(papers)
        shifts = self.identify_methodology_shift(papers)
        authors = self.identify_key_authors(papers)

        # 趋势叙述
        narrative = self._generate_narrative(papers, shifts, cat_counter)

        result = {
            "status": "ok",
            "species": species_name,
            "papers": len(papers),
            "year_range": year_range,
            "category_distribution": dict(cat_counter.most_common()),
            "era_distribution": era_dist,
            "methodology_shifts": shifts,
            "key_authors": authors,
            "narrative": narrative,
        }
        return result

    def _generate_narrative(self, papers: List[dict],
                            shifts: List[dict],
                            cat_counter: Counter) -> str:
        """生成趋势叙述文本"""
        lines = []

        # 总览
        total = len(papers)
        lines.append(f"📊 共 {total} 篇论文")

        # 方向分布
        lines.append(f"📂 研究方向分布:")
        for cat, count in cat_counter.most_common():
            pct = round(count / total * 100)
            lines.append(f"   {cat:15s} {count:3d}篇 ({pct}%)")

        # 方法学跃迁
        lines.append(f"\n🔄 方法学跃迁路径:")
        for shift in shifts:
            methods_str = " → ".join(
                [f"{m['name']}({m['pct']}%)" for m in shift["top_methods"][:3]]
            )
            lines.append(f"   {shift['era']:10s} ({shift['total_papers']}篇): {methods_str}")

        # 识别趋势信号
        if len(shifts) >= 2:
            latest = shifts[-1]
            prev = shifts[-2]
            latest_methods = {m["name"] for m in latest["top_methods"]}
            prev_methods = {m["name"] for m in prev["top_methods"]}
            emerging = latest_methods - prev_methods
            declining = prev_methods - latest_methods

            if emerging:
                lines.append(f"\n🚀 新兴方向: {', '.join(emerging)}")
            if declining:
                lines.append(f"📉 衰退方向: {', '.join(declining)}")

        return "\n".join(lines)


# ═══════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════

def format_report(result: dict) -> str:
    """格式化为可读报告"""
    if result.get("status") != "ok":
        return f"⚠️ 未找到: {result.get('species', '?')}"

    lines = []
    lines.append("=" * 60)
    lines.append(f"📈 研究趋势分析报告")
    lines.append(f"   物种: {result['species']}")
    lines.append(f"   文献: {result['papers']} 篇 ({result['year_range']})")
    lines.append("=" * 60)
    lines.append("")
    lines.append(result.get("narrative", ""))

    # 核心作者
    authors = result.get("key_authors", [])
    if authors:
        lines.append(f"\n🏆 核心研究者:")
        for i, a in enumerate(authors[:8], 1):
            lines.append(f"   {i:2d}. {a['name']:30s} ({a['papers']}篇)")

    # 年份桶
    lines.append(f"\n📅 年代分布:")
    for shift in result.get("methodology_shifts", []):
        era = shift["era"]
        total = shift["total_papers"]
        top = ", ".join([f"{m['name']}({m['count']})" for m in shift["top_methods"][:3]])
        bar = "█" * min(total, 20) + "░" * max(0, 20 - total)
        lines.append(f"   {era:10s} │{bar}│ {total}篇  {top}")

    return "\n".join(lines)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="研究趋势分析引擎")
    parser.add_argument("species", nargs="?", default="珠星三块鱼", help="物种名")
    parser.add_argument("--json", action="store_true", help="JSON 输出")
    parser.add_argument("--compare", default="", help="对比物种")
    args = parser.parse_args()

    analyzer = TrendAnalyzer()
    result = analyzer.analyze(args.species)

    if args.compare:
        result2 = analyzer.analyze(args.compare)
        if result.get("status") == "ok" and result2.get("status") == "ok":
            print("=" * 60)
            print(f"📊 双物种趋势对比")
            print(f"   {args.species} ({result['papers']}篇) vs {args.compare} ({result2['papers']}篇)")
            print("=" * 60)
            print()
            for label, r in [(args.species, result), (args.compare, result2)]:
                print(f"  [{label}]")
                print(f"  文献范围: {r['year_range']}")
                cats = ", ".join([f"{c}({n})" for c, n in r["category_distribution"].items()])
                print(f"  研究方向: {cats}")
                print()

            # 对比空白
            cats1 = set(result["category_distribution"].keys())
            cats2 = set(result2["category_distribution"].keys())
            only1 = cats1 - cats2
            only2 = cats2 - cats1
            if only1:
                print(f"  🔵 {args.species} 独有方向: {', '.join(only1)}")
            if only2:
                print(f"  🟠 {args.compare} 独有方向: {', '.join(only2)}")
    elif args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print(format_report(result))


if __name__ == "__main__":
    main()
