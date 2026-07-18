#!/usr/bin/env python3
"""gap_analyzer.py — 研究空白识别引擎 (v1.0)

独立可执行脚本 + 可导入模块。
分析知识库论文覆盖度，识别研究盲区和空白，包括:
  §8.1 研究方向空白 — 哪些方向论文数量显著不足
  §8.2 地理空白 — 哪些地区/国家缺乏研究
  §8.3 时间空白 — 哪些年代有断层
  §8.4 方法论空白 — 哪些现代方法尚未应用
  §8.5 近亲对比空白 — 与近亲物种相比缺失了什么

对应: LIT_SEARCH_PROTOCOL.md §8 研究空白识别

用法:
  python scripts/gap_analyzer.py "珠星三块鱼"
  python scripts/gap_analyzer.py "珠星三块鱼" --json
  python scripts/gap_analyzer.py "珠星三块鱼" --compare "滩头三块鱼"

  from scripts.gap_analyzer import GapAnalyzer
  analyzer = GapAnalyzer()
  gaps = analyzer.analyze("珠星三块鱼")
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

_WORKSPACE = Path(__file__).resolve().parent.parent
KB_PATH = _WORKSPACE / "config" / "root_config" / "species_kb.yaml"
GRAPH_PATH = _WORKSPACE / "config" / "root_config" / "species_graph.yaml"

# 基准研究方向 — 鱼类研究应有的覆盖方向
BENCHMARK_DIRECTIONS = {
    "taxonomy": {"label": "🔖 分类学", "expected": ">5篇", "weight": 0.10},
    "morphology": {"label": "🔬 形态学", "expected": ">5篇", "weight": 0.10},
    "genetics_population": {"label": "🧬 种群遗传学", "expected": ">5篇", "weight": 0.15},
    "phylogeography": {"label": "🌍 生物地理学", "expected": ">3篇", "weight": 0.10},
    "genomics": {"label": "🧬 基因组学", "expected": ">3篇", "weight": 0.15},
    "ecology": {"label": "🌿 生态学", "expected": ">10篇", "weight": 0.15},
    "physiology": {"label": "⚡ 生理学", "expected": ">5篇", "weight": 0.10},
    "toxicology": {"label": "☣️ 毒理学", "expected": ">3篇", "weight": 0.05},
    "conservation": {"label": "🛡️ 保护生物学", "expected": ">3篇", "weight": 0.05},
    "aquaculture": {"label": "🐟 养殖与繁育", "expected": ">5篇", "weight": 0.05},
}

# 现代方法基准
MODERN_METHODS = [
    "环境DNA (eDNA)", "全基因组关联分析 (GWAS)", "CRISPR基因编辑",
    "群体基因组学 (Population Genomics)", "单细胞转录组",
    "稳定同位素分析", "耳石微化学", "气候模型耦合",
    "机器学习物种识别", "宏基因组学",
]

# 地区
REGIONS = ["中国", "日本", "韩国", "俄罗斯", "其他"]


class GapAnalyzer:
    """研究空白识别引擎"""

    def __init__(self, kb_path: Path = KB_PATH, graph_path: Path = GRAPH_PATH):
        self.kb_path = kb_path
        self.graph_path = graph_path

    def load_papers(self, species: str) -> List[dict]:
        """从统一数据源加载论文"""
        from scripts.kb_loader import get_papers
        return get_papers(species)

    def load_related_species(self, species: str, related_name: str) -> List[dict]:
        """加载近亲物种的论文"""
        # 先从 KB 加载
        papers = self.load_papers(related_name)
        if papers:
            return papers
        # 再从图谱加载
        if not self.graph_path.exists():
            return []
        with open(self.graph_path, "r", encoding="utf-8") as f:
            graph = yaml.safe_load(f) or {}
        result = []
        for p in graph.get("graph", {}).get("papers", []):
            sp_list = p.get("species", [])
            sp_names = [s.lower() for s in sp_list] if isinstance(sp_list, list) else [str(sp_list).lower()]
            if related_name.lower() in sp_names:
                result.append(p)
        return result

    def analyze_direction_gaps(self, papers: List[dict]) -> List[dict]:
        """研究方向空白分析"""
        cat_counter = Counter()
        for p in papers:
            cat_counter[p.get("category", "others")] += 1

        # 映射 KB 类别到基准方向
        kb_to_benchmark = {
            "genomics": "genomics",
            "genetics": "genetics_population",
            "ecology": "ecology",
            "morphology": "morphology",
            "physiology": "physiology",
            "parasitology": "ecology",  # 寄生虫学归入生态学
            "toxicology": "toxicology",
            "others": "conservation",
        }

        gaps = []
        for key, info in BENCHMARK_DIRECTIONS.items():
            # 找匹配的论文
            matched_cats = [k for k, v in kb_to_benchmark.items() if v == key]
            count = sum(cat_counter.get(c, 0) for c in matched_cats)
            expected_raw = int(info["expected"].replace(">", "").replace("篇", ""))
            if count < expected_raw:
                severity = "🔴" if count == 0 else "🟠" if count < expected_raw * 0.5 else "🟡"
                gaps.append({
                    "direction": info["label"],
                    "key": key,
                    "count": count,
                    "expected": info["expected"],
                    "severity": severity,
                    "severity_label": "缺失" if count == 0 else "严重不足" if count < expected_raw * 0.5 else "不足",
                })
        return gaps

    def analyze_geographic_gaps(self, papers: List[dict],
                                distribution: Optional[dict] = None) -> List[dict]:
        """地理空白分析"""
        # 从标题/作者推断地理来源
        region_keywords = {
            "中国": ["China", "Chinese", "绥芬河", "Suifen River", "图们江", "Tumen River",
                     "长江", "Yangtze", "Jianghan University", "Wuhan"],
            "日本": ["Japan", "Japanese", "Hokkaido", "Honshu", "Ugui", "Watarase",
                     "Fukushima", "Seomjin"],
            "韩国": ["Korea", "Korean", "Seomjin River"],
            "俄罗斯": ["Russia", "Russian", "Primorsky", "Sakhalin", "TINRO",
                       "Far Eastern", "Vladivostok"],
        }

        region_counter = Counter()
        unknown = 0
        for p in papers:
            text = (p.get("title", "") + " " + " ".join(p.get("authors", []))).lower()
            found = False
            for region, kws in region_keywords.items():
                for kw in kws:
                    if kw.lower() in text:
                        region_counter[region] += 1
                        found = True
                        break
                if found:
                    break
            if not found:
                unknown += 1

        gaps = []
        if region_counter.get("中国", 0) < 3:
            gaps.append({
                "region": "🇨🇳 中国",
                "count": region_counter.get("中国", 0),
                "note": "绥芬河是中国唯一分布区，但研究极少" if region_counter.get("中国", 0) < 3 else "",
                "severity": "🔴" if region_counter.get("中国", 0) == 0 else "🟠",
            })
        return gaps

    def analyze_temporal_gaps(self, papers: List[dict]) -> List[dict]:
        """时间空白分析"""
        years = [int(p["year"]) for p in papers
                 if isinstance(p.get("year"), (int, float)) and p["year"] > 0]
        if not years:
            return []

        year_min, year_max = min(years), max(years)
        decade_gaps = []
        for decade_start in range(year_min // 10 * 10, year_max + 1, 10):
            decade_end = decade_start + 9
            count = sum(1 for y in years if decade_start <= y <= decade_end)
            if count == 0:
                decade_gaps.append({
                    "period": f"{decade_start}-{decade_end}",
                    "count": 0,
                    "severity": "🔴",
                    "note": "完全空白",
                })
            elif count < 3 and decade_start >= 1980:
                decade_gaps.append({
                    "period": f"{decade_start}-{decade_end}",
                    "count": count,
                    "severity": "🟡",
                    "note": "论文极少",
                })

        # 最近5年活跃度
        recent = sum(1 for y in years if y >= 2020)
        recent_gap = None
        if recent < 5:
            recent_gap = {
                "period": "2020-2026",
                "count": recent,
                "severity": "🟠",
                "note": "近期活跃度低",
            }

        return decade_gaps + ([recent_gap] if recent_gap else [])

    def analyze_methodology_gaps(self, papers: List[dict]) -> List[dict]:
        """方法论空白分析"""
        text_all = " ".join([
            p.get("title", "") + " " + p.get("journal", "")
            for p in papers
        ]).lower()

        gaps = []
        for method in MODERN_METHODS:
            # 检查该方法是否已被使用
            keywords = method.lower().split()
            found = any(kw.rstrip("(),").strip() in text_all for kw in keywords if len(kw) > 3)
            if not found:
                gaps.append({
                    "method": method,
                    "severity": "🔵",
                    "note": "该方法尚未应用于该物种",
                    "potential": self._estimate_potential(method),
                })
        return gaps

    def _estimate_potential(self, method: str) -> str:
        """估计方法应用潜力"""
        high = ["eDNA", "稳定同位素", "耳石微化学", "群体基因组学"]
        med = ["气候模型", "GWAS", "机器学习"]
        for h in high:
            if h in method:
                return "高"
        for m in med:
            if m in method:
                return "中"
        return "低"

    def compare_with_relative(self, papers: List[dict],
                              relative_papers: List[dict]) -> List[dict]:
        """与近亲物种对比空白"""
        if not relative_papers:
            return []

        self_cats = set(p.get("category", "others") for p in papers)
        rel_cats = set(p.get("category", "others") for p in relative_papers)

        # 近亲有但我们没有的方向
        missing = rel_cats - self_cats
        # 近亲更多但我们也有的方向
        self_count = Counter(p.get("category", "others") for p in papers)
        rel_count = Counter(p.get("category", "others") for p in relative_papers)

        gaps = []
        for cat in missing:
            gaps.append({
                "direction": cat,
                "self_count": 0,
                "relative_count": rel_count[cat],
                "severity": "🔴",
                "note": f"近亲有 {rel_count[cat]} 篇但本物种完全空白",
            })

        for cat in self_cats & rel_cats:
            if self_count[cat] < rel_count[cat] * 0.3 and rel_count[cat] >= 5:
                gaps.append({
                    "direction": cat,
                    "self_count": self_count[cat],
                    "relative_count": rel_count[cat],
                    "severity": "🟠",
                    "note": f"近亲 {rel_count[cat]} 篇 vs 本物种 {self_count[cat]} 篇，差距显著",
                })

        return gaps

    def analyze(self, species: str, relative_species: str = "") -> dict:
        """全量空白分析"""
        papers = self.load_papers(species)
        if not papers:
            return {"status": "not_found", "species": species, "papers": 0}

        direction_gaps = self.analyze_direction_gaps(papers)
        geo_gaps = self.analyze_geographic_gaps(papers)
        temporal_gaps = self.analyze_temporal_gaps(papers)
        method_gaps = self.analyze_methodology_gaps(papers)

        relative_gaps = []
        if relative_species:
            rel_papers = self.load_related_species(species, relative_species)
            if rel_papers:
                relative_gaps = self.compare_with_relative(papers, rel_papers)

        # 汇总
        all_gaps = {
            "direction_gaps": direction_gaps,
            "geographic_gaps": geo_gaps,
            "temporal_gaps": temporal_gaps,
            "methodology_gaps": method_gaps,
            "relative_gaps": relative_gaps,
            "summary": self._generate_summary(direction_gaps, geo_gaps,
                                              temporal_gaps, method_gaps,
                                              relative_gaps),
        }
        return all_gaps

    def _generate_summary(self, *gap_lists) -> str:
        """生成空白摘要"""
        lines = []
        total = sum(len(g) for g in gap_lists)
        lines.append(f"🔍 共识别 {total} 处研究空白")

        severe = sum(1 for gl in gap_lists for g in gl if g.get("severity") == "🔴")
        if severe:
            lines.append(f"🔴 严重空白 {severe} 处 — 完全缺失的方向/年代/方法")
        for gl in gap_lists:
            for g in gl:
                if g.get("severity") == "🔴":
                    title = g.get("direction") or g.get("region") or g.get("period") or g.get("method") or "?"
                    lines.append(f"   🔴 {title}: {g.get('note', '')}")
        for gl in gap_lists:
            for g in gl:
                if g.get("severity") == "🟠" or g.get("severity") == "🟡":
                    title = g.get("direction") or g.get("region") or g.get("period") or g.get("method") or "?"
                    lines.append(f"   {g['severity']} {title}: {g.get('note', '')}")
        return "\n".join(lines)


# ═══════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════

def format_gaps(result: dict) -> str:
    """格式化为可读报告"""
    lines = []
    lines.append("=" * 60)
    lines.append("🔍 研究空白分析报告")
    lines.append("=" * 60)
    lines.append("")
    lines.append(result.get("summary", ""))

    direction = result.get("direction_gaps", [])
    if direction:
        lines.append(f"\n📂 研究方向空白:")
        for g in direction:
            lines.append(f"   {g['severity']} {g['direction']}: 仅 {g['count']} 篇 (期望 {g['expected']})")

    geo = result.get("geographic_gaps", [])
    if geo:
        lines.append(f"\n🌍 地理空白:")
        for g in geo:
            lines.append(f"   {g['severity']} {g['region']}: {g['count']} 篇 — {g.get('note', '')}")

    temporal = result.get("temporal_gaps", [])
    if temporal:
        lines.append(f"\n📅 时间空白:")
        for g in temporal:
            lines.append(f"   {g['severity']} {g['period']}: {g['count']} 篇 — {g.get('note', '')}")

    method = result.get("methodology_gaps", [])
    if method:
        lines.append(f"\n🔬 方法空白:")
        for g in method[:8]:
            potential = g.get('potential', '?')
            lines.append(f"   {g['severity']} {g['method']} (潜力: {potential})")

    relative = result.get("relative_gaps", [])
    if relative:
        lines.append(f"\n🔄 近亲对比空白:")
        for g in relative:
            lines.append(f"   {g['severity']} {g['direction']}: {g.get('note', '')}")

    return "\n".join(lines)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="研究空白识别引擎")
    parser.add_argument("species", nargs="?", default="珠星三块鱼", help="物种名")
    parser.add_argument("--json", action="store_true", help="JSON 输出")
    parser.add_argument("--compare", default="", help="对比近亲物种")
    args = parser.parse_args()

    analyzer = GapAnalyzer()
    result = analyzer.analyze(args.species, args.compare)

    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print(format_gaps(result))


if __name__ == "__main__":
    main()
