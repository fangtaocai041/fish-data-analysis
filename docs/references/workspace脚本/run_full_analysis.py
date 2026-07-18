#!/usr/bin/env python3
"""
run_full_analysis.py — 物种全景分析编排脚本 (v1.0)

一键执行所有分析管线，适用于任意物种。
自动检测物种是否存在，加载数据，依次运行:
  Phase 0: 物种画像 (KB查询)
  Phase 1: 文献统计 (年份/方向/数量)
  Phase 2: 趋势分析 (年代跃迁)
  Phase 3: 空白识别 (方向/方法/地理)
  Phase 4: 跨物种涌现 (与近亲对比)

不依赖任何物种特定的硬编码字符串，全部从数据实时计算。

用法:
  python scripts/run_full_analysis.py "珠星三块鱼"
  python scripts/run_full_analysis.py "珠星三块鱼" "Tribolodon brandti"
  python scripts/run_full_analysis.py "Pseudaspius hakonensis"
  python scripts/run_full_analysis.py --list-species      # 列出所有可用物种
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

_WORKSPACE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_WORKSPACE))
KB_PATH = _WORKSPACE / "config" / "root_config" / "species_kb.yaml"
GRAPH_PATH = _WORKSPACE / "config" / "root_config" / "species_graph.yaml"


def load_kb_species() -> List[dict]:
    """加载KB中的所有物种"""
    if not KB_PATH.exists():
        return []
    with open(KB_PATH, "r", encoding="utf-8") as f:
        kb = yaml.safe_load(f) or {}
    return kb.get("species", [])


def load_graph_species() -> List[str]:
    """加载图谱中的所有物种"""
    if not GRAPH_PATH.exists():
        return []
    with open(GRAPH_PATH, "r", encoding="utf-8") as f:
        graph = yaml.safe_load(f) or {}
    species_set = set()
    for p in graph.get("graph", {}).get("papers", []):
        sp = p.get("species", [])
        if isinstance(sp, list):
            species_set.update(sp)
        else:
            species_set.add(str(sp))
    return sorted(species_set)


def format_banner(title: str) -> str:
    return f"\n{'='*60}\n{title}\n{'='*60}"


def run_phase_portrait(species_name: str) -> str:
    """Phase 0: 物种画像"""
    lines = [format_banner("Phase 0: 物种画像")]
    for sp in load_kb_species():
        if (species_name.lower() in sp.get("name", "").lower()
                or species_name.lower() in sp.get("scientific", "").lower()
                or species_name.lower() in [a.lower() for a in sp.get("aliases", [])]):
            papers = sp.get("literature", [])
            lines.append(f"  名称: {sp.get('name', '?')}")
            lines.append(f"  学名: {sp.get('scientific', '?')}")
            lines.append(f"  科: {sp.get('family', '?')}")
            lines.append(f"  生态: {sp.get('ecology', '?')}")
            lines.append(f"  论文数: {len(papers)} 篇")
            return "\n".join(lines)
    lines.append(f"  ℹ️ 图谱加载中...")
    return "\n".join(lines)


def run_phase_stats(species_name: str) -> str:
    """Phase 1: 文献统计"""
    lines = [format_banner("Phase 1: 文献统计")]
    try:
        from scripts.cross_synthesis import CrossSynthesis
        engine = CrossSynthesis()
        papers = engine.load_papers(species_name)
    except ImportError:
        papers = []

    if not papers:
        lines.append(f"  ⚠️ 未找到 '{species_name}' 的论文")
        return "\n".join(lines)

    from collections import Counter
    cats = Counter()
    years = []
    journals = Counter()
    for p in papers:
        y = p.get("year", 0)
        if isinstance(y, (int, float)) and y > 0:
            years.append(int(y))
        cats[p.get("category", "others")] += 1
        j = p.get("journal", "")
        if j:
            journals[j] += 1

    lines.append(f"  论文总数: {len(papers)} 篇")
    if years:
        lines.append(f"  时间跨度: {min(years)}-{max(years)}")
    lines.append(f"  研究方向:")
    for cat, n in cats.most_common():
        pct = round(n / len(papers) * 100)
        bar = "█" * max(1, pct // 5) + "░" * max(0, 20 - pct // 5)
        lines.append(f"    {cat:15s} {bar} {n:3d}篇 ({pct:2d}%)")
    if journals:
        top3 = journals.most_common(3)
        lines.append(f"  主要期刊: {' | '.join([f'{j}({n})' for j, n in top3])}")
    return "\n".join(lines)


def run_phase_trend(species_name: str) -> str:
    """Phase 2: 趋势分析"""
    lines = [format_banner("Phase 2: 研究趋势")]
    try:
        from scripts.trend_analyzer import TrendAnalyzer
        ta = TrendAnalyzer()
        result = ta.analyze(species_name, verbose=False)
        if result.get("status") == "ok":
            lines.append(f"  文献跨度: {result['year_range']} ({result['papers']}篇)")
            lines.append(f"  方法学跃迁:")
            for s in result.get("methodology_shifts", []):
                top = " → ".join([f"{m['name']}({m['pct']}%)" for m in s["top_methods"][:3]])
                lines.append(f"    {s['era']:12s} ({s['total_papers']}篇): {top}")
        else:
            lines.append(f"  ⚠️ 趋势分析无结果")
    except ImportError:
        lines.append(f"  ⚠️ trend_analyzer 不可用")
    except Exception as e:
        lines.append(f"  ⚠️ {e}")
    return "\n".join(lines)


def run_phase_gaps(species_name: str) -> str:
    """Phase 3: 空白识别"""
    lines = [format_banner("Phase 3: 研究空白")]
    try:
        from scripts.gap_analyzer import GapAnalyzer
        ga = GapAnalyzer()
        gaps = ga.analyze(species_name)
        dir_gaps = gaps.get("direction_gaps", [])
        if dir_gaps:
            lines.append(f"  研究方向空白:")
            for g in dir_gaps:
                lines.append(f"    {g['severity']} {g['direction']}: 仅 {g['count']} 篇 (期望 {g['expected']})")
        geo = gaps.get("geographic_gaps", [])
        if geo:
            lines.append(f"  地理空白:")
            for g in geo:
                lines.append(f"    {g['severity']} {g['region']}: {g.get('note', '')}")
        method = gaps.get("methodology_gaps", [])
        high = [g for g in method if g.get("potential") == "高"]
        if high:
            lines.append(f"  高潜力未用方法: {' | '.join([g['method'] for g in high[:4]])}")
    except ImportError:
        lines.append(f"  ⚠️ gap_analyzer 不可用")
    except Exception as e:
        lines.append(f"  ⚠️ {e}")
    return "\n".join(lines)


def run_phase_synthesis(sp1: str, sp2: str) -> str:
    """Phase 4: 跨物种涌现"""
    lines = [format_banner(f"Phase 4: 跨物种涌现 ({sp1} ↔ {sp2})")]
    try:
        from scripts.cross_synthesis import CrossSynthesis
        ce = CrossSynthesis()
        result = ce.synthesize(sp1, sp2)
        if result.get("evidence"):
            lines.append(f"  检测器输出 ({len(result['evidence'])} 条):")
            for ev in result["evidence"]:
                icon = "🔬" if ev["confidence"] >= 0.7 else "🧪"
                doi_str = f" (DOI: {', '.join(ev['dois'][:3])})" if ev["dois"] else ""
                lines.append(f"    {icon} [{ev['detector']}] {ev['pattern'][:60]}{doi_str}")
        if result.get("hypotheses"):
            lines.append(f"\n  涌现假说 ({len(result['hypotheses'])} 条):")
            for h in result["hypotheses"]:
                icon = {"高": "🔬", "中": "🧪", "探索": "💭"}.get(h["confidence"], "📌")
                lines.append(f"    {icon} {h['title'][:55]}")
                lines.append(f"       DOI引用: {h['doi_count']}篇")
    except ImportError:
        lines.append(f"  ⚠️ cross_synthesis 不可用")
    except Exception as e:
        lines.append(f"  ⚠️ {e}")
    return "\n".join(lines)


def run_phase_reasoning(sp1: str, sp2: str) -> str:
    """Phase 5: 生态假说推理"""
    lines = [format_banner(f"Phase 5: 生态假说推理 ({sp1} ↔ {sp2})")]
    try:
        from scripts.reasoning_engine import EcologyReasoner
        engine = EcologyReasoner()
        hypotheses = engine.reason(sp1, sp2)

        if hypotheses:
            lines.append(f"  生态假说 ({len(hypotheses)} 条):")
            for h in hypotheses:
                icon = {"高": "🔬", "中": "🧪", "探索": "💭"}.get(h["confidence"], "📌")
                lines.append(f"  {icon} {h['title']} [{h['type']}]")
                lines.append(f"     置信度: {h['confidence']}")
                for step in h.get("reasoning_chain", [])[:3]:
                    lines.append(f"     ├─ {step[:70]}")
                lines.append(f"     └─ 预测: {h.get('prediction', '')[:60]}")
                lines.append(f"        验证: {h.get('test_method', '')[:60]}")
        else:
            lines.append(f"  未产生生态假说（数据不足以支持推理）")
    except Exception as e:
        lines.append(f"  ⚠️ {e}")
    return "\n".join(lines)


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="物种全景分析编排脚本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python scripts/run_full_analysis.py "珠星三块鱼"
  python scripts/run_full_analysis.py "珠星三块鱼" "Tribolodon brandti"
  python scripts/run_full_analysis.py --list-species
        """
    )
    parser.add_argument("species", nargs="?", default="", help="物种名")
    parser.add_argument("compare", nargs="?", default="", help="对比物种")
    parser.add_argument("--list-species", action="store_true", help="列出可用物种")
    parser.add_argument("--json", action="store_true", help="JSON 输出")
    args = parser.parse_args()

    if args.list_species:
        print("=== KB 物种 ===")
        for sp in load_kb_species():
            print(f"  {sp.get('name', '?'):20s} {sp.get('scientific', '?')}")
        print(f"\n=== 图谱物种 ===")
        for s in load_graph_species():
            if s.lower() != "test":
                print(f"  {s}")
        sys.exit(0)

    if not args.species:
        parser.print_help()
        sys.exit(1)

    species = args.species
    compare = args.compare

    # 自动查找近亲物种（图谱中同属）
    if not compare:
        graph_species = [s for s in load_graph_species() if s.lower() != species.lower() and s.lower() != "test"]
        if graph_species:
            compare = graph_species[0]  # 默认取第一个

    # 运行所有阶段
    report = {
        "species": species,
        "compare": compare or "",
        "phases": {}
    }

    report["phases"]["portrait"] = run_phase_portrait(species)
    report["phases"]["stats"] = run_phase_stats(species)
    report["phases"]["trend"] = run_phase_trend(species)
    report["phases"]["gaps"] = run_phase_gaps(species)
    if compare:
        report["phases"]["synthesis"] = run_phase_synthesis(species, compare)
        report["phases"]["reasoning"] = run_phase_reasoning(species, compare)
    else:
        report["phases"]["synthesis"] = "\n(未指定对比物种，跳过涌现)"
        report["phases"]["reasoning"] = ""

    if args.json:
        import json as j
        print(j.dumps(report, ensure_ascii=False, indent=2))
    else:
        for phase_name, phase_output in report["phases"].items():
            print(phase_output)
        print(f"\n{'='*60}")
        print(f"✅ 分析完成: {species}")
        if compare:
            print(f"   对比: {compare}")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
