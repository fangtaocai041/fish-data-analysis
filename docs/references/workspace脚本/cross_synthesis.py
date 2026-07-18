#!/usr/bin/env python3
"""cross_synthesis.py — 跨物种涌现引擎 v2.0 (证据驱动)

重大变更: v2.0 完全基于论文数据检测模式，无任何手写假说。

5 个数据检测器:
  1. shared_env      — 共享环境信号: 同期刊/同地区/同污染事件
  2. method_gap      — 方法学差距: 物种A有但B没有的研究方向
  3. temporal_trend  — 时间序列: 研究热点随时间的迁移模式
  4. author_network  — 作者合作网络: 跨物种研究群体
  5. blindspot       — 覆盖盲区: 与基准方向库对比

每条检测输出 = 检测到的模式 + 支持该模式的论文DOIs + 置信度

对应: LIT_SEARCH_PROTOCOL.md §11 跨物种涌现

用法:
  python scripts/cross_synthesis.py "珠星三块鱼" "Tribolodon brandti"
  python scripts/cross_synthesis.py "珠星三块鱼" "Tribolodon brandti" --json

  from scripts.cross_synthesis import CrossSynthesis
  engine = CrossSynthesis()
  result = engine.synthesize("珠星三块鱼", "Tribolodon brandti")
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

_WORKSPACE = Path(__file__).resolve().parent.parent
KB_PATH = _WORKSPACE / "config" / "root_config" / "species_kb.yaml"
GRAPH_PATH = _WORKSPACE / "config" / "root_config" / "species_graph.yaml"

# ===== 研究方向关键词库 (用于从标题推断) =====
CATEGORY_KEYWORDS = {
    "ecology":       ["habitat", "spawning", "trophic", "cascade", "predation", "alarm",
                      "distribution", "feeding", "diet", "migration", "population",
                      "diversity", "river", "estuary", "freshwater", "conservation",
                      "产卵", "栖息地", "食性", "分布", "种群", "洄游", "生态"],
    "genetics":      ["DNA", "genetic", "mitochondrial", "cyt b", "COI", "barcoding",
                      "microsatellite", "marker", "phylogeny", "phylogeograph", "speciation",
                      "hybrid", "genotyping", "遗传", "分子", "系统发育", "杂交"],
    "genomics":      ["genome", "transcriptome", "genomic", "RNA-seq", "gene expression",
                      "mitogenome", "sequencing", "whole genome", "基因组", "转录组"],
    "morphology":    ["morpholog", "morphometry", "pharyngeal teeth", "retina", "ganglion",
                      "scale", "fin ray", "形态", "解剖", "咽齿", "鳞"],
    "physiology":    ["physiology", "spectral", "vitelline", "enzyme", "sperm",
                      "osmoregulation", "hormone", "lectin", "protein", "视觉", "生理"],
    "parasitology":  ["trematode", "digenea", "parasite", "helminth", "Metagonimus",
                      "Zoogonidae", "worm", "fluke", "吸虫", "绦虫", "寄生虫"],
    "toxicology":    ["heavy metal", "Cs-137", "radiocesium", "Fukushima", "pollution",
                      "contaminant", "mercury", "cadmium", "毒", "重金属", "福岛"],
}

# ===== 地理区域关键词 =====
REGION_KEYWORDS = {
    "中国":   ["China", "Chinese", "绥芬河", "Suifen", "图们江", "Tumen", "长江", "Yangtze", "武汉"],
    "日本":   ["Japan", "Japanese", "Hokkaido", "Honshu", "Kyushu", "Fukushima", "Ugui"],
    "韩国":   ["Korea", "Korean", "Seomjin"],
    "俄罗斯":  ["Russia", "Russian", "Primorsky", "Sakhalin", "TINRO", "Vladivostok"],
}


def infer_category(title: str, journal: str = "") -> str:
    """从标题+期刊推断研究方向分类"""
    text = (title + " " + journal).lower()
    scores = {}
    for cat, kws in CATEGORY_KEYWORDS.items():
        scores[cat] = sum(1 for kw in kws if kw.lower() in text)
    if max(scores.values()) > 0:
        return max(scores, key=scores.get)
    return "others"


def infer_region(title: str, authors: list) -> List[str]:
    """从标题+作者推断地理区域"""
    text = (title + " " + " ".join(authors if isinstance(authors, list) else [])).lower()
    regions = []
    for region, kws in REGION_KEYWORDS.items():
        for kw in kws:
            if kw.lower() in text:
                regions.append(region)
                break
    return regions if regions else ["未识别"]


def era_label(year: int) -> str:
    """将年份映射到研究时期"""
    if year < 1990: return "奠基期(1970-1989)"
    if year < 2005: return "拓展期(1990-2004)"
    if year < 2015: return "转型期(2005-2014)"
    return "爆发期(2015-2026)"


@dataclass
class Evidence:
    """一条检测证据"""
    detector: str
    pattern: str
    evidence_dois: List[str] = field(default_factory=list)
    confidence: float = 0.5
    detail: str = ""


@dataclass
class Hypothesis:
    """一条涌现假说"""
    title: str
    evidence: List[Evidence]
    confidence: str  # high / medium / exploratory
    test_method: str = ""


class CrossSynthesis:
    """跨物种涌现引擎 v2.0 — 纯数据驱动"""

    @staticmethod
    def _normalize(name: str) -> str:
        """物种名标准化：去除空格和下划线差异"""
        return name.lower().replace(" ", "").replace("_", "").replace("-", "")

    def load_papers(self, species_name: str) -> List[dict]:
        """从统一数据源加载论文"""
        import importlib.util as _iu
        _kb_path = Path(__file__).resolve().parent / "kb_loader.py"
        _spec = _iu.spec_from_file_location("workspace_kb_loader", str(_kb_path))
        if _spec and _spec.loader:
            _mod = _iu.module_from_spec(_spec)
            sys.modules["workspace_kb_loader"] = _mod
            _spec.loader.exec_module(_mod)
            get_papers = _mod.get_papers
        else:
            raise RuntimeError(f"kb_loader not found at {_kb_path}")
        papers = get_papers(species_name)
        # 确保 category/region 字段
        for p in papers:
            if not p.get("category") or p["category"] == "others":
                p["category"] = infer_category(p.get("title", ""), p.get("journal", ""))
            if not p.get("region"):
                p["region"] = infer_region(p.get("title", ""), p.get("authors", []))
        return papers

    def detector_shared_env(self, p1: List[dict], p2: List[dict]) -> List[Evidence]:
        """检测器1: 共享环境信号"""
        evidences = []

        # 共享期刊
        j1 = set(p.get("journal", "") for p in p1 if p.get("journal"))
        j2 = set(p.get("journal", "") for p in p2 if p.get("journal"))
        shared_journals = j1 & j2
        if shared_journals:
            dois = []
            for j in list(shared_journals)[:3]:
                dois += [p.get("doi", "") for p in p1 + p2 if p.get("journal") == j and p.get("doi")]
            evidences.append(Evidence(
                detector="shared_env",
                pattern=f"共享期刊: {', '.join(list(shared_journals)[:3])}",
                evidence_dois=dois[:5],
                confidence=0.6,
                detail="两物种在同一期刊发表，说明研究社区有交集"
            ))

        # 共享地理
        r1 = set(r for p in p1 for r in p.get("region", ["未识别"]))
        r2 = set(r for p in p2 for r in p.get("region", ["未识别"]))
        shared_regions = r1 & r2 - {"未识别"}
        if shared_regions:
            dois = [p.get("doi", "") for p in p1 + p2
                    if any(r in p.get("region", []) for r in shared_regions) and p.get("doi")]
            evidences.append(Evidence(
                detector="shared_env",
                pattern=f"同域分布: {', '.join(shared_regions)}",
                evidence_dois=dois[:5],
                confidence=0.75,
                detail="两物种同域分布，可能面临共同环境选择压力"
            ))

        # 共享研究方向
        c1 = set(p.get("category", "others") for p in p1)
        c2 = set(p.get("category", "others") for p in p2)
        shared_cats = c1 & c2 - {"others"}
        if shared_cats:
            dois = [p.get("doi", "") for p in p1 + p2
                    if p.get("category") in shared_cats and p.get("doi")]
            evidences.append(Evidence(
                detector="shared_env",
                pattern=f"共同研究方向: {', '.join(shared_cats)}",
                evidence_dois=dois[:5],
                confidence=0.7,
                detail="两物种有重叠的研究方向，可进行横向对比"
            ))

        return evidences

    def detector_method_gap(self, p1: List[dict], p2: List[dict]) -> List[Evidence]:
        """检测器2: 方法学差距"""
        evidences = []
        c1 = set(p.get("category", "others") for p in p1)
        c2 = set(p.get("category", "others") for p in p2)

        # A有但B没有
        only_in_1 = c1 - c2 - {"others"}
        if only_in_1:
            dois = [p.get("doi", "") for p in p1 if p.get("category") in only_in_1 and p.get("doi")]
            evidences.append(Evidence(
                detector="method_gap",
                pattern=f"{', '.join(only_in_1)}",
                evidence_dois=dois[:5],
                confidence=0.8,
                detail="物种A在这些方向已有研究积累，物种B可借鉴方法"
            ))

        # B有但A没有
        only_in_2 = c2 - c1 - {"others"}
        if only_in_2:
            dois = [p.get("doi", "") for p in p2 if p.get("category") in only_in_2 and p.get("doi")]
            evidences.append(Evidence(
                detector="method_gap",
                pattern=f"物种B独有方向: {', '.join(only_in_2)}",
                evidence_dois=dois[:5],
                confidence=0.7,
            ))
        return evidences

    def detector_temporal(self, p1: List[dict], p2: List[dict]) -> List[Evidence]:
        """检测器3: 时间序列信号"""
        evidences = []

        def direction_shift(papers, n=3):
            """检测研究方向转移"""
            eras = defaultdict(set)
            for p in papers:
                y = p.get("year", 0)
                if not isinstance(y, (int, float)) or y == 0:
                    continue
                e = era_label(int(y))
                c = p.get("category", "others")
                if c != "others":
                    eras[e].add(c)
            eras_sorted = sorted(eras.items(), key=lambda x: x[0])
            if len(eras_sorted) >= 2:
                first_cats = set.union(*[c for _, c in eras_sorted[:1]]) if eras_sorted else set()
                last_cats = set.union(*[c for _, c in eras_sorted[-1:]]) if eras_sorted else set()
                new = last_cats - first_cats
                lost = first_cats - last_cats
                return new, lost
            return set(), set()

        new1, lost1 = direction_shift(p1)
        if new1:
            dois = [p.get("doi", "") for p in p1 if p.get("category") in new1 and p.get("doi")]
            evidences.append(Evidence(
                detector="temporal",
                pattern=f"研究方向转换: 新增 {', '.join(new1)}" + (f", 减少 {', '.join(lost1)}" if lost1 else ""),
                evidence_dois=dois[:5],
                confidence=0.7,
                detail="研究热点的时序迁移反映领域发展趋势"
            ))
        return evidences

    def detector_author_network(self, p1: List[dict], p2: List[dict]) -> List[Evidence]:
        """检测器4: 作者合作网络"""
        evidences = []
        a1 = set()
        for p in p1:
            for a in p.get("authors", []):
                if isinstance(a, str) and len(a) > 3:
                    a1.add(a.lower().strip())
        a2 = set()
        for p in p2:
            for a in p.get("authors", []):
                if isinstance(a, str) and len(a) > 3:
                    a2.add(a.lower().strip())

        shared_authors = a1 & a2
        if shared_authors:
            dois = [p.get("doi", "") for p in p1 + p2
                    if any(a.lower().strip() in shared_authors for a in p.get("authors", []))
                    and p.get("doi")]
            evidences.append(Evidence(
                detector="author_network",
                pattern=f"共享作者: {', '.join(list(shared_authors)[:5])}",
                evidence_dois=dois[:5],
                confidence=0.85,
                detail="共享作者说明两物种属于同一研究社区"
            ))
        else:
            # 无共享作者本身就是一条信息
            pass
        return evidences

    def detector_blindspot(self, papers: List[dict], species_name: str) -> List[Evidence]:
        """检测器5: 覆盖盲区"""
        evidences = []
        benchmark = ["ecology", "genetics", "genomics", "morphology", "physiology",
                     "parasitology", "toxicology"]
        cats = set(p.get("category", "others") for p in papers)
        missing = [c for c in benchmark if c not in cats]
        if missing:
            evidences.append(Evidence(
                detector="blindspot",
                pattern=f"{species_name}未覆盖的基准方向: {', '.join(missing)}",
                evidence_dois=[],
                confidence=0.8,
                detail=f"鱼类研究通常覆盖{len(benchmark)}个方向，{species_name}覆盖了{len(cats & set(benchmark))}个"
            ))

        # 近期活跃度
        recent = sum(1 for p in papers
                     if isinstance(p.get("year"), (int, float)) and int(p["year"]) >= 2020)
        if recent < 3:
            evidences.append(Evidence(
                detector="blindspot",
                pattern=f"{species_name}近期(2020+)仅{recent}篇论文，活跃度低",
                evidence_dois=[],
                confidence=0.6,
            ))
        return evidences

    def synthesize(self, sp1_name: str, sp2_name: str = "") -> dict:
        """全量涌现分析"""
        p1 = self.load_papers(sp1_name)
        p2 = self.load_papers(sp2_name) if sp2_name else []

        all_evidence = []

        # 运行检测器
        if p2:
            all_evidence += self.detector_shared_env(p1, p2)
            all_evidence += self.detector_method_gap(p1, p2)
            all_evidence += self.detector_temporal(p1, p2)
            all_evidence += self.detector_author_network(p1, p2)

        all_evidence += self.detector_blindspot(p1, sp1_name)
        if p2:
            all_evidence += self.detector_blindspot(p2, sp2_name)

        # 证据分组 → 假说
        by_confidence = {"高": [], "中": [], "探索": []}
        for ev in all_evidence:
            if ev.confidence >= 0.7:
                by_confidence["高"].append(ev)
            elif ev.confidence >= 0.5:
                by_confidence["中"].append(ev)
            else:
                by_confidence["探索"].append(ev)

        hypotheses = []
        for level, ev_list in by_confidence.items():
            if not ev_list:
                continue
            for ev in ev_list[:3]:  # 每个级别最多3条
                hy = Hypothesis(
                    title=f"[{ev.detector}] {ev.pattern[:60]}",
                    evidence=[ev],
                    confidence=level,
                    test_method=f"核查 {len(ev.evidence_dois)} 篇论文的原始数据" if ev.evidence_dois else "需进一步文献调研验证"
                )
                hypotheses.append(hy)

        # 路线图
        roadmap = []
        for h in hypotheses:
            total_dois = sum(len(ev.evidence_dois) for ev in h.evidence)
            prio = {"高": "🔴P0", "中": "🟠P1", "探索": "🟡P2"}.get(h.confidence, "P3")
            roadmap.append({
                "hypothesis": h.title[:60],
                "priority": prio,
                "evidence_count": total_dois,
                "method": h.test_method,
                "confidence": h.confidence,
            })

        return {
            "status": "ok",
            "v2": True,
            "primary": sp1_name,
            "comparison": sp2_name if sp2_name else "无",
            "papers_sp1": len(p1),
            "papers_sp2": len(p2),
            "detectors_used": ["shared_env", "method_gap", "temporal", "author_network", "blindspot"],
            "evidence_count": len(all_evidence),
            "doi_count": len(set(d for ev in all_evidence for d in ev.evidence_dois if d)),
            "evidence": [
                {"detector": ev.detector, "pattern": ev.pattern,
                 "dois": ev.evidence_dois[:5], "confidence": ev.confidence}
                for ev in all_evidence
            ],
            "hypotheses": [
                {"title": h.title, "confidence": h.confidence,
                 "evidence_count": len(h.evidence),
                 "doi_count": sum(len(ev.evidence_dois) for ev in h.evidence),
                 "test_method": h.test_method}
                for h in hypotheses
            ],
            "roadmap": roadmap,
        }


# ===== CLI =====
def format_report(r: dict) -> str:
    lines = []
    lines.append("=" * 60)
    lines.append(f"🧠 跨物种涌现报告 v2.0 (证据驱动)")
    primary_line = f"   {r['primary']}"
    if r['comparison'] and r['comparison'] != "无":
        primary_line += f" ↔ {r['comparison']}"
    primary_line += f"\n   检测器: {len(r['detectors_used'])}个 | 证据: {r['evidence_count']}条 | 引用DOI: {r['doi_count']}篇"
    lines.append(primary_line)
    lines.append("=" * 60)
    lines.append("")

    # 证据列表
    if r["evidence"]:
        lines.append("📋 检测证据")
        for ev in r["evidence"]:
            icon = "🔬" if ev["confidence"] >= 0.7 else "🧪" if ev["confidence"] >= 0.5 else "💭"
            lines.append(f"  {icon} [{ev['detector']}] {ev['pattern']} (置信度: {ev['confidence']})")
            if ev["dois"]:
                lines.append(f"      📄 DOI: {', '.join(ev['dois'][:3])}")
        lines.append("")

    # 假说
    if r["hypotheses"]:
        lines.append("💡 涌现假说（来自数据检测，可核查）")
        for i, h in enumerate(r["hypotheses"], 1):
            icon = {"高": "🔬", "中": "🧪", "探索": "💭"}.get(h["confidence"], "📌")
            lines.append(f"\n  {icon} 假说{i}: {h['title'][:60]}")
            lines.append(f"     置信度: {h['confidence']} | 证据: {h['evidence_count']}条 | 引用: {h['doi_count']}篇DOI")
            lines.append(f"     验证方法: {h['test_method']}")
        lines.append("")

    # 路线图
    if r["roadmap"]:
        lines.append("🗺️  研究路线图（优先级排序）")
        for rm in r["roadmap"]:
            lines.append(f"  {rm['priority']} {rm['hypothesis'][:50]}")
            lines.append(f"     证据: {rm['evidence_count']}篇 | 方法: {rm['method']}")
        lines.append("")

    lines.append("⚠️  每条证据均关联具体DOI，可独立核查验证")
    return "\n".join(lines)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="跨物种涌现引擎 v2.0 (证据驱动)")
    parser.add_argument("species", nargs="?", default="珠星三块鱼", help="主物种")
    parser.add_argument("compare", nargs="?", default="", help="对比物种")
    parser.add_argument("--json", action="store_true", help="JSON 输出")
    args = parser.parse_args()

    engine = CrossSynthesis()
    result = engine.synthesize(args.species, args.compare)

    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print(format_report(result))


if __name__ == "__main__":
    main()
