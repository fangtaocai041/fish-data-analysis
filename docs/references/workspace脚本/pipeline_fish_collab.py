#!/usr/bin/env python3
"""
pipeline_fish_collab.py — f项目自主运行 → 多项目协作管线

调用链:
  Phase 1: fish-ecology-assistant (自主运行) → 物种知识库查询
  Phase 2: cognitive-search-engine      → 文献搜索
  Phase 3: conflict_verdict (来自 lookup_species, 不重复计算)
  Phase 4: porpoise-agent / coilia-agent → 领域评估
  Phase 5: 全栈汇总

用法:
  from pipeline_fish_collab import run_fish_pipeline
  result = run_fish_pipeline("鳤")
  print(result["summary"])

或者命令行:
  python workspace/pipeline_fish_collab.py 鳤
"""

from __future__ import annotations

import logging
import os
import sys
import time
from typing import Any, Dict, List, Optional

# Force UTF-8 for stdout to prevent GBK encoding errors with Chinese text
os.environ["PYTHONIOENCODING"] = "utf-8"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger("pipeline")


def run_fish_pipeline(
    species_name: str,
    enable_literature: bool = True,
    enable_conflict: bool = True,
    enable_assessment: bool = False,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    f项目自主运行 → 多项目协作管线。

    Args:
        species_name: 物种中文名或学名
        enable_literature: 是否启用文献搜索协同 (Phase 2)
        enable_conflict: 是否启用冲突仲裁 (Phase 3)
        enable_assessment: 是否启用领域评估 (Phase 4)
        verbose: 是否输出阶段日志

    Returns:
        结构化管线结果 dict
    """
    _log = logger.info if verbose else lambda _: None

    # 确定工作区根目录 (D:\Reasonix)
    _WORKSPACE_ROOT = _get_workspace_root()

    _log("\n" + "=" * 60)
    _log(f"  🔥 管线启动 — {species_name}")
    _log("=" * 60)

    # ── Phase 0: 导入依赖 ──
    _log("\n▸ Phase 0/5: 初始化项目适配器")
    try:
        from workspace import (
            lookup_species, search_species, assess_conflict,
            assess_conservation, assess_species,
        )
        _log("  ✓ 适配器就绪")
    except ImportError as exc:
        return {"status": "error", "error": f"适配器导入失败: {exc}"}

    result: Dict[str, Any] = {
        "species": species_name,
        "phases": {},
        "status": "ok",
        "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    # ══════════════════════════════════════════
    # Phase 1: f项目自主运行 — 知识库查询
    # ══════════════════════════════════════════
    _log(f"\n▸ Phase 1/5: {_phase_tag('1', 'fish-ecology-assistant')} 自主运行")
    _log("  └─ 物种知识库查询...")

    phase1_start = time.time()
    try:
        profile = lookup_species(species_name)
        result["phases"]["1_fish"] = profile
        elapsed = time.time() - phase1_start

        if profile.get("known_species"):
            sd = profile.get("species_data", {}) or {}
            _log(f"     ✓ 命中知识库")
            _log(f"       学名:     {profile.get('scientific_name', '?')}")
            _log(f"       中文名:   {profile.get('chinese_name', '?')}")
            _log(f"       科:       {sd.get('family', sd.get('family', '?'))}")
            if sd.get("distribution"):
                dist = sd["distribution"]
                _log(f"       大陆:     {', '.join(dist.get('continents', []))}")
                _log(f"       国家:     {', '.join(dist.get('countries', []))}")
                _log(f"       流域:     {', '.join(dist.get('basins', []))}")
            _log(f"     ⏱ {elapsed:.2f}s")
        else:
            _log(f"     ⚠️ 知识库未收录此物种")
            _log(f"     ⏱ {elapsed:.2f}s")
    except Exception as exc:
        _log(f"     ❌ 查询失败: {exc}")
        result["phases"]["1_fish"] = {"status": "error", "error": str(exc)}

    # ══════════════════════════════════════════
    # Phase 2: 协同 — 文献搜索 + 三角验证评分
    # ══════════════════════════════════════════
    if enable_literature:
        _log(f"\n▸ Phase 2/5: {_phase_tag('2', 'cognitive-search-engine')} 协同")
        _log("  └─ 文献搜索 (lit-search v3.1 + credibility_score)...")

        phase2_start = time.time()
        try:
            # Read species_graph directly (0 token)
            import yaml
            _graph_path = _WORKSPACE_ROOT / "cognitive-search-engine" / "config" / "species_graph.yaml"
            if _graph_path.exists():
                with open(_graph_path, encoding="utf-8") as f:
                    _graph = yaml.safe_load(f)

                # Find species
                _q = species_name.lower().strip()
                _sid = None
                for _sp in _graph["graph"]["species"]:
                    _spid = _sp.get("id", "").lower().replace("_", " ")
                    if _q in _spid or _q in _sp.get("name", "").lower() or _q in _sp.get("chinese", "").lower():
                        _sid = _sp.get("id", "")
                        _sci_name = _sp.get("name", "")
                        _cn_name = _sp.get("chinese", "")
                        break

                if _sid:
                    _papers = [p for p in _graph["graph"]["papers"] if _sid in p.get("species", [])]
                    # Credibility scoring
                    try:
                        _SCRIPTS = str(_WORKSPACE_ROOT / "fish-ecology-assistant" / "scripts")
                        if _SCRIPTS not in sys.path:
                            sys.path.insert(0, _SCRIPTS)
                        from credibility_scorer import score_papers, format_credibility
                        _scored = score_papers(_papers, species_name=_sci_name)
                    except ImportError:
                        _scored = _papers

                    _high = sum(1 for p in _scored if p.get("_credibility_label") == "高")
                    _mid = sum(1 for p in _scored if p.get("_credibility_label") == "中")
                    _low = sum(1 for p in _scored if p.get("_credibility_label") == "低")

                    result["phases"]["2_cognitive"] = {
                        "total_papers": len(_papers),
                        "scientific_name": _sci_name,
                        "chinese_name": _cn_name,
                        "mode": "graph_lookup",
                        "credibility": {"high": _high, "mid": _mid, "low": _low},
                        "papers_top5": [
                            {"title": p.get("title","")[:60], "score": p.get("_credibility_score",0),
                             "flag": p.get("_credibility_flag","")}
                            for p in sorted(_scored, key=lambda x: x.get("_credibility_score",0), reverse=True)[:5]
                        ],
                    }
                    elapsed = time.time() - phase2_start
                    total = result["phases"]["2_cognitive"]["total_papers"]
                    _log(f"     ✓ 图谱: {_cn_name} — {total} 篇论文")
                    _log(f"     ✓ 置信: 🟢{_high} 🟡{_mid} 🟠{_low}")
                    _log(f"     ⏱ {elapsed:.2f}s")

                    # Writeback: update KB papers_in_graph count
                    try:
                        _kb_path = _WORKSPACE_ROOT / "fish-ecology-assistant" / "config" / "fish_species_kb.yaml"
                        if _kb_path.exists():
                            with open(_kb_path, encoding="utf-8") as f:
                                _kb_content = f.read()

                            # Find tribolodon_brandti entry and update papers_in_graph
                            _old = f"papers_in_graph:"
                            _new = f"papers_in_graph: {total}"
                            if _old in _kb_content:
                                # Replace just the first occurrence after the species section
                                _pos = _kb_content.find(_old)
                                if _pos > 0:
                                    _line_end = _kb_content.find("\n", _pos)
                                    _old_line = _kb_content[_pos:_line_end]
                                    _kb_content = _kb_content.replace(_old_line, f"papers_in_graph: {total}", 1)
                                    with open(_kb_path, "w", encoding="utf-8") as f:
                                        f.write(_kb_content)
                                    _log(f"     ✓ KB写回: papers_in_graph {total}")
                    except Exception:
                        pass  # best effort

                    # Self-evolve log
                    try:
                        _se_path = _WORKSPACE_ROOT / "fish-ecology-assistant" / "scripts" / "self_evolve.py"
                        if _se_path.exists():
                            sys.path.insert(0, str(_se_path.parent))
                            from self_evolve import log_search
                            log_search(_sid, {
                                "mode": "graph_lookup",
                                "layers_activated": ["L1"],
                                "layers_producing": {"L1": total},
                                "known_papers": total,
                                "new_papers": 0,
                                "total_papers": total,
                                "tokens_estimated": 0,
                                "mode_auto": False,
                            })
                    except Exception:
                        pass
                else:
                    # Fallback: coordinated_search
                    lit_result = search_species(species_name)
                    result["phases"]["2_cognitive"] = {
                        "total_papers": getattr(lit_result, "total_papers", 0),
                        "mode": "coordinated_search",
                    }
                    _log(f"     ✓ 引擎检索: {getattr(lit_result, 'total_papers', 0)} 篇")
            else:
                _log("     ⚠️ species_graph.yaml not found, skipping")
                result["phases"]["2_cognitive"] = {"status": "skipped"}
        except Exception as exc:
            _log(f"     ❌ 文献搜索失败: {exc}")
            result["phases"]["2_cognitive"] = {"status": "error", "error": str(exc)}
    else:
        result["phases"]["2_cognitive"] = {"status": "skipped"}

    # ══════════════════════════════════════════
    # Phase 3: 冲突裁决 (已由 lookup_species 内置计算)
    # ══════════════════════════════════════════
    if enable_conflict:
        _log(f"\n▸ Phase 3/5: {_phase_tag('3', 'conflict-arbiter')} 协同")
        _log("  └─ 冲突裁决 (已内嵌在 lookup_species 中)...")
        phase3_start = time.time()
        try:
            # 从 Phase 1 结果读取已计算的 conflict_verdict
            p1 = result.get("phases", {}).get("1_fish", {})
            cv = p1.get("conflict_verdict", {})
            if cv:
                result["phases"]["3_conflict"] = {
                    "conflict_level": cv.get("conflict_level"),
                    "consensus": cv.get("consensus"),
                    "verdict": cv.get("verdict"),
                    "source": "lookup_species (内嵌)",
                }
                elapsed = time.time() - phase3_start
                cl = cv.get("conflict_level", "?")
                _log(f"     ✓ 裁决: {cv.get('verdict','')}")
                _log(f"     ⏱ {elapsed:.2f}s")
            else:
                _log(f"     ⚠️ 无保护等级数据，跳过仲裁")
                result["phases"]["3_conflict"] = {"status": "skipped", "reason": "无保护等级数据"}
        except Exception as exc:
            _log(f"     ❌ 读取冲突裁决失败: {exc}")
            result["phases"]["3_conflict"] = {"status": "error", "error": str(exc)}
    else:
        result["phases"]["3_conflict"] = {"status": "skipped"}
# ══════════════════════════════════════════
    # Phase 4: 协同 — 领域评估
    # ══════════════════════════════════════════
    if enable_assessment:
        _log(f"\n▸ Phase 4/5: {_phase_tag('4', 'porpoise/coilia')} 协同")
        _log("  └─ 领域评估...")
        phase4_start = time.time()
        try:
            # 自动路由: 江豚→porpoise, 刀鲚→coilia, 其他→porpoise 兜底
            if "江豚" in species_name or "neophocaena" in species_name.lower():
                domain_result = assess_conservation(species_name)
                domain_source = "porpoise-agent"
            elif "刀鲚" in species_name or "coilia" in species_name.lower():
                domain_result = assess_species(species_name)
                domain_source = "coilia-agent"
            else:
                domain_result = assess_conservation(species_name)
                domain_source = "porpoise-agent"

            result["phases"]["4_domain"] = {
                "source": domain_source,
                "data": domain_result,
            }
            elapsed = time.time() - phase4_start
            status = domain_result.get("status", "ok")
            _log(f"     ✓ {domain_source} 评估完成 (status={status})")
            _log(f"     ⏱ {elapsed:.2f}s")
        except Exception as exc:
            _log(f"     ❌ 领域评估失败: {exc}")
            result["phases"]["4_domain"] = {"status": "error", "error": str(exc)}
    else:
        result["phases"]["4_domain"] = {"status": "skipped"}

    # ══════════════════════════════════════════
    # Phase 5: 全栈汇总
    # ══════════════════════════════════════════
    _log(f"\n▸ Phase 5/5: 全栈汇总")

    summary_lines = [f"  {_status_icon(result['phases'].get('1_fish', {}))} 知识库",]

    if enable_literature:
        p2 = result["phases"].get("2_cognitive", {})
        total = p2.get("total_papers", 0)
        summary_lines.append(f"  {_status_icon(p2)} 文献 ({total} 篇)")

    if enable_conflict:
        p3 = result["phases"].get("3_conflict", {})
        cl = p3.get("conflict_level", "?")
        summary_lines.append(f"  {_status_icon(p3)} 冲突仲裁 (等级 {cl})")

    if enable_assessment:
        p4 = result["phases"].get("4_domain", {})
        summary_lines.append(f"  {_status_icon(p4)} 领域评估")

    _log("\n管线汇总:")
    for line in summary_lines:
        _log(line)

    _log("\n" + "=" * 60)
    _log(f"  管线完成 ✓")
    _log("=" * 60)

    result["completed_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    result["summary"] = "\n".join(summary_lines)
    return result


# ══════════════════════════════════════════════
# 辅助函数
# ══════════════════════════════════════════════

def _get_workspace_root() -> Path:
    """返回 Reasonix 工作区根目录 D:/Reasonix。

    无论从哪个目录运行此脚本都能正确解析。
    """
    from pathlib import Path
    # pipeline_fish_collab.py 位于 workspace/ 目录下
    # workspace 目录的父目录就是 D:\Reasonix
    return Path(__file__).resolve().parent.parent


def _phase_tag(phase: str, project: str) -> str:
    """生成阶段标签。"""
    tags = {
        "fish-ecology-assistant": "f项目",
        "cognitive-search-engine": "木",
        "conflict-arbiter": "火",
        "porpoise/coilia": "金/水",
    }
    tag = tags.get(project, project)
    return f"[Phase {phase} · {tag}]"


def _status_icon(phase_data: Dict[str, Any]) -> str:
    """根据阶段数据返回状态图标。"""
    status = phase_data.get("status", "")
    if status == "error":
        return "❌"
    if status == "skipped":
        return "⏭️"
    if phase_data.get("known_species") or phase_data.get("total_papers", 0) > 0:
        return "✅"
    return "⚠️"


def _build_conflict_sources(profile: Dict[str, Any]) -> List[Dict[str, Any]]:
    """从知识库查询结果中提取保护等级，构建冲突检测输入。"""
    sources = []
    sd = profile.get("species_data", {}) or profile

    # IUCN 状态
    iucn = sd.get("iucn", sd.get("iucn_status", ""))
    if iucn:
        sources.append({"source": "iucn", "iucn": iucn})

    # 中国保护等级
    prot = sd.get("protection_level", "")
    if prot:
        sources.append({"source": "chinese_red_list", "protection_level": prot})

    # 本地 conservation 字段
    cons = sd.get("conservation", "")
    if cons and cons not in ("", "无"):
        sources.append({"source": "provincial_protection", "protection_level": cons})

    return sources


def _fallback_conflict_sources(profile: Dict[str, Any]) -> List[Dict[str, Any]]:
    """当 KB 无显式保护等级时，从 category 推导冲突检测输入。"""
    sources = []
    sd = profile.get("species_data", {}) or profile
    category = sd.get("category", "")

    # category → 推定保护等级映射
    CATEGORY_MAP = {
        "protected_recorded":     {"iucn": "CR", "protection_level": "国家一级"},
        "protected_missing":      {"iucn": "EX", "protection_level": "国家一级"},
        "endangered_in_graph":    {"protection_level": "濒危"},
        "mammal":                 {"iucn": "CR", "protection_level": "国家一级"},
        "dominant":               {"iucn": "LC", "conservation": "无"},
        "diadromous":             {"conservation": "无"},
    }
    mapping = CATEGORY_MAP.get(category, {})
    if mapping:
        if "iucn" in mapping:
            sources.append({"source": "iucn", "iucn": mapping["iucn"]})
        if "protection_level" in mapping:
            sources.append({"source": "chinese_red_list", "protection_level": mapping["protection_level"]})
        if "conservation" in mapping:
            sources.append({"source": "provincial_protection", "protection_level": mapping["conservation"]})
    return sources


def _clean_verdict(verdict: str) -> str:
    """去除裁决文本中的 emoji 以便终端显示。"""
    import re
    return re.sub(r'[\U0001F300-\U0001FFFF\U00002000-\U00002BFF]', '', verdict)


# ══════════════════════════════════════════════
# 命令行入口
# ══════════════════════════════════════════════

if __name__ == "__main__":
    species = sys.argv[1] if len(sys.argv) > 1 else "鳤"
    result = run_fish_pipeline(
        species,
        enable_literature=True,
        enable_conflict=True,
        enable_assessment="--full" in sys.argv,
    )
    sys.exit(0 if result.get("status") == "ok" else 1)
