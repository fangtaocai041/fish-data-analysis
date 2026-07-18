"""
鱼类/水生生物物种综合查询工具 — 聚合 GBIF + Wikipedia

功能：
  1. get_species_profile(scientific_name) — 获取物种完整档案
  2. search_species(common_name) — 按俗名/学名搜索

用法示例：
  python scripts/species_tools.py profile "Coilia nasus"
  python scripts/species_tools.py search "江豚"
  python scripts/species_tools.py profile "Neophocaena asiaeorientalis" --lang en

输出：结构化 JSON，包含分类、保护状态、分布、描述等
"""

from __future__ import annotations

import json
import os
import sys
import urllib.parse
import urllib.request
from typing import Any

USER_AGENT = "Reasonix-Species-Tool/1.0 (fish-ecology-agent)"
_DEBUG = os.environ.get("SPECIES_TOOL_DEBUG", "0") == "1"


def _log_warning(msg: str) -> None:
    """输出警告到 stderr。"""
    if _DEBUG:
        import traceback
        traceback.print_exc()
    print(f"[warn] {msg}", file=sys.stderr)


# ── 网络/代理 ──────────────────────────────────────────


def _make_opener() -> urllib.request.OpenerDirector:
    """智能选择 opener：代理可达时用系统代理，否则直连。"""
    import socket

    for var in ["HTTPS_PROXY", "HTTP_PROXY", "https_proxy", "http_proxy"]:
        proxy_url = os.environ.get(var, "")
        if proxy_url:
            try:
                raw = proxy_url.replace("http://", "").replace("https://", "")
                host, port_str = raw.split(":")
                s = socket.create_connection((host, int(port_str)), timeout=1)
                s.close()
                return urllib.request.build_opener()  # 使用系统默认代理
            except (OSError, ValueError, IndexError):
                continue
    return urllib.request.build_opener(urllib.request.ProxyHandler({}))


def _fetch_json(url: str) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with _make_opener().open(req, timeout=15) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _fetch_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with _make_opener().open(req, timeout=15) as resp:
        return resp.read().decode("utf-8")


# ── Wikipedia 查询 ────────────────────────────────────


def _wiki_summary(title: str, lang: str = "zh") -> dict[str, Any]:
    """查询 Wikipedia 摘要。"""
    url = f"https://{lang}.wikipedia.org/api/rest_v1/page/summary/{urllib.parse.quote(title.replace(' ', '_'))}"
    try:
        data = _fetch_json(url)
        return {
            "title": data.get("title", title),
            "extract": data.get("extract", ""),
            "url": data.get("content_urls", {}).get("desktop", {}).get("page", ""),
            "image": data.get("originalimage", {}).get("source", ""),
        }
    except Exception:
        _log_warning(f"Wikipedia summary '{title}' failed (lang={lang})")
        return {}


# ── GBIF 查询 ─────────────────────────────────────────


def _gbif_search(name: str) -> list[dict[str, Any]]:
    """通过 GBIF 搜索物种。"""
    url = f"https://api.gbif.org/v1/species/match?verbose=true&name={urllib.parse.quote(name)}"
    try:
        data = _fetch_json(url)
        if data.get("matchType", "") == "NONE":
            return []
        return [{
            "scientific_name": data.get("scientificName", name),
            "canonical_name": data.get("canonicalName", name),
            "rank": data.get("rank", ""),
            "status": data.get("status", ""),
            "taxon_key": data.get("speciesKey", data.get("usageKey")),
            "kingdom": data.get("kingdom", ""),
            "phylum": data.get("phylum", ""),
            "class": data.get("class", ""),
            "order": data.get("order", ""),
            "family": data.get("family", ""),
            "genus": data.get("genus", ""),
        }]
    except Exception:
        _log_warning(f"GBIF species match '{name}' failed")
        return []


def _gbif_iucn(taxon_key: int) -> str:
    """查询 GBIF 的 IUCN 状态。"""
    url = f"https://api.gbif.org/v1/species/{taxon_key}/iucn"
    try:
        data = _fetch_json(url)
        return data.get("category", "未知")
    except Exception:
        _log_warning(f"GBIF IUCN query for taxon {taxon_key} failed")
        return "未知"


def _gbif_occurrences(taxon_key: int) -> dict[str, Any]:
    """获取 GBIF 分布统计。"""
    url = f"https://api.gbif.org/v1/occurrence/counts?taxonKey={taxon_key}"
    try:
        data = _fetch_json(url)
        return {"total": data.get("total", 0)}
    except Exception:
        _log_warning(f"GBIF occurrence count for taxon {taxon_key} failed")
        return {"total": 0}


# ── 公开函数 ──────────────────────────────────────────


def search_species(query: str, limit: int = 5) -> list[dict[str, Any]]:
    """
    按名称搜索物种（同时查 GBIF + 中文维基）。

    参数:
        query: 物种名称（学名或俗名）
        limit: 结果上限

    返回:
        [{"name": ..., "scientific_name": ..., "rank": ..., ...}, ...]
    """
    results: list[dict[str, Any]] = []

    # GBIF 搜索
    gbif_url = f"https://api.gbif.org/v1/species/search?q={urllib.parse.quote(query)}&limit={limit}"
    try:
        data = _fetch_json(gbif_url)
        for r in data.get("results", []):
            results.append({
                "name": r.get("scientificName", ""),
                "canonical_name": r.get("canonicalName", ""),
                "rank": r.get("rank", ""),
                "taxon_key": r.get("key"),
                "kingdom": r.get("kingdom", ""),
                "phylum": r.get("phylum", ""),
                "class": r.get("class", ""),
                "order": r.get("order", ""),
                "family": r.get("family", ""),
                "genus": r.get("genus", ""),
                "source": "GBIF",
            })
    except Exception:
        _log_warning(f"GBIF species search '{query}' failed")

    # Wikipedia 搜索（中文）
    wiki_url = f"https://zh.wikipedia.org/w/api.php?action=query&list=search&srsearch={urllib.parse.quote(query)}&srlimit={limit}&format=json&origin=*"
    try:
        data = _fetch_json(wiki_url)
        for r in data.get("query", {}).get("search", []):
            results.append({
                "name": r.get("title", ""),
                "page_id": r.get("pageid"),
                "snippet": r.get("snippet", "").replace("<span class=\"searchmatch\">", "").replace("</span>", ""),
                "source": "Wikipedia(zh)",
            })
    except Exception:
        _log_warning(f"Wikipedia search '{query}' failed")

    return results


def get_species_profile(scientific_name: str, lang: str = "zh") -> dict[str, Any]:
    """
    获取物种完整档案（聚合 GBIF + Wikipedia）。

    参数:
        scientific_name: 学名（如 "Coilia nasus"）
        lang: Wikipedia 语言

    返回:
        {"taxonomy": ..., "conservation": ..., "description": ..., "distribution": ...}
    """
    profile: dict[str, Any] = {
        "scientific_name": scientific_name,
        "sources": [],
    }

    # 1. GBIF 分类 + IUCN
    gbif_matches = _gbif_search(scientific_name)
    if gbif_matches:
        tax = gbif_matches[0]
        profile["taxonomy"] = {
            "kingdom": tax["kingdom"],
            "phylum": tax["phylum"],
            "class": tax["class"],
            "order": tax["order"],
            "family": tax["family"],
            "genus": tax["genus"],
            "scientific_name": tax["scientific_name"],
            "rank": tax["rank"],
            "status": tax["status"],
        }
        profile["sources"].append("GBIF")
        if tax["taxon_key"]:
            profile["conservation"] = {"iucn": _gbif_iucn(tax["taxon_key"])}
            occ = _gbif_occurrences(tax["taxon_key"])
            profile["occurrence_count"] = occ.get("total", 0)

    # 2. Wikipedia 描述
    wiki_info = _wiki_summary(scientific_name, lang="en")
    if wiki_info:
        profile["description_en"] = wiki_info.get("extract", "")
        profile["image"] = wiki_info.get("image", "")
        profile["wiki_url"] = wiki_info.get("url", "")
        profile["sources"].append("Wikipedia(en)")

    # 中文维基
    if lang == "zh":
        wiki_zh = _wiki_summary(urllib.parse.unquote(scientific_name), lang="zh")
        if wiki_zh and wiki_zh.get("extract"):
            profile["description_zh"] = wiki_zh.get("extract", "")
            if not profile.get("image"):
                profile["image"] = wiki_zh.get("image", "")
            if not profile.get("wiki_url"):
                profile["wiki_url"] = wiki_zh.get("url", "")
            profile["sources"].append("Wikipedia(zh)")

    # 中文常见名搜索
    if lang == "zh":
        wiki_search_url = f"https://zh.wikipedia.org/w/api.php?action=query&list=search&srsearch={urllib.parse.quote(scientific_name)}&srlimit=3&format=json&origin=*"
        try:
            data = _fetch_json(wiki_search_url)
            for r in data.get("query", {}).get("search", []):
                title = r.get("title", "")
                if title and title != scientific_name:
                    wiki_alt = _wiki_summary(title, lang="zh")
                    if wiki_alt and wiki_alt.get("extract"):
                        profile["common_name"] = title
                        if not profile.get("description_zh"):
                            profile["description_zh"] = wiki_alt.get("extract", "")
                        break
        except Exception:
            _log_warning(f"Wikipedia common name search '{scientific_name}' failed")

    return profile


# ── CLI ───────────────────────────────────────────────


def main():
    if len(sys.argv) < 2:
        print("用法: python species_tools.py <命令> [参数...]")
        print("")
        print("命令:")
        print("  search  <关键词>  [--limit 5]    — 搜索物种")
        print("  profile <学名>    [--lang zh]    — 获取完整档案")
        sys.exit(1)

    command = sys.argv[1]
    lang = "zh"
    limit = 5
    args = sys.argv[2:]
    positional: list[str] = []

    i = 0
    while i < len(args):
        if args[i] == "--lang" and i + 1 < len(args):
            lang = args[i + 1]
            i += 2
        elif args[i] == "--limit" and i + 1 < len(args):
            limit = int(args[i + 1])
            i += 2
        else:
            positional.append(args[i])
            i += 1

    query = " ".join(positional) if positional else ""

    if command == "search":
        if not query:
            print("错误: 需要搜索关键词"); sys.exit(1)
        results = search_species(query, limit=limit)
        print(json.dumps(results, ensure_ascii=False, indent=2))

    elif command == "profile":
        if not query:
            print("错误: 需要学名"); sys.exit(1)
        profile = get_species_profile(query, lang=lang)
        print(json.dumps(profile, ensure_ascii=False, indent=2))

    else:
        print(f"未知命令: {command}")
        sys.exit(1)


if __name__ == "__main__":
    main()
