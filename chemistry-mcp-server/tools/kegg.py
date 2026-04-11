"""KEGG tools — pathways, compounds, genes, diseases."""

import asyncio
import json

import httpx

from server import mcp

KEGG_API = "https://rest.kegg.jp"


def _kegg_get(endpoint: str) -> str | None:
    r = httpx.get(f"{KEGG_API}/{endpoint}", timeout=30.0)
    if r.status_code != 200:
        return None
    return r.text


def _parse_kegg_list(text: str) -> list[dict]:
    """Parse KEGG list format (tab-separated key-value pairs)."""
    results = []
    for line in text.strip().splitlines():
        if "\t" not in line:
            continue
        parts = line.split("\t", 1)
        results.append({"id": parts[0], "name": parts[1] if len(parts) > 1 else ""})
    return results


def _parse_kegg_entry(text: str) -> dict:
    """Parse KEGG flat-file entry into a dict."""
    result = {}
    current_key = None
    for line in text.splitlines():
        if not line:
            continue
        if line[0] != " ":
            key, _, value = line.partition(" ")
            current_key = key
            if current_key in result:
                if isinstance(result[current_key], list):
                    result[current_key].append(value.strip())
                else:
                    result[current_key] = [result[current_key], value.strip()]
            else:
                result[current_key] = value.strip()
        elif current_key:
            if isinstance(result[current_key], list):
                result[current_key][-1] += " " + line.strip()
            else:
                result[current_key] += " " + line.strip()
    return result


# ─── Search ───────────────────────────────────────────────────

@mcp.tool()
async def kegg_search(
    database: str,
    query: str,
    max_results: int = 20,
) -> str:
    """Search KEGG database for entries.

    Args:
        database: KEGG database to search: pathway, compound, drug, disease, enzyme, genes, organism
        query: Search query string
        max_results: Max results (1-50, default 20)

    Returns:
        JSON array of {id, name} entries
    """
    def _search():
        text = _kegg_get(f"find/{database}/{query}")
        if text is None:
            return []
        results = _parse_kegg_list(text)
        return results[:max(1, min(max_results, 50))]

    results = await asyncio.to_thread(_search)
    return json.dumps(results, indent=2, ensure_ascii=False)


# ─── Get entry details ────────────────────────────────────────

@mcp.tool()
async def kegg_get_entry(entry_id: str) -> str:
    """Get detailed information about a KEGG entry.

    Args:
        entry_id: KEGG entry ID (e.g. 'cpd:C00031' for glucose, 'hsa:10458' for a gene,
                  'path:hsa00010' for glycolysis, 'D00075' for a drug)

    Returns:
        JSON dictionary with entry details (name, definition, pathway, genes, etc.)
    """
    def _get():
        text = _kegg_get(f"get/{entry_id}")
        if text is None:
            raise ValueError(f"KEGG entry not found: {entry_id}")
        parsed = _parse_kegg_entry(text)
        return parsed

    result = await asyncio.to_thread(_get)
    return json.dumps(result, indent=2, ensure_ascii=False)


# ─── List entries ─────────────────────────────────────────────

@mcp.tool()
async def kegg_list(
    database: str,
    organism: str | None = None,
) -> str:
    """List entries in a KEGG database.

    Args:
        database: KEGG database: pathway, compound, drug, disease, enzyme, module, reaction
        organism: Optional organism code (e.g. 'hsa' for human, 'mmu' for mouse)

    Returns:
        JSON array of {id, name} entries
    """
    def _list():
        endpoint = f"list/{database}"
        if organism:
            endpoint += f"/{organism}"
        text = _kegg_get(endpoint)
        if text is None:
            return []
        return _parse_kegg_list(text)

    results = await asyncio.to_thread(_list)
    return json.dumps(results, indent=2, ensure_ascii=False)


# ─── Link (cross-references) ──────────────────────────────────

@mcp.tool()
async def kegg_link(
    target_db: str,
    source_id: str,
) -> str:
    """Find cross-references between KEGG entries (e.g. genes in a pathway, compounds in a reaction).

    Args:
        target_db: Target database (e.g. pathway, compound, gene, enzyme)
        source_id: Source entry ID or database (e.g. 'hsa:10458', 'path:hsa00010', 'hsa')

    Returns:
        JSON array of {source_id, target_id} pairs
    """
    def _link():
        text = _kegg_get(f"link/{target_db}/{source_id}")
        if text is None:
            return []
        pairs = []
        for line in text.strip().splitlines():
            if "\t" not in line:
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                pairs.append({"source_id": parts[0], "target_id": parts[1]})
        return pairs

    results = await asyncio.to_thread(_link)
    return json.dumps(results, indent=2, ensure_ascii=False)


# ─── Pathway compounds ────────────────────────────────────────

@mcp.tool()
async def kegg_pathway_compounds(pathway_id: str) -> str:
    """Get all compounds involved in a KEGG pathway.

    Args:
        pathway_id: KEGG pathway ID (e.g. 'hsa00010' for glycolysis)

    Returns:
        JSON array of compound IDs
    """
    def _get():
        text = _kegg_get(f"link/cpd/{pathway_id}")
        if text is None:
            return []
        compounds = []
        for line in text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 2 and "cpd:" in parts[1]:
                compounds.append(parts[1])
        return compounds

    results = await asyncio.to_thread(_get)
    return json.dumps(results, indent=2, ensure_ascii=False)
