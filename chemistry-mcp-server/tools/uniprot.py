"""UniProt tools — protein sequences, functions, metadata."""

import asyncio
import json

import httpx

from server import mcp

UNIPROT_REST = "https://rest.uniprot.org"


def _uniprot_get(endpoint: str, params: dict | None = None, accept: str = "json") -> dict | None:
    r = httpx.get(f"{UNIPROT_REST}/{endpoint}", params=params, timeout=30.0,
                  headers={"Accept": f"application/{accept}"})
    if r.status_code != 200:
        return None
    return r.json()


def _clean_entry(entry: dict) -> dict:
    """Extract key fields from a UniProt entry."""
    acc = entry.get("primaryAccession", "")
    organism = entry.get("organism", {})
    comments = entry.get("comments", [])
    keywords = entry.get("keywords", [])

    # Extract function description
    function = ""
    for c in comments:
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts", [])
            function = texts[0].get("value", "") if texts else ""
            break

    # Extract catalytic activity
    catalytic = ""
    for c in comments:
        if c.get("commentType") == "CATALYTIC_ACTIVITY":
            reaction = c.get("reaction", {})
            catalytic = reaction.get("name", "")
            break

    # Extract subcellular location
    locations = []
    for c in comments:
        if c.get("commentType") == "SUBCELLULAR_LOCATION":
            for loc in c.get("subcellularLocations", []):
                loc_name = loc.get("location", {}).get("value", "")
                if loc_name:
                    locations.append(loc_name)

    # Extract gene names
    genes = entry.get("genes", [])
    gene_name = ""
    gene_synonyms = []
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")
        gene_synonyms = [s.get("value", "") for s in genes[0].get("synonyms", [])]

    # Extract database cross-references
    db_refs = {}
    for ref in entry.get("uniProtKBCrossReferences", []):
        db = ref.get("database", "")
        if db in ("PDB", "KEGG", "ChEMBL", "DrugBank", "GeneID", "Ensembl"):
            db_refs.setdefault(db, []).append(ref.get("id", ""))

    # Sequence
    sequence = entry.get("sequence", {})

    return {
        "accession": acc,
        "id": entry.get("uniProtKBEntryId", acc),
        "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
        "gene_name": gene_name,
        "gene_synonyms": gene_synonyms[:5],
        "organism": organism.get("scientificName", ""),
        "organism_id": str(organism.get("taxonId", "")),
        "sequence_length": sequence.get("length", 0),
        "sequence": sequence.get("value", ""),
        "molecular_weight": sequence.get("molWeight", 0),
        "function": function,
        "catalytic_activity": catalytic,
        "subcellular_locations": locations[:5],
        "keywords": [k.get("name", "") for k in keywords][:10],
        "reviewed": entry.get("entryType", "") == "UniProtKB-reviewed (Swiss-Prot)",
        "database_cross_references": db_refs,
    }


# ─── Search ───────────────────────────────────────────────────

@mcp.tool()
async def uniprot_search(
    query: str,
    organism: str | None = None,
    reviewed: bool = True,
    max_results: int = 10,
) -> str:
    """Search UniProt for proteins by name, gene, or keyword.

    Args:
        query: Search term (protein name, gene name, keyword)
        organism: Optional organism filter (e.g. 'human', 'mouse')
        reviewed: Only reviewed (Swiss-Prot) entries (default True)
        max_results: Max results (1-25, default 10)

    Returns:
        JSON array of protein entries with accession, name, gene, organism, sequence
    """
    def _search():
        lucene_parts = [query]
        if reviewed:
            lucene_parts.append("reviewed:true")
        if organism:
            lucene_parts.append(f"organism_name:{organism}")

        full_query = " AND ".join(lucene_parts)
        r = httpx.get(
            f"{UNIPROT_REST}/uniprotkb/search",
            params={
                "query": full_query,
                "size": max(1, min(max_results, 25)),
                "fields": "accession,id,protein_name,gene_names,organism_name,length",
            },
            headers={"Accept": "application/json"},
            timeout=30.0,
        )
        if r.status_code != 200:
            return []
        data = r.json()

        results = []
        for entry in data.get("results", []):
            # Lightweight version for search results
            genes = entry.get("genes", [])
            gene_name = genes[0].get("geneName", {}).get("value", "") if genes else ""
            organism_name = entry.get("organism", {}).get("scientificName", "")

            desc = entry.get("proteinDescription", {})
            prot_name = desc.get("recommendedName", {}).get("fullName", {}).get("value", "")

            seq = entry.get("sequence", {})
            results.append({
                "accession": entry.get("primaryAccession", ""),
                "id": entry.get("uniProtKBEntryId", ""),
                "protein_name": prot_name,
                "gene_name": gene_name,
                "organism": organism_name,
                "sequence_length": seq.get("length", 0),
            })
        return results

    results = await asyncio.to_thread(_search)
    return json.dumps(results, indent=2, ensure_ascii=False)


# ─── Get protein details ──────────────────────────────────────

@mcp.tool()
async def uniprot_get_protein(accession: str) -> str:
    """Get detailed protein information from UniProt by accession ID.

    Args:
        accession: UniProt accession (e.g. 'P00533' for EGFR)

    Returns:
        JSON with protein name, gene, organism, sequence, function, locations, cross-references
    """
    def _get():
        data = _uniprot_get(f"uniprotkb/{accession}")
        if data is None:
            raise ValueError(f"Protein not found: {accession}")
        return _clean_entry(data)

    result = await asyncio.to_thread(_get)
    return json.dumps(result, indent=2, ensure_ascii=False)


# ─── ID mapping ───────────────────────────────────────────────

@mcp.tool()
async def uniprot_map_ids(
    ids: list[str],
    from_type: str = "UniProtKB_AC-ID",
    to_type: str = "GeneID",
) -> str:
    """Map protein identifiers between databases (UniProt, GeneID, Ensembl, PDB, KEGG, etc.).

    Args:
        ids: List of identifiers to map (max 100)
        from_type: Source database type (e.g. UniProtKB_AC-ID, GeneID, Ensembl, PDB)
        to_type: Target database type (e.g. GeneID, UniProtKB_AC-ID, Ensembl, ChEMBL)

    Returns:
        JSON array of {from_id, to_id} mapping pairs
    """
    def _map():
        ids_limited = ids[:100]
        # Submit mapping job
        r = httpx.post(
            f"{UNIPROT_REST}/idmapping/run",
            data={"from": from_type, "to": to_type, "ids": " ".join(ids_limited)},
            timeout=30.0,
        )
        if r.status_code != 200:
            raise ValueError(f"ID mapping failed: {r.text[:200]}")
        job_id = r.json()["jobId"]

        # Poll for completion
        for _ in range(30):
            status_r = httpx.get(f"{UNIPROT_REST}/idmapping/status/{job_id}", timeout=10.0)
            status = status_r.json()
            if status.get("jobStatus") == "FINISHED":
                break
            if status.get("jobStatus") == "FAILURE":
                raise ValueError(f"Mapping job failed: {status.get('failureReason', '')}")
            import time
            time.sleep(1)
        else:
            raise ValueError("ID mapping timed out")

        # Fetch results
        results_r = httpx.get(f"{UNIPROT_REST}/idmapping/stream/{job_id}", timeout=30.0)
        if results_r.status_code != 200:
            return []
        data = results_r.json()
        return [{"from": r.get("from", ""), "to": r.get("to", "")} for r in data.get("results", [])]

    results = await asyncio.to_thread(_map)
    return json.dumps(results, indent=2, ensure_ascii=False)
