"""ChEMBL tools — bioactivity data, drug targets, molecule search."""

import asyncio
import json

import httpx

from server import mcp

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"


def _get(endpoint: str, params: dict | None = None) -> dict | list | None:
    params = params or {}
    params["format"] = "json"
    r = httpx.get(f"{CHEMBL_API}/{endpoint}", params=params, timeout=30.0)
    if r.status_code != 200:
        return None
    return r.json()


# ─── Molecule search ──────────────────────────────────────────

@mcp.tool()
async def chembl_search_molecule(
    query: str,
    search_type: str = "name",
    max_results: int = 10,
) -> str:
    """Search ChEMBL for molecules by name, SMILES, substructure, or similarity.

    Args:
        query: Molecule name, SMILES, or ChEMBL ID
        search_type: One of: name, substructure, similarity, chembl_id
        max_results: Max results to return (1-20, default 10)

    Returns:
        JSON array of molecules with ChEMBL ID, name, SMILES, properties
    """
    def _search():
        if search_type == "chembl_id":
            data = _get(f"molecule/{query}")
            if data is None:
                return []
            return [_clean_molecule(data)]

        if search_type == "name":
            data = _get("molecule/search", {"q": query, "limit": max(1, min(max_results, 20))})
        elif search_type == "substructure":
            data = _get("substructure", {"smiles": query, "limit": max(1, min(max_results, 20))})
        elif search_type == "similarity":
            data = _get("similarity", {"smiles": query, "limit": max(1, min(max_results, 20))})
        else:
            raise ValueError(f"Invalid search_type: {search_type!r}. Use: name, substructure, similarity, chembl_id")

        if data is None:
            return []
        molecules = data.get("molecules", [])
        return [_clean_molecule(m) for m in molecules]

    results = await asyncio.to_thread(_search)
    return json.dumps(results, indent=2, ensure_ascii=False)


def _clean_molecule(m: dict) -> dict:
    """Extract key fields from a ChEMBL molecule record."""
    props = m.get("molecule_properties", {}) or {}
    struct = m.get("molecule_structures", {}) or {}
    return {
        "chembl_id": m.get("molecule_chembl_id"),
        "name": m.get("pref_name"),
        "max_phase": m.get("max_phase"),
        "molecule_type": m.get("molecule_type"),
        "smiles": struct.get("canonical_smiles"),
        "inchi": struct.get("standard_inchi"),
        "inchikey": struct.get("standard_inchi_key"),
        "mw": props.get("full_mwt"),
        "logp": props.get("alogp"),
        "psa": props.get("psa"),
        "hbd": props.get("hbd"),
        "hba": props.get("hba"),
        "ro3_pass": props.get("num_ro5_violations"),
    }


# ─── Bioactivity ──────────────────────────────────────────────

@mcp.tool()
async def chembl_get_bioactivity(
    chembl_id: str,
    target_chembl_id: str | None = None,
    activity_type: str | None = None,
    max_results: int = 20,
) -> str:
    """Get bioactivity data (IC50, Ki, EC50 etc.) for a molecule from ChEMBL.

    Args:
        chembl_id: ChEMBL molecule ID (e.g. 'CHEMBL25')
        target_chembl_id: Optional target ChEMBL ID to filter results
        activity_type: Optional filter: IC50, Ki, EC50, Kd, etc.
        max_results: Max results (1-50, default 20)

    Returns:
        JSON array of activity records with target, value, units, assay info
    """
    def _get_activity():
        params = {
            "molecule_chembl_id": chembl_id,
            "limit": max(1, min(max_results, 50)),
        }
        if target_chembl_id:
            params["target_chembl_id"] = target_chembl_id
        if activity_type:
            params["standard_type"] = activity_type

        data = _get("activity", params)
        if data is None:
            return []
        return [_clean_activity(a) for a in data.get("activities", [])]

    results = await asyncio.to_thread(_get_activity)
    return json.dumps(results, indent=2, ensure_ascii=False)


def _clean_activity(a: dict) -> dict:
    return {
        "activity_id": a.get("activity_id"),
        "type": a.get("standard_type"),
        "value": a.get("standard_value"),
        "units": a.get("standard_units"),
        "relation": a.get("standard_relation"),
        "target_chembl_id": a.get("target_chembl_id"),
        "target_name": a.get("target_pref_name"),
        "target_type": a.get("target_type"),
        "assay_chembl_id": a.get("assay_chembl_id"),
        "assay_description": a.get("assay_description"),
        "pchembl_value": a.get("pchembl_value"),
    }


# ─── Targets ──────────────────────────────────────────────────

@mcp.tool()
async def chembl_search_target(
    query: str,
    target_type: str | None = None,
    max_results: int = 10,
) -> str:
    """Search ChEMBL for biological targets (proteins, organisms).

    Args:
        query: Target name, gene name, or ChEMBL target ID
        target_type: Optional filter: single_protein, protein_family, organism, etc.
        max_results: Max results (1-20, default 10)

    Returns:
        JSON array of targets with names, types, organism info
    """
    def _search():
        params = {"q": query, "limit": max(1, min(max_results, 20))}
        if target_type:
            params["target_type"] = target_type

        data = _get("target/search", params)
        if data is None:
            return []
        return [_clean_target(t) for t in data.get("targets", [])]

    results = await asyncio.to_thread(_search)
    return json.dumps(results, indent=2, ensure_ascii=False)


def _clean_target(t: dict) -> dict:
    comps = t.get("target_components", [])
    accessions = []
    gene_names = []
    for c in comps:
        accessions.extend(c.get("target_component_synonyms", []))
        if c.get("accession"):
            gene_names.append(c.get("accession"))

    return {
        "chembl_id": t.get("target_chembl_id"),
        "name": t.get("pref_name"),
        "type": t.get("target_type"),
        "organism": t.get("organism"),
        "accession": gene_names[0] if gene_names else None,
    }


# ─── Drug indications ─────────────────────────────────────────

@mcp.tool()
async def chembl_get_drug_indications(chembl_id: str, max_results: int = 20) -> str:
    """Get drug indications (diseases/conditions a drug is approved for).

    Args:
        chembl_id: ChEMBL molecule ID (e.g. 'CHEMBL25' for aspirin)
        max_results: Max results (1-50, default 20)

    Returns:
        JSON array of indication records with disease terms and mesh IDs
    """
    def _get():
        data = _get("drug_indication", {
            "molecule_chembl_id": chembl_id,
            "limit": max(1, min(max_results, 50)),
        })
        if data is None:
            return []
        return [
            {
                "indication": d.get("indication"),
                "mesh_heading": d.get("mesh_heading"),
                "mesh_id": d.get("mesh_id"),
                "efo_term": d.get("efo_term"),
                "efo_id": d.get("efo_id"),
                "max_phase_for_ind": d.get("max_phase_for_ind"),
            }
            for d in data.get("drug_indications", [])
        ]

    results = await asyncio.to_thread(_get)
    return json.dumps(results, indent=2, ensure_ascii=False)
