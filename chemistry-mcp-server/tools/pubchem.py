"""PubChem search tools — compound search, properties, synonyms."""

import asyncio
import json

import pubchempy as pcp

from server import mcp
from helpers import compound_to_dict, resolve_namespace, parse_cid


@mcp.tool()
async def search_compound(
    query: str,
    query_type: str = "name",
    max_results: int = 5,
) -> str:
    """Search PubChem for chemical compounds by name, SMILES, InChI, InChIKey, formula, or CID.

    Args:
        query: Search query string (or numeric CID)
        query_type: One of: name, smiles, inchi, inchikey, formula, cid
        max_results: Max results to return (1-20, default 5)

    Returns:
        JSON array of compound records with CID, SMILES, InChI, molecular weight, formula, etc.
    """
    resolve_namespace(query_type)  # validate before async
    capped = max(1, min(max_results, 20))

    def _search():
        if query_type == "cid":
            comp = pcp.Compound.from_cid(parse_cid(query))
            return [compound_to_dict(comp)]
        else:
            comps = pcp.get_compounds(query, query_type, max_records=capped)
            return [compound_to_dict(c) for c in comps]

    results = await asyncio.to_thread(_search)
    return json.dumps(results, indent=2, ensure_ascii=False)


@mcp.tool()
async def get_compound_properties(
    query: str,
    query_type: str = "name",
    properties: list[str] | None = None,
) -> str:
    """Get specific properties for compounds from PubChem API.

    Args:
        query: Compound identifier
        query_type: One of: name, smiles, inchi, inchikey, formula, cid
        properties: Property names to retrieve. Available: MolecularFormula, MolecularWeight,
            CanonicalSMILES, IsomericSMILES, InChI, InChIKey, IUPACName, XLogP, ExactMass,
            MonoisotopicMass, TPSA, Complexity, Charge, HBondDonorCount, HBondAcceptorCount,
            RotatableBondCount, HeavyAtomCount, IsotopeAtomCount, AtomStereoCount,
            DefinedAtomStereoCount, UndefinedAtomStereoCount, BondStereoCount,
            CovalentUnitCount, Volume3D.  If omitted, returns common properties.

    Returns:
        JSON array of property dictionaries
    """
    ns = resolve_namespace(query_type)
    ident = parse_cid(query) if query_type == "cid" else query

    def _get():
        if properties is None:
            props = [
                "MolecularFormula", "MolecularWeight", "CanonicalSMILES",
                "IsomericSMILES", "InChI", "InChIKey", "IUPACName",
                "XLogP", "ExactMass", "TPSA", "Complexity",
                "HBondDonorCount", "HBondAcceptorCount",
                "RotatableBondCount", "HeavyAtomCount",
            ]
        else:
            props = properties

        results = pcp.get_properties(props, ident, ns)
        return results if results else []

    results = await asyncio.to_thread(_get)
    return json.dumps(results, indent=2, ensure_ascii=False)


@mcp.tool()
async def get_synonyms(query: str, query_type: str = "name") -> str:
    """Get synonyms (alternative names) for a compound from PubChem.

    Args:
        query: Compound identifier
        query_type: One of: name, smiles, inchi, inchikey, formula, cid

    Returns:
        JSON list of synonym strings
    """
    ns = resolve_namespace(query_type)
    ident = parse_cid(query) if query_type == "cid" else query

    def _get():
        results = pcp.get_synonyms(ident, ns)
        if results and len(results) > 0:
            return results[0].get("Synonym", [])[:50]
        return []

    results = await asyncio.to_thread(_get)
    return json.dumps(results, indent=2, ensure_ascii=False)
