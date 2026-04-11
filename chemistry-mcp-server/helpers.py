"""Shared helpers for chemistry MCP tools."""

from rdkit import Chem


def mol_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def compound_to_dict(c) -> dict:
    """PubChemPy Compound -> clean dict."""
    d = {
        "cid": c.cid,
        "molecular_formula": getattr(c, "molecular_formula", None),
        "molecular_weight": getattr(c, "molecular_weight", None),
        "canonical_smiles": getattr(c, "canonical_smiles", None),
        "isomeric_smiles": getattr(c, "isomeric_smiles", None),
        "iupac_name": getattr(c, "iupac_name", None),
        "inchi": getattr(c, "inchi", None),
        "inchikey": getattr(c, "inchikey", None),
        "xlogp": getattr(c, "xlogp", None),
        "exact_mass": getattr(c, "exact_mass", None),
        "tpsa": getattr(c, "tpsa", None),
        "complexity": getattr(c, "complexity", None),
        "charge": getattr(c, "charge", None),
        "h_bond_donor_count": getattr(c, "h_bond_donor_count", None),
        "h_bond_acceptor_count": getattr(c, "h_bond_acceptor_count", None),
        "rotatable_bond_count": getattr(c, "rotatable_bond_count", None),
        "heavy_atom_count": getattr(c, "heavy_atom_count", None),
    }
    try:
        syns = c.synonyms
        if syns:
            d["synonyms"] = syns[:20]
    except Exception:
        pass
    return d


def resolve_namespace(query_type: str) -> str:
    return "cid" if query_type == "cid" else query_type
