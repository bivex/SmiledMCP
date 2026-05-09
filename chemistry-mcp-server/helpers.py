"""Shared helpers for chemistry MCP tools."""

from typing import Optional

from rdkit import Chem


def _fix_aromatic_nitrogen(smiles: str) -> Optional[str]:
    """Try to fix bare aromatic nitrogen ('n') by replacing with [nH].

    RDKit's kekulizer cannot always determine where to place implicit hydrogens
    on aromatic nitrogen atoms in fused heterocycles (e.g. purines like adenine).
    This helper attempts to add explicit [nH] to bare aromatic N atoms so that
    kekulization succeeds.

    Returns a corrected SMILES string, or None if no fix was found.
    """
    # Find positions of bare lowercase 'n' not inside bracket atoms
    bare_n_positions = []
    in_bracket = False
    for i, c in enumerate(smiles):
        if c == "[":
            in_bracket = True
        elif c == "]":
            in_bracket = False
        elif c == "n" and not in_bracket:
            bare_n_positions.append(i)

    if not bare_n_positions:
        return None

    # Try replacing each bare 'n' with '[nH]' one at a time
    for pos in bare_n_positions:
        modified = smiles[:pos] + "[nH]" + smiles[pos + 1 :]
        if Chem.MolFromSmiles(modified) is not None:
            return modified

    # Try replacing all bare 'n' with '[nH]' at once
    modified = smiles
    offset = 0
    for pos in bare_n_positions:
        p = pos + offset
        modified = modified[:p] + "[nH]" + modified[p + 1 :]
        offset += 3  # '[nH]' is 3 chars longer than 'n'
    if Chem.MolFromSmiles(modified) is not None:
        return modified

    return None


def mol_from_smiles(smiles: str) -> Chem.Mol:
    """Parse SMILES with fallback aromatic nitrogen correction.

    Tries standard parsing first.  If it fails because the kekulizer
    cannot assign implicit hydrogens to aromatic nitrogens (common in
    fused heterocycles like purines), automatically attempts to insert
    explicit ``[nH]`` atoms and re-parses.

    Raises:
        ValueError: If the SMILES cannot be parsed even after correction.
    """
    if not smiles or not smiles.strip():
        raise ValueError(f"Invalid SMILES: {smiles!r} (empty)")

    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and mol.GetNumAtoms() > 0:
        return mol

    # Attempt to fix bare aromatic nitrogen atoms
    fixed = _fix_aromatic_nitrogen(smiles)
    if fixed is not None:
        mol = Chem.MolFromSmiles(fixed)
        if mol is not None and mol.GetNumAtoms() > 0:
            return mol

    # Fallback with sanitize=False to give a more informative error
    mol_nosani = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol_nosani is not None and mol_nosani.GetNumAtoms() > 0:
        raise ValueError(
            f"SMILES parsed without sanitization but failed full validation: {smiles!r}. "
            f"This usually means the structure has invalid valence or kekulization issues. "
            f"For aromatic nitrogen heterocycles, try using explicit [nH] notation."
        )

    raise ValueError(f"Invalid SMILES: {smiles!r} (RDKit could not parse this structure)")


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


VALID_QUERY_TYPES = {"name", "smiles", "inchi", "inchikey", "formula", "cid"}


def resolve_namespace(query_type: str) -> str:
    if query_type not in VALID_QUERY_TYPES:
        raise ValueError(
            f"Invalid query_type: {query_type!r}. "
            f"Must be one of: {sorted(VALID_QUERY_TYPES)}"
        )
    return "cid" if query_type == "cid" else query_type


def parse_cid(query: str) -> int:
    try:
        return int(query)
    except (ValueError, TypeError):
        raise ValueError(f"CID must be a number, got: {query!r}")
