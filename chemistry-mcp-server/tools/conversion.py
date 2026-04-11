"""Format conversion tools — SMILES <-> InChI."""

import json

from rdkit import Chem

from server import mcp
from helpers import mol_from_smiles


@mcp.tool()
def convert_format(smiles: str) -> str:
    """Convert SMILES to canonical SMILES, InChI, and InChIKey.

    Args:
        smiles: Input SMILES string

    Returns:
        JSON with canonical_smiles, inchi, inchikey
    """
    mol = mol_from_smiles(smiles)
    canonical = Chem.MolToSmiles(mol)
    inchi = Chem.MolToInchi(mol)
    if inchi is None:
        raise ValueError(f"Could not generate InChI for: {smiles}")
    result = {
        "input_smiles": smiles,
        "canonical_smiles": canonical,
        "inchi": inchi,
        "inchikey": Chem.MolToInchiKey(mol),
    }
    return json.dumps(result, indent=2, ensure_ascii=False)


@mcp.tool()
def inchi_to_smiles(inchi: str) -> str:
    """Convert InChI string to SMILES format.

    Args:
        inchi: InChI string

    Returns:
        JSON with smiles, inchikey
    """
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        raise ValueError(f"Invalid InChI: {inchi}")
    result = {
        "inchi": inchi,
        "smiles": Chem.MolToSmiles(mol),
        "inchikey": Chem.MolToInchiKey(mol),
    }
    return json.dumps(result, indent=2, ensure_ascii=False)
