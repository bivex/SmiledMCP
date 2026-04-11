"""Structure analysis tools — substructure, similarity, scaffold, fragmentation."""

import json

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

from server import mcp
from helpers import mol_from_smiles


@mcp.tool()
def check_substructure(smiles: str, smarts: str) -> str:
    """Check if a molecule contains a substructure defined by SMARTS pattern.

    Args:
        smiles: SMILES of the target molecule
        smarts: SMARTS pattern to search for (e.g. "[OX2H]" for hydroxyl, "c1ccccc1" for benzene)

    Returns:
        JSON with has_substructure (bool), match_count, and matching atom indices
    """
    mol = mol_from_smiles(smiles)
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        raise ValueError(f"Invalid SMARTS: {smarts}")

    has_match = mol.HasSubstructMatch(query)
    result = {
        "smiles": smiles,
        "smarts": smarts,
        "has_substructure": has_match,
    }

    if has_match:
        matches = mol.GetSubstructMatches(query)
        result["match_count"] = len(matches)
        result["matching_atoms"] = [list(m) for m in matches[:20]]

    return json.dumps(result, indent=2, ensure_ascii=False)


@mcp.tool()
def calculate_similarity(
    smiles1: str,
    smiles2: str,
    radius: int = 2,
    n_bits: int = 2048,
) -> str:
    """Calculate Tanimoto similarity between two molecules using Morgan fingerprints.

    Args:
        smiles1: SMILES of first molecule
        smiles2: SMILES of second molecule
        radius: Morgan fingerprint radius (default 2)
        n_bits: Fingerprint bit count (default 2048)

    Returns:
        JSON with tanimoto_similarity (0.0 to 1.0, where 1.0 = identical)
    """
    mol1 = mol_from_smiles(smiles1)
    mol2 = mol_from_smiles(smiles2)

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=n_bits)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=n_bits)

    sim = DataStructs.TanimotoSimilarity(fp1, fp2)

    return json.dumps({
        "smiles1": smiles1,
        "smiles2": smiles2,
        "tanimoto_similarity": round(sim, 4),
    }, indent=2)


@mcp.tool()
def get_scaffold(smiles: str, generic: bool = False) -> str:
    """Extract Murcko scaffold from a molecule.

    Args:
        smiles: SMILES string
        generic: If true, convert all atoms to C and bonds to single (generic scaffold)

    Returns:
        JSON with scaffold SMILES (and generic scaffold if requested)
    """
    mol = mol_from_smiles(smiles)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smi = Chem.MolToSmiles(scaffold)

    result = {"smiles": smiles, "scaffold": scaffold_smi}

    if generic:
        gen = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        result["generic_scaffold"] = Chem.MolToSmiles(gen)

    return json.dumps(result, indent=2)


@mcp.tool()
def fragment_molecule(smiles: str, max_cuts: int = 3) -> str:
    """Fragment a molecule for Matched Molecular Pair Analysis (MMPA).

    Args:
        smiles: SMILES string
        max_cuts: Maximum number of bond cuts (1-5, default 3)

    Returns:
        JSON with list of (core, side_chain) fragment pairs
    """
    from rdkit.Chem import rdMMPA

    mol = mol_from_smiles(smiles)
    max_cuts = max(1, min(max_cuts, 5))
    fragments = rdMMPA.FragmentMol(mol, maxCuts=max_cuts)

    frags = []
    for core, side in fragments:
        if core is not None and side is not None:
            frags.append({
                "core": Chem.MolToSmiles(core),
                "side_chains": Chem.MolToSmiles(side),
            })

    return json.dumps({
        "smiles": smiles,
        "fragment_count": len(frags),
        "fragments": frags[:50],
    }, indent=2, ensure_ascii=False)
