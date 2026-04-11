"""Molecular property tools — RDKit descriptors and Lipinski check."""

import json

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from server import mcp
from helpers import mol_from_smiles


@mcp.tool()
def molecular_info(smiles: str) -> str:
    """Calculate comprehensive molecular properties from SMILES using RDKit.

    Returns: molecular weight, formula, LogP, TPSA, H-bond donors/acceptors,
    rotatable bonds, ring counts, fraction CSP3, Lipinski rule check, and more.

    Args:
        smiles: SMILES string of the molecule

    Returns:
        JSON dictionary with all calculated molecular properties
    """
    mol = mol_from_smiles(smiles)

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    info = {
        "smiles": Chem.MolToSmiles(mol),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(mw, 4),
        "exact_molecular_weight": round(Descriptors.ExactMolWt(mol), 6),
        "logp": round(logp, 4),
        "tpsa": round(Descriptors.TPSA(mol), 4),
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba,
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "aromatic_rings": Descriptors.NumAromaticRings(mol),
        "aliphatic_rings": Descriptors.NumAliphaticRings(mol),
        "total_rings": Descriptors.RingCount(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "atom_count": mol.GetNumAtoms(),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 4),
        "num_valence_electrons": Descriptors.NumValenceElectrons(mol),
        "num_radical_electrons": Descriptors.NumRadicalElectrons(mol),
        "heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "amide_bonds": rdMolDescriptors.CalcNumAmideBonds(mol),
        "lipinski_rule_of_5": {
            "mw_le_500": mw <= 500,
            "logp_le_5": logp <= 5,
            "hbd_le_5": hbd <= 5,
            "hba_le_10": hba <= 10,
            "violations": sum([mw > 500, logp > 5, hbd > 5, hba > 10]),
            "passes": all([mw <= 500, logp <= 5, hbd <= 5, hba <= 10]),
        },
    }

    return json.dumps(info, indent=2, ensure_ascii=False)


DESCRIPTOR_MAP = {
    "molecular_weight": lambda m: round(Descriptors.MolWt(m), 4),
    "exact_molecular_weight": lambda m: round(Descriptors.ExactMolWt(m), 6),
    "logp": lambda m: round(Descriptors.MolLogP(m), 4),
    "tpsa": lambda m: round(Descriptors.TPSA(m), 4),
    "h_bond_donors": lambda m: Descriptors.NumHDonors(m),
    "h_bond_acceptors": lambda m: Descriptors.NumHAcceptors(m),
    "rotatable_bonds": lambda m: Descriptors.NumRotatableBonds(m),
    "aromatic_rings": lambda m: Descriptors.NumAromaticRings(m),
    "aliphatic_rings": lambda m: Descriptors.NumAliphaticRings(m),
    "total_rings": lambda m: Descriptors.RingCount(m),
    "heavy_atoms": lambda m: m.GetNumHeavyAtoms(),
    "atom_count": lambda m: m.GetNumAtoms(),
    "fraction_csp3": lambda m: round(rdMolDescriptors.CalcFractionCSP3(m), 4),
    "heteroatoms": lambda m: rdMolDescriptors.CalcNumHeteroatoms(m),
    "num_valence_electrons": lambda m: Descriptors.NumValenceElectrons(m),
    "num_radical_electrons": lambda m: Descriptors.NumRadicalElectrons(m),
    "amide_bonds": lambda m: rdMolDescriptors.CalcNumAmideBonds(m),
    "lipinski_hba": lambda m: rdMolDescriptors.CalcNumLipinskiHBA(m),
    "lipinski_hbd": lambda m: rdMolDescriptors.CalcNumLipinskiHBD(m),
    "labute_asa": lambda m: round(rdMolDescriptors.CalcLabuteASA(m), 4),
    "chi0v": lambda m: round(rdMolDescriptors.CalcChi0v(m), 4),
    "chi1v": lambda m: round(rdMolDescriptors.CalcChi1v(m), 4),
    "chi2v": lambda m: round(rdMolDescriptors.CalcChi2v(m), 4),
    "chi3v": lambda m: round(rdMolDescriptors.CalcChi3v(m), 4),
    "chi4v": lambda m: round(rdMolDescriptors.CalcChi4v(m), 4),
    "kappa1": lambda m: round(rdMolDescriptors.CalcKappa1(m), 4),
    "kappa2": lambda m: round(rdMolDescriptors.CalcKappa2(m), 4),
    "kappa3": lambda m: round(rdMolDescriptors.CalcKappa3(m), 4),
}


@mcp.tool()
def compute_descriptors(
    smiles_list: list[str],
    descriptor_names: list[str],
) -> str:
    """Compute specific molecular descriptors for multiple molecules.

    Available descriptor names:
        molecular_weight, exact_molecular_weight, logp, tpsa,
        h_bond_donors, h_bond_acceptors, rotatable_bonds,
        aromatic_rings, aliphatic_rings, total_rings,
        heavy_atoms, num_atoms, fraction_csp3, heteroatoms,
        num_valence_electrons, num_radical_electrons, amide_bonds,
        lipinski_hba, lipinski_hbd, labute_asa,
        chi0v, chi1v, chi2v, chi3v, chi4v,
        kappa1, kappa2, kappa3

    Args:
        smiles_list: List of SMILES strings (max 50)
        descriptor_names: List of descriptor names to compute

    Returns:
        JSON array — each element is a dict mapping descriptor name -> value for one molecule
    """
    unknown = [d for d in descriptor_names if d not in DESCRIPTOR_MAP]
    if unknown:
        raise ValueError(
            f"Unknown descriptors: {unknown}. "
            f"Available: {list(DESCRIPTOR_MAP.keys())}"
        )

    smiles_list = smiles_list[:50]
    funcs = {name: DESCRIPTOR_MAP[name] for name in descriptor_names}

    results = []
    for smi in smiles_list:
        mol = mol_from_smiles(smi)
        row = {"smiles": smi}
        for name, func in funcs.items():
            try:
                row[name] = func(mol)
            except Exception as e:
                row[name] = f"error: {e}"
        results.append(row)

    return json.dumps(results, indent=2, ensure_ascii=False)
