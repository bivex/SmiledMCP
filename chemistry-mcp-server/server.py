"""
Chemistry MCP Server — unified chemistry tools for LLMs.

Combines PubChem search (pubchempy), molecular analysis (RDKit),
equation balancing (chempy), format conversion, structure analysis,
and visualization into one clean MCP server.
"""

import asyncio
import base64
import io
import json
from typing import Optional

from mcp.server.fastmcp import FastMCP
import pubchempy as pcp
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw, AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

mcp = FastMCP("chemistry-mcp-server")

# Optional: chempy for equation balancing
try:
    from chempy import balance_stoichiometry
    HAS_CHEMPY = True
except ImportError:
    HAS_CHEMPY = False


# ════════════════════════════════════════════════════════════
#  Helpers
# ════════════════════════════════════════════════════════════

def _mol_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def _compound_to_dict(c) -> dict:
    """PubChemPy Compound → clean dict."""
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


def _resolve_namespace(query_type: str) -> str:
    return "cid" if query_type == "cid" else query_type


# ════════════════════════════════════════════════════════════
#  1. PubChem Search
# ════════════════════════════════════════════════════════════

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
    capped = max(1, min(max_results, 20))

    def _search():
        if query_type == "cid":
            comp = pcp.Compound.from_cid(int(query))
            return [_compound_to_dict(comp)]
        else:
            comps = pcp.get_compounds(query, query_type, max_records=capped)
            return [_compound_to_dict(c) for c in comps]

    results = await asyncio.to_thread(_search)
    return json.dumps(results, indent=2, ensure_ascii=False)


@mcp.tool()
async def get_compound_properties(
    query: str,
    query_type: str = "name",
    properties: Optional[list[str]] = None,
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

        ns = _resolve_namespace(query_type)
        ident = int(query) if query_type == "cid" else query
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
    def _get():
        ns = _resolve_namespace(query_type)
        ident = int(query) if query_type == "cid" else query
        results = pcp.get_synonyms(ident, ns)
        if results and len(results) > 0:
            return results[0].get("Synonym", [])[:50]
        return []

    results = await asyncio.to_thread(_get)
    return json.dumps(results, indent=2, ensure_ascii=False)


# ════════════════════════════════════════════════════════════
#  2. Molecular Properties (RDKit)
# ════════════════════════════════════════════════════════════

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
    mol = _mol_from_smiles(smiles)

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
        "num_atoms": mol.GetNumAtoms(),
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
        JSON array — each element is a dict mapping descriptor name → value for one molecule
    """
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
        "num_atoms": lambda m: m.GetNumAtoms(),
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
        mol = _mol_from_smiles(smi)
        row = {"smiles": smi}
        for name, func in funcs.items():
            try:
                row[name] = func(mol)
            except Exception as e:
                row[name] = f"error: {e}"
        results.append(row)

    return json.dumps(results, indent=2, ensure_ascii=False)


# ════════════════════════════════════════════════════════════
#  3. Format Conversion
# ════════════════════════════════════════════════════════════

@mcp.tool()
def convert_format(smiles: str) -> str:
    """Convert SMILES to canonical SMILES, InChI, and InChIKey.

    Args:
        smiles: Input SMILES string

    Returns:
        JSON with canonical_smiles, inchi, inchikey
    """
    mol = _mol_from_smiles(smiles)
    result = {
        "input_smiles": smiles,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "inchi": Chem.MolToInchi(mol),
        "inchikey": Chem.MolToInchiKey(mol),
    }
    return json.dumps(result, indent=2, ensure_ascii=False)


@mcp.tool()
def inchi_to_smiles(inchi: str) -> str:
    """Convert InChI string to SMILES format.

    Args:
        inchi: InChI string

    Returns:
        JSON with smiles, canonical_smiles, inchikey
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


# ════════════════════════════════════════════════════════════
#  4. Equation Balancing
# ════════════════════════════════════════════════════════════

@mcp.tool()
def balance_equation(equation: str) -> str:
    """Balance a chemical equation.

    Uses chempy for stoichiometric balancing. Examples:
      "H2 + O2 -> H2O"
      "Fe + O2 -> Fe2O3"
      "CH4 + O2 -> CO2 + H2O"

    Args:
        equation: Chemical equation with '->' or '=' separator

    Returns:
        JSON with balanced equation and stoichiometric coefficients
    """
    if not HAS_CHEMPY:
        return json.dumps({
            "error": "chempy is not installed. Run: pip install chempy"
        })

    try:
        if "->" in equation:
            left_str, right_str = equation.split("->", 1)
        elif "=" in equation:
            left_str, right_str = equation.split("=", 1)
        else:
            return json.dumps({"error": "Use '->' or '=' to separate reactants and products"})

        def parse_substances(s):
            return [part.strip() for part in s.split("+") if part.strip()]

        reactants = parse_substances(left_str)
        products = parse_substances(right_str)

        balanced_r, balanced_p = balance_stoichiometry(reactants, products)

        def format_side(coeffs):
            parts = []
            for substance, coeff in coeffs.items():
                parts.append(f"{coeff} {substance}" if coeff != 1 else str(substance))
            return " + ".join(parts)

        balanced_eq = f"{format_side(balanced_r)} -> {format_side(balanced_p)}"

        result = {
            "original": equation,
            "balanced": balanced_eq,
            "reactant_coefficients": {str(k): int(v) for k, v in balanced_r.items()},
            "product_coefficients": {str(k): int(v) for k, v in balanced_p.items()},
        }
        return json.dumps(result, indent=2, ensure_ascii=False)

    except Exception as e:
        return json.dumps({"error": f"Failed to balance: {str(e)}"})


# ════════════════════════════════════════════════════════════
#  5. Structure Analysis
# ════════════════════════════════════════════════════════════

@mcp.tool()
def check_substructure(smiles: str, smarts: str) -> str:
    """Check if a molecule contains a substructure defined by SMARTS pattern.

    Args:
        smiles: SMILES of the target molecule
        smarts: SMARTS pattern to search for (e.g. "[OX2H]" for hydroxyl, "c1ccccc1" for benzene)

    Returns:
        JSON with has_substructure (bool), match_count, and matching atom indices
    """
    mol = _mol_from_smiles(smiles)
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
    mol1 = _mol_from_smiles(smiles1)
    mol2 = _mol_from_smiles(smiles2)

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
    mol = _mol_from_smiles(smiles)
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

    mol = _mol_from_smiles(smiles)
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


# ════════════════════════════════════════════════════════════
#  6. Visualization & 3D
# ════════════════════════════════════════════════════════════

@mcp.tool()
def draw_molecule(
    smiles: str,
    width: int = 400,
    height: int = 400,
    image_format: str = "png",
) -> str:
    """Draw a 2D molecular structure. Returns base64-encoded image data.

    Args:
        smiles: SMILES string
        width: Image width in pixels (default 400)
        height: Image height in pixels (default 400)
        image_format: "png" or "svg"

    Returns:
        Base64-encoded image data (decode and save to view the structure)
    """
    mol = _mol_from_smiles(smiles)

    if image_format.lower() == "svg":
        from rdkit.Chem.Draw import rdMolDraw2D
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_data = drawer.GetDrawingText()
        b64 = base64.b64encode(svg_data.encode()).decode()
        return json.dumps({
            "format": "svg",
            "data": b64,
            "note": "Decode base64 to get SVG string",
        })
    else:
        img = Draw.MolToImage(mol, size=(width, height))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode()
        return json.dumps({
            "format": "png",
            "data": b64,
            "note": "Decode base64 and save as .png to view the structure",
        })


@mcp.tool()
def draw_molecule_grid(
    smiles_list: list[str],
    legends: Optional[list[str]] = None,
    sub_img_size: int = 300,
    image_format: str = "png",
) -> str:
    """Draw a grid of 2D molecular structures. Returns base64-encoded image.

    Args:
        smiles_list: List of SMILES strings (max 20)
        legends: Optional list of legend strings for each molecule
        sub_img_size: Size of each sub-image in pixels (default 300)
        image_format: "png" or "svg"

    Returns:
        Base64-encoded image data of the molecule grid
    """
    smiles_list = smiles_list[:20]
    mols = []
    valid_legends = []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)
            if legends and i < len(legends):
                valid_legends.append(legends[i])
            else:
                valid_legends.append(smi)

    if not mols:
        raise ValueError("No valid molecules found")

    if image_format.lower() == "svg":
        from rdkit.Chem.Draw import rdMolDraw2D
        n_cols = min(len(mols), 4)
        n_rows = (len(mols) + n_cols - 1) // n_cols
        drawer = rdMolDraw2D.MolDraw2DSVG(
            sub_img_size * n_cols, sub_img_size * n_rows,
            sub_img_size, sub_img_size,
        )
        drawer.DrawMolecules(mols, legends=valid_legends)
        drawer.FinishDrawing()
        svg_data = drawer.GetDrawingText()
        b64 = base64.b64encode(svg_data.encode()).decode()
        return json.dumps({"format": "svg", "data": b64})
    else:
        img = Draw.MolsToGridImage(
            mols, molsPerRow=4, subImgSize=(sub_img_size, sub_img_size),
            legends=valid_legends,
        )
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode()
        return json.dumps({"format": "png", "data": b64})


@mcp.tool()
def generate_3d_structure(smiles: str, num_conformers: int = 1, random_seed: int = 42) -> str:
    """Generate 3D molecular coordinates for a molecule.

    Uses RDKit's ETKDGv3 algorithm for 3D conformer generation.

    Args:
        smiles: SMILES string
        num_conformers: Number of conformers to generate (1-50, default 1)
        random_seed: Random seed for reproducibility

    Returns:
        JSON with 3D SDF (mol block) for each conformer
    """
    mol = _mol_from_smiles(smiles)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed

    num_conformers = max(1, min(num_conformers, 50))

    if num_conformers == 1:
        conf_id = AllChem.EmbedMolecule(mol, params)
        if conf_id == -1:
            return json.dumps({"error": "Failed to generate 3D coordinates"})
        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
        conformers = [{"id": conf_id, "sdf": Chem.MolToMolBlock(mol, confId=conf_id)}]
    else:
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)
        conformers = []
        for cid in conf_ids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid)
            conformers.append({
                "id": int(cid),
                "sdf": Chem.MolToMolBlock(mol, confId=cid),
            })

    return json.dumps({
        "smiles": smiles,
        "conformer_count": len(conformers),
        "conformers": conformers,
    }, indent=2, ensure_ascii=False)


# ════════════════════════════════════════════════════════════
#  Entry point
# ════════════════════════════════════════════════════════════

def main():
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
