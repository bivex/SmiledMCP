"""Visualization & 3D tools — drawing molecules and generating conformers."""

import base64
import io
import json

from rdkit import Chem
from rdkit.Chem import Draw, AllChem

from server import mcp
from helpers import mol_from_smiles


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
    mol = mol_from_smiles(smiles)

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
    legends: list[str] | None = None,
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
    mol = mol_from_smiles(smiles)
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
