# Chemistry MCP Server

Unified MCP server for chemistry tasks — PubChem search, molecular analysis, equation balancing, format conversion, structure analysis, and visualization in one server.

Built by combining the best parts of [PubChem-MCP-Server](https://github.com/JackKuo666/PubChem-MCP-Server), [rdkit-mcp-server](https://github.com/tandemai-inc/rdkit-mcp-server), and [PubChempy_MCPserver](https://github.com/thinktraveller/PubChempy_MCPserver) into 15 clean tools.

## Stack

| Library | Purpose |
|---------|---------|
| `mcp` (FastMCP) | MCP server SDK |
| `rdkit` | Molecules, SMILES, descriptors, 2D/3D, drawing |
| `pubchempy` | PubChem REST API — search, properties, synonyms |
| `chempy` | Equation balancing, stoichiometry |
| `pillow` | Image generation |

## Installation

```bash
pip install mcp[cli] rdkit pubchempy chempy pillow
```

## Running

```bash
python server.py
```

The server runs on **stdio** transport by default.

## MCP Client Configuration

### Claude Desktop

Add to `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "chemistry": {
      "command": "python3",
      "args": ["/absolute/path/to/chemistry-mcp-server/server.py"]
    }
  }
}
```

### Cursor / Windsurf / Cline

```json
{
  "mcpServers": {
    "chemistry": {
      "command": "python3",
      "args": ["/absolute/path/to/chemistry-mcp-server/server.py"]
    }
  }
}
```

## Tools

### 1. PubChem Search

#### `search_compound`
Search PubChem by name, SMILES, InChI, InChIKey, formula, or CID.

```
query: "aspirin"  |  query_type: "name"|"smiles"|"inchi"|"inchikey"|"formula"|"cid"  |  max_results: 1-20
```

Returns: CID, molecular formula, weight, SMILES, InChI, InChIKey, LogP, TPSA, synonyms, etc.

#### `get_compound_properties`
Get specific properties from PubChem API.

```
query: "aspirin"  |  query_type: "name"  |  properties: ["MolecularWeight", "XLogP", "TPSA"]
```

Available properties: `MolecularFormula`, `MolecularWeight`, `CanonicalSMILES`, `IsomericSMILES`, `InChI`, `InChIKey`, `IUPACName`, `XLogP`, `ExactMass`, `TPSA`, `Complexity`, `HBondDonorCount`, `HBondAcceptorCount`, `RotatableBondCount`, `HeavyAtomCount`, `Volume3D`, etc.

#### `get_synonyms`
Get alternative names for a compound.

```
query: "aspirin"  |  query_type: "name"
```

### 2. Molecular Properties (RDKit)

#### `molecular_info`
Comprehensive molecular analysis from a SMILES string.

```
smiles: "CC(=O)Oc1ccccc1C(=O)O"
```

Returns: molecular weight, formula, LogP, TPSA, H-bond donors/acceptors, rotatable bonds, ring counts, fraction CSP3, and **Lipinski Rule of 5** check.

#### `compute_descriptors`
Compute named descriptors for multiple molecules in batch.

```
smiles_list: ["CCO", "c1ccccc1"]  |  descriptor_names: ["molecular_weight", "logp", "tpsa"]
```

Available descriptors (27): `molecular_weight`, `exact_molecular_weight`, `logp`, `tpsa`, `h_bond_donors`, `h_bond_acceptors`, `rotatable_bonds`, `aromatic_rings`, `aliphatic_rings`, `total_rings`, `heavy_atoms`, `num_atoms`, `fraction_csp3`, `heteroatoms`, `num_valence_electrons`, `num_radical_electrons`, `amide_bonds`, `lipinski_hba`, `lipinski_hbd`, `labute_asa`, `chi0v`–`chi4v`, `kappa1`–`kappa3`.

### 3. Format Conversion

#### `convert_format`
SMILES → canonical SMILES, InChI, InChIKey.

```
smiles: "CCO"
```

#### `inchi_to_smiles`
InChI → SMILES.

```
inchi: "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
```

### 4. Equation Balancing

#### `balance_equation`
Balance chemical equations via chempy.

```
equation: "H2 + O2 -> H2O"
equation: "CH4 + O2 -> CO2 + H2O"
equation: "Fe + O2 -> Fe2O3"
```

Returns: balanced equation string + reactant/product coefficients.

### 5. Structure Analysis

#### `check_substructure`
Check if a molecule contains a SMARTS pattern.

```
smiles: "CCO"  |  smarts: "[OX2H]"         → hydroxyl check
smiles: "c1ccccc1C(=O)O"  |  smarts: "c1ccccc1"  → benzene ring check
```

#### `calculate_similarity`
Tanimoto similarity between two molecules (Morgan fingerprints).

```
smiles1: "CCO"  |  smiles2: "CCCO"  |  radius: 2  |  n_bits: 2048
```

Returns: similarity score 0.0–1.0 (1.0 = identical).

#### `get_scaffold`
Extract Murcko scaffold.

```
smiles: "CC(=O)Oc1ccccc1C(=O)O"  |  generic: false
```

Set `generic: true` to get the generic scaffold (all atoms → C, all bonds → single).

#### `fragment_molecule`
Fragment molecule for Matched Molecular Pair Analysis (MMPA).

```
smiles: "c1ccc(CC(=O)O)cc1"  |  max_cuts: 3
```

### 6. Visualization & 3D

#### `draw_molecule`
Draw 2D structure → base64-encoded PNG or SVG.

```
smiles: "c1ccccc1"  |  width: 400  |  height: 400  |  image_format: "png"|"svg"
```

#### `draw_molecule_grid`
Draw a grid of molecules → base64-encoded image.

```
smiles_list: ["CCO", "c1ccccc1", "CC(=O)O"]  |  legends: ["Ethanol", "Benzene", "Acetic acid"]
```

#### `generate_3d_structure`
Generate 3D conformers (ETKDGv3 + MMFF optimization).

```
smiles: "CCO"  |  num_conformers: 5  |  random_seed: 42
```

Returns: SDF mol block for each conformer.

## Examples

```
> Search for aspirin on PubChem
search_compound(query="aspirin", query_type="name")

> Get molecular properties of ethanol
molecular_info(smiles="CCO")

> Balance combustion of methane
balance_equation(equation="CH4 + O2 -> CO2 + H2O")

> Convert SMILES to InChI
convert_format(smiles="c1ccccc1")

> Check if aspirin has a benzene ring
check_substructure(smiles="CC(=O)Oc1ccccc1C(=O)O", smarts="c1ccccc1")

> Compare two molecules
calculate_similarity(smiles1="CCO", smiles2="CCCO")

> Draw a molecule
draw_molecule(smiles="CC(=O)Oc1ccccc1C(=O)O")

> Generate 3D structure
generate_3d_structure(smiles="CCCCCC", num_conformers=5)
```

## Testing

```bash
pip install pytest
python3 -m pytest tests/ -v
```

60 tests covering all tools including live PubChem API calls.

## Comparison with original servers

| | PubChem-MCP-Server | rdkit-mcp-server | PubChempy_MCPserver | **This server** |
|---|---|---|---|---|
| Tools | 4 | ~72 | 43 | **15** |
| PubChem search | Yes | No | Yes | **Yes** |
| RDKit descriptors | No | Yes | No | **Yes** |
| Equation balancing | No | No | No | **Yes** |
| Format conversion | No | No | Partial | **Yes** |
| Lipinski check | No | No | No | **Yes** |
| 2D drawing | No | Yes | No | **Yes** |
| 3D conformers | No | Yes | No | **Yes** |
| Input format | Mixed | PickledMol | Plain | **SMILES/names** |

## License

MIT
