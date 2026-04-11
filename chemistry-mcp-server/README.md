# Chemistry & Biology MCP Server

Unified MCP server for chemistry and biology tasks — PubChem/ChEMBL search, molecular analysis, protein data, metabolic pathways, equation balancing, format conversion, structure analysis, and visualization in one server.

Built by combining the best parts of [PubChem-MCP-Server](https://github.com/JackKuo666/PubChem-MCP-Server), [rdkit-mcp-server](https://github.com/tandemai-inc/rdkit-mcp-server), and [PubChempy_MCPserver](https://github.com/thinktraveller/PubChempy_MCPserver), plus **ChEMBL**, **UniProt**, and **KEGG** integration — 27 tools in one server.

## Stack

| Library | Purpose |
|---------|---------|
| `mcp` (FastMCP) | MCP server SDK |
| `rdkit` | Molecules, SMILES, descriptors, 2D/3D, drawing |
| `pubchempy` | PubChem REST API — search, properties, synonyms |
| `httpx` | ChEMBL, UniProt, KEGG REST APIs |
| `chempy` | Equation balancing, stoichiometry |
| `pillow` | Image generation |

## Installation

```bash
pip install mcp[cli] rdkit pubchempy chempy pillow httpx
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

### Claude Code (global)

```bash
claude mcp add --transport stdio --scope user chemistry -- python3 /absolute/path/to/chemistry-mcp-server/server.py
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

## Tools (27)

### 1. PubChem Search

#### `search_compound`
Search PubChem by name, SMILES, InChI, InChIKey, formula, or CID.

```
query: "aspirin"  |  query_type: "name"|"smiles"|"inchi"|"inchikey"|"formula"|"cid"  |  max_results: 1-20
```

#### `get_compound_properties`
Get specific properties from PubChem API.

```
query: "aspirin"  |  query_type: "name"  |  properties: ["MolecularWeight", "XLogP", "TPSA"]
```

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

28 available descriptors: `molecular_weight`, `logp`, `tpsa`, `h_bond_donors`, `h_bond_acceptors`, `rotatable_bonds`, `aromatic_rings`, `kappa1-3`, `chi0v-chi4v`, etc.

### 3. Format Conversion

#### `convert_format`
SMILES to canonical SMILES, InChI, and InChIKey.

```
smiles: "CCO"
```

#### `inchi_to_smiles`
InChI to SMILES.

```
inchi: "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
```

### 4. Equation Balancing

#### `balance_equation`
Balance chemical equations via chempy.

```
equation: "H2 + O2 -> H2O"
equation: "CH4 + O2 -> CO2 + H2O"
```

### 5. Structure Analysis

#### `check_substructure`
Check if a molecule contains a SMARTS pattern.

```
smiles: "CCO"  |  smarts: "[OX2H]"
```

#### `calculate_similarity`
Tanimoto similarity between two molecules (Morgan fingerprints).

```
smiles1: "CCO"  |  smiles2: "CCCO"  |  radius: 2  |  n_bits: 2048
```

#### `get_scaffold`
Extract Murcko scaffold.

```
smiles: "CC(=O)Oc1ccccc1C(=O)O"  |  generic: false
```

#### `fragment_molecule`
Fragment molecule for Matched Molecular Pair Analysis (MMPA).

```
smiles: "c1ccc(CC(=O)O)cc1"  |  max_cuts: 3
```

### 6. Visualization & 3D

#### `draw_molecule`
Draw 2D structure, returns base64-encoded PNG or SVG.

```
smiles: "c1ccccc1"  |  width: 400  |  height: 400  |  image_format: "png"|"svg"
```

#### `draw_molecule_grid`
Draw a grid of molecules, returns base64-encoded image.

```
smiles_list: ["CCO", "c1ccccc1", "CC(=O)O"]  |  legends: ["Ethanol", "Benzene", "Acetic acid"]
```

#### `generate_3d_structure`
Generate 3D conformers (ETKDGv3 + MMFF optimization).

```
smiles: "CCO"  |  num_conformers: 5  |  random_seed: 42
```

### 7. ChEMBL — Bioactivity & Drug Data

#### `chembl_search_molecule`
Search ChEMBL for molecules by name, SMILES, substructure, or similarity.

```
query: "aspirin"  |  search_type: "name"|"substructure"|"similarity"|"chembl_id"  |  max_results: 10
```

Returns: ChEMBL ID, name, SMILES, max phase, molecular properties.

#### `chembl_get_bioactivity`
Get bioactivity data (IC50, Ki, EC50, etc.) for a molecule.

```
chembl_id: "CHEMBL25"  |  target_chembl_id: optional  |  activity_type: "IC50"  |  max_results: 20
```

Returns: activity values, target info, assay descriptions, pChEMBL values.

#### `chembl_search_target`
Search ChEMBL for biological targets (proteins, organisms).

```
query: "EGFR"  |  target_type: optional  |  max_results: 10
```

#### `chembl_get_drug_indications`
Get approved disease indications for a drug.

```
chembl_id: "CHEMBL25"  |  max_results: 20
```

### 8. UniProt — Protein Data

#### `uniprot_search`
Search UniProt for proteins by name, gene, or keyword.

```
query: "EGFR"  |  organism: "human"  |  reviewed: true  |  max_results: 10
```

Returns: accession, protein name, gene name, organism, sequence length.

#### `uniprot_get_protein`
Get detailed protein information by accession ID.

```
accession: "P00533"
```

Returns: full sequence, function description, catalytic activity, subcellular locations, GO terms, cross-references (PDB, KEGG, ChEMBL, DrugBank).

#### `uniprot_map_ids`
Map protein identifiers between databases.

```
ids: ["P00533", "P04626"]  |  from_type: "UniProtKB_AC-ID"  |  to_type: "GeneID"
```

Supported: UniProtKB, GeneID, Ensembl, PDB, KEGGG, ChEMBL, DrugBank, RefSeq, and 100+ more.

### 9. KEGG — Pathways & Metabolism

#### `kegg_search`
Search KEGG databases for entries.

```
database: "compound"|"pathway"|"drug"|"disease"|"enzyme"  |  query: "glucose"  |  max_results: 20
```

#### `kegg_get_entry`
Get detailed information about a KEGG entry.

```
entry_id: "cpd:C00031"  |  "path:hsa00010"  |  "hsa:10458"  |  "D00075"
```

#### `kegg_list`
List entries in a KEGG database.

```
database: "pathway"  |  organism: "hsa"
```

#### `kegg_link`
Find cross-references between KEGG entries.

```
target_db: "pathway"  |  source_id: "hsa:10458"
```

#### `kegg_pathway_compounds`
Get all compounds involved in a KEGG pathway.

```
pathway_id: "hsa00010"
```

## Examples

```
> Search for aspirin on PubChem
search_compound(query="aspirin", query_type="name")

> Get molecular properties of ethanol
molecular_info(smiles="CCO")

> Balance combustion of methane
balance_equation(equation="CH4 + O2 -> CO2 + H2O")

> Find drugs targeting EGFR on ChEMBL
chembl_get_bioactivity(chembl_id="CHEMBL25", target_chembl_id="CHEMBL203")

> Get human EGFR protein details
uniprot_get_protein(accession="P00533")

> Find glucose in KEGG
kegg_search(database="compound", query="glucose")

> Get all compounds in glycolysis pathway
kegg_pathway_compounds(pathway_id="hsa00010")

> Map UniProt IDs to GeneIDs
uniprot_map_ids(ids=["P00533", "P04626"], from_type="UniProtKB_AC-ID", to_type="GeneID")

> Draw a molecule
draw_molecule(smiles="CC(=O)Oc1ccccc1C(=O)O")

> Generate 3D structure
generate_3d_structure(smiles="CCCCCC", num_conformers=5)
```

## Project Structure

```
chemistry-mcp-server/
├── server.py              # FastMCP instance + entry point
├── helpers.py             # Shared utilities (SMILES parsing, validation)
├── tools/
│   ├── __init__.py        # Imports all modules to register tools
│   ├── pubchem.py         # PubChem search (3 tools)
│   ├── properties.py      # RDKit descriptors + Lipinski (2 tools)
│   ├── conversion.py      # SMILES <-> InChI (2 tools)
│   ├── equations.py       # Equation balancing (1 tool)
│   ├── structure.py       # Substructure, similarity, scaffold, fragments (4 tools)
│   ├── visualization.py   # 2D drawing + 3D conformers (3 tools)
│   ├── chembl.py          # ChEMBL bioactivity & drug data (4 tools)
│   ├── uniprot.py         # UniProt protein search & metadata (3 tools)
│   └── kegg.py            # KEGG pathways & metabolism (5 tools)
├── tests/
│   ├── test_server.py     # 60 core tests
│   └── test_bugs.py       # 131 edge case + bug regression tests
├── pyproject.toml
└── README.md
```

## Testing

```bash
pip install pytest
python3 -m pytest tests/ -v
```

191 tests covering all tools including live PubChem, ChEMBL, UniProt, and KEGG API calls.

## Reference Servers

| Repository | Used for |
|-----------|----------|
| [PubChem-MCP-Server](https://github.com/JackKuo666/PubChem-MCP-Server) | PubChem search patterns |
| [rdkit-mcp-server](https://github.com/tandemai-inc/rdkit-mcp-server) | RDKit descriptor & drawing patterns |
| [PubChempy_MCPserver](https://github.com/thinktraveller/PubChempy_MCPserver) | PubChem + chempy patterns |
| [chembl_webresource_client](https://github.com/chembl/chembl_webresource_client) | ChEMBL API reference |
| [uniprot](https://github.com/boscoh/uniprot) | UniProt API reference |
| [KEGGRESTpy](https://github.com/guokai8/KEGGRESTpy) | KEGG API reference |

## License

MIT
