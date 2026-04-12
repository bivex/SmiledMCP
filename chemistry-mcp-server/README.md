<h1 align="center">Chemistry & Biology MCP Server</h1>

<p align="center">
  <strong>27 tools. One server. Every chemistry & biology task your AI assistant needs.</strong>
</p>

<p align="center">
  Search PubChem & ChEMBL. Analyze molecules. Pull protein data. Walk metabolic pathways.<br>
  Balance equations. Generate 3D structures. All from your AI chat.
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+">
  <img src="https://img.shields.io/badge/tools-27-green" alt="27 Tools">
  <img src="https://img.shields.io/badge/tests-191%20passing-brightgreen" alt="191 Tests">
  <img src="https://img.shields.io/badge/license-MIT-yellow" alt="MIT License">
</p>

---

## What It Does

This is a unified [MCP server](https://modelcontextprotocol.io) that gives AI assistants like Claude, Cursor, and Windsurf direct access to chemistry and biology tools. Instead of copying data between websites and your AI, the AI calls these tools natively.

**Drug discovery** — Search ChEMBL for bioactivity data. Find protein targets. Check drug indications.

**Molecular analysis** — Compute descriptors. Check Lipinski's Rule of 5. Compare structures with Tanimoto similarity. Extract Murcko scaffolds.

**Chemical search** — Look up any compound on PubChem by name, SMILES, InChI, formula, or CID. Get properties and synonyms.

**Equations** — Balance any chemical equation instantly.

**Proteins** — Search UniProt. Get sequences, functions, subcellular locations. Map IDs across databases.

**Pathways** — Browse KEGG pathways. Find compounds in glycolysis. Link genes to reactions.

**Visualization** — Draw 2D structures. Generate 3D conformers. Export as PNG or SVG.

---

## Quick Start

```bash
pip install mcp[cli] rdkit pubchempy chempy pillow httpx
```

```bash
python server.py
```

That's it. The server starts on **stdio** transport by default.

---

## Connect to Your AI

### Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS)
or `%APPDATA%\Claude\claude_desktop_config.json` (Windows):

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

### Claude Code

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

---

## All 27 Tools

### PubChem — Compound Search

| Tool | What it does |
|------|-------------|
| `search_compound` | Search PubChem by name, SMILES, InChI, InChIKey, formula, or CID |
| `get_compound_properties` | Get specific properties — molecular weight, LogP, TPSA, and more |
| `get_synonyms` | Get all alternative names for a compound |

```python
search_compound(query="aspirin", query_type="name")
get_compound_properties(query="aspirin", properties=["MolecularWeight", "XLogP", "TPSA"])
get_synonyms(query="aspirin")
```

### RDKit — Molecular Properties

| Tool | What it does |
|------|-------------|
| `molecular_info` | Full molecular profile — weight, formula, LogP, TPSA, H-bonds, Lipinski Rule of 5 |
| `compute_descriptors` | Batch-compute named descriptors across multiple molecules |

```python
molecular_info(smiles="CC(=O)Oc1ccccc1C(=O)O")
compute_descriptors(smiles_list=["CCO", "c1ccccc1"], descriptor_names=["molecular_weight", "logp", "tpsa"])
```

28 available descriptors: `molecular_weight`, `logp`, `tpsa`, `h_bond_donors`, `h_bond_acceptors`, `rotatable_bonds`, `aromatic_rings`, `kappa1-3`, `chi0v-chi4v`, and more.

### Format Conversion

| Tool | What it does |
|------|-------------|
| `convert_format` | SMILES → canonical SMILES, InChI, InChIKey |
| `inchi_to_smiles` | InChI → SMILES |

```python
convert_format(smiles="CCO")
inchi_to_smiles(inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
```

### Equation Balancing

| Tool | What it does |
|------|-------------|
| `balance_equation` | Balance any chemical equation |

```python
balance_equation(equation="H2 + O2 -> H2O")
balance_equation(equation="CH4 + O2 -> CO2 + H2O")
```

### Structure Analysis

| Tool | What it does |
|------|-------------|
| `check_substructure` | Check if a molecule matches a SMARTS pattern |
| `calculate_similarity` | Tanimoto similarity between two molecules (Morgan fingerprints) |
| `get_scaffold` | Extract Murcko scaffold from a molecule |
| `fragment_molecule` | Fragment for Matched Molecular Pair Analysis (MMPA) |

```python
check_substructure(smiles="CCO", smarts="[OX2H]")
calculate_similarity(smiles1="CCO", smiles2="CCCO")
get_scaffold(smiles="CC(=O)Oc1ccccc1C(=O)O")
fragment_molecule(smiles="c1ccc(CC(=O)O)cc1", max_cuts=3)
```

### Visualization & 3D

| Tool | What it does |
|------|-------------|
| `draw_molecule` | Draw a 2D structure — returns PNG or SVG |
| `draw_molecule_grid` | Draw multiple molecules side-by-side |
| `generate_3d_structure` | Generate 3D conformers (ETKDGv3 + MMFF optimization) |

```python
draw_molecule(smiles="c1ccccc1", image_format="png")
draw_molecule_grid(smiles_list=["CCO", "c1ccccc1"], legends=["Ethanol", "Benzene"])
generate_3d_structure(smiles="CCCCCC", num_conformers=5)
```

### ChEMBL — Bioactivity & Drug Data

| Tool | What it does |
|------|-------------|
| `chembl_search_molecule` | Search molecules by name, SMILES, substructure, or similarity |
| `chembl_get_bioactivity` | Get IC50, Ki, EC50 data for a molecule against targets |
| `chembl_search_target` | Search for biological targets — proteins, organisms |
| `chembl_get_drug_indications` | Get approved disease indications for a drug |

```python
chembl_search_molecule(query="aspirin", search_type="name")
chembl_get_bioactivity(chembl_id="CHEMBL25", activity_type="IC50")
chembl_search_target(query="EGFR")
chembl_get_drug_indications(chembl_id="CHEMBL25")
```

### UniProt — Protein Data

| Tool | What it does |
|------|-------------|
| `uniprot_search` | Search proteins by name, gene, or keyword |
| `uniprot_get_protein` | Get full details — sequence, function, locations, cross-references |
| `uniprot_map_ids` | Map identifiers across 100+ databases (GeneID, Ensembl, PDB, KEGG, ...) |

```python
uniprot_search(query="EGFR", organism="human")
uniprot_get_protein(accession="P00533")
uniprot_map_ids(ids=["P00533", "P04626"], from_type="UniProtKB_AC-ID", to_type="GeneID")
```

### KEGG — Pathways & Metabolism

| Tool | What it does |
|------|-------------|
| `kegg_search` | Search KEGG databases — compounds, pathways, drugs, diseases, enzymes |
| `kegg_get_entry` | Get full details for any KEGG entry |
| `kegg_list` | List all entries in a KEGG database |
| `kegg_link` | Find cross-references between KEGG entries |
| `kegg_pathway_compounds` | Get all compounds in a metabolic pathway |

```python
kegg_search(database="compound", query="glucose")
kegg_get_entry(entry_id="cpd:C00031")
kegg_list(database="pathway", organism="hsa")
kegg_link(target_db="pathway", source_id="hsa:10458")
kegg_pathway_compounds(pathway_id="hsa00010")
```

---

## Real-World Examples

Ask your AI assistant something like:

> "Find all drugs approved for hypertension and show me their molecular properties"

Behind the scenes, it chains `chembl_search_molecule` → `chembl_get_drug_indications` → `molecular_info`.

> "What compounds are involved in human glycolysis?"

Calls `kegg_pathway_compounds` for `hsa00010`.

> "Compare the Tanimoto similarity between aspirin and ibuprofen"

Calls `search_compound` for both → `calculate_similarity` on the SMILES.

> "Balance the combustion of octane"

Calls `balance_equation("C8H18 + O2 -> CO2 + H2O")`.

> "Get the full sequence and function of human EGFR"

Calls `uniprot_get_protein(accession="P00533")`.

---

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

---

## Testing

```bash
pip install pytest
python3 -m pytest tests/ -v
```

191 tests covering every tool — including live API calls to PubChem, ChEMBL, UniProt, and KEGG.

---

## Tech Stack

| Library | Role |
|---------|------|
| `mcp` (FastMCP) | MCP server SDK |
| `rdkit` | Molecules, SMILES, descriptors, 2D/3D, drawing |
| `pubchempy` | PubChem REST API — search, properties, synonyms |
| `httpx` | ChEMBL, UniProt, KEGG REST APIs |
| `chempy` | Equation balancing and stoichiometry |
| `pillow` | Image generation |

---

## Credits

Built by combining the best patterns from:

- [PubChem-MCP-Server](https://github.com/JackKuo666/PubChem-MCP-Server) — PubChem search
- [rdkit-mcp-server](https://github.com/tandemai-inc/rdkit-mcp-server) — RDKit descriptors & drawing
- [PubChempy_MCPserver](https://github.com/thinktraveller/PubChempy_MCPserver) — PubChem + chempy
- [chembl_webresource_client](https://github.com/chembl/chembl_webresource_client) — ChEMBL API reference
- [uniprot](https://github.com/boscoh/uniprot) — UniProt API reference
- [KEGGRESTpy](https://github.com/guokai8/KEGGRESTpy) — KEGG API reference

## License

MIT
