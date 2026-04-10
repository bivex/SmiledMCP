"""Tests for chemistry MCP server."""

import json
import sys
import os
import asyncio

import pytest

# Make sure server.py is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from server import (
    search_compound,
    get_compound_properties,
    get_synonyms,
    molecular_info,
    compute_descriptors,
    convert_format,
    inchi_to_smiles,
    balance_equation,
    check_substructure,
    calculate_similarity,
    get_scaffold,
    fragment_molecule,
    draw_molecule,
    draw_molecule_grid,
    generate_3d_structure,
    _mol_from_smiles,
)


def _json(result: str) -> dict | list:
    """Parse JSON string returned by tools."""
    return json.loads(result)


# ════════════════════════════════════════════════════════════
#  Helpers
# ════════════════════════════════════════════════════════════

class TestHelpers:
    def test_mol_from_smiles_valid(self):
        mol = _mol_from_smiles("CCO")
        assert mol is not None
        assert mol.GetNumHeavyAtoms() == 3

    def test_mol_from_smiles_invalid(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            _mol_from_smiles("not_a_molecule_XYZZZ")


# ════════════════════════════════════════════════════════════
#  1. PubChem Search (async — need to await)
# ════════════════════════════════════════════════════════════

class TestSearchCompound:
    """PubChem search tools — make real API calls."""

    def test_search_by_name(self):
        result = _json(asyncio.run(search_compound("aspirin", "name", 1)))
        assert isinstance(result, list)
        assert len(result) >= 1
        comp = result[0]
        assert comp["cid"] is not None
        assert comp["canonical_smiles"] is not None
        assert comp["molecular_formula"] is not None

    def test_search_by_name_multiple(self):
        result = _json(asyncio.run(search_compound("benzene", "name", 3)))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_search_by_cid(self):
        result = _json(asyncio.run(search_compound("2244", "cid")))
        assert isinstance(result, list)
        assert len(result) == 1
        assert result[0]["cid"] == 2244  # Aspirin

    def test_search_by_smiles(self):
        result = _json(asyncio.run(search_compound("CC(=O)Oc1ccccc1C(=O)O", "smiles", 1)))
        assert isinstance(result, list)
        assert len(result) >= 1
        assert result[0]["cid"] == 2244

    def test_search_not_found(self):
        result = _json(asyncio.run(search_compound("XYZNOTAREALCOMPOUND12345", "name", 1)))
        assert isinstance(result, list)
        assert len(result) == 0


class TestGetCompoundProperties:
    def test_default_properties(self):
        result = _json(asyncio.run(get_compound_properties("aspirin", "name")))
        assert isinstance(result, list)
        assert len(result) >= 1
        prop = result[0]
        assert "MolecularFormula" in prop
        assert "MolecularWeight" in prop
        assert "CID" in prop

    def test_specific_properties(self):
        result = _json(asyncio.run(
            get_compound_properties("aspirin", "name", ["MolecularWeight", "XLogP"])
        ))
        assert isinstance(result, list)
        assert len(result) >= 1
        prop = result[0]
        assert "MolecularWeight" in prop
        assert "XLogP" in prop

    def test_by_cid(self):
        result = _json(asyncio.run(get_compound_properties("2244", "cid")))
        assert isinstance(result, list)
        assert len(result) >= 1


class TestGetSynonyms:
    def test_synonyms_aspirin(self):
        result = _json(asyncio.run(get_synonyms("aspirin", "name")))
        assert isinstance(result, list)
        assert len(result) > 0
        # Aspirin should have many synonyms
        assert any("aspirin" in s.lower() for s in result)

    def test_synonyms_by_cid(self):
        result = _json(asyncio.run(get_synonyms("2244", "cid")))
        assert isinstance(result, list)
        assert len(result) > 0


# ════════════════════════════════════════════════════════════
#  2. Molecular Properties (RDKit)
# ════════════════════════════════════════════════════════════

class TestMolecularInfo:
    def test_aspirin(self):
        info = _json(molecular_info("CC(=O)Oc1ccccc1C(=O)O"))
        assert info["molecular_formula"] == "C9H8O4"
        assert abs(info["molecular_weight"] - 180.157) < 0.1
        assert info["h_bond_donors"] == 1
        assert info["h_bond_acceptors"] == 3
        assert info["aromatic_rings"] == 1
        assert info["lipinski_rule_of_5"]["passes"] is True

    def test_ethanol(self):
        info = _json(molecular_info("CCO"))
        assert info["molecular_formula"] == "C2H6O"
        assert abs(info["molecular_weight"] - 46.069) < 0.1
        assert info["h_bond_donors"] == 1
        assert info["h_bond_acceptors"] == 1
        assert info["lipinski_rule_of_5"]["passes"] is True

    def test_large_molecule_lipinski_fail(self):
        # A large molecule that likely fails Lipinski
        info = _json(molecular_info("C1CCCCC1C(CCCCCCCCCCCC)CCCCCCCCCCCCCCCCCCCCCCCC"))
        lipinski = info["lipinski_rule_of_5"]
        assert lipinski["violations"] > 0
        assert lipinski["passes"] is False

    def test_benzene(self):
        info = _json(molecular_info("c1ccccc1"))
        assert info["molecular_formula"] == "C6H6"
        assert info["aromatic_rings"] == 1
        assert info["total_rings"] == 1
        assert info["rotatable_bonds"] == 0

    def test_invalid_smiles(self):
        with pytest.raises(ValueError):
            molecular_info("INVALID_SMILES")


class TestComputeDescriptors:
    def test_single_molecule(self):
        result = _json(compute_descriptors(
            ["CCO"],
            ["molecular_weight", "logp", "h_bond_donors"],
        ))
        assert isinstance(result, list)
        assert len(result) == 1
        row = result[0]
        assert row["smiles"] == "CCO"
        assert "molecular_weight" in row
        assert "logp" in row
        assert "h_bond_donors" in row
        assert row["h_bond_donors"] == 1

    def test_multiple_molecules(self):
        result = _json(compute_descriptors(
            ["CCO", "c1ccccc1", "CC(=O)O"],
            ["molecular_weight", "tpsa"],
        ))
        assert len(result) == 3
        assert all("molecular_weight" in r for r in result)
        assert all("tpsa" in r for r in result)

    def test_unknown_descriptor(self):
        with pytest.raises(ValueError, match="Unknown descriptors"):
            compute_descriptors(["CCO"], ["nonexistent_descriptor"])

    def test_max_50_molecules(self):
        smiles = ["CCO"] * 55
        result = _json(compute_descriptors(smiles, ["molecular_weight"]))
        assert len(result) == 50  # capped

    def test_all_descriptor_types(self):
        result = _json(compute_descriptors(
            ["CCO"],
            ["molecular_weight", "exact_molecular_weight", "logp", "tpsa",
             "h_bond_donors", "h_bond_acceptors", "rotatable_bonds",
             "aromatic_rings", "total_rings", "heavy_atoms", "num_atoms",
             "fraction_csp3", "heteroatoms", "kappa1", "kappa2", "kappa3"],
        ))
        assert len(result) == 1
        row = result[0]
        assert all(v != f"error:" for k, v in row.items() if k != "smiles")


# ════════════════════════════════════════════════════════════
#  3. Format Conversion
# ════════════════════════════════════════════════════════════

class TestConvertFormat:
    def test_smiles_to_inchi(self):
        result = _json(convert_format("CCO"))
        assert result["canonical_smiles"] == "CCO"
        assert result["inchi"].startswith("InChI=1S/")
        assert len(result["inchikey"]) == 27  # InChIKey format: XXXXXXXXXXXXXX-XXXXXXXXXX-X

    def test_canonicalization(self):
        # Different SMILES for the same molecule should give same canonical form
        r1 = _json(convert_format("C(C)O"))
        r2 = _json(convert_format("OCC"))
        r3 = _json(convert_format("CCO"))
        assert r1["canonical_smiles"] == r2["canonical_smiles"] == r3["canonical_smiles"]

    def test_inchikey_format(self):
        result = _json(convert_format("c1ccccc1"))
        inchikey = result["inchikey"]
        assert "-" in inchikey
        parts = inchikey.split("-")
        assert len(parts[0]) == 14


class TestInchiToSmiles:
    def test_roundtrip(self):
        # SMILES -> InChI -> SMILES should give the same molecule
        r1 = _json(convert_format("CCO"))
        inchi = r1["inchi"]
        r2 = _json(inchi_to_smiles(inchi))
        assert r2["smiles"] == "CCO"

    def test_invalid_inchi(self):
        with pytest.raises(ValueError, match="Invalid InChI"):
            inchi_to_smiles("not_a_valid_inchi")


# ════════════════════════════════════════════════════════════
#  4. Equation Balancing
# ════════════════════════════════════════════════════════════

class TestBalanceEquation:
    def test_simple_water(self):
        result = _json(balance_equation("H2 + O2 -> H2O"))
        assert "balanced" in result
        assert result["reactant_coefficients"]["H2"] == 2
        assert result["reactant_coefficients"]["O2"] == 1
        assert result["product_coefficients"]["H2O"] == 2

    def test_combustion(self):
        result = _json(balance_equation("CH4 + O2 -> CO2 + H2O"))
        assert result["reactant_coefficients"]["CH4"] == 1
        assert result["reactant_coefficients"]["O2"] == 2
        assert result["product_coefficients"]["CO2"] == 1
        assert result["product_coefficients"]["H2O"] == 2

    def test_iron_oxide(self):
        result = _json(balance_equation("Fe + O2 -> Fe2O3"))
        assert result["reactant_coefficients"]["Fe"] == 4
        assert result["reactant_coefficients"]["O2"] == 3
        assert result["product_coefficients"]["Fe2O3"] == 2

    def test_equals_separator(self):
        result = _json(balance_equation("H2 + O2 = H2O"))
        assert "balanced" in result

    def test_no_separator(self):
        result = _json(balance_equation("H2 + O2 H2O"))
        assert "error" in result


# ════════════════════════════════════════════════════════════
#  5. Structure Analysis
# ════════════════════════════════════════════════════════════

class TestCheckSubstructure:
    def test_hydroxyl_in_ethanol(self):
        result = _json(check_substructure("CCO", "[OX2H]"))
        assert result["has_substructure"] is True
        assert result["match_count"] >= 1

    def test_no_benzene_in_ethanol(self):
        result = _json(check_substructure("CCO", "c1ccccc1"))
        assert result["has_substructure"] is False

    def test_benzene_in_aspirin(self):
        result = _json(check_substructure("CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1"))
        assert result["has_substructure"] is True

    def test_carboxylic_acid(self):
        result = _json(check_substructure("CC(=O)O", "[CX3](=O)[OX2H1]"))
        assert result["has_substructure"] is True

    def test_invalid_smarts(self):
        with pytest.raises(ValueError, match="Invalid SMARTS"):
            check_substructure("CCO", "not_valid_smarts%%%")

    def test_multiple_matches(self):
        # Two hydroxyl groups in ethylene glycol
        result = _json(check_substructure("OCCO", "[OX2H]"))
        assert result["match_count"] == 2


class TestCalculateSimilarity:
    def test_identical_molecules(self):
        result = _json(calculate_similarity("CCO", "CCO"))
        assert result["tanimoto_similarity"] == 1.0

    def test_different_molecules(self):
        result = _json(calculate_similarity("CCO", "c1ccccc1"))
        assert 0.0 <= result["tanimoto_similarity"] < 1.0

    def test_custom_params(self):
        result = _json(calculate_similarity("CCO", "CCO", radius=3, n_bits=1024))
        assert result["tanimoto_similarity"] == 1.0

    def test_similarity_range(self):
        result = _json(calculate_similarity("CCO", "CCCCCC"))
        sim = result["tanimoto_similarity"]
        assert 0.0 <= sim <= 1.0


class TestGetScaffold:
    def test_aspirin_scaffold(self):
        result = _json(get_scaffold("CC(=O)Oc1ccccc1C(=O)O"))
        assert "scaffold" in result
        assert result["scaffold"] != ""

    def test_generic_scaffold(self):
        result = _json(get_scaffold("CC(=O)Oc1ccccc1C(=O)O", generic=True))
        assert "scaffold" in result
        assert "generic_scaffold" in result
        # Generic scaffold should only contain carbon atoms
        gen = result["generic_scaffold"]
        assert "C" in gen

    def test_no_scaffold(self):
        # Simple acyclic molecule — scaffold should be empty or minimal
        result = _json(get_scaffold("CCO"))
        assert "scaffold" in result


class TestFragmentMolecule:
    def test_simple_fragmentation(self):
        result = _json(fragment_molecule("c1ccc(CC(=O)O)cc1", max_cuts=1))
        assert "fragments" in result
        assert result["fragment_count"] >= 0

    def test_max_cuts_param(self):
        r1 = _json(fragment_molecule("c1ccc(CC(=O)O)cc1", max_cuts=1))
        r3 = _json(fragment_molecule("c1ccc(CC(=O)O)cc1", max_cuts=3))
        # More cuts should generally give more fragments
        assert r3["fragment_count"] >= r1["fragment_count"]


# ════════════════════════════════════════════════════════════
#  6. Visualization & 3D
# ════════════════════════════════════════════════════════════

class TestDrawMolecule:
    def test_png(self):
        result = _json(draw_molecule("CCO"))
        assert result["format"] == "png"
        assert len(result["data"]) > 0
        # Valid base64 should decode without error
        import base64
        decoded = base64.b64decode(result["data"])
        assert decoded[:4] == b"\x89PNG"

    def test_svg(self):
        result = _json(draw_molecule("CCO", image_format="svg"))
        assert result["format"] == "svg"
        import base64
        decoded = base64.b64decode(result["data"]).decode()
        assert "<svg" in decoded

    def test_custom_size(self):
        result = _json(draw_molecule("c1ccccc1", width=600, height=300))
        assert result["format"] == "png"
        assert len(result["data"]) > 0


class TestDrawMoleculeGrid:
    def test_grid_png(self):
        result = _json(draw_molecule_grid(["CCO", "c1ccccc1", "CC(=O)O"]))
        assert result["format"] == "png"
        assert len(result["data"]) > 0

    def test_grid_with_legends(self):
        result = _json(draw_molecule_grid(
            ["CCO", "c1ccccc1"],
            legends=["Ethanol", "Benzene"],
        ))
        assert result["format"] == "png"
        assert len(result["data"]) > 0

    def test_grid_svg(self):
        result = _json(draw_molecule_grid(["CCO"], image_format="svg"))
        assert result["format"] == "svg"

    def test_grid_max_20(self):
        smiles = ["CCO"] * 25
        result = _json(draw_molecule_grid(smiles))
        assert result["format"] == "png"

    def test_grid_invalid_smiles_skipped(self):
        result = _json(draw_molecule_grid(["CCO", "INVALID_XYZ", "c1ccccc1"]))
        assert result["format"] == "png"
        # Should still produce image with valid molecules

    def test_grid_all_invalid(self):
        with pytest.raises(ValueError, match="No valid molecules"):
            draw_molecule_grid(["INVALID1", "INVALID2"])


class TestGenerate3DStructure:
    def test_single_conformer(self):
        result = _json(generate_3d_structure("CCO"))
        assert result["conformer_count"] == 1
        assert len(result["conformers"]) == 1
        assert "sdf" in result["conformers"][0]
        # SDF should contain 3D coordinates
        assert "V2000" in result["conformers"][0]["sdf"]

    def test_multiple_conformers(self):
        result = _json(generate_3d_structure("CCCCCC", num_conformers=5))
        assert result["conformer_count"] == 5
        assert len(result["conformers"]) == 5

    def test_ethane(self):
        result = _json(generate_3d_structure("CC"))
        assert result["smiles"] == "CC"
        assert result["conformer_count"] >= 1

    def test_reproducible_with_seed(self):
        r1 = _json(generate_3d_structure("CCO", random_seed=42))
        r2 = _json(generate_3d_structure("CCO", random_seed=42))
        assert r1["conformers"][0]["sdf"] == r2["conformers"][0]["sdf"]
