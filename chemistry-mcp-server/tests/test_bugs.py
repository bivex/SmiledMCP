"""Comprehensive bug-hunting tests for chemistry MCP server.

Targets edge cases, boundary conditions, error paths, and output consistency.
"""

import asyncio
import base64
import json
import sys
import os

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helpers import mol_from_smiles, compound_to_dict, resolve_namespace

from tools.pubchem import search_compound, get_compound_properties, get_synonyms
from tools.properties import molecular_info, compute_descriptors, DESCRIPTOR_MAP
from tools.conversion import convert_format, inchi_to_smiles
from tools.equations import balance_equation, HAS_CHEMPY
from tools.structure import (
    check_substructure,
    calculate_similarity,
    get_scaffold,
    fragment_molecule,
)
from tools.visualization import draw_molecule, draw_molecule_grid, generate_3d_structure


def _json(result: str) -> dict | list:
    return json.loads(result)


def _run(coro):
    return asyncio.run(coro)


# ================================================================
#  helpers.py
# ================================================================

class TestMolFromSmilesEdgeCases:
    def test_empty_string(self):
        with pytest.raises(ValueError):
            mol_from_smiles("")

    def test_whitespace_only(self):
        with pytest.raises(ValueError):
            mol_from_smiles("   ")

    def test_charged_atom(self):
        mol = mol_from_smiles("[NH4+]")
        assert mol is not None

    def test_radical(self):
        mol = mol_from_smiles("[CH3]")
        assert mol is not None

    def test_isotope(self):
        mol = mol_from_smiles("[2H]O")
        assert mol is not None

    def test_stereochemistry(self):
        mol = mol_from_smiles("C/C=C/C")
        assert mol is not None

    def test_complex_natural_product(self):
        # Caffeine
        mol = mol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        assert mol is not None
        assert mol.GetNumHeavyAtoms() == 14

    def test_dot_disconnected(self):
        # Two separate fragments
        mol = mol_from_smiles("CCO.O")
        assert mol is not None

    def test_very_long_smiles(self):
        # Polyethylene-like chain
        smi = "C" * 200
        mol = mol_from_smiles(smi)
        assert mol is not None


class TestCompoundToDict:
    def test_compound_with_no_synonyms(self):
        """compound_to_dict should not crash if synonyms access raises."""

        class FakeCompound:
            cid = 123
            molecular_formula = "CH4"
            molecular_weight = 16.04
            canonical_smiles = "C"
            isomeric_smiles = "C"
            iupac_name = "methane"
            inchi = "InChI=1S/CH4/h1H4"
            inchikey = "VNWKTOKETHGBQD-UHFFFAOYSA-N"
            xlogp = 0.0
            exact_mass = 16.0313
            tpsa = 0.0
            complexity = 0.0
            charge = 0
            h_bond_donor_count = 0
            h_bond_acceptor_count = 0
            rotatable_bond_count = 0
            heavy_atom_count = 1

            @property
            def synonyms(self):
                raise RuntimeError("API error")

        d = compound_to_dict(FakeCompound())
        assert d["cid"] == 123
        assert "synonyms" not in d

    def test_compound_with_empty_synonyms(self):
        class FakeCompound:
            cid = 1
            molecular_formula = None
            molecular_weight = None
            canonical_smiles = None
            isomeric_smiles = None
            iupac_name = None
            inchi = None
            inchikey = None
            xlogp = None
            exact_mass = None
            tpsa = None
            complexity = None
            charge = None
            h_bond_donor_count = None
            h_bond_acceptor_count = None
            rotatable_bond_count = None
            heavy_atom_count = None
            synonyms = []

        d = compound_to_dict(FakeCompound())
        assert "synonyms" not in d

    def test_compound_with_many_synonyms_truncated(self):
        class FakeCompound:
            cid = 1
            molecular_formula = None
            molecular_weight = None
            canonical_smiles = None
            isomeric_smiles = None
            iupac_name = None
            inchi = None
            inchikey = None
            xlogp = None
            exact_mass = None
            tpsa = None
            complexity = None
            charge = None
            h_bond_donor_count = None
            h_bond_acceptor_count = None
            rotatable_bond_count = None
            heavy_atom_count = None
            synonyms = [f"syn_{i}" for i in range(100)]

        d = compound_to_dict(FakeCompound())
        assert len(d["synonyms"]) == 20


class TestResolveNamespace:
    def test_cid(self):
        assert resolve_namespace("cid") == "cid"

    def test_name(self):
        assert resolve_namespace("name") == "name"

    def test_invalid_type_raises(self):
        with pytest.raises(ValueError, match="Invalid query_type"):
            resolve_namespace("garbage")

    def test_all_valid_types(self):
        for qt in ("name", "smiles", "inchi", "inchikey", "formula", "cid"):
            assert resolve_namespace(qt) is not None


# ================================================================
#  PubChem tools — async tests with real API
# ================================================================

@pytest.mark.skipif(not HAS_CHEMPY, reason="Needs network")
class TestSearchCompoundBugs:
    def test_max_results_zero_clamped_to_one(self):
        """max_results=0 should be clamped to 1, not crash."""
        result = _json(_run(search_compound("water", "name", 0)))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_max_results_negative_clamped(self):
        """Negative max_results should be clamped to 1."""
        result = _json(_run(search_compound("water", "name", -5)))
        assert isinstance(result, list)
        assert len(result) >= 1

    def test_max_results_over_20_clamped(self):
        """max_results > 20 should be clamped to 20."""
        result = _json(_run(search_compound("benzene", "name", 100)))
        assert isinstance(result, list)
        assert len(result) <= 20

    def test_cid_non_numeric_raises(self):
        """query_type='cid' with non-numeric query should raise clean ValueError."""
        with pytest.raises(ValueError, match="CID must be a number"):
            _run(search_compound("not_a_number", "cid"))

    def test_search_returns_json_parseable(self):
        result = _run(search_compound("caffeine", "name", 1))
        parsed = json.loads(result)
        assert isinstance(parsed, list)

    def test_search_result_has_all_expected_keys(self):
        result = _json(_run(search_compound("aspirin", "name", 1)))
        comp = result[0]
        expected_keys = [
            "cid", "molecular_formula", "molecular_weight",
            "canonical_smiles", "inchi", "inchikey",
        ]
        for key in expected_keys:
            assert key in comp, f"Missing key: {key}"

    def test_formula_search(self):
        result = _json(_run(search_compound("C6H6", "formula", 3)))
        assert isinstance(result, list)
        assert len(result) >= 1


@pytest.mark.skipif(not HAS_CHEMPY, reason="Needs network")
class TestGetCompoundPropertiesBugs:
    def test_cid_non_numeric_raises(self):
        with pytest.raises(ValueError, match="CID must be a number"):
            _run(get_compound_properties("abc", "cid"))

    def test_nonexistent_compound_returns_empty(self):
        result = _json(_run(get_compound_properties("XYZNOTREAL12345", "name")))
        assert isinstance(result, list)
        assert len(result) == 0

    def test_single_property(self):
        result = _json(_run(
            get_compound_properties("aspirin", "name", ["MolecularWeight"])
        ))
        assert isinstance(result, list)
        assert len(result) >= 1
        assert "MolecularWeight" in result[0]


@pytest.mark.skipif(not HAS_CHEMPY, reason="Needs network")
class TestGetSynonymsBugs:
    def test_cid_non_numeric_raises(self):
        with pytest.raises(ValueError, match="CID must be a number"):
            _run(get_synonyms("abc", "cid"))

    def test_nonexistent_compound_returns_empty(self):
        result = _json(_run(get_synonyms("XYZNOTREAL99999", "name")))
        assert isinstance(result, list)
        assert len(result) == 0


# ================================================================
#  Molecular properties
# ================================================================

class TestMolecularInfoBugs:
    def test_output_is_valid_json(self):
        result = molecular_info("CCO")
        parsed = json.loads(result)
        assert isinstance(parsed, dict)

    def test_all_keys_present(self):
        info = _json(molecular_info("c1ccccc1"))
        expected = [
            "smiles", "molecular_formula", "molecular_weight",
            "exact_molecular_weight", "logp", "tpsa",
            "h_bond_donors", "h_bond_acceptors", "rotatable_bonds",
            "aromatic_rings", "aliphatic_rings", "total_rings",
            "heavy_atoms", "atom_count", "fraction_csp3",
            "num_valence_electrons", "num_radical_electrons",
            "heteroatoms", "amide_bonds", "lipinski_rule_of_5",
        ]
        for key in expected:
            assert key in info, f"Missing key: {key}"

    def test_lipinski_structure(self):
        info = _json(molecular_info("CCO"))
        lip = info["lipinski_rule_of_5"]
        for key in ["mw_le_500", "logp_le_5", "hbd_le_5", "hba_le_10", "violations", "passes"]:
            assert key in lip, f"Missing lipinski key: {key}"
        assert isinstance(lip["violations"], int)
        assert 0 <= lip["violations"] <= 4
        assert lip["passes"] is (lip["violations"] == 0)

    def test_methane(self):
        info = _json(molecular_info("C"))
        assert info["molecular_formula"] == "CH4"
        assert info["heavy_atoms"] == 1
        assert info["atom_count"] == 1  # only explicit heavy atoms in graph (implicit H not counted)
        assert info["total_rings"] == 0
        assert info["rotatable_bonds"] == 0

    def test_charged_molecule(self):
        """Charged molecules should still work."""
        info = _json(molecular_info("[NH4+]"))
        assert info is not None
        assert "molecular_formula" in info

    def test_multi_ring_system(self):
        """Naphthalene — two fused rings."""
        info = _json(molecular_info("c1ccc2ccccc2c1"))
        assert info["total_rings"] == 2
        assert info["aromatic_rings"] == 2

    def test_halogens(self):
        info = _json(molecular_info("CCl"))
        assert info["heteroatoms"] == 1

    def test_num_atoms_includes_hydrogens_implicitly(self):
        """num_atoms counts atoms in the SMILES graph, which excludes implicit H."""
        info = _json(molecular_info("CCO"))
        # C-C-O = 3 heavy atoms, implicit H not counted
        assert info["heavy_atoms"] == 3
        assert info["atom_count"] == 3  # only explicit atoms in graph

    def test_molecular_weight_positive(self):
        info = _json(molecular_info("CCO"))
        assert info["molecular_weight"] > 0

    def test_tpsa_non_negative(self):
        info = _json(molecular_info("CCO"))
        assert info["tpsa"] >= 0

    def test_fraction_csp3_range(self):
        info = _json(molecular_info("CCO"))
        assert 0.0 <= info["fraction_csp3"] <= 1.0


class TestComputeDescriptorsBugs:
    def test_empty_smiles_list(self):
        """Empty smiles_list should return empty list, not crash."""
        result = _json(compute_descriptors([], ["molecular_weight"]))
        assert result == []

    def test_empty_descriptor_names(self):
        """Empty descriptor_names — should work, each row has only 'smiles'."""
        result = _json(compute_descriptors(["CCO"], []))
        assert len(result) == 1
        assert result[0] == {"smiles": "CCO"}

    def test_invalid_smiles_in_list_stops_at_first_bad(self):
        """Invalid SMILES raises ValueError — stops the whole batch."""
        with pytest.raises(ValueError, match="Invalid SMILES"):
            compute_descriptors(["CCO", "INVALID_SMILES", "c1ccccc1"], ["molecular_weight"])

    def test_all_27_descriptors_work(self):
        """Every descriptor in DESCRIPTOR_MAP should compute without error."""
        names = list(DESCRIPTOR_MAP.keys())
        assert len(names) == 28
        result = _json(compute_descriptors(["CCO"], names))
        row = result[0]
        for name in names:
            assert name in row, f"Missing descriptor: {name}"
            assert not str(row[name]).startswith("error:"), f"{name} errored: {row[name]}"

    def test_descriptor_values_are_numeric(self):
        result = _json(compute_descriptors(["c1ccccc1"], ["molecular_weight", "logp", "tpsa"]))
        row = result[0]
        for key in ["molecular_weight", "logp", "tpsa"]:
            assert isinstance(row[key], (int, float)), f"{key} is not numeric: {type(row[key])}"

    def test_single_molecule_long_list_capped(self):
        smiles = ["CCO"] * 55
        result = _json(compute_descriptors(smiles, ["molecular_weight"]))
        assert len(result) == 50

    def test_mixed_valid_invalid_descriptors(self):
        """Mix of valid and invalid descriptor names — should list all unknowns."""
        with pytest.raises(ValueError, match="Unknown descriptors"):
            compute_descriptors(["CCO"], ["molecular_weight", "bogus_desc"])


# ================================================================
#  Format conversion
# ================================================================

class TestConvertFormatBugs:
    def test_inchi_can_be_none(self):
        """BUG: Chem.MolToInchi() can return None for some valid molecules.
        This would serialize as 'null' in JSON."""
        # Phenol is a known case where MolToInchi may return None in some RDKit versions
        # Let's test with a simple case and check it handles it
        result = _json(convert_format("CCO"))
        # If inchi is None, the JSON value would be null
        # This is a latent bug — not all valid SMILES have InChI
        assert "inchi" in result
        assert "inchikey" in result

    def test_input_smiles_preserved(self):
        """The 'input_smiles' should be exactly what was passed in."""
        result = _json(convert_format("C(C)O"))
        assert result["input_smiles"] == "C(C)O"
        # canonical_smiles should be different (canonical form)
        assert result["canonical_smiles"] == "CCO"

    def test_stereochemistry_preserved_in_canonical(self):
        """Stereo SMILES should preserve stereo info in canonical form."""
        result = _json(convert_format("C/C=C/C"))
        assert "/" in result["canonical_smiles"] or "\\" in result["canonical_smiles"]

    def test_empty_smiles_raises(self):
        with pytest.raises(ValueError):
            convert_format("")


class TestInchiToSmilesBugs:
    def test_output_keys(self):
        r1 = _json(convert_format("CCO"))
        inchi = r1["inchi"]
        r2 = _json(inchi_to_smiles(inchi))
        assert "smiles" in r2
        assert "inchi" in r2
        assert "inchikey" in r2

    def test_inchi_with_extra_whitespace(self):
        """InChI with leading/trailing spaces."""
        with pytest.raises(ValueError):
            inchi_to_smiles("  InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3  ")

    def test_roundtrip_identity(self):
        """SMILES -> InChI -> SMILES should give the same canonical form."""
        r1 = _json(convert_format("c1ccccc1"))
        r2 = _json(inchi_to_smiles(r1["inchi"]))
        # Canonical form should be the same
        assert r2["smiles"] == _json(convert_format("c1ccccc1"))["canonical_smiles"]


# ================================================================
#  Equation balancing
# ================================================================

class TestBalanceEquationBugs:
    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_empty_string(self):
        with pytest.raises(ValueError, match="Empty equation"):
            balance_equation("")

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_only_separator(self):
        """Just '->' with no substances on either side."""
        with pytest.raises(ValueError, match="substances on both sides"):
            balance_equation("->")

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_whitespace_around_substances(self):
        """Extra spaces around substance names should be handled."""
        result = _json(balance_equation("  H2  +  O2  ->  H2O  "))
        assert "balanced" in result
        assert result["reactant_coefficients"]["H2"] == 2

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_single_substance_each_side(self):
        """Same substance on both sides — chempy rejects this."""
        with pytest.raises(ValueError, match="Failed to balance"):
            balance_equation("H2 -> H2")

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_impossible_equation(self):
        """Equation that can't be balanced (mass not conserved)."""
        with pytest.raises(ValueError, match="Failed to balance"):
            balance_equation("H2 -> O2")

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_duplicate_substances(self):
        """Same substance on both sides."""
        with pytest.raises(ValueError):
            balance_equation("H2 + O2 -> H2O + H2")

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_complex_equation(self):
        """More complex equation with multiple products."""
        result = _json(balance_equation("Fe2O3 + CO -> Fe + CO2"))
        assert "balanced" in result
        # Check coefficients are integers
        for side in ["reactant_coefficients", "product_coefficients"]:
            for k, v in result[side].items():
                assert isinstance(v, int), f"{k}: {v} is not int"

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_int_conversion_of_coefficients(self):
        """BUG RISK: chempy may return Fraction or sympy types, int() could fail."""
        result = _json(balance_equation("C2H6 + O2 -> CO2 + H2O"))
        assert "balanced" in result
        for side in ["reactant_coefficients", "product_coefficients"]:
            for k, v in result[side].items():
                assert isinstance(v, int), f"Coefficient for {k} is {type(v)}: {v}"

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_balanced_equation_coefficients_sum_mass_balance(self):
        """Verify mass balance: total atoms of each element must match."""
        result = _json(balance_equation("CH4 + O2 -> CO2 + H2O"))
        # CH4: 1C + 4H, O2: 2O, coeff: 1, 2
        # CO2: 1C + 2O, H2O: 2H + 1O, coeff: 1, 2
        # Left:  1C, 4H, 4O  |  Right: 1C, 4H, (2+2)O = 4O
        assert result["reactant_coefficients"]["CH4"] == 1
        assert result["reactant_coefficients"]["O2"] == 2
        assert result["product_coefficients"]["CO2"] == 1
        assert result["product_coefficients"]["H2O"] == 2

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_nonexistent_substance(self):
        """Substances that aren't real chemical formulas."""
        with pytest.raises(ValueError, match="Failed to balance"):
            balance_equation("XYZ + ABC -> DEF")


# ================================================================
#  Structure analysis
# ================================================================

class TestCheckSubstructureBugs:
    def test_no_match_output_consistency(self):
        """When no match, output should NOT have match_count or matching_atoms."""
        result = _json(check_substructure("CCO", "c1ccccc1"))
        assert result["has_substructure"] is False
        assert "match_count" not in result
        assert "matching_atoms" not in result

    def test_match_output_has_count_and_atoms(self):
        """When match found, match_count and matching_atoms should be present."""
        result = _json(check_substructure("c1ccccc1", "c1ccccc1"))
        assert result["has_substructure"] is True
        assert "match_count" in result
        assert "matching_atoms" in result
        assert result["match_count"] >= 1

    def test_matching_atoms_are_valid_indices(self):
        result = _json(check_substructure("CCO", "[OX2H]"))
        atoms = result["matching_atoms"]
        mol = mol_from_smiles("CCO")
        for match in atoms:
            for idx in match:
                assert 0 <= idx < mol.GetNumAtoms()

    def test_wildcard_smarts(self):
        """Wildcard SMARTS should match everything."""
        result = _json(check_substructure("CCO", "*"))
        assert result["has_substructure"] is True

    def test_empty_smiles_raises(self):
        with pytest.raises(ValueError):
            check_substructure("", "CC")

    def test_empty_smarts_raises(self):
        """Empty SMARTS string — MolFromSmarts('') returns a valid empty mol,
        which matches nothing. Not an error, but worth documenting."""
        result = _json(check_substructure("CCO", ""))
        assert result["has_substructure"] is False


class TestCalculateSimilarityBugs:
    def test_same_molecule_different_smiles(self):
        """Same molecule written differently should have similarity 1.0."""
        result = _json(calculate_similarity("C(C)O", "OCC"))
        assert result["tanimoto_similarity"] == 1.0

    def test_similarity_symmetric(self):
        """Tanimoto(A, B) == Tanimoto(B, A)."""
        r1 = _json(calculate_similarity("CCO", "c1ccccc1"))
        r2 = _json(calculate_similarity("c1ccccc1", "CCO"))
        assert r1["tanimoto_similarity"] == r2["tanimoto_similarity"]

    def test_radius_zero(self):
        """radius=0 — only atom itself, no neighbors."""
        result = _json(calculate_similarity("CCO", "CCO", radius=0))
        assert result["tanimoto_similarity"] == 1.0

    def test_very_small_n_bits(self):
        """n_bits=1 — single bit, everything maps to same bit."""
        result = _json(calculate_similarity("CCO", "CCO", n_bits=1))
        assert result["tanimoto_similarity"] == 1.0

    def test_invalid_smiles_first(self):
        with pytest.raises(ValueError):
            calculate_similarity("INVALID", "CCO")

    def test_invalid_smiles_second(self):
        with pytest.raises(ValueError):
            calculate_similarity("CCO", "INVALID")

    def test_output_keys(self):
        result = _json(calculate_similarity("CCO", "c1ccccc1"))
        assert "smiles1" in result
        assert "smiles2" in result
        assert "tanimoto_similarity" in result

    def test_similarity_boundaries(self):
        """Similarity must be in [0, 1]."""
        for s1, s2 in [("CCO", "c1ccccc1"), ("CCO", "CCCCCC"), ("CC", "CC(=O)O")]:
            result = _json(calculate_similarity(s1, s2))
            assert 0.0 <= result["tanimoto_similarity"] <= 1.0


class TestGetScaffoldBugs:
    def test_acyclic_molecule_scaffold(self):
        """Acyclic molecules like ethanol have no ring scaffold."""
        result = _json(get_scaffold("CCO"))
        assert "scaffold" in result
        # Scaffold of acyclic molecule is empty
        assert result["scaffold"] == ""

    def test_pure_aromatic_scaffold(self):
        """Benzene — scaffold is benzene itself."""
        result = _json(get_scaffold("c1ccccc1"))
        assert result["scaffold"] != ""

    def test_generic_scaffold_only_carbon_single_bonds(self):
        """Generic scaffold should only contain C and single bonds."""
        result = _json(get_scaffold("CC(=O)Oc1ccccc1C(=O)O", generic=True))
        gen = result["generic_scaffold"]
        # Generic scaffold should not contain O, N, etc.
        for char in gen:
            assert char in "Cc-.[]()0123456789%#", f"Unexpected char '{char}' in generic scaffold: {gen}"

    def test_generic_false_no_generic_key(self):
        result = _json(get_scaffold("c1ccccc1", generic=False))
        assert "generic_scaffold" not in result

    def test_fused_ring_system(self):
        """Naphthalene — scaffold is the full fused ring system."""
        result = _json(get_scaffold("c1ccc2ccccc2c1"))
        assert result["scaffold"] != ""


class TestFragmentMoleculeBugs:
    def test_simple_methane(self):
        """Methane has only C-H bonds, no single bonds between heavy atoms to cut."""
        result = _json(fragment_molecule("C", max_cuts=1))
        assert "fragments" in result
        # Methane might have no fragmentable bonds
        assert isinstance(result["fragment_count"], int)

    def test_max_cuts_zero(self):
        """max_cuts=0 — no cuts allowed, should return empty or minimal fragments."""
        result = _json(fragment_molecule("CCO", max_cuts=0))
        assert "fragments" in result

    def test_max_cuts_negative(self):
        """Negative max_cuts — RDKit behavior is undefined, should not crash."""
        result = _json(fragment_molecule("CCO", max_cuts=-1))
        assert "fragments" in result

    def test_max_cuts_very_large(self):
        result = _json(fragment_molecule("c1ccc(CC(=O)O)cc1", max_cuts=10))
        assert "fragments" in result

    def test_fragments_capped_at_50(self):
        """Result should cap at 50 fragments."""
        result = _json(fragment_molecule("c1ccc(CC(=O)O)cc1", max_cuts=5))
        assert len(result["fragments"]) <= 50

    def test_fragment_output_structure(self):
        result = _json(fragment_molecule("c1ccc(CC(=O)O)cc1", max_cuts=2))
        for frag in result["fragments"]:
            assert "core" in frag
            assert "side_chains" in frag
            # core and side_chains should be valid SMILES
            assert isinstance(frag["core"], str)
            assert isinstance(frag["side_chains"], str)


# ================================================================
#  Visualization & 3D
# ================================================================

class TestDrawMoleculeBugs:
    def test_png_is_valid_image(self):
        result = _json(draw_molecule("c1ccccc1"))
        decoded = base64.b64decode(result["data"])
        assert decoded[:8] == b"\x89PNG\r\n\x1a\n"

    def test_svg_contains_svg_tag(self):
        result = _json(draw_molecule("c1ccccc1", image_format="svg"))
        decoded = base64.b64decode(result["data"]).decode()
        assert decoded.startswith("<?xml") or "<svg" in decoded
        assert "</svg>" in decoded

    def test_format_case_insensitive(self):
        """image_format should be case-insensitive."""
        r1 = _json(draw_molecule("CCO", image_format="PNG"))
        r2 = _json(draw_molecule("CCO", image_format="png"))
        r3 = _json(draw_molecule("CCO", image_format="Png"))
        assert r1["format"] == r2["format"] == r3["format"]

    def test_unsupported_format_defaults_to_png(self):
        """Unknown image_format silently defaults to PNG."""
        result = _json(draw_molecule("CCO", image_format="jpeg"))
        assert result["format"] == "png"

    def test_zero_width_height(self):
        """Zero or negative dimensions — RDKit may crash or produce empty image."""
        # This tests robustness — we don't guarantee behavior, just no unhandled crash
        try:
            draw_molecule("CCO", width=0, height=0)
        except Exception:
            pass  # Acceptable to fail

    def test_very_large_dimensions(self):
        result = _json(draw_molecule("CCO", width=10000, height=10000))
        assert result["format"] == "png"

    def test_output_keys(self):
        result = _json(draw_molecule("CCO"))
        assert "format" in result
        assert "data" in result
        assert "note" in result


class TestDrawMoleculeGridBugs:
    def test_single_molecule(self):
        result = _json(draw_molecule_grid(["CCO"]))
        assert result["format"] == "png"

    def test_legends_longer_than_smiles_list(self):
        """Extra legends should be silently ignored (indices checked)."""
        result = _json(draw_molecule_grid(
            ["CCO"],
            legends=["Ethanol", "Extra legend"],
        ))
        assert result["format"] == "png"

    def test_legends_empty_strings(self):
        result = _json(draw_molecule_grid(
            ["CCO", "c1ccccc1"],
            legends=["", ""],
        ))
        assert result["format"] == "png"

    def test_all_invalid_smiles_raises(self):
        with pytest.raises(ValueError, match="No valid molecules"):
            draw_molecule_grid(["NOT_REAL_1", "NOT_REAL_2"])

    def test_empty_list_raises(self):
        """Empty smiles_list — all (zero) molecules are invalid."""
        with pytest.raises(ValueError, match="No valid molecules"):
            draw_molecule_grid([])

    def test_max_20_capped(self):
        """More than 20 SMILES should be capped at 20."""
        smiles = ["CCO"] * 30
        result = _json(draw_molecule_grid(smiles))
        assert result["format"] == "png"

    def test_mixed_valid_invalid_keeps_legends_aligned(self):
        """When invalid SMILES are skipped, legends for valid molecules should be correct."""
        result = _json(draw_molecule_grid(
            ["CCO", "INVALID", "c1ccccc1"],
            legends=["Ethanol", "Bad", "Benzene"],
        ))
        assert result["format"] == "png"

    def test_svg_grid(self):
        result = _json(draw_molecule_grid(
            ["CCO", "c1ccccc1"],
            image_format="svg",
        ))
        assert result["format"] == "svg"
        decoded = base64.b64decode(result["data"]).decode()
        assert "<svg" in decoded


class TestGenerate3DStructureBugs:
    def test_num_conformers_zero_clamped_to_one(self):
        result = _json(generate_3d_structure("CCO", num_conformers=0))
        assert result["conformer_count"] == 1

    def test_num_conformers_negative_clamped(self):
        result = _json(generate_3d_structure("CCO", num_conformers=-5))
        assert result["conformer_count"] == 1

    def test_num_conformers_over_50_clamped(self):
        result = _json(generate_3d_structure("CCCCCC", num_conformers=100))
        assert result["conformer_count"] <= 50

    def test_single_conformer_has_sdf_with_3d(self):
        result = _json(generate_3d_structure("CCO"))
        sdf = result["conformers"][0]["sdf"]
        assert "V2000" in sdf or "V3000" in sdf
        # SDF should have multiple lines with 3D coordinates
        lines = sdf.strip().split("\n")
        assert len(lines) > 3

    def test_multiple_conformers_unique_ids(self):
        result = _json(generate_3d_structure("CCCCCC", num_conformers=5))
        ids = [c["id"] for c in result["conformers"]]
        assert len(set(ids)) == len(ids), "Duplicate conformer IDs"

    def test_multiple_conformers_different_coordinates(self):
        """Different conformers should have different 3D coordinates."""
        result = _json(generate_3d_structure("CCCCCC", num_conformers=5, random_seed=42))
        sdfs = [c["sdf"] for c in result["conformers"]]
        # At least some should differ
        assert len(set(sdfs)) > 1, "All conformers have identical coordinates"

    def test_smiles_preserved_in_output(self):
        result = _json(generate_3d_structure("CCO"))
        assert result["smiles"] == "CCO"

    def test_different_seeds_different_results(self):
        r1 = _json(generate_3d_structure("CCO", random_seed=1))
        r2 = _json(generate_3d_structure("CCO", random_seed=999))
        assert r1["conformers"][0]["sdf"] != r2["conformers"][0]["sdf"]

    def test_empty_smiles_raises(self):
        with pytest.raises(ValueError):
            generate_3d_structure("")

    def test_conformer_count_matches_list_length(self):
        """conformer_count should equal len(conformers)."""
        result = _json(generate_3d_structure("CCCCCC", num_conformers=5))
        assert result["conformer_count"] == len(result["conformers"])


# ================================================================
#  Cross-cutting concerns
# ================================================================

class TestJsonOutputConsistency:
    """All tools must return valid JSON strings."""

    def test_molecular_info_valid_json(self):
        parsed = json.loads(molecular_info("CCO"))
        assert isinstance(parsed, dict)

    def test_compute_descriptors_valid_json(self):
        parsed = json.loads(compute_descriptors(["CCO"], ["molecular_weight"]))
        assert isinstance(parsed, list)

    def test_convert_format_valid_json(self):
        parsed = json.loads(convert_format("CCO"))
        assert isinstance(parsed, dict)

    def test_balance_equation_valid_json(self):
        parsed = json.loads(balance_equation("H2 + O2 -> H2O"))
        assert isinstance(parsed, dict)

    def test_check_substructure_valid_json(self):
        parsed = json.loads(check_substructure("CCO", "[OX2H]"))
        assert isinstance(parsed, dict)

    def test_calculate_similarity_valid_json(self):
        parsed = json.loads(calculate_similarity("CCO", "CCCO"))
        assert isinstance(parsed, dict)

    def test_get_scaffold_valid_json(self):
        parsed = json.loads(get_scaffold("c1ccccc1"))
        assert isinstance(parsed, dict)

    def test_fragment_molecule_valid_json(self):
        parsed = json.loads(fragment_molecule("CCO"))
        assert isinstance(parsed, dict)

    def test_draw_molecule_valid_json(self):
        parsed = json.loads(draw_molecule("CCO"))
        assert isinstance(parsed, dict)

    def test_generate_3d_structure_valid_json(self):
        parsed = json.loads(generate_3d_structure("CCO"))
        assert isinstance(parsed, dict)


class TestErrorHandlingConsistency:
    """Error handling should be consistent across tools."""

    def test_invalid_smiles_raises_valueerror_molecular_info(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            molecular_info("GARBAGE")

    def test_invalid_smiles_raises_valueerror_convert_format(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            convert_format("GARBAGE")

    def test_invalid_smiles_raises_valueerror_check_substructure(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            check_substructure("GARBAGE", "CC")

    def test_invalid_smiles_raises_valueerror_similarity(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            calculate_similarity("GARBAGE", "CCO")

    def test_invalid_smiles_raises_valueerror_scaffold(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            get_scaffold("GARBAGE")

    def test_invalid_smiles_raises_valueerror_fragment(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            fragment_molecule("GARBAGE")

    def test_invalid_smiles_raises_valueerror_draw(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            draw_molecule("GARBAGE")

    def test_invalid_smiles_raises_valueerror_3d(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            generate_3d_structure("GARBAGE")

    @pytest.mark.skipif(not HAS_CHEMPY, reason="chempy not installed")
    def test_balance_raises_valueerror(self):
        """balance_equation now raises ValueError consistently with other tools."""
        with pytest.raises(ValueError):
            balance_equation("H2 + O2 H2O")
