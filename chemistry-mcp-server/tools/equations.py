"""Equation balancing tools — chempy stoichiometry."""

import re
import json

from server import mcp

try:
    from chempy import balance_stoichiometry
    HAS_CHEMPY = True
except ImportError:
    HAS_CHEMPY = False


def _parse_int_coeff(v) -> int:
    """Safely convert chempy coefficient to int (handles Fraction, sympy, etc.)."""
    return int(v)


# Patterns that suggest a biological polymer rather than a small molecule
_POLYMER_INDICATORS = [
    re.compile(r'\((?:[A-Z][a-z]?[0-9]*[A-Z][a-z]?[0-9]*){30,}\)'),  # long AA/NA sequence
    re.compile(r'(?:Ala|Arg|Asn|Asp|Cys|Glu|Gln|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr){3,}'),
    re.compile(r'(?:A|C|G|T|U){10,}'),  # long nucleotide runs (not chemical formulas)
    re.compile(r'~[A-Z][a-z]?[0-9]*~'),  # repeating unit notation
    re.compile(r'\(\d+,\d+\)'),          # large coordinate pairs
]


def _looks_like_polymer(formula: str) -> bool:
    """Detect if a substance string looks like a biological polymer."""
    for pat in _POLYMER_INDICATORS:
        if pat.search(formula):
            return True
    # Check for very long chemical formulas (polymers have huge atom counts)
    # A reasonable small-molecule formula has < 50 heavy atom symbols
    stripped = re.sub(r'[+\-()]', '', formula)
    try:
        total_atoms = sum(
            int(m.group(1) or 1)
            for m in re.finditer(r'[A-Z][a-z]?(\d*)', stripped)
        )
    except (ValueError, AttributeError):
        total_atoms = 0
    if total_atoms > 200:
        return True
    # Check for peptide/RNA-style repeating patterns
    if re.search(r'(?:NH?-?CH?R?-?CO?)', formula) or re.search(r'(?:ribose|deoxyribose|phosphate)', formula, re.IGNORECASE):
        return True
    return False


@mcp.tool()
def balance_equation(equation: str) -> str:
    """Balance a chemical equation.

    Uses chempy for stoichiometric balancing. Works well for small molecules
    and inorganic reactions.

    Examples:
      "H2 + O2 -> H2O"
      "Fe + O2 -> Fe2O3"
      "CH4 + O2 -> CO2 + H2O"

    Note:
      This tool is designed for **small-molecule** stoichiometry. It does not
      handle biological polymers (proteins, RNA, DNA, polysaccharides) because
      chempy uses integer linear programming on atom counts, and polymer
      repeat units produce degenerate or trivially large solutions. For
      polymer balancing, use a dedicated polymer chemistry tool or balance
      the *monomer* reaction instead.

    Args:
        equation: Chemical equation with '->' or '=' separator

    Returns:
        JSON with balanced equation and stoichiometric coefficients

    Raises:
        ValueError: If the equation format is invalid or balancing fails
    """
    if not HAS_CHEMPY:
        raise ValueError("chempy is not installed. Run: pip install chempy")

    if not equation or not equation.strip():
        raise ValueError("Empty equation")

    if "->" in equation:
        left_str, right_str = equation.split("->", 1)
    elif "=" in equation:
        left_str, right_str = equation.split("=", 1)
    else:
        raise ValueError("Use '->' or '=' to separate reactants and products")

    def parse_substances(s):
        return [part.strip() for part in s.split("+") if part.strip()]

    reactants = parse_substances(left_str)
    products = parse_substances(right_str)

    if not reactants or not products:
        raise ValueError("Equation must have substances on both sides")

    # Warn if any substance looks like a polymer
    for r in reactants:
        if _looks_like_polymer(r):
            return json.dumps({
                "error": "polymer_detected",
                "message": (
                    f"Reactant '{r}' appears to be a biological polymer. "
                    f"Equation balancing is designed for small molecules (typically < 50 atoms). "
                    f"For polymer reactions, balance the monomer unit instead, e.g.:\n"
                    f"  Instead of: AminoAcid_n + ... -> Protein + ...\n"
                    f"  Try:        AminoAcid + ... -> Dipeptide + H2O\n"
                    f"  Then note: n AminoAcids → (AminoAcid)_n + (n-1) H2O"
                ),
                "detected_polymer": r,
            }, indent=2, ensure_ascii=False)

    for p in products:
        if _looks_like_polymer(p):
            return json.dumps({
                "error": "polymer_detected",
                "message": (
                    f"Product '{p}' appears to be a biological polymer. "
                    f"Equation balancing is designed for small molecules (typically < 50 atoms). "
                    f"For polymer reactions, balance the monomer unit instead."
                ),
                "detected_polymer": p,
            }, indent=2, ensure_ascii=False)

    try:
        balanced_r, balanced_p = balance_stoichiometry(reactants, products)
    except Exception as e:
        # Provide helpful suggestions for common failures
        error_msg = str(e)
        if "not found" in error_msg.lower() or "component" in error_msg.lower():
            raise ValueError(
                f"Failed to balance: {error_msg}. "
                f"Please check that all chemical formulas are valid and "
                f"that elements on both sides of the equation match."
            ) from e
        raise ValueError(f"Failed to balance: {error_msg}") from e

    def format_side(coeffs):
        parts = []
        for substance, coeff in coeffs.items():
            parts.append(f"{coeff} {substance}" if coeff != 1 else str(substance))
        return " + ".join(parts)

    balanced_eq = f"{format_side(balanced_r)} -> {format_side(balanced_p)}"

    result = {
        "original": equation,
        "balanced": balanced_eq,
        "reactant_coefficients": {str(k): _parse_int_coeff(v) for k, v in balanced_r.items()},
        "product_coefficients": {str(k): _parse_int_coeff(v) for k, v in balanced_p.items()},
    }
    return json.dumps(result, indent=2, ensure_ascii=False)
