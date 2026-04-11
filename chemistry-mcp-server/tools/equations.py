"""Equation balancing tools — chempy stoichiometry."""

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

    try:
        balanced_r, balanced_p = balance_stoichiometry(reactants, products)
    except Exception as e:
        raise ValueError(f"Failed to balance: {e}") from e

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
