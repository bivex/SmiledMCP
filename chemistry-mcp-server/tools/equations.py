"""Equation balancing tools — chempy stoichiometry."""

import json

from server import mcp

try:
    from chempy import balance_stoichiometry
    HAS_CHEMPY = True
except ImportError:
    HAS_CHEMPY = False


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
