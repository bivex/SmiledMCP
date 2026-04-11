"""Chemistry MCP Server — FastMCP instance and entry point."""

import sys

from mcp.server.fastmcp import FastMCP

mcp = FastMCP("chemistry-mcp-server")

# When running as `python server.py`, the module is __main__ but tool modules
# import `from server import mcp`.  Without this alias, Python loads server.py
# a second time under the name "server", creating a *different* FastMCP instance
# — tools register on that ghost instance while .run() uses the original (empty) one.
sys.modules.setdefault("server", sys.modules[__name__])


def main():
    # Import tools to trigger registration
    import tools  # noqa: F401
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
