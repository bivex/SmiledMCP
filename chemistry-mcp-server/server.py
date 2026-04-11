"""Chemistry MCP Server — FastMCP instance and entry point."""

from mcp.server.fastmcp import FastMCP

mcp = FastMCP("chemistry-mcp-server")


def main():
    # Import tools to trigger registration
    import tools  # noqa: F401
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
