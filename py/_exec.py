"""Rendering helpers for the executed documentation examples.

Not part of the ``pilotr`` package. The docs pages run their code at build time via
markdown-exec; ``table`` prints a ``pilotr`` result as a Markdown table and ``show`` prints a
matplotlib figure as inline SVG. Kept dependency-light: matplotlib is a docs-build dependency
only.
"""

import matplotlib
matplotlib.use("Agg")  # headless backend for the CI build
import matplotlib.pyplot as plt
from io import StringIO

# The pilotr palette, shared with the R package website.
BLUE = "#2c6fb0"
RED = "#b0402c"
GREEN = "#2e8b57"

plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "font.size": 11})


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.3g}"
    return str(v)


# markdown-exec injects its own ``print`` into each executed block, so these helpers must
# *return* their rendered text (the block then prints it) rather than printing directly.

def table(rows, columns=None):
    """Return a list of row dicts as a GitHub-flavoured Markdown table."""
    rows = list(rows)
    if columns is None:
        columns = list(rows[0].keys())
    lines = ["", "| " + " | ".join(columns) + " |",
             "| " + " | ".join("---" for _ in columns) + " |"]
    for r in rows:
        lines.append("| " + " | ".join(_fmt(r.get(c, "")) for c in columns) + " |")
    return "\n".join(lines)


def show(fig):
    """Return a matplotlib figure as an inline SVG string (use with ``html="true"``)."""
    buf = StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()
