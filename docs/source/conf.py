# Configuration file for the Sphinx documentation builder.
# LignoForge v0.1.0

import os
import sys

# Let Sphinx find the lignoforge package for autodoc
_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, _root)

# ── Project info ──────────────────────────────────────────────────────────────
project   = "LignoForge"
copyright = "2025, DAIMON Team"
author    = "DAIMON Team"
version   = "0.1"
release   = "0.1.0"

# ── Sphinx extensions ─────────────────────────────────────────────────────────
extensions = [
    "sphinx.ext.autodoc",       # extract docstrings automatically
    "sphinx.ext.autosummary",   # class/function summary tables
    "sphinx.ext.napoleon",      # Google + NumPy docstring styles
    "sphinx.ext.intersphinx",   # cross-link to Python / NumPy docs
    "sphinx.ext.mathjax",       # LaTeX math rendering
    "sphinx.ext.viewcode",      # [source] links next to API entries
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "sphinx_design",            # grid cards, tabs, badges
]

# autodoc
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members":           True,
    "undoc-members":     False,
    "show-inheritance":  True,
    "special-members":   "__init__",
}

# napoleon
napoleon_google_docstring = True
napoleon_numpy_docstring  = True
napoleon_include_init_with_doc = True
napoleon_use_ivar  = True
napoleon_use_param = True
napoleon_use_rtype = True

# ── General ───────────────────────────────────────────────────────────────────
templates_path   = ["_templates"]
source_suffix    = ".rst"
master_doc       = "index"
language         = "en"
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
pygments_style   = "friendly"

# ── HTML output ───────────────────────────────────────────────────────────────
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth":    4,
    "collapse_navigation": False,
    "sticky_navigation":   True,
    "titles_only":         False,
}
html_static_path = ["_static"]
html_css_files   = ["css/custom.css"]
html_title       = "LignoForge Documentation"
html_short_title = "LignoForge"

# ── intersphinx ───────────────────────────────────────────────────────────────
intersphinx_mapping = {
    "python":   ("https://docs.python.org/3", None),
    "numpy":    ("https://numpy.org/doc/stable/", None),
    "networkx": ("https://networkx.org/documentation/stable/", None),
}

todo_include_todos = True
