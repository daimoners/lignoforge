"""
Interactive 3D visualization helpers for coarse-grained topology JSON.
"""

from __future__ import annotations

import json
from typing import Optional

# Monomer residue codes → bead colour
_BEAD_COLORS: dict = {
    "HPU": "#e74c3c",   # H – p-Hydroxyphenyl – red
    "GYU": "#2ecc71",   # G – Guaiacyl          – green
    "SYU": "#3498db",   # S – Syringyl           – blue
    "UNK": "#95a5a6",   # unknown                – grey
}

# Linkage type → bond colour  (C-O types blue/teal; C-C cross-links warm)
_LINK_COLORS: dict = {
    "beta-O-4":  "#1abc9c",
    "alpha-O-4": "#16a085",
    "4-O-5":     "#8e44ad",
    "5-5":       "#e74c3c",
    "beta-5":    "#e67e22",
    "beta-beta": "#f39c12",
    "beta-1":    "#d35400",
    "unknown":   "#95a5a6",
}

# C-C cross-links and ring-closing bonds are rendered dashed; C-O backbone solid
_LINK_DASH: dict = {
    "beta-O-4":  "solid",
    "alpha-O-4": "solid",
    "4-O-5":     "dash",
    "5-5":       "dash",
    "beta-5":    "dash",
    "beta-beta": "dash",
    "beta-1":    "dot",
    "unknown":   "solid",
}

_MONOMER_LABEL: dict = {
    "HPU": "H (p-Hydroxyphenyl)",
    "GYU": "G (Guaiacyl)",
    "SYU": "S (Syringyl)",
}


def write_cg_topology_viewer_html(
    cg_json_path: str,
    html_output_path: str,
    title: str = "LignoForge Coarse-Grained Topology",
) -> str:
    """
    Read a coarse-grained topology JSON and generate an interactive 3D HTML viewer.

    * Beads are colour-coded by monomer type: H (red), G (green), S (blue).
    * Inter-monomer bonds are colour-coded by linkage type.
    * C-C cross-links (5-5, beta-5, beta-beta, beta-1, 4-O-5) are dashed lines;
      backbone C-O bonds (beta-O-4, alpha-O-4) are solid.
    """
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError(
            "plotly is required for interactive CG visualization. "
            "Install with `pip install plotly`."
        ) from exc

    with open(cg_json_path, "r") as f:
        data = json.load(f)

    chains = data.get("chains", [data])

    fig = go.Figure()

    # Track which legend labels have already been registered (shared across chains)
    legend_beads_registered: set = set()
    legend_links_registered: set = set()

    for chain in chains:
        chain_id = chain.get("chain_id", "A")
        beads = chain.get("beads", [])
        links = chain.get("links", [])

        if not beads:
            continue

        bead_xyz: dict = {b["bead_id"]: (b["x"], b["y"], b["z"]) for b in beads}

        # ── Beads: one trace per monomer type for a clean legend ──────────────
        type_groups: dict = {}
        for b in beads:
            btype = b.get("bead_type") or "UNK"
            type_groups.setdefault(btype, []).append(b)

        for btype, bead_list in sorted(type_groups.items()):
            color = _BEAD_COLORS.get(btype, "#95a5a6")
            show_in_legend = btype not in legend_beads_registered
            legend_beads_registered.add(btype)

            fig.add_trace(go.Scatter3d(
                x=[b["x"] for b in bead_list],
                y=[b["y"] for b in bead_list],
                z=[b["z"] for b in bead_list],
                mode="markers",
                marker=dict(
                    size=12,
                    color=color,
                    line=dict(width=1, color="white"),
                    opacity=0.9,
                ),
                text=[
                    f"Chain {chain_id} | Bead {b['bead_id']} | {btype}"
                    for b in bead_list
                ],
                hovertemplate="%{text}<extra></extra>",
                name=_MONOMER_LABEL.get(btype, btype),
                legendgroup=btype,
                showlegend=show_in_legend,
            ))

        # ── Bonds: one trace per linkage type ─────────────────────────────────
        link_groups: dict = {}
        for lk in links:
            ltype = lk.get("linkage_type") or "unknown"
            link_groups.setdefault(ltype, []).append(lk)

        for ltype, link_list in sorted(link_groups.items()):
            x_edges: list = []
            y_edges: list = []
            z_edges: list = []
            for lk in link_list:
                s = lk["source_bead_id"]
                t = lk["target_bead_id"]
                if s not in bead_xyz or t not in bead_xyz:
                    continue
                x_edges += [bead_xyz[s][0], bead_xyz[t][0], None]
                y_edges += [bead_xyz[s][1], bead_xyz[t][1], None]
                z_edges += [bead_xyz[s][2], bead_xyz[t][2], None]

            if not x_edges:
                continue

            color = _LINK_COLORS.get(ltype, "#95a5a6")
            dash  = _LINK_DASH.get(ltype, "solid")
            show_in_legend = ltype not in legend_links_registered
            legend_links_registered.add(ltype)

            fig.add_trace(go.Scatter3d(
                x=x_edges,
                y=y_edges,
                z=z_edges,
                mode="lines",
                line=dict(width=4, color=color, dash=dash),
                name=ltype,
                legendgroup=ltype,
                showlegend=show_in_legend,
                hoverinfo="none",
            ))

    fig.update_layout(
        title=dict(text=title, font=dict(size=16)),
        scene=dict(
            xaxis_title="X (Å)",
            yaxis_title="Y (Å)",
            zaxis_title="Z (Å)",
            xaxis=dict(showgrid=True, gridcolor="#ecf0f1"),
            yaxis=dict(showgrid=True, gridcolor="#ecf0f1"),
            zaxis=dict(showgrid=True, gridcolor="#ecf0f1"),
            bgcolor="#f8f9fa",
        ),
        legend=dict(
            title=dict(text="Monomer / Linkage"),
            itemsizing="constant",
        ),
        margin=dict(l=0, r=0, t=50, b=0),
    )

    fig.write_html(html_output_path, include_plotlyjs="cdn")
    return html_output_path
