#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
graph_metrics_v2.py — Figuras de diversidad y composición (solo Python/conda)

Genera:
  - alfa_diversity_violin.html/png  (Shannon, Simpson, Observed, Chao1)   [ya existente]
  - alpha_observed_chao1.html/png   [NUEVO: Observed y Chao1 en una sola figura]
  - alpha_shannon_simpson.html/png  [NUEVO: Shannon y Simpson en otra figura]
  - braycurtis_heatmap_clust.html/png  (y jaccard_heatmap_clust)
  - pcoa_braycurtis.html/png  (y pcoa_jaccard) con ellipses por grupo (si hay metadata)
  - top_genus_stackedbar.html/png
  - top_genus_heatmap_clust.html/png
  - cooccurrence_network.html/png  (opcional, si hay suficientes taxones y networkx)

Requisitos (conda-forge):
  conda install -c conda-forge pandas numpy plotly scipy networkx kaleido

Usa kaleido si está instalado para imágenes estáticas. Siempre escribe HTML.
"""

from __future__ import annotations
import argparse
import os, re, sys
from pathlib import Path
from typing import Optional, Dict, List, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# SciPy para ordenar heatmaps por clustering
try:
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

# networkx para red de co-ocurrencia (inspirado en SpiecEasi)
try:
    import networkx as nx
    _HAS_NX = True
except Exception:
    _HAS_NX = False

# export estático
_HAS_KALEIDO = False
try:
    import kaleido  # noqa: F401
    _HAS_KALEIDO = True
except Exception:
    pass


def eprint(*a, **k): print(*a, file=sys.stderr, **k)


def _ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def _read_tsv(path: Path, **kw) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except Exception as e:
        eprint(f"[WARN] No pude leer {path}: {e}")
        return None


def _maybe_features_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Asegura orientación: filas=features, columnas=muestras."""
    idx_like_sample = pd.Series(df.index.astype(str)).str.match(r"^(SRR|ERR|DRR)\d+").any()
    if idx_like_sample:
        return df.T
    return df


def _save(fig: go.Figure, html_path: Path, png_path: Optional[Path] = None, width=1100, height=700):
    fig.update_layout(template="plotly_white", width=width, height=height)
    fig.write_html(str(html_path))
    if png_path and _HAS_KALEIDO:
        ext = png_path.suffix.lower()
        if ext in [".png", ".svg", ".pdf"]:
            fig.write_image(str(png_path))
    elif png_path:
        eprint(f"[INFO] Para exportar estático instala 'kaleido'. Me quedo solo con HTML: {html_path}")


def load_inputs(results_dir: Path):
    paths = {
        "alpha": results_dir / "alpha_diversity.tsv",
        "rel": results_dir / "feature_table_relabund.tsv",
        "tax": results_dir / "taxonomy.tsv",
        "bc_dm": results_dir / "braycurtis_dm.tsv",
        "ja_dm": results_dir / "jaccard_dm.tsv",
        "bc_coords": results_dir / "pcoa_braycurtis_coords.tsv",
        "bc_var": results_dir / "pcoa_braycurtis_variance.tsv",
        "ja_coords": results_dir / "pcoa_jaccard_coords.tsv",
        "ja_var": results_dir / "pcoa_jaccard_variance.tsv",
    }
    data = {k: _read_tsv(p, index_col=0 if k.endswith("_dm") or k.endswith("_coords") or k in ["rel"] else None)
            if p.exists() else None
            for k, p in paths.items()}
    return data, paths


def load_metadata(meta_path: Optional[Path], id_col: str) -> Optional[pd.DataFrame]:
    if not meta_path:
        return None
    df = _read_tsv(meta_path)
    if df is None:
        return None
    if id_col not in df.columns:
        cand = [c for c in df.columns if re.search(r"sample|id", c, re.I)]
        if cand:
            df.rename(columns={cand[0]: id_col}, inplace=True)
        else:
            eprint(f"[WARN] Metadata no tiene columna '{id_col}'. Continuo sin agrupar.")
            return None
    df[id_col] = df[id_col].astype(str)
    return df


# ---------- ALFA: violin + jitter (figura general con 4 métricas) ----------
def plot_alpha(alpha_df: pd.DataFrame, meta: Optional[pd.DataFrame], group_col: Optional[str],
               figdir: Path):
    if alpha_df is None or alpha_df.empty:
        return
    alpha = alpha_df.copy()
    if "sample" in alpha.columns:
        alpha.set_index("sample", inplace=True)
    alpha.index.name = "sample"
    alpha.reset_index(inplace=True)

    # largo
    al = alpha.melt(id_vars="sample", var_name="metric", value_name="value")

    if meta is not None and group_col and group_col in meta.columns:
        al = al.merge(meta[[group_col]].rename_axis("sample").reset_index(), on="sample", how="left")
    else:
        al[group_col or "Group"] = "All"

    fig = px.violin(
        al, x="metric", y="value", color=group_col or "Group", box=True, points="all",
        hover_data=["sample"], title="Alpha diversidad (violin + puntos)"
    )
    _save(fig, figdir / "alfa_diversity_violin.html", figdir / "alfa_diversity_violin.png")


# ---------- NUEVO: pares Observed+Chao1 y Shannon+Simpson en figuras separadas ----------
def _alpha_two_metric_violin(alpha_df: pd.DataFrame, m1: str, m2: str, title: str) -> go.Figure:
    """Crea una sola figura con dos violines superpuestos (m1 y m2)."""
    df = alpha_df.copy()
    # estandariza índice y nombres
    if "sample" in df.columns:
        df = df.set_index("sample")
    df.columns = [c.lower() for c in df.columns]
    need = {m1, m2}
    if not need.issubset(df.columns):
        faltan = sorted(list(need.difference(df.columns)))
        raise ValueError(f"Faltan columnas en alpha_diversity.tsv para este plot: {faltan}")
    y1 = pd.to_numeric(df[m1], errors="coerce").dropna()
    y2 = pd.to_numeric(df[m2], errors="coerce").dropna()
    fig = go.Figure()
    fig.add_trace(go.Violin(
        y=y1, name=m1.capitalize(),
        side="negative", box_visible=True, meanline_visible=True,
        points="all", jitter=0.2, scalemode="count", spanmode="hard"
    ))
    fig.add_trace(go.Violin(
        y=y2, name=m2.capitalize(),
        side="positive", box_visible=True, meanline_visible=True,
        points="all", jitter=0.2, scalemode="count", spanmode="hard"
    ))
    fig.update_traces(marker=dict(opacity=0.6, size=5))
    fig.update_layout(
        title=title,
        xaxis_title="Métrica",
        yaxis_title="Valor",
        violingap=0, violinmode="overlay",
        legend_title="Métrica"
    )
    return fig


def plot_alpha_pairs(alpha_df: pd.DataFrame, figdir: Path):
    """Genera:
       - alpha_observed_chao1.(html/png)
       - alpha_shannon_simpson.(html/png)
    """
    if alpha_df is None or alpha_df.empty:
        return
    # 1) Observed + Chao1
    try:
        fig1 = _alpha_two_metric_violin(alpha_df, "observed", "chao1",
                                        "Riqueza: Observed vs Chao1")
        _save(fig1, figdir / "alpha_observed_chao1.html", figdir / "alpha_observed_chao1.png")
    except Exception as e:
        eprint(f"[WARN] No se pudo crear alpha_observed_chao1: {e}")

    # 2) Shannon + Simpson
    try:
        fig2 = _alpha_two_metric_violin(alpha_df, "shannon", "simpson",
                                        "Diversidad: Shannon vs Simpson")
        _save(fig2, figdir / "alpha_shannon_simpson.html", figdir / "alpha_shannon_simpson.png")
    except Exception as e:
        eprint(f"[WARN] No se pudo crear alpha_shannon_simpson: {e}")


# ---------- HEATMAP DM con clustering ----------
def _heatmap_with_clustering(D: pd.DataFrame, title: str) -> go.Figure:
    ordered = D
    if _HAS_SCIPY and D.shape[0] >= 3:
        try:
            A = D.values
            # Si es una matriz de distancia cuadrada (simétrica y diagonal≈0), usar squareform -> linkage
            if D.shape[0] == D.shape[1] and np.allclose(A, A.T, atol=1e-12) and np.allclose(np.diag(A), 0.0, atol=1e-12):
                Z = linkage(squareform(A, checks=False), method="average")
                order = leaves_list(Z)
                ordered = D.iloc[order, :].iloc[:, order]
            else:
                # Si no es DM, tratamos filas como observaciones
                Z = linkage(A, method="average")
                order = leaves_list(Z)
                ordered = D.iloc[order, :].iloc[:, order]
        except Exception as e:
            eprint(f"[WARN] No se pudo clusterizar heatmap: {e}")
            ordered = D

    fig = go.Figure(
        data=go.Heatmap(
            z=ordered.values,
            x=ordered.columns.astype(str),
            y=ordered.index.astype(str),
            colorbar=dict(title="distancia")
        )
    )
    fig.update_layout(title=title, xaxis_title="Muestras", yaxis_title="Muestras")
    return fig


def plot_heatmaps(dm_bc: Optional[pd.DataFrame], dm_ja: Optional[pd.DataFrame], figdir: Path):
    if dm_bc is not None and not dm_bc.empty:
        fig = _heatmap_with_clustering(dm_bc, "Bray–Curtis — heatmap/clust")
        _save(fig, figdir / "braycurtis_heatmap_clust.html", figdir / "braycurtis_heatmap_clust.png")
    if dm_ja is not None and not dm_ja.empty:
        fig = _heatmap_with_clustering(dm_ja, "Jaccard — heatmap/clust")
        _save(fig, figdir / "jaccard_heatmap_clust.html", figdir / "jaccard_heatmap_clust.png")


# ---------- PCoA con “ellipses” por grupo ----------
def _group_ellipses(df2d: pd.DataFrame, group_col: str) -> List[go.Scatter]:
    traces = []
    for g, sub in df2d.groupby(group_col):
        if sub.shape[0] < 3:
            continue
        xy = sub[["PC1", "PC2"]].values
        mu = xy.mean(axis=0)
        cov = np.cov(xy, rowvar=False)
        try:
            vals, vecs = np.linalg.eigh(cov)
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:, order]
            theta = np.linspace(0, 2*np.pi, 100)
            circle = np.array([np.cos(theta), np.sin(theta)])
            scale = 2.0 * np.sqrt(vals)  # ~95%
            ellipse = (vecs @ (scale[:, None] * circle)) + mu[:, None]
            traces.append(go.Scatter(
                x=ellipse[0], y=ellipse[1], mode="lines", name=f"{g} ellipse",
                line=dict(dash="dash")
            ))
        except Exception:
            pass
    return traces


def plot_pcoa(coords: Optional[pd.DataFrame], var: Optional[pd.DataFrame],
              meta: Optional[pd.DataFrame], group_col: Optional[str],
              prefix: str, figdir: Path):
    if coords is None or coords.empty:
        return
    cs = coords.copy()
    if "PC1" not in cs.columns or "PC2" not in cs.columns:
        cs.columns = [f"PC{i+1}" for i in range(cs.shape[1])]
    df = cs.reset_index().rename(columns={"index": "sample"})
    if meta is not None and group_col and group_col in meta.columns:
        df = df.merge(meta[[group_col]].rename_axis("sample").reset_index(), on="sample", how="left")
    else:
        df[group_col or "Group"] = "All"

    pc1lab = "PC1"
    pc2lab = "PC2"
    if var is not None and "proportion_explained" in var.columns and len(var) >= 2:
        ve = var["proportion_explained"].values
        pc1lab += f" ({ve[0]*100:.1f}%)"
        pc2lab += f" ({ve[1]*100:.1f}%)"

    fig = px.scatter(
        df, x="PC1", y="PC2", color=group_col or "Group",
        hover_name="sample", title=f"PCoA — {prefix}"
    )
    if group_col and group_col in df.columns:
        for tr in _group_ellipses(df, group_col):
            fig.add_trace(tr)
    fig.update_layout(xaxis_title=pc1lab, yaxis_title=pc2lab)
    _save(fig, figdir / f"pcoa_{prefix}.html", figdir / f"pcoa_{prefix}.png")


# ---------- Composición al nivel género ----------
def _parse_genus_from_qiime_tax(s: str) -> str:
    m = re.search(r"g__([^;]+)", str(s))
    g = m.group(1) if m else "Unassigned"
    g = "Unassigned" if g in ["", "uncultured", "metagenome"] else g
    return g

def genus_tables(rel: pd.DataFrame, tax: pd.DataFrame, topn: int = 12) -> Tuple[pd.DataFrame, List[str]]:
    if rel is None or rel.empty or tax is None or tax.empty:
        return pd.DataFrame(), []
    rel0 = _maybe_features_rows(rel.copy())
    rel0 = rel0.fillna(0.0)
    if "feature_id" in tax.columns:
        tx = tax.set_index("feature_id")
    else:
        tx = tax.copy()
        tx.index.name = "feature_id"
    if "taxonomy" not in tx.columns:
        raise RuntimeError("taxonomy.tsv debe tener columna 'taxonomy' con cadena QIIME.")

    genus = tx["taxonomy"].map(_parse_genus_from_qiime_tax)
    common = rel0.index.intersection(genus.index)
    if len(common) == 0:
        rel0 = rel0.T
        common = rel0.index.intersection(genus.index)
    rel0 = rel0.loc[common]
    genus = genus.loc[common]

    G = rel0.groupby(genus).sum()
    colsum = G.sum(axis=0)
    colsum[colsum == 0] = 1.0
    G = G / colsum

    top = (G.mean(axis=1)).sort_values(ascending=False).head(topn).index.tolist()
    Gtop = G.loc[top].copy()
    if len(G.index.difference(top)) > 0:
        others = G.loc[G.index.difference(top)].sum()
        Gtop.loc["Others"] = others
    return Gtop, top


def plot_genus_stacked(G: pd.DataFrame, figdir: Path):
    if G is None or G.empty:
        return
    df = G.T.reset_index().rename(columns={"index": "sample"})
    df_long = df.melt(id_vars="sample", var_name="genus", value_name="rel_abund")
    fig = px.bar(
        df_long, x="sample", y="rel_abund", color="genus",
        title="Abundancia relativa — top géneros (apilado)",
        hover_data={"rel_abund":":.3f"}
    )
    fig.update_layout(yaxis_title="Proporción", xaxis_tickangle=45, legend_title="Género")
    _save(fig, figdir / "top_genus_stackedbar.html", figdir / "top_genus_stackedbar.png", height=800)


def plot_genus_heatmap(G: pd.DataFrame, figdir: Path):
    if G is None or G.empty:
        return
    M = G.copy()
    if _HAS_SCIPY and M.shape[0] >= 3 and M.shape[1] >= 3:
        try:
            rorder = leaves_list(linkage(M.values, method="average"))
            corder = leaves_list(linkage(M.values.T, method="average"))
            M = M.iloc[rorder, :].iloc[:, corder]
        except Exception:
            pass
    fig = go.Figure(
        data=go.Heatmap(
            z=M.values, x=M.columns.astype(str), y=M.index.astype(str),
            colorbar=dict(title="rel. abundancia")
        )
    )
    fig.update_layout(title="Heatmap — top géneros", xaxis_title="Muestras", yaxis_title="Género")
    _save(fig, figdir / "top_genus_heatmap_clust.html", figdir / "top_genus_heatmap_clust.png")


# ---------- Red de co-ocurrencia (Spearman) ----------
def plot_cooccurrence_network(G: pd.DataFrame, figdir: Path,
                              corr_thr: float = 0.5, top_edges: int = 200):
    if not _HAS_NX or G is None or G.empty:
        return
    X = G.values
    corr = pd.DataFrame(np.corrcoef(X), index=G.index, columns=G.index)
    corr = G.T.corr(method="spearman")

    A = corr.where(np.triu(np.ones_like(corr, dtype=bool), 1))
    edges = (
        A.stack()
         .rename("rho")
         .reindex(A.columns, level=0)
         .dropna()
         .sort_values(key=lambda s: s.abs(), ascending=False)
    )
    if edges.empty:
        return
    edges = edges.iloc[:min(len(edges), top_edges)]
    edges_keep = edges[edges.abs() >= corr_thr]

    if edges_keep.empty:
        edges_keep = edges.iloc[:min(len(edges), 50)]

    Gx = nx.Graph()
    for g in G.index:
        Gx.add_node(g, size=float(G.loc[g].mean()))
    for (g1, g2), r in edges_keep.items():
        Gx.add_edge(g1, g2, weight=float(abs(r)), sign=np.sign(r))

    pos = nx.spring_layout(Gx, seed=42, k=0.5)

    edgelines = []
    for (u, v, d) in Gx.edges(data=True):
        x0, y0 = pos[u]; x1, y1 = pos[v]
        dash = "solid" if d["sign"] >= 0 else "dot"
        edgelines.append(go.Scatter(x=[x0, x1], y=[y0, y1],
                                    mode="lines", line=dict(width=1+3*d["weight"], dash=dash),
                                    hoverinfo="skip", showlegend=False))

    nodelines = go.Scatter(
        x=[pos[n][0] for n in Gx.nodes()],
        y=[pos[n][1] for n in Gx.nodes()],
        mode="markers+text",
        text=[str(n) for n in Gx.nodes()],
        textposition="top center",
        marker=dict(size=[8 + 40*Gx.nodes[n]["size"] for n in Gx.nodes()]),
        hovertext=[f"{n}<br>mean rel.abund={Gx.nodes[n]['size']:.3f}" for n in Gx.nodes()],
        hoverinfo="text",
        name="géneros"
    )

    fig = go.Figure(edgelines + [nodelines])
    fig.update_layout(title="Red de co-ocurrencia (ρ Spearman) — inspirado en SPIEC-EASI",
                      xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                      yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
    _save(fig, figdir / "cooccurrence_network.html", figdir / "cooccurrence_network.png")


# ---------- CLI / MAIN ----------
def build_parser():
    p = argparse.ArgumentParser(
        description="Graficador de métricas a partir de salidas TSV (EMU pipeline)."
    )
    p.add_argument("--results-dir", required=True, type=Path,
                   help="Carpeta que contiene alpha_diversity.tsv, feature_table_relabund.tsv, taxonomy.tsv, etc.")
    p.add_argument("--metadata", type=Path, default=None,
                   help="TSV opcional con metadata por muestra (debe contener columna --id-col).")
    p.add_argument("--id-col", type=str, default="sample", help="Nombre de la columna de IDs de muestra en metadata.")
    p.add_argument("--group-col", type=str, default=None, help="Columna de metadata para colorear/ellipses (p.ej. Grupo).")
    p.add_argument("--topn-genus", type=int, default=12, help="Número de géneros top para plots de composición.")
    p.add_argument("--figdir", type=Path, default=None, help="Carpeta de figuras (por defecto results-dir/figures).")
    return p

def main():
    args = build_parser().parse_args()
    results_dir = args.results_dir.resolve()
    if not results_dir.exists():
        sys.exit(f"[ERROR] No existe --results-dir: {results_dir}")

    figdir = _ensure_dir(args.figdir or (results_dir / "figures"))

    data, paths = load_inputs(results_dir)
    meta = load_metadata(args.metadata, id_col=args.id_col)

    # ALFA (figura general con 4 métricas)
    if data["alpha"] is not None:
        plot_alpha(data["alpha"], meta, args.group_col, figdir)
        # NUEVO: pares Observed+Chao1 y Shannon+Simpson
        plot_alpha_pairs(data["alpha"], figdir)
    else:
        eprint("[INFO] Sin alpha_diversity.tsv — salto.")

    # HEATMAPS DM
    if data["bc_dm"] is not None or data["ja_dm"] is not None:
        plot_heatmaps(data["bc_dm"], data["ja_dm"], figdir)
    else:
        eprint("[INFO] Sin matrices de distancia — salto heatmaps.")

    # PCoA
    if data["bc_coords"] is not None:
        plot_pcoa(data["bc_coords"], data["bc_var"], meta, args.group_col, "braycurtis", figdir)
    if data["ja_coords"] is not None:
        plot_pcoa(data["ja_coords"], data["ja_var"], meta, args.group_col, "jaccard", figdir)

    # COMPOSICIÓN
    if data["rel"] is not None and data["tax"] is not None:
        Gtop, _ = genus_tables(data["rel"], data["tax"], topn=args.topn_genus)
        if not Gtop.empty:
            plot_genus_stacked(Gtop, figdir)
            plot_genus_heatmap(Gtop, figdir)
            try:
                plot_cooccurrence_network(Gtop, figdir)
            except Exception as e:
                eprint(f"[WARN] No se pudo construir la red de co-ocurrencia: {e}")
    else:
        eprint("[INFO] Faltan feature_table_relabund.tsv o taxonomy.tsv — salto composición.")

    eprint(f"[OK] Figuras en: {figdir}")

if __name__ == "__main__":
    main()

