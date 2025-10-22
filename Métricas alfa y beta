#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
metricas.py — Calcula diversidad alfa/beta + PCoA y genera figuras a partir de salidas de EMU.

Entradas esperadas en --emu-outdir:
  - feature_table_counts.tsv           (preferida; features x muestras)
  - feature_table_relabund.tsv         (respaldo si no hay counts; features x muestras)
  - taxonomy.tsv                       (feature_id, taxonomy estilo QIIME; opcional para barras por género)
  - metadata.tsv                       (opcional vía --metadata; filas=muestras)

Salidas escritas en --emu-outdir:
  - alpha_diversity.tsv
  - braycurtis_dm.tsv, jaccard_dm.tsv
  - pcoa_braycurtis_coords.tsv / _eigvals.tsv / _variance.tsv (+ pcoa_braycurtis.png / .html)
  - pcoa_jaccard_coords.tsv  / _eigvals.tsv  / _variance.tsv  (+ pcoa_jaccard.png  / .html)
  - plots_static/bars_topN_genus.png (si taxonomy y --plots-static)
  - plots_interactive/bars_topN_genus.html (si taxonomy y --plots-interactive)
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Matplotlib (figuras estáticas)
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

# Plotly (figuras interactivas)
try:
    import plotly.express as px
    HAS_PLOTLY = True
except Exception:
    HAS_PLOTLY = False


# ----------------------------- utilidades básicas -----------------------------

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def _read_tsv(path: Path, index_col: Optional[int | str] = 0) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=index_col)
    return df


def _to_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.fillna(0.0)
    return out


def _looks_like_sample(name: str) -> bool:
    return bool(re.match(r"^(SRR|ERR|DRR)\d+", str(name)))


def _ensure_features_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Garantiza orientación: filas=features, columnas=muestras.
    Heurística:
      - Si columnas parecen muestras (SRR/ERR/DRR), se deja como está.
      - Si el índice parece muestras, se transpone.
      - Si no hay pista, se deja como está.
    """
    cols_samples = sum(_looks_like_sample(c) for c in df.columns)
    idx_samples = sum(_looks_like_sample(i) for i in df.index)
    if idx_samples > cols_samples:
        return df.T
    return df


def _align_with_metadata(df: pd.DataFrame, meta: Optional[pd.DataFrame]) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    if meta is None:
        return df, None
    common = [c for c in df.columns if c in meta.index]
    if len(common) < 2:
        eprint("[WARN] Menos de 2 muestras en común entre tabla y metadata; se ignorará metadata.")
        return df, None
    df2 = df[common].copy()
    meta2 = meta.loc[common].copy()
    return df2, meta2


def _drop_all_zero_features(df: pd.DataFrame) -> pd.DataFrame:
    mask = df.sum(axis=1) > 0
    dropped = (~mask).sum()
    if dropped > 0:
        eprint(f"[INFO] Eliminadas {dropped} features con suma total 0.")
    return df[mask]


def _drop_zero_samples(df: pd.DataFrame) -> pd.DataFrame:
    mask = df.sum(axis=0) > 0
    dropped = (~mask).sum()
    if dropped > 0:
        eprint(f"[INFO] Eliminadas {dropped} muestras con suma total 0.")
    return df.loc[:, mask]


# ------------------------------- carga de datos -------------------------------

def load_feature_tables(emu_outdir: Path) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    counts_p = emu_outdir / "feature_table_counts.tsv"
    rel_p    = emu_outdir / "feature_table_relabund.tsv"

    counts_df = None
    rel_df = None

    if counts_p.exists():
        df = _read_tsv(counts_p, index_col=0)
        df = _ensure_features_rows(df)
        df = _to_numeric_df(df)
        counts_df = df
        eprint(f"[OK] Cargada tabla de conteos: {counts_p}")
    else:
        eprint("[WARN] No se encontró feature_table_counts.tsv")

    if rel_p.exists():
        df = _read_tsv(rel_p, index_col=0)
        df = _ensure_features_rows(df)
        df = _to_numeric_df(df)
        rel_df = df
        eprint(f"[OK] Cargada tabla de abundancias relativas: {rel_p}")
    else:
        eprint("[WARN] No se encontró feature_table_relabund.tsv")

    if counts_df is None and rel_df is None:
        raise FileNotFoundError("No se encontraron ni 'feature_table_counts.tsv' ni 'feature_table_relabund.tsv'.")

    return counts_df, rel_df


def load_taxonomy(emu_outdir: Path) -> Optional[pd.DataFrame]:
    tax_p = emu_outdir / "taxonomy.tsv"
    if not tax_p.exists():
        eprint("[WARN] taxonomy.tsv no encontrado; se omitirán barras por género.")
        return None
    df = pd.read_csv(tax_p, sep="\t")
    if "feature_id" not in df.columns or "taxonomy" not in df.columns:
        eprint("[WARN] taxonomy.tsv no tiene columnas esperadas; se omiten barras por género.")
        return None
    df = df[["feature_id", "taxonomy"]].dropna()
    return df


def _load_metadata(meta_path: Optional[Path]) -> Optional[pd.DataFrame]:
    if meta_path is None:
        return None
    if not meta_path.exists():
        eprint(f"[WARN] Metadata no encontrada: {meta_path}. Se continúa sin metadata.")
        return None
    df = pd.read_csv(meta_path, sep="\t")
    if df.empty:
        eprint("[WARN] Metadata vacía; se ignora.")
        return None
    # La primera columna que parezca ID de muestra se usa como índice; de lo contrario, usa la primera
    if "sample" in df.columns:
        df = df.set_index("sample")
    else:
        df = df.set_index(df.columns[0])
    return df


# ------------------------------ alfa diversidad -------------------------------

def chao1(counts: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=float)
    Sobs = np.sum(counts > 0)
    f1 = np.sum(counts == 1)
    f2 = np.sum(counts == 2)
    if f2 > 0:
        return float(Sobs + (f1 * f1) / (2.0 * f2))
    return float(Sobs + (f1 * (f1 - 1.0)) / 2.0)


def shannon_entropy(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p))) if p.size else 0.0


def simpson_index(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    return float(1.0 - np.sum(p * p))


def alpha_from_table(counts_df: Optional[pd.DataFrame], rel_df: Optional[pd.DataFrame]) -> pd.DataFrame:
    """
    Prefiere counts para Observed/Chao1.
    Si no hay counts, usa relativas; Observed=features>0, Chao1≈Observed.
    """
    if counts_df is not None:
        table = _drop_all_zero_features(_drop_zero_samples(counts_df))
        if table.shape[1] == 0:
            eprint("[WARN] Todas las muestras de counts tienen suma 0; se intentará con relativas.")
            table = None
    else:
        table = None

    if table is None:
        assert rel_df is not None
        table = _drop_all_zero_features(_drop_zero_samples(rel_df))
        use_counts = False
    else:
        use_counts = True

    rows = []
    for sample in table.columns:
        x = table[sample].values.astype(float)
        if use_counts:
            total = x.sum()
            p = x / total if total > 0 else x
            obs = int(np.sum(x > 0))
            c1 = chao1(x)
        else:
            p = x
            obs = int(np.sum(p > 0))
            c1 = float(obs)
        H = shannon_entropy(p)
        S = simpson_index(p)
        rows.append({"sample": sample, "observed": obs, "chao1": c1, "shannon": H, "simpson": S})

    return pd.DataFrame(rows).set_index("sample")


# ------------------------------ beta + PCoA -----------------------------------

def _bray_curtis_dm(X: np.ndarray) -> np.ndarray:
    n = X.shape[0]
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        xi = X[i]
        for j in range(i + 1, n):
            xj = X[j]
            num = np.abs(xi - xj).sum()
            den = (xi + xj).sum()
            D[i, j] = D[j, i] = 0.0 if den == 0 else num / den
    return D


def _jaccard_dm_bin(Xbin: np.ndarray) -> np.ndarray:
    n = Xbin.shape[0]
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        ai = Xbin[i] > 0
        for j in range(i + 1, n):
            aj = Xbin[j] > 0
            inter = np.logical_and(ai, aj).sum()
            union = np.logical_or(ai, aj).sum()
            D[i, j] = D[j, i] = 1.0 if union == 0 else 1.0 - (inter / union)
    return D


def _classical_pcoa(D: np.ndarray):
    n = D.shape[0]
    if n < 2:
        return np.zeros((n, 0)), np.array([]), np.array([])
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pos = eigvals > 1e-12
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]
    coords = eigvecs_pos * np.sqrt(eigvals_pos)
    var = eigvals_pos / eigvals_pos.sum() if eigvals_pos.size else np.array([])
    return coords, eigvals_pos, var


def _save_dm_and_pcoa(D: np.ndarray, sample_ids: List[str], outdir: Path, prefix: str,
                      make_png: bool, make_html: bool, meta: Optional[pd.DataFrame], group_col: Optional[str]):
    ensure_dir(outdir)
    # DM
    dm_path = outdir / f"{prefix}_dm.tsv"
    pd.DataFrame(D, index=sample_ids, columns=sample_ids).to_csv(dm_path, sep="\t")

    # PCoA
    coords, eigvals, var = _classical_pcoa(D)
    k = min(3, coords.shape[1]) if coords.size else 0

    coords_path = outdir / f"pcoa_{prefix}_coords.tsv"
    eig_path    = outdir / f"pcoa_{prefix}_eigvals.tsv"
    var_path    = outdir / f"pcoa_{prefix}_variance.tsv"

    pd.DataFrame(coords[:, :k], index=sample_ids,
                 columns=[f"PC{i+1}" for i in range(k)]).to_csv(coords_path, sep="\t")
    pd.DataFrame({"eigenvalue": eigvals}).to_csv(eig_path, sep="\t", index=False)
    pd.DataFrame({"proportion_explained": var}).to_csv(var_path, sep="\t", index=False)

    # Figuras
    if k >= 2:
        pc1 = f"{(var[0] * 100):.1f}%" if var.size > 0 else ""
        pc2 = f"{(var[1] * 100):.1f}%" if var.size > 1 else ""
        df_plot = pd.DataFrame({
            "sample": sample_ids,
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
        }).set_index("sample")

        color_series = None
        if meta is not None:
            # Elegir columna de agrupación
            col = group_col
            if col and col in meta.columns:
                color_series = meta[col]
            else:
                # Heurística: primera columna categórica con <=20 niveles
                for c in meta.columns:
                    if meta[c].dtype == "O" or meta[c].dtype.name.startswith("category"):
                        if meta[c].nunique() <= 20:
                            color_series = meta[c]
                            break
            if color_series is not None:
                # Alinear por índice de df_plot
                color_series = color_series.reindex(df_plot.index)

        # PNG
        if make_png and HAS_MPL:
            plt.figure()
            if color_series is None:
                plt.scatter(df_plot["PC1"], df_plot["PC2"])
            else:
                groups = color_series.fillna("NA").astype(str).values
                uniq = pd.unique(groups)
                for g in uniq:
                    sel = (groups == g)
                    plt.scatter(df_plot["PC1"][sel], df_plot["PC2"][sel], label=str(g))
                plt.legend(frameon=False, fontsize=8)
            for s, (x, y) in df_plot[["PC1", "PC2"]].iterrows():
                plt.text(x, y, str(s), fontsize=7)
            plt.xlabel(f"PC1 ({pc1})"); plt.ylabel(f"PC2 ({pc2})")
            plt.title(f"PCoA — {prefix}")
            plt.tight_layout()
            fig_dir = outdir / "plots_static"
            ensure_dir(fig_dir)
            plt.savefig(fig_dir / f"pcoa_{prefix}.png", dpi=150)
            plt.close()

        # HTML
        if make_html and HAS_PLOTLY:
            df_html = df_plot.reset_index()
            if color_series is not None:
                df_html["group"] = color_series.values
                fig = px.scatter(df_html, x="PC1", y="PC2", text="sample", color="group",
                                 title=f"PCoA — {prefix} (PC1 {pc1}, PC2 {pc2})")
            else:
                fig = px.scatter(df_html, x="PC1", y="PC2", text="sample",
                                 title=f"PCoA — {prefix} (PC1 {pc1}, PC2 {pc2})")
            fig.update_traces(textposition="top center")
            html_dir = outdir / "plots_interactive"
            ensure_dir(html_dir)
            fig.write_html(str(html_dir / f"pcoa_{prefix}.html"), include_plotlyjs="cdn")


def beta_and_pcoa(counts_df: Optional[pd.DataFrame],
                  rel_df: Optional[pd.DataFrame],
                  outdir: Path,
                  make_png: bool,
                  make_html: bool,
                  meta: Optional[pd.DataFrame],
                  group_col: Optional[str]):
    # Preparar matriz de trabajo
    if counts_df is not None:
        tab = _drop_all_zero_features(_drop_zero_samples(counts_df))
        if tab.shape[1] < 2:
            eprint("[WARN] <2 muestras con counts útiles; se intentará con relativas.")
            tab = None
    else:
        tab = None

    if tab is None:
        if rel_df is None:
            eprint("[WARN] No hay datos para beta-diversidad.")
            return
        tab = _drop_all_zero_features(_drop_zero_samples(rel_df))

    # Alinear con metadata
    tab, meta2 = _align_with_metadata(tab, meta)

    if tab.shape[1] < 2:
        eprint("[WARN] <2 muestras tras alineación; se omite beta/PCoA.")
        return

    # Bray–Curtis (en relativas)
    X = tab.values.T.astype(float)  # samples x features
    # Normalizar a relativas por fila si vienen de counts
    row_sums = X.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    X_rel = X / row_sums
    D_bray = _bray_curtis_dm(X_rel)
    _save_dm_and_pcoa(D_bray, list(tab.columns), outdir, "braycurtis",
                      make_png, make_html, meta2, group_col)

    # Jaccard (binario)
    Xbin = (X > 0).astype(int)
    D_jac = _jaccard_dm_bin(Xbin)
    _save_dm_and_pcoa(D_jac, list(tab.columns), outdir, "jaccard",
                      make_png, make_html, meta2, group_col)


# ------------------------------ barras por género -----------------------------

def _extract_genus_from_taxstr(tax: str) -> str:
    # tax string estilo "k__...; p__...; ...; g__Genus; s__Species"
    parts = [t.strip() for t in str(tax).split(";")]
    for p in parts:
        if p.startswith("g__"):
            g = p[3:].strip()
            return g if g else "Unassigned"
    return "Unassigned"


def plot_bars_topn_genus(rel_df: pd.DataFrame,
                         taxonomy_df: pd.DataFrame,
                         outdir: Path,
                         topn: int = 10,
                         make_png: bool = True,
                         make_html: bool = False,
                         meta: Optional[pd.DataFrame] = None,
                         group_col: Optional[str] = None):
    # Necesitamos relativas por feature x muestra
    if rel_df is None or rel_df.empty:
        eprint("[WARN] No hay tabla de relativas para barras por género.")
        return
    # Alinear muestras con metadata (opcional)
    rel_df2, meta2 = _align_with_metadata(rel_df, meta)

    # Mapear feature -> genus
    tax_map = taxonomy_df.set_index("feature_id")["taxonomy"].map(_extract_genus_from_taxstr)
    # Intersectar features
    feats_common = rel_df2.index.intersection(tax_map.index)
    if len(feats_common) == 0:
        eprint("[WARN] taxonomy.tsv no coincide con feature_table; se omiten barras por género.")
        return
    rel_use = rel_df2.loc[feats_common].copy()
    genus = tax_map.loc[feats_common].values

    # Colapsar a género
    rel_use["__genus__"] = genus
    rel_genus = rel_use.groupby("__genus__").sum(numeric_only=True)

    # TopN por media global
    mean_abund = rel_genus.mean(axis=1).sort_values(ascending=False)
    top_gen = list(mean_abund.head(topn).index)
    rel_top = rel_genus.loc[top_gen]
    # "Otros"
    if rel_genus.shape[0] > len(top_gen):
        others = rel_genus.drop(index=top_gen).sum(axis=0)
        rel_top.loc["Otros"] = others

    # Normalizar a 1 por muestra (por si venían de counts)
    col_sums = rel_top.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    rel_top = rel_top / col_sums

    # --- PNG (matplotlib)
    if make_png and HAS_MPL:
        fig_dir = outdir / "plots_static"
        ensure_dir(fig_dir)
        rel_top.T.plot(kind="bar", stacked=True, figsize=(12, 5))
        plt.ylabel("Abundancia relativa")
        plt.xlabel("Muestras")
        plt.title(f"Top {topn} géneros (+ Otros)")
        plt.tight_layout()
        plt.savefig(fig_dir / "bars_topN_genus.png", dpi=150)
        plt.close()

    # --- HTML (plotly)
    if make_html and HAS_PLOTLY:
        html_dir = outdir / "plots_interactive"
        ensure_dir(html_dir)
        df_melt = rel_top.T.reset_index().melt(id_vars="index", var_name="Genus", value_name="Relative Abundance")
        df_melt = df_melt.rename(columns={"index": "Sample"})
        color = "Genus"
        facet_col = None
        if meta2 is not None:
            col = group_col
            if not col or col not in meta2.columns:
                # Heurística: primera categórica con <= 20 niveles
                for c in meta2.columns:
                    if (meta2[c].dtype == "O" or str(meta2[c].dtype).startswith("category")) and meta2[c].nunique() <= 20:
                        col = c
                        break
            if col and col in meta2.columns:
                df_melt["Group"] = meta2[col].reindex(df_melt["Sample"]).values
                facet_col = "Group"
        fig = px.bar(df_melt, x="Sample", y="Relative Abundance", color=color, facet_row=None,
                     facet_col=facet_col, title=f"Top {topn} géneros (+ Otros)", barmode="stack")
        fig.write_html(str(html_dir / "bars_topN_genus.html"), include_plotlyjs="cdn")


# ----------------------------------- CLI --------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Calcula diversidad alfa/beta + PCoA y genera figuras desde salidas de EMU.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--emu-outdir", required=True, type=Path, help="Carpeta con feature_table_*.tsv y taxonomy.tsv")
    p.add_argument("--metadata", type=Path, default=None, help="TSV con filas=muestras (opcional)")
    p.add_argument("--group-col", type=str, default=None, help="Columna de metadata para colorear/ facetar")
    p.add_argument("--plots-static", action="store_true", help="Guardar figuras PNG")
    p.add_argument("--plots-interactive", action="store_true", help="Guardar figuras HTML interactivas")
    p.add_argument("--topn-genus", type=int, default=10, help="Top N géneros para barras apiladas")
    return p


# ---------------------------------- main --------------------------------------

def main() -> int:
    args = build_parser().parse_args()
    outdir = args.emu_outdir
    ensure_dir(outdir)

    # Cargar tablas base
    counts_df, rel_df = load_feature_tables(outdir)

    # Metadata
    meta = _load_metadata(args.metadata)

    # Calcular ALFA
    eprint("[STEP] Alfa diversidad…")
    alpha_df = alpha_from_table(counts_df, rel_df)
    alpha_path = outdir / "alpha_diversity.tsv"
    alpha_df.to_csv(alpha_path, sep="\t")
    eprint(f"[OK] Escrito: {alpha_path}")

    # Calcular BETA + PCoA
    eprint("[STEP] Beta diversidad + PCoA…")
    beta_and_pcoa(counts_df, rel_df, outdir, args.plots_static, args.plots_interactive, meta, args.group_col)
    eprint("[OK] Beta + PCoA listos.")

    # Barras apiladas por género (requiere taxonomy + relativas)
    tax_df = load_taxonomy(outdir)
    if tax_df is not None:
        # Asegurar relativas (si no había, normalizar desde counts)
        if rel_df is None and counts_df is not None:
            tab = _drop_all_zero_features(_drop_zero_samples(counts_df))
            if tab.shape[1] > 0:
                rel_df2 = tab / tab.sum(axis=0).replace(0, 1.0)
            else:
                rel_df2 = None
        else:
            rel_df2 = rel_df
        if rel_df2 is not None and rel_df2.shape[1] > 0:
            eprint("[STEP] Barras apiladas por género…")
            plot_bars_topn_genus(rel_df2, tax_df, outdir, topn=args.topn_genus,
                                 make_png=args.plots_static, make_html=args.plots_interactive,
                                 meta=meta, group_col=args.group_col)
            eprint("[OK] Barras por género listas.")
        else:
            eprint("[WARN] No hay relativas utilizables; se omiten barras por género.")

    eprint("Listo. Métricas y figuras generadas.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
