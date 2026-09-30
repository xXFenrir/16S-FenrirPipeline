#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
metricas.py
-----------
Calcula métricas de diversidad alfa y beta (Bray–Curtis, Jaccard y PCoA)
a partir de salidas de EMU. Agrega de forma robusta tablas de conteos y
abundancias relativas, y cae automáticamente a relativas si los conteos
están vacíos.

Entradas esperadas dentro de --emu-outdir:
- (Opcional) feature_table_counts.tsv
- (Opcional) feature_table_relabund.tsv
- (Opcional) taxonomy.tsv
- O bien, TSVs por muestra con '*counts*.tsv' y '*rel*abundance*.tsv'

Salidas (en --outdir, por defecto = --emu-outdir):
- feature_table_counts.tsv (agregada o existente)
- feature_table_relabund.tsv (agregada o existente)
- alpha_diversity.tsv
- braycurtis_dm.tsv, jaccard_dm.tsv
- pcoa_braycurtis_coords.tsv, pcoa_braycurtis_variance.tsv
- pcoa_jaccard_coords.tsv,   pcoa_jaccard_variance.tsv
- plots_py/  (PNG)
- plots_interactive/ (HTML) si se pide --plots-interactive
"""

import argparse
import os
import sys
import re
import glob
import math
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd

# For static plots
import matplotlib
# Si corres en servidor sin display, forzamos backend no interactivo
if os.environ.get("DISPLAY", "") == "":
    matplotlib.use("Agg")

import matplotlib.pyplot as plt

# Plotly (solo si piden interactivo)
try:
    import plotly.express as px
except Exception:
    px = None


# ----------------------- Utilidades de E/S -----------------------

def _mkdir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def _is_tsv(path: str) -> bool:
    return path.lower().endswith(".tsv")


def _find_sample_id_from_path(path: str) -> str:
    """Extrae ID de muestra (SRR/ERR/DRR...) del path o usa el stem del archivo."""
    m = re.search(r'(SRR\d+|ERR\d+|DRR\d+)', path)
    if m:
        return m.group(1)
    return Path(path).stem


def _read_tsv(p: str) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t")


def _list_files_recursive(base: str, pattern: str) -> List[str]:
    return glob.glob(os.path.join(base, "**", pattern), recursive=True)


def _save_tsv(df: pd.DataFrame, out_path: str, index: bool = True) -> None:
    df.to_csv(out_path, sep="\t", index=index)


# ----------------------- Agregación robusta -----------------------

COUNT_COL_CANDIDATES = {"count", "counts", "read_count", "reads", "n"}
REL_COL_CANDIDATES = {"relative_abundance", "rel_abundance", "relabundance", "abundance", "rel", "proportion"}

def _pick_col(df: pd.DataFrame, candidates: set, human_label: str, path: str) -> str:
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in low:
            return low[cand]
    # último recurso: si hay una sola numérica plausible
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if len(numeric_cols) == 1:
        return numeric_cols[0]
    raise ValueError(f"No encontré columna de {human_label} en {path}. "
                     f"Disponibles: {list(df.columns)}")


def _normalize_rank_name(df: pd.DataFrame, rank: str, path: str) -> str:
    """Devuelve el nombre de columna de taxonomía 'rank' con case-insensitive."""
    if rank in df.columns:
        return rank
    low = {c.lower(): c for c in df.columns}
    if rank.lower() in low:
        return low[rank.lower()]
    raise ValueError(f"No encontré la columna taxonómica '{rank}' en {path}.")


def _aggregate_counts_from_per_sample(emu_outdir: str, rank: str) -> Optional[pd.DataFrame]:
    """Busca TSVs por muestra que contengan counts y agrega a tabla wide (filas=taxa, cols=muestras)."""
    files = _list_files_recursive(emu_outdir, "*count*.tsv")
    if not files:
        return None

    series = []
    for p in files:
        try:
            df = _read_tsv(p)
            rank_col = _normalize_rank_name(df, rank, p)
            count_col = _pick_col(df, COUNT_COL_CANDIDATES, "conteos", p)
            tmp = df[[rank_col, count_col]].dropna()
            tmp[count_col] = pd.to_numeric(tmp[count_col], errors="coerce").fillna(0)
            tmp = tmp.groupby(rank_col, as_index=False)[count_col].sum()
            s = tmp.set_index(rank_col)[count_col]
            s.name = _find_sample_id_from_path(p)
            series.append(s)
        except Exception as e:
            print(f"[WARN] Saltando counts de {p}: {e}", file=sys.stderr)

    if not series:
        return None

    counts = pd.concat(series, axis=1).fillna(0)
    # Asegurar tipo entero si tiene sentido
    try:
        counts = counts.round().astype(int)
    except Exception:
        counts = counts.astype(float)
    return counts


def _aggregate_rel_from_per_sample(emu_outdir: str, rank: str) -> Optional[pd.DataFrame]:
    """Busca TSVs por muestra con abundancia relativa y agrega wide."""
    files = _list_files_recursive(emu_outdir, "*rel*abundance*.tsv")
    if not files:
        return None

    series = []
    for p in files:
        try:
            df = _read_tsv(p)
            rank_col = _normalize_rank_name(df, rank, p)
            rel_col = _pick_col(df, REL_COL_CANDIDATES, "abundancias relativas", p)
            tmp = df[[rank_col, rel_col]].dropna()
            tmp[rel_col] = pd.to_numeric(tmp[rel_col], errors="coerce").fillna(0.0)
            tmp = tmp.groupby(rank_col, as_index=False)[rel_col].sum()
            s = tmp.set_index(rank_col)[rel_col]
            s.name = _find_sample_id_from_path(p)
            series.append(s)
        except Exception as e:
            print(f"[WARN] Saltando relativas de {p}: {e}", file=sys.stderr)

    if not series:
        return None

    rel = pd.concat(series, axis=1).fillna(0.0)
    # Normaliza cada muestra a 1 (por seguridad)
    colsum = rel.sum(axis=0).replace(0, np.nan)
    rel = rel.div(colsum, axis=1).fillna(0.0)
    return rel


def _ensure_aggregated_tables(emu_outdir: str, rank: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Devuelve (counts, rel). Si no existen agregadas, intenta construirlas desde TSVs por muestra.
    """
    counts_path = os.path.join(emu_outdir, "feature_table_counts.tsv")
    rel_path = os.path.join(emu_outdir, "feature_table_relabund.tsv")

    counts = pd.read_csv(counts_path, sep="\t", index_col=0) if os.path.exists(counts_path) else None
    rel = pd.read_csv(rel_path, sep="\t", index_col=0) if os.path.exists(rel_path) else None

    if counts is None:
        counts = _aggregate_counts_from_per_sample(emu_outdir, rank)
        if counts is not None:
            _save_tsv(counts, counts_path)

    if rel is None:
        rel = _aggregate_rel_from_per_sample(emu_outdir, rank)
        if rel is not None:
            _save_tsv(rel, rel_path)

    if counts is None and rel is None:
        raise RuntimeError("No encontré ni tablas agregadas ni TSVs por muestra para counts/relabund.")

    # Asegurar orientación (filas=taxa, columnas=muestras)
    if counts is not None and counts.shape[0] < counts.shape[1] and counts.columns.str.match(r'^(SRR|ERR|DRR|Sample|S\d+)', na=False).any():
        pass  # parece correcto (cols son muestras)
    if rel is not None and rel.shape[0] < rel.shape[1] and rel.columns.str.match(r'^(SRR|ERR|DRR|Sample|S\d+)', na=False).any():
        pass

    return counts, rel


# ----------------------- Métricas de diversidad -----------------------

def _as_relative_from_counts(counts: pd.DataFrame) -> pd.DataFrame:
    counts = counts.clip(lower=0)
    sums = counts.sum(axis=0).replace(0, np.nan)
    rel = counts.div(sums, axis=1).fillna(0.0)
    return rel


def alpha_diversity(counts: Optional[pd.DataFrame], rel: Optional[pd.DataFrame]) -> pd.DataFrame:
    """
    Calcula Observed, Shannon, Simpson; Chao1 si hay conteos > 0.
    Si counts es None o todo 0, usa relativas para Shannon/Simpson/Observed (Chao1=NaN).
    """
    use_rel = False
    if counts is None or counts.values.sum() == 0:
        if rel is None:
            raise RuntimeError("No hay counts ni relativas para calcular alfa.")
        use_rel = True
        X = rel.copy()
        print("[WARN] Tabla de conteos vacía; alfa con relativas (Chao1=NaN).", file=sys.stderr)
    else:
        X = _as_relative_from_counts(counts)

    samples = X.columns
    obs = (X > 0).sum(axis=0)
    shannon = -(X.replace(0, np.nan) * np.log(X.replace(0, np.nan))).sum(axis=0).fillna(0.0)
    simpson = 1.0 - (X**2).sum(axis=0)

    chao1_vals = []
    if not use_rel:
        for s in samples:
            col = counts[s]
            if col.sum() == 0:
                chao1_vals.append(np.nan)
                continue
            f1 = (col == 1).sum()
            f2 = (col == 2).sum()
            sobs = (col > 0).sum()
            if f2 == 0:
                chao = np.nan  # o sobs + f1*(f1-1)/2 para corrección
            else:
                chao = sobs + (f1*f1)/(2.0*f2)
            chao1_vals.append(chao)
    else:
        chao1_vals = [np.nan] * len(samples)

    df = pd.DataFrame({
        "sample": samples,
        "observed": obs.values.astype(float),
        "shannon": shannon.values.astype(float),
        "simpson": simpson.values.astype(float),
        "chao1": chao1_vals
    })
    return df


def bray_curtis_and_jaccard(counts: Optional[pd.DataFrame], rel: Optional[pd.DataFrame]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calcula Bray–Curtis y Jaccard. Si counts está vacío, usa relativas.
    Bray–Curtis sobre relativas: 0.5 * L1.
    Jaccard sobre presencia/ausencia (a partir de la misma X usada).
    """
    if counts is not None and counts.values.sum() > 0:
        X = _as_relative_from_counts(counts)
    elif rel is not None:
        print("[WARN] Tabla de conteos vacía; beta y PCoA desde relativas.", file=sys.stderr)
        X = rel.copy()
    else:
        raise RuntimeError("No hay matriz válida para beta diversidad.")

    X = X.clip(lower=0.0)
    # normaliza columnas a 1 (por seguridad)
    X = X.div(X.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)

    samples = list(X.columns)
    n = len(samples)

    # Bray–Curtis
    bray = pd.DataFrame(0.0, index=samples, columns=samples)
    for i in range(n):
        xi = X.iloc[:, i].values
        for j in range(i+1, n):
            xj = X.iloc[:, j].values
            d = 0.5 * np.abs(xi - xj).sum()
            bray.iat[i, j] = bray.iat[j, i] = d

    # Jaccard (presencia/ausencia)
    PA = (X > 0).astype(int)
    jac = pd.DataFrame(0.0, index=samples, columns=samples)
    for i in range(n):
        ai = PA.iloc[:, i].values
        for j in range(i+1, n):
            aj = PA.iloc[:, j].values
            inter = int((ai & aj).sum())
            union = int((ai | aj).sum())
            d = 1.0 if union == 0 else 1.0 - (inter / union)
            jac.iat[i, j] = jac.iat[j, i] = d

    return bray, jac


def pcoa(distance_matrix: pd.DataFrame, k: int = 3) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """PCoA clásica (eigendecomposition del double-centered D^2)."""
    labels = distance_matrix.index.tolist()
    D = distance_matrix.values.astype(float)
    D2 = D ** 2
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1]
    w = w[idx]
    V = V[:, idx]
    pos = w > 0
    w_pos = w[pos][:k]
    V_pos = V[:, pos][:, :k]
    coords = V_pos * np.sqrt(w_pos)
    var = (w_pos / w[pos].sum()) if pos.any() else np.array([])
    coords_df = pd.DataFrame(coords, index=labels, columns=[f"PC{i+1}" for i in range(coords.shape[1])])
    var_df = pd.DataFrame({"axis": [f"PC{i+1}" for i in range(len(var))],
                           "explained_variance": var})
    return coords_df, var_df


# ----------------------- Plotting -----------------------

def _load_metadata(meta_path: Optional[str]) -> Optional[pd.DataFrame]:
    if not meta_path:
        return None
    df = pd.read_csv(meta_path, sep="\t")
    # se espera columna 'sample'; 'group' es opcional
    if "sample" not in df.columns:
        low = {c.lower(): c for c in df.columns}
        if "sample" in low:
            df.rename(columns={low["sample"]: "sample"}, inplace=True)
        else:
            raise ValueError("Metadata debe tener una columna 'sample'.")
    return df


def plot_alpha(df_alpha: pd.DataFrame, outdir: str, metadata: Optional[pd.DataFrame]) -> None:
    _mkdir(outdir)
    for metric in ["shannon", "simpson", "observed"]:
        plt.figure(figsize=(8, 4))
        if metadata is not None and "group" in metadata.columns:
            merged = df_alpha.merge(metadata[["sample", "group"]], on="sample", how="left")
            groups = merged["group"].fillna("NA").unique()
            for g in groups:
                sub = merged[merged["group"].fillna("NA") == g]
                plt.scatter(np.arange(len(sub)), sub[metric], label=str(g), alpha=0.8)
            plt.legend(title="group", fontsize=9)
            plt.xticks([])
        else:
            plt.bar(df_alpha["sample"], df_alpha[metric])
            plt.xticks(rotation=90, fontsize=7)
        plt.title(f"Alpha — {metric}")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"alpha_{metric}.png"), dpi=150)
        plt.close()


def plot_pcoa(coords: pd.DataFrame, var: pd.DataFrame, outdir: str, name: str,
              metadata: Optional[pd.DataFrame], interactive: bool = False) -> None:
    _mkdir(outdir)
    xlab = f"PC1 ({(var['explained_variance'].iloc[0]*100):.1f}%)" if len(var) > 0 else "PC1"
    ylab = f"PC2 ({(var['explained_variance'].iloc[1]*100):.1f}%)" if len(var) > 1 else "PC2"

    # Static
    plt.figure(figsize=(6, 5))
    if metadata is not None and "group" in metadata.columns:
        merged = coords.merge(metadata[["sample", "group"]], left_index=True, right_on="sample", how="left")
        for g, sub in merged.groupby(merged["group"].fillna("NA")):
            plt.scatter(sub["PC1"], sub["PC2"], label=str(g), alpha=0.9)
        plt.legend(title="group", fontsize=9)
    else:
        plt.scatter(coords["PC1"], coords["PC2"], alpha=0.9)
        for s, row in coords.iterrows():
            plt.text(row["PC1"], row["PC2"], s, fontsize=7)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(f"PCoA — {name}")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"pcoa_{name}.png"), dpi=150)
    plt.close()

    # Interactive
    if interactive and px is not None:
        html_dir = os.path.join(os.path.dirname(outdir), "plots_interactive")
        _mkdir(html_dir)
        plot_df = coords.copy()
        plot_df["sample"] = plot_df.index
        color = None
        if metadata is not None and "group" in metadata.columns:
            plot_df = plot_df.merge(metadata[["sample", "group"]], on="sample", how="left")
            color = "group"
        fig = px.scatter(plot_df, x="PC1", y="PC2", hover_name="sample", color=color,
                         title=f"PCoA — {name} ({xlab} vs {ylab})")
        fig.write_html(os.path.join(html_dir, f"pcoa_{name}.html"), include_plotlyjs="cdn")


def plot_distance_heatmap(D: pd.DataFrame, outdir: str, name: str) -> None:
    _mkdir(outdir)
    plt.figure(figsize=(6.5, 5.5))
    im = plt.imshow(D.values, aspect="auto", interpolation="nearest")
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.xticks(ticks=np.arange(D.shape[1]), labels=D.columns, rotation=90, fontsize=7)
    plt.yticks(ticks=np.arange(D.shape[0]), labels=D.index, fontsize=7)
    plt.title(f"Distance heatmap — {name}")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"heatmap_{name}.png"), dpi=150)
    plt.close()


def plot_composition_stacked(rel: pd.DataFrame, outdir: str, topn: int = 12) -> None:
    _mkdir(outdir)
    # top-N por abundancia media
    mean_ab = rel.mean(axis=1).sort_values(ascending=False)
    top_taxa = mean_ab.index[:topn]
    rel_top = rel.loc[top_taxa].copy()
    other = rel.drop(index=top_taxa, errors="ignore").sum(axis=0)
    rel_top.loc["Other"] = other

    # Re-normaliza a 1 por muestra
    rel_top = rel_top.div(rel_top.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)

    # Plot
    samples = list(rel_top.columns)
    bottoms = np.zeros(len(samples))
    plt.figure(figsize=(max(8, len(samples)*0.35), 5))
    for taxon in rel_top.index:
        vals = rel_top.loc[taxon, samples].values
        plt.bar(samples, vals, bottom=bottoms, label=taxon, width=0.8)
        bottoms += vals
    plt.xticks(rotation=90, fontsize=7)
    plt.ylabel("Relative abundance")
    plt.title(f"Composition (top {topn})")
    plt.legend(fontsize=7, ncol=2, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"stacked_bar_top{topn}.png"), dpi=150)
    plt.close()


# ----------------------- Main -----------------------

def parse_args():
    p = argparse.ArgumentParser(description="Diversidad alfa/beta y PCoA desde salidas EMU.")
    p.add_argument("--emu-outdir", required=True, help="Carpeta con resultados de EMU (por muestra o agregados).")
    p.add_argument("--outdir", default=None, help="Carpeta de salida. Por defecto = --emu-outdir.")
    p.add_argument("--rank", default="species", help="Rango taxonómico para indexar/agregar (default: species).")
    p.add_argument("--metadata", default=None, help="TSV con columna 'sample' y opcional 'group'.")
    p.add_argument("--plots-static", action="store_true", help="Guardar PNG de figuras.")
    p.add_argument("--plots-interactive", action="store_true", help="Guardar HTML interactivo (Plotly).")
    p.add_argument("--topn-genus", type=int, default=12, help="Top-N para barras apiladas (usa índices tal como estén).")
    return p.parse_args()


def main():
    args = parse_args()
    emu_outdir = os.path.abspath(args.emu_outdir)
    outdir = os.path.abspath(args.outdir or emu_outdir)
    _mkdir(outdir)

    print(f"[INFO] EMU outdir: {emu_outdir}")
    print(f"[INFO] Outdir    : {outdir}")
    print(f"[INFO] Rank      : {args.rank}")

    # 1) Asegurar tablas agregadas
    counts, rel = _ensure_aggregated_tables(emu_outdir, args.rank)

    if counts is None:
        print("[WARN] No se obtuvo tabla de conteos agregada.", file=sys.stderr)
    else:
        # persistimos (por si la acabamos de construir)
        _save_tsv(counts, os.path.join(outdir, "feature_table_counts.tsv"))

    if rel is None and counts is not None and counts.values.sum() > 0:
        rel = _as_relative_from_counts(counts)
        _save_tsv(rel, os.path.join(outdir, "feature_table_relabund.tsv"))
    elif rel is not None:
        _save_tsv(rel, os.path.join(outdir, "feature_table_relabund.tsv"))

    # 2) Alfa diversidad
    df_alpha = alpha_diversity(counts, rel)
    _save_tsv(df_alpha.set_index("sample"), os.path.join(outdir, "alpha_diversity.tsv"))

    # 3) Beta diversidad
    bray, jac = bray_curtis_and_jaccard(counts, rel)
    _save_tsv(bray, os.path.join(outdir, "braycurtis_dm.tsv"))
    _save_tsv(jac,  os.path.join(outdir, "jaccard_dm.tsv"))

    # 4) PCoA
    pcoa_b_coords, pcoa_b_var = pcoa(bray)
    pcoa_j_coords, pcoa_j_var = pcoa(jac)
    _save_tsv(pcoa_b_coords, os.path.join(outdir, "pcoa_braycurtis_coords.tsv"))
    pcoa_b_var.to_csv(os.path.join(outdir, "pcoa_braycurtis_variance.tsv"), sep="\t", index=False)
    _save_tsv(pcoa_j_coords, os.path.join(outdir, "pcoa_jaccard_coords.tsv"))
    pcoa_j_var.to_csv(os.path.join(outdir, "pcoa_jaccard_variance.tsv"), sep="\t", index=False)

    # 5) Plots
    meta = _load_metadata(args.metadata) if args.metadata else None
    if args.plots_static:
        plots_dir = os.path.join(outdir, "plots_py")
        _mkdir(plots_dir)
        plot_alpha(df_alpha, plots_dir, meta)
        plot_distance_heatmap(bray, plots_dir, "braycurtis")
        plot_distance_heatmap(jac, plots_dir, "jaccard")
        plot_composition_stacked(rel if rel is not None else _as_relative_from_counts(counts),
                                 plots_dir, topn=args.topn_genus)
        # PCoA
        pcoa_dir = plots_dir
        plot_pcoa(pcoa_b_coords, pcoa_b_var, pcoa_dir, "braycurtis", meta, interactive=False)
        plot_pcoa(pcoa_j_coords, pcoa_j_var, pcoa_dir, "jaccard", meta, interactive=False)

    if args.plots_interactive:
        if px is None:
            print("[WARN] Plotly no disponible; omitiendo HTML interactivos.", file=sys.stderr)
        else:
            html_dir = os.path.join(outdir, "plots_interactive")
            _mkdir(html_dir)
            # Solo PCoA interactivo (lo más útil)
            plot_pcoa(pcoa_b_coords, pcoa_b_var, os.path.join(outdir, "plots_py"),
                      "braycurtis", meta, interactive=True)
            plot_pcoa(pcoa_j_coords, pcoa_j_var, os.path.join(outdir, "plots_py"),
                      "jaccard", meta, interactive=True)

    print("[OK] Métricas y figuras listas en:", outdir)


if __name__ == "__main__":
    main()

