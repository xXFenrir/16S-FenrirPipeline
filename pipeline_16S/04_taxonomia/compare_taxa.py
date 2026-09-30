#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compare_taxa.py — Comparación de composiciones taxonómicas (tus resultados vs. artículo)

Permite cargar tablas agregadas o per-muestra (relativas) y comparar en un rango
taxonómico específico (p. ej., genus). Genera tablas de comparación y figuras.

Requisitos (conda-forge):
  conda install -c conda-forge pandas numpy plotly kaleido

Ejemplos de uso:

# 1) Mis resultados ya agregados (relativas + taxonomy) VS artículo (counts+taxonomy)
python compare_taxa.py \
  --my-results-dir /ruta/EMU_out \
  --article-counts /ruta/art/emu_16S_counts_2020.csv \
  --article-taxonomy /ruta/art/emu_16S_taxonomy_2020.csv \
  --rank genus

# 2) Mis per-muestra relativas (glob) VS artículo per-muestra relativas (glob)
python compare_taxa.py \
  --my-relabund-glob "/ruta/EMU_out/samples/*/*_rel-abundance.tsv" \
  --article-relabund-glob "/ruta/art/rel_abund/barcode*.t_rel-abundance*.tsv" \
  --rank genus \
  --outdir /ruta/compare_out

Salidas (en --outdir, por defecto <my-results-dir>/compare_taxa):
  - mine_{rank}_rel.tsv            (tabla rank x muestras, relativas)
  - article_{rank}_rel.tsv
  - compare_{rank}_means.tsv       (medias por taxón y deltas)
  - compare_{rank}_correlations.tsv (Pearson/Spearman sobre taxones comunes)
  - compare_{rank}_scatter.html/png (dispersión medias: mío vs artículo)
  - compare_{rank}_top{N}_bars.html/png (barras lado a lado top taxa)
"""

from __future__ import annotations
import argparse
import sys, os, re, glob
from pathlib import Path
from typing import Optional, Tuple, Dict, List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# export estático
_HAS_KALEIDO = False
try:
    import kaleido  # noqa: F401
    _HAS_KALEIDO = True
except Exception:
    pass

TAX_RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


# -------------------------- utilidades generales --------------------------

def eprint(*a, **k): print(*a, file=sys.stderr, **k)

def _ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def _read_table_any(path: Path, index_col: Optional[int|str] = None) -> pd.DataFrame:
    """Lee CSV o TSV (detecta sep automáticamente)."""
    return pd.read_csv(path, sep=None, engine="python", index_col=index_col)

def _coerce_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.fillna(0.0)

def _is_sample_like_index(idx: pd.Index) -> bool:
    s = pd.Series(idx.astype(str))
    return s.str.match(r"^(SRR|ERR|DRR)\d+").any() or s.str.contains(r"^barcode", case=False, na=False).any()

def _maybe_features_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Asegura orientación: filas=features, columnas=muestras."""
    if _is_sample_like_index(df.index):
        return df.T
    return df

def _save_fig(fig: go.Figure, html_path: Path, png_path: Optional[Path] = None,
              width=1100, height=750, title: Optional[str] = None):
    if title:
        fig.update_layout(title=title)
    fig.update_layout(template="plotly_white", width=width, height=height)
    fig.write_html(str(html_path))
    if png_path and _HAS_KALEIDO:
        fig.write_image(str(png_path))
    elif png_path:
        eprint(f"[INFO] Instala 'kaleido' para PNG/SVG/PDF. Guardé HTML: {html_path}")

# ---------------------- parsing taxonómico (cadenas QIIME) ----------------------

_QIIME_PATTERNS = {
    "superkingdom": re.compile(r"k__([^;]+)"),
    "phylum": re.compile(r"p__([^;]+)"),
    "class": re.compile(r"c__([^;]+)"),
    "order": re.compile(r"o__([^;]+)"),
    "family": re.compile(r"f__([^;]+)"),
    "genus": re.compile(r"g__([^;]+)"),
    "species": re.compile(r"s__([^;]+)"),
}

def _from_qiime(s: str, rank: str) -> str:
    if not isinstance(s, str):
        return "Unassigned"
    m = _QIIME_PATTERNS[rank].search(s)
    if not m:
        return "Unassigned"
    val = m.group(1).strip()
    if val in ["", "uncultured", "metagenome"]:
        return "Unassigned"
    return val

def taxonomy_to_rank_map(tax_df: pd.DataFrame, rank: str) -> pd.Series:
    """Devuelve Serie index=feature_id, values=nombre del rank."""
    df = tax_df.copy()
    # estandariza nombre de id
    if "feature_id" in df.columns:
        df = df.set_index("feature_id")
    # detecta si ya trae columna del rank
    if rank in df.columns:
        s = df[rank].astype(str)
    elif "taxonomy" in df.columns:
        s = df["taxonomy"].map(lambda x: _from_qiime(str(x), rank))
    else:
        # intenta reconstruir cadena taxonomy si tiene columnas de ranks
        cols_lower = {c.lower(): c for c in df.columns}
        if all(r in cols_lower for r in TAX_RANKS):
            joined = (
                "k__" + df[cols_lower["superkingdom"]].astype(str) + "; " +
                "p__" + df[cols_lower["phylum"]].astype(str) + "; " +
                "c__" + df[cols_lower["class"]].astype(str) + "; " +
                "o__" + df[cols_lower["order"]].astype(str) + "; " +
                "f__" + df[cols_lower["family"]].astype(str) + "; " +
                "g__" + df[cols_lower["genus"]].astype(str) + "; " +
                "s__" + df[cols_lower["species"]].astype(str)
            )
            s = joined.map(lambda x: _from_qiime(str(x), rank))
        else:
            raise RuntimeError("taxonomy no tiene ni 'taxonomy' (QIIME) ni columna del rank solicitado.")
    s.name = "rank_name"
    return s

# ------------------------- detectar columna de relativas -------------------------

def detect_rel_col(df: pd.DataFrame) -> Optional[str]:
    # candidata que sume ~1
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]):
            s = pd.to_numeric(df[c], errors="coerce").fillna(0).sum()
            if 0.9 <= s <= 1.1:
                return c
    # por nombre
    cand = [c for c in df.columns if re.search(r"abund", c, flags=re.I)]
    for c in cand:
        s = pd.to_numeric(df[c], errors="coerce").fillna(0).sum()
        if 0.9 <= s <= 1.1:
            return c
    return None

# --------------------------- construcción de tablas rank ---------------------------

def counts_to_relative(counts: pd.DataFrame) -> pd.DataFrame:
    C = _coerce_numeric_df(counts)
    if "feature_id" in C.columns:
        C = C.set_index("feature_id")
    C = _maybe_features_rows(C)
    colsum = C.sum(axis=0)
    colsum[colsum == 0] = 1.0
    return C / colsum

def aggregate_by_rank_rel(rel: pd.DataFrame, taxmap: pd.Series) -> pd.DataFrame:
    """Entrada: rel (features x samples), taxmap (index=features -> rank_name)."""
    R = _maybe_features_rows(rel.copy()).fillna(0.0)
    # alinear
    common = R.index.intersection(taxmap.index)
    if common.empty:
        # quizá rel trae features en columnas
        R = R.T
        common = R.index.intersection(taxmap.index)
    if common.empty:
        raise RuntimeError("No coinciden los feature_id entre tabla y taxonomía.")
    R = R.loc[common]
    taxmap = taxmap.loc[common]
    G = R.groupby(taxmap).sum()
    # normaliza por muestra a proporciones
    colsum = G.sum(axis=0)
    colsum[colsum == 0] = 1.0
    return G / colsum

def load_my(
    my_results_dir: Optional[Path],
    rank: str,
    my_counts: Optional[Path],
    my_taxonomy: Optional[Path],
    my_relabund_glob: Optional[str],
) -> pd.DataFrame:
    """
    Devuelve tabla rank x muestras (relativas) para tus datos.
    Orden de preferencia:
      1) counts + taxonomy  -> relativas + agregación por rank
      2) results_dir con feature_table_relabund.tsv + taxonomy.tsv
      3) glob de *_rel-abundance.tsv (requiere tener rank o taxonomy en cada TSV,
         o proporcionar --my-taxonomy para mapear).
    """
    # 1) counts + taxonomy
    if my_counts and my_taxonomy:
        C = _read_table_any(Path(my_counts), index_col=0)
        T = _read_table_any(Path(my_taxonomy))
        taxmap = taxonomy_to_rank_map(T, rank=rank)
        R = counts_to_relative(C)
        return aggregate_by_rank_rel(R, taxmap)

    # 2) results_dir
    if my_results_dir:
        rel_p = my_results_dir / "feature_table_relabund.tsv"
        tax_p = my_results_dir / "taxonomy.tsv"
        if rel_p.exists() and tax_p.exists():
            Rel = _read_table_any(rel_p, index_col=0)
            T = _read_table_any(tax_p)
            taxmap = taxonomy_to_rank_map(T, rank=rank)
            return aggregate_by_rank_rel(Rel, taxmap)

    # 3) glob per-muestra
    if my_relabund_glob:
        paths = glob.glob(my_relabund_glob)
        if not paths:
            raise RuntimeError(f"No encontré TSV con el patrón: {my_relabund_glob}")
        tabs = []
        need_tax = False
        for p in paths:
            df = _read_table_any(Path(p))
            # detectar rel col
            relc = detect_rel_col(df)
            if relc is None:
                raise RuntimeError(f"No detecté columna de relativas en: {p}")
            # rank directo o desde taxonomy
            if rank in df.columns:
                df_rank = df[[rank, relc]].groupby(rank).sum()
            elif "taxonomy" in df.columns:
                df["_rank_"] = df["taxonomy"].map(lambda s: _from_qiime(str(s), rank))
                df_rank = df[["_rank_", relc]].groupby("_rank_").sum()
                df_rank.index.name = rank
            else:
                need_tax = True
                tabs = []  # invalida acumulado
                break
            sname = _sample_name_from_path(Path(p))
            tabs.append(df_rank.rename(columns={relc: sname}))
        if need_tax:
            if not my_taxonomy:
                raise RuntimeError("Tus TSV no traen columna rank ni 'taxonomy'. Pasa --my-taxonomy para mapear.")
            # construir rels por feature -> luego agregar por rank
            # asumimos que TSV por muestra traen filas=features y una columna de relativas
            rels = []
            for p in paths:
                df = _read_table_any(Path(p))
                relc = detect_rel_col(df)
                if relc is None:
                    continue
                if "feature_id" in df.columns:
                    df = df.set_index("feature_id")
                elif "id" in df.columns:
                    df = df.set_index("id")
                else:
                    # si no trae id explícito, lo intentamos con una columna que parezca id
                    raise RuntimeError(f"{p} no tiene 'feature_id' ni 'id' para alinear con taxonomy.")
                sname = _sample_name_from_path(Path(p))
                rels.append(df[[relc]].rename(columns={relc: sname}))
            if not rels:
                raise RuntimeError("No pude construir relativas por feature desde tus TSV.")
            Rel = pd.concat(rels, axis=1).fillna(0.0)
            T = _read_table_any(Path(my_taxonomy))
            taxmap = taxonomy_to_rank_map(T, rank=rank)
            return aggregate_by_rank_rel(Rel, taxmap)
        if not tabs:
            raise RuntimeError("No pude extraer información de relativas por muestra (mis datos).")
        G = pd.concat(tabs, axis=1).fillna(0.0)
        # normaliza columnas a 1 (por seguridad)
        colsum = G.sum(axis=0); colsum[colsum == 0] = 1.0
        return G / colsum

    raise RuntimeError(
        "No encontré tablas agregadas ni pude construirlas para 'mis datos'. "
        "Prueba alguna opción: "
        "--my-results-dir (relabund+taxonomy), "
        "--my-counts + --my-taxonomy, "
        "o --my-relabund-glob."
    )

def load_article(
    rank: str,
    article_counts: Optional[Path],
    article_taxonomy: Optional[Path],
    article_relabund_glob: Optional[str],
) -> pd.DataFrame:
    """Devuelve rank x muestras (relativas) para el artículo."""
    # A) counts + taxonomy
    if article_counts and article_taxonomy:
        C = _read_table_any(Path(article_counts), index_col=0)
        T = _read_table_any(Path(article_taxonomy))
        taxmap = taxonomy_to_rank_map(T, rank=rank)
        R = counts_to_relative(C)
        return aggregate_by_rank_rel(R, taxmap)

    # B) glob per-muestra de relativas
    if article_relabund_glob:
        paths = glob.glob(article_relabund_glob)
        if not paths:
            raise RuntimeError(f"No encontré TSV del artículo con el patrón: {article_relabund_glob}")
        tabs = []
        need_tax = False
        for p in paths:
            df = _read_table_any(Path(p))
            relc = detect_rel_col(df)
            if relc is None:
                raise RuntimeError(f"Artículo: no detecté columna de relativas en: {p}")
            if rank in df.columns:
                df_rank = df[[rank, relc]].groupby(rank).sum()
            elif "taxonomy" in df.columns:
                df["_rank_"] = df["taxonomy"].map(lambda s: _from_qiime(str(s), rank))
                df_rank = df[["_rank_", relc]].groupby("_rank_").sum()
                df_rank.index.name = rank
            else:
                need_tax = True
                tabs = []
                break
            sname = _sample_name_from_path(Path(p))
            tabs.append(df_rank.rename(columns={relc: sname}))
        if need_tax:
            if not article_taxonomy:
                raise RuntimeError("Los TSV del artículo no traen columna rank ni 'taxonomy'. Pasa --article-taxonomy.")
            rels = []
            for p in paths:
                df = _read_table_any(Path(p))
                relc = detect_rel_col(df)
                if relc is None:
                    continue
                if "feature_id" in df.columns:
                    df = df.set_index("feature_id")
                elif "id" in df.columns:
                    df = df.set_index("id")
                else:
                    raise RuntimeError(f"{p} no tiene 'feature_id' ni 'id' para alinear con taxonomy (artículo).")
                sname = _sample_name_from_path(Path(p))
                rels.append(df[[relc]].rename(columns={relc: sname}))
            if not rels:
                raise RuntimeError("No pude construir relativas por feature desde los TSV del artículo.")
            Rel = pd.concat(rels, axis=1).fillna(0.0)
            T = _read_table_any(Path(article_taxonomy))
            taxmap = taxonomy_to_rank_map(T, rank=rank)
            return aggregate_by_rank_rel(Rel, taxmap)
        if not tabs:
            raise RuntimeError("No pude extraer relativas por muestra (artículo).")
        G = pd.concat(tabs, axis=1).fillna(0.0)
        colsum = G.sum(axis=0); colsum[colsum == 0] = 1.0
        return G / colsum

    raise RuntimeError(
        "Debes pasar --article-relabund-glob o --article-counts + --article-taxonomy para el artículo."
    )

def _sample_name_from_path(p: Path) -> str:
    """Nombre de muestra estable a partir de archivo rel-abundance TSV."""
    # intenta capturar 'barcodeNN'
    m = re.search(r"(barcode\d+)", p.name, flags=re.I)
    if m:
        return m.group(1)
    # si viene como SRR/ERR/DRR
    m = re.search(r"(SRR|ERR|DRR)\d{5,}", p.name)
    if m:
        return m.group(0)
    # si el stem tiene sufijo tipo '.t_rel-abundance'
    stem = p.stem
    stem = re.sub(r"\.t_?rel-?abundance.*$", "", stem, flags=re.I)
    return stem

# ------------------------------- comparación / plots -------------------------------

def compare_means(mine: pd.DataFrame, art: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Retorna (df_means, df_corrs). df_means: mine_mean, art_mean, delta, abs_delta."""
    mine_mean = mine.mean(axis=1) if mine.shape[1] > 0 else pd.Series(dtype=float)
    art_mean  = art.mean(axis=1) if art.shape[1] > 0 else pd.Series(dtype=float)
    # alinear por taxón (intersección)
    common_taxa = mine_mean.index.intersection(art_mean.index)
    if common_taxa.empty:
        raise RuntimeError("No hay taxones en común para comparar promedios.")
    mm = pd.DataFrame({
        "mine_mean": mine_mean.loc[common_taxa],
        "art_mean":  art_mean.loc[common_taxa]
    })
    mm["delta"] = mm["mine_mean"] - mm["art_mean"]
    mm["abs_delta"] = mm["delta"].abs()

    # Correlaciones sobre comunes (valores por taxón, tomando medias)
    from scipy.stats import pearsonr, spearmanr
    pear = pearsonr(mm["mine_mean"].values, mm["art_mean"].values)
    spear = spearmanr(mm["mine_mean"].values, mm["art_mean"].values)
    dfc = pd.DataFrame({
        "metric": ["pearson_r", "pearson_p", "spearman_rho", "spearman_p"],
        "value": [pear.statistic, pear.pvalue, spear.statistic, spear.pvalue],
    })
    return mm.sort_values("mine_mean", ascending=False), dfc

def plot_scatter_means(mm: pd.DataFrame, outdir: Path, rank: str, pseudocount: float = 1e-6):
    """Dispersión de medias (mío vs artículo), escala log10 opcional con pseudocuenta."""
    x = mm["mine_mean"].clip(lower=0) + pseudocount
    y = mm["art_mean"].clip(lower=0) + pseudocount
    fig = px.scatter(mm.reset_index(), x=x, y=y, hover_name="index")
    fig.update_traces(marker=dict(size=10))
    # línea y=x
    minv = float(min(x.min(), y.min()))
    maxv = float(max(x.max(), y.max()))
    line = go.Scatter(x=[minv, maxv], y=[minv, maxv], mode="lines", name="y=x", line=dict(dash="dash"))
    fig.add_trace(line)
    fig.update_layout(
        xaxis_title="Mis medias (rel. abund.) + 1e-6",
        yaxis_title="Artículo medias (rel. abund.) + 1e-6",
        xaxis_type="log", yaxis_type="log",
        title=f"Comparación de medias por {rank} (log10)"
    )
    _save_fig(fig, outdir / f"compare_{rank}_scatter.html", outdir / f"compare_{rank}_scatter.png")

def plot_top_bars(mine: pd.DataFrame, art: pd.DataFrame, outdir: Path, rank: str, topn: int = 20):
    """Barras lado a lado de los top taxones por media (en cualquiera de los dos)."""
    mmean = mine.mean(axis=1)
    amean = art.mean(axis=1)
    taxa_union = mmean.index.union(amean.index)
    comb = pd.DataFrame({"mine": mmean.reindex(taxa_union).fillna(0.0),
                         "art":  amean.reindex(taxa_union).fillna(0.0)})
    # top por suma de ambas medias
    top = comb.assign(sum=comb["mine"] + comb["art"]).sort_values("sum", ascending=False).head(topn)
    df_long = (top.drop(columns=["sum"])
                   .rename_axis(rank)
                   .reset_index()
                   .melt(id_vars=rank, var_name="dataset", value_name="rel_abund"))
    fig = px.bar(df_long, x=rank, y="rel_abund", color="dataset", barmode="group",
                 title=f"Top {topn} {rank} — medias relativas (mío vs artículo)")
    fig.update_layout(xaxis_tickangle=45)
    _save_fig(fig, outdir / f"compare_{rank}_top{topn}_bars.html", outdir / f"compare_{rank}_top{topn}_bars.png")

# ---------------------------------- CLI / MAIN ----------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Compara composiciones taxonómicas (tus resultados vs. artículo) a un rango dado."
    )
    # Mis datos
    p.add_argument("--my-results-dir", type=Path, default=None,
                   help="Carpeta con feature_table_relabund.tsv y taxonomy.tsv.")
    p.add_argument("--my-counts", type=Path, default=None,
                   help="Tabla de conteos (CSV/TSV) con filas=feature_id, columnas=muestras.")
    p.add_argument("--my-taxonomy", type=Path, default=None,
                   help="Taxonomía (CSV/TSV) con columna 'feature_id' y 'taxonomy' o columnas por rank.")
    p.add_argument("--my-relabund-glob", type=str, default=None,
                   help="Glob para tus *_rel-abundance.tsv por muestra.")

    # Artículo
    p.add_argument("--article-counts", type=Path, default=None,
                   help="Tabla de conteos del artículo (CSV/TSV).")
    p.add_argument("--article-taxonomy", type=Path, default=None,
                   help="Taxonomía del artículo (CSV/TSV).")
    p.add_argument("--article-relabund-glob", type=str, default=None,
                   help="Glob a *_rel-abundance.tsv del artículo.")

    # Parámetros
    p.add_argument("--rank", type=str, default="genus", choices=TAX_RANKS,
                   help="Nivel taxonómico para la comparación.")
    p.add_argument("--outdir", type=Path, default=None,
                   help="Carpeta de salida. Por defecto <my-results-dir>/compare_taxa o ./compare_taxa.")
    p.add_argument("--topn", type=int, default=20, help="Top N taxones para barras lado a lado.")
    return p

def main():
    args = build_parser().parse_args()

    # outdir por defecto
    if args.outdir is None:
        base = args.my_results_dir if args.my_results_dir else Path(".")
        args.outdir = base / "compare_taxa"
    outdir = _ensure_dir(args.outdir.resolve())

    # Cargar mis datos
    mine = load_my(
        my_results_dir=args.my_results_dir.resolve() if args.my_results_dir else None,
        rank=args.rank,
        my_counts=args.my_counts,
        my_taxonomy=args.my_taxonomy,
        my_relabund_glob=args.my_relabund_glob,
    )
    eprint(f"[OK] Mis datos: tabla {args.rank} x {mine.shape[1]} muestras (relativas).")

    # Cargar artículo
    art = load_article(
        rank=args.rank,
        article_counts=args.article_counts,
        article_taxonomy=args.article_taxonomy,
        article_relabund_glob=args.article_relabund_glob,
    )
    eprint(f"[OK] Artículo: tabla {args.rank} x {art.shape[1]} muestras (relativas).")

    # Guardar tablas rank
    mine_path = outdir / f"mine_{args.rank}_rel.tsv"
    art_path  = outdir / f"article_{args.rank}_rel.tsv"
    mine.to_csv(mine_path, sep="\t")
    art.to_csv(art_path, sep="\t")

    # Comparar promedios / correlaciones
    means_df, corrs_df = compare_means(mine, art)
    means_path = outdir / f"compare_{args.rank}_means.tsv"
    corr_path  = outdir / f"compare_{args.rank}_correlations.tsv"
    means_df.to_csv(means_path, sep="\t", index=True, header=True)
    corrs_df.to_csv(corr_path, sep="\t", index=False)
    eprint("[OK] Medias y correlaciones guardadas.")

    # Figuras
    try:
        plot_scatter_means(means_df, outdir, args.rank)
    except Exception as e:
        eprint(f"[WARN] No se pudo generar scatter: {e}")
    try:
        plot_top_bars(mine, art, outdir, args.rank, topn=args.topn)
    except Exception as e:
        eprint(f"[WARN] No se pudo generar barras top: {e}")

    eprint(f"[DONE] Resultados en: {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

