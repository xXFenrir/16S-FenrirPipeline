#!/usr/bin/env python3

# Compara las abundancias relativas medias por taxón entre mis resultados de EMU y los del
# artículo guía, y saca la correlación de Pearson (validación del Objetivo 1).
# Ejemplo:
#   python compare_taxa.py --my-results-dir EMU_out \
#     --article-counts emu_16S_counts_2020.csv --article-taxonomy emu_16S_taxonomy_2020.csv \
#     --rank family genus species
# Para PNG hace falta kaleido, si no solo guarda los HTML.

from __future__ import annotations
import argparse
import sys, os, re, glob
from pathlib import Path
from typing import Optional, Tuple, Dict, List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

_HAS_KALEIDO = False
try:
    import kaleido  # noqa: F401
    _HAS_KALEIDO = True
except Exception:
    pass

TAX_RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_ES = {"superkingdom": "Dominio", "phylum": "Filo", "class": "Clase", "order": "Orden",
           "family": "Familia", "genus": "Género", "species": "Especie"}
ETQ_MIO, ETQ_ART = "Este estudio", "Artículo guía"
COLORES = {ETQ_MIO: "#1f77b4", ETQ_ART: "#d62728"}


def eprint(*a, **k): print(*a, file=sys.stderr, **k)

def _ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def _read_table_any(path: Path, index_col: Optional[int|str] = None) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", index_col=index_col)

def _coerce_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.fillna(0.0)

def _is_sample_like_index(idx: pd.Index) -> bool:
    s = pd.Series(idx.astype(str))
    return s.str.match(r"^(SRR|ERR|DRR)\d+").any() or s.str.contains(r"^barcode", case=False, na=False).any()

# filas = taxones, columnas = muestras
def _maybe_features_rows(df: pd.DataFrame) -> pd.DataFrame:
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
        fig.write_image(str(png_path), scale=3)
    elif png_path:
        eprint(f"[INFO] Instala 'kaleido' para PNG/SVG/PDF. Guardé HTML: {html_path}")


_QIIME_PATTERNS = {
    "superkingdom": re.compile(r"k__([^;]+)"),
    "phylum": re.compile(r"p__([^;]+)"),
    "class": re.compile(r"c__([^;]+)"),
    "order": re.compile(r"o__([^;]+)"),
    "family": re.compile(r"f__([^;]+)"),
    "genus": re.compile(r"g__([^;]+)"),
    "species": re.compile(r"s__([^;]+)"),
}

# estas no son taxones; se quitan y se renormaliza, igual que pasa con las celdas vacías
# de los TSV del artículo (groupby las descarta)
_NO_ASIGNADO = {"", "nan", "none", "na", "unassigned", "uncultured", "metagenome"}

def _drop_unassigned(G: pd.DataFrame) -> pd.DataFrame:
    keep = ~G.index.astype(str).str.strip().str.lower().isin(_NO_ASIGNADO)
    return G.loc[keep]

def _from_qiime(s: str, rank: str) -> str:
    if not isinstance(s, str):
        return "Unassigned"
    m = _QIIME_PATTERNS[rank].search(s)
    if not m:
        return "Unassigned"
    val = m.group(1).strip()
    if val.lower() in _NO_ASIGNADO:
        return "Unassigned"
    return val

def taxonomy_to_rank_map(tax_df: pd.DataFrame, rank: str) -> pd.Series:
    df = tax_df.copy()
    if "feature_id" in df.columns:
        df = df.set_index("feature_id")
    if rank in df.columns:
        s = df[rank].astype(str)
    elif "taxonomy" in df.columns:
        s = df["taxonomy"].map(lambda x: _from_qiime(str(x), rank))
    else:
        # columnas por rango, se arma la cadena tipo QIIME
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


def detect_rel_col(df: pd.DataFrame) -> Optional[str]:
    # la que sume ~1
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]):
            s = pd.to_numeric(df[c], errors="coerce").fillna(0).sum()
            if 0.9 <= s <= 1.1:
                return c
    cand = [c for c in df.columns if re.search(r"abund", c, flags=re.I)]
    for c in cand:
        s = pd.to_numeric(df[c], errors="coerce").fillna(0).sum()
        if 0.9 <= s <= 1.1:
            return c
    return None


def counts_to_relative(counts: pd.DataFrame) -> pd.DataFrame:
    C = _coerce_numeric_df(counts)
    if "feature_id" in C.columns:
        C = C.set_index("feature_id")
    C = _maybe_features_rows(C)
    colsum = C.sum(axis=0)
    colsum[colsum == 0] = 1.0
    return C / colsum

def aggregate_by_rank_rel(rel: pd.DataFrame, taxmap: pd.Series) -> pd.DataFrame:
    R = _maybe_features_rows(rel.copy()).fillna(0.0)
    common = R.index.intersection(taxmap.index)
    if common.empty:
        # puede venir traspuesta
        R = R.T
        common = R.index.intersection(taxmap.index)
    if common.empty:
        raise RuntimeError("No coinciden los feature_id entre tabla y taxonomía.")
    R = R.loc[common]
    taxmap = taxmap.loc[common]
    G = _drop_unassigned(R.groupby(taxmap).sum())
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
    # se prueba en orden: counts + taxonomy, carpeta de resultados, TSV por muestra
    if my_counts and my_taxonomy:
        C = _read_table_any(Path(my_counts), index_col=0)
        T = _read_table_any(Path(my_taxonomy))
        taxmap = taxonomy_to_rank_map(T, rank=rank)
        R = counts_to_relative(C)
        return aggregate_by_rank_rel(R, taxmap)

    if my_results_dir:
        rel_p = my_results_dir / "feature_table_relabund.tsv"
        tax_p = my_results_dir / "taxonomy.tsv"
        if rel_p.exists() and tax_p.exists():
            Rel = _read_table_any(rel_p, index_col=0)
            T = _read_table_any(tax_p)
            taxmap = taxonomy_to_rank_map(T, rank=rank)
            return aggregate_by_rank_rel(Rel, taxmap)

    if my_relabund_glob:
        paths = glob.glob(my_relabund_glob)
        if not paths:
            raise RuntimeError(f"No encontré TSV con el patrón: {my_relabund_glob}")
        tabs = []
        need_tax = False
        for p in paths:
            df = _read_table_any(Path(p))
            relc = detect_rel_col(df)
            if relc is None:
                raise RuntimeError(f"No detecté columna de relativas en: {p}")
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
            if not my_taxonomy:
                raise RuntimeError("Tus TSV no traen columna rank ni 'taxonomy'. Pasa --my-taxonomy para mapear.")
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
        G = _drop_unassigned(pd.concat(tabs, axis=1).fillna(0.0))
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
    if article_counts and article_taxonomy:
        C = _read_table_any(Path(article_counts), index_col=0)
        T = _read_table_any(Path(article_taxonomy))
        taxmap = taxonomy_to_rank_map(T, rank=rank)
        R = counts_to_relative(C)
        return aggregate_by_rank_rel(R, taxmap)

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
        G = _drop_unassigned(pd.concat(tabs, axis=1).fillna(0.0))
        colsum = G.sum(axis=0); colsum[colsum == 0] = 1.0
        return G / colsum

    raise RuntimeError(
        "Debes pasar --article-relabund-glob o --article-counts + --article-taxonomy para el artículo."
    )

def _sample_name_from_path(p: Path) -> str:
    m = re.search(r"(barcode\d+)", p.name, flags=re.I)
    if m:
        return m.group(1)
    # nombres tipo SRR/ERR/DRR
    m = re.search(r"(SRR|ERR|DRR)\d{5,}", p.name)
    if m:
        return m.group(0)
    stem = p.stem
    stem = re.sub(r"\.t_?rel-?abundance.*$", "", stem, flags=re.I)
    return stem


def compare_means(mine: pd.DataFrame, art: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    mine_mean = mine.mean(axis=1) if mine.shape[1] > 0 else pd.Series(dtype=float)
    art_mean  = art.mean(axis=1) if art.shape[1] > 0 else pd.Series(dtype=float)
    # solo taxones en común
    common_taxa = mine_mean.index.intersection(art_mean.index)
    if common_taxa.empty:
        raise RuntimeError("No hay taxones en común para comparar promedios.")
    mm = pd.DataFrame({
        "mine_mean": mine_mean.loc[common_taxa],
        "art_mean":  art_mean.loc[common_taxa]
    })
    mm["delta"] = mm["mine_mean"] - mm["art_mean"]
    mm["abs_delta"] = mm["delta"].abs()

    # Pearson sobre las medias por taxón
    from scipy.stats import pearsonr
    r, p = pearsonr(mm["mine_mean"].values, mm["art_mean"].values)
    dfc = pd.DataFrame({
        "metric": ["pearson_r", "pearson_p", "n_taxones_comunes"],
        "value": [r, p, len(mm)],
    })
    return mm.sort_values("mine_mean", ascending=False), dfc

def plot_scatter_means(mm: pd.DataFrame, outdir: Path, rank: str, r: Optional[float] = None,
                       pseudocount: float = 1e-6):
    # pseudocuenta para poder usar escala log
    x = (mm["mine_mean"].clip(lower=0) + pseudocount) * 100
    y = (mm["art_mean"].clip(lower=0) + pseudocount) * 100
    fig = px.scatter(mm.reset_index(), x=x, y=y, hover_name=mm.index)
    fig.update_traces(marker=dict(size=9, color=COLORES[ETQ_MIO], opacity=0.75), name=RANK_ES[rank],
                      showlegend=False)
    minv = float(min(x.min(), y.min()))
    maxv = float(max(x.max(), y.max()))
    line = go.Scatter(x=[minv, maxv], y=[minv, maxv], mode="lines", name="Concordancia perfecta (y = x)",
                      line=dict(dash="dash", color="gray"))
    fig.add_trace(line)
    fig.update_layout(
        xaxis_title=f"Abundancia relativa media — {ETQ_MIO.lower()} (%)",
        yaxis_title=f"Abundancia relativa media — {ETQ_ART.lower()} (%)",
        xaxis_type="log", yaxis_type="log",
        title=f"Abundancia relativa media por {RANK_ES[rank].lower()}: {ETQ_MIO.lower()} vs. {ETQ_ART.lower()}",
        legend=dict(yanchor="top", y=0.98, xanchor="left", x=0.02),
    )
    if r is not None:
        fig.add_annotation(xref="paper", yref="paper", x=0.98, y=0.04, showarrow=False,
                           text=f"Pearson r = {r:.3f}", font=dict(size=16))
    _save_fig(fig, outdir / f"compare_{rank}_scatter.html", outdir / f"compare_{rank}_scatter.png")

def plot_top_bars(mine: pd.DataFrame, art: pd.DataFrame, outdir: Path, rank: str, topn: int = 20):
    mmean = mine.mean(axis=1)
    amean = art.mean(axis=1)
    taxa_union = mmean.index.union(amean.index)
    comb = pd.DataFrame({"mine": mmean.reindex(taxa_union).fillna(0.0),
                         "art":  amean.reindex(taxa_union).fillna(0.0)})
    # top por la suma de las dos medias
    top = comb.assign(sum=comb["mine"] + comb["art"]).sort_values("sum", ascending=False).head(topn)
    col_rank = RANK_ES[rank]
    df_long = (top.drop(columns=["sum"])
                   .rename(columns={"mine": ETQ_MIO, "art": ETQ_ART}) * 100)
    df_long = (df_long.rename_axis(col_rank)
                   .reset_index()
                   .melt(id_vars=col_rank, var_name="Conjunto de datos",
                         value_name="Abundancia relativa media (%)"))
    fig = px.bar(df_long, x=col_rank, y="Abundancia relativa media (%)", color="Conjunto de datos",
                 barmode="group", color_discrete_map=COLORES,
                 title=f"{topn} taxones más abundantes a nivel de {col_rank.lower()}: "
                       f"{ETQ_MIO.lower()} vs. {ETQ_ART.lower()}")
    fig.update_layout(xaxis_tickangle=45)
    _save_fig(fig, outdir / f"compare_{rank}_top{topn}_bars.html", outdir / f"compare_{rank}_top{topn}_bars.png")


def build_parser():
    p = argparse.ArgumentParser(
        description="Compara composiciones taxonómicas (tus resultados vs. artículo) a un rango dado."
    )
    p.add_argument("--my-results-dir", type=Path, default=None,
                   help="Carpeta con feature_table_relabund.tsv y taxonomy.tsv.")
    p.add_argument("--my-counts", type=Path, default=None,
                   help="Tabla de conteos (CSV/TSV) con filas=feature_id, columnas=muestras.")
    p.add_argument("--my-taxonomy", type=Path, default=None,
                   help="Taxonomía (CSV/TSV) con columna 'feature_id' y 'taxonomy' o columnas por rank.")
    p.add_argument("--my-relabund-glob", type=str, default=None,
                   help="Glob para tus *_rel-abundance.tsv por muestra.")

    p.add_argument("--article-counts", type=Path, default=None,
                   help="Tabla de conteos del artículo (CSV/TSV).")
    p.add_argument("--article-taxonomy", type=Path, default=None,
                   help="Taxonomía del artículo (CSV/TSV).")
    p.add_argument("--article-relabund-glob", type=str, default=None,
                   help="Glob a *_rel-abundance.tsv del artículo.")

    p.add_argument("--rank", type=str, nargs="+", default=["genus"], choices=TAX_RANKS,
                   help="Uno o varios niveles taxonómicos (p. ej. family genus species).")
    p.add_argument("--outdir", type=Path, default=None,
                   help="Carpeta de salida. Por defecto <my-results-dir>/compare_taxa o ./compare_taxa.")
    p.add_argument("--topn", type=int, default=20, help="Top N taxones para barras lado a lado.")
    return p

def main():
    args = build_parser().parse_args()

    if args.outdir is None:
        base = args.my_results_dir if args.my_results_dir else Path(".")
        args.outdir = base / "compare_taxa"
    outdir = _ensure_dir(args.outdir.resolve())

    resumen = {}
    for rank in args.rank:
        eprint(f"\n===== Rango: {rank} =====")
        mine = load_my(
            my_results_dir=args.my_results_dir.resolve() if args.my_results_dir else None,
            rank=rank,
            my_counts=args.my_counts,
            my_taxonomy=args.my_taxonomy,
            my_relabund_glob=args.my_relabund_glob,
        )
        eprint(f"[OK] Mis datos: tabla {rank} x {mine.shape[1]} muestras (relativas).")

        art = load_article(
            rank=rank,
            article_counts=args.article_counts,
            article_taxonomy=args.article_taxonomy,
            article_relabund_glob=args.article_relabund_glob,
        )
        eprint(f"[OK] Artículo: tabla {rank} x {art.shape[1]} muestras (relativas).")

        mine.to_csv(outdir / f"mine_{rank}_rel.tsv", sep="\t")
        art.to_csv(outdir / f"article_{rank}_rel.tsv", sep="\t")

        means_df, corrs_df = compare_means(mine, art)
        means_df.to_csv(outdir / f"compare_{rank}_means.tsv", sep="\t", index=True, header=True)
        corrs_df.to_csv(outdir / f"compare_{rank}_correlations.tsv", sep="\t", index=False)
        r = float(corrs_df.set_index("metric").loc["pearson_r", "value"])
        resumen[RANK_ES[rank]] = r
        eprint(f"[OK] Pearson r ({rank}) = {r:.4f}")

        try:
            plot_scatter_means(means_df, outdir, rank, r=r)
        except Exception as e:
            eprint(f"[WARN] No se pudo generar scatter: {e}")
        try:
            plot_top_bars(mine, art, outdir, rank, topn=args.topn)
        except Exception as e:
            eprint(f"[WARN] No se pudo generar barras top: {e}")

    # tabla de Pearson con una columna por rango
    tabla = pd.DataFrame([resumen], index=["Pearson"])
    tabla.to_csv(outdir / "correlacion_pearson.tsv", sep="\t")
    try:
        tabla.round(4).to_excel(outdir / "correlacion_pearson.xlsx")
    except Exception as e:
        eprint(f"[WARN] No se pudo escribir XLSX: {e}")
    eprint("\n" + tabla.round(4).to_string())

    eprint(f"[DONE] Resultados en: {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

