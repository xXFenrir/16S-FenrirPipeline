#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
taxonomy_profiling.py — Gráficas de perfil taxonómico (barras apiladas de
abundancia relativa) a nivel de especie, género o familia, mostrando los
Top N taxones más abundantes + 'Otros'. Las muestras se agrupan por Sistema
Agrícola (Empresarial / Campesina / Agroecológica), igual que en rarefaccion.py.

Entradas:
  - tabla_abundancia_relativa.tsv/.xlsx  (id_taxon x muestras, id_taxon
    formato 'rango|tax_id|nombre', normalmente a nivel de especie)
  - taxonomia.tsv (id_taxon, taxonomia) con linaje QIIME-like:
    'k__..; p__..; c__..; o__..; f__..; g__..; s__..'
  - Mapa_Barcodes_Microbioma.csv (barcode->ID) + Sistemas Agrícolas y
    Muestras.xlsx (ID->Sistema), para etiquetar y agrupar las muestras.

Uso típico:
  python3 taxonomy_profiling.py --model hac --dataset todo --rank genus --top-n 10
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches


def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)


RANK_PREFIX = {
    "superkingdom": "k__", "phylum": "p__", "class": "c__",
    "order": "o__", "family": "f__", "genus": "g__", "species": "s__",
}
RANK_ORDER = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


def load_relabund(path: Path) -> pd.DataFrame:
    ext = path.suffix.lower()
    if ext == ".xlsx":
        df = pd.read_excel(path, index_col=0)
    else:
        df = pd.read_csv(path, sep="\t", index_col=0)
    resumen_cols = [c for c in df.columns if str(c).startswith(("Total_counts", "Frecuencia_"))]
    if resumen_cols:
        df = df.drop(columns=resumen_cols)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def load_taxonomy(path: Path) -> dict:
    """Devuelve {id_taxon: {'family': ..., 'genus': ..., 'species': ...}}."""
    df = pd.read_csv(path, sep="\t")
    if "id_taxon" not in df.columns or "taxonomia" not in df.columns:
        raise SystemExit(f"ERROR: {path} debe tener columnas 'id_taxon' y 'taxonomia'.")
    tax_map = {}
    for _, row in df.iterrows():
        fid = str(row["id_taxon"])
        parts = [p.strip() for p in str(row["taxonomia"]).split(";")]
        levels = {}
        for rank in RANK_ORDER:
            prefix = RANK_PREFIX[rank]
            match = next((p[len(prefix):] for p in parts if p.startswith(prefix)), None)
            levels[rank] = match if match and match.lower() not in ("", "unassigned", "nan") else "Unclassified"
        tax_map[fid] = levels
    return tax_map


def rank_label_from_feature_id(fid: str, tax_map: dict, rank: str) -> str:
    """
    Determina la etiqueta de 'rank' para un feature_id. Si el feature_id ya
    ES ese rango (p. ej. viene de una tabla a nivel de especie y rank='species'),
    usa directamente el nombre incluido en el feature_id. Si no, consulta
    taxonomy.tsv para subir/bajar de nivel.
    """
    if fid in tax_map:
        val = tax_map[fid].get(rank)
        if val:
            return val
    # Respaldo: parsear directamente el feature_id 'rango|tax_id|nombre'
    parts = fid.split("|")
    if len(parts) >= 3 and parts[0].lower() == rank:
        return parts[-1]
    return "Unclassified"


def aggregate_by_rank(relabund: pd.DataFrame, tax_map: dict, rank: str) -> pd.DataFrame:
    labels = pd.Series(
        {fid: rank_label_from_feature_id(fid, tax_map, rank) for fid in relabund.index}
    )
    grouped = relabund.groupby(labels).sum()
    return grouped


def top_n_plus_others(grouped: pd.DataFrame, top_n: int) -> pd.DataFrame:
    mean_abund = grouped.mean(axis=1).sort_values(ascending=False)
    top_taxa = [t for t in mean_abund.index if t != "Unclassified"][:top_n]
    keep = top_taxa.copy()
    if "Unclassified" in grouped.index:
        keep.append("Unclassified")
    others = grouped.drop(index=keep, errors="ignore").sum(axis=0)
    out = grouped.loc[[t for t in keep if t in grouped.index]].copy()
    if others.sum() > 0:
        out.loc["Otros"] = others
    # Orden: taxones por abundancia media descendente, 'Unclassified' y 'Otros' al final
    order = [t for t in top_taxa if t in out.index]
    if "Unclassified" in out.index:
        order.append("Unclassified")
    if "Otros" in out.index:
        order.append("Otros")
    return out.loc[order]


def load_double_mapping(bridge_path: Path, meta_path: Path) -> dict:
    """Barcode(int) -> {'id':..., 'sistema':...}. Igual que en rarefaccion.py."""
    mapping = {}
    try:
        df_bridge = pd.read_csv(bridge_path, sep=';')
        if not {"ID", "Barcode"}.issubset(df_bridge.columns):
            df_bridge = pd.read_csv(bridge_path, sep=',')
        barcode_to_id = {}
        for _, row in df_bridge.iterrows():
            try:
                barcode_to_id[int(row["Barcode"])] = str(row["ID"]).strip()
            except (ValueError, TypeError):
                continue
        df_meta = pd.read_excel(meta_path, sheet_name=0, engine='openpyxl').dropna(subset=['ID'])
        id_to_sistema = dict(zip(df_meta['ID'].astype(str).str.strip(), df_meta['Sistema'].astype(str).str.strip()))
        for bc_num, finca_id in barcode_to_id.items():
            mapping[bc_num] = {"id": finca_id, "sistema": id_to_sistema.get(finca_id, "Desconocido")}
    except Exception as e:
        eprint(f"[ERROR] Fallo leyendo puente/metadata: {e}")
    return mapping


def sample_info(col_name: str, mapping: dict):
    m = re.search(r'barcode0*(\d+)', col_name.lower())
    bc_num = int(m.group(1)) if m else None
    info = mapping.get(bc_num) if bc_num is not None else None
    if info:
        return info["id"], info["sistema"]
    return (f"BC{bc_num:02d}" if bc_num is not None else col_name), "Sin_metadata"


def build_color_map(taxa: list) -> dict:
    fixed = {"Otros": "#B0B0B0", "Unclassified": "#7A7A7A"}
    variable = [t for t in taxa if t not in fixed]
    n = len(variable)
    cmap = cm.get_cmap('tab20', max(n, 1)) if n <= 20 else cm.get_cmap('nipy_spectral', n)
    colors = {t: cmap(i) for i, t in enumerate(variable)}
    colors.update({t: fixed[t] for t in fixed if t in taxa})
    return colors


def plot_profile(top_df: pd.DataFrame, mapping: dict, rank: str, top_n: int, model: str, dataset: str, out_path: Path):
    samples = list(top_df.columns)
    info = {s: sample_info(s, mapping) for s in samples}
    order = sorted(samples, key=lambda s: (info[s][1], info[s][0]))
    sistemas = list(dict.fromkeys(info[s][1] for s in order))  # orden de aparición

    color_map = build_color_map(list(top_df.index))

    fig, axes = plt.subplots(
        1, len(sistemas), figsize=(max(10, 0.5 * len(order) + 3), 7),
        sharey=True, dpi=300,
        gridspec_kw={"width_ratios": [sum(1 for s in order if info[s][1] == g) for g in sistemas]}
    )
    if len(sistemas) == 1:
        axes = [axes]

    for ax, sistema in zip(axes, sistemas):
        grupo_samples = [s for s in order if info[s][1] == sistema]
        bottoms = np.zeros(len(grupo_samples))
        for taxon in top_df.index:
            vals = (top_df.loc[taxon, grupo_samples].values * 100.0)
            ax.bar(range(len(grupo_samples)), vals, bottom=bottoms, color=color_map[taxon],
                   edgecolor='white', linewidth=0.3, width=0.85)
            bottoms += vals
        ax.set_xticks(range(len(grupo_samples)))
        ax.set_xticklabels([info[s][0] for s in grupo_samples], rotation=90, fontsize=8)
        ax.set_title(sistema, fontsize=11, fontweight='bold')
        ax.set_ylim(0, 100)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].set_ylabel("Abundancia relativa (%)", fontsize=12, fontweight='bold')
    fig.suptitle(f"Perfil taxonómico ({rank.capitalize()}, Top {top_n}) — {model.upper()} {dataset}",
                 fontsize=14, fontweight='bold')

    handles = [mpatches.Patch(color=color_map[t], label=t) for t in top_df.index]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9)

    plt.tight_layout(rect=[0, 0, 0.85, 0.95])
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    eprint(f"[OK] Gráfica guardada: {out_path}")


def main():
    base = Path("/home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/EMU_propio")

    ap = argparse.ArgumentParser(description="Perfiles taxonómicos (barras apiladas, Top N + Otros).")
    ap.add_argument("--model", choices=["hac", "sup"], default="hac")
    ap.add_argument("--dataset", choices=["todo", "results"], default="todo")
    ap.add_argument("--rank", choices=["family", "genus", "species"], default="genus",
                     help="Nivel taxonómico a graficar.")
    ap.add_argument("--top-n", type=int, default=10, help="Número de taxones más abundantes a mostrar.")
    ap.add_argument("--relabund", type=Path, default=None,
                     help="Por defecto: EMU{model}_{dataset}/tabla_abundancia_relativa.tsv")
    ap.add_argument("--taxonomy", type=Path, default=None,
                     help="Por defecto: EMU{model}_{dataset}/taxonomia.tsv")
    ap.add_argument("--meta", type=Path, default=base / "Diversidad" / "Sistemas Agrícolas y Muestras.xlsx")
    ap.add_argument("--bridge", type=Path,
                     default=Path("/home/fenrir/Documentos/Tesis/datos_gulupa/Mapa Barcodes Microbioma.csv"))
    ap.add_argument("--outdir", type=Path, default=None,
                     help="Por defecto: EMU{model}_{dataset}/figures")
    args = ap.parse_args()

    emu_dir = base / f"EMU{args.model}_{args.dataset}"
    if args.relabund is None:
        args.relabund = emu_dir / "tabla_abundancia_relativa.tsv"
    if args.taxonomy is None:
        args.taxonomy = emu_dir / "taxonomia.tsv"
    if args.outdir is None:
        args.outdir = emu_dir / "figures"
    args.outdir.mkdir(parents=True, exist_ok=True)

    eprint(f"[INFO] Modelo: {args.model} | dataset: {args.dataset} | rank: {args.rank} | top_n: {args.top_n}")
    eprint(f"[INFO] relabund: {args.relabund}")
    eprint(f"[INFO] taxonomy: {args.taxonomy}")

    relabund = load_relabund(args.relabund)
    tax_map = load_taxonomy(args.taxonomy)
    mapping = load_double_mapping(args.bridge, args.meta)
    eprint(f"[INFO] Muestras con ID/Sistema resueltos: {len(mapping)}")

    grouped = aggregate_by_rank(relabund, tax_map, args.rank)
    top_df = top_n_plus_others(grouped, args.top_n)

    out_png = args.outdir / f"taxonomy_profile_{args.rank}_top{args.top_n}_{args.model}_{args.dataset}.png"
    plot_profile(top_df, mapping, args.rank, args.top_n, args.model, args.dataset, out_png)

    out_tsv = args.outdir / f"taxonomy_profile_{args.rank}_top{args.top_n}_{args.model}_{args.dataset}.tsv"
    (top_df * 100.0).round(4).to_csv(out_tsv, sep="\t")
    eprint(f"[OK] Tabla guardada: {out_tsv}")


if __name__ == "__main__":
    main()
