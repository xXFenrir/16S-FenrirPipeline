#!/usr/bin/env python3

# Compara, muestra por muestra, los resultados del pipeline con los del artículo guía
# (Erlandson et al., 2024; PRJNA1020132) en cada paso: datos iniciales, pre y post limpieza,
# lecturas asignadas por EMU, composición, diversidad alfa y beta, y PERMANOVA.
# Del artículo se usan los archivos de su repositorio (github.com/serlandson/sterile_sentinels):
# las salidas de EMU por barcode, su tabla de conteos final, los datos de muestra y sus
# resultados de adonis. Para 2020 no publicaron estadísticas por etapa, así que las lecturas
# de cada paso se comparan contra sus lecturas asignadas por EMU.

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import braycurtis, pdist, squareform
from scipy.stats import pearsonr
from skbio import DistanceMatrix
from skbio.stats.distance import mantel

YO, ART, GRIS, TINTA, TINTA2 = "#2a78d6", "#eb6834", "#9a9993", "#0b0b0b", "#52514e"


def leer_estadisticas(ruta: Path) -> pd.DataFrame:
    # sirve con los excel de los anexos 3 y 4 (encabezado corrido) o con la salida de stats_fastq.py
    if ruta.suffix.lower() == ".xlsx":
        crudo = pd.read_excel(ruta, header=None)
    else:
        crudo = pd.read_csv(ruta, sep="\t", header=None)
    fila = next(i for i, r in crudo.iterrows() if r.astype(str).str.strip().eq("Muestra").any())
    df = crudo.iloc[fila + 1:].copy()
    df.columns = crudo.iloc[fila].astype(str).str.strip()
    df = df.loc[:, ~df.columns.isin(["nan"])]
    df["Run"] = df["Muestra"].astype(str).str.extract(r"(SRR\d+)", expand=False)
    df = df[df.Run.notna()]
    for c in ("lecturas", "bases"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df[["Run", "lecturas", "bases"]]


def leer_metadatos(rutas) -> pd.DataFrame:
    # el SraRunTable.csv y el excel del Anexo 2 no traen exactamente las mismas corridas
    # (al csv le falta S1A2 y al excel S1A8), así que se unen los que se pasen
    tablas = []
    for ruta in rutas:
        if ruta.suffix.lower() == ".xlsx":
            crudo = pd.read_excel(ruta, header=None)
            fila = next(i for i, r in crudo.iterrows() if r.astype(str).str.strip().eq("Run").any())
            df = crudo.iloc[fila + 1:].copy()
            df.columns = crudo.iloc[fila].astype(str).str.strip()
        else:
            df = pd.read_csv(ruta)
        df = df[df["Run"].astype(str).str.startswith("SRR")]
        tablas.append(df[["Run", "Sample Name", "Bases"]])
    meta = pd.concat(tablas).drop_duplicates("Run")
    meta["Bases"] = pd.to_numeric(meta["Bases"])
    return meta.reset_index(drop=True)


def shannon(tabla: pd.DataFrame) -> pd.Series:
    p = tabla.div(tabla.sum(axis=1), axis=0)
    return -(p * np.log(p.where(p > 0, 1))).sum(axis=1)


def clr(tabla: pd.DataFrame) -> pd.DataFrame:
    # +1 igual que el artículo: decostand(comm + 1, method = "clr")
    x = np.log(tabla + 1)
    return x.sub(x.mean(axis=1), axis=0)


def distancias(tabla: pd.DataFrame, metrica: str = "euclidean") -> DistanceMatrix:
    return DistanceMatrix(squareform(pdist(tabla.values, metrica)), ids=list(tabla.index))


def permanova_vegan(dm: DistanceMatrix, grupos, permutaciones: int = 999, semilla: int = 42):
    # PERMANOVA de una vía como adonis2 de vegan. Con 8 muestras solo hay 35 formas de partirlas
    # en dos grupos, así que muchas permutaciones dan el mismo F que el observado; vegan compara
    # con una tolerancia (EPS) y sin ella el valor p cambia según el redondeo en cada corrida
    d2 = dm.data ** 2
    g = np.asarray(grupos)
    n, k = len(g), len(np.unique(g))
    sst = d2[np.triu_indices(n, 1)].sum() / n

    def ss_dentro(etq):
        total = 0.0
        for e in np.unique(etq):
            idx = np.where(etq == e)[0]
            total += d2[np.ix_(idx, idx)][np.triu_indices(len(idx), 1)].sum() / len(idx)
        return total

    def pseudo_f(etq):
        ssw = ss_dentro(etq)
        return ((sst - ssw) / (k - 1)) / (ssw / (n - k)), (sst - ssw) / sst

    F, r2 = pseudo_f(g)
    rng = np.random.default_rng(semilla)
    f_perm = np.array([pseudo_f(rng.permutation(g))[0] for _ in range(permutaciones)])
    eps = np.sqrt(np.finfo(float).eps)
    p = ((f_perm >= F - eps).sum() + 1) / (permutaciones + 1)
    return F, r2, p


def cargar_datos(args):
    art = args.articulo / "data"
    sra = leer_metadatos(args.sra)
    sra["sample_id"] = sra["Sample Name"].str.replace("_bac", "", regex=False)
    muestras_art = pd.read_csv(art / "16S_sample_data_2020.csv")
    m = sra.merge(muestras_art[["sample_id", "bacteria_barcode", "type", "time", "rotation"]], on="sample_id", how="inner")

    pre = leer_estadisticas(args.pre).rename(columns={"lecturas": "lect_iniciales", "bases": "bases_iniciales"})
    post = leer_estadisticas(args.post).rename(columns={"lecturas": "lect_limpias"})[["Run", "lect_limpias"]]
    m = m.merge(pre, on="Run").merge(post, on="Run")

    conteos = pd.read_csv(args.conteos, sep="\t", index_col=0)
    conteos.index = [str(i).split("|")[1] if "|" in str(i) else str(i) for i in conteos.index]
    conteos = conteos.groupby(level=0).sum()
    m = m[m.Run.isin(conteos.columns)].copy()

    # salidas de EMU del artículo por barcode, sin los "unassigned"
    por_muestra = {}
    for _, fila in m.iterrows():
        d = pd.read_csv(art / "bacteria raw data 2020" / f"{fila.bacteria_barcode}.t_rel-abundance.tsv", sep="\t")
        d = d[d.tax_id.astype(str) != "unassigned"]
        por_muestra[fila.Run] = d.groupby(d.tax_id.astype(str))["estimated counts"].sum()
    A = pd.DataFrame(por_muestra).fillna(0)
    U = conteos[m.Run]
    taxones = A.index.union(U.index)
    A, U = A.reindex(taxones, fill_value=0), U.reindex(taxones, fill_value=0)

    m["emu_yo"] = m.Run.map(U.sum())
    m["emu_art"] = m.Run.map(A.sum())
    m = m.reset_index(drop=True)
    return m, U, A


def analizar(args, m, U, A):
    res = {}
    Ur, Ar = U / U.sum(), A / A.sum()
    runs = list(m.Run)

    # composición: misma muestra en los dos análisis vs muestras distintas
    res["bc_misma"] = np.array([braycurtis(Ur[s], Ar[s]) for s in runs])
    res["bc_distintas"] = np.array([braycurtis(Ur[a], Ar[b]) for a in runs for b in runs if a != b])

    # tabla final del artículo y la mía con su mismo filtro (>4 lecturas y presente en >=3 muestras)
    tabla_art = pd.read_csv(args.articulo / "data" / "emu_16S_counts_2020.csv", index_col=0).loc[m.bacteria_barcode]
    tabla_art.index = runs
    tabla_yo = U.T.copy()
    tabla_yo = tabla_yo.loc[:, tabla_yo.sum() > 4]
    tabla_yo = np.ceil(tabla_yo)
    tabla_yo = tabla_yo.loc[:, (tabla_yo > 0).sum() >= 3]

    res["sh_yo"], res["sh_art"] = shannon(tabla_yo), shannon(tabla_art)
    res["obs_yo"], res["obs_art"] = (tabla_yo > 0).sum(axis=1), (tabla_art > 0).sum(axis=1)

    res["ait_yo"], res["ait_art"] = distancias(clr(tabla_yo)), distancias(clr(tabla_art))
    res["mantel_ait"] = mantel(res["ait_yo"], res["ait_art"], permutations=999, seed=42)
    res["mantel_bc"] = mantel(distancias(Ur.T, "braycurtis"), distancias(Ar.T, "braycurtis"), permutations=999, seed=42)

    # PERMANOVA por tipo y tiempo, como en permANOVA_functions.R del artículo (CLR + euclidiana)
    adonis = pd.read_csv(args.articulo / "data" / "bacteria_2020_adonis_results.csv", index_col=0)
    filas = []
    for _, r in adonis.iterrows():
        sel = m[(m.type == r.type) & (m.time == r.time)]
        fila = dict(Tipo=r.type, Semana=r.week, n=len(sel), R2_publicado=r.rotation_R2, p_publicado=r.rotation_p_value)
        for etiqueta, tabla in (("este_trabajo", tabla_yo), ("tabla_articulo", tabla_art)):
            dm = distancias(clr(tabla.loc[sel.Run]))
            _, r2, p = permanova_vegan(dm, sel.rotation.values)
            fila[f"R2_{etiqueta}"] = r2
            fila[f"p_{etiqueta}"] = p
        filas.append(fila)
    res["permanova"] = pd.DataFrame(filas)
    return res


def figura(m, res, salida: Path):
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA,
                         "xtick.color": TINTA2, "ytick.color": TINTA2, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                         "grid.linewidth": 0.6, "axes.axisbelow": True})

    def diagonal(ax, lo, hi):
        ax.plot([lo, hi], [lo, hi], ls="--", lw=1, color=GRIS, zorder=1)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)

    def nota(ax, texto, abajo=False):
        x, y, ha, va = (0.97, 0.03, "right", "bottom") if abajo else (0.03, 0.97, "left", "top")
        ax.text(x, y, texto, transform=ax.transAxes, ha=ha, va=va, fontsize=8.5, color=TINTA,
                bbox=dict(fc="white", ec="none", alpha=0.85))

    def puntos(ax, x, y):
        ax.scatter(x, y, s=22, color=YO, ec="white", lw=0.6, zorder=3)

    fig, axs = plt.subplots(2, 4, figsize=(17.5, 9))
    axs = axs.ravel()

    ax = axs[0]
    puntos(ax, m.Bases / 1e6, m.bases_iniciales / 1e6)
    diagonal(ax, 15, 65)
    ax.set_xlabel("Bases según el BioProject (Mb)"); ax.set_ylabel("Bases en los FASTQ descargados (Mb)")
    ax.set_title("A. Datos iniciales", loc="left", fontweight="bold")
    nota(ax, f"{int((m.Bases == m.bases_iniciales).sum())} de {len(m)} muestras idénticas")

    # B, C y D contra las lecturas asignadas por EMU del artículo, que es lo único por muestra que publicaron
    def etapa(ax, col, titulo, ylab, extra, misma_etapa):
        x, y = m.emu_art / 1e3, m[col] / 1e3
        puntos(ax, x, y)
        hi = max(x.max(), y.max()) * 1.08
        k = 1.0
        if misma_etapa:
            diagonal(ax, 0, hi)
        else:
            k = (x * y).sum() / (x * x).sum()
            ax.plot([0, hi], [0, k * hi], color=GRIS, lw=1.2, zorder=1)
            ax.set_xlim(0, hi); ax.set_ylim(0, hi)
            extra = f"pendiente = {k:.2f} (línea gris, ajuste por el origen)\n" + extra
        ax.set_xlabel("Artículo: lecturas asignadas por EMU (miles)"); ax.set_ylabel(ylab)
        ax.set_title(titulo, loc="left", fontweight="bold")
        # si la nube queda por encima de la diagonal, la nota va abajo para no tapar puntos
        nota(ax, f"r = {pearsonr(m.emu_art, m[col])[0]:.3f}\n" + extra, abajo=k > 1)

    ret = m.lect_limpias / m.lect_iniciales
    etapa(axs[1], "lect_iniciales", "B. Pre-limpieza", "Este trabajo: lecturas iniciales (miles)", "", False)
    etapa(axs[2], "lect_limpias", "C. Post-limpieza", "Este trabajo: lecturas limpias (miles)",
          f"retención: {100 * ret.mean():.0f} % ({100 * ret.min():.0f}–{100 * ret.max():.0f} %)", False)
    etapa(axs[3], "emu_yo", "D. Lecturas asignadas por EMU", "Este trabajo: lecturas asignadas (miles)",
          f"de las lecturas iniciales: este trabajo {100 * (m.emu_yo / m.lect_iniciales).mean():.0f} %,\n"
          f"artículo {100 * (m.emu_art / m.lect_iniciales).mean():.0f} %", True)

    ax = axs[4]
    bins = np.linspace(0, 1, 41)
    ax.hist(res["bc_distintas"], bins=bins, density=True, color=GRIS, alpha=0.55, label="Muestras distintas", ec="white", lw=0.5)
    ax.hist(res["bc_misma"], bins=bins, density=True, color=YO, alpha=0.9, label="Misma muestra", ec="white", lw=0.5)
    ax.set_xlabel("Disimilitud de Bray-Curtis (especie)"); ax.set_ylabel("Densidad")
    ax.set_title("E. Composición", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="center right", fontsize=8)
    nota(ax, f"mediana misma muestra: {np.median(res['bc_misma']):.2f}\n"
             f"mediana muestras distintas: {np.median(res['bc_distintas']):.2f}")

    ax = axs[5]
    puntos(ax, res["sh_art"], res["sh_yo"])
    lo = min(res["sh_art"].min(), res["sh_yo"].min()) - 0.1
    hi = max(res["sh_art"].max(), res["sh_yo"].max()) + 0.1
    diagonal(ax, lo, hi)
    ax.set_xlabel("Shannon, artículo"); ax.set_ylabel("Shannon, este trabajo")
    ax.set_title("F. Diversidad alfa por muestra", loc="left", fontweight="bold")
    nota(ax, f"r = {pearsonr(res['sh_art'], res['sh_yo'])[0]:.3f}\n"
             f"riqueza observada: r = {pearsonr(res['obs_art'], res['obs_yo'])[0]:.3f}")

    ax = axs[6]
    da, du = res["ait_art"].condensed_form(), res["ait_yo"].condensed_form()
    ax.hexbin(da, du, gridsize=45, cmap="Blues", mincnt=1, linewidths=0)
    diagonal(ax, min(da.min(), du.min()), max(da.max(), du.max()))
    ax.set_xlabel("Distancia de Aitchison, artículo"); ax.set_ylabel("Distancia de Aitchison, este trabajo")
    ax.set_title("G. Diversidad beta (pares de muestras)", loc="left", fontweight="bold")
    r_ait, p_ait, _ = res["mantel_ait"]
    nota(ax, f"Mantel r = {r_ait:.3f} (p = {p_ait:.3f})\nBray-Curtis: Mantel r = {res['mantel_bc'][0]:.3f}")

    ax = axs[7]
    P = res["permanova"].copy()
    P["etiqueta"] = P.Tipo.map({"bulk": "Suelo", "bag": "Bolsa"}) + " s" + P.Semana.astype(str)
    yy = np.arange(len(P))[::-1]
    for i, (_, r) in enumerate(P.iterrows()):
        ax.plot([r.R2_publicado, r.R2_este_trabajo], [yy[i]] * 2, color="#d4d3cd", lw=1.5, zorder=1)
    for col, c, mk in (("publicado", ART, "o"), ("este_trabajo", YO, "D")):
        sig = (P[f"p_{col}"] < 0.05).values
        ax.scatter(P[f"R2_{col}"][~sig], yy[~sig], s=40, facecolor="white", edgecolor=c, lw=1.6, marker=mk, zorder=3)
        ax.scatter(P[f"R2_{col}"][sig], yy[sig], s=40, color=c, marker=mk, zorder=3)
    ax.set_yticks(yy); ax.set_yticklabels(P.etiqueta); ax.grid(axis="y", visible=False)
    ax.set_ylim(-0.7, len(P) + 0.6)
    ax.set_xlabel("R² del sistema de rotación (PERMANOVA)")
    ax.set_title("H. Efecto de la rotación (PERMANOVA)", loc="left", fontweight="bold")
    leyenda = [Line2D([], [], marker="o", ls="", color=ART, label="Publicado (artículo)"),
               Line2D([], [], marker="D", ls="", color=YO, label="Este trabajo"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=TINTA2, label="relleno: p < 0.05; vacío: p ≥ 0.05")]
    ax.legend(handles=leyenda, frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=2)
    iguales = ((P.p_publicado < 0.05) == (P.p_este_trabajo < 0.05)).sum()
    ax.text(0.97, 0.99, f"misma conclusión (α = 0.05)\nen {iguales} de {len(P)} subconjuntos", transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, color=TINTA)

    fig.suptitle(f"Comparación paso a paso con el artículo guía (PRJNA1020132, {len(m)} muestras)",
                 x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(salida, dpi=300)
    plt.close(fig)
    return iguales


def tablas(m, res, iguales, salida: Path):
    P = res["permanova"]
    ret = m.lect_limpias / m.lect_iniciales
    r = lambda a, b: pearsonr(m[a], m[b])[0]
    resumen = pd.DataFrame([
        ["1. Datos iniciales", "Bases por muestra: BioProject vs FASTQ descargados",
         f"{int((m.Bases == m.bases_iniciales).sum())}/{len(m)} idénticas"],
        ["2. Pre-limpieza", "Lecturas iniciales por muestra vs lecturas en EMU del artículo (Pearson)",
         f"r = {r('emu_art', 'lect_iniciales'):.4f}"],
        ["3. Post-limpieza", "Lecturas limpias por muestra vs lecturas en EMU del artículo (Pearson)",
         f"r = {r('emu_art', 'lect_limpias'):.4f}"],
        ["3. Post-limpieza", "Retención de la limpieza (media, mín–máx)",
         f"{100 * ret.mean():.1f} % ({100 * ret.min():.1f}–{100 * ret.max():.1f} %)"],
        ["4. Clasificación EMU", "Lecturas asignadas por muestra (Pearson)", f"r = {r('emu_art', 'emu_yo'):.4f}"],
        ["4. Clasificación EMU", "Lecturas asignadas / lecturas iniciales",
         f"este trabajo {100 * (m.emu_yo / m.lect_iniciales).mean():.1f} %; "
         f"artículo {100 * (m.emu_art / m.lect_iniciales).mean():.1f} %"],
        ["5. Composición", "Bray-Curtis misma muestra vs muestras distintas (mediana)",
         f"{np.median(res['bc_misma']):.3f} vs {np.median(res['bc_distintas']):.3f}"],
        ["6. Diversidad alfa", "Shannon por muestra (Pearson)", f"r = {pearsonr(res['sh_art'], res['sh_yo'])[0]:.4f}"],
        ["6. Diversidad alfa", "Riqueza observada por muestra (Pearson)",
         f"r = {pearsonr(res['obs_art'], res['obs_yo'])[0]:.4f}"],
        ["7. Diversidad beta", "Mantel, distancias de Aitchison",
         f"r = {res['mantel_ait'][0]:.4f}, p = {res['mantel_ait'][1]:.3f}"],
        ["7. Diversidad beta", "Mantel, distancias de Bray-Curtis",
         f"r = {res['mantel_bc'][0]:.4f}, p = {res['mantel_bc'][1]:.3f}"],
        ["8. PERMANOVA", "Subconjuntos con la misma conclusión (α = 0.05)", f"{iguales} de {len(P)}"],
    ], columns=["Paso", "Comparación", "Resultado"])

    por_muestra = pd.DataFrame({
        "SRR": m.Run, "Muestra": m["Sample Name"], "Barcode artículo": m.bacteria_barcode,
        "Tipo": m.type, "Tiempo": m.time, "Rotación": m.rotation,
        "Bases BioProject": m.Bases, "Bases descargadas": m.bases_iniciales,
        "Lecturas iniciales": m.lect_iniciales, "Lecturas limpias (este trabajo)": m.lect_limpias,
        "Lecturas a EMU (este trabajo)": m.emu_yo, "Lecturas a EMU (artículo)": m.emu_art,
        "Shannon (este trabajo)": m.Run.map(res["sh_yo"]), "Shannon (artículo)": m.Run.map(res["sh_art"]),
        "Bray-Curtis misma muestra": res["bc_misma"],
    })
    perm = P.rename(columns={"R2_publicado": "R² publicado", "p_publicado": "p publicado",
                             "R2_este_trabajo": "R² este trabajo", "p_este_trabajo": "p este trabajo",
                             "R2_tabla_articulo": "R² tabla del artículo (recalculado)",
                             "p_tabla_articulo": "p tabla del artículo (recalculado)"})
    with pd.ExcelWriter(salida) as w:
        resumen.to_excel(w, sheet_name="Resumen", index=False)
        por_muestra.to_excel(w, sheet_name="Por_muestra", index=False)
        perm.to_excel(w, sheet_name="PERMANOVA", index=False)
        # ancho de columnas según el contenido, para que se lea sin ajustar a mano
        for hoja in w.sheets.values():
            for col in hoja.columns:
                largo = max(len(str(c.value)) if c.value is not None else 0 for c in col)
                hoja.column_dimensions[col[0].column_letter].width = min(largo + 2, 60)
    return resumen


def main():
    ap = argparse.ArgumentParser(description="Comparación paso a paso con el artículo guía (PRJNA1020132).")
    ap.add_argument("--sra", type=Path, nargs="+", required=True,
                    help="Metadatos del BioProject (Anexo 2): SraRunTable.csv y/o el excel; si se pasan los dos se unen.")
    ap.add_argument("--pre", type=Path, required=True, help="Estadísticas antes de la limpieza (Anexo 3 o salida de stats_fastq.py).")
    ap.add_argument("--post", type=Path, required=True, help="Estadísticas después de la limpieza (Anexo 4 o salida de stats_fastq.py).")
    ap.add_argument("--conteos", type=Path, required=True, help="Tabla de conteos de EMU por especie (tabla_conteos.tsv).")
    ap.add_argument("--articulo", type=Path, required=True, help="Carpeta clonada de github.com/serlandson/sterile_sentinels.")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("comparacion_articulo"), help="Carpeta de salida.")
    args = ap.parse_args()

    if not (args.articulo / "data" / "bacteria raw data 2020").is_dir():
        sys.exit(f"[!] No encuentro 'data/bacteria raw data 2020' dentro de {args.articulo}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    m, U, A = cargar_datos(args)
    print(f"[OK] Muestras con datos en los dos análisis: {len(m)}")
    res = analizar(args, m, U, A)
    iguales = figura(m, res, args.outdir / "comparacion_articulo_pasos.png")
    resumen = tablas(m, res, iguales, args.outdir / "comparacion_articulo_pasos.xlsx")
    print(resumen.to_string(index=False))
    print(f"\n[OK] Resultados en: {args.outdir}")


if __name__ == "__main__":
    main()
