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

    # Pearson entre los perfiles de la misma muestra y especies compartidas. Solo cuentan las especies
    # presentes en alguna de las dos versiones de esa muestra; si no, los ceros compartidos inflan la r
    pearson, compartidas, solo_yo, solo_art, ab_yo, ab_art = [], [], [], [], [], []
    for s in runs:
        yo, ar = Ur[s] > 0, Ar[s] > 0
        alguna = yo | ar
        pearson.append(pearsonr(Ur[s][alguna], Ar[s][alguna])[0])
        compartidas.append(int((yo & ar).sum()))
        solo_yo.append(int((yo & ~ar).sum()))
        solo_art.append(int((ar & ~yo).sum()))
        # qué fracción de las lecturas de cada versión cae en las especies compartidas
        ab_yo.append(Ur[s][yo & ar].sum())
        ab_art.append(Ar[s][yo & ar].sum())
    res["pearson_misma"] = np.array(pearson)
    res["compartidas"], res["solo_yo"], res["solo_art"] = np.array(compartidas), np.array(solo_yo), np.array(solo_art)
    res["ab_yo"], res["ab_art"] = np.array(ab_yo), np.array(ab_art)

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
    # 2 columnas x 4 filas al ancho útil de una hoja carta (16,5 cm), para que en el documento
    # no haya que reducirla; el título y las explicaciones van en la leyenda de la figura
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
                         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                         "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA,
                         "xtick.color": TINTA2, "ytick.color": TINTA2, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                         "grid.linewidth": 0.5, "axes.axisbelow": True, "axes.linewidth": 0.6})

    def diagonal(ax, lo, hi):
        ax.plot([lo, hi], [lo, hi], ls="--", lw=0.8, color=GRIS, zorder=1)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)

    def nota(ax, texto, donde="arriba_izq"):
        x, y, ha, va = {"arriba_izq": (0.03, 0.97, "left", "top"), "arriba_der": (0.97, 0.97, "right", "top"),
                        "abajo_der": (0.97, 0.03, "right", "bottom")}[donde]
        ax.text(x, y, texto, transform=ax.transAxes, ha=ha, va=va, fontsize=6.5, color=TINTA,
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))

    def puntos(ax, x, y):
        ax.scatter(x, y, s=9, color=YO, ec="white", lw=0.3, zorder=3)

    def titulo(ax, texto):
        ax.set_title(texto, loc="left", fontweight="bold", pad=4)

    fig, axs = plt.subplots(4, 2, figsize=(6.5, 7.6))
    axs = axs.ravel()

    ax = axs[0]
    puntos(ax, m.Bases / 1e6, m.bases_iniciales / 1e6)
    diagonal(ax, 15, 65)
    ax.set_xlabel("Bases según el BioProject (Mb)"); ax.set_ylabel("Bases descargadas (Mb)")
    titulo(ax, "A. Datos iniciales")
    nota(ax, f"{int((m.Bases == m.bases_iniciales).sum())} de {len(m)} muestras idénticas")

    # B, C y D contra las lecturas asignadas por EMU del artículo, que es lo único por muestra que publicaron
    def etapa(ax, col, texto_titulo, ylab, extra, misma_etapa):
        x, y = m.emu_art / 1e3, m[col] / 1e3
        puntos(ax, x, y)
        hi = max(x.max(), y.max()) * 1.08
        k = 1.0
        if misma_etapa:
            diagonal(ax, 0, hi)
        else:
            k = (x * y).sum() / (x * x).sum()
            ax.plot([0, hi], [0, k * hi], color=GRIS, lw=0.9, zorder=1)
            ax.set_xlim(0, hi); ax.set_ylim(0, hi)
            extra = f"pendiente = {k:.2f}\n" + extra
        ax.set_xlabel("Artículo: lecturas asignadas por EMU (miles)"); ax.set_ylabel(ylab)
        titulo(ax, texto_titulo)
        # si la nube queda por encima de la diagonal, la nota va abajo para no tapar puntos
        nota(ax, (f"r = {pearsonr(m.emu_art, m[col])[0]:.3f}\n" + extra).strip(),
             "abajo_der" if k > 1 else "arriba_izq")

    ret = m.lect_limpias / m.lect_iniciales
    etapa(axs[1], "lect_iniciales", "B. Pre-limpieza", "Lecturas iniciales (miles)", "", False)
    etapa(axs[2], "lect_limpias", "C. Post-limpieza", "Lecturas limpias (miles)",
          f"retención: {100 * ret.mean():.0f} % ({100 * ret.min():.0f}–{100 * ret.max():.0f} %)", False)
    etapa(axs[3], "emu_yo", "D. Lecturas asignadas por EMU", "Lecturas asignadas (miles)",
          f"asignadas/iniciales: {100 * (m.emu_yo / m.lect_iniciales).mean():.0f} %\n"
          f"(artículo: {100 * (m.emu_art / m.lect_iniciales).mean():.0f} %)", True)

    ax = axs[4]
    bins = np.linspace(0, 1, 41)
    ax.hist(res["bc_distintas"], bins=bins, density=True, color=GRIS, alpha=0.55, label="Muestras distintas", ec="white", lw=0.3)
    ax.hist(res["bc_misma"], bins=bins, density=True, color=YO, alpha=0.9, label="Misma muestra", ec="white", lw=0.3)
    ax.set_xlabel("Disimilitud de Bray-Curtis (especie)"); ax.set_ylabel("Densidad")
    titulo(ax, "E. Composición")
    ax.legend(frameon=False, loc="center right", bbox_to_anchor=(1, 0.3), fontsize=6.5, handlelength=1.2)
    pct = 100 * res["compartidas"] / (res["compartidas"] + res["solo_yo"] + res["solo_art"])
    # la nota va arriba a la derecha, donde las barras son bajas (todo son medianas por muestra)
    nota(ax, f"Bray-Curtis misma muestra: {np.median(res['bc_misma']):.2f}\n"
             f"Bray-Curtis muestras distintas: {np.median(res['bc_distintas']):.2f}\n"
             f"Pearson misma muestra: r = {np.median(res['pearson_misma']):.3f}\n"
             f"especies compartidas: {np.median(res['compartidas']):.0f} ({np.median(pct):.0f} %)\n"
             f"lecturas en ellas: {100 * np.median(res['ab_yo']):.0f} % / {100 * np.median(res['ab_art']):.0f} %",
         "arriba_der")

    ax = axs[5]
    puntos(ax, res["sh_art"], res["sh_yo"])
    lo = min(res["sh_art"].min(), res["sh_yo"].min()) - 0.1
    hi = max(res["sh_art"].max(), res["sh_yo"].max()) + 0.1
    diagonal(ax, lo, hi)
    ax.set_xlabel("Shannon, artículo"); ax.set_ylabel("Shannon, este trabajo")
    titulo(ax, "F. Diversidad alfa")
    nota(ax, f"Shannon: r = {pearsonr(res['sh_art'], res['sh_yo'])[0]:.3f}\n"
             f"riqueza observada: r = {pearsonr(res['obs_art'], res['obs_yo'])[0]:.3f}")

    ax = axs[6]
    da, du = res["ait_art"].condensed_form(), res["ait_yo"].condensed_form()
    ax.hexbin(da, du, gridsize=35, cmap="Blues", mincnt=1, linewidths=0)
    diagonal(ax, min(da.min(), du.min()), max(da.max(), du.max()))
    ax.set_xlabel("Aitchison, artículo"); ax.set_ylabel("Aitchison, este trabajo")
    titulo(ax, "G. Diversidad beta (pares de muestras)")
    r_ait, p_ait, _ = res["mantel_ait"]
    nota(ax, f"Mantel r = {r_ait:.3f} (p = {p_ait:.3f})\nBray-Curtis: Mantel r = {res['mantel_bc'][0]:.3f}")

    ax = axs[7]
    P = res["permanova"].copy()
    P["etiqueta"] = P.Tipo.map({"bulk": "Suelo", "bag": "Bolsa"}) + " s" + P.Semana.astype(str)
    yy = np.arange(len(P))[::-1]
    for i, (_, r) in enumerate(P.iterrows()):
        ax.plot([r.R2_publicado, r.R2_este_trabajo], [yy[i]] * 2, color="#d4d3cd", lw=1, zorder=1)
    for col, c, mk in (("publicado", ART, "o"), ("este_trabajo", YO, "D")):
        sig = (P[f"p_{col}"] < 0.05).values
        ax.scatter(P[f"R2_{col}"][~sig], yy[~sig], s=16, facecolor="white", edgecolor=c, lw=1, marker=mk, zorder=3)
        ax.scatter(P[f"R2_{col}"][sig], yy[sig], s=16, color=c, marker=mk, zorder=3)
    ax.set_yticks(yy); ax.set_yticklabels(P.etiqueta, fontsize=6); ax.grid(axis="y", visible=False)
    ax.set_ylim(-0.7, len(P) - 0.3)
    ax.set_xlabel("R² del sistema de rotación (PERMANOVA)")
    titulo(ax, "H. Efecto de la rotación")
    iguales = ((P.p_publicado < 0.05) == (P.p_este_trabajo < 0.05)).sum()
    nota(ax, f"misma conclusión en\n{iguales} de {len(P)} subconjuntos", "arriba_der")
    leyenda = [Line2D([], [], marker="o", ls="", color=ART, markersize=4, label="Publicado"),
               Line2D([], [], marker="D", ls="", color=YO, markersize=4, label="Este trabajo"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=TINTA2, markersize=4, label="vacío: p ≥ 0.05")]
    ax.legend(handles=leyenda, frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.42, -0.3),
              ncol=3, handletextpad=0.3, columnspacing=0.8)

    fig.tight_layout(h_pad=1.0, w_pad=1.5)
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
        ["5. Composición", "Pearson entre perfiles de la misma muestra (mediana, mín–máx)",
         f"r = {np.median(res['pearson_misma']):.3f} ({res['pearson_misma'].min():.3f}–{res['pearson_misma'].max():.3f})"],
        ["5. Composición", "Especies compartidas por muestra (mediana)",
         f"{np.median(res['compartidas']):.0f} de {np.median(res['compartidas'] + res['solo_yo'] + res['solo_art']):.0f} "
         f"({np.median(100 * res['compartidas'] / (res['compartidas'] + res['solo_yo'] + res['solo_art'])):.1f} %)"],
        ["5. Composición", "Lecturas en especies compartidas (mediana)",
         f"este trabajo {100 * np.median(res['ab_yo']):.1f} %; artículo {100 * np.median(res['ab_art']):.1f} %"],
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
        "Riqueza (este trabajo)": m.Run.map(res["obs_yo"]), "Riqueza (artículo)": m.Run.map(res["obs_art"]),
        "Bray-Curtis misma muestra": res["bc_misma"], "Pearson misma muestra": res["pearson_misma"],
        "Especies compartidas": res["compartidas"], "Especies solo en este trabajo": res["solo_yo"],
        "Especies solo en el artículo": res["solo_art"],
        "% lecturas en compartidas (este trabajo)": 100 * res["ab_yo"],
        "% lecturas en compartidas (artículo)": 100 * res["ab_art"],
    })
    # medias por tipo de muestra (bolsa o suelo)
    cols = ["Bray-Curtis misma muestra", "Pearson misma muestra", "Especies compartidas",
            "Shannon (este trabajo)", "Shannon (artículo)", "Riqueza (este trabajo)", "Riqueza (artículo)"]
    tipo = por_muestra["Tipo"].map({"bag": "Bolsa", "bulk": "Suelo"})
    medias = por_muestra[cols].groupby(tipo).agg(["mean", "std"]).round(3)
    medias.columns = [f"{c} ({'media' if e == 'mean' else 'desv. est.'})" for c, e in medias.columns]
    medias = medias.reset_index().rename(columns={"Tipo": "Tipo de muestra"})
    perm = P.rename(columns={"R2_publicado": "R² publicado", "p_publicado": "p publicado",
                             "R2_este_trabajo": "R² este trabajo", "p_este_trabajo": "p este trabajo",
                             "R2_tabla_articulo": "R² tabla del artículo (recalculado)",
                             "p_tabla_articulo": "p tabla del artículo (recalculado)"})
    with pd.ExcelWriter(salida) as w:
        resumen.to_excel(w, sheet_name="Resumen", index=False)
        por_muestra.to_excel(w, sheet_name="Por_muestra", index=False)
        medias.to_excel(w, sheet_name="Medias_por_tipo", index=False)
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
