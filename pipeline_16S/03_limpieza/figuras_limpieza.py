#!/usr/bin/env python3

# Figuras 10 y 11 de la tesis: lecturas y calidad antes y después de la limpieza (PRJNA1020132).
# La 10 y el panel B de la 11 salen de las tablas de estadísticas (Anexos 3 y 4). El panel A de
# la 11 necesita el QScore de cada lectura, así que solo se hace si se pasan las carpetas de FASTQ.
# Ejemplo:
#   python figuras_limpieza.py --pre Anexo_03_Estadisticas_pre_limpieza.xlsx \
#     --post Anexo_04_Estadisticas_post_limpieza.xlsx \
#     --fastq-antes Muestras_16S --fastq-despues Limpieza -o figuras

import argparse
import gzip
import random
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator

AZUL, NARANJA, GRIS = "#0072B2", "#E69F00", "#9a9a9a"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e5e5e5", "grid.linewidth": 0.8, "axes.axisbelow": True})


# números con punto de miles y coma decimal, como en el texto de la tesis
def miles(x, _=None):
    return f"{x:,.0f}".replace(",", ".")


def decimal(n):
    return FuncFormatter(lambda x, _=None: f"{x:.{n}f}".replace(".", ","))


def coma(x, n=1):
    return f"{x:.{n}f}".replace(".", ",")


def leer_estadisticas(ruta: Path) -> pd.DataFrame:
    # sirve con los excel de los anexos (encabezado corrido) o con la salida de stats_fastq.py
    crudo = pd.read_excel(ruta, header=None) if ruta.suffix.lower() == ".xlsx" else pd.read_csv(ruta, sep="\t", header=None)
    fila = next(i for i, r in crudo.iterrows() if r.astype(str).str.strip().eq("Muestra").any())
    df = crudo.iloc[fila + 1:].copy()
    df.columns = crudo.iloc[fila].astype(str).str.strip()
    df["Run"] = df["Muestra"].astype(str).str.extract(r"(SRR\d+)", expand=False)
    df = df[df.Run.notna()].set_index("Run")
    return df[["lecturas", "QScore promedio"]].apply(pd.to_numeric)


def figura_10(pre, post, salida: Path):
    ant, des = pre.lecturas, post.lecturas.reindex(pre.index)
    ret = 100 * des / ant
    med = ret.median()
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.8))

    tope = max(ant.max(), des.max()) * 1.07
    a.plot([0, tope], [0, tope], ls="--", color=GRIS, lw=1.3, label="Sin pérdida (100 %)")
    a.plot([0, tope], [0, tope * med / 100], color=NARANJA, lw=1.8, label=f"Retención mediana ({coma(med)} %)")
    a.scatter(ant, des, s=34, color=AZUL, ec="white", lw=0.7, zorder=3, label=f"Muestras (n = {len(ant)})")
    a.set_xlim(0, tope); a.set_ylim(0, tope)
    a.xaxis.set_major_formatter(FuncFormatter(miles)); a.yaxis.set_major_formatter(FuncFormatter(miles))
    a.set_xlabel("Lecturas antes del filtrado"); a.set_ylabel("Lecturas después del filtrado")
    a.legend(frameon=False, loc="upper left")

    bordes = np.arange(np.floor(ret.min() * 2) / 2, np.ceil(ret.max() * 2) / 2 + 0.5, 0.5)
    conteo, _, _ = b.hist(ret, bins=bordes, color=AZUL, alpha=0.85, ec="white", lw=1)
    b.plot([med, med], [0, conteo.max() * 1.08], color=NARANJA, lw=2)
    b.set_xlim(bordes[0] - 0.8, bordes[-1] + 0.8)
    b.xaxis.set_major_locator(MultipleLocator(1))
    b.yaxis.set_major_locator(MultipleLocator(3)); b.set_ylim(0, conteo.max() * 1.5)
    b.set_xlabel("Lecturas retenidas por muestra (%)"); b.set_ylabel("Número de muestras")
    b.text(0.03, 0.96, f"Mediana: {coma(med)} %\nRango: {coma(ret.min())}–{coma(ret.max())} %\n"
                       f"Total: {miles(ant.sum())} → {miles(des.sum())} lecturas",
           transform=b.transAxes, va="top", bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))

    for ax, letra in ((a, "A"), (b, "B")):
        ax.text(-0.1, 1.04, letra, transform=ax.transAxes, fontsize=14, fontweight="bold")
    fig.tight_layout(); fig.savefig(salida, dpi=300); plt.close(fig)


def caja_y_puntos(ax, pre, post):
    antes, despues = pre["QScore promedio"], post["QScore promedio"].reindex(pre.index)
    rng = np.random.default_rng(42)
    xa = 0.17 + rng.uniform(-0.04, 0.04, len(antes))
    xd = 0.83 + rng.uniform(-0.04, 0.04, len(despues))
    for x0, y0, x1, y1 in zip(xa, antes, xd, despues):
        ax.plot([x0, x1], [y0, y1], color="#c8c8c8", lw=0.7, zorder=1)
    ax.scatter(xa, antes, s=22, color=AZUL, ec="white", lw=0.4, zorder=3)
    ax.scatter(xd, despues, s=22, color=NARANJA, ec="white", lw=0.4, zorder=3)
    for pos, datos, color in ((0.0, antes, AZUL), (1.0, despues, NARANJA)):
        bp = ax.boxplot(datos, positions=[pos], widths=0.12, patch_artist=True, showfliers=False)
        bp["boxes"][0].set(facecolor=color, alpha=0.55)
        bp["medians"][0].set(color="black", lw=1.5)
    ax.set_xticks([0.17, 0.83], ["Antes", "Después"]); ax.set_xlim(-0.15, 1.15)
    ax.yaxis.set_major_formatter(decimal(1)); ax.grid(axis="x", visible=False)
    ax.set_ylabel("QScore promedio por muestra")
    cambio = (despues - antes).median()
    ax.set_title(f"n = {len(antes)} muestras · cambio mediano: +{coma(cambio, 2)}", fontsize=10, color="#444444")


# QScore medio de cada lectura (promedio de los Phred, igual que stats_fastq.py)
def qscores_por_lectura(ruta: Path, n: int, semilla: int) -> np.ndarray:
    abrir = gzip.open if ruta.suffix == ".gz" else open
    elegidas, i = [], 0
    rng = random.Random(semilla)
    with abrir(ruta, "rt") as fh:
        while True:
            if not fh.readline():
                break
            fh.readline(); fh.readline()
            cal = fh.readline().strip()
            q = (np.frombuffer(cal.encode(), dtype=np.uint8) - 33).mean()
            # muestreo de reservorio: n lecturas al azar sin cargar el archivo completo
            if i < n:
                elegidas.append(q)
            else:
                j = rng.randint(0, i)
                if j < n:
                    elegidas[j] = q
            i += 1
    return np.array(elegidas)


def buscar_fastq(carpeta: Path, runs) -> dict:
    archivos = {}
    for f in sorted(carpeta.rglob("*")):
        if re.search(r"\.(fastq|fq)(\.gz)?$", f.name):
            m = re.search(r"(SRR\d+)", f.name)
            if m and m.group(1) in runs and m.group(1) not in archivos:
                archivos[m.group(1)] = f
    return archivos


def figura_11(pre, post, salida: Path, fastq_antes=None, fastq_despues=None, total=120000):
    if fastq_antes and fastq_despues:
        fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.8))
        # hasta 120.000 lecturas por grupo, repartidas por igual entre las muestras
        por_muestra = total // len(pre)
        grupos = {}
        for nombre, carpeta in (("Antes", fastq_antes), ("Después", fastq_despues)):
            archivos = buscar_fastq(carpeta, set(pre.index))
            faltan = sorted(set(pre.index) - set(archivos))
            if faltan:
                print(f"[!] {nombre}: no encontré FASTQ de {len(faltan)} muestras: {faltan[:5]}")
            grupos[nombre] = np.concatenate([qscores_por_lectura(f, por_muestra, k) for k, f in enumerate(archivos.values())])
        bordes = np.arange(8, 26.25, 0.25)
        for nombre, color in (("Antes", AZUL), ("Después", NARANJA)):
            q = grupos[nombre]
            peso = np.full(len(q), 100 / len(q))
            mediana = np.median(q)
            a.hist(q, bins=bordes, weights=peso, histtype="stepfilled", color=color, alpha=0.15)
            a.hist(q, bins=bordes, weights=peso, histtype="step", color=color, lw=1.8,
                   label=f"{nombre} (mediana {coma(mediana)})")
            a.axvline(mediana, color=color, ls="--", lw=1.4)
        a.set_xlim(7, 27); a.xaxis.set_major_locator(MultipleLocator(2))
        a.yaxis.set_major_formatter(decimal(1))
        a.set_xlabel("QScore medio por lectura"); a.set_ylabel("Lecturas (%)")
        a.legend(frameon=False, loc="upper left")
        caja_y_puntos(b, pre, post)
        for ax, letra in ((a, "A"), (b, "B")):
            ax.text(-0.1, 1.04, letra, transform=ax.transAxes, fontsize=14, fontweight="bold")
    else:
        # sin FASTQ solo se puede hacer el panel B
        fig, b = plt.subplots(figsize=(5.75, 4.8))
        caja_y_puntos(b, pre, post)
        b.text(-0.12, 1.04, "B", transform=b.transAxes, fontsize=14, fontweight="bold")
    fig.tight_layout(); fig.savefig(salida, dpi=300); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Figuras 10 y 11: lecturas y calidad antes y después de la limpieza.")
    ap.add_argument("--pre", type=Path, required=True, help="Estadísticas antes de la limpieza (Anexo 3).")
    ap.add_argument("--post", type=Path, required=True, help="Estadísticas después de la limpieza (Anexo 4).")
    ap.add_argument("--fastq-antes", type=Path, default=None, help="Carpeta con los FASTQ descargados (para el panel 11A).")
    ap.add_argument("--fastq-despues", type=Path, default=None, help="Carpeta con los FASTQ limpios (para el panel 11A).")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("figuras_limpieza"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    pre, post = leer_estadisticas(args.pre), leer_estadisticas(args.post)
    comunes = pre.index.intersection(post.index)
    pre, post = pre.loc[comunes], post.loc[comunes]
    print(f"[OK] Muestras: {len(pre)}")

    figura_10(pre, post, args.outdir / "figura_10_retencion.png")
    nombre_11 = "figura_11_calidad.png" if args.fastq_antes and args.fastq_despues else "figura_11B_calidad_por_muestra.png"
    figura_11(pre, post, args.outdir / nombre_11, args.fastq_antes, args.fastq_despues)
    print(f"[OK] Resultados en: {args.outdir}")


if __name__ == "__main__":
    main()
