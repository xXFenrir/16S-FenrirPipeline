#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Comparación gráfica ANTES vs. DESPUÉS (FASTQ por muestra + per-read).

Por defecto (--estilo nuevo) genera dos figuras con TODAS las muestras:
- 01_lecturas_retencion.png: (A) lecturas antes vs. después por muestra; (B) % de retención por muestra.
- 02_calidad_qscore.png: (A) distribución del QScore medio por lectura (requiere FASTQ);
  (B) QScore promedio por muestra antes y después (sale de las tablas de stats_fastq.py).
Con --estilo anterior (o ambos) genera además las figuras originales:
- Mancuernas (# lecturas y # bases) con %remoción por muestra (+ submuestreo).
  * Matplotlib: pares correctos (línea por muestra), colores/forma distintos (Antes vs. Después), etiquetas rem fuera.
  * Plotly (--use-plotly): HTML interactivo con segmentos por muestra (None separator) + PNG si hay 'kaleido'.
- Longitudes por lectura: Violín estándar o Raincloud (--raincloud) y opción log10 (--lengths-log10).
- N50 por muestra: Violín estándar o Raincloud.
- Scatter de QScore por lectura (no promedios), con grupos cercanos y jitter suave.
- CSV de retención/remoción.
- PDF mosaico de una sola página con las figuras PNG/JPG.
- Ruta de resultados configurable con --outdir o alias --results-dir.
"""

from __future__ import annotations
import argparse
import gzip
import re
import sys
from glob import glob
from pathlib import Path
from typing import Iterable, Tuple, List, Optional
from fnmatch import fnmatch

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.image import imread

# ====== Opcionales (se usan si están disponibles) ======
# Seaborn
try:
    import seaborn as sns
    _HAS_SNS = True
    sns.set_context("talk")
    sns.set_style("whitegrid")
except Exception:
    _HAS_SNS = False

# Plotly
try:
    import plotly.graph_objects as go
    _HAS_PLOTLY = True
except Exception:
    _HAS_PLOTLY = False

# Raincloud (ptitprince)
try:
    import ptitprince as pt
    _HAS_PTIP = True
except Exception:
    _HAS_PTIP = False


# ------------------------
# Columnas esperadas
# ------------------------

HEADERS_ES = [
    "Muestra","lecturas","bases","longitud promedio","longitud mínima","longitud máxima",
    "N50","GC%","QScore promedio"
]

NUM_COLS = {
    "lecturas": float,
    "bases": float,
    "longitud promedio": float,
    "longitud mínima": float,
    "longitud máxima": float,
    "N50": float,
    "GC%": float,
    "QScore promedio": float,
}

# Mínimos para empatar/graficar lo requerido
REQUIRED_FOR_PLOTS = [
    "Muestra","lecturas","bases","longitud promedio","N50","QScore promedio"
]


# ------------------------
# Lectura de tablas
# ------------------------

def _parse_pretty_txt(path: Path) -> pd.DataFrame:
    """Lee un TXT 'bonito' con barras verticales."""
    rows = []
    headers = None
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.rstrip('\n')
            if not s or set(s) == {'-'}:
                continue
            t = s.strip()
            if t.startswith("|") and t.endswith("|"):
                parts = [c.strip() for c in t[1:-1].split("|")]
                if headers is None:
                    headers = parts
                else:
                    rows.append(parts)
    if headers is None:
        raise ValueError(f"No se encontraron cabeceras en: {path}")
    df = pd.DataFrame(rows, columns=headers)
    for col in NUM_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col].replace({"": np.nan}), errors="coerce")
    return df


def load_table(path: Path) -> pd.DataFrame:
    ext = path.suffix.lower()
    if path.is_dir():
        raise FileNotFoundError(f"Se recibió un directorio en lugar de archivo: {path}")
    if ext == ".xlsx":
        df = pd.read_excel(path)
    elif ext in (".tsv", ".tab"):
        df = pd.read_csv(path, sep="\t")
    elif ext == ".txt":
        df = _parse_pretty_txt(path)
    else:
        # Intentos flexibles
        try:
            df = pd.read_excel(path)
        except Exception:
            try:
                df = pd.read_csv(path, sep="\t")
            except Exception:
                df = _parse_pretty_txt(path)
    if "Muestra" not in df.columns:
        raise ValueError(f"Falta columna 'Muestra' en {path}")
    for col in NUM_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# ------------------------
# Normalización y unión
# ------------------------

def normalize_sample_name(s: str) -> str:
    """Normaliza nombre de muestra para empatar 'antes'/'después'."""
    if s is None:
        return ""
    base = str(s).strip().split("/")[-1].lower()
    base = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
    suffixes = ("clean","filtered","trim","trimmed","dada2","denoise","dentrim","oriented","reoriented","final")
    changed = True
    while changed:
        changed = False
        for sep in ("_", "-"):
            for suf in suffixes:
                pat = rf"{sep}{suf}$"
                if re.search(pat, base, flags=re.IGNORECASE):
                    base = re.sub(pat, "", base, flags=re.IGNORECASE)
                    changed = True
    return base.strip()


def join_before_after(df_before: pd.DataFrame, df_after: pd.DataFrame) -> pd.DataFrame:
    b = df_before.copy()
    a = df_after.copy()
    b["__key__"] = b["Muestra"].apply(normalize_sample_name)
    a["__key__"] = a["Muestra"].apply(normalize_sample_name)

    set_b = set(b["__key__"]); set_a = set(a["__key__"])
    inter = set_b & set_a
    print(f"[INFO] Muestras ANTES: {len(set_b)} | DESPUÉS: {len(set_a)} | Intersección: {len(inter)}", file=sys.stderr)

    only_b = sorted(set_b - set_a)
    only_a = sorted(set_a - set_b)
    if only_b:
        print(f"[WARN] Solo en ANTES (tras normalizar): {only_b[:10]}{' ...' if len(only_b)>10 else ''}", file=sys.stderr)
    if only_a:
        print(f"[WARN] Solo en DESPUÉS (tras normalizar): {only_a[:10]}{' ...' if len(only_a)>10 else ''}", file=sys.stderr)

    df = b.merge(a, on="__key__", how="inner", suffixes=("_before", "_after"))
    if df.empty:
        print("[ERROR] Tras normalizar, no hay coincidencias de muestras.", file=sys.stderr)
        return df

    same = (df["Muestra_before"].str.lower() == df["Muestra_after"].str.lower())
    df["Muestra"] = np.where(same, df["Muestra_before"], df["Muestra_after"])

    # métricas de retención/remoción
    for m in ("lecturas","bases"):
        cb, ca = f"{m}_before", f"{m}_after"
        if cb in df.columns and ca in df.columns:
            with np.errstate(divide="ignore", invalid="ignore"):
                df[f"{m}_retencion_%"] = 100.0 * df[ca] / df[cb]
                df[f"{m}_remocion_%"] = 100.0 - df[f"{m}_retencion_%"]

    cols = ["Muestra"] + [c for c in df.columns if c not in {"Muestra","Muestra_before","Muestra_after","__key__"}]
    df = df[cols]
    return df


# ------------------------
# Helpers de guardado
# ------------------------

def ensure_outdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def savefig(fname: Path):
    plt.tight_layout()
    plt.savefig(fname, dpi=220)
    plt.close()


# ------------------------
# Subconjunto de muestras para gráficos
# ------------------------

def choose_subset(df: pd.DataFrame, metric: str, max_samples: int, subset_mode: str) -> pd.DataFrame:
    """Devuelve un subconjunto ordenado para gráficos legibles."""
    if max_samples is None or max_samples <= 0 or len(df) <= max_samples:
        return df
    cb, ca = f"{metric}_before", f"{metric}_after"
    if subset_mode == "highest_removal" and f"{metric}_remocion_%" in df.columns:
        take = df.dropna(subset=[f"{metric}_remocion_%"]).sort_values(by=f"{metric}_remocion_%", ascending=False).head(max_samples)
    elif subset_mode == "largest_before" and cb in df.columns:
        take = df.dropna(subset=[cb]).sort_values(by=cb, ascending=False).head(max_samples)
    elif subset_mode == "largest_after" and ca in df.columns:
        take = df.dropna(subset=[ca]).sort_values(by=ca, ascending=False).head(max_samples)
    else:
        take = df.sample(n=max_samples, random_state=7)
    return take


# ------------------------
# 1) Dumbbell (lecturas / bases) — Matplotlib (pareado + colores)
# ------------------------

def dumbbell_plot(df: pd.DataFrame, metric: str, outpath: Path,
                  max_samples: Optional[int] = None, subset_mode: str = "highest_removal",
                  annotate_top: int = 0, fig_scale: float = 1.0, lineheight: float = 0.38):
    """
    Mancuerna ANTES vs DESPUÉS con pares por muestra.
    Antes: azul / círculo; Después: ámbar / cuadrado.
    """
    import matplotlib.lines as mlines

    req = [f"{metric}_before", f"{metric}_after", "Muestra"]
    if any(r not in df.columns for r in req):
        print(f"[WARN] Faltan columnas para dumbbell {metric}.", file=sys.stderr); return

    data = df.dropna(subset=[f"{metric}_before", f"{metric}_after"]).copy()
    data = choose_subset(data, metric, max_samples, subset_mode)
    if data.empty:
        print(f"[WARN] No hay datos válidos para dumbbell {metric}.", file=sys.stderr); return

    order = data.sort_values(by=f"{metric}_after", ascending=True).reset_index(drop=True)
    y = np.arange(len(order))

    with np.errstate(divide="ignore", invalid="ignore"):
        ret = 100.0 * order[f"{metric}_after"] / order[f"{metric}_before"]
        rem = 100.0 - ret

    xs = np.concatenate([order[f"{metric}_before"].values, order[f"{metric}_after"].values])
    xs = xs[np.isfinite(xs)]
    x_min, x_max = (np.nanmin(xs), np.nanmax(xs)) if xs.size else (0, 1)
    pad = 0.18 * (x_max - x_min + 1e-9)  # espacio para etiquetas a la derecha

    h = max(5.0*fig_scale, len(order) * lineheight * fig_scale)
    w = 12.5 * fig_scale
    plt.figure(figsize=(w, h))
    ax = plt.gca()

    # top-N anotado por %remoción
    idxs = np.argsort(rem.values)[::-1] if np.isfinite(rem).any() else np.arange(len(order))
    annotate_mask = np.zeros(len(order), dtype=bool)
    if annotate_top and annotate_top > 0:
        annotate_mask[idxs[:min(annotate_top, len(order))]] = True
    else:
        annotate_mask[:] = True

    # Colores/estética
    c_before = "#1f77b4"  # azul
    c_after  = "#ff7f0e"  # ámbar

    for i, row in order.iterrows():
        xb = row[f"{metric}_before"]; xa = row[f"{metric}_after"]
        # línea por muestra
        ax.plot([xb, xa], [i, i], color="0.45", linewidth=2, zorder=1)
        # puntos
        ax.scatter([xb], [i], color=c_before, s=28, zorder=3, marker='o')
        ax.scatter([xa], [i], color=c_after,  s=28, zorder=3, marker='s')
        # etiqueta de remoción fuera de la línea
        if annotate_mask[i] and np.isfinite(rem.iloc[i]):
            x_anchor = max(xb, xa)
            dy = 0.10 if (i % 2 == 0) else -0.10
            ax.text(x_anchor + 0.02*(x_max-x_min), i+dy, f"rem={rem.iloc[i]:.1f}%",
                    va="center", ha="left", fontsize=9)

    ax.set_yticks(y)
    ax.set_yticklabels(order["Muestra"], fontsize=9)
    ax.set_xlabel(metric)
    ax.set_title(f"Dumbbell: {metric} (antes vs. después)", fontsize=13, pad=10)
    ax.grid(True, axis='x', linewidth=0.5, alpha=0.5)
    ax.set_xlim(x_min, x_max + pad)

    # Leyenda clara
    import matplotlib.lines as mlines
    h_before = mlines.Line2D([], [], color=c_before, marker='o', linestyle='None', markersize=6, label='Antes')
    h_after  = mlines.Line2D([], [], color=c_after,  marker='s', linestyle='None', markersize=6, label='Después')
    ax.legend(handles=[h_before, h_after], loc='best', fontsize=9, frameon=True)

    plt.tight_layout()
    plt.savefig(outpath, dpi=220)
    plt.close()


# ------------------------
# 1b) Dumbbell Plotly (interactivo pareado + colores)
# ------------------------

def dumbbell_plot_plotly(df: pd.DataFrame, metric: str, out_html: Path, out_png: Path | None = None,
                         max_samples: Optional[int] = None, subset_mode: str = "highest_removal",
                         annotate_top: int = 0):
    """
    Mancuerna interactiva correcta:
    - segmentos por muestra con separadores None (no zig-zag)
    - colores/forma distintos (Antes=azul/◯, Después=ámbar/■)
    - etiquetas rem=% a la derecha
    """
    if not _HAS_PLOTLY:
        print("[WARN] Plotly no está disponible. Instala 'plotly' y 'kaleido' para PNG.", file=sys.stderr)
        return

    req = [f"{metric}_before", f"{metric}_after", "Muestra"]
    if any(r not in df.columns for r in req):
        print(f"[WARN] Faltan columnas para dumbbell {metric}.", file=sys.stderr)
        return

    data = df.dropna(subset=[f"{metric}_before", f"{metric}_after"]).copy()
    data = choose_subset(data, metric, max_samples, subset_mode)
    if data.empty:
        print(f"[WARN] No hay datos válidos para dumbbell {metric}.", file=sys.stderr)
        return

    order = data.sort_values(by=f"{metric}_after", ascending=True).reset_index(drop=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ret = 100.0 * order[f"{metric}_after"] / order[f"{metric}_before"]
        rem = 100.0 - ret

    # Top-N anotado
    annotate_mask = np.zeros(len(order), dtype=bool)
    if annotate_top and np.isfinite(rem).any():
        idxs = np.argsort(rem.values)[::-1][:min(annotate_top, len(order))]
        annotate_mask[idxs] = True
    else:
        annotate_mask[:] = True

    ylabels = order["Muestra"].tolist()
    y = np.arange(len(order))

    # 1) Segmentos por muestra usando None como separador
    seg_x = []
    seg_y = []
    for i, row in order.iterrows():
        xb = row[f"{metric}_before"]; xa = row[f"{metric}_after"]
        seg_x.extend([xb, xa, None])
        seg_y.extend([i,  i,  None])

    fig = go.Figure()
    fig.add_scatter(
        x=seg_x, y=seg_y, mode="lines",
        line=dict(color="rgba(80,80,80,0.65)", width=2),
        hoverinfo="skip", showlegend=False
    )

    # 2) Puntos “Antes”
    fig.add_scatter(
        x=order[f"{metric}_before"], y=y,
        mode="markers", name="Antes",
        marker=dict(size=8, color="#1f77b4", symbol="circle"),
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            f"{metric} antes: %{{x}}<br>"
            f"{metric} después: %{{customdata[1]}}<br>"
            "rem=%{customdata[2]:.1f}%<extra></extra>"
        ),
        customdata=np.c_[order["Muestra"], order[f"{metric}_after"], rem]
    )

    # 3) Puntos “Después”
    fig.add_scatter(
        x=order[f"{metric}_after"], y=y,
        mode="markers", name="Después",
        marker=dict(size=8, color="#ff7f0e", symbol="square"),
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            f"{metric} después: %{{x}}<br>"
            f"{metric} antes: %{{customdata[1]}}<br>"
            "rem=%{customdata[2]:.1f}%<extra></extra>"
        ),
        customdata=np.c_[order["Muestra"], order[f"{metric}_before"], rem]
    )

    # 4) Etiquetas rem=% a la derecha
    x_max = np.nanmax([order[f"{metric}_before"].values, order[f"{metric}_after"].values])
    pad = 0.02 * (x_max + 1e-9)
    annot_x = np.maximum(order[f"{metric}_before"].values, order[f"{metric}_after"].values) + pad
    annot_y = y + np.where((np.arange(len(y)) % 2) == 0, 0.12, -0.12)
    labels = [f"rem={v:.1f}%" if (np.isfinite(v) and annotate_mask[i]) else "" for i, v in enumerate(rem.values)]
    fig.add_scatter(
        x=annot_x, y=annot_y, mode="text", text=labels,
        textposition="middle left", textfont=dict(size=11),
        showlegend=False, hoverinfo="skip"
    )

    fig.update_yaxes(tickmode="array", tickvals=y, ticktext=ylabels, automargin=True)
    fig.update_xaxes(showgrid=True, gridcolor="rgba(0,0,0,0.2)", zeroline=False)
    fig.update_layout(
        title=f"Dumbbell: {metric} (antes vs. después)",
        margin=dict(l=10, r=20, t=50, b=10),
        height=max(500, int(30*len(y))),
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )

    out_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out_html))
    if out_png is not None:
        try:
            # requiere kaleido: pip install kaleido
            fig.write_image(str(out_png), scale=2)
        except Exception:
            print("[INFO] Para exportar PNG con Plotly instala 'kaleido'. Se generó el HTML.", file=sys.stderr)


# ------------------------
# 2) Lectura FASTQ per-read (longitudes y QScore)
# ------------------------

def _iter_fastq_records(path: str) -> Iterable[Tuple[str, str]]:
    """Genera (seq, qual) por lectura. Phred+33; soporta .gz."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', errors='ignore') as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            fh.readline()  # plus
            qual = fh.readline().rstrip('\n')
            if not seq or not qual:
                continue
            yield (seq, qual)


def _expand_fastq(pat: str) -> list[str]:
    """
    Devuelve una lista de FASTQ (*.fastq, *.fq, *.fastq.gz, *.fq.gz).
    - Si 'pat' es carpeta: búsqueda RECURSIVA en subdirectorios.
    - Si 'pat' es patrón glob: usa glob(recursive=True).
    """
    p = Path(pat)
    exts = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    if p.is_dir():
        files = []
        for e in exts:
            files += [str(x) for x in p.rglob(e)]
        return sorted(files)
    return sorted(glob(pat, recursive=True))


def _qmean_read(qual: str, modo: str = "aritmetica") -> float:
    """QScore medio de una lectura.
    - 'aritmetica': promedio simple de los Phred (mismo criterio que stats_fastq.py y las tablas de la tesis).
    - 'probabilidad': -10*log10(promedio de las probabilidades de error), convención de Dorado/NanoPlot.
    """
    if not qual:
        return np.nan
    qs = np.frombuffer(qual.encode("ascii", "ignore"), dtype=np.uint8).astype(float) - 33.0
    if modo == "probabilidad":
        return float(-10.0 * np.log10(np.mean(10.0 ** (-qs / 10.0))))
    return float(qs.mean())


def sample_lengths_qscores(paths: List[str], reads_cap: int, reads_per_file_cap: int, seed: int = 13,
                           q_modo: str = "aritmetica", balanceado: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Lee longitudes y QScore medio por lectura.
    Con balanceado=True, el tope total se reparte por igual entre archivos, de modo que todas las
    muestras aporten lecturas (antes, con 120.000 de tope y 20.000 por archivo, solo entraban las
    primeras 6 muestras en orden alfabético)."""
    lengths = []
    qmeans = []
    total = 0
    rng = np.random.default_rng(seed)
    if balanceado and paths:
        reads_per_file_cap = max(1, min(reads_per_file_cap, int(np.ceil(reads_cap / len(paths)))))
    for p in paths:
        if total >= reads_cap:
            break
        taken = 0
        try:
            for seq, qual in _iter_fastq_records(p):
                L = len(seq)
                qmean = _qmean_read(qual, q_modo)
                lengths.append(L)
                qmeans.append(qmean)
                taken += 1
                total += 1
                if taken >= reads_per_file_cap or total >= reads_cap:
                    break
        except Exception as e:
            print(f"[WARN] No se pudo leer {p}: {e}", file=sys.stderr)
    return np.array(lengths, dtype=float), np.array(qmeans, dtype=float)


# ------------------------
# 3) Longitudes y QScore per-read (Violín/Raincloud + Scatter)
# ------------------------

def violin_lengths(before_lengths: np.ndarray, after_lengths: np.ndarray, outpath: Path,
                   lengths_log10: bool = False, fig_scale: float = 1.0, use_raincloud: bool = False):
    a = before_lengths[np.isfinite(before_lengths)]
    b = after_lengths[np.isfinite(after_lengths)]
    if a.size==0 and b.size==0:
        print("[WARN] Violín/Raincloud longitudes vacío.", file=sys.stderr); return

    if lengths_log10:
        a = np.log10(a[a>0]); b = np.log10(b[b>0])
        ylab = "log10(Longitud por lectura)"
    else:
        ylab = "Longitud por lectura (nt)"

    df_long = pd.DataFrame({"valor": np.concatenate([a, b]),
                            "estado": np.array(["Antes"]*a.size + ["Después"]*b.size)})

    # Tamaño base más grande
    plt.figure(figsize=(12.0*fig_scale, 7.2*fig_scale))
    if use_raincloud and _HAS_PTIP:
        ax = pt.RainCloud(x="estado", y="valor", data=df_long,
                          width_viol=0.8, width_box=0.18, point_size=2.0,
                          move=0.18, orient="v", bw=0.3, alpha=0.75)
    elif _HAS_SNS:
        ax = sns.violinplot(data=df_long, x="estado", y="valor", cut=0, inner="quartile")
        sns.stripplot(data=df_long.sample(min(4000, len(df_long)), random_state=7),
                      x="estado", y="valor", size=2.2, alpha=0.22, color="k")
    else:
        plt.violinplot(dataset=[a,b], showmeans=True, showextrema=True, showmedians=True)
        plt.xticks([1,2], ["Antes","Después"])
        ax = plt.gca()

    ax.set_ylabel(ylab)
    ax.set_xlabel("")
    ax.set_title(("Raincloud: " if (use_raincloud and _HAS_PTIP) else "Diagrama de violín: ")
                  + "longitudes de lecturas (Antes vs. Después)", pad=12)
    ax.tick_params(axis='both', labelsize=11)
    plt.tight_layout()
    plt.savefig(outpath, dpi=240)
    plt.close()


def scatter_qscores(before_q: np.ndarray, after_q: np.ndarray, outpath: Path, fig_scale: float = 1.0):
    a = before_q[np.isfinite(before_q)]
    b = after_q[np.isfinite(after_q)]
    if a.size==0 and b.size==0:
        print("[WARN] Scatter QScore vacío.", file=sys.stderr); return

    # Tamaño base más grande
    plt.figure(figsize=(11.0*fig_scale, 7.0*fig_scale))
    # grupos cercanos: 1.00 y 1.15
    jitter_before = (np.random.rand(a.size)-0.5)*0.08 + 1.00
    jitter_after  = (np.random.rand(b.size)-0.5)*0.08 + 1.15
    plt.scatter(jitter_before, a, s=8, alpha=0.35, edgecolors='none', label="Antes", color="#1f77b4")
    plt.scatter(jitter_after,  b, s=8, alpha=0.35, edgecolors='none', label="Después", color="#ff7f0e")
    plt.xticks([1.00, 1.15], ["Antes","Después"])
    plt.xlim(0.88, 1.27)
    plt.ylabel("QScore por lectura (media por lectura)")
    plt.title("Scatter: QScore por lectura (Antes vs. Después)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=240)
    plt.close()


def violin_n50(df: pd.DataFrame, outpath: Path, fig_scale: float = 1.0, use_raincloud: bool = False):
    a = df.get("N50_before"); b = df.get("N50_after")
    if a is None or b is None:
        print("[WARN] Faltan columnas N50 para violín/raincloud.", file=sys.stderr); return
    a = pd.Series(a).dropna().values; b = pd.Series(b).dropna().values
    if len(a)==0 and len(b)==0:
        print("[WARN] Violín/Raincloud N50 vacío.", file=sys.stderr); return

    df_long = pd.DataFrame({"N50": np.concatenate([a, b]),
                            "estado": np.array(["Antes"]*len(a) + ["Después"]*len(b))})

    # Tamaño base más grande
    plt.figure(figsize=(12.0*fig_scale, 7.0*fig_scale))
    if use_raincloud and _HAS_PTIP:
        ax = pt.RainCloud(x="estado", y="N50", data=df_long,
                          width_viol=0.8, width_box=0.18, point_size=2.0,
                          move=0.18, orient="v", bw=0.3, alpha=0.75)
    elif _HAS_SNS:
        ax = sns.violinplot(data=df_long, x="estado", y="N50", cut=0, inner="quartile")
        sns.stripplot(data=df_long, x="estado", y="N50", size=2.2, alpha=0.22, color="k")
    else:
        plt.violinplot(dataset=[a,b], showmeans=True, showextrema=True, showmedians=True)
        plt.xticks([1,2], ["Antes","Después"])
        ax = plt.gca()

    ax.set_ylabel("N50 (nt)")
    ax.set_xlabel("")
    ax.set_title(("Raincloud: " if (use_raincloud and _HAS_PTIP) else "Diagrama de violín: ")
                  + "N50 por muestra (Antes vs. Después)", pad=12)
    ax.tick_params(axis='both', labelsize=11)
    plt.tight_layout()
    plt.savefig(outpath, dpi=240)
    plt.close()


# ------------------------
# 4) Figuras nuevas (--estilo nuevo): lecturas y calidad, con todas las muestras
# ------------------------

C_ANTES = "#0072B2"    # azul (Okabe-Ito, apto para daltonismo)
C_DESPUES = "#E69F00"  # naranja (Okabe-Ito)


def _num_es(x: float, dec: int = 0) -> str:
    """Número con formato español: miles con punto y decimales con coma."""
    s = f"{x:,.{dec}f}"
    return s.replace(",", "§").replace(".", ",").replace("§", ".")


_RC_TESIS = {
    "font.size": 11, "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 9.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "0.3", "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": "0.85", "grid.linewidth": 0.6,
}


def _estilo_tesis():
    """Estilo limpio para la tesis, independiente del estilo de seaborn que usan las figuras anteriores."""
    from contextlib import ExitStack
    stack = ExitStack()
    stack.enter_context(plt.style.context("default"))
    stack.enter_context(plt.rc_context(_RC_TESIS))
    return stack


def _letra_panel(ax, letra: str):
    ax.text(-0.13, 1.04, letra, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def _con_estilo(func):
    import functools
    @functools.wraps(func)
    def envoltura(*a, **k):
        with _estilo_tesis():
            return func(*a, **k)
    return envoltura


@_con_estilo
def figura_lecturas(df: pd.DataFrame, outpath: Path, fig_scale: float = 1.0):
    """(A) Lecturas antes vs. después por muestra, con todas las muestras.
       (B) Distribución del porcentaje de retención por muestra."""
    d = df.dropna(subset=["lecturas_before", "lecturas_after"]).copy()
    if d.empty:
        print("[WARN] Sin datos de lecturas para la figura nueva.", file=sys.stderr); return
    ret = 100.0 * d["lecturas_after"] / d["lecturas_before"]
    med = float(np.median(ret))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0 * fig_scale, 4.6 * fig_scale),
                                   gridspec_kw={"width_ratios": [1.15, 1]})
    fmt_miles = plt.FuncFormatter(lambda v, _: _num_es(v))

    # (A) dispersión antes vs. después
    xmax = float(d["lecturas_before"].max()) * 1.06
    xs = np.array([0, xmax])
    ax1.plot(xs, xs, ls="--", lw=1.1, color="0.55", zorder=1, label="Sin pérdida (100 %)")
    ax1.plot(xs, xs * med / 100.0, lw=1.4, color=C_DESPUES, zorder=2,
             label=f"Retención mediana ({_num_es(med, 1)} %)")
    ax1.scatter(d["lecturas_before"], d["lecturas_after"], s=30, color=C_ANTES, alpha=0.85,
                edgecolor="white", linewidth=0.5, zorder=3, label=f"Muestras (n = {len(d)})")
    ax1.set_xlim(0, xmax); ax1.set_ylim(0, xmax)
    ax1.set_aspect("equal", adjustable="box")
    ax1.xaxis.set_major_formatter(fmt_miles); ax1.yaxis.set_major_formatter(fmt_miles)
    ax1.set_xlabel("Lecturas antes del filtrado")
    ax1.set_ylabel("Lecturas después del filtrado")
    ax1.legend(loc="upper left", frameon=False)
    _letra_panel(ax1, "A")

    # (B) porcentaje de retención por muestra
    lo, hi = np.floor(ret.min()), np.ceil(ret.max())
    bins = np.arange(lo, hi + 0.5, 0.5)
    cuentas, _, _ = ax2.hist(ret, bins=bins, color=C_ANTES, alpha=0.85, edgecolor="white", linewidth=0.8)
    ax2.vlines(med, 0, cuentas.max() * 1.08, color=C_DESPUES, lw=1.8)
    ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: _num_es(v, 0)))
    ax2.set_xlabel("Lecturas retenidas por muestra (%)")
    ax2.set_ylabel("Número de muestras")
    from matplotlib.ticker import MaxNLocator
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.set_ylim(0, ax2.get_ylim()[1] * 1.45)
    ax2.text(0.03, 0.97, f"Mediana: {_num_es(med, 1)} %\n"
                         f"Rango: {_num_es(ret.min(), 1)}–{_num_es(ret.max(), 1)} %\n"
                         f"Total: {_num_es(d['lecturas_before'].sum())} → {_num_es(d['lecturas_after'].sum())} lecturas",
             transform=ax2.transAxes, ha="left", va="top", fontsize=9.5,
             bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.8"))
    _letra_panel(ax2, "B")

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


@_con_estilo
def figura_calidad(df: pd.DataFrame, before_q: Optional[np.ndarray], after_q: Optional[np.ndarray],
                   outpath: Path, q_umbral: Optional[float] = None, fig_scale: float = 1.0):
    """(A) Distribución del QScore medio por lectura antes y después (requiere FASTQ).
       (B) QScore promedio por muestra, antes y después, con una línea por muestra (usa las tablas)."""
    tiene_a = before_q is not None and after_q is not None and np.isfinite(before_q).any() and np.isfinite(after_q).any()
    tiene_b = "QScore promedio_before" in df.columns and "QScore promedio_after" in df.columns
    if not (tiene_a or tiene_b):
        print("[WARN] Sin datos de calidad para la figura nueva.", file=sys.stderr); return
    n = int(tiene_a) + int(tiene_b)
    fig, axes = plt.subplots(1, n, figsize=((5.6 * n) * fig_scale, 4.6 * fig_scale), squeeze=False)
    axes = axes[0]
    i = 0
    if tiene_a:
        ax = axes[i]; i += 1
        a = before_q[np.isfinite(before_q)]; b = after_q[np.isfinite(after_q)]
        bins = np.arange(np.floor(min(a.min(), b.min())), np.ceil(max(a.max(), b.max())) + 0.25, 0.25)
        for v, c, lab in ((a, C_ANTES, "Antes"), (b, C_DESPUES, "Después")):
            w = np.full(v.size, 100.0 / v.size)
            ax.hist(v, bins=bins, weights=w, histtype="stepfilled", color=c, alpha=0.16, linewidth=0)
            ax.hist(v, bins=bins, weights=w, histtype="step", color=c, lw=1.5,
                    label=f"{lab} (mediana {_num_es(np.median(v), 1)})")
            ax.axvline(np.median(v), color=c, ls="--", lw=1.2)
        if q_umbral is not None:
            ax.axvline(q_umbral, color="0.25", lw=1.1, ls=":")
            ax.text(q_umbral, ax.get_ylim()[1] * 0.98, f" umbral Q{_num_es(q_umbral, 0)}", fontsize=9, va="top")
        from matplotlib.ticker import MultipleLocator
        ax.xaxis.set_major_locator(MultipleLocator(2))  # ticks enteros: evita que 17,5 se rotule como "18"
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: _num_es(v, 0)))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: _num_es(v, 1)))
        ax.set_xlabel("QScore medio por lectura")
        ax.set_ylabel("Lecturas (%)")
        ax.legend(loc="upper left", frameon=False)
        _letra_panel(ax, "A")
    if tiene_b:
        ax = axes[i]
        d = df.dropna(subset=["QScore promedio_before", "QScore promedio_after"])
        qb = d["QScore promedio_before"].values; qa = d["QScore promedio_after"].values
        for y0, y1 in zip(qb, qa):
            ax.plot([0, 1], [y0, y1], color="0.6", lw=0.7, alpha=0.5, zorder=1)
        rng = np.random.default_rng(7)
        ax.scatter(rng.uniform(-0.04, 0.04, qb.size), qb, s=16, color=C_ANTES, zorder=3, edgecolor="white", linewidth=0.4)
        ax.scatter(1 + rng.uniform(-0.04, 0.04, qa.size), qa, s=16, color=C_DESPUES, zorder=3, edgecolor="white", linewidth=0.4)
        bp = ax.boxplot([qb, qa], positions=[-0.22, 1.22], widths=0.14, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=1.2))
        for patch, c in zip(bp["boxes"], (C_ANTES, C_DESPUES)):
            patch.set_facecolor(c); patch.set_alpha(0.55)
        delta = float(np.median(qa - qb))
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Antes", "Después"])
        ax.set_xlim(-0.45, 1.45)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: _num_es(v, 1)))
        ax.set_ylabel("QScore promedio por muestra")
        ax.grid(axis="x", visible=False)
        signo = "+" if delta >= 0 else "−"
        ax.text(0.5, 1.01, f"n = {len(d)} muestras · cambio mediano: {signo}{_num_es(abs(delta), 2)}",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=9.5, color="0.2")
        _letra_panel(ax, "B" if tiene_a else "A")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


# ------------------------
# PDF (una sola página, mosaico de imágenes PNG/JPG)
# ------------------------

def build_onepage_pdf(image_paths: List[Path], pdf_path: Path):
    """Crea un PDF de una página con mosaico de todas las imágenes (PNG/JPG)."""
    if not image_paths:
        print("[WARN] No hay imágenes para el PDF.", file=sys.stderr); return

    # Orden por prefijo numérico si existe
    def _key(p: Path):
        m = re.match(r"(\d+)_", p.name)
        return (int(m.group(1)) if m else 9999, p.name)
    image_paths = [p for p in image_paths if p.suffix.lower() in (".png", ".jpg", ".jpeg")]
    image_paths = sorted(image_paths, key=_key)
    if not image_paths:
        print("[WARN] No hay PNG/JPG para el PDF.", file=sys.stderr); return

    n = len(image_paths)
    if n <= 4:
        cols = 2
    elif n <= 9:
        cols = 3
    else:
        cols = 4
    rows = int(np.ceil(n/cols))

    fig_w = max(10, cols*4)
    fig_h = max(7, rows*3.5)
    fig, axes = plt.subplots(rows, cols, figsize=(fig_w, fig_h))
    axes = np.atleast_2d(axes).reshape(rows, cols)

    for i, p in enumerate(image_paths):
        r, c = divmod(i, cols)
        ax = axes[r, c]
        img = imread(str(p))
        ax.imshow(img)
        ax.set_title(p.name, fontsize=8)
        ax.axis('off')

    for j in range(n, rows*cols):
        r, c = divmod(j, cols)
        axes[r, c].axis('off')

    fig.tight_layout()
    fig.savefig(pdf_path, dpi=200)
    plt.close(fig)


# ------------------------
# Fastq → métricas por muestra (si no hay tablas)
# ------------------------

def _safe_name_from_path(p: str) -> str:
    base = Path(p).name
    base = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
    return base


def compute_metrics_from_fastq(paths: list[str], reads_per_file_cap: int | None = None) -> pd.DataFrame:
    rows = []
    for p in paths:
        name = normalize_sample_name(_safe_name_from_path(p))
        lengths = []
        q_all = []
        gc_bases = 0
        total_bases = 0
        taken = 0
        try:
            for seq, qual in _iter_fastq_records(p):
                L = len(seq)
                lengths.append(L)
                total_bases += L
                gc_bases += sum(1 for c in seq if c in "GgCc")
                for ch in qual:
                    q_all.append(ord(ch)-33)
                taken += 1
                if reads_per_file_cap and taken >= reads_per_file_cap:
                    break
        except Exception as e:
            print(f"[WARN] No se pudo leer {p}: {e}", file=sys.stderr)
            continue

        if not lengths:
            continue
        arr = np.array(lengths, dtype=float)
        lecturas = float(len(arr))
        bases = float(arr.sum())
        lprom = float(arr.mean())
        lmin = float(arr.min())
        lmax = float(arr.max())
        s = np.sort(arr)
        csum = np.cumsum(s)
        half = bases/2.0
        idx = np.searchsorted(csum, half)
        n50 = float(s[idx]) if idx < s.size else float('nan')
        gc_pct = float(gc_bases/total_bases*100.0) if total_bases>0 else float('nan')
        qscore_prom = float(np.mean(q_all)) if q_all else float('nan')

        rows.append({
            'Muestra': name,
            'lecturas': lecturas,
            'bases': bases,
            'longitud promedio': lprom,
            'longitud mínima': lmin,
            'longitud máxima': lmax,
            'N50': n50,
            'GC%': gc_pct,
            'QScore promedio': qscore_prom
        })
    return pd.DataFrame(rows)


# ------------------------
# Main
# ------------------------

def main():
    ap = argparse.ArgumentParser(description="Comparación gráfica de métricas antes vs. después (FASTQ).")

    # Tablas / FASTQ
    ap.add_argument("--before", required=False, default=None, help="Ruta a métricas ANTES (xlsx/tsv/txt bonito). Opcional si usas --derive-metrics.")
    ap.add_argument("--after", required=False, default=None, help="Ruta a métricas DESPUÉS (xlsx/tsv/txt bonito). Opcional si usas --derive-metrics.")
    ap.add_argument("--derive-metrics", action="store_true", help="Deriva métricas por muestra leyendo los FASTQ (antes/después)")
    ap.add_argument("--metrics-reads-per-file-cap", type=int, default=None, help="Cap de lecturas por archivo al derivar métricas")

    # FASTQ per-read (longitudes / QScore)
    ap.add_argument("--fastq-before", default=None, help="Glob o carpeta con FASTQ ANTES (e.g., '/raw/**/*.fastq.gz')")
    ap.add_argument("--fastq-after",  default=None, help="Glob o carpeta con FASTQ DESPUÉS (e.g., '/clean/**/*.fastq.gz')")
    ap.add_argument("--reads-cap", type=int, default=120000, help="Máx. lecturas totales a leer por grupo (antes/después)")
    ap.add_argument("--reads-per-file-cap", type=int, default=20000, help="Máx. lecturas por archivo")
    ap.add_argument("--after-filter", default="*_final.fastq*", help="Patrón de nombre para FASTQ del grupo DESPUÉS (ej. '*_final.fastq*'). Usa vacío para no filtrar.")

    # Salidas / rutas
    ap.add_argument("--outdir", required=False, default=None, help="Carpeta de salida para archivos (PNGs, CSV, etc.)")
    ap.add_argument("--results-dir", default=None, help="Alias de --outdir. Si se define, tiene prioridad.")
    ap.add_argument("--prefix", default="", help="Prefijo para archivos de salida (p. ej., 'run1_').")
    ap.add_argument("--pdf-name", default="figuras_antes_despues.pdf", help="Nombre del PDF que se generará.")
    ap.add_argument("--pdf-onepage", action="store_true", help="Si se activa, crea un PDF con todas las figuras en una sola página.")

    # Mancuernas
    ap.add_argument("--max-samples", type=int, default=25, help="Máximo de muestras a mostrar en dumbbells (0=sin límite).")
    ap.add_argument("--subset-mode", choices=["highest_removal","largest_before","largest_after","random"], default="highest_removal",
                    help="Criterio para elegir el subconjunto en dumbbells")
    ap.add_argument("--annotate-top", type=int, default=0, help="Anotar solo Top-N por %% remoción en dumbbells (0=todas).")
    ap.add_argument("--use-plotly", action="store_true", help="Usa Plotly para las mancuernas (HTML interactivo; PNG si kaleido está instalado).")

    # Escalas / estética
    ap.add_argument("--fig-scale", type=float, default=1.0, help="Factor global para escalar tamaños de figura.")
    ap.add_argument("--lineheight", type=float, default=0.38, help="Alto por muestra en mancuernas (pulgadas por fila).")
    ap.add_argument("--lengths-log10", action="store_true", help="Usa escala log10 para el violín de longitudes por lectura.")
    ap.add_argument("--raincloud", action="store_true", help="Usa gráficos tipo raincloud (ptitprince). Si no está, cae a violín estándar.")

    # Figuras nuevas
    ap.add_argument("--estilo", choices=["nuevo", "anterior", "ambos"], default="nuevo",
                    help="'nuevo': 2 figuras (lecturas y calidad, todas las muestras). "
                         "'anterior': mancuernas, violines y scatter originales. 'ambos': todo.")
    ap.add_argument("--q-por-lectura", choices=["aritmetica", "probabilidad"], default="aritmetica",
                    help="Cómo calcular el QScore medio de cada lectura. 'aritmetica' coincide con stats_fastq.py.")
    ap.add_argument("--q-umbral", type=float, default=None,
                    help="Dibuja una línea vertical en este QScore en la figura de calidad (opcional).")

    args = ap.parse_args()

    # Resolver carpeta de salida
    if args.results_dir:
        args.outdir = args.results_dir
    if not args.outdir:
        raise SystemExit("Debes pasar --outdir o --results-dir")
    outdir = Path(args.outdir).expanduser().resolve()
    ensure_outdir(outdir)

    # Construcción de tablas de métricas
    if args.derive_metrics:
        if not args.fastq_before or not args.fastq_after:
            raise SystemExit("--derive-metrics requiere --fastq-before y --fastq-after")

        paths_b = _expand_fastq(args.fastq_before)
        paths_a = _expand_fastq(args.fastq_after)
        if args.after_filter:
            paths_a = [p for p in paths_a if fnmatch(Path(p).name, args.after_filter)]

        if not paths_b or not paths_a:
            raise SystemExit("No se encontraron FASTQ para derivar métricas.")

        print(f"[INFO] FASTQ antes: {len(paths_b)} | FASTQ después (filtrados): {len(paths_a)}", file=sys.stderr)

        print("[INFO] Derivando métricas ANTES desde FASTQ…", file=sys.stderr)
        before = compute_metrics_from_fastq(paths_b, reads_per_file_cap=args.metrics_reads_per_file_cap)

        print("[INFO] Derivando métricas DESPUÉS desde FASTQ…", file=sys.stderr)
        after  = compute_metrics_from_fastq(paths_a, reads_per_file_cap=args.metrics_reads_per_file_cap)

    else:
        if not args.before or not args.after:
            raise SystemExit("Debes pasar --before y --after, o usar --derive-metrics.")
        before = load_table(Path(args.before))
        after  = load_table(Path(args.after))

    # Verifica columnas mínimas
    for col in REQUIRED_FOR_PLOTS:
        if col not in before.columns:
            raise ValueError(f"Falta columna '{col}' en BEFORE")
        if col not in after.columns:
            raise ValueError(f"Falta columna '{col}' en AFTER")

    df = join_before_after(before, after)
    if df.empty:
        sys.exit(1)

    # Helper para nombres con prefijo
    def OUT(name: str) -> Path:
        return outdir / f"{args.prefix}{name}"

    produced: List[Path] = []
    viejo = args.estilo in ("anterior", "ambos")
    nuevo = args.estilo in ("nuevo", "ambos")

    # QScore por lectura (se lee una sola vez y se usa en las figuras que lo necesiten)
    b_len = b_q = a_len = a_q = None
    if args.fastq_before and args.fastq_after:
        paths_b = _expand_fastq(args.fastq_before)
        paths_a = _expand_fastq(args.fastq_after)
        if args.after_filter:
            paths_a = [p for p in paths_a if fnmatch(Path(p).name, args.after_filter)]
        if not paths_b:
            print(f"[WARN] No FASTQ en --fastq-before: {args.fastq_before}", file=sys.stderr)
        if not paths_a:
            print(f"[WARN] No FASTQ en --fastq-after (tras filtro): {args.fastq_after}", file=sys.stderr)
        if paths_b and paths_a:
            b_len, b_q = sample_lengths_qscores(paths_b, args.reads_cap, args.reads_per_file_cap, q_modo=args.q_por_lectura)
            a_len, a_q = sample_lengths_qscores(paths_a, args.reads_cap, args.reads_per_file_cap, q_modo=args.q_por_lectura)
    else:
        print("[INFO] Sin FASTQ: la figura de calidad solo tendrá el panel por muestra.", file=sys.stderr)

    # ====== 0) Figuras nuevas ======
    if nuevo:
        fA = OUT("01_lecturas_retencion.png")
        fB = OUT("02_calidad_qscore.png")
        figura_lecturas(df, fA, fig_scale=args.fig_scale)
        figura_calidad(df, b_q, a_q, fB, q_umbral=args.q_umbral, fig_scale=args.fig_scale)
        for p in (fA, fB):
            if p.exists():
                produced.append(p)

    # ====== 1) Mancuernas ======
    if not viejo:
        pass
    elif args.use_plotly:
        html1 = OUT("01_dumbbell_lecturas.html")
        html2 = OUT("02_dumbbell_bases.html")
        png1  = OUT("01_dumbbell_lecturas.png")
        png2  = OUT("02_dumbbell_bases.png")
        dumbbell_plot_plotly(df, "lecturas", html1, png1, max_samples=args.max_samples,
                             subset_mode=args.subset_mode, annotate_top=args.annotate_top)
        dumbbell_plot_plotly(df, "bases", html2, png2, max_samples=args.max_samples,
                             subset_mode=args.subset_mode, annotate_top=args.annotate_top)
        for p in [png1, png2, html1, html2]:
            if p.exists():
                produced.append(p)
    else:
        f1 = OUT("01_dumbbell_lecturas.png")
        f2 = OUT("02_dumbbell_bases.png")
        dumbbell_plot(df, "lecturas", f1, max_samples=args.max_samples, subset_mode=args.subset_mode,
                      annotate_top=args.annotate_top, fig_scale=args.fig_scale, lineheight=args.lineheight)
        dumbbell_plot(df, "bases",    f2, max_samples=args.max_samples, subset_mode=args.subset_mode,
                      annotate_top=args.annotate_top, fig_scale=args.fig_scale, lineheight=args.lineheight)
        for p in [f1, f2]:
            if p.exists():
                produced.append(p)

    # ====== 2) Longitudes por lectura + Scatter QScore (si hay FASTQ) ======
    if viejo and b_q is not None and a_q is not None:
        f3 = OUT("03_violin_read_lengths.png")
        f4 = OUT("04_scatter_qscore_per_read.png")
        violin_lengths(b_len, a_len, f3, lengths_log10=args.lengths_log10,
                       fig_scale=args.fig_scale, use_raincloud=args.raincloud)
        scatter_qscores(b_q, a_q, f4, fig_scale=args.fig_scale)
        for p in [f3, f4]:
            if p.exists():
                produced.append(p)

    # ====== 3) N50 por muestra ======
    if viejo:
        f5 = OUT("05_violin_n50.png")
        violin_n50(df, f5, fig_scale=args.fig_scale, use_raincloud=args.raincloud)
        if f5.exists():
            produced.append(f5)

    # ====== 4) CSV resumen (retención/remoción) ======
    try:
        sel = df[[
            "Muestra",
            "lecturas_before","lecturas_after","lecturas_retencion_%","lecturas_remocion_%",
            "bases_before","bases_after","bases_retencion_%","bases_remocion_%"
        ]].copy()
        csv_out = OUT("00_resumen_retencion_remocion.csv")
        sel.to_csv(csv_out, index=False)
        produced.append(csv_out)
    except Exception as e:
        print(f"[WARN] No se exportó resumen retención/remocion: {e}", file=sys.stderr)

    # ====== 5) PDF una sola página (mosaico) ======
    if args.pdf_onepage:
        pdf_path = OUT(args.pdf_name)
        imgs = [p for p in produced if p.suffix.lower() in (".png", ".jpg", ".jpeg")]
        build_onepage_pdf(imgs, pdf_path)
        if pdf_path.exists():
            print(f"[OK] PDF (una página) creado: {pdf_path.resolve()}")

    print(f"[OK] Resultados guardados en: {outdir.resolve()}")


if __name__ == "__main__":
    main()
