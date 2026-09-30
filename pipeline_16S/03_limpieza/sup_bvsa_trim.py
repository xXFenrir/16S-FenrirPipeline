#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Comparación gráfica ANTES vs. DESPUÉS usando Seaborn (Formato Vertical / Lista de Barcodes).
Métricas soportadas: Lecturas, Bases y QScore.
"""

import argparse
import re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Configuración estética global de Seaborn
sns.set_theme(style="whitegrid", palette="muted")
sns.set_context("paper", font_scale=1.1)

def load_data(path: Path) -> pd.DataFrame:
    """Carga archivos TSV, CSV o XLSX."""
    ext = path.suffix.lower()
    if ext == ".xlsx":
        df = pd.read_excel(path)
    elif ext in (".tsv", ".tab"):
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    return df

def simplify_barcode_name(name: str) -> str:
    """Convierte cadenas como 'barcode01' o 'pass_barcode01_...' a 'B01'."""
    match = re.search(r'barcode(\d+)', str(name), re.IGNORECASE)
    if match:
        num = int(match.group(1))
        return f"B{num:02d}"
    return str(name)

def prepare_tidy_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Estandariza columnas y transforma la tabla a formato 'tidy' (long format).
    Abrevia las etiquetas de los barcodes a B01, B02, etc.
    """
    sample_col = "barcode" if "barcode" in df.columns else ("Muestra" if "Muestra" in df.columns else df.columns[0])
    
    records = []
    for idx, row in df.iterrows():
        sample = simplify_barcode_name(row[sample_col])
        
        # 1. Lecturas
        if 'lecturas raw' in df.columns:
            records.append({'Muestra': sample, 'Métrica': 'Lecturas', 'Estado': 'Antes', 'Valor': row['lecturas raw']})
        if 'lecturas after filtlong' in df.columns:
            records.append({'Muestra': sample, 'Métrica': 'Lecturas', 'Estado': 'Después', 'Valor': row['lecturas after filtlong']})
        elif 'lecturas after pychopper' in df.columns:
            records.append({'Muestra': sample, 'Métrica': 'Lecturas', 'Estado': 'Después', 'Valor': row['lecturas after pychopper']})
            
        # 2. Bases
        if 'bases raw' in df.columns:
            records.append({'Muestra': sample, 'Métrica': 'Bases', 'Estado': 'Antes', 'Valor': row['bases raw']})
        if 'bases after' in df.columns:
            records.append({'Muestra': sample, 'Métrica': 'Bases', 'Estado': 'Después', 'Valor': row['bases after']})
            
        # 3. QScore
        if 'QScore promedio raw' in df.columns:
            records.append({'Muestra': sample, 'Métrica': 'QScore', 'Estado': 'Antes', 'Valor': row['QScore promedio raw']})
        if 'QScore promedio after' in df.columns:
            records.append({'Muestra': sample, 'Métrica': 'QScore', 'Estado': 'Después', 'Valor': row['QScore promedio after']})
            
    return pd.DataFrame(records)

def plot_vertical_paired_bars(df_metric: pd.DataFrame, metric_name: str, outpath: Path):
    """Grafica la comparación en formato vertical (Y=Muestra, X=Valor)."""
    
    n_samples = df_metric['Muestra'].nunique()
    
    # Altura dinámica: ~0.28 pulgadas por barcode para que no se encimen las barras
    fig_height = max(8, n_samples * 0.28)
    plt.figure(figsize=(10, fig_height))
    
    ax = sns.barplot(
        data=df_metric,
        y='Muestra', x='Valor', hue='Estado',
        palette={'Antes': '#1f77b4', 'Después': '#ff7f0e'},
        edgecolor='black', linewidth=0.5, alpha=0.9
    )
    
    ax.set_title(f'Comparación por Barcode: {metric_name} (Antes vs. Después)', fontsize=14, pad=12)
    ax.set_ylabel('Barcode', fontsize=12)
    ax.set_xlabel(metric_name, fontsize=12)
    
    # Formateo numérico para miles/millones en el eje horizontal
    if metric_name in ['Lecturas', 'Bases']:
        ax.xaxis.set_major_formatter('{x:,.0f}')
        
    ax.legend(title="Estado", loc='lower right')
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

def plot_distribution_summary(tidy_df: pd.DataFrame, outpath: Path):
    """Grafica la distribución global (Boxplot + Stripplot) para las 3 métricas."""
    g = sns.catplot(
        data=tidy_df,
        x='Estado', y='Valor', hue='Estado',
        col='Métrica', sharey=False, 
        kind='box',
        height=4.5, aspect=0.9, boxprops=dict(alpha=0.6)
    )
    
    g.map_dataframe(
        sns.stripplot,
        x='Estado', y='Valor', hue='Estado',
        dodge=False, jitter=0.2, alpha=0.7, size=4, palette={'Antes': '#1f77b4', 'Después': '#ff7f0e'}
    )
    
    g.fig.subplots_adjust(top=0.82)
    g.fig.suptitle('Resumen Global de Filtrado: Distribución de Métricas', fontsize=14)
    g.set_axis_labels("", "Valor")
    
    plt.savefig(outpath, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Graficar comparaciones Antes vs. Después en formato vertical.")
    parser.add_argument("--input", required=True, help="Ruta al archivo de resumen (TSV, CSV, XLSX).")
    parser.add_argument("--outdir", required=True, help="Carpeta de salida para las figuras.")
    
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    raw_df = load_data(Path(args.input))
    tidy_df = prepare_tidy_data(raw_df)
    
    if tidy_df.empty:
        print("[ERROR] No se pudieron mapear las columnas requeridas (Lecturas, Bases, QScore).")
        return

    # Generar barras pareadas verticales por cada métrica
    for metric in ['Lecturas', 'Bases', 'QScore']:
        df_sub = tidy_df[tidy_df['Métrica'] == metric]
        if not df_sub.empty:
            out_png = outdir / f"comparacion_{metric.lower()}.png"
            plot_vertical_paired_bars(df_sub, metric, out_png)
            print(f"[OK] Gráfica guardada (Formato Vertical): {out_png}")
            
    # Resumen de distribución global
    out_dist = outdir / "resumen_distribucion_global.png"
    plot_distribution_summary(tidy_df, out_dist)
    print(f"[OK] Gráfica global guardada: {out_dist}")

if __name__ == "__main__":
    main()
