#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Comparación gráfica POST-LIMPIEZA entre modelos de Basecalling (HAC vs. SUP).
Formatos soportados: Vertical por Barcode y Resumen Distributivo.
"""

import argparse
import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Configuración estética global
sns.set_theme(style="whitegrid", palette="muted")
sns.set_context("paper", font_scale=1.1)

def load_data(path: Path) -> pd.DataFrame:
    """Carga archivos TSV, CSV o XLSX verificando previamente su existencia."""
    if not path.exists():
        raise FileNotFoundError(f"[ERROR] No se encontró el archivo en la ruta especificadas:\n  -> {path}")
        
    ext = path.suffix.lower()
    if ext == ".xlsx":
        return pd.read_excel(path)
    elif ext in (".tsv", ".tab"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)

def simplify_barcode_name(name: str) -> str:
    """Abrevia nombres de barcodes a formato B01, B02, etc."""
    match = re.search(r'barcode(\d+)', str(name), re.IGNORECASE)
    if match:
        num = int(match.group(1))
        return f"B{num:02d}"
    return str(name)

def extract_clean_metrics(df: pd.DataFrame, method_label: str) -> pd.DataFrame:
    """Extrae únicamente los datos recuperados POST-LIMPIEZA."""
    sample_col = "barcode" if "barcode" in df.columns else ("Muestra" if "Muestra" in df.columns else df.columns[0])
    
    records = []
    for idx, row in df.iterrows():
        sample = simplify_barcode_name(row[sample_col])
        
        # Lecturas limpias (filtlong o pychopper)
        lecturas_col = 'lecturas after filtlong' if 'lecturas after filtlong' in df.columns else 'lecturas after pychopper'
        if lecturas_col in df.columns:
            records.append({'Muestra': sample, 'Método': method_label, 'Métrica': 'Lecturas Limpias', 'Valor': row[lecturas_col]})
            
        # QScore promedio post-filtro
        if 'QScore promedio after' in df.columns:
            records.append({'Muestra': sample, 'Método': method_label, 'Métrica': 'QScore Limpio', 'Valor': row['QScore promedio after']})
            
    return pd.DataFrame(records)

def plot_vertical_hac_vs_sup(df_metric: pd.DataFrame, metric_name: str, outpath: Path):
    """Genera gráfico vertical con pares HAC vs. SUP por cada barcode."""
    n_samples = df_metric['Muestra'].nunique()
    
    # Altura dinámica: ~0.28 pulgadas por muestra
    fig_height = max(8, n_samples * 0.28)
    plt.figure(figsize=(10, fig_height))
    
    ax = sns.barplot(
        data=df_metric,
        y='Muestra', x='Valor', hue='Método',
        palette={'HAC': '#2b5c8f', 'SUP': '#d95f02'}, # Azul para HAC, Naranja/Rojo para SUP
        edgecolor='black', linewidth=0.5, alpha=0.9
    )
    
    ax.set_title(f'Comparación Post-Limpieza: {metric_name} (HAC vs. SUP)', fontsize=14, pad=12)
    ax.set_ylabel('Barcode', fontsize=12)
    ax.set_xlabel(metric_name, fontsize=12)
    
    # Formato numérico horizontal
    if 'Lecturas' in metric_name:
        ax.xaxis.set_major_formatter('{x:,.0f}')
        
    ax.legend(title="Modelo Basecalling", loc='lower right')
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

def plot_global_comparison(tidy_df: pd.DataFrame, outpath: Path):
    """Resumen de caja y puntos (Boxplot) comparando el desempeño global HAC vs. SUP."""
    g = sns.catplot(
        data=tidy_df,
        x='Método', y='Valor', hue='Método',
        col='Métrica', sharey=False, 
        kind='box',
        palette={'HAC': '#2b5c8f', 'SUP': '#d95f02'},
        height=4.5, aspect=0.9, boxprops=dict(alpha=0.6)
    )
    
    g.map_dataframe(
        sns.stripplot,
        x='Método', y='Valor', hue='Método',
        dodge=False, jitter=0.2, alpha=0.7, size=4,
        palette={'HAC': '#2b5c8f', 'SUP': '#d95f02'}
    )
    
    g.fig.subplots_adjust(top=0.82)
    g.fig.suptitle('Comparativa Global Post-Limpieza: HAC vs. SUP', fontsize=14)
    g.set_axis_labels("", "")
    # Etiqueta del eje Y con el nombre de cada métrica (en vez de "Valor")
    ylabels = {'Lecturas Limpias': 'Lecturas', 'QScore Limpio': 'QScore'}
    for metric, ax in g.axes_dict.items():
        ax.set_ylabel(ylabels.get(metric, metric))
    
    plt.savefig(outpath, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Contrastar métricas limpias entre HAC y SUP.")
    parser.add_argument("--hac", required=True, help="Ruta al archivo de resumen HAC.")
    parser.add_argument("--sup", required=True, help="Ruta al archivo de resumen SUP.")
    parser.add_argument("--outdir", required=True, help="Carpeta de salida.")
    
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Cargar datos de ambos métodos
    df_hac_raw = load_data(Path(args.hac))
    df_sup_raw = load_data(Path(args.sup))
    
    # Extraer métricas limpias
    df_hac_clean = extract_clean_metrics(df_hac_raw, 'HAC')
    df_sup_clean = extract_clean_metrics(df_sup_raw, 'SUP')
    
    # Unir datasets
    tidy_df = pd.concat([df_hac_clean, df_sup_clean], ignore_index=True)
    
    if tidy_df.empty:
        print("[ERROR] No se encontraron columnas compatibles en los archivos.")
        return

    # Generar gráficos verticales por cada métrica
    for metric in tidy_df['Métrica'].unique():
        df_sub = tidy_df[tidy_df['Métrica'] == metric]
        out_png = outdir / f"comparativo_hac_vs_sup_{metric.lower().replace(' ', '_')}.png"
        plot_vertical_hac_vs_sup(df_sub, metric, out_png)
        print(f"[OK] Gráfica vertical generada: {out_png}")
            
    # Generar gráfico de caja de distribución global
    out_dist = outdir / "resumen_global_hac_vs_sup.png"
    plot_global_comparison(tidy_df, out_dist)
    print(f"[OK] Gráfica de distribución global generada: {out_dist}")

if __name__ == "__main__":
    main()
