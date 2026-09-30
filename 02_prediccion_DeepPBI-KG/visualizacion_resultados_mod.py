import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import plotly.graph_objects as go
import os
import argparse

def main():
    # =====================================================================
    # 1. RUTAS Y CONFIGURACIÓN (ahora por parámetro, no fijas)
    # =====================================================================
    parser = argparse.ArgumentParser(description="Genera gráficas de densidad y embudo a partir de result.csv de DeepPBI-KG.")
    parser.add_argument('--dir_base', type=str, required=True,
                         help="Carpeta que contiene result.csv (normalmente <combo>/output).")
    parser.add_argument('--umbral', type=float, default=0.85,
                         help="Umbral de probabilidad para las líneas/cuadrantes (por defecto: 0.85).")
    parser.add_argument('--etiqueta', type=str, default="",
                         help="Etiqueta opcional para incluir en los títulos (ej. 'HAC - todo').")
    args = parser.parse_args()

    DIR_BASE = args.dir_base
    ARCHIVO_RESULTADOS = os.path.join(DIR_BASE, "result.csv")
    DIR_OUTPUT = os.path.join(DIR_BASE, "graficas_tesis")
    os.makedirs(DIR_OUTPUT, exist_ok=True)
    UMBRAL = args.umbral
    sufijo_titulo = f" ({args.etiqueta})" if args.etiqueta else ""

    print(f"Cargando matriz teórica desde: {ARCHIVO_RESULTADOS}...")
    try:
        df = pd.read_csv(ARCHIVO_RESULTADOS)
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo {ARCHIVO_RESULTADOS}.")
        return

    col_kg, col_wgs = "key_gene_output", "wgs_output"

    # LIMPIEZA: Rellenar vacíos y limitar a 1.0 para evitar errores matemáticos
    df[col_kg] = df[col_kg].fillna(0.0).clip(lower=0.0, upper=1.0)
    df[col_wgs] = df[col_wgs].fillna(0.0).clip(lower=0.0, upper=1.0)

    # Probabilidad compuesta
    col_comp = "weight_output"
    df[col_comp] = (df[col_kg] + df[col_wgs]) / 2.0
    total_interacciones = len(df)

    print(f"¡Matriz cargada! Procesando {total_interacciones} interacciones...\n")

    # =====================================================================
    # 2. GRÁFICO DE DENSIDAD 2D
    # =====================================================================
    print("1/2 Generando Gráfico de Densidad 2D...")
    sns.set_theme(style="whitegrid")
    fig1, ax1 = plt.subplots(figsize=(10, 8))

    h = ax1.hist2d(df[col_wgs], df[col_kg], bins=100, cmap="mako", cmin=1, norm=mcolors.LogNorm())
    cbar = fig1.colorbar(h[3], ax=ax1)

    cbar.set_label('Densidad de Interacciones', rotation=270, labelpad=20, fontweight='bold')

    # Líneas de umbral
    ax1.axvline(UMBRAL, color='#e63946', linestyle='--', linewidth=2)
    ax1.axhline(UMBRAL, color='#e63946', linestyle='--', linewidth=2)

    # Sombreado sutil del cuadrante ideal
    ax1.axhspan(UMBRAL, 1.0, xmin=UMBRAL, xmax=1.0, color='#e63946', alpha=0.15)

    ax1.set_xlabel("Probabilidad WGS", fontsize=12, fontweight='bold')
    ax1.set_ylabel("Probabilidad Genes Clave", fontsize=12, fontweight='bold')
    ax1.set_title(f"Densidad de interacciones{sufijo_titulo}", fontsize=13, fontweight='bold')
    ax1.set_xlim(0, 1.0)
    ax1.set_ylim(0, 1.0)

    plt.tight_layout()
    fig1.savefig(os.path.join(DIR_OUTPUT, "1_Densidad_Cuadrantes.png"), dpi=300)
    plt.close(fig1)

    # =====================================================================
    # 3. GRÁFICO DE EMBUDO (FUNNEL)
    # =====================================================================
    print("2/2 Generando Diagrama de Embudo...")
    etapas = [
        "Interacciones totales",
        f"Interacción genes clave ≥ {int(UMBRAL*100)}%",
        f"Interacción WGS ≥ {int(UMBRAL*100)}%",
        f"Interacción compuesta (Promedio) ≥ {int(UMBRAL*100)}%"
    ]

    set_kg = set(df[df[col_kg] >= UMBRAL].index)
    set_wgs = set(df[df[col_wgs] >= UMBRAL].index)
    set_comp = set(df[df[col_comp] >= UMBRAL].index)
    valores = [total_interacciones, len(set_kg), len(set_wgs), len(set_comp)]

    total_inicial = valores[0]
    textos_personalizados = [
        f"{f'{v:,}'.replace(',', '.')}<br>{(v/total_inicial)*100:.1f}%"
        for v in valores
    ]

    colores_vibrantes = ["#00BFFF", "#32CD32", "#FF3333", "#FFA500"]

    fig3 = go.Figure(go.Funnel(
        y=etapas,
        x=valores,
        text=textos_personalizados,
        textinfo="text",
        textposition="auto",
        textfont=dict(size=14, color="black", family="Arial"),
        marker=dict(color=colores_vibrantes)
    ))

    fig3.update_layout(
        title=dict(text=f"Reducción de Dimensionalidad Computacional{sufijo_titulo}", font=dict(size=16), x=0.5),
        plot_bgcolor="white", paper_bgcolor="white", showlegend=False, margin=dict(l=20, r=20, t=60, b=20)
    )

    fig3.write_image(os.path.join(DIR_OUTPUT, "3_Funnel_Reduccion.png"), scale=3)

    print("\n" + "="*60)
    print(f"¡Proceso finalizado! Las gráficas se guardaron en:\n{DIR_OUTPUT}")
    print("="*60)

if __name__ == "__main__":
    main()
