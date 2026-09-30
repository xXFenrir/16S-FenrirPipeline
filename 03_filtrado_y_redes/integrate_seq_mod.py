import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3
import plotly.graph_objects as go
import os

def cargar_y_estandarizar(ruta, nombre_columna_prob):
    """Carga un CSV y estandariza los nombres de las columnas para poder unificarlos."""
    df = pd.read_csv(ruta)
    columnas_originales = list(df.columns)
    # Asumimos que col 0 = Fago, col 1 = Bacteria, col 2 = Probabilidad
    df = df.rename(columns={
        columnas_originales[0]: 'Bacteriofago', 
        columnas_originales[1]: 'Bacteria', 
        columnas_originales[2]: nombre_columna_prob
    })
    return df[['Bacteriofago', 'Bacteria', nombre_columna_prob]]

def main():
    # =====================================================================
    # 1. RUTAS Y DIRECTORIOS
    # =====================================================================
    DIR_BASE = "/home/fenrir/Documentos/Tesis/Algoritmo/comprobacion/resultados_prediccion"
    
    archivo_wgs = os.path.join(DIR_BASE, "interacciones_wgs.csv")
    archivo_kg = os.path.join(DIR_BASE, "interacciones_key_gene.csv")
    archivo_comp = os.path.join(DIR_BASE, "interacciones_comp.csv")
    
    DIR_OUTPUT = os.path.join(DIR_BASE, "graficas_tesis")
    os.makedirs(DIR_OUTPUT, exist_ok=True)
    
    UMBRAL = 0.85

    # =====================================================================
    # 2. CARGA Y UNIFICACIÓN DE LAS TRES MATRICES
    # =====================================================================
    print("Cargando y unificando matrices de predicción...")
    try:
        df_wgs = cargar_y_estandarizar(archivo_wgs, 'Prob_WGS')
        df_kg = cargar_y_estandarizar(archivo_kg, 'Prob_KeyGene')
        df_comp = cargar_y_estandarizar(archivo_comp, 'Prob_Compuesta')
    except FileNotFoundError as e:
        print(f"Error cargando archivos: {e}")
        return

    # Unir las tablas usando el Fago y la Bacteria como llaves (Outer join para no perder nada)
    df_merged = df_kg.merge(df_wgs, on=['Bacteriofago', 'Bacteria'], how='outer')
    df_merged = df_merged.merge(df_comp, on=['Bacteriofago', 'Bacteria'], how='outer')
    
    # Rellenar valores vacíos con 0 (cruces que aparecieron en una tabla pero no en otra)
    df_merged = df_merged.fillna(0)
    total_interacciones = len(df_merged)
    print(f"¡Matriz unificada con éxito! Procesando {total_interacciones} interacciones teóricas...\n")

    # =====================================================================
    # 3. GRÁFICO DE DENSIDAD 2D CON CUADRANTES
    # =====================================================================
    print("1/4 Generando Gráfico de Densidad 2D...")
    sns.set_theme(style="whitegrid")
    fig1, ax1 = plt.subplots(figsize=(10, 8))
    
    sns.histplot(
        data=df_merged, x='Prob_WGS', y='Prob_KeyGene', 
        bins=100, cmap="mako", cbar=True, cbar_kws={'label': 'Densidad de Interacciones'}, ax=ax1
    )
    
    # Líneas de cuadrantes
    ax1.axvline(UMBRAL, color='#e63946', linestyle='--', linewidth=2)
    ax1.axhline(UMBRAL, color='#e63946', linestyle='--', linewidth=2)
    
    # Sombreado del cuadrante ideal
    ax1.axhspan(UMBRAL, 1.0, xmin=UMBRAL, xmax=1.0, color='#e63946', alpha=0.15)
    ax1.text(UMBRAL + 0.01, 0.95, 'Interactoma Óptimo\n(Alta Lisis + Alto WGS)', 
             color='#d62828', fontsize=11, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    ax1.set_xlabel("Probabilidad WGS (Viabilidad Genómica Global)", fontsize=12, fontweight='bold')
    ax1.set_ylabel("Probabilidad Key Gene (Genes de Infección Lítica)", fontsize=12, fontweight='bold')
    ax1.set_xlim(0, 1.0)
    ax1.set_ylim(0, 1.0)
    
    plt.tight_layout()
    fig1.savefig(f"{DIR_OUTPUT}/1_Densidad_Cuadrantes.png", dpi=300)
    plt.close(fig1)

    # =====================================================================
    # 4. DIAGRAMA DE VENN DE INTERACCIONES
    # =====================================================================
    print("2/4 Generando Diagrama de Venn...")
    # Crear sets con los índices de las interacciones que superan el umbral
    set_kg = set(df_merged[df_merged['Prob_KeyGene'] >= UMBRAL].index)
    set_wgs = set(df_merged[df_merged['Prob_WGS'] >= UMBRAL].index)
    set_comp = set(df_merged[df_merged['Prob_Compuesta'] >= UMBRAL].index)

    fig2, ax2 = plt.subplots(figsize=(8, 8))
    venn_out = venn3(
        [set_kg, set_wgs, set_comp], 
        set_labels=(f'Key Gene $\ge$ {UMBRAL}', f'WGS $\ge$ {UMBRAL}', f'Compuesto $\ge$ {UMBRAL}'),
        set_colors=('#457b9d', '#e63946', '#2a9d8f'), alpha=0.7, ax=ax2
    )
    plt.title("Intersección de Viabilidad Genómica y Ecológica", fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    fig2.savefig(f"{DIR_OUTPUT}/2_Diagrama_Venn.png", dpi=300)
    plt.close(fig2)

    # =====================================================================
    # 5. GRÁFICO DE EMBUDO (FUNNEL)
    # =====================================================================
    print("3/4 Generando Diagrama de Embudo...")
    etapas = [
        "Total de cruces teóricos", 
        f"Filtrado Key Gene ($\ge${UMBRAL})", 
        f"Filtrado WGS ($\ge${UMBRAL})", 
        f"Interactoma Compuesto ($\ge${UMBRAL})"
    ]
    valores = [total_interacciones, len(set_kg), len(set_wgs), len(set_comp)]

    fig3 = go.Figure(go.Funnel(
        y = etapas,
        x = valores,
        textinfo = "value+percent initial",
        textposition = "inside",
        textfont = dict(size=14, color="white", family="Arial"),
        marker = dict(color = ["#1d3557", "#457b9d", "#e63946", "#2a9d8f"])
    ))

    fig3.update_layout(
        title = dict(text="Reducción de Dimensionalidad Computacional", font=dict(size=16), x=0.5),
        plot_bgcolor="white", paper_bgcolor="white", showlegend=False, margin=dict(l=20, r=20, t=60, b=20)
    )
    
    fig3.write_image(f"{DIR_OUTPUT}/3_Funnel_Reduccion.png", scale=3)

    # =====================================================================
    # 6. TABLA DEL "TOP 5" INTERACCIONES
    # =====================================================================
    print("4/4 Exportando Top 5 Interacciones...")
    # Ordenar por el score compuesto
    top5_df = df_merged.sort_values(by='Prob_Compuesta', ascending=False).head(5)
    
    for col in ['Prob_KeyGene', 'Prob_WGS', 'Prob_Compuesta']:
        top5_df[col] = top5_df[col].apply(lambda x: round(x, 4))
        
    # Renombrar columnas para exportación limpia
    top5_df.columns = ["Bacteriófago", "Bacteria Hospedera", "Probabilidad Genes Clave", "Probabilidad WGS", "Probabilidad Compuesta"]
    top5_df.to_csv(f"{DIR_OUTPUT}/4_Top5_Interacciones.csv", index=False)
    
    print("\n" + "="*60)
    print(f"¡Proceso finalizado! Las 4 gráficas se guardaron en:\n{DIR_OUTPUT}")
    print("="*60)

if __name__ == "__main__":
    main()
