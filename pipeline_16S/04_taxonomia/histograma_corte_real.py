import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

def main():
    # =====================================================================
    # 1. CARGA Y PROCESAMIENTO DE DATOS REALES
    # =====================================================================
    input_file = '/home/fenrir/Documentos/Tesis/results_dentrim_Q9/EMU_Q9/feature_table_relabund.tsv'
    
    print(f"Cargando datos desde: {input_file}")
    # Cargar matriz (taxones en filas, muestras en columnas)
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    
    # Calcular la abundancia relativa promedio de cada taxón en todas las muestras
    abundancias_promedio = df.mean(axis=1).values
    
    # Excluir ceros absolutos (si un taxón tiene 0.0 en todas las muestras) 
    abundancias_promedio = abundancias_promedio[abundancias_promedio > 0]
    
    umbral = 0.001
    taxones_retenidos = np.sum(abundancias_promedio >= umbral)
    taxones_excluidos = np.sum(abundancias_promedio < umbral)
    total_taxones = len(abundancias_promedio)

    print(f"Total taxones procesados: {total_taxones}")
    print(f"Taxones retenidos (>= 0.001): {taxones_retenidos}")
    print(f"Taxones excluidos (< 0.001): {taxones_excluidos}")

    # =====================================================================
    # 2. CONFIGURACIÓN DEL GRÁFICO (Estilo Académico APA)
    # =====================================================================
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Bins en escala logarítmica
    bins = np.logspace(np.log10(abundancias_promedio.min()), np.log10(abundancias_promedio.max()), 55)
    
    # Dibujar histograma
    n, bins_edges, patches = ax.hist(
        abundancias_promedio, bins=bins, edgecolor='black', linewidth=0.5
    )
    
    # Colorear según la posición respecto al umbral
    for i in range(len(patches)):
        if bins_edges[i] >= umbral:
            patches[i].set_facecolor('#457b9d') # Azul (Retenidos)
            patches[i].set_alpha(0.9)
        else:
            patches[i].set_facecolor('#cccccc') # Gris (Excluidos)
            patches[i].set_alpha(0.7)
            
    # Trazar la línea de corte
    ax.axvline(x=umbral, color='#e63946', linestyle='--', linewidth=2.5)
    
    # Ajuste dinámico del texto del umbral para que no choque con las barras
    y_max = ax.get_ylim()[1]
    ax.text(umbral * 1.2, y_max * 0.85, 'Umbral de Corte\n(0.001)', 
            color='#e63946', fontsize=11, fontweight='bold')
    
    ax.set_xscale('log')
    
    # Leyendas con los datos extraídos del archivo
    legend_elements = [
        Patch(facecolor='#cccccc', edgecolor='black', alpha=0.7, 
              label=f'Taxones poco abundantes ({taxones_excluidos} taxones)'),
        Patch(facecolor='#457b9d', edgecolor='black', alpha=0.9, 
              label=f'Taxones más dominantes ({taxones_retenidos} taxones)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=11)
    
    # Ejes
    ax.set_xlabel("Abundancia Relativa Promedio (Escala Logarítmica)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Número de Taxones", fontsize=12, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=11)
    
    # Guardar en la misma ruta de EMU
    output_img = "/home/fenrir/Documentos/Tesis/results_dentrim_Q9/EMU_Q9/histograma_abundancia_real.png"
    plt.tight_layout()
    plt.savefig(output_img, dpi=300)
    print(f"\nGráfico guardado con éxito en: {output_img}")
    plt.close()

if __name__ == "__main__":
    main()
