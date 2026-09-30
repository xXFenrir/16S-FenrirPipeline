import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def main():
    sns.set_theme(style="whitegrid")
    
    # Datos extraídos de la tabla del artículo y de tu CSV
    categorias = ['Artículo Guía (2020)', 'Mis Secuencias (Tesis)']
    raw_reads = [3575783, 2365456]
    filtered_reads = [1850995, 1223953]
    
    # Cálculo de retención
    retencion = [(f / r) * 100 for f, r in zip(filtered_reads, raw_reads)]

    x = np.arange(len(categorias))
    width = 0.35  # Ancho de las barras

    fig, ax1 = plt.subplots(figsize=(9, 6))

    # Barras de Lecturas (Eje Y principal)
    bar1 = ax1.bar(x - width/2, raw_reads, width, label='Lecturas Iniciales (Crudas)', color='#1d3557', edgecolor='black')
    bar2 = ax1.bar(x + width/2, filtered_reads, width, label='Lecturas Post-Filtro (Limpias)', color='#2a9d8f', edgecolor='black')

    ax1.set_ylabel('Millones de Lecturas', fontsize=12, fontweight='bold')
    ax1.set_title('Eficacia del Filtrado de Calidad: Comparativa de Retención', fontsize=14, fontweight='bold', pad=15)
    ax1.set_xticks(x)
    ax1.set_xticklabels(categorias, fontsize=12, fontweight='bold')
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: f"{x*1e-6:.1f}M"))

    # Añadir números sobre las barras
    for bar in bar1 + bar2:
        height = bar.get_height()
        ax1.annotate(f'{int(height):,}'.replace(',', '.'),
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=10)

    # Línea de Porcentaje de Retención (Eje Y secundario)
    ax2 = ax1.twinx()
    ax2.plot(x, retencion, color='#e63946', marker='o', markersize=10, linewidth=3, linestyle='-', label='% de Retención')
    ax2.set_ylabel('Porcentaje de Retención (%)', fontsize=12, fontweight='bold', color='#e63946')
    ax2.tick_params(axis='y', labelcolor='#e63946')
    
    # Forzar el límite del eje secundario para que la línea se vea estable (50-55%)
    ax2.set_ylim(40, 60)

    # Anotar los porcentajes exactos en los puntos de la línea
    for i, txt in enumerate(retencion):
        ax2.annotate(f'{txt:.2f}%', (x[i], retencion[i]), 
                     textcoords="offset points", xytext=(0, 10), ha='center', 
                     fontsize=11, fontweight='bold', color='#e63946',
                     bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#e63946", lw=1))

    # Leyendas combinadas
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right', bbox_to_anchor=(1, 0.95))

    plt.tight_layout()
    plt.savefig('Comparativa_Retencion_Proporcional.png', dpi=300)
    print("¡Gráfico generado exitosamente como 'Comparativa_Retencion_Proporcional.png'!")

if __name__ == "__main__":
    main()
