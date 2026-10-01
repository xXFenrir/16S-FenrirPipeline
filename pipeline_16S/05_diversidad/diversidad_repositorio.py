#!/usr/bin/env python3
import argparse
import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*global interpreter lock.*")

from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix
from skbio.diversity import alpha_diversity

def main():
    parser = argparse.ArgumentParser(description="Script de Tesis: Diversidad Alfa y Beta (Gráficos Optimizados)")
    parser.add_argument(
        '-i', '--input', 
        default='/home/fenrir/Documentos/Tesis/results_dentrim_Q9/EMU_Q9/feature_table_relabund.tsv', 
        help="Ruta al archivo de la tabla de abundancias relativas."
    )
    parser.add_argument(
        '-m', '--metadata',
        default='/home/fenrir/Documentos/Tesis/Muestras_16S/SraRunTable.csv',
        help="SraRunTable del BioProject (columnas 'Run' y 'crop_rotation')."
    )
    parser.add_argument(
        '-o', '--output_dir', 
        default='/home/fenrir/Documentos/Tesis/results_dentrim_Q9/EMU_Q9/estad', 
        help="Directorio de salida de resultados."
    )
    parser.add_argument('-p', '--prefix', default="Tesis_2020", help="Prefijo de salida.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    OUTPUT_PREFIX = os.path.join(args.output_dir, args.prefix)

    # 1. CARGA DE TABLA DE ABUNDANCIA DE EMU Y TRANSPOSICIÓN
    print(f"Cargando tabla de abundancia de EMU: {args.input}")
    raw_data = pd.read_csv(args.input, index_col=0, sep='\t')
    abundance = raw_data.T
    counts = abundance.copy()
    srr_samples = list(counts.index)

    # 2. CARGA DE METADATOS (SraRunTable)
    print(f"Cargando metadatos base desde: {args.metadata}")
    df_meta = pd.read_csv(args.metadata)
    df_meta['Run'] = df_meta['Run'].astype(str).str.strip()
    df_meta['crop_rotation'] = df_meta['crop_rotation'].astype(str).str.strip()
    run_a_rotacion = dict(zip(df_meta['Run'], df_meta['crop_rotation']))

    # 3. ASIGNACIÓN DE TRATAMIENTOS POR IDENTIFICADOR SRR (no por posición)
    print("Asignando tratamientos por identificador SRR...")
    NOMBRES_GRUPO = {"CSCS": "CS", "CS": "CS", "CSSWP": "CSSwP", "NOT APPLICABLE": "Control"}
    GROUPS = {}
    for srr_id in srr_samples:
        rotation_val = run_a_rotacion.get(srr_id)
        GROUPS[srr_id] = NOMBRES_GRUPO.get(rotation_val.upper(), rotation_val) if rotation_val else "Sin metadata"

    print("\n--- Conteo de muestras asignadas por grupo ---")
    counts_series = pd.Series(GROUPS.values())
    print(counts_series.value_counts())
    print("----------------------------------------------\n")

    # Guardar matriz de metadatos consolidada para skbio
    metadata = pd.DataFrame(index=abundance.index)
    metadata['Group'] = metadata.index.map(GROUPS)

    # =========================================================================
    # MULTI-ANÁLISIS 1: ALFA DIVERSIDAD (SHANNON Y OBSERVADOS)
    # =========================================================================
    print("Calculando índices de Alfa Diversidad...")
    shannon = alpha_diversity('shannon', counts.values, ids=counts.index)
    observados = (counts.values > 0).sum(axis=1)

    df_alpha = pd.DataFrame({
        'Sample': counts.index,
        'Shannon': shannon.values,
        'Observados': observados,
        'Group': [GROUPS[s] for s in counts.index]
    })
    df_alpha.to_csv(f'{OUTPUT_PREFIX}_Alfa_Diversidad_Indices.tsv', sep='\t', index=False)

    df_w = df_alpha[df_alpha['Group'].isin(['CS', 'CSSwP'])]
    g1 = df_w[df_w['Group'] == 'CS']
    g2 = df_w[df_w['Group'] == 'CSSwP']

    wilcoxon_results = {}
    if len(g1) > 1 and len(g2) > 1:
        print("Calculando pruebas de Wilcoxon para Alfa Diversidad...")
        for metric in ['Shannon', 'Observados']:
            stat, p_val = mannwhitneyu(g1[metric], g2[metric], alternative='two-sided')
            wilcoxon_results[metric] = {'U_statistic': stat, 'p-value': p_val}
        
        df_wilcoxon = pd.DataFrame(wilcoxon_results).T
        df_wilcoxon.to_csv(f'{OUTPUT_PREFIX}_Alfa_Wilcoxon_Resultados.tsv', sep='\t')

    # Boxplots de Alfa Diversidad
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    sns.set_theme(style="whitegrid")
    for i, metric in enumerate(['Shannon', 'Observados']):
        sns.boxplot(data=df_w, x='Group', y=metric, hue='Group', palette='Set1', ax=axes[i], legend=False)
        sns.stripplot(data=df_w, x='Group', y=metric, color='black', alpha=0.6, size=6, ax=axes[i])
        p_str = f"p={round(wilcoxon_results[metric]['p-value'], 4)}" if metric in wilcoxon_results else "p=N/A"
        axes[i].set_title(f'{metric} ({p_str})')
        axes[i].set_xlabel('Rotación de cultivos')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_PREFIX}_Alfa_Diversidad_Boxplots.png', dpi=300)
    plt.close()

    # =========================================================================
    # MULTI-ANÁLISIS 2: DISTANCIA BETA Y HEATMAPS (LETRA PEQUEÑA CORREGIDA)
    # =========================================================================
    print("Calculando matrices de distancia Beta...")
    bray_distances = pdist(abundance.values, metric='braycurtis')
    jaccard_distances = pdist(abundance.values, metric='jaccard')
    
    bray_dm = DistanceMatrix(squareform(bray_distances), ids=abundance.index)
    jaccard_dm = DistanceMatrix(squareform(jaccard_distances), ids=abundance.index)

    bray_df = pd.DataFrame(bray_dm.data, index=abundance.index, columns=abundance.index)
    jaccard_df = pd.DataFrame(jaccard_dm.data, index=abundance.index, columns=abundance.index)

    bray_df.to_csv(f'{OUTPUT_PREFIX}_braycurtis.tsv', sep='\t')
    jaccard_df.to_csv(f'{OUTPUT_PREFIX}_jaccard.tsv', sep='\t')

    for matrix, name in zip([bray_df, jaccard_df], ['BrayCurtis', 'Jaccard']):
        fig, ax = plt.subplots(figsize=(12, 10))
        # Ajustamos cbar_kws para que la barra de color sea estéticamente proporcional
        sns.heatmap(matrix, cmap='viridis', xticklabels=True, yticklabels=True, ax=ax,
                    cbar_kws={'label': f'Distancia de {name}'})
        
        # AJUSTE CRUCIAL: Reducción del tamaño de la fuente para las 88 muestras
        ax.tick_params(axis='both', which='major', labelsize=5.5)
        
        plt.title(f'Matriz de Distancia - {name}', fontsize=14, pad=15)
        plt.tight_layout()
        plt.savefig(f'{OUTPUT_PREFIX}_{name}_heatmap.png', dpi=300)
        plt.close()

    # =========================================================================
    # MULTI-ANÁLISIS 3: PCoA CON ELIPSES COMPLETAS Y MÁRGENES EXPANDIDOS
    # =========================================================================
    print("Calculando PCoA...")
    pcoa_results = pcoa(bray_dm, number_of_dimensions=2)
    coords = pcoa_results.samples.iloc[:, 0:2]
    coords.columns = ['PC1', 'PC2']
    coords['Group'] = coords.index.map(GROUPS)
    coords.to_csv(f'{OUTPUT_PREFIX}_PCoA.tsv', sep='\t')

    coords_plot = coords[coords['Group'].isin(['CS', 'CSSwP'])]

    fig, ax = plt.subplots(figsize=(9, 7))
    sns.scatterplot(
        data=coords_plot, x='PC1', y='PC2', hue='Group', 
        palette='Set1', s=100, edgecolor='black', zorder=5, ax=ax
    )

    import matplotlib.patches as patches
    handles, labels = ax.get_legend_handles_labels()
    color_map = {label: handle.get_color() for handle, label in zip(handles, labels)}

    # Listas para rastrear la extensión total de las elipses y evitar que se corten
    all_ellipse_x = []
    all_ellipse_y = []

    for group_name, group_data in coords_plot.groupby('Group'):
        if len(group_data) > 2:
            x = group_data['PC1'].values
            y = group_data['PC2'].values
            
            # Centro de la elipse
            mean_x, mean_y = np.mean(x), np.mean(y)
            
            # Calcular matriz de covarianza y valores/vectores propios
            cov = np.cov(x, y)
            vals, vecs = np.linalg.eig(cov)
            
            # Ordenar de mayor a menor componente principal
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:, order]
            
            # Ángulo de rotación en grados
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            
            # Radios de la elipse usando el intervalo de confianza del 95% (1.96 * sqrt(val))
            width, height = 2 * 1.96 * np.sqrt(vals)
            
            ellipse = patches.Ellipse(
                xy=(mean_x, mean_y),
                width=width, height=height, angle=theta,
                alpha=0.15,
                facecolor=color_map.get(group_name, 'gray'),
                edgecolor=color_map.get(group_name, 'black'),
                linestyle='--', linewidth=1.5, zorder=2
            )
            ax.add_patch(ellipse)
            
            # Aproximar los límites de la elipse para ajustar los ejes dinámicamente
            all_ellipse_x.extend([mean_x - width/2, mean_x + width/2])
            all_ellipse_y.extend([mean_y - height/2, mean_y + height/2])

    # AJUSTE CRUCIAL: Expandir los límites del gráfico un 20% más allá de las elipses para que no se corten
    all_ellipse_x.extend(coords_plot['PC1'].tolist())
    all_ellipse_y.extend(coords_plot['PC2'].tolist())
    if all_ellipse_x and all_ellipse_y:
        x_min, x_max = min(all_ellipse_x), max(all_ellipse_x)
        y_min, y_max = min(all_ellipse_y), max(all_ellipse_y)
        x_padding = (x_max - x_min) * 0.20
        y_padding = (y_max - y_min) * 0.20
        
        ax.set_xlim(x_min - x_padding, x_max + x_padding)
        ax.set_ylim(y_min - y_padding, y_max + y_padding)

    pc1_var = round(pcoa_results.proportion_explained.iloc[0] * 100, 2)
    pc2_var = round(pcoa_results.proportion_explained.iloc[1] * 100, 2)

    plt.xlabel(f"PC1 ({pc1_var}%)", fontsize=11)
    plt.ylabel(f"PC2 ({pc2_var}%)", fontsize=11)
    plt.title('Diversidad Beta: Ordenación PCoA Bray-Curtis', fontsize=13, pad=15)
    ax.legend(title='Rotación de cultivos')
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_PREFIX}_PCoA.png', dpi=300)
    plt.close()

    # =========================================================================
    # MULTI-ANÁLISIS 4: PERMANOVA
    # =========================================================================
    metadata_filtered = metadata[metadata['Group'].isin(['CS', 'CSSwP'])]
    
    print("Calculando PERMANOVA...")
    bray_dm_filtered = bray_dm.filter(metadata_filtered.index)
    permanova_results = permanova(
        distance_matrix=bray_dm_filtered,
        grouping=metadata_filtered['Group'],
        permutations=999,
        seed=42
    )
    with open(f'{OUTPUT_PREFIX}_PERMANOVA.txt', 'w') as f:
        f.write(str(permanova_results))
    
    print("\n--- RESULTADOS DEL PERMANOVA ---")
    print(permanova_results)
    print(f"\n¡Imágenes optimizadas generadas con éxito en: {args.output_dir}")

if __name__ == '__main__':
    main()
