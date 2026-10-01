#!/usr/bin/env python3
import argparse
import os
import re
import sys
import warnings
import itertools
from pathlib import Path
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

def imprimir_error(*a, **k):
    print(*a, file=sys.stderr, **k)

def cargar_mapeo_doble(ruta_puente: Path, ruta_meta: Path, columna_grupo: str) -> dict:
    """
    Barcode(int) -> Sistema (o el valor de columna_grupo que corresponda), usando
    el mismo puente CSV (Mapa Barcodes Microbioma.csv) + maestro (Sistemas
    Agrícolas y Muestras.xlsx) que rarefaccion.py / taxonomy_profiling.py.
    """
    mapeo = {}
    df_puente = pd.read_csv(ruta_puente, sep=';')
    if not {"ID", "Barcode"}.issubset(df_puente.columns):
        df_puente = pd.read_csv(ruta_puente, sep=',')
    barcode_a_id = {}
    for _, fila in df_puente.iterrows():
        try:
            barcode_a_id[int(fila["Barcode"])] = str(fila["ID"]).strip()
        except (ValueError, TypeError):
            continue

    df_meta = pd.read_excel(ruta_meta, sheet_name=0, engine='openpyxl').dropna(subset=['ID'])
    id_a_grupo = dict(zip(df_meta['ID'].astype(str).str.strip(), df_meta[columna_grupo].astype(str).str.strip()))

    for num_bc, id_finca in barcode_a_id.items():
        mapeo[num_bc] = id_a_grupo.get(id_finca, "Desconocido")
    return mapeo

def num_barcode_desde_muestra(muestra: str):
    m = re.search(r'barcode0*(\d+)', muestra.lower())
    return int(m.group(1)) if m else None

def main():
    base = Path("/home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8")
    emu_propio = base / "EMU_propio"

    parser = argparse.ArgumentParser(description="Script de Tesis: Diversidad con Traducción Doble (Barcode -> ID -> Sistema)")
    parser.add_argument('--model', choices=['hac', 'sup'], default='hac',
                        help="Modelo de basecalling: define los valores por defecto de --input y --output_dir.")
    parser.add_argument('--dataset', choices=['todo', 'results'], default='todo',
                        help="'todo' -> EMU{model}_todo; 'results' -> EMU{model}_results (>500 lecturas).")
    parser.add_argument('-i', '--input', default=None, help="Ruta a la tabla de abundancias relativas (TSV). Por defecto: EMU{model}_{dataset}/feature_table_relabund_{model}_{dataset}.tsv")
    parser.add_argument('-m', '--metadata', default=emu_propio / "Diversidad" / "Sistemas Agrícolas y Muestras.xlsx",
                        help="Ruta al archivo de metadatos (Sistemas Agrícolas y Muestras).")
    parser.add_argument('-b', '--barcode_map', default=base / "Mapa Barcodes Microbioma.csv",
                        help="Ruta al archivo puente barcode->ID (Mapa Barcodes Microbioma.csv).")
    parser.add_argument('-o', '--output_dir', default=None, help="Directorio de salida. Por defecto: Diversidad/{model}_{dataset}")
    parser.add_argument('-p', '--prefix', default=None, help="Prefijo de salida. Por defecto: Tesis_{model}_{dataset}")
    parser.add_argument('-c', '--group_col', default="Sistema", help="Columna de tratamientos (por defecto: Sistema).")
    parser.add_argument('-g', '--groups', nargs='+', default=["Empresarial", "Campesina", "Agroecológica"],
                        help="Grupos a comparar (por defecto, los 3 Sistemas).")
    args = parser.parse_args()

    carpeta_emu = emu_propio / f"EMU{args.model}_{args.dataset}"
    if args.input is None:
        args.input = carpeta_emu / f"feature_table_relabund_{args.model}_{args.dataset}.tsv"
    if args.output_dir is None:
        args.output_dir = emu_propio / "Diversidad" / f"{args.model}_{args.dataset}"
    if args.prefix is None:
        args.prefix = f"Tesis_{args.model}_{args.dataset}"

    imprimir_error(f"[INFO] Modelo: {args.model} | dataset: {args.dataset}")
    imprimir_error(f"[INFO] input: {args.input}")
    imprimir_error(f"[INFO] output_dir: {args.output_dir}")

    os.makedirs(args.output_dir, exist_ok=True)
    PREFIJO_SALIDA = os.path.join(args.output_dir, args.prefix)

    # 1. CARGA DE TABLA DE ABUNDANCIA
    print(f"Cargando tabla de abundancia: {args.input}")
    datos_crudos = pd.read_csv(args.input, index_col=0, sep='\t')
    columnas_resumen = [c for c in datos_crudos.columns if str(c).startswith(("Total_counts", "Frecuencia_"))]
    if columnas_resumen:
        datos_crudos = datos_crudos.drop(columns=columnas_resumen)
    abundancia = datos_crudos.T
    conteos = abundancia.copy()
    nombres_muestras = list(conteos.index)

    # 2 y 3. PUENTE (Barcode -> ID) + METADATOS (ID -> Sistema), en un solo paso
    print(f"Cargando puente ({args.barcode_map}) y metadatos ({args.metadata})...")
    barcode_a_grupo = cargar_mapeo_doble(Path(args.barcode_map), Path(args.metadata), args.group_col)

    # 4. TRADUCCIÓN POR MUESTRA (Barcode -> Sistema)
    print("Traduciendo Barcodes a Sistemas...")
    GRUPOS = {}
    sin_mapeo = []
    for muestra in nombres_muestras:
        num_bc = num_barcode_desde_muestra(muestra)
        grupo = barcode_a_grupo.get(num_bc) if num_bc is not None else None
        if grupo is None:
            grupo = "Desconocido"
            sin_mapeo.append(muestra)
        GRUPOS[muestra] = grupo

    if sin_mapeo:
        imprimir_error(f"[WARN] {len(sin_mapeo)} muestra(s) sin Sistema resuelto: {sin_mapeo}")

    print("\n--- Conteo verificado de muestras asignadas por grupo ---")
    serie_conteos = pd.Series(GRUPOS.values())
    print(serie_conteos.value_counts())
    print("----------------------------------------------------------\n")

    metadatos = pd.DataFrame(index=abundancia.index)
    metadatos['Group'] = metadatos.index.map(GRUPOS)

    # =========================================================================
    # MULTI-ANÁLISIS 1: ALFA DIVERSIDAD (SHANNON Y OBSERVADOS)
    # =========================================================================
    print("Calculando índices de Alfa Diversidad...")
    shannon = alpha_diversity('shannon', conteos.values, ids=conteos.index)
    observados = (conteos.values > 0).sum(axis=1)

    df_alfa = pd.DataFrame({
        'Sample': conteos.index,
        'Shannon': shannon.values,
        'Observados': observados,
        'Group': [GRUPOS[s] for s in conteos.index]
    })
    df_alfa.to_csv(f'{PREFIJO_SALIDA}_Alfa_Diversidad_Indices.tsv', sep='\t', index=False)

    df_alfa_filtrado = df_alfa[df_alfa['Group'].isin(args.groups)]

    # Pruebas estadísticas automáticas
    resultados_wilcoxon = {}
    print(f"Calculando pruebas de Wilcoxon para las combinaciones de {len(args.groups)} grupos...")
    for metrica in ['Shannon', 'Observados']:
        resultados_metrica = {}
        for nombre_g1, nombre_g2 in itertools.combinations(args.groups, 2):
            grupo1 = df_alfa_filtrado[df_alfa_filtrado['Group'] == nombre_g1]
            grupo2 = df_alfa_filtrado[df_alfa_filtrado['Group'] == nombre_g2]
            if len(grupo1) > 1 and len(grupo2) > 1:
                estadistico, valor_p = mannwhitneyu(grupo1[metrica], grupo2[metrica], alternative='two-sided')
                resultados_metrica[f"{nombre_g1} vs {nombre_g2}"] = valor_p
        resultados_wilcoxon[metrica] = resultados_metrica

    df_valores_p = pd.DataFrame(resultados_wilcoxon)
    df_valores_p.to_csv(f'{PREFIJO_SALIDA}_Alfa_Wilcoxon_pvalues.tsv', sep='\t')

    # Boxplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    sns.set_theme(style="whitegrid")
    for i, metrica in enumerate(['Shannon', 'Observados']):
        sns.boxplot(data=df_alfa_filtrado, x='Group', y=metrica, hue='Group', palette='Set1', ax=axes[i], legend=False)
        sns.stripplot(data=df_alfa_filtrado, x='Group', y=metrica, color='black', alpha=0.6, size=6, ax=axes[i])
        axes[i].set_title(f'{metrica}')
        axes[i].set_xlabel('Sistema de manejo')
        axes[i].tick_params(axis='x', rotation=15)
    plt.tight_layout()
    plt.savefig(f'{PREFIJO_SALIDA}_Alfa_Diversidad_Boxplots.png', dpi=300)
    plt.close()

    # =========================================================================
    # MULTI-ANÁLISIS 2: DISTANCIA BETA Y HEATMAPS
    # =========================================================================
    print("Calculando matrices de distancia Beta...")
    distancias_bray = pdist(abundancia.values, metric='braycurtis')
    distancias_jaccard = pdist(abundancia.values, metric='jaccard')

    dm_bray = DistanceMatrix(squareform(distancias_bray), ids=abundancia.index)
    dm_jaccard = DistanceMatrix(squareform(distancias_jaccard), ids=abundancia.index)

    df_bray = pd.DataFrame(dm_bray.data, index=abundancia.index, columns=abundancia.index)
    df_jaccard = pd.DataFrame(dm_jaccard.data, index=abundancia.index, columns=abundancia.index)

    df_bray.to_csv(f'{PREFIJO_SALIDA}_braycurtis.tsv', sep='\t')
    df_jaccard.to_csv(f'{PREFIJO_SALIDA}_jaccard.tsv', sep='\t')

    for matriz, nombre in zip([df_bray, df_jaccard], ['BrayCurtis', 'Jaccard']):
        fig, ax = plt.subplots(figsize=(12, 10))
        sns.heatmap(matriz, cmap='viridis', xticklabels=True, yticklabels=True, ax=ax,
                    cbar_kws={'label': f'Distancia de {nombre}'})
        ax.tick_params(axis='both', which='major', labelsize=6)
        plt.title(f'Matriz de Distancia - {nombre}', fontsize=14, pad=15)
        plt.tight_layout()
        plt.savefig(f'{PREFIJO_SALIDA}_{nombre}_heatmap.png', dpi=300)
        plt.close()

    # =========================================================================
    # MULTI-ANÁLISIS 3: PCoA
    # =========================================================================
    print("Calculando PCoA...")
    resultado_pcoa = pcoa(dm_bray, number_of_dimensions=2)
    coordenadas = resultado_pcoa.samples.iloc[:, 0:2]
    coordenadas.columns = ['PC1', 'PC2']
    coordenadas['Group'] = coordenadas.index.map(GRUPOS)
    coordenadas.to_csv(f'{PREFIJO_SALIDA}_PCoA.tsv', sep='\t')

    coordenadas_grafico = coordenadas[coordenadas['Group'].isin(args.groups)]

    fig, ax = plt.subplots(figsize=(10, 7))
    sns.scatterplot(
        data=coordenadas_grafico, x='PC1', y='PC2', hue='Group',
        palette='Set1', s=100, edgecolor='black', zorder=5, ax=ax
    )

    import matplotlib.patches as patches
    manejadores, etiquetas = ax.get_legend_handles_labels()
    mapa_color = {etiqueta: manejador.get_color() for manejador, etiqueta in zip(manejadores, etiquetas)}

    todos_elipse_x, todos_elipse_y = [], []

    for nombre_grupo, datos_grupo in coordenadas_grafico.groupby('Group'):
        if len(datos_grupo) > 2:
            x, y = datos_grupo['PC1'].values, datos_grupo['PC2'].values
            media_x, media_y = np.mean(x), np.mean(y)
            covarianza = np.cov(x, y)
            valores_propios, vectores_propios = np.linalg.eig(covarianza)
            orden = valores_propios.argsort()[::-1]
            valores_propios, vectores_propios = valores_propios[orden], vectores_propios[:, orden]
            theta = np.degrees(np.arctan2(*vectores_propios[:, 0][::-1]))
            ancho, alto = 2 * 1.96 * np.sqrt(valores_propios)

            elipse = patches.Ellipse(
                xy=(media_x, media_y), width=ancho, height=alto, angle=theta,
                alpha=0.15, facecolor=mapa_color.get(nombre_grupo, 'gray'),
                edgecolor=mapa_color.get(nombre_grupo, 'black'),
                linestyle='--', linewidth=1.5, zorder=2
            )
            ax.add_patch(elipse)
            todos_elipse_x.extend([media_x - ancho/2, media_x + ancho/2])
            todos_elipse_y.extend([media_y - alto/2, media_y + alto/2])

    todos_elipse_x.extend(coordenadas_grafico['PC1'].tolist())
    todos_elipse_y.extend(coordenadas_grafico['PC2'].tolist())
    if todos_elipse_x and todos_elipse_y:
        x_min, x_max = min(todos_elipse_x), max(todos_elipse_x)
        y_min, y_max = min(todos_elipse_y), max(todos_elipse_y)
        relleno_x, relleno_y = (x_max - x_min) * 0.20, (y_max - y_min) * 0.20
        ax.set_xlim(x_min - relleno_x, x_max + relleno_x)
        ax.set_ylim(y_min - relleno_y, y_max + relleno_y)

    var_pc1 = round(resultado_pcoa.proportion_explained.iloc[0] * 100, 2)
    var_pc2 = round(resultado_pcoa.proportion_explained.iloc[1] * 100, 2)

    plt.xlabel(f"PC1 ({var_pc1}%)", fontsize=11)
    plt.ylabel(f"PC2 ({var_pc2}%)", fontsize=11)
    plt.title('Diversidad Beta: Ordenación PCoA Bray-Curtis', fontsize=13, pad=15)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Sistema de manejo')
    plt.tight_layout()
    plt.savefig(f'{PREFIJO_SALIDA}_PCoA.png', dpi=300)
    plt.close()

    # =========================================================================
    # MULTI-ANÁLISIS 4: PERMANOVA
    # =========================================================================
    metadatos_filtrados = metadatos[metadatos['Group'].isin(args.groups)]

    print("Calculando PERMANOVA...")
    dm_bray_filtrada = dm_bray.filter(metadatos_filtrados.index)
    resultado_permanova = permanova(
        distance_matrix=dm_bray_filtrada,
        grouping=metadatos_filtrados['Group'],
        permutations=999,
        seed=42
    )
    with open(f'{PREFIJO_SALIDA}_PERMANOVA.txt', 'w') as f:
        f.write(str(resultado_permanova))

    print("\n--- RESULTADOS DEL PERMANOVA ---")
    print(resultado_permanova)
    print(f"\n¡Métricas generadas con éxito en: {args.output_dir}")

if __name__ == '__main__':
    main()
