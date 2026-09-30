#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    parser = argparse.ArgumentParser(description="Estadísticas Generales y Top 20 Taxones de EMU (HAC vs SUP)")
    parser.add_argument('--dir_hac', required=True, help="Directorio con los resultados de HAC")
    parser.add_argument('--dir_sup', required=True, help="Directorio con los resultados de SUP")
    parser.add_argument('-o', '--out_dir', default='.', help="Directorio de salida para gráficas y tablas")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # =====================================================================
    # 1. BÚSQUEDA Y CARGA DE ARCHIVOS
    # =====================================================================
    def get_filepath(directory, base_name, suffix):
        path_suffix = os.path.join(directory, f"{base_name}_{suffix}.tsv")
        path_normal = os.path.join(directory, f"{base_name}.tsv")
        if os.path.exists(path_suffix):
            return path_suffix
        return path_normal

    file_hac_rel = get_filepath(args.dir_hac, 'tabla_abundancia_relativa', 'hac')
    file_sup_rel = get_filepath(args.dir_sup, 'tabla_abundancia_relativa', 'sup')
    file_hac_counts = get_filepath(args.dir_hac, 'tabla_conteos', 'hac')
    file_sup_counts = get_filepath(args.dir_sup, 'tabla_conteos', 'sup')
    
    print("Cargando archivos...")
    try:
        df_hac_rel = pd.read_csv(file_hac_rel, sep='\t', index_col=0)
        df_sup_rel = pd.read_csv(file_sup_rel, sep='\t', index_col=0)
        df_hac_counts = pd.read_csv(file_hac_counts, sep='\t', index_col=0)
        df_sup_counts = pd.read_csv(file_sup_counts, sep='\t', index_col=0)
    except FileNotFoundError as e:
        print(f"Error cargando archivos: {e}")
        return

    # =====================================================================
    # 2. CÁLCULO DE ESTADÍSTICAS GENERALES (TABLA SIMPLIFICADA)
    # =====================================================================
    def calc_stats(df_rel, df_counts, name):
        n_taxa = df_rel.shape[0]
        n_samples = df_rel.shape[1]
        
        total_reads = df_counts.sum().sum()
        mean_reads_per_sample = df_counts.sum(axis=0).mean()
        
        # Taxones con abundancia relativa promedio mayor al 1%
        mean_rel = df_rel.mean(axis=1)
        if df_rel.max().max() <= 1.1:
            abundant_taxa = (mean_rel > 0.01).sum()
        else:
            abundant_taxa = (mean_rel > 1.0).sum()
        
        return {
            "Modelo": name,
            "Taxones Detectados": n_taxa,
            "Muestras Analizadas": n_samples,
            "Total Lecturas Asignadas": f"{int(total_reads):,}".replace(",", "."),
            "Promedio Lecturas por Muestra": f"{int(mean_reads_per_sample):,}".replace(",", "."),
            "Taxones Dominantes (>1%)": abundant_taxa
        }

    stats = [
        calc_stats(df_hac_rel, df_hac_counts, "EMU - HAC"),
        calc_stats(df_sup_rel, df_sup_counts, "EMU - SUP")
    ]
    
    df_stats = pd.DataFrame(stats)
    
    # Exportar tabla lista para Excel
    output_csv = os.path.join(args.out_dir, "estadisticas_generales_modelos.csv")
    df_stats.to_csv(output_csv, sep=';', index=False, encoding='utf-8-sig')
    print(f"¡Tabla de estadísticas guardada con éxito como '{output_csv}'!\n")

    # =====================================================================
    # 3. GRÁFICOS DE BARRAS (TOP 20 TAXONES)
    # =====================================================================
    def plot_top20(df_rel, title, filename, color):
        mean_abund = df_rel.mean(axis=1).sort_values(ascending=False).head(20)
        
        if df_rel.max().max() <= 1.1:
            mean_abund = mean_abund * 100
            
        plt.figure(figsize=(10, 7))
        sns.barplot(x=mean_abund.values, y=mean_abund.index, color=color)
        plt.title(f"Top 20 Taxones Más Abundantes - {title}", fontsize=14, fontweight='bold', pad=15)
        plt.xlabel("Abundancia Relativa Promedio (%)", fontsize=12, fontweight='bold')
        plt.ylabel("", fontsize=12)
        
        clean_labels = []
        for label in mean_abund.index:
            # 1. Separar por si viene con el linaje completo (con ;)
            parts = str(label).split(';')
            last_part = parts[-1].strip() if len(parts) > 0 else str(label)
            
            # 2. Separar por el pipe (|) para quitar "species|1234|" y dejar solo el nombre
            sub_parts = last_part.split('|')
            clean_name = sub_parts[-1].strip()
            
            clean_labels.append(clean_name)
        
        plt.yticks(ticks=range(len(clean_labels)), labels=clean_labels, fontsize=11, fontstyle='italic')
        
        for i, v in enumerate(mean_abund.values):
            plt.text(v + 0.1, i + 0.15, f"{v:.1f}%", color='black', fontsize=10)
            
        plt.tight_layout()
        plt.savefig(os.path.join(args.out_dir, filename), dpi=300)
        plt.close()

    print("Generando gráficos de barras con nombres limpios...")
    plot_top20(df_hac_rel, "Modelo HAC", "top20_hac.png", "#457b9d")
    plot_top20(df_sup_rel, "Modelo SUP", "top20_sup.png", "#e63946")
    print(f"¡Gráficos guardados exitosamente en el directorio: {args.out_dir}")

if __name__ == "__main__":
    main()
