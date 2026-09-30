import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Filtro triple de predicciones DeepPBI-KG para Gephi.")
    parser.add_argument('--input', type=str, required=True, help="Ruta del archivo CSV original.")
    parser.add_argument('--output_dir', type=str, required=True, help="Carpeta donde se guardarán los 3 archivos CSV.")
    parser.add_argument('--threshold', type=float, default=0.90, help="Umbral de probabilidad (ej. 0.90).")

    args = parser.parse_args()

    print(f"Cargando datos desde: {args.input}")
    df = pd.read_csv(args.input)
    
    # 1. Limpieza de datos: Convertir a numérico y reemplazar valores vacíos (NaN) con 0
    # Esto es vital para que el promedio matemático no falle si falta un gen clave
    df['wgs_output'] = pd.to_numeric(df['wgs_output'], errors='coerce').fillna(0)
    df['key_gene_output'] = pd.to_numeric(df['key_gene_output'], errors='coerce').fillna(0)

    # ---------------------------------------------------------
    # TABLA 1: INTERACCIONES KEY GENE
    # ---------------------------------------------------------
    df_kg = df[df['key_gene_output'] >= args.threshold].copy()
    df_kg.rename(columns={'phage': 'Source', 'bacterium': 'Target', 'key_gene_output': 'Weight'}, inplace=True)
    df_kg[['Source', 'Target', 'Weight']].to_csv(os.path.join(args.output_dir, 'interacciones_key_gene.csv'), index=False)
    print(f"[Key Gene] Generado: {len(df_kg)} interacciones.")

    # ---------------------------------------------------------
    # TABLA 2: INTERACCIONES WGS
    # ---------------------------------------------------------
    df_wgs = df[df['wgs_output'] >= args.threshold].copy()
    df_wgs.rename(columns={'phage': 'Source', 'bacterium': 'Target', 'wgs_output': 'Weight'}, inplace=True)
    df_wgs[['Source', 'Target', 'Weight']].to_csv(os.path.join(args.output_dir, 'interacciones_wgs.csv'), index=False)
    print(f"[WGS] Generado: {len(df_wgs)} interacciones.")

    # ---------------------------------------------------------
    # TABLA 3: INTERACCIONES COMPUESTAS (PROMEDIO)
    # ---------------------------------------------------------
    # Calculamos el promedio de ambas métricas
    df['Weight'] = (df['wgs_output'] + df['key_gene_output']) / 2.0
    
    df_comp = df[df['Weight'] >= args.threshold].copy()
    df_comp.rename(columns={'phage': 'Source', 'bacterium': 'Target'}, inplace=True)
    
    # Conservar el detalle individual más el promedio final
    cols_comp = ['Source', 'Target', 'key_gene_output', 'wgs_output', 'Weight']
    df_comp[cols_comp].to_csv(os.path.join(args.output_dir, 'interacciones_comp.csv'), index=False)
    print(f"[Compuesta/Promedio] Generado: {len(df_comp)} interacciones líticas ecológicamente viables.")

if __name__ == '__main__':
    main()
