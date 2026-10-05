# Filtra las predicciones de DeepPBI-KG por umbral y saca tres tablas para Gephi:
# genes clave, genoma completo y el promedio de las dos

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
    
    # vacíos a 0, si no el promedio falla cuando no hay genes clave
    df['wgs_output'] = pd.to_numeric(df['wgs_output'], errors='coerce').fillna(0)
    df['key_gene_output'] = pd.to_numeric(df['key_gene_output'], errors='coerce').fillna(0)

    # genes clave
    df_kg = df[df['key_gene_output'] >= args.threshold].copy()
    df_kg.rename(columns={'phage': 'Source', 'bacterium': 'Target', 'key_gene_output': 'Weight'}, inplace=True)
    df_kg[['Source', 'Target', 'Weight']].to_csv(os.path.join(args.output_dir, 'interacciones_key_gene.csv'), index=False)
    print(f"[Key Gene] Generado: {len(df_kg)} interacciones.")

    # genoma completo
    df_wgs = df[df['wgs_output'] >= args.threshold].copy()
    df_wgs.rename(columns={'phage': 'Source', 'bacterium': 'Target', 'wgs_output': 'Weight'}, inplace=True)
    df_wgs[['Source', 'Target', 'Weight']].to_csv(os.path.join(args.output_dir, 'interacciones_wgs.csv'), index=False)
    print(f"[WGS] Generado: {len(df_wgs)} interacciones.")

    # compuesta (promedio de las dos)
    df['Weight'] = (df['wgs_output'] + df['key_gene_output']) / 2.0
    
    df_comp = df[df['Weight'] >= args.threshold].copy()
    df_comp.rename(columns={'phage': 'Source', 'bacterium': 'Target'}, inplace=True)
    
    # se dejan también las dos salidas por separado
    cols_comp = ['Source', 'Target', 'key_gene_output', 'wgs_output', 'Weight']
    df_comp[cols_comp].to_csv(os.path.join(args.output_dir, 'interacciones_comp.csv'), index=False)
    print(f"[Compuesta/Promedio] Generado: {len(df_comp)} interacciones líticas ecológicamente viables.")

if __name__ == '__main__':
    main()
