import pandas as pd

def clean_labels(text):
    if pd.isna(text):
        return "Unknown"
    # Limpieza profunda: puntos por pipes, guiones bajos por espacios y limpieza de comillas
    clean = str(text).replace('|', '.').replace(' ', '_').replace('-', '_').replace('"', '')
    return clean

def process_thesis_data(counts_path, tax_path, output_prefix):
    print(f"--- Procesando datos de EMU: {counts_path} ---")
    
    # 1. Cargar tablas (usando sep='\t' porque son .tsv)
    counts = pd.read_csv(counts_path, sep='\t', index_col=0)
    tax = pd.read_csv(tax_path, sep='\t', index_col=0)

    # 2. Limpiar los IDs de los taxones en ambas tablas
    counts.index = [clean_labels(x) for x in counts.index]
    tax.index = [clean_labels(x) for x in tax.index]

    # 3. Limpiar el contenido de las categorías taxonómicas
    for col in tax.columns:
        tax[col] = tax[col].apply(clean_labels)

    # 4. Asegurar que no haya taxones huerfanos
    common_ids = counts.index.intersection(tax.index)
    counts_final = counts.loc[common_ids]
    tax_final = tax.loc[common_ids]

    # 5. Exportar a CSV (formato más amigable para R)
    counts_final.to_csv(f"{output_prefix}_counts_cleaned.csv")
    tax_final.to_csv(f"{output_prefix}_tax_cleaned.csv")
    
    print(f"¡Listo! {len(common_ids)} taxones procesados.")
    print(f"Nuevos archivos: {output_prefix}_counts_cleaned.csv y {output_prefix}_tax_cleaned.csv")

# --- EJECUCIÓN ---
process_thesis_data(
    counts_path='feature_table_counts.tsv', 
    tax_path='taxonomy.tsv', 
    output_prefix='Tesis_2020'
)
