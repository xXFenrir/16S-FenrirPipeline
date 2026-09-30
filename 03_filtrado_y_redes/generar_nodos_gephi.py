import pandas as pd
import argparse
import os
import re
import glob

def build_phage_name_map(fasta_source):
    """Recorre uno o varios fasta con headers tipo:
    '>GQ303259.1 Mycobacterium phage Colbert, complete genome'
    y devuelve un diccionario: accession_sin_version -> nombre descriptivo."""
    name_map = {}
    if os.path.isdir(fasta_source):
        files = glob.glob(os.path.join(fasta_source, "*.fasta")) + glob.glob(os.path.join(fasta_source, "*.fna"))
    else:
        files = [fasta_source]

    header_re = re.compile(r'^>(\S+)\s+(.*)$')
    for fpath in files:
        with open(fpath, 'r', errors='ignore') as f:
            for line in f:
                if line.startswith('>'):
                    m = header_re.match(line.strip())
                    if not m:
                        continue
                    accession_full, desc = m.groups()
                    accession = accession_full.split('.')[0]  # quita el .1 de versión
                    desc_clean = re.split(
                        r',\s*(complete genome|complete sequence|partial genome|partial sequence|genome)',
                        desc, flags=re.IGNORECASE
                    )[0]
                    name_map[accession] = desc_clean.strip()
    return name_map

def clean_bacteria_name(bact_id):
    """Convierte 'taxid_1358_Lactococcus_lactis' -> 'Lactococcus lactis'
    'GENERO_Pseudonocardia_882449' -> 'Pseudonocardia'
    'name_Algo_asi' -> 'Algo asi'"""
    m = re.match(r'^taxid_\d+_(.+)$', bact_id)
    if m:
        return m.group(1).replace('_', ' ')
    m = re.match(r'^GENERO_([A-Za-z0-9]+)_\d+$', bact_id)
    if m:
        return m.group(1)
    m = re.match(r'^name_(.+)$', bact_id)
    if m:
        return m.group(1).replace('_', ' ')
    return bact_id

def main():
    parser = argparse.ArgumentParser(description="Genera un CSV Source/Target/Weight con nombres legibles, listo para Gephi.")
    parser.add_argument('--edges', type=str, required=True,
                         help="CSV de aristas original (Source=fago, Target=bacteria, Weight=probabilidad).")
    parser.add_argument('--phage_fasta', type=str, required=True,
                         help="Ruta a all_phage_seq.fasta O a la carpeta fagos_fna con archivos individuales.")
    parser.add_argument('--output', type=str, required=True, help="Ruta de salida del CSV final.")
    args = parser.parse_args()

    df = pd.read_csv(args.edges)
    phage_names = build_phage_name_map(args.phage_fasta)

    df_final = pd.DataFrame()
    df_final['Source'] = df['Source'].apply(lambda x: phage_names.get(x, x))
    df_final['Target'] = df['Target'].apply(clean_bacteria_name)
    df_final['Weight'] = df['Weight']

    n_antes = len(df_final)

    # Si dos IDs distintos terminan con el mismo nombre (ej. dos cepas de la misma especie),
    # se fusionarían en un solo nodo en Gephi. Avisamos y nos quedamos con el peso máximo.
    duplicados = df_final.duplicated(subset=['Source', 'Target'], keep=False)
    if duplicados.any():
        n_dup_grupos = df_final[duplicados].groupby(['Source', 'Target']).ngroups
        print(f"[!] Aviso: {n_dup_grupos} pares Source-Target quedaron duplicados tras renombrar "
              f"(probablemente organismos distintos con el mismo nombre legible). "
              f"Se conserva el Weight máximo de cada grupo.")

    df_final = df_final.groupby(['Source', 'Target'], as_index=False)['Weight'].max()

    df_final.to_csv(args.output, index=False)

    sin_nombre = df[~df['Source'].isin(phage_names.keys())]['Source'].unique()
    print(f"Generado {args.output}: {n_antes} aristas originales -> {len(df_final)} aristas finales "
          f"({df_final['Source'].nunique()} fagos, {df_final['Target'].nunique()} bacterias).")
    if len(sin_nombre) > 0:
        print(f"[!] {len(sin_nombre)} fagos no encontraron nombre en el fasta (quedaron con su ID). "
              f"Ejemplos: {list(sin_nombre)[:10]}")

if __name__ == "__main__":
    main()
