import pandas as pd
import argparse
import os
import re
import glob


def build_phage_name_map(fasta_source):
    """Recorre uno o varios fasta con headers tipo:
    '>NC_024148.1 Mycobacterium phage Phantastic, complete genome'
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
                    accession = accession_full.split('.')[0]  # quita el .1 de version
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


def normalizar(nombre):
    """Normaliza para comparar sin importar mayusculas ni el typo weisella/weissella."""
    n = nombre.strip().lower()
    n = n.replace("weisella", "weissella")
    return n


def main():
    parser = argparse.ArgumentParser(
        description="Limpia nombres de fago/bacteria de una tabla de interacciones de comprobacion "
                    "(interacciones_comp.csv, interacciones_key_gene.csv o interacciones_wgs.csv) "
                    "y la filtra a solo ciertos taxones bacterianos, como filtrar_top_taxones.py "
                    "pero operando sobre IDs crudos (taxid_/GENERO_/accession) en vez de nombres ya legibles."
    )
    parser.add_argument('--input', type=str, required=True,
                         help="CSV de entrada con columnas Source (accession de fago), Target (bact_id crudo) y Weight.")
    parser.add_argument('--output', type=str, required=True, help="CSV de salida filtrado y con nombres legibles.")
    parser.add_argument('--taxones', type=str, required=True,
                         help="Lista de taxones separados por coma, ej: 'Lactococcus lactis,Weissella soli'")
    parser.add_argument('--phage_fasta', type=str, default=None,
                         help="Carpeta con fasta de fagos (ej. comprobacion/datos_phage) para reemplazar el "
                              "accession del Source por su nombre descriptivo. Si se omite, se deja el accession tal cual.")
    args = parser.parse_args()

    taxones_normalizados = set(normalizar(t) for t in args.taxones.split(','))

    df = pd.read_csv(args.input)

    df['Target'] = df['Target'].apply(clean_bacteria_name)

    if args.phage_fasta:
        phage_names = build_phage_name_map(args.phage_fasta)
        sin_nombre = df[~df['Source'].isin(phage_names.keys())]['Source'].unique()
        df['Source'] = df['Source'].apply(lambda x: phage_names.get(x, x))
        if len(sin_nombre) > 0:
            print(f"  [!] {len(sin_nombre)} fagos no encontraron nombre en el fasta (quedaron con su accession). "
                  f"Ejemplos: {list(sin_nombre)[:10]}")

    # Si dos accessions distintos terminan con el mismo nombre de fago, o el mismo bact_id
    # se limpia igual, podrian fusionarse en un solo par Source-Target. Nos quedamos con el
    # Weight maximo de cada grupo.
    duplicados = df.duplicated(subset=['Source', 'Target'], keep=False)
    if duplicados.any():
        n_dup_grupos = df[duplicados].groupby(['Source', 'Target']).ngroups
        print(f"  [!] Aviso: {n_dup_grupos} pares Source-Target quedaron duplicados tras renombrar. "
              f"Se conserva la fila con el Weight maximo de cada grupo.")
        df = df.sort_values('Weight', ascending=False).drop_duplicates(subset=['Source', 'Target'], keep='first')

    df_filtrado = df[df['Target'].apply(normalizar).isin(taxones_normalizados)].copy()
    df_filtrado.to_csv(args.output, index=False)

    print(f"Filtrado {args.input} -> {args.output}")
    print(f"  {len(df)} aristas totales -> {len(df_filtrado)} aristas retenidas")
    print(f"  Taxones encontrados en los datos: {sorted(df_filtrado['Target'].unique())}")

    faltantes = taxones_normalizados - set(df_filtrado['Target'].apply(normalizar).unique())
    if faltantes:
        print(f"  [!] Taxones pedidos que NO aparecieron en esta combinacion: {sorted(faltantes)}")


if __name__ == "__main__":
    main()
