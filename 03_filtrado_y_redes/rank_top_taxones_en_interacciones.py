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
                    accession = accession_full.split('.')[0]
                    desc_clean = re.split(
                        r',\s*(complete genome|complete sequence|partial genome|partial sequence|genome)',
                        desc, flags=re.IGNORECASE
                    )[0]
                    name_map[accession] = desc_clean.strip()
    return name_map


def parse_bact_id(bact_id):
    """Devuelve (nombre_limpio, nivel) donde nivel es 'species' o 'genus'."""
    m = re.match(r'^taxid_\d+_(.+)$', bact_id)
    if m:
        return m.group(1).replace('_', ' '), 'species'
    m = re.match(r'^GENERO_([A-Za-z0-9]+)_\d+$', bact_id)
    if m:
        return m.group(1), 'genus'
    m = re.match(r'^name_(.+)$', bact_id)
    if m:
        return m.group(1).replace('_', ' '), 'species'
    return bact_id, 'species'


def normalizar(nombre):
    n = nombre.strip().lower()
    n = n.replace("weisella", "weissella")
    return n


def genero_de(especie_normalizada):
    """Extrae el genero del primer token, saltando 'candidatus' si aparece."""
    tokens = especie_normalizada.split()
    if not tokens:
        return especie_normalizada
    if tokens[0] == 'candidatus' and len(tokens) > 1:
        return tokens[1]
    return tokens[0]


def main():
    parser = argparse.ArgumentParser(
        description="Dentro de una tabla de interacciones (Source=fago, Target=bact_id crudo) ya filtrada por "
                    "umbral, selecciona los N taxones bacterianos MAS ABUNDANTES de entre los que SI aparecen en "
                    "esa tabla (usando la abundancia relativa de una feature_table de EMU), y filtra la tabla a "
                    "solo esos N taxones."
    )
    parser.add_argument('--interacciones', type=str, required=True,
                         help="CSV de interacciones (Source, Target=bact_id crudo, Weight, ...).")
    parser.add_argument('--feature_table', type=str, required=True,
                         help="feature_table_relabund.tsv de EMU (index tipo 'species|taxid|Genero especie').")
    parser.add_argument('--top_n', type=int, default=4, help="Cuantos taxones tomar (default 4).")
    parser.add_argument('--output', type=str, required=True, help="CSV de salida filtrado.")
    parser.add_argument('--phage_fasta', type=str, default=None,
                         help="Carpeta con fasta de fagos para reemplazar el accession del Source por su nombre.")
    args = parser.parse_args()

    # --- Abundancias desde la feature table (nivel especie y agregado a nivel genero) ---
    ft = pd.read_csv(args.feature_table, sep='\t', index_col=0)
    ft['mean_abundance'] = ft.mean(axis=1)

    especie_abund = {}
    genero_abund = {}
    for feature_id, row in ft.iterrows():
        parts = feature_id.split('|')
        nombre = parts[2] if len(parts) >= 3 else feature_id
        nombre_norm = normalizar(nombre)
        especie_abund[nombre_norm] = especie_abund.get(nombre_norm, 0.0) + row['mean_abundance']
        genero_abund[genero_de(nombre_norm)] = genero_abund.get(genero_de(nombre_norm), 0.0) + row['mean_abundance']

    # --- Cargar y limpiar la tabla de interacciones ---
    df = pd.read_csv(args.interacciones)
    parsed = df['Target'].apply(parse_bact_id)
    df['Target'] = parsed.apply(lambda x: x[0])
    df['_nivel'] = parsed.apply(lambda x: x[1])

    if args.phage_fasta:
        phage_names = build_phage_name_map(args.phage_fasta)
        df['Source'] = df['Source'].apply(lambda x: phage_names.get(x, x))

    duplicados = df.duplicated(subset=['Source', 'Target'], keep=False)
    if duplicados.any():
        n_dup_grupos = df[duplicados].groupby(['Source', 'Target']).ngroups
        print(f"  [!] Aviso: {n_dup_grupos} pares Source-Target quedaron duplicados tras renombrar. "
              f"Se conserva la fila con el Weight maximo de cada grupo.")
        df = df.sort_values('Weight', ascending=False).drop_duplicates(subset=['Source', 'Target'], keep='first')

    # --- Rankear los taxones que SI aparecen en la tabla de interacciones por su abundancia real ---
    taxones_en_tabla = df[['Target', '_nivel']].drop_duplicates()
    filas = []
    for _, r in taxones_en_tabla.iterrows():
        nombre_norm = normalizar(r['Target'])
        if r['_nivel'] == 'genus':
            abund = genero_abund.get(nombre_norm, 0.0)
        else:
            abund = especie_abund.get(nombre_norm)
            if abund is None:
                # Nombre a nivel especie que no aparece tal cual en la feature table:
                # probamos como genero (ej. si la feature table solo tiene otras especies del mismo genero).
                abund = 0.0
        filas.append({'Target': r['Target'], 'nivel': r['_nivel'], 'abundancia_media': abund})

    ranking = pd.DataFrame(filas).sort_values('abundancia_media', ascending=False)
    print("Ranking de abundancia de los taxones presentes en la tabla de interacciones:")
    print(ranking.to_string(index=False))

    top = ranking.head(args.top_n)
    taxones_top = set(normalizar(t) for t in top['Target'])
    print(f"\nTop {args.top_n} taxones seleccionados: {sorted(top['Target'])}")

    no_detectados = top[top['abundancia_media'] == 0.0]
    if len(no_detectados) > 0:
        print(f"  [!] Aviso: {len(no_detectados)} de los seleccionados tienen abundancia 0 en la feature table "
              f"(no se encontro coincidencia de nombre): {list(no_detectados['Target'])}")

    df_filtrado = df[df['Target'].apply(normalizar).isin(taxones_top)].drop(columns=['_nivel']).copy()
    df_filtrado.to_csv(args.output, index=False)

    print(f"\nFiltrado {args.interacciones} -> {args.output}")
    print(f"  {len(df)} aristas totales -> {len(df_filtrado)} aristas retenidas")


if __name__ == "__main__":
    main()
