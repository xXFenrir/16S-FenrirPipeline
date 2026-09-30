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
                    # Quita coletillas tipo ", complete genome" / ", partial sequence"
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
    return bact_id  # si no matchea ningún patrón, deja el id tal cual

def main():
    parser = argparse.ArgumentParser(description="Genera nodes.csv (Id, Label, Type) para Gephi con nombres legibles.")
    parser.add_argument('--edges', type=str, required=True,
                         help="CSV de aristas (columnas Source=fago, Target=bacteria).")
    parser.add_argument('--phage_fasta', type=str, required=True,
                         help="Ruta a all_phage_seq.fasta O a la carpeta fagos_fna con archivos individuales.")
    parser.add_argument('--output', type=str, required=True, help="Ruta de salida para nodes.csv.")
    args = parser.parse_args()

    df = pd.read_csv(args.edges)
    phage_names = build_phage_name_map(args.phage_fasta)

    fagos_ids = df['Source'].unique()
    fagos = pd.DataFrame({'Id': fagos_ids})
    fagos['Label'] = fagos['Id'].apply(lambda x: phage_names.get(x, x))
    fagos['Type'] = 'Fago'

    bacterias_ids = df['Target'].unique()
    bacterias = pd.DataFrame({'Id': bacterias_ids})
    bacterias['Label'] = bacterias['Id'].apply(clean_bacteria_name)
    bacterias['Type'] = 'Bacteria'

    nodos = pd.concat([fagos, bacterias], ignore_index=True)
    nodos.to_csv(args.output, index=False)

    sin_nombre = fagos[fagos['Label'] == fagos['Id']]
    print(f"Generado {args.output}: {len(fagos)} fagos, {len(bacterias)} bacterias, {len(nodos)} nodos totales.")
    if len(sin_nombre) > 0:
        ejemplos = list(sin_nombre['Id'])[:10]
        print(f"[!] {len(sin_nombre)} fagos no encontraron nombre descriptivo en el fasta (quedaron con su ID). Ejemplos: {ejemplos}")

if __name__ == "__main__":
    main()
