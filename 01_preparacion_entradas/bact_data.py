# Versión anterior de descarga_genomas.py (la usada en el Objetivo 2): las N bacterias más abundantes, con rutas fijas

import pandas as pd
import subprocess
import os
import glob
import shutil
import time
import argparse

PATH_ABUNDANCIA = '/home/fenrir/Documentos/Tesis/results_dentrim_Q9/EMU_Q9/feature_table_relabund.tsv'
PHAGE_SOURCE = '/home/fenrir/Documentos/Tesis/Algoritmo/datos_phage/phage_raw_data'
WORK_DIR = '/home/fenrir/Documentos/Tesis/Algoritmo/DeepPBI-KG_scripts'

HOST_RAW = os.path.join(WORK_DIR, 'host_raw_data')
PHAGE_RAW = os.path.join(WORK_DIR, 'phage_raw_data')

def setup_workspace():
    print(f"[*] Configurando directorios en {WORK_DIR}...")
    dirs = [HOST_RAW, PHAGE_RAW, os.path.join(WORK_DIR, 'all_host_db')]
    for d in dirs:
        os.makedirs(d, exist_ok=True)

def descargar_ncbidatasets(query, filename_prefix, desc_log):
    # primero genomas completos, si no hay, chromosome/scaffold
    niveles = ["complete", "chromosome,scaffold"]
    for nivel in niveles:
        archivo_zip = f"temp_{filename_prefix}.zip"
        cmd = [
            "datasets", "download", "genome", "taxon", str(query),
            "--assembly-level", nivel, "--filename", archivo_zip
        ]
        result = subprocess.run(cmd, capture_output=True)
        
        if result.returncode == 0 and os.path.exists(archivo_zip) and os.path.getsize(archivo_zip) > 1500:
            extract_path = f"temp_extract_{filename_prefix}"
            subprocess.run(["unzip", "-q", "-o", archivo_zip, "-d", extract_path])
            fna_files = glob.glob(f"{extract_path}/ncbi_dataset/data/**/*.fna", recursive=True)
            if fna_files:
                destino = os.path.join(HOST_RAW, f"{filename_prefix}.fna")
                shutil.move(fna_files[0], destino)
                shutil.rmtree(extract_path); os.remove(archivo_zip)
                print(f"   [OK] {desc_log} ({nivel}) guardado.")
                return True
            if os.path.exists(extract_path): shutil.rmtree(extract_path)
        if os.path.exists(archivo_zip): os.remove(archivo_zip)
    return False

def procesar_bacteria(feature_id, num_actual, total):
    parts = feature_id.split('|')
    if len(parts) < 3: return False
    
    taxid_especie = parts[1]
    nombre_completo = parts[2]
    genero = nombre_completo.split(' ')[0]
    nombre_archivo = nombre_completo.replace(" ", "_").replace(".", "")

    print(f"\n[{num_actual}/{total}] Intentando: {nombre_completo} (TaxID: {taxid_especie})")
    
    # 1. por TaxID de la especie
    if descargar_ncbidatasets(taxid_especie, f"taxid_{taxid_especie}_{nombre_archivo}", "Especie exacta"):
        return True
    
    # 2. por nombre
    print(f"   [!] No hallado por TaxID. Intentando por nombre científico...")
    if descargar_ncbidatasets(nombre_completo, f"name_{nombre_archivo}", "Búsqueda por nombre"):
        return True

    # 3. cualquier genoma del género
    print(f"   [!] No hallado por nombre. Buscando cualquier representante del género: {genero}...")
    if descargar_ncbidatasets(genero, f"GENERO_{genero}_{taxid_especie}", f"Representante de {genero}"):
        return True

    print(f"   [!] Sin éxito para {nombre_completo} en ninguna categoría.")
    return False

def vincular_fagos():
    print(f"\n[*] Vinculando fagos...")
    fagos = glob.glob(os.path.join(PHAGE_SOURCE, "*.fasta")) + glob.glob(os.path.join(PHAGE_SOURCE, "*.fna"))
    for f in fagos:
        destino = os.path.join(PHAGE_RAW, os.path.basename(f))
        if not os.path.exists(destino):
            try: os.symlink(f, destino)
            except: shutil.copy(f, destino)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num', type=int, default=15)
    args = parser.parse_args()

    setup_workspace()
    
    df = pd.read_csv(PATH_ABUNDANCIA, sep='\t', index_col=0)
    df['mean_abundance'] = df.mean(axis=1)
    top_n = df.sort_values(by='mean_abundance', ascending=False).head(args.num)

    exitos = 0
    for i, (feature_id, row) in enumerate(top_n.iterrows(), 1):
        if procesar_bacteria(feature_id, i, args.num):
            exitos += 1
        time.sleep(1)

    vincular_fagos()
    print(f"\n[FIN] Total obtenido: {exitos}/{args.num} genomas en {HOST_RAW}")
