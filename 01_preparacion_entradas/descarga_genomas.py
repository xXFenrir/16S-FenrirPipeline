# Descarga de NCBI los genomas de las bacterias con abundancia promedio >= --min_abund en la tabla de EMU

import pandas as pd
import subprocess
import os
import glob
import shutil
import time
import argparse

MAX_ZIP_BYTES = 300 * 1024 * 1024  # 300MB, uno o pocos genomas no deberían pasar de esto

def descargar_ncbidatasets(query, filename_prefix, desc_log, out_bact, reference=False):
    # primero genomas completos, si no hay, chromosome/scaffold
    niveles = ["complete", "chromosome,scaffold"]
    for nivel in niveles:
        archivo_zip = f"temp_{filename_prefix}.zip"
        if os.path.exists(archivo_zip):
            os.remove(archivo_zip)
        cmd = [
            "datasets", "download", "genome", "taxon", str(query),
            "--assembly-level", nivel, "--filename", archivo_zip
        ]
        if reference:
            # sin --reference por género baja todos los ensamblajes del género (decenas de GB)
            cmd.append("--reference")

        # se vigila el tamaño del zip: algunas especies (p. ej. Klebsiella pneumoniae) tienen
        # miles de ensamblajes y la descarga puede ser de decenas de GB
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        aborted = False
        while proc.poll() is None:
            if os.path.exists(archivo_zip) and os.path.getsize(archivo_zip) > MAX_ZIP_BYTES:
                proc.kill()
                proc.wait()
                aborted = True
                print(f"   [!] Descarga abortada: supera {MAX_ZIP_BYTES // (1024*1024)}MB (probable descarga masiva no intencionada).")
                break
            time.sleep(1)
        returncode = -1 if aborted else proc.returncode

        if not aborted and returncode == 0 and os.path.exists(archivo_zip) and os.path.getsize(archivo_zip) > 1500:
            extract_path = f"temp_extract_{filename_prefix}"
            subprocess.run(["unzip", "-q", "-o", archivo_zip, "-d", extract_path])
            fna_files = glob.glob(f"{extract_path}/ncbi_dataset/data/**/*.fna", recursive=True)
            
            if fna_files:
                destino = os.path.join(out_bact, f"{filename_prefix}.fna")
                shutil.move(fna_files[0], destino)
                shutil.rmtree(extract_path)
                os.remove(archivo_zip)
                print(f"   [OK] {desc_log} ({nivel}) guardado en {destino}.")
                return True
                
            if os.path.exists(extract_path): 
                shutil.rmtree(extract_path)
                
        if os.path.exists(archivo_zip): 
            os.remove(archivo_zip)
            
    return False

def procesar_bacteria(feature_id, num_actual, total, out_bact):
    parts = feature_id.split('|')
    if len(parts) < 3: 
        return False
    
    taxid_especie = parts[1]
    nombre_completo = parts[2]
    genero = nombre_completo.split(' ')[0]
    nombre_archivo = nombre_completo.replace(" ", "_").replace(".", "")

    print(f"\n[{num_actual}/{total}] Intentando: {nombre_completo} (TaxID: {taxid_especie})")
    
    # 1. por TaxID de la especie
    if descargar_ncbidatasets(taxid_especie, f"taxid_{taxid_especie}_{nombre_archivo}", "Especie exacta", out_bact):
        return True
    
    # 2. por nombre
    print(f"   [!] No hallado por TaxID. Intentando por nombre científico...")
    if descargar_ncbidatasets(nombre_completo, f"name_{nombre_archivo}", "Búsqueda por nombre", out_bact):
        return True

    # 3. un representante del género
    print(f"   [!] No hallado por nombre. Buscando representante del género: {genero}...")
    if descargar_ncbidatasets(genero, f"GENERO_{genero}_{taxid_especie}", f"Representante de {genero}", out_bact, reference=True):
        return True

    print(f"   [x] Sin éxito para {nombre_completo} en ninguna categoría.")
    return False

def vincular_fagos(phage_source, out_phage):
    if not phage_source or not out_phage:
        return
        
    print(f"\n[*] Vinculando fagos desde {phage_source} a {out_phage}...")
    os.makedirs(out_phage, exist_ok=True)
    fagos = glob.glob(os.path.join(phage_source, "*.fasta")) + glob.glob(os.path.join(phage_source, "*.fna"))
    
    for f in fagos:
        destino = os.path.join(out_phage, os.path.basename(f))
        if not os.path.exists(destino):
            try: 
                os.symlink(f, destino)
            except OSError: 
                shutil.copy(f, destino)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filtra abundancias relativas y descarga genomas WGS desde NCBI.")
    
    parser.add_argument('-t', '--tsv', required=True, type=str, help='Ruta al archivo TSV de abundancias relativas de EMU.')
    parser.add_argument('-ob', '--out_bact', required=True, type=str, help='Ruta donde se guardarán los FASTA/FNA bacterianos.')
    
    parser.add_argument('-ps', '--phage_source', type=str, help='Ruta de origen de los fagos (opcional).')
    parser.add_argument('-op', '--out_phage', type=str, help='Ruta donde se vincularán los fagos (opcional).')
    parser.add_argument('-min', '--min_abund', type=float, default=0.001, help='Abundancia relativa promedio mínima (por defecto: 0.001).')
    
    args = parser.parse_args()

    os.makedirs(args.out_bact, exist_ok=True)
    
    print(f"[*] Leyendo tabla de abundancias: {args.tsv}")
    df = pd.read_csv(args.tsv, sep='\t', index_col=0)
    
    df['mean_abundance'] = df.mean(axis=1)
    
    taxones_filtrados = df[df['mean_abundance'] >= args.min_abund].sort_values(by='mean_abundance', ascending=False)
    total_taxones = len(taxones_filtrados)
    
    print(f"[*] Se encontraron {total_taxones} taxones con abundancia promedio >= {args.min_abund} ({(args.min_abund * 100):.3f}%)")

    exitos = 0
    for i, (feature_id, row) in enumerate(taxones_filtrados.iterrows(), 1):
        if procesar_bacteria(feature_id, i, total_taxones, args.out_bact):
            exitos += 1
        time.sleep(1) # límite de peticiones de NCBI

    if args.phage_source and args.out_phage:
        vincular_fagos(args.phage_source, args.out_phage)

    print(f"\n[FIN] Total obtenido: {exitos}/{total_taxones} genomas bacterianos en {args.out_bact}")
