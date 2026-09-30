import os
import subprocess
import glob

BASE_DIR = "/home/fenrir/Documentos/Tesis/Algoritmo/DeepPBI_Workspace"
HOST_IN = os.path.join(BASE_DIR, "host_raw_data")
PHAGE_IN = os.path.join(BASE_DIR, "phage_raw_data")
HOST_OUT = os.path.join(BASE_DIR, "host_annotation")
PHAGE_OUT = os.path.join(BASE_DIR, "phage_annotation")

PROKKA_BIN = subprocess.getoutput("which prokka")

def run_prokka(input_folder, output_base, kingdom):
    print(f"\n[*] Iniciando anotación de {kingdom}...")
    files = glob.glob(os.path.join(input_folder, "*.fna")) + \
            glob.glob(os.path.join(input_folder, "*.fasta"))
    
    if not files:
        print(f"[!] No se encontraron archivos en {input_folder}")
        return

    for f in files:
        sample_name = os.path.basename(f).split('.')[0]
        out_dir = os.path.join(output_base, f"prokka_{sample_name}")
        
        if os.path.exists(out_dir):
            continue
            
        print(f"[>] Anotando: {sample_name}")
        cmd = [
            PROKKA_BIN, f,
            "--outdir", out_dir,
            "--prefix", sample_name,
            "--kingdom", kingdom,
            "--cpus", "4",
            "--quiet",
            "--notbl2asn"  # Ahora sí funcionará con tu versión 1.13
        ]
        subprocess.run(cmd)

if __name__ == "__main__":
    os.makedirs(HOST_OUT, exist_ok=True)
    os.makedirs(PHAGE_OUT, exist_ok=True)
    run_prokka(HOST_IN, HOST_OUT, "Bacteria")
    run_prokka(PHAGE_IN, PHAGE_OUT, "Viruses")
    print("\n[FIN] Proceso terminado.")
