# 16S-FenrirPipeline
Taxa ans diversity pipeline for ONT 16S
Para el procesamiento de muestras de 16S con ONT, fue necesario hacer una revisión bibliográfica de otros estudios, con el fin de identificar aquellas herramientas que son útiles para diferentes estancias en la construcción del pipeline. Por lo que se identificaron los siguientes pasos:
- Basecalling
- Demultiplexing
- Denoising y Trimming
- Taxonomía
- Filogenía (opcional)
- Diversidad

La elaboración de este pipeline se hace desde un entorno Linux con Anaconda.

# Basecalling
Este es el primer paso del pipeline, pues aquí se busca convertir los archivos fast5 a fastq. Los archivos fast5 son el output de la secuenciación con MinION (Oxford Nanopore Technologies). En este caso, ONT proporciona Dorado como una herramienta para este paso. 
Para usar esta herramienta es necesario hacer la instalación de la misma para entorno Linux. Esta descarga se hizo por medio de línea de comando desde el Terminal, pues `Dorado` se encuentra como un paquete de `Anaconda` ([Dorado](https://anaconda.org/HCC/dorado)).

# Demultiplexing
Para realizar este paso, se encontraron que algunos estudios usaron herramientas especializadas como `Porechop`, sin embargo, se decidió usar `Dorado` porque esta herramienta es la recomendada por ONT.

# Denoising y Trimming
Aquí se busca depurar las muestras y reducir el ruido presente en las secuencias del gen 16S. Una de las herramientas encontradas en bibliografía para la limpieza de las muestras es `Chopper`. Para proceder con la instalación, existen dos métodos por línea de comando desde la terminal de Linux.

Desde [GitHub](https://github.com/epi2me-labs/pychopper):
```
conda install -c nanoporetech -c conda-forge -c bioconda "nanoporetech::pychopper"
```

Desde [Anaconda](https://anaconda.org/bioconda/pychopper):
```
conda install bioconda::pychopper 
```
Posteriormente, se debe revisar el comando base para usar `Chopper`, este puede verse con:
```
pychopper -h
```
Mostrando así la composición básica del comando.
```
usage: pychopper [-h] [-b primers] [-g phmm_file] [-c config_file] [-k {PCS109,PCS110,PCS111,PCS114,LSK114,PCB111,PCB114}] [-q cutoff] [-Q min_qual] [-z min_len] [-r report_pdf]
                 [-u unclass_output] [-l len_fail_output] [-w rescue_output] [-S stats_output] [-K qc_fail_output] [-Y autotune_nr] [-L autotune_samples] [-A scores_output]
                 [-m method] [-x rescue] [-p] [-t threads] [-B batch_size] [-D read stats] [-y] [-U]
                 input_fastx [output_fastx]
```
En este caso se debe crear una carpeta y su respectiva ruta, para que allí sean descargados los resultados después de la limpieza de la muestra. Adicionalmente, que se tenga una archivo `.fasta` con los primers y otro `.txt` con la orientación de los mismos. De esta manera, ya es posible hacer uso del comando base de `Pychopper`.
```
mkdir -p resultados

pychopper -m edlib \ #Seleccionar el método
  -b resources/primers_stesen.fasta \ #Identifica la ruta y el archivo donde está la secuencia de los primers
  -c resources/primer_stesen.txt \ #Identica la ruta y el archivo con la orientación de los mismos
  -Q 9 -z 1200 -t 4 \ #`-Q 9` QScore mínimo, `-z 1200` longitud mínima en pb, `-t` cantidad de núcleos a usar
  -r results/sample1_report.pdf \ #Reporte gráfico de los resultados en formato PDF
  -u results/sample1_unclassified.fastq \ #fastq de las lecturas que se excluyeron
  -w results/sample1_rescued.fastq \ #fastq de las lecturas recuperadas
  "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels/SRR26147154.fastq.gz" \ #Dirección de la muestra a tratar
  "results/sample1_oriented_trimmed.fastq" 
```
Sin embargo, `Pychopper` es una herramienta que filtra según la longitud, pero solo permite considerar largo mínimo. Por lo que, se usó de forma complementaria la herramienta `Filtlong`, pues esta permite establecer el rango ideal de 1300pb a 1700pb. Para instalar esta herramienta se usó el comando:
```
 conda install bioconda::filtlong
```
En este caso, los archivos fueron filtrados inicialmente con `Filtlong` y posteriormente con `Pychopper`. De esta manera, los archivos resultantes serán filtrados según las condiciones que se desee establecer. Adicionalmente, se realizó el código para no tener que limpiar cada archivo por separado, si no que, todos las secuencais son filtradas y Los archivos `fastq` resultantes son dispuestos en su carpeta correspondiente.
```
mkdir -p results && \ #Crear carpeta
filtlong --min_length 1300 --max_length 1700 --keep_percent 90 \ #Condiciones de filtrado para Filtlong
  "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels/SRR26147165.fastq.gz" \ #Dirección de los archivos originales
| pychopper -m edlib \ 
  -b "/home/fenrir/resources/primers_stesen.fasta" \ #Secuencia de los primers
  -c "/home/fenrir/resources/primers_stesen.txt" \ #Dirección de los primers
  -Q 9 -z 1300 -t 8 \ #Condiciones de filtrado para Pychopper
  -Y 0 -q 0.52 \ #`-Y` cuántos reads muestrea `q`, `-q` es para definir que tan exigente es en el alineamiento de los primers
  -r "results/SRR26147165_report.pdf" \ #Genera métricas y gráficos de la limpieza por Pychopper
  -S "results/SRR26147165_stats.tsv" \ #Estadísticas globales tabuladas
  -A "results/SRR26147165_scores.tsv" \ #Score por lectura
  -K "results/SRR26147165_qc_fail.fastq" \ #Lecturas que no pasaron el filtro de calidad
  -l "results/SRR26147165_len_fail.fastq" \ #Lecturas que quedan por debajo de la longitud mínima (debería ser 0)
  -u "results/SRR26147165_unclassified.fastq" \ #Lecturas no clasificadas
  -w "results/SRR26147165_rescued.fastq" \ #Lecturas rescatadas
  /dev/stdin "results/SRR26147165_oriented_trimmed.fastq" #Obtiene el archivo que saca Fitlong y lo reemplaza con las muestras filtradas y orientadas
```
```
mkdir -p results && \
filtlong --min_length 1300 --max_length 1700 \
  "$IN" \
| pychopper -m edlib \
    -b "$PRIMERS" -c "$PCONFIG" \
    -Q 9 -z 1300 -t 8 -Y 0 -q 0.52 \
    -r "results/${SAMPLE}_report.pdf" \
    -S "results/${SAMPLE}_stats.tsv" \
    -A "results/${SAMPLE}_scores.tsv" \
    -K "results/${SAMPLE}_qc_fail.fastq" \
    -l "results/${SAMPLE}_len_fail.fastq" \
    -u "results/${SAMPLE}_unclassified.fastq" \
    -w "results/${SAMPLE}_rescued.fastq" \
    /dev/stdin "results/${SAMPLE}_oriented_trimmed.fastq"
```
# Taxonomía
EMU
```
conda create -n emu -c conda-forge -c bioconda emu
```
Librerías de EMU
```
# Definir dónde guardar la base
export EMU_DATABASE_DIR="/home/fenrir/emu_db"
mkdir -p "$EMU_DATABASE_DIR"
cd "$EMU_DATABASE_DIR"

# Instalar cliente para descargar desde OSF
conda install -c conda-forge osfclient

# Descargar la base de datos
osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar

# Extraer contenido
tar -xvf emu.tar
```
Graficas
```
conda install -c conda-forge graphviz pygraphviz
```
exd
```
python - <<'PY'
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# 1) Cargar la tabla de EMU
df = pd.read_csv("results_emu/SRR26147165_rel-abundance.tsv", sep="\t")

# 2) Quedarnos SOLO con la clasificación (sin abundancias ni diversidad)
tax_cols = ["superkingdom","phylum","class","order","family","genus","species"]
df_tax = df[tax_cols].fillna("Unassigned")

# 3) Guardar una tabla con solo la clasificación
out_csv = "results_emu/SRR26147165_taxonomy_only.csv"
df_tax.to_csv(out_csv, index=False)
print(f"Tabla de clasificación guardada en: {out_csv}")

# 4) Construir un grafo jerárquico (superkingdom -> ... -> species)
G = nx.DiGraph()
for _, row in df_tax.iterrows():
    path = [row[c] for c in tax_cols if row[c] != "Unassigned"]
    for i in range(len(path)-1):
        G.add_edge(path[i], path[i+1])

# 5) Dibujar el grafo (solo clasificación, sin abundancias)
plt.figure(figsize=(14, 10))
pos = nx.spring_layout(G, k=0.8, seed=42)
nx.draw(G, pos, with_labels=True, node_size=1200, font_size=7, node_color="lightblue", edge_color="gray")
plt.title("EMU — Clasificación taxonómica (estructura jerárquica)")
plt.tight_layout()
out_png = "results_emu/SRR26147165_taxonomy_graph.png"
plt.savefig(out_png, dpi=150)
print(f"Grafo guardado en: {out_png}")
PY
```
Crear TSV de EMU
```
emu abundance \
  --db "/home/fenrir/emu_db" \
  --threads 8 \
  --output-dir results_emu \
  --output-basename SRR26147165 \
  results/SRR26147165/SRR26147165_oriented_trimmed.fastq
```
Para todas las carpetas
```
set -euo pipefail

DB="/home/fenrir/emu_db"        # Base de EMU
THREADS=8
OUTDIR="results_emu"
mkdir -p "$OUTDIR"

echo "=== Iniciando pipeline completo ==="

# --- EMU para todas las SRR filtradas ---
found=0
for fq in results/*/SRR*_oriented_trimmed*.fastq*; do
  [ -e "$fq" ] || continue
  found=1
  
  # Detectar el SRR
  sample="$(basename "$fq" | grep -o 'SRR[0-9]\{6,\}' || true)"
  if [ -z "$sample" ]; then
    sample="$(basename "$(dirname "$fq")")"
  fi

  # Crear subcarpeta para la muestra
  sample_dir="$OUTDIR/$sample"
  mkdir -p "$sample_dir"

  echo ">>> Procesando muestra: $sample"
  echo "    FASTQ: $fq"
  echo "    Carpeta de salida: $sample_dir"

  # Ejecutar EMU
  emu abundance \
    --db "$DB" \
    --threads "$THREADS" \
    --output-dir "$sample_dir" \
    --output-basename "$sample" \
    --keep-read-assignments --keep-counts \
    "$fq"
done

if [ "$found" -eq 0 ]; then
  echo "No se encontraron archivos FASTQ finales en results/*/SRR*_oriented_trimmed*.fastq*"
  exit 0
fi

# --- Postproceso: solo taxonomía (tabla + 2 grafos) ---
python - <<'PY'
import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
try:
    from networkx.drawing.nx_agraph import graphviz_layout
    HAS_DOT = True
except Exception:
    HAS_DOT = False

BASE_DIR = "results_emu"
tax_cols = ["superkingdom","phylum","class","order","family","genus","species"]

# Buscar TSV en todas las subcarpetas
for root, _, files in os.walk(BASE_DIR):
    for f in files:
        if not f.endswith(".tsv"):
            continue
        if "rel-abundance" not in f:
            continue
        
        tsv_path = os.path.join(root, f)
        sample = os.path.splitext(f)[0]
        print(f"\nPostproceso para {sample}")
        
        try:
            df = pd.read_csv(tsv_path, sep="\t")
        except Exception as e:
            print(f"  ⚠ Error leyendo {tsv_path}: {e}")
            continue

        # Verificar columnas taxonómicas
        col_low = {c.lower(): c for c in df.columns}
        present = {c: col_low.get(c) for c in tax_cols}
        if sum(v is not None for v in present.values()) < 3:
            print(f"  ⚠ {tsv_path} no parece tabla taxonómica.")
            continue

        # Tabla solo taxonomía
        cols_src = []
        for c in tax_cols:
            src = present[c]
            if src is None:
                df[c] = "Unassigned"
                cols_src.append(c)
            else:
                cols_src.append(src)

        df_tax = df[cols_src].copy()
        df_tax.columns = tax_cols
        df_tax = df_tax.fillna("Unassigned")

        out_csv = os.path.join(root, f"{sample}_taxonomy_only.csv")
        df_tax.to_csv(out_csv, index=False)
        print(f"  ✔ Tabla taxonómica: {out_csv}")

        # Grafo
        G = nx.DiGraph()
        for _, row in df_tax.iterrows():
            path = [row[c] for c in tax_cols if row[c] != "Unassigned"]
            for i in range(len(path)-1):
                G.add_edge(path[i], path[i+1])

        # Grafo spring
        plt.figure(figsize=(14, 10), constrained_layout=True)
        pos_spring = nx.spring_layout(G, k=0.8, seed=42)
        nx.draw(G, pos_spring, with_labels=True,
                node_size=1200, font_size=7,
                node_color="lightblue", edge_color="gray")
        plt.title(f"EMU — Jerarquía taxonómica (spring) — {sample}")
        out_png1 = os.path.join(root, f"{sample}_taxonomy_graph_spring.png")
        plt.savefig(out_png1, dpi=150)
        plt.close()
        print(f"  ✔ Grafo (spring): {out_png1}")

        # Grafo jerárquico (Graphviz)
        if HAS_DOT:
            try:
                plt.figure(figsize=(12, 14), constrained_layout=True)
                pos_dot = graphviz_layout(G, prog="dot")
                nx.draw(G, pos_dot, with_labels=True,
                        node_size=1200, font_size=7,
                        node_color="lightblue", edge_color="gray")
                plt.title(f"EMU — Jerarquía taxonómica (jerárquico) — {sample}")
                out_png2 = os.path.join(root, f"{sample}_taxonomy_graph_hier.png")
                plt.savefig(out_png2, dpi=150)
                plt.close()
                print(f"  ✔ Grafo (jerárquico): {out_png2}")
            except Exception as e:
                print(f"  ⚠ Error con Graphviz para {sample}: {e}")
        else:
            print("  ⚠ Graphviz no disponible, omitiendo grafo jerárquico.")
PY
```

# Métricas de diversidad
