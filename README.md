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
Este es el primer paso del pipeline, pues aquí se busca convertir los archivos fast5 a `.fastq`. Los archivos `.fast5` son el output de la secuenciación con MinION (Oxford Nanopore Technologies). En este caso, ONT proporciona Dorado como una herramienta para este paso. 
Para usar esta herramienta es necesario hacer la instalación de la misma para entorno Linux. Esta descarga se hizo por medio de línea de comando desde el Terminal, pues `Dorado` se encuentra como un paquete de `Anaconda` ([Dorado](https://anaconda.org/HCC/dorado)). Se descarga el binario precompilado con:
```
curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.1.1-linux-x64.tar.gz" -o dorado-1.1.1-linux-x64.tar.gz
```

Luego, se debe extraer el archivo con:
```
tar -xzf dorado-1.1.1-linux-x64.tar.gz
```

Sin embargo, en las versiones actuales de `Dorado` no se utilizan como entradas los archivos en formato `.fast5` si no archivos `.POD5`, por lo que antes de aplicar el comando básido de la herrmaienta se deben convertir los archivos de ser necesario.
```
pod5 convert fast5 /ruta/fast5/*.fast5 --output pod5_out/

dorado basecaller hac /ruta/pod5/ \
  --emit-fastq \
  --no-trim \
  -x auto > basecalls.fastq
```

- `hac` es para usar el modelo de alta precicsión.
- `--emit-fastq` hace que la salida sea en formato `.fastq`
- `--no-trim` es para evitar que recorte barcodes y adaptadores.
- `-x auto` decide automáticamente si usar GPU o CPU.

# Demultiplexing
Para realizar este paso, se encontraron que algunos estudios usaron herramientas especializadas como `Porechop`, sin embargo, se decidió usar `Dorado` porque esta herramienta es la recomendada por ONT.

Generalmente, este paso se hace al mismo tiempo que el basecalling, por lo que se modifica un poco el comando de basecalling anteriormente puesto. Para usar esta herramienta se usó el siguiente código:
```
dorado basecaller hac /ruta/pod5/ \
  --kit-name SQK-16S114-24 \
  --trim all \
  --primer-sequences primers_27F_1492R.fasta \
  -x auto > calls.bam

dorado demux --no-classify --emit-fastq -o demux_fastq calls.bam
```

- `--kit name` nombre del kit de barcodes usados.
- `-- trim all` quitar barcodes y primers. En este caso no afecta, pues dorado garantiza que el recorte no interfiere con la demultiplexación.
- `--primer-sequences` archivo `.fasta` von las secuencias de los primers.
- `-x auto` automáticamente elige si usar la GPU o CPU.
- `--no-clasiffy` no vuelve a clasificar, simplemente lee la asignación de barcode que ya quedó guardada en `calls.bam` durante el basecalling, y divide las lecturas en archivos separados por barcode.
- `--emit-fastq` la salida pasa a ser en formato `.fastq`.
- `-o` carpeta donde estaran las salidas.

# Reporte de estadísticas básicas

Antes de proceder con este paso, es necesario saber si las muestras necesitan o no ser limpiados, pues en el caso de los repositorios algunos ya vienen filtrados. Por lo que, se deben tabular la estadísticas básicas de los archivos `.fastq`, de tal manera que se pueda observar si vale la pena filtrar o estos ya están limpios. Por lo que, se usó el código presente en el archivo de [Reporte de estadísticas](https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/Reporte%20de%20estad%C3%ADsticas)

Los paquetes que se importanron fueron:
- `argparse` evita que se tenga que editar el código si se cambia el formato de la entradaa tanto con una carpeta como con un archivo.
- `gzip` lee `.fastq.gz` sin descomprimir.
- `os` para saber tamaño del archivo.
- `pathlib.Path` reconoce espacios, ~, etc.
- `re` permite reconocer solo las bases nitrogenadas en archivos FASTA.
- `sys` permite imprimir avisos.
- `time` permite conocer el tiempo que demoró en ejecutarse.
- `shutil` permite ver si es posible convertir un archivo de TSV a XLSX
- `subprocess` permite separa cada columna en el archivo XLSX.
- `pandas` crea un XLSX con las columnas establecidas.

Con el parámetro `input` puede ser un archivo o una carpeta, lo que permite procesar una sola muestra o una carpeta con varios de estas. Además, con `--recursive` se busca también en subcarpetas, ignorando así las jerarquías. Luego, con `--output` se define la ruta y el nombre del archivo, y con `--tsv` permite generar un archivo separado por tabulaciones. Usando `--phred` se puede conocer la calidad de las lecturas. Finalmente, con `--max-reads` se procesa solo las primeras N lecturas de cada archivo y se detiene.

Para abrir `.fastq.gz` se usó `gzip.open` pues evita descomprimir y mantienes un flujo de lectura constante. Además, con `encoding="ascii"` y `errors="ignore"` se evita que cuando aparezca un carácter extraño, sean ignorados y sacar métricas del resto sin detenerse.

`iter_fastq_reads` lee cada cuatro líneas del archivo `.fastq`, pues corresponden al encabezado, la secuencia y la calidad. Como salida entrega únicamente la secuencia y su calidad, deteniéndose si encontraba el fin del archivo.

Adicionalmente, se incluyeron funciones para encontrar presencia de primers y quimeras. Por lo que, para saber si todavía quedan primers en los extremos de las lecturas, se debe pasar un archivo `.fasta` con las secuencias de tus primers usando `--primers`, y el programa examina solo los extremos de cada lectura. Para tolerar errores de secuenciación, permite con `--primer-max-mismatches` desajustes y además prueba tanto la secuencia del primer como su reverse. Mientras que, la estimación de quimeras se activa con `--chimera-denovo` y requiere tener `vsearch` instalado. En este caso, se toma una muestra de lecturas por archivo, se convierten a `.fasta` temporal y luego correr `vsearch --uchime3_denovo` para detectar quimeras sin referencia.

Para la parte estadística, se usó `n50_from_lengths` para optener el valor N50 como indicador de la longitud típica por bases de las lecturas. Otras de las métricas son el número de lecturas, bases totales, estadísticos de longitud (media, mediana, mínimo, máximo y N50), composición (GC% y N%) y calidad (Q media por base y % de bases ≥Q20/≥Q30). En este caso, la mediana se incluye para hacer más sensible el código a valores anormales, GC% y N% ayudan a detectar contaminación y revela ambigüedades o errores de lectura, la calidad media por base, y los umbrales Q20/Q30 como referencias comparativas. Finalmente, el script asume Phred+33, aunque, se expone `--phred` en caso de datos atípicos.

La búsqueda de archivos con `find_fastqs` acepta `.fastq` y `.fq`. Adicionalmente, en `main`, si se termina en `.tsv` se fuerzan tabs aunque olvide. Por último, se imprime `[INFO] con el número de archivos detectados y `[OK]` con la ruta final.

La escritura del reporte usa `csv.DictWriter` para mantener un orden estable de columnas y evitar errores de formato. Se abre el archivo con `newline=""` para prevenir líneas en blanco extra en algunos sistemas. Si al procesar un archivo ocurre una excepción, se imprime un `[ERROR]` en consola y se continúa; el diccionario que devuelve process_fastq puede incluir un campo error, pero la fila se escribe con los campos conocidos. Esta tolerancia controlada es intencional en trabajos por lote: te permite terminar la corrida y luego revisar con calma los casos problemáticos.

Se usó `csv.DictWriter` para garantizar un orden fijo de columnas y un formato consistente sin pelear con separadores. Luego, se abre el archivo con `newline=""` para evitar líneas en blanco extra. Si al procesar un `.fastq` ocurre una excepción, se imprime un `[ERROR] en consola y el script continúa con los demás archivos, donde la fila problemática se escribe solo con los campos conocidos, pero con las métricas vacías.

Una vez creado el script con el código, se ejecúta:
```
python3 "/home/fenrir/scriptsbioinf/fastq_estads.py" \
  "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels" \
  -r \
  --primers "/home/fenrir/resources/primers_stesen.fasta" \
  -o "/home/fenrir/scriptsbioinf/estadisticas_pretty.txt" \
  --xlsx-out "/home/fenrir/scriptsbioinf/estadisticas.xlsx"
```
(SE DEBE AGREGAR LOS CAMBIOS QUE SE LE HICIERON AL CÓDIGO, QUE ES EL # DE BASES Y LA SALIDA EN FORMATO .XSLX)

# Denoising y Trimming

Es uno de los pasos más importantes en un análisis bioinformático de datos de secuenciación 16S rRNA es la depuración de lecturas crudas. Pues el objetivo aquí es mejorar la calidad de los datos y asegurar que únicamente las lecturas confiables y relevantes pasen a la etapa de taxonomía. En este caso se quiere:

- Filtrar lecturas por longitud y calidad, eliminando aquellas demasiado cortas, largas o de baja calidad que puedan corresponder a artefactos de secuenciación o fragmentos incompletos.
- Orientar y recortar las lecturas basándonos en la ubicación de los primers, asegurando que todas las secuencias tengan la misma dirección y contenido correcto.

Una de las herramientas encontradas en bibliografía para la limpieza de las muestras es `Pychopper`. Para proceder con la instalación, existen dos métodos por línea de comando desde la terminal de Linux.

Desde [GitHub](https://github.com/epi2me-labs/pychopper):
```
conda install -c nanoporetech -c conda-forge -c bioconda "nanoporetech::pychopper"
```

Desde [Anaconda](https://anaconda.org/bioconda/pychopper):
```
conda install bioconda::pychopper 
```

Pychopper identifica los primers en cada lectura, determina su orientación y recorta el amplicón para dejar únicamente la región de interés. Además, genera reportes que muestran cuántas lecturas fueron clasificadas, rechazadas o rescatadas. Adicionalmente, es una herramienta que permite filtrar aquellas lecturas que poseen un QScore por debajo de 9. 

Sin embargo, es necesario crear un archivo `.fasta` con la secuencia de los primers y otro `.txt` con la orientación de los mismos. De esta manera, ya es posible hacer uso del comando base de `Pychopper`.
```
mkdir -p resultados

pychopper -m edlib \ 
  -b "/home/fenrir/resources/primers_stesen.fasta" \ 
  -c "/home/fenrir/resources/primers_stesen.txt" \ 
  -Q 9 -z 1300 -t 8 \
  -Y 0 -q 0.52 \ 
  "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels/SRR26147154.fastq.gz" \ 
  "results/sample1_oriented_trimmed.fastq" 
```

- `mkdir -p` crear una carpeta.
- `-m edlib` habilita un alineador para encontrar los primers.
- `-b` ruta del archivo con la secuencia de los primers.
- `-c` ruta del archivo con la orientación de los primers.
- `-Q` filtro por QScore.
- `-z`filtro por longitud mínima.
- `-t` núcleos de CPU a usar.
- `-Y` desactiva la función de muestreo automático para hacer el proceso reproducible.
- `-q` ajusta la exigencia del alineamiento a los primers.
- `-r` reporte de los resultados en PDF.
- `-S` estadísticas globales tabuladas.
- `-A` score por lectura.
- `-K` lecturas que no pasaron el filtro de calidad.
- `-l` lecturas por debajo de la longitud mínima.
- `-u` lecturas excluidas.
- `-w` lecturas guardadas.

Las dos últimas líneas de código permiten buscar el archivo `.fastq` a limpiar. Luego, se define la ruta en que se va a guardar el archivo `.fastq` limpiado.

Sin embargo, `Pychopper` es una herramienta que filtra por longitud mínima y no por rango. Por lo que, se usó de forma complementaria la herramienta `Filtlong`, pues esta permite establecer el rango ideal de 1300pb a 1700pb. 

Para instalar esta herramienta se usó el comando:
```
 conda install bioconda::filtlong
```

En este caso, los archivos fueron filtrados inicialmente con `Filtlong` y posteriormente con `Pychopper`. De esta manera, los archivos resultantes serán filtrados según las condiciones que se desee establecer. Adicionalmente, se realizó el código para no tener que limpiar cada archivo por separado, si no que, todos las secuencais son filtradas y los archivos `.fastq` resultantes son dispuestos en su carpeta correspondiente.
```
mkdir -p results_2 && \ #Crear carpeta
filtlong --min_length 1300 --max_length 1700 \ 
  "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels/SRR26147165.fastq.gz" \ 
| pychopper -m edlib \ 
  -b "/home/fenrir/resources/primers_stesen.fasta" \ 
  -c "/home/fenrir/resources/primers_stesen.txt" \ 
  -Q 9 -z 1300 -t 8 \ 
  -Y 0 -q 0.52 \
  /dev/stdin "results/SRR26147165_oriented_trimmed.fastq" 
```

`Filtlong` tiene un comando más sencillo de usar, pues es una herramienta especializada en filtrar por rangom de longitud. Por lo que, `--min_length` define el rango mínimo y `--max_length` define el rango máximo. Adicionalmente, se debe inidcar la dirección del archivo a filtrar. Posteriormente, se usa `|` para que los archivos que salen de `Filtlong` entren directamente a `Pychopper`, y se usa `/dev/stdin` para que `Pychopper` lea como entrada lo que entra por `|`. Finalmente, el archivo entregado se espera que sea un archivo `.fastq` filtrado con los parámetros definidos en ambas herramientas.

Ahora, para que este código funcione para todos los archivos del estudio se usó:
```
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob  

BASE="/home/fenrir/Documentos/Muestras 16S/sterile_sentinels"
PRIMERS="/home/fenrir/resources/primers_stesen.fasta"
PCONFIG="/home/fenrir/resources/primers_stesen.txt"
OUTDIR="results3int"

mkdir -p "$OUTDIR"

inputs=( "$BASE"/*.fastq.gz "$BASE"/*.fq.gz "$BASE"/*.fastq "$BASE"/*.fq )
if ((${#inputs[@]}==0)); then
  echo "No se encontraron FASTQ en: $BASE" >&2
  exit 1
fi

for fq in "${inputs[@]}"; do
  SAMPLE=$(basename "$fq")
  SAMPLE=${SAMPLE%.gz}
  SAMPLE=${SAMPLE%.fastq}
  SAMPLE=${SAMPLE%.fq}

  sdir="$OUTDIR/$SAMPLE"
  mkdir -p "$sdir"
  echo "Procesando $SAMPLE …"

  filtlong --min_length 1300 --max_length 1700 "$fq" \
  | pychopper -m edlib \
      -b "$PRIMERS" -c "$PCONFIG" \
      -Q 9 -z 1300 -t 8 -Y 0 -q 0.52 \
      -r "$sdir/${SAMPLE}_report.pdf" \
      -S "$sdir/${SAMPLE}_stats.tsv" \
      -A "$sdir/${SAMPLE}_scores.tsv" \
      -K "$sdir/${SAMPLE}_qc_fail.fastq" \
      -l "$sdir/${SAMPLE}_len_fail.fastq" \
      -u "$sdir/${SAMPLE}_unclassified.fastq" \
      -w "$sdir/${SAMPLE}_rescued.fastq" \
      - "$sdir/${SAMPLE}_oriented_trimmed.fastq"
done
```
# CAAMBIOS (EL SIGUIENTE ES EL NUEVO)
```
cat > dentrim.py <<'PY'
#!/usr/bin/env python3
"""
dentrim.py — Filtra y orienta lecturas 16S (filtlong + pychopper) y produce SOLO el FASTQ final.

Provee:
  - run_dentrim(input_dir, outdir, minlen, maxlen, qscore, threads, primers, pconfig, verbose=True)
  - CLI si se ejecuta como script.

Requisitos en PATH: filtlong, pychopper, bash
"""

import argparse
import glob
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List


def _must_exist(path: str, kind: str) -> None:
    if kind == "dir" and not os.path.isdir(path):
        sys.exit(f"ERROR: No existe directorio: {path}")
    if kind == "file" and not os.path.isfile(path):
        sys.exit(f"ERROR: No existe archivo: {path}")


def _which_or_die(cmd: str) -> None:
    if shutil.which(cmd) is None:
        sys.exit(f"ERROR: No se encontró '{cmd}' en PATH")


def _sample_name(f: str) -> str:
    name = Path(f).name
    for suf in (".gz", ".fastq", ".fq"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name


def run_dentrim(
    input_dir: str,
    outdir: str,
    minlen: int,
    maxlen: int,
    qscore: int,
    threads: int,
    primers: str,
    pconfig: str,
    verbose: bool = True,
) -> Dict[str, List[str]]:
    """
    Ejecuta el pipeline sobre todos los FASTQ/FQ del directorio.

    Devuelve:
      {
        "processed": [lista de muestras procesadas],
        "skipped":   [lista de muestras saltadas],
        "outdir":    ruta de salida (str)
      }
    """
    # Validaciones de entrada
    _must_exist(input_dir, "dir")
    _must_exist(primers, "file")
    _must_exist(pconfig, "file")
    if minlen > maxlen:
        sys.exit(f"ERROR: minlen ({minlen}) > maxlen ({maxlen})")

    # Dependencias
    _which_or_die("filtlong")
    _which_or_die("pychopper")
    _which_or_die("bash")

    os.makedirs(outdir, exist_ok=True)

    # Recolectar entradas
    patterns = ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]
    inputs: List[str] = []
    for pat in patterns:
        inputs.extend(glob.glob(os.path.join(input_dir, pat)))
    if not inputs:
        sys.exit(f"No se encontraron FASTQ en: {input_dir}")

    if verbose:
        print("== Parámetros ==")
        print(f"Input:   {input_dir}")
        print(f"Output:  {outdir}")
        print(f"MinLen:  {minlen}")
        print(f"MaxLen:  {maxlen}")
        print(f"QScore:  {qscore}")
        print(f"Threads: {threads}")
        print(f"Primers: {primers}")
        print(f"PConfig: {pconfig}\n")

    processed: List[str] = []
    skipped: List[str] = []

    for fq in inputs:
        sample = _sample_name(fq)
        sdir = os.path.join(outdir, sample)
        os.makedirs(sdir, exist_ok=True)
        out_fastq = os.path.join(sdir, f"{sample}_oriented_trimmed.fastq")
        if verbose:
            print(f"Procesando {sample} …")

        # Quoted paths (maneja espacios)
        fq_q = shlex.quote(fq)
        primers_q = shlex.quote(primers)
        pconfig_q = shlex.quote(pconfig)
        out_q = shlex.quote(out_fastq)

        # filtlong lee el archivo directamente; pychopper desde stdin.
        cmd = (
            "set -o pipefail; "
            f"filtlong --min_length {minlen} --max_length {maxlen} {fq_q}"
            f" | pychopper -m edlib -b {primers_q} -c {pconfig_q}"
            f"     -Q {qscore} -z {minlen} -t {threads} -Y 0 -q 0.52"
            f"     - {out_q}"
        )
        proc = subprocess.run(cmd, shell=True, executable="/bin/bash")

        if proc.returncode != 0:
            # Pudo quedar en 0 lecturas o error interno en pychopper; limpiamos y marcamos como saltado
            try:
                os.remove(out_fastq)
            except FileNotFoundError:
                pass
            skipped.append(sample)
            if verbose:
                print(f"⚠️  {sample}: sin lecturas tras filtrado u otro error. Saltando.", file=sys.stderr)
            continue

        processed.append(sample)

    if verbose:
        print(f"\n✅ Listo. Resultados en: {outdir}")
        if skipped:
            print("⚠️  Muestras saltadas:", ", ".join(skipped))

    return {"processed": processed, "skipped": skipped, "outdir": outdir}


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Filtra y orienta lecturas 16S (filtlong+pychopper) y deja SOLO el FASTQ final."
    )
    p.add_argument("-i", "--input", required=True, help="Carpeta con FASTQ/FQ (.gz opcional)")
    p.add_argument("-o", "--outdir", required=True, help="Carpeta de salida")
    p.add_argument("--minlen", required=True, type=int, help="Longitud mínima (filtlong, pychopper -z)")
    p.add_argument("--maxlen", required=True, type=int, help="Longitud máxima (filtlong)")
    p.add_argument("-Q", "--qscore", required=True, type=int, help="QScore mínimo (pychopper -Q)")
    p.add_argument("-t", "--threads", required=True, type=int, help="Núcleos/hilos (pychopper -t)")
    p.add_argument("--primers", required=True, help="FASTA de primers")
    p.add_argument("--pconfig", required=True, help="Config de primers (TXT de pychopper)")
    return p


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    run_dentrim(
        input_dir=args.input,
        outdir=args.outdir,
        minlen=args.minlen,
        maxlen=args.maxlen,
        qscore=args.qscore,
        threads=args.threads,
        primers=args.primers,
        pconfig=args.pconfig,
        verbose=True,
    )


if __name__ == "__main__":
    main()
PY

chmod +x dentrim.py
```
# ASÍ SE LLAMA
```
/home/fenrir/scriptsbioinf/dentrim.py \
  -i "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels" \
  -o "dentrim_stesen" \
  --minlen 1000 --maxlen 1700 \
  -Q 12 -t 8 \
  --primers "/home/fenrir/resources/primers_stesen.fasta" \
  --pconfig "/home/fenrir/resources/primers_stesen.txt"
```
# MENCIONAR CAMBIOS
- `set -euo pipefail` evita que se ejecute pacialmente a causa de fallos.
- `shopt -s nullglob` si no hay coincidencias es igual a 0
- `BASE` carpeta donde esten las lecturas `.fastq`.
- `PRIMERS` archivo `.fasta` con las secuencias de primers.
- `PCONFIG` archivo con la orientación de primers.
- `OUTDIR` carpeta para todos los resultados de cada muestra.
- `inputs` un array para que sea capaz de reconocer cualquier forma en la que se pueda encontrar el archivo `.fastq`.

Con `if ((${#inputs[@]}==0))` se busca que cuando el array no tenga un sufijo de `.fastq` se detenga. Luego, con `for fq in "${inputs[@]}"; do` hace la iteración por cada archivo `.fastq` reconocido. Finalmente, con `SAMPLE=` se busca que poco a poco se quiten los sufijos del archivo hasta quedar únicamente con el directorio y este es el que se usa para nombrar la carpeta en la que se agruparan los resultados.













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

Conflictos con paquetes de R y Bioconductor.

```
cat > metricasd.py <<'PY'
#!/usr/bin/env python3
import os, sys, csv, math, argparse, gzip
from collections import defaultdict, OrderedDict

TAX_COLS = ["superkingdom","phylum","class","order","family","genus","species"]

# ---------- Utilidades de E/S ----------
def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def read_tsv(path):
    with open_any(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = [{(k.strip().lower() if k else k): (v.strip() if isinstance(v,str) else v)
                 for k,v in r.items()} for r in reader]
    return rows

def write_csv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header); w.writerows(rows)

# ---------- Lectura EMU ----------
def feature_id_from_row(r):
    sp = r.get("species","")
    if sp and sp.lower() != "unassigned":
        return sp
    parts = [r.get(c,"") for c in TAX_COLS]
    parts = [p for p in parts if p and p.lower()!="unassigned"]
    return ";".join(parts) if parts else "Unassigned"

def load_sample(sample_dir):
    rel = None; cnt = None
    for fn in os.listdir(sample_dir):
        if fn.endswith("rel-abundance.tsv"): rel = os.path.join(sample_dir, fn)
        elif fn.endswith("counts.tsv"):      cnt = os.path.join(sample_dir, fn)
    if not rel: return None

    rel_rows = read_tsv(rel)
    counts_by_key = {}
    if cnt:
        cnt_rows = read_tsv(cnt)
        def key_from_row(r): return "|".join([r.get(c,"") for c in TAX_COLS])
        counts_by_key = { key_from_row(r): int((r.get("count","0") or "0")) for r in cnt_rows }

    sample_id = os.path.basename(sample_dir)
    rel_map, count_map = {}, {}

    for r in rel_rows:
        fid = feature_id_from_row(r)
        try: abu = float(r.get("abundance","0") or 0.0)
        except: abu = 0.0
        rel_map[fid] = rel_map.get(fid, 0.0) + abu
        if counts_by_key:
            key = "|".join([r.get(c,"") for c in TAX_COLS])
            cval = counts_by_key.get(key, 0)
            count_map[fid] = count_map.get(fid, 0) + cval

    return sample_id, rel_map, (count_map if count_map else None)

def union_features(samples):
    feats = set()
    for _, rel_map, _ in samples: feats.update(rel_map.keys())
    for _, _, cnt_map in samples:
        if cnt_map: feats.update(cnt_map.keys())
    return sorted(feats)

# ---------- Métricas ----------
def bray_curtis(a, b):
    num = 0.0; den = 0.0
    for x,y in zip(a,b): num += abs(x-y); den += (x+y)
    return (num/den) if den>0 else 0.0

def jaccard_binary(a, b, eps=0.0):
    A = sum(1 for x in a if x>eps); B = sum(1 for y in b if y>eps)
    I = sum(1 for x,y in zip(a,b) if x>eps and y>eps)
    U = A + B - I
    return (1 - I / U) if U>0 else 0.0

def hellinger_euclidean(a, b):
    s = 0.0
    for x,y in zip(a,b):
        dx = math.sqrt(x) - math.sqrt(y); s += dx*dx
    return math.sqrt(s)

def shannon(p):
    s = 0.0
    for x in p:
        if x>0: s -= x*math.log(x)
    return s

def simpson_1D(p): return 1.0 - sum(x*x for x in p)
def pielou_evenness(H, S): return (H / math.log(S)) if S>1 and H>0 else 0.0
def chao1(counts):
    f1 = sum(1 for c in counts if c==1); f2 = sum(1 for c in counts if c==2)
    S_obs = sum(1 for c in counts if c>0)
    if S_obs==0: return 0.0
    return S_obs + (f1*f1)/(2.0*f2) if f2>0 else S_obs + (f1*(f1-1))/2.0

# ---------- Visual (SVG) ----------
def color_palette(n):
    cols=[]; 
    for i in range(n): cols.append(hsl_to_rgb((i*1.0/n), 0.6, 0.55))
    return cols
def hsl_to_rgb(h,s,l):
    def hue2rgb(p,q,t):
        if t<0: t+=1
        if t>1: t-=1
        if t<1/6: return p+(q-p)*6*t
        if t<1/2: return q
        if t<2/3: return p+(q-p)*(2/3 - t)*6
        return p
    q = l*(1+s) if l<0.5 else l+s - l*s; p = 2*l - q
    r = hue2rgb(p,q,h+1/3); g = hue2rgb(p,q,h); b = hue2rgb(p,q,h-1/3)
    return "#%02x%02x%02x" % (int(r*255), int(g*255), int(b*255))

def svg_bar_chart(title, labels, values, width=900, height=300, margin=50):
    maxv = max(values) if values else 1.0
    bar_w = (width-2*margin) / max(1,len(values))
    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    svg.append(f'<text x="{margin}" y="25" font-size="16" font-family="sans-serif">{title}</text>')
    svg.append(f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="#333"/>')
    for i,(lab,val) in enumerate(zip(labels,values)):
        h = 0 if maxv==0 else (val/maxv)*(height-2*margin)
        x = margin + i*bar_w + 4; y = height - margin - h
        svg.append(f'<rect x="{x}" y="{y}" width="{bar_w-8}" height="{h}" fill="#6a93d6"/>')
        svg.append(f'<text x="{x+bar_w/2-8}" y="{height-margin+14}" font-size="10" font-family="sans-serif" transform="rotate(45 {x+bar_w/2-8},{height-margin+14})">{lab}</text>')
    for k in range(5):
        y = height - margin - (k/4)*(height-2*margin); val = (k/4)*maxv
        svg.append(f'<line x1="{margin-5}" y1="{y}" x2="{width-margin}" y2="{y}" stroke="#eee"/>')
        svg.append(f'<text x="10" y="{y+4}" font-size="10" font-family="sans-serif">{val:.2f}</text>')
    svg.append('</svg>'); return "\n".join(svg)

def svg_heatmap(title, samples, M, width=900, height=450, margin=90):
    n = len(samples); cell = (width-2*margin)/max(1,n)
    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    svg.append(f'<text x="{margin}" y="25" font-size="16" font-family="sans-serif">{title}</text>')
    for i in range(n):
        for j in range(n):
            v = M[i][j]; c = int(255 - min(max(v,0.0),1.0)*180); col = f"rgb({c},{c+20},{255})"
            x = margin + j*cell; y = margin + i*cell
            svg.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{col}" stroke="white" stroke-width="0.5"><title>{samples[i]} vs {samples[j]}: {v:.3f}</title></rect>')
    for i,s in enumerate(samples):
        x = margin + i*cell + cell/2
        svg.append(f'<text x="{x}" y="{margin-10}" font-size="10" font-family="sans-serif" text-anchor="end" transform="rotate(-45 {x},{margin-10})">{s}</text>')
        svg.append(f'<text x="{margin-10}" y="{margin + i*cell + cell/2}" font-size="10" font-family="sans-serif" text-anchor="end">{s}</text>')
    svg.append('</svg>'); return "\n".join(svg)

def svg_stacked_bars(title, samples, series_dict, width=900, height=350, margin=60):
    order = list(series_dict.keys()); palette = color_palette(len(order))
    bar_w = (width-2*margin) / max(1,len(samples))
    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    svg.append(f'<text x="{margin}" y="25" font-size="16" font-family="sans-serif">{title}</text>')
    svg.append(f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="#333"/>')
    for i, sid in enumerate(samples):
        acc = 0.0; x = margin + i*bar_w + 4
        for k, name in enumerate(order):
            val = series_dict[name][i]; h = val*(height-2*margin); y = height - margin - acc - h
            svg.append(f'<rect x="{x}" y="{y}" width="{bar_w-8}" height="{h}" fill="{palette[k]}"><title>{name}: {val:.3f}</title></rect>')
            acc += h
        svg.append(f'<text x="{x+bar_w/2-8}" y="{height-margin+14}" font-size="10" font-family="sans-serif" transform="rotate(45 {x+bar_w/2-8},{height-margin+14})">{sid}</text>')
    lx, ly = width - margin - 220, 40
    for k, name in enumerate(order[:18]):
        svg.append(f'<rect x="{lx}" y="{ly + 18*k}" width="12" height="12" fill="{palette[k]}"/>')
        svg.append(f'<text x="{lx+18}" y="{ly + 18*k + 10}" font-size="11" font-family="sans-serif">{name}</text>')
    svg.append('</svg>'); return "\n".join(svg)

# ---------- Composición por rank (desde EMU) ----------
def value_for_rank(r, rank):
    rank = rank.lower()
    if rank in ("genus","family","species"):
        val = r.get(rank,"") or r.get(rank.capitalize(),"")
        if val: return val
    # heurísticas de respaldo
    if rank == "genus":
        sp = r.get("species","")
        if sp: return sp.strip().split()[0]
        parts = [r.get(c,"") for c in TAX_COLS]
        parts = [p for p in parts if p]
        if len(parts)>=6 and parts[5]: return parts[5]
        return "Unassigned"
    if rank == "family":
        parts = [r.get(c,"") for c in TAX_COLS]
        parts = [p for p in parts if p]
        if len(parts)>=5 and parts[4]: return parts[4]
        return "Unassigned"
    if rank == "species":
        sp = r.get("species","")
        if sp: return sp
        # si no hay species, intenta genus + " sp."
        g = value_for_rank(r, "genus")
        return f"{g} sp." if g else "Unassigned"
    return "Unassigned"

def collect_rank_matrix(emu_dir, rank):
    sample_dirs = sorted([os.path.join(emu_dir,d) for d in os.listdir(emu_dir)
                          if os.path.isdir(os.path.join(emu_dir, d))])
    samples = []; per_sample=[]
    for d in sample_dirs:
        rel = None
        for fn in os.listdir(d):
            if fn.endswith("rel-abundance.tsv"): rel = os.path.join(d, fn); break
        if not rel: continue
        rows = read_tsv(rel); r_map = defaultdict(float)
        for r in rows:
            try: abu = float(r.get("abundance","0") or 0.0)
            except: abu = 0.0
            key = value_for_rank(r, rank)
            r_map[key] += abu
        s = sum(r_map.values())
        if s>0:
            for k in list(r_map.keys()): r_map[k] = r_map[k]/s
        samples.append(os.path.basename(d)); per_sample.append(r_map)
    all_keys = sorted(set().union(*[set(m.keys()) for m in per_sample])) if per_sample else []
    rank_abund = {k: [per_sample[j].get(k,0.0) for j in range(len(samples))] for k in all_keys}
    return samples, rank_abund

# ---------- Pipeline principal ----------
def main():
    ap = argparse.ArgumentParser(description="Abundancia, alfa, beta + reporte HTML (sin R)")
    ap.add_argument("--emu-dir", default="results_emu", help="Carpeta con subcarpetas por muestra (EMU)")
    ap.add_argument("--out", default="mp_outputs", help="Carpeta de salida")
    ap.add_argument("--top", type=int, default=10, help="Top N categorías para barras apiladas")
    ap.add_argument("--rank", choices=["species","genus","family"], default="species",
                    help="Nivel taxonómico para composición apilada (default: species)")
    args = ap.parse_args()

    if not os.path.isdir(args.emu_dir):
        sys.stderr.write("No existe la carpeta %s\n" % args.emu_dir); sys.exit(1)

    # --- Cargar muestras de EMU ---
    sample_dirs = sorted([os.path.join(args.emu_dir, d) for d in os.listdir(args.emu_dir)
                          if os.path.isdir(os.path.join(args.emu_dir, d))])
    samples = [load_sample(d) for d in sample_dirs]
    samples = [s for s in samples if s]
    if not samples:
        sys.stderr.write("No se encontraron rel-abundance.tsv en %s\n" % args.emu_dir); sys.exit(1)

    sample_ids = [s[0] for s in samples]
    feats = union_features(samples)

    # --- Matrices de abundancia relativa y (si hay) conteos ---
    rel_mat, cnt_mat = [], []
    have_counts = all(s[2] is not None for s in samples)
    for fid in feats:
        row_rel, row_cnt = [], []
        for _, rel_map, cnt_map in samples:
            row_rel.append(rel_map.get(fid, 0.0))
            if have_counts: row_cnt.append(cnt_map.get(fid, 0))
        rel_mat.append(row_rel)
        if have_counts: cnt_mat.append(row_cnt)

    # Normaliza columnas a suma 1
    col_sums = [0.0]*len(sample_ids)
    for row in rel_mat:
        for j, x in enumerate(row): col_sums[j] += x
    for j, s in enumerate(col_sums):
        if s>0:
            for i in range(len(rel_mat)): rel_mat[i][j] = rel_mat[i][j] / s

    # --- Exporta abundancias por feature ---
    rows = [[fid] + ["%.10f" % rel_mat[i][j] for j in range(len(sample_ids))]
            for i, fid in enumerate(feats)]
    write_csv(os.path.join(args.out, "abundance_relative_by_feature.csv"),
              ["feature_id"] + sample_ids, rows)

    # --- Alfa-diversidad ---
    alpha_rows = []
    for j, sid in enumerate(sample_ids):
        p = [rel_mat[i][j] for i in range(len(feats))]
        S_obs = sum(1 for x in p if x>0)
        H = shannon(p); sim1D = simpson_1D(p); J = pielou_evenness(H, S_obs)
        ch1 = ""
        if have_counts:
            counts = [cnt_mat[i][j] for i in range(len(feats))]
            ch1 = "%.6f" % chao1(counts)
        alpha_rows.append([sid, str(S_obs), "%.6f" % H, "%.6f" % sim1D, "%.6f" % J, ch1])
    write_csv(os.path.join(args.out, "alpha_diversity_metrics.csv"),
              ["sample_id","Observed","Shannon","Simpson_1D","Pielou","Chao1"], alpha_rows)

    # --- Distancias (beta) ---
    cols_rel = [[rel_mat[i][j] for i in range(len(feats))] for j in range(len(sample_ids))]
    def matrix_from(dist_fn):
        return [[dist_fn(cols_rel[a], cols_rel[b]) for b in range(len(sample_ids))]
                for a in range(len(sample_ids))]
    M_bray = matrix_from(bray_curtis)
    M_jacc = matrix_from(jaccard_binary)
    M_hell = matrix_from(hellinger_euclidean)
    write_csv(os.path.join(args.out, "beta_distance_bray.csv"), [""]+sample_ids,
              [[sample_ids[i]]+["%.6f"%v for v in row] for i,row in enumerate(M_bray)])
    write_csv(os.path.join(args.out, "beta_distance_jaccard.csv"), [""]+sample_ids,
              [[sample_ids[i]]+["%.6f"%v for v in row] for i,row in enumerate(M_jacc)])
    write_csv(os.path.join(args.out, "beta_distance_hellinger_euclidean.csv"), [""]+sample_ids,
              [[sample_ids[i]]+["%.6f"%v for v in row] for i,row in enumerate(M_hell)])

    # --- Composición por rank solicitado y HTML ---
    comp_samples, rank_abund = collect_rank_matrix(args.emu_dir, args.rank)
    means = [(k, sum(v)/max(1,len(v))) for k,v in rank_abund.items()]
    means.sort(key=lambda x: x[1], reverse=True)
    top = [k for k,_ in means[:args.top]]
    series = OrderedDict((k, rank_abund[k]) for k in top)
    other = [0.0]*len(comp_samples)
    for k,vals in rank_abund.items():
        if k in series: continue
        for i,v in enumerate(vals): other[i]+=v
    series["Other"] = other

    # SVGs
    shannon_vals = [float(x[2]) for x in alpha_rows]
    observed_vals = [float(x[1]) for x in alpha_rows]
    svg1 = svg_bar_chart("Alpha — Shannon", sample_ids, shannon_vals)
    svg2 = svg_bar_chart("Alpha — Observed (riqueza)", sample_ids, observed_vals)
    svg3 = svg_heatmap("Beta — Bray–Curtis (heatmap)", sample_ids, M_bray)
    title_rank = {"species":"especie","genus":"género","family":"familia"}[args.rank]
    svg4 = svg_stacked_bars(f"Abundancia relativa por {title_rank} (Top {args.top} + Other)", comp_samples, series)

    html = f"""<!DOCTYPE html>
<html lang="es"><head>
<meta charset="utf-8"/>
<title>Informe Microbiota (sin R)</title>
<style>
 body{{font-family:system-ui,-apple-system,Segoe UI,Roboto,Ubuntu,"Helvetica Neue",Arial,sans-serif;margin:20px;background:#fafafa;color:#222}}
 .card{{background:#fff;border:1px solid #eee;border-radius:12px;box-shadow:0 2px 8px rgba(0,0,0,.04);padding:16px;margin:16px auto;max-width:980px}}
 h2{{margin:8px 0 12px 0}}
 .note{{font-size:13px;color:#666}}
 code{{background:#f2f2f2;padding:2px 6px;border-radius:6px}}
</style>
</head><body>
<div class="card"><h2>Alpha diversidad</h2><div>{svg1}</div><div style="height:12px"></div><div>{svg2}</div>
<p class="note">Las barras están escaladas al máximo observado.</p></div>
<div class="card"><h2>Beta diversidad</h2><div>{svg3}</div>
<p class="note">Escala 0→1 (menor→mayor distancia).</p></div>
<div class="card"><h2>Composición ({title_rank})</h2><div>{svg4}</div>
<p class="note">Calculado desde *_rel-abundance.tsv de EMU; columnas normalizadas por muestra.</p></div>
<div class="card"><h2>Archivos</h2>
<ul>
<li><code>{os.path.abspath(os.path.join(args.out, "abundance_relative_by_feature.csv"))}</code></li>
<li><code>{os.path.abspath(os.path.join(args.out, "alpha_diversity_metrics.csv"))}</code></li>
<li><code>{os.path.abspath(os.path.join(args.out, "beta_distance_bray.csv"))}</code></li>
<li><code>{os.path.abspath(os.path.join(args.out, "beta_distance_jaccard.csv"))}</code></li>
<li><code>{os.path.abspath(os.path.join(args.out, "beta_distance_hellinger_euclidean.csv"))}</code></li>
</ul></div>
</body></html>"""
    os.makedirs(args.out, exist_ok=True)
    out_html = os.path.join(args.out, "report.html")
    with open(out_html, "w", encoding="utf-8") as f: f.write(html)

    print("Metricas y reporte generados en:", os.path.abspath(args.out))
    print("Abra:", out_html)

if __name__ == "__main__":
    main()
PY
```
```
chmod +x metricasd.py
./metricasd.py --emu-dir results_emu --out mp_outputs --rank species --top 12
```
