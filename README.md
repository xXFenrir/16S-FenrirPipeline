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

# BASECALLING

Este es el primer paso del pipeline, pues aquí se busca convertir los archivos `FAST5` a `FASTQ`. Los archivos `FAST5` son el output de la secuenciación con MinION (Oxford Nanopore Technologies). En este caso, ONT proporciona Dorado como una herramienta para este paso. 
Para usar esta herramienta es necesario hacer la instalación de la misma para entorno Linux. Esta descarga se hizo por medio de línea de comando desde el Terminal, pues `Dorado` se encuentra como un paquete de `Anaconda` ([Dorado](https://anaconda.org/HCC/dorado)). Se descarga el binario precompilado con:
```
curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.1.1-linux-x64.tar.gz" -o dorado-1.1.1-linux-x64.tar.gz
```

Luego, se debe extraer el archivo con:
```
tar -xzf dorado-1.1.1-linux-x64.tar.gz
```

Sin embargo, en las versiones actuales de `Dorado` no se utilizan como entradas los archivos en formato `FAST5` si no archivos `POD5`, por lo que antes de aplicar el comando básido de la herrmaienta se deben convertir los archivos de ser necesario.
```
pod5 convert fast5 /ruta/fast5/*.fast5 --output pod5_out/

dorado basecaller hac /ruta/pod5/ \
  --emit-fastq \
  --no-trim \
  -x auto > basecalls.fastq
```

- `hac` es para usar el modelo de alta precicsión.
- `--emit-fastq` hace que la salida sea en formato `FASTQ`.
- `--no-trim` es para evitar que recorte barcodes y adaptadores.
- `-x auto` decide automáticamente si usar GPU o CPU.

# DEMULTIPLEXING

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
- `-- trim all` quitar barcodes, adaptadores y primers. En este caso no afecta, pues dorado garantiza que el recorte no interfiere con la demultiplexación.
- `--primer-sequences` archivo `.fasta` von las secuencias de los primers.
- `-x auto` automáticamente elige si usar la GPU o CPU.
- `--no-clasiffy` no vuelve a clasificar, simplemente lee la asignación de barcode que ya quedó guardada en `calls.bam` durante el basecalling, y divide las lecturas en archivos separados por barcode.
- `--emit-fastq` la salida pasa a ser en formato `FASTQ`.
- `-o` carpeta donde estaran las salidas.

# REPORTE DE ESTADÍSTICAS BÁSICAS

Es necesario saber si las muestras necesitan o no ser limpiados, pues en el caso de algunos repositorios, estos ya viene filtrados. Por lo que, se deben tabular la estadísticas básicas de los archivos `FASTQ`, de tal manera que se pueda observar si vale la pena filtrar o no. Adicvionalmente, se podría llevar un registo sobre la calidad en la que estan las secuencias antes y después de filtradas. Por lo que, se usó el código que se encuentra en el archivo de [Reporte de estadísticas](https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/Reporte%20de%20estad%C3%ADsticas). Este código permite generar un archivo `TXT` y `XLSX` donde se genera una tabla con el reporte de estadísticas básicas de cada una de las secuencais.

## Paquetes importados

Los paquetes que se importanron fueron:
- `argparse` evita que se tenga que editar el código si se cambia el formato de la entradaa tanto con una carpeta como con un archivo.
- `gzip` lee `FASTQ.GZ` sin descomprimir.
- `os` para saber tamaño del archivo.
- `pathlib.Path` reconoce espacios, ~, etc.
- `re` permite reconocer solo las bases nitrogenadas en archivos `FASTA`.
- `sys` permite imprimir avisos.
- `time` permite conocer el tiempo que demoró en ejecutarse.
- `shutil` permite ver si es posible convertir un archivo de `TSV` a `XLSX`.
- `subprocess` permite separa cada columna en el archivo `XLSX`.
- `pandas` crea un `XLSX` con las columnas establecidas.

## Funciones definidas

- `is_gzip()`: Detecta si el archivo a analizar posee el sufijo `.GZ`, pues convierte el string en booleano y hace la lectura de este.
- `open_maybe_gzip()`: Abre un archivo de texto, si es `.GZ`, usa `gzip.open()` para leer el archivo sin descomprimir, si no, usa `open()`.
- `iter_fastq_reads()`: Itera continuamente sobre un archivo `FASTQ` guardando únicamente la secuencia y su calidad.
- `n50_from_lengths()`: Ordena las longitudes de mayor a menor y acumula hasta alcanzar al menos el 50 % del total de bases, y esta longitud corresponde al N50. En caso de que no haya datos lo representa como 0.
- `read_fasta_seqs()`: Carga un archivo `FASTA` y devuelve solo las secuencias. Para cada `FASTA`, concatena líneas de secuencia y aplica una limpieza con `regex` para eliminar cualquier carácter que no sea ACGTN. De esta manera, se estandariza la entrada para la detección de primers.
- `revcomp()`: Genera la secuencia reverse de ADN. Primero traduce, y luego invierte la cadena, de tal manera, permite buscar primers en ambas hebras.
- `hamming_leq_k()`: Con una matriz, comprueba que la longitud del fragmento de la secuencia corresponda a la del primer. Además, se usa para saber si hay secuencias N, permitiendo la lectura aún cuando  haya errores de secuenciación.
- `any_primer_in_window()`: Comprueba si en la lectura forward o reverse aparece algún primer. Donde, si encuentra una coincidencia es `True` y `False` si no hay ninguna.
- `process_fastq()`: Procesa un archivo `FASTQ` y calcula métricas de número de lecturas, bases totales, longitudes, GC%, QScore promedio y porcentajes de bases ≥Q20, ≥Q30 y % primers.
- `_fmt()`: Si el valor es float, fija el número de decimales, si es entero, lo convierte a string. Permitiendo que el formato sea el mismo en todas las columnas.
- `row_to_display()`: Transforma el conjunto de métricas en una lista de strings y en el orden exacto de columnas para el `TXT`.
- `row_to_excel()`: Convierte el mismo conjunto de métricas en int/float/None, lo que, permite que `pandas` u `openpyxl` escriban un `XLSX`.
- `write_pretty_table()`: Calcula el ancho máximo de cada columna, rellena con `ljust` e inserta una línea de guiones del mismo largo para cada columna y cada fila.
- `write_tsv()`: Genera un `TSV` separado por tabulaciones usando los valores de `row_to_excel`.
- `try_write_xlsx()`: Intenta escribir el `XLSX`. Primero intenta con `pandas`, si falla, usa `openpyxl`, y si tampoco es posible, utiliza LibreOffice headless convirtiendo un `TSV` de `write_tsv()` a `XLSX`. Finalmente, indica `True` si logró crear el Excel y muestra el método usado.
- `find_fastqs()`: Permite comprobar la ruta sin considerar jerarquías, pues si es un archivo, valida que tenga extensión `FASTQ` o `FQ`. Si es carpeta, usa `glob` o `rglob` según si está o no comprimido. Por último, devuelve las rutas ordenadas y sin duplicados.
- `parse_args()`: Define la interfaz de línea de comandos, como la entrada, la salida y ruta de los primers.
- `main()`: Parsea argumentos y anuncia inicio, carga primers si se proporcionan, localiza los `FASTQ`, informa cuántos encontró, procesa cada archivo con `process_fastq`, escribe el `TXT` con `write_pretty_table` y luego intenta el `XLSX` con `try_write_xlsx`. Finalmente, imprime un resumen con la ruta de salida y el tiempo de ejecución.

## Ejecución del código

Una vez creado el script con el código, se ejecúta:
```
python3 "/home/fenrir/scriptsbioinf/fastq_estads.py" \
  "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels" \
  -r \
  --primers "/home/fenrir/resources/primers_stesen.fasta" \
  -o "/home/fenrir/og_stats/estadisticas_og_stesen.txt" \
  --xlsx-out "/home/fenrir/og_stats/estadisticas_og_stesen.xlsx"
```
1. Primero se escribe la ruta en la que se encuentra el script.
2. `-r` busca en subcarpetas.
3. `--primers` es la ruta del archivo `FASTA` de los primers.
4. `-o` ruta del archivo `TXT` de salida.
5. `--xlsx-out` rutas del archivo `XLSX` de salida. 

# Denoising y Trimming

Es uno de los pasos más importantes en un análisis bioinformático de datos de secuenciación 16S rRNA es la depuración de lecturas crudas. Pues el objetivo aquí es mejorar la calidad de los datos y asegurar que únicamente las lecturas confiables y relevantes pasen a la etapa de taxonomía. En este caso se quiere:

- Filtrar lecturas por longitud y calidad, eliminando aquellas demasiado cortas, largas o de baja calidad que puedan corresponder a artefactos de secuenciación o fragmentos incompletos.
- Orientar y recortar las lecturas basándonos en la ubicación de los primers, asegurando que todas las secuencias tengan la misma dirección y contenido correcto.

## Pychopper

Una de las herramientas encontradas en bibliografía para la limpieza de las muestras es `Pychopper`. La instalación se hace por línea de código desde desde [GitHub](https://github.com/epi2me-labs/pychopper) con:
```
conda install -c nanoporetech -c conda-forge -c bioconda "nanoporetech::pychopper"
```

Pychopper identifica los primers en cada lectura, determina su orientación y recorta el amplicón para dejar únicamente la región de interés. Además, genera reportes que muestran cuántas lecturas fueron clasificadas, rechazadas o rescatadas. Adicionalmente, es una herramienta que permite filtrar aquellas lecturas que poseen un QScore por debajo de 9. 

Sin embargo, es necesario crear un archivo `FASTA` con la secuencia de los primers y otro `TXT` con la orientación de los mismos. De esta manera, ya es posible hacer uso del comando base de `Pychopper`.
```
pychopper -m edlib \ 
  -b "/home/fenrir/resources/primers_stesen.fasta" \ 
  -c "/home/fenrir/resources/primers_stesen.txt" \ 
  -Q 9 -z 1300 -t 8 \
  -Y 0 -q 0.52 \ 
  "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels/SRR26147154.fastq.gz" \ 
  "results/sample1_oriented_trimmed.fastq" 
```

- `-m edlib` habilita un alineador para encontrar los primers.
- `-b` ruta del archivo con la secuencia de los primers.
- `-c` ruta del archivo con la orientación de los primers.
- `-Q` filtro por QScore.
- `-z`filtro por longitud mínima.
- `-t` núcleos de CPU a usar.
- `-Y` desactiva la función de muestreo automático para hacer el proceso reproducible.
- `-q` ajusta la exigencia del alineamiento a los primers.

Las dos últimas líneas de código permiten buscar el archivo `FASTQ` a limpiar. Luego, se define la ruta en que se va a guardar el archivo `FASTQ` limpiado.

## Filtlong

Sin embargo, `Pychopper` es una herramienta que filtra por longitud mínima y no por rango. Por lo que, se usó de forma complementaria la herramienta `Filtlong`, pues esta permite establecer el rango ideal en pb. 

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

## Descripción del código

El script para el [Denosing y trimming](https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/Denoising%20y%20trimming) permite filtrar las secuencias `FASTQ` crudas por longitud, QScore y remover primers.

### Paquetes importados

Los paquetes importados fueron:
- `argparse`: Permite construir la interfaz de línea de comandos definiendo las entradas y variables.
- `glob`: Sirve para encontrar archivos con sufijo `FASTQ` independientemente de si está o no comprimido. En caso de no encontrar nada, el programa avisa y termina. 
- `os`: Verifica la existencia de los archivos de entrada, crea carpetas de salida por muestra, compone rutas portables y elimina salidas parciales si la ejecución falla.
- `shlex`: Permite que las rutas o archivos que estén compuestos con espacios sean tomados con una sola palabra con `shlex.quote()`.
- `shutil.which`: Comprueba si falta `filtlong`, `pychopper` y `bash` están disponibles antes de arrancar, generando un mensaje si falta alguno de estos.
- `subprocess`: Arma un proceso `filtlong … | pychopper …` y lo corre con `subprocess.run()`. Se antepone `set -o pipefail` para que cualquier fallo en el proceso se refleje sea informado. Luego, se revisa `returncode` para decidir si se salta una muestra y se limpia la salida parcial.
- `sys`: Proporciona control sobre la salida y la terminación del programa por medio de avisos.
- `pathlib.Path`: Obtiene el nombre base del archivo quitando sus sufijos.
- `typing`: Permite que se informe de errores antes de correr el script.

### Funciones definidas

- `_must_exist(path, kind)`: Valida que una ruta exista antes de correr el pipeline. Si `kind` es `dir`, exige que la entrada sea una carpeta, si es `file`, exige que sea un archivo. Cuando la comprobación falla, termina el programa con un mensaje claro usando `sys.exit.
- `_which_or_die(cmd)`: Verificar la presencia de `filtlong`, `pychopper` y `bash`. Si alguna no está instalada, finaliza con un mensaje explicando qué falta.
- `_sample_name(f)`: Toma la ruta de un archivo `FASTQ` comprimido o no, y obtiene el nombre base con `Path(f).name` y elimina, los sufijos que le acompañen. De este modo, se usa para crear el subdirectorio de salida y el nombre del archivo resultante, manteniendo una nomenclatura consistente.
- `run_dentrim()`: Valida entradas y dependencias, localiza todos los archivos de lectura en `input_dir`. Pues, ejecuta `filtlong` leyendo directamente el archivo para filtrar por longitud entre `minlen` y `maxlen`. Luego, envía su salida por `stdin` a `pychopper`, que recibe y orienta usando los primers con `-b primers` y `-c pconfig`, el umbral de calidad qscore, y el número de hilos threads. Además, se usa `-z minlen` para mantener coherencia con el mínimo de longitud. Antes de su ejecución se activa `set -o pipefail` para que cualquier fallo intermedio se refleje en la terminal.
- `_build_parser()`: Define todas las variables de entrada como obligatorias, para definir facilmente desde la interfaz de línea de comandos.
- `main()`: Actúa como punto de entrada cuando el archivo se ejecuta directamente. Donde, interpreta los argumentos proporcionados por el usuario y da paso a `run_dentrim()`, activando mensajes de salida que muestran los parámetros y el avance.

### Ejecución del código

Para ejecutar el código se usa el siguiente comando:
```
/home/fenrir/scriptsbioinf/dentrim.py \
  -i "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels" \
  -o "dentrim_stesen" \
  --minlen 1000 --maxlen 1700 \
  -Q 12 -t 8 \
  --primers "/home/fenrir/resources/primers_stesen.fasta" \
  --pconfig "/home/fenrir/resources/primers_stesen.txt"
```

1. Se define la ruta del script.
2. `-i` es la ruta de los datos de entrada.
3. `-o` nombre de la carpeta de salida.
4. `--minlen` y `--maxlen` definen el rango de longitud aceptada.
5. `-Q` es el QScore mínimo aceptado.
6. `-t` es el número de núcleos para ejecutar la herramienta.
7. `--primers` es la ruta con el `FASTA` con la secuencia de los primers.
8. `--pconfig` es la ruta con la dirección de los primers.

# Taxonomía
EMU
```
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EMU pipeline wrapper + aggregator (+ opcional Bray–Curtis/Jaccard + PCoA)

Novedades:
- --do-pcoa: calcula matrices Bray–Curtis/Jaccard y PCoA usando feature_table_counts.tsv
- --pcoa-outdir: directorio para los outputs de PCoA (default: <outdir>/BETA_PCOA)
- --pcoa-relative: normaliza a abundancias relativas para Bray–Curtis

Requisitos:
- Python 3.8+
- pandas, numpy
- (opcional) biom-format para exportar BIOM
- (para PCoA) scikit-bio, matplotlib
- EMU instalado en PATH o especificar --emu-cmd

Ejemplo:
python EMU.py \
  --db /home/fenrir/emu_db \
  --input-dir /home/fenrir/results_dentrim_Q10 \
  --pattern '*_final.fastq' \
  --outdir /home/fenrir/results_dentrim_Q10/EMU_Q10 \
  --threads 8 \
  --rank species \
  --keep-counts \
  --keep-assignments \
  --emu-cmd /home/fenrir/anaconda3/envs/pipelinefenrir/bin/emu \
  --do-pcoa --pcoa-relative
"""
import argparse
import os
import sys
import glob
import subprocess
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd

# ----------------------------- Utilidades básicas -----------------------------

TAX_COLS_CANON = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
RANKS_ALLOWED = set(TAX_COLS_CANON)

def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def which(cmd: str) -> Optional[str]:
    from shutil import which as _which
    return _which(cmd)

def safe_basename_noext(p: Path) -> str:
    """Obtén un nombre de muestra razonable desde el FASTQ."""
    name = p.name
    import re
    m = re.search(r"(SRR\d{6,})", name, flags=re.IGNORECASE)
    if m:
        return m.group(1)
    stem = p.stem
    if stem.endswith(".fastq"):
        stem = Path(stem).stem
    parent = p.parent.name
    return f"{parent}_{stem}"

def read_table_maybe(path: Path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as e:
        eprint(f"[WARN] No pude leer {path}: {e}")
        return None

def to_qiime_tax_string(row: pd.Series) -> str:
    parts = []
    prefixes = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
    for col, pref in zip(TAX_COLS_CANON, prefixes):
        val = str(row.get(col, "") or "").strip()
        if not val or val.lower() == "unassigned":
            parts.append(pref)
        else:
            parts.append(pref + val)
    return ";".join(parts)

def chao1(counts: np.ndarray) -> float:
    counts = counts[counts > 0]
    S_obs = (counts > 0).sum()
    if counts.size == 0:
        return 0.0
    f1 = (counts == 1).sum()
    f2 = (counts == 2).sum()
    if f2 == 0:
        return float(S_obs + (f1 * (f1 - 1)) / 2.0)
    return float(S_obs + (f1 * f1) / (2.0 * f2))

def shannon_entropy(p: np.ndarray) -> float:
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    return float(-(p * np.log(p)).sum())

def simpson_index(p: np.ndarray) -> float:
    if p.size == 0:
        return 0.0
    return float(1.0 - (p * p).sum())

# ----------------------------- EMU ejecución ---------------------------------

def run_emu_for_sample(
    emu_cmd: str,
    db: Path,
    fastq: Path,
    outdir_sample: Path,
    sample: str,
    threads: int,
    keep_counts: bool,
    keep_assignments: bool,
) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    ensure_dir(outdir_sample)
    cmd = [
        emu_cmd, "abundance",
        "--db", str(db),
        "--threads", str(threads),
        "--output-dir", str(outdir_sample),
        "--output-basename", sample,
    ]
    if keep_counts:
        cmd.append("--keep-counts")
    if keep_assignments:
        cmd.append("--keep-read-assignments")
    cmd.append(str(fastq))

    eprint(f"[EMU] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    rel_abund = next(outdir_sample.glob(f"{sample}*rel-abundance.tsv"), None)
    counts = next(outdir_sample.glob(f"{sample}*counts.tsv"), None) if keep_counts else None
    assigns = next(outdir_sample.glob(f"{sample}*read-assignments.tsv"), None) if keep_assignments else None
    return rel_abund, counts, assigns

# ----------------------------- Agregación ------------------------------------

def select_rank(df: pd.DataFrame, rank: str) -> pd.DataFrame:
    lower = {c.lower(): c for c in df.columns}
    if rank not in lower:
        return df.copy()
    col = lower[rank]
    df2 = df.copy()
    df2 = df2[df2[col].notna() & (df2[col].astype(str).str.strip() != "")]
    return df2

def build_feature_id(row: pd.Series, rank: str) -> str:
    tax_id = str(row.get("tax_id", "")).strip()
    name = ""
    if rank in row and pd.notna(row[rank]) and str(row[rank]).strip():
        name = str(row[rank]).strip()
    elif "species" in row and pd.notna(row["species"]) and str(row["species"]).strip():
        name = str(row["species"]).strip()
    elif "tax_name" in row and pd.notna(row["tax_name"]) and str(row["tax_name"]).strip():
        name = str(row["tax_name"]).strip()
    name = name.replace("|", "_")
    tax_id = tax_id.replace("|", "_")
    base = f"{rank}|{tax_id}|{name}" if tax_id or name else f"{rank}|NA|NA"
    return base

def aggregate_tables(
    per_sample_counts: Dict[str, Path],
    per_sample_rel: Dict[str, Path],
    outdir: Path,
    rank: str,
    min_abundance: float = 0.0
) -> Tuple[Path, Path, Path, Path]:
    features = set()
    sample_order = sorted(set(per_sample_counts.keys()) | set(per_sample_rel.keys()))
    rows_tax_meta: Dict[str, Dict[str, str]] = {}
    counts_map: Dict[str, Dict[str, int]] = {s: {} for s in sample_order}
    rel_map: Dict[str, Dict[str, float]] = {s: {} for s in sample_order}

    def process_one_tsv(tsv_path: Path, sample: str, kind: str):
        df = read_table_maybe(tsv_path)
        if df is None or df.empty:
            return
        col_low = {c.lower(): c for c in df.columns}
        for col in TAX_COLS_CANON:
            if col not in col_low:
                df[col] = np.nan
            else:
                df[col] = df[col_low[col]]

        df = select_rank(df, rank)

        if kind == "rel":
            rel_col = None
            for cand in ["relative_abundance", "rel_abundance", "rel-abundance", "abundance", "fraction_total_reads", "fraction"]:
                if cand in df.columns:
                    rel_col = cand
                    break
            if rel_col is None:
                eprint(f"[WARN] No hallé columna de abundancia relativa en {tsv_path}.")
            else:
                if min_abundance > 0:
                    df = df[df[rel_col] >= min_abundance]
        return df

    for sample, tsv in per_sample_rel.items():
        df = process_one_tsv(tsv, sample, "rel")
        if df is None or df.empty:
            continue
        rel_col = None
        for cand in ["relative_abundance", "rel_abundance", "rel-abundance", "abundance", "fraction_total_reads", "fraction"]:
            if cand in df.columns:
                rel_col = cand
                break
        if rel_col is None:
            continue

        for _, row in df.iterrows():
            fid = build_feature_id(row, rank)
            features.add(fid)
            rel_map[sample][fid] = float(row[rel_col]) if pd.notna(row[rel_col]) else 0.0
            rows_tax_meta.setdefault(fid, {c: "" for c in TAX_COLS_CANON})
            for c in TAX_COLS_CANON:
                val = row.get(c, "")
                rows_tax_meta[fid][c] = "" if pd.isna(val) else str(val)

    for sample, tsv in per_sample_counts.items():
        df = process_one_tsv(tsv, sample, "counts")
        if df is None or df.empty:
            continue
        count_col = None
        for cand in ["count", "counts", "read_count", "reads"]:
            if cand in df.columns:
                count_col = cand
                break
        if count_col is None:
            eprint(f"[WARN] No hallé columna de conteo en {tsv}.")
            continue

        for _, row in df.iterrows():
            fid = build_feature_id(row, rank)
            features.add(fid)
            val = row[count_col]
            try:
                val = int(val)
            except Exception:
                try:
                    val = int(float(val))
                except Exception:
                    val = 0
            counts_map[sample][fid] = val
            rows_tax_meta.setdefault(fid, {c: "" for c in TAX_COLS_CANON})
            for c in TAX_COLS_CANON:
                valx = row.get(c, "")
                if rows_tax_meta[fid].get(c, "") == "":
                    rows_tax_meta[fid][c] = "" if pd.isna(valx) else str(valx)

    features = sorted(features)
    mat_counts = np.zeros((len(features), len(sample_order)), dtype=int)
    for j, s in enumerate(sample_order):
        for i, fid in enumerate(features):
            mat_counts[i, j] = int(counts_map[s].get(fid, 0))

    mat_rel = np.zeros((len(features), len(sample_order)), dtype=float)
    for j, s in enumerate(sample_order):
        for i, fid in enumerate(features):
            mat_rel[i, j] = float(rel_map[s].get(fid, 0.0))

    df_counts = pd.DataFrame(mat_counts, index=features, columns=sample_order)
    df_rel = pd.DataFrame(mat_rel, index=features, columns=sample_order)

    tax_rows = []
    for fid in features:
        meta = rows_tax_meta.get(fid, {c: "" for c in TAX_COLS_CANON})
        tax_rows.append({"feature_id": fid, "taxonomy": to_qiime_tax_string(pd.Series(meta))})
    df_tax = pd.DataFrame(tax_rows, columns=["feature_id", "taxonomy"])

    alpha_rows = []
    for s in sample_order:
        counts_vec = df_counts[s].to_numpy()
        rel_vec = df_rel[s].to_numpy()
        observed = int((counts_vec > 0).sum()) if counts_vec.sum() > 0 else int((rel_vec > 0).sum())
        chao = chao1(counts_vec) if counts_vec.sum() > 0 else float(observed)
        if counts_vec.sum() > 0:
            p = counts_vec / counts_vec.sum()
        else:
            total_rel = rel_vec.sum()
            p = rel_vec / total_rel if total_rel > 0 else rel_vec
        shan = shannon_entropy(p)
        simp = simpson_index(p)
        alpha_rows.append({"sample": s, "observed": observed, "chao1": chao, "shannon": shan, "simpson": simp})
    df_alpha = pd.DataFrame(alpha_rows, columns=["sample", "observed", "chao1", "shannon", "simpson"])

    counts_tsv = outdir / "feature_table_counts.tsv"
    rel_tsv = outdir / "feature_table_relabund.tsv"
    tax_tsv = outdir / "taxonomy.tsv"
    alpha_tsv = outdir / "alpha_diversity.tsv"
    df_counts.to_csv(counts_tsv, sep="\t", index=True, header=True)
    df_rel.to_csv(rel_tsv, sep="\t", index=True, header=True)
    df_tax.to_csv(tax_tsv, sep="\t", index=False)
    df_alpha.to_csv(alpha_tsv, sep="\t", index=False)
    eprint(f"[OK] Matriz de conteos: {counts_tsv}")
    eprint(f"[OK] Matriz de abundancia relativa: {rel_tsv}")
    eprint(f"[OK] Taxonomía: {tax_tsv}")
    eprint(f"[OK] Alfa diversidad: {alpha_tsv}")

    try:
        import biom
        from biom.table import Table
        biom_out = outdir / "feature_table.biom"
        table = Table(df_counts.to_numpy(), observation_ids=df_counts.index.tolist(), sample_ids=df_counts.columns.tolist())
        with biom.util.biom_open(str(biom_out), 'w') as f:
            table.to_hdf5(f, "emu_pipeline_pack")
        eprint(f"[OK] BIOM: {biom_out}")
    except Exception as e:
        eprint(f"[INFO] Saltando export BIOM (instala 'biom-format' si lo necesitas). Motivo: {e}")

    return counts_tsv, rel_tsv, tax_tsv, alpha_tsv

# ----------------------------- Beta + PCoA -----------------------------------

def compute_beta_pcoa(counts_tsv: Path, outdir: Path, use_relative: bool = False):
    """
    Calcula Bray–Curtis y Jaccard + PCoA a partir de feature_table_counts.tsv
    Escribe:
      - braycurtis_dm.tsv, jaccard_dm.tsv
      - pcoa_braycurtis_coords.tsv / eigvals.tsv / variance.tsv (+ PNG)
      - pcoa_jaccard_coords.tsv / eigvals.tsv / variance.tsv (+ PNG)
    """
    ensure_dir(outdir)
    eprint(f"[PCOA] Cargando tabla: {counts_tsv}")
    df = pd.read_csv(counts_tsv, sep="\t", index_col=0).fillna(0)
    # filtrar vacíos
    df = df[(df.sum(axis=1) > 0)]
    df = df.loc[:, (df.sum(axis=0) > 0)]
    if df.shape[1] < 2:
        raise ValueError("Se necesitan al menos 2 muestras con conteos > 0 para PCoA.")

    if use_relative:
        df_bray = df / df.sum(axis=0).replace(0, np.nan)
        df_bray = df_bray.fillna(0.0)
    else:
        df_bray = df.copy()

    df_jacc = (df > 0).astype(int)

    try:
        from skbio.diversity import beta_diversity
        from skbio.stats.ordination import pcoa
    except Exception as e:
        raise RuntimeError("Falta scikit-bio. Instala con: conda install -c conda-forge scikit-bio") from e

    X_bray = df_bray.T.values
    X_jacc = df_jacc.T.values
    sample_ids = df_bray.columns.astype(str).tolist()

    eprint("[PCOA] Bray–Curtis…")
    dm_bray = beta_diversity(metric="braycurtis", counts=X_bray, ids=sample_ids)
    pd.DataFrame(dm_bray.data, index=dm_bray.ids, columns=dm_bray.ids).to_csv(outdir / "braycurtis_dm.tsv", sep="\t")

    eprint("[PCOA] Jaccard…")
    dm_jacc = beta_diversity(metric="jaccard", counts=X_jacc, ids=sample_ids)
    pd.DataFrame(dm_jacc.data, index=dm_jacc.ids, columns=dm_jacc.ids).to_csv(outdir / "jaccard_dm.tsv", sep="\t")

    eprint("[PCOA] Ordination Bray–Curtis…")
    pcoa_bray = pcoa(dm_bray)
    _write_pcoa_outputs(pcoa_bray, "braycurtis", outdir)

    eprint("[PCOA] Ordination Jaccard…")
    pcoa_jacc = pcoa(dm_jacc)
    _write_pcoa_outputs(pcoa_jacc, "jaccard", outdir)

    # figuras opcionales
    try:
        import matplotlib.pyplot as plt
        _plot_pcoa(pcoa_bray, "braycurtis", outdir)
        _plot_pcoa(pcoa_jacc, "jaccard", outdir)
    except Exception as e:
        eprint(f"[INFO] No se generarán PNG (matplotlib no disponible): {e}")

    eprint(f"[PCOA] Listo. Resultados en: {outdir}")

def _write_pcoa_outputs(ord_res, prefix: str, outdir: Path):
    coords = ord_res.samples.copy()
    coords.index.name = "sample"
    coords.to_csv(outdir / f"pcoa_{prefix}_coords.tsv", sep="\t")
    ev = pd.Series(ord_res.eigvals, name="eigenvalue")
    ev.to_csv(outdir / f"pcoa_{prefix}_eigvals.tsv", sep="\t", header=True)
    var = pd.Series(ord_res.proportion_explained, name="proportion_explained")
    var.to_csv(outdir / f"pcoa_{prefix}_variance.tsv", sep="\t", header=True)

def _plot_pcoa(ord_res, prefix: str, outdir: Path):
    import matplotlib.pyplot as plt
    coords = ord_res.samples
    if coords.shape[1] < 2:
        eprint("[INFO] PCoA con <2 ejes; omito gráfica.")
        return
    x, y = coords.iloc[:, 0], coords.iloc[:, 1]
    var = ord_res.proportion_explained
    xlab = f"PC1 ({var.iloc[0]*100:.1f}%)"
    ylab = f"PC2 ({var.iloc[1]*100:.1f}%)"
    plt.figure(figsize=(6, 5))
    plt.scatter(x, y)
    for sid, xi, yi in zip(coords.index, x, y):
        plt.text(xi, yi, str(sid), fontsize=8, ha="center", va="bottom")
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(f"PCoA — {prefix}")
    plt.tight_layout()
    out_png = outdir / f"pcoa_{prefix}.png"
    plt.savefig(out_png, dpi=150)
    plt.close()
    eprint(f"[OK] Figura: {out_png}")

# ----------------------------- Descubrimiento --------------------------------

def discover_fastqs(args) -> List[Path]:
    paths: List[Path] = []
    if args.input_glob:
        for pat in args.input_glob:
            for p in glob.glob(pat):
                if os.path.isfile(p):
                    paths.append(Path(p))
    if args.input_dir:
        base = Path(args.input_dir)
        pat = args.pattern or "*.fastq*"
        paths.extend(base.rglob(pat))
    uniq = sorted(set([p.resolve() for p in paths]))
    return uniq

def merge_assignments(assign_paths: Dict[str, Path], outdir: Path) -> Optional[Path]:
    if not assign_paths:
        return None
    rows = []
    nfiles = len(assign_paths)
    eprint(f"[INFO] Fusionando lecturas asignadas de {nfiles} muestras… (puede ser pesado)")
    for sample, p in assign_paths.items():
        df = read_table_maybe(p)
        if df is None or df.empty:
            continue
        colmap = {c.lower(): c for c in df.columns}
        read_col = colmap.get("read_id", None)
        taxid_col = colmap.get("tax_id", None)
        tname_col = colmap.get("tax_name", None)
        if read_col is None or taxid_col is None:
            eprint(f"[WARN] {p} no contiene columnas esperadas (read_id/tax_id). Lo omito.")
            continue
        sub = pd.DataFrame({
            "sample": sample,
            "read_id": df[read_col].astype(str),
            "tax_id": df[taxid_col].astype(str),
            "tax_name": df[tname_col].astype(str) if tname_col else ""
        })
        rows.append(sub)
    if not rows:
        return None
    big = pd.concat(rows, ignore_index=True)
    outcsv = outdir / "read_assignments_merged.csv"
    big.to_csv(outcsv, index=False)
    eprint(f"[OK] Lecturas asignadas fusionadas: {outcsv} (filas={len(big)})")
    return outcsv

# ----------------------------- Main ------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Pipeline EMU + agregación + (opcional) Bray/Jaccard + PCoA.")
    ap.add_argument("--db", required=True, help="Ruta a la base de datos de EMU")
    ap.add_argument("--outdir", required=True, help="Directorio de salida")
    ap.add_argument("--threads", type=int, default=8, help="Hilos para EMU (default: 8)")
    ap.add_argument("--emu-cmd", default="emu", help="Comando EMU (default: 'emu')")
    ap.add_argument("--input-glob", nargs="+", help="Uno o más patrones glob (ej: 'results/*/*.fastq.gz')")
    ap.add_argument("--input-dir", help="Carpeta raíz para buscar FASTQ(s)")
    ap.add_argument("--pattern", help="Patrón (default: '*.fastq*') para --input-dir")
    ap.add_argument("--keep-counts", action="store_true", help="Guardar tablas de conteos por muestra")
    ap.add_argument("--keep-assignments", action="store_true", help="Guardar lecturas asignadas por muestra")
    ap.add_argument("--rank", default="species", choices=list(RANKS_ALLOWED), help="Nivel taxonómico a consolidar (default: species)")
    ap.add_argument("--min-abundance", type=float, default=0.0, help="Filtro mínimo de abundancia relativa para incluir features (default: 0.0)")
    ap.add_argument("--merge-assignments", action="store_true", help="Fusionar lecturas asignadas en un solo CSV maestro")
    # PCoA
    ap.add_argument("--do-pcoa", action="store_true", help="Calcular Bray–Curtis/Jaccard y PCoA a partir de la tabla generada")
    ap.add_argument("--pcoa-outdir", help="Directorio de salida para PCoA (default: <outdir>/BETA_PCOA)")
    ap.add_argument("--pcoa-relative", action="store_true", help="Usar abundancias relativas para Bray–Curtis en PCoA")
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    ensure_dir(outdir)

    emu_bin = args.emu_cmd
    if os.path.sep not in emu_bin:
        wb = which(emu_bin)
        if wb is None:
            eprint(f"[ERROR] No se encontró '{emu_bin}' en PATH. Especifica --emu-cmd o ajusta tu entorno.")
            sys.exit(1)
        emu_bin = wb

    db = Path(args.db).resolve()
    if not db.exists():
        eprint(f"[ERROR] Base de EMU no existe: {db}")
        sys.exit(1)

    fastqs = discover_fastqs(args)
    if not fastqs:
        eprint("[ERROR] No se encontraron FASTQ(s). Usa --input-glob o --input-dir/--pattern.")
        sys.exit(1)

    eprint(f"[INFO] FASTQ(s) detectados: {len(fastqs)}")
    manifest_rows = []

    per_sample_rel: Dict[str, Path] = {}
    per_sample_counts: Dict[str, Path] = {}
    per_sample_assigns: Dict[str, Path] = {}

    for fq in fastqs:
        sample = safe_basename_noext(fq)
        sdir = outdir / sample
        ensure_dir(sdir)
        try:
            rel, cnt, asg = run_emu_for_sample(
                emu_cmd=emu_bin, db=db, fastq=fq, outdir_sample=sdir, sample=sample,
                threads=args.threads, keep_counts=args.keep_counts, keep_assignments=args.keep_assignments
            )
        except subprocess.CalledProcessError as e:
            eprint(f"[ERROR] EMU falló para {fq}: {e}")
            continue

        if rel and rel.exists():
            per_sample_rel[sample] = rel
        if cnt and cnt.exists():
            per_sample_counts[sample] = cnt
        if asg and asg.exists():
            per_sample_assigns[sample] = asg

        manifest_rows.append({
            "sample": sample,
            "fastq": str(fq),
            "outdir_sample": str(sdir),
            "rel_abundance_tsv": str(rel) if rel else "",
            "counts_tsv": str(cnt) if cnt else "",
            "assignments_tsv": str(asg) if asg else ""
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = outdir / "manifest.tsv"
    if not manifest.empty:
        manifest.to_csv(manifest_path, sep="\t", index=False)
        eprint(f"[OK] Manifiesto: {manifest_path}")
    else:
        eprint("[ERROR] No hay resultados para agregar. Revisa logs anteriores.")
        sys.exit(1)

    if not per_sample_rel:
        eprint("[ERROR] No se hallaron archivos *rel-abundance.tsv*. ¿EMU generó salidas?")
        sys.exit(1)

    counts_tsv, rel_tsv, tax_tsv, alpha_tsv = aggregate_tables(
        per_sample_counts=per_sample_counts,
        per_sample_rel=per_sample_rel,
        outdir=outdir,
        rank=args.rank,
        min_abundance=args.min_abundance
    )

    if args.merge_assignments and per_sample_assigns:
        merge_assignments(per_sample_assigns, outdir)

    # ---------- PCoA opcional ----------
    if args.do_pcoa:
        pcoa_dir = Path(args.pcoa_outdir).resolve() if args.pcoa_outdir else (outdir / "BETA_PCOA")
        ensure_dir(pcoa_dir)
        try:
            compute_beta_pcoa(counts_tsv=counts_tsv, outdir=pcoa_dir, use_relative=args.pcoa_relative)
        except Exception as e:
            eprint(f"[ERROR] Falló cálculo de PCoA: {e}")
            # no abortamos el pipeline principal
    # -----------------------------------

    eprint("\n=== Listo. Archivos principales ===")
    eprint(f"- {manifest_path}")
    eprint(f"- {counts_tsv}")
    eprint(f"- {rel_tsv}")
    eprint(f"- {tax_tsv}")
    eprint(f"- {alpha_tsv}")
    if args.do_pcoa:
        eprint(f"- Carpeta PCoA: {pcoa_dir}")
    eprint("Para UniFrac/PD necesitas un árbol filogenético externo (QIIME2/phylogeny, etc.).")
    return 0

if __name__ == "__main__":
    sys.exit(main())
```
```
python /home/fenrir/scriptsbioinf/EMU.py \
  --db /home/fenrir/emu_db \
  --input-dir /home/fenrir/results_dentrim_Q10 \
  --pattern '*_final.fastq' \
  --outdir /home/fenrir/results_dentrim_Q10/EMU_Q10 \
  --threads 8 \
  --rank species \
  --keep-counts \
  --keep-assignments \
  --emu-cmd /home/fenrir/anaconda3/envs/pipelinefenrir/bin/emu \
  --do-pcoa --pcoa-relative
```
