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

# Basecalling y Demultiplexing

El basecalling es el primer paso obligatorio del pipeline: convierte las señales eléctricas crudas capturadas por la plataforma Oxford Nanopore (formato `POD5`, sucesor de `FAST5`) en lecturas de secuencia. Sin este paso no existen datos de secuencia sobre los que trabajar. El demultiplexing, por su parte, separa las lecturas según el barcode de cada muestra para poder procesarlas de forma independiente en las etapas siguientes.

Ambos pasos se resuelven con **Dorado**, la herramienta de basecalling recomendada actualmente por Oxford Nanopore Technologies (sucesora de Guppy/Albacore). El binario precompilado se instala así:
```
curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-2.1.0-linux-x64.tar.gz" -o dorado-2.1.0-linux-x64.tar.gz
tar -xzf dorado-2.1.0-linux-x64.tar.gz
dorado-2.1.0-linux-x64/bin/dorado --version
```

Si la secuenciación quedó en formato `FAST5` en vez de `POD5`, se convierte antes de usar Dorado:
```
pod5 convert fast5 /ruta/fast5/*.fast5 --output pod5_out/
```

Siguiendo la recomendación del fabricante, basecalling y demultiplexing se ejecutan en un solo comando: se indica el kit de barcodes con `--kit-name` para que Dorado anote el barcode de cada lectura directamente en el BAM de salida, en vez de correr un demultiplexado por separado con otra herramienta (como `Porechop`). También se aplica un filtro de calidad mínima (`--min-qscore`) y se usa el modelo de alta precisión **HAC** (High Accuracy).

Comando final usado:
```
~/dorado-2.1.0-linux-x64/bin/dorado basecaller hac \
/home/fenrir/Documentos/Tesis/datos_gulupa/20260715_1807_MN30942_FAZ24575_f8347d8d/pod5 \
--kit-name SQK-NBD114-96 \
--min-qscore 8 \
--emit-summary \
--output-dir /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/hac_8 \
> /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/hac8_calls.bam
```

- `hac` usa el modelo de alta precisión.
- `--kit-name` nombre del kit de barcodes usado; con esto Dorado hace también el demultiplexing durante el basecalling.
- `--min-qscore` calidad mínima de las lecturas conservadas.
- `--emit-summary` genera un `sequencing_summary.txt` con métricas por lectura.
- `--output-dir` carpeta donde Dorado escribe un BAM por barcode (`barcode01/`, `barcode02/`, ...).

**Resultado esperado:** un archivo BAM por carpeta de barcode, con las lecturas ya llamadas y etiquetadas por muestra, más el resumen de secuenciación. Estas carpetas por barcode son la entrada directa del paso de Denoising y Trimming.

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
  "/home/fenrir/Documentos/Muestras_16S/sterile_sentinels" \
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

El objetivo de este paso es depurar las lecturas crudas de 16S para que solo lleguen a la clasificación taxonómica fragmentos confiables: se descartan lecturas incompletas o con una longitud que no corresponde al gen 16S, y se orientan/recortan según la posición de los primers para que todas queden en la misma dirección y contengan únicamente el amplicón de interés. Es necesario porque el basecalling por sí solo no filtra por calidad ni corrige la orientación de las lecturas.

El script [`dentrim_bam.py`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/Denoising%20y%20trimming) (entorno `dentrim_env`) automatiza este filtrado por lotes, un barcode a la vez, con dos herramientas complementarias:

## Pychopper

**Pychopper** detecta los primers en cada lectura, determina su orientación, recorta el amplicón y genera reportes de lecturas clasificadas, rechazadas o rescatadas. También filtra por un QScore mínimo. Se instala con:
```
conda install -c nanoporetech -c conda-forge -c bioconda "nanoporetech::pychopper"
```
Necesita un `FASTA` con la secuencia de los primers y un `TXT` con su orientación.

## Filtlong

**Filtlong** filtra por longitud mínima y máxima (algo que Pychopper no hace, pues solo fija un mínimo), por lo que se ejecuta antes para descartar de entrada las lecturas fuera del rango del gen 16S. Se instala con:
```
conda install bioconda::filtlong
```

## Cómo lo hace el script

Para cada barcode, `dentrim_bam.py`:
1. Convierte el BAM de Dorado a FASTQ con `samtools fastq`.
2. Orienta y recorta con `Pychopper` (backend `edlib` o `hmmer`), aplicando el umbral de calidad `-Q`.
3. Filtra por longitud y calidad media con `Filtlong` (`--min_length`/`--max_length`/`--min_mean_q`).
4. Calcula lecturas, bases, longitud media/mediana, N50 y QScore promedio antes y después del filtrado (cálculo propio a partir del FASTQ, sin dependencias externas), para cuantificar cuánto se removió en cada muestra.

**Paquetes/herramientas usados:** `samtools`, `pychopper` (requiere el módulo `edlib` en el mismo intérprete de Python), `filtlong`, y de la librería estándar `argparse`, `subprocess`, `csv`, `math`, `shutil`.

Comando:
```
python3 dentrim_bam.py \
  -i /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/hac_8 \
  -o /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/Limpieza/hac_8_trim_edlib \
  --primers /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/Limpieza/primers_gulupa/primers.fasta \
  --pconfig /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/Limpieza/primers_gulupa/primers_config.txt \
  --max-barcode 96 \
  --post-minlen 1000 --post-maxlen 1700 \
  -Q 8 -q 0.52 --fl-min-mean-q 12 \
  -m edlib -t 20 -v
```

- `-i` carpeta con las subcarpetas `barcodeXX` que salieron de Dorado.
- `-o` carpeta de salida.
- `--primers` / `--pconfig` archivos de primers que necesita Pychopper.
- `--max-barcode` número máximo de barcodes a procesar.
- `--post-minlen` / `--post-maxlen` rango de longitud aceptado (pb).
- `-Q` QScore mínimo de Pychopper; `--fl-min-mean-q` calidad media mínima para Filtlong.
- `-m` motor de alineamiento de primers (`edlib` o `hmmer`).
- `-t` número de hilos.
- `-v` imprime cada comando ejecutado.

**Resultado esperado:** un FASTQ limpio por barcode (`<barcode>_clean.fastq`) con las lecturas orientadas, recortadas y filtradas por longitud/calidad, más un `dentrim_summary.tsv` con las métricas antes/después y el porcentaje de remoción por muestra (lecturas, bases, longitud, N50, QScore).

# Taxonomía

Este paso asigna identidad taxonómica (a nivel de especie) a cada lectura 16S depurada, para saber qué bacterias componen cada muestra de suelo. Al trabajar con lecturas de longitud completa del gen (a diferencia de amplicones cortos), no es necesario agrupar en OTU/ASV: se puede clasificar directamente contra una base de referencia.

La clasificación se hace con **EMU**, que alinea las lecturas de cada muestra contra una base de referencia 16S (rrnDB + NCBI 16S RefSeq) usando **minimap2**, y luego aplica un algoritmo de **Expectation–Maximization (EM)** para resolver lecturas que mapean a varios taxones cercanos y estimar la abundancia relativa real de cada especie.

## Cómo lo hace el script

[`EMU_propio.py`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/Taxonom%C3%ADa%20y%20diversidad) orquesta EMU sobre lotes de muestras:
- Filtra de entrada los FASTQ con menos lecturas que `--min-reads-input`, para no gastar cómputo en muestras demasiado pequeñas para ser representativas.
- Ejecuta `emu abundance` por muestra en un directorio temporal y solo lo mueve a su ubicación final si terminó sin errores (ejecución atómica), para que una muestra fallida no deje resultados a medias.
- Agrega las tablas por muestra en dos tablas globales (`feature_table_relabund.tsv` y `taxonomy.tsv`), normalizando nombres de columnas taxonómicas y construyendo un ID único por taxón (`rank|tax_id|nombre`).

Como EMU entrega abundancias **relativas** y sus conteos por lectura no siempre reflejan el total real de lecturas limpias de cada muestra, [`rebuild_counts.py`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/Taxonom%C3%ADa%20y%20diversidad) reconstruye una tabla de conteos enteros confiable: multiplica cada abundancia relativa por el total real de lecturas de esa muestra (tomado del `dentrim_summary.tsv` del paso anterior, con las asignaciones de EMU o el FASTQ original como respaldo) y redondea preservando la suma total (`feature_table_counts.tsv`).

Con esa tabla de conteos, `rarefaccion.py` construye curvas de rarefacción: para cada muestra, submuestrea aleatoriamente distintas profundidades de lectura y mide cuántos taxones distintos se observan en promedio (con un paso adaptativo por muestra, proporcional a su propia profundidad, para no comparar curvas con distinta resolución). El objetivo es evaluar si la profundidad de secuenciación alcanzada fue suficiente para capturar la riqueza real del rizobioma, es decir, si la curva alcanza una meseta.

**Paquetes/herramientas usados:** `emu` (CLI, basada en `minimap2`), `pandas`, `numpy`, `matplotlib`.

## Ejecución del código
```
# 1) Clasificación taxonómica por muestra + agregación global
python3 EMU_propio.py \
  --db /home/fenrir/Documentos/Tesis/emu_db \
  --input-dir /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/Limpieza/hac_8_trim_edlib \
  --pattern "*_clean.fastq" \
  --outdir /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/EMU_propio/EMUhac_results \
  --threads 20 --rank species --min-reads-input 500 \
  --keep-counts --keep-assignments --force

# 2) Reconstrucción de la tabla de conteos enteros
python3 rebuild_counts.py --model hac --dataset results --rank species

# 3) Curva de rarefacción
python3 rarefaccion.py --model hac --dataset results --n-points 30 --iterations 10
```

- `--db` ruta a la base de datos de EMU.
- `--input-dir` / `--pattern` carpeta y patrón para localizar los FASTQ ya limpios.
- `--rank` nivel taxonómico de agregación (por defecto especie).
- `--min-reads-input` umbral mínimo de lecturas para procesar una muestra.
- `--keep-counts` / `--keep-assignments` piden a EMU que conserve conteos y asignaciones por lectura.
- `--model` / `--dataset` en `rebuild_counts.py` y `rarefaccion.py` seleccionan, por convención de carpetas, qué corrida de EMU (modelo de basecalling / subconjunto de muestras) usar.

**Resultado esperado:** `feature_table_relabund.tsv` (abundancia relativa por taxón y muestra), `feature_table_counts.tsv` (conteos enteros reconstruidos) y `taxonomy.tsv` (cadena taxonómica completa por ID de taxón), más una curva de rarefacción en PNG por muestra que indica si la profundidad de secuenciación fue suficiente (curva en meseta) o si haría falta secuenciar más.

# Diversidad

Este es el análisis que responde la pregunta central del proyecto: si la estructura y composición de la comunidad microbiana del suelo (rizobioma) difiere entre los sistemas de manejo agrícola (Empresarial, Campesino, Agroecológico). Se hace porque una tabla de abundancias por sí sola no dice si esas diferencias son reales o solo ruido de muestreo.

## Cómo lo hace el script

[`diversidad_mod.py`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/Taxonom%C3%ADa%20y%20diversidad) toma la tabla de abundancias relativas de EMU y traduce cada barcode a su sistema agrícola (barcode → ID de finca → Sistema, usando un CSV puente y un Excel de metadatos), y luego corre cuatro análisis con **scikit-bio**, **scipy** y **seaborn/matplotlib**:

- **Diversidad alfa** (estructura interna de cada muestra): índice de Simpson y riqueza observada (`skbio.diversity.alpha_diversity`), con prueba de Mann-Whitney U (Wilcoxon de dos muestras) entre cada par de sistemas agrícolas y boxplots comparativos.
- **Diversidad beta** (disimilitud entre muestras): matrices de distancia de **Bray-Curtis** (cuantitativa, sensible a cambios de abundancia) y **Jaccard** (cualitativa, presencia/ausencia), calculadas con `scipy.spatial.distance.pdist` y visualizadas como heatmaps.
- **Ordenación**: un Análisis de Coordenadas Principales (**PCoA**, `skbio.stats.ordination.pcoa`) sobre la matriz de Bray-Curtis, graficado con elipses de confianza al 95% por sistema agrícola, para ver visualmente si los grupos se solapan o se separan.
- **Contraste estadístico**: **PERMANOVA** (`skbio.stats.distance.permanova`, 999 permutaciones) para determinar si la separación observada en el PCoA es estadísticamente significativa (estadístico pseudo-F y valor p).

**Paquetes/herramientas usados:** `pandas`, `numpy`, `scikit-bio` (`alpha_diversity`, `DistanceMatrix`, `pcoa`, `permanova`), `scipy` (`pdist`, `mannwhitneyu`), `seaborn`, `matplotlib`.

## Ejecución del código
```
python3 diversidad_mod.py \
  -i /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/EMU_propio/EMUhac_results/feature_table_relabund.tsv \
  -o /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/EMU_propio/Diversidad/hac_results \
  -p Tesis_hac_results \
  -c Sistema \
  -g Empresarial Campesina Agroecológica
```

- `-i` tabla de abundancias relativas producida en Taxonomía.
- `-o` carpeta de salida; `-p` prefijo de los archivos generados.
- `-c` columna de metadatos que define los grupos a comparar (por defecto `Sistema`).
- `-g` lista de grupos a contrastar entre sí.

**Resultado esperado:** tabla de índices de diversidad alfa por muestra + p-valores de Wilcoxon y sus boxplots (PNG); matrices de distancia Bray-Curtis/Jaccard (TSV) y sus heatmaps (PNG); coordenadas y gráfico de PCoA con elipses de confianza (TSV + PNG); y un archivo de texto con el resultado de PERMANOVA (pseudo-F y valor p) que indica si los sistemas agrícolas albergan comunidades microbianas estadísticamente distintas.
