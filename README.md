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
  -i "/home/fenrir/Documentos/Muestras_16S/sterile_sentinels" \
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
Una de la herramientas más utilizadas según la bibliografía revisada para la clasificación taxonómica, es EMU. Pues es una herramienta basada en minimap 2 (otro clasificador taxonómico), pero que ya integra las bases de datos de referencia y permite extraer la abundancia relativa de las muestras. El objetivo de este paso es obtener tablas conteos, abundancias relativas, asignación de lecturas y mapeo taxonómico, que esten listas para análisis de diversidad alfa y beta. Por lo que, el script integra un módulo opcional para calcular Bray–Curtis/Jaccard y su PCoA directamente desde la tabla de conteos.

## Descripción del código
El script para la [Taxonomía](https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/Taxonom%C3%ADa) recibe como entrada los archivos FASTQ filtrados. Donde ejecuta `emu abundance` por muestra usando la base de datos que la documentación proporciona ([rrnDB v5.6](https://rrndb.umms.med.umich.edu/?_x_tr_sl=en&_x_tr_tl=fr&_x_tr_hl=fr&_x_tr_pto=sc) y [NCBI 16S RefSeq](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/)). Además, guarda abundancias relativas, conteos, lecturas asignadas, genera un mapeo de la clasificación taxonómica, y calcula métricas de diversidad alfa (Observed, Chao1, Shannon, Simpson). Usando `--do-pcoa`, construye distancias Bray–Curtis y Jaccard, realiza PCoA y escribe coordenadas, varianza explicada y figuras.

## Paquetes importados

Los paquetes importados fueron:
- `argparse ()`: Permite construir la interfaz de línea de comandos definiendo las entradas y variables.
- `os ()`: Verifica la existencia de los archivos de entrada, crea carpetas de salida por muestra, compone rutas portables y elimina salidas parciales si la ejecución falla.
- `sys ()`: Proporciona control sobre la salida y la terminación del programa por medio de avisos.
- `glob ()`: Sirve para encontrar archivos con sufijo `FASTQ` independientemente de si está o no comprimido. En caso de no encontrar nada, el programa avisa y termina.
- `subprocess ()`: Permite que se ejecuten subprocesos, como dirigir cada salida de EMU (`emu abundance`, etc), captura fallos y continúa con otras muestras sin cortar el ciclo.
- `pathlib.Path ()`: Obtiene el nombre base del archivo quitando sus sufijos.
- `typing ()`: Permite que se informe de errores antes de correr el script.
- `numpy`: Usa las métricas para calcular alfa diversidad (Shannon, Simpson, Chao1), normalizaciones y operaciones rápidas sobre matrices.
- `pandas`: Esencial para las tablas, pues lee los TSV de EMU, armoniza columnas taxonómicas y construye matrices de conteos.
- `biom ()`: Exporta la tabla de conteos a BIOM para interoperar con phyloseq. (Aún en revisión)
- `scikit-bio ()`: Obtiene la beta-diversidad con Bray–Curtis y Jaccard, y realiza la PCoA, escribiendo las distancias y coordenadas. (Aún en revisión)
- `matplotlib.pyplot ()`: Genera las figuras de los PCoA. (Aún en revisión)
- `re ()`: Reconoce una expresión regular para detectar espécíficamente el archivo deseado en un grupo de datos.
- `shutil.which ()`: Localiza la ruta de EMU.

## Funciones definidas

Se definieron las siguientes funciones:
- `eprint()`: Imprime mensajes en la salida en caso de error del procesamiento de una muestra o de todo el script.
- `ensure_dir(path)`: Crea el directorio indicado si no existe.
- `which(cmd)`: Permite validar y ubicar el binario de EMU, para evitar conflictos en el entorno.
- `safe_basename_noext(path)`: Garantiza nombres coherentes para subcarpetas y prefijos de salida, para conservar el ID de la muestra.
- `read_table_maybe(path)`: Intenta leer un TSV con `pandas` y devuelve `None` si falla, registrando una advertencia.
- `to_qiime_tax_string(row)`: Construye una cadena taxonómica separado por columnas. 
- `chao1(counts)`: Calcula el estimador de riqueza Chao1 usando singletons y doubletons.
- `shannon_entropy(p)`: Calcula el índice de Shannon sobre proporciones (p), ignorando ceros.
- `simpson_index(p)`: Calcula el índice de Simpson como `(1-\sum p^2)`.
- `run_emu_for_sample()`: Llama a `emu abundance` por muestra, verifica el retorno y envía las salidas (`rel-abundance.tsv`, `counts.tsv` y `read-assignments.tsv`) a la ruta solicitada.
- `select_rank()`: Filtra la tabla de EMU al nivel taxonómico elegido.
- `build_feature_id(row, rank)`: Genera un identificador por taxón con el patrón `rank|tax_id|tax_name`, evitando confuciones entre nombres iguales con distintos `tax_id`.
- `aggregate_tables()`: Pricipal encargado de generar tablas, pues organiza columnas taxonómicas, aplica `select_rank` y el umbral `min_abundance`, construye matrices multi-muestra de conteos y relativas, emite `taxonomy.tsv` y calcula alfa (Observed/Chao1/Shannon/Simpson).
- `compute_beta_pcoa(counts_tsv, outdir, use_relative)`: A partir de la tabla de conteos, calcula distancias Bray–Curtis y Jaccard, ejecuta PCoA y escribe matrices de distancia, coordenadas, eigenvalores, varianza explicada y, si puede, figuras.
- `_write_pcoa_outputs(ord_res, prefix, outdir)`: Estandariza las salidas de PCoA en TSV. 
- `_plot_pcoa(ord_res, prefix, outdir)`: Genera un scatter PC1 vs PC2 con las proporciones de varianza en etiquetas y anota las muestras.
- `discover_fastqs(args)`: Localiza archivos de entrada usando `--input-glob` y `--input-dir`+`--pattern`, desduplica y ordena rutas absolutas. 
- `merge_assignments(assign_paths, outdir)`: Une las tablas `read-assignments` de todas las muestras en un CSV con columna `sample`.
- `main()`: Parsea argumentos, valida el entorno, descubre FASTQ, corre EMU por muestra, construye el `manifest.tsv`, agrega resultados vía `aggregate_tables`, fusiona asignaciones y ejecuta `compute_beta_pcoa`, y reporta un resumen final.

## Ejecución del código
```
python /home/fenrir/scriptsbioinf/EMU.py \
  --db /home/fenrir/emu_db \
  --input-dir /home/fenrir/results_dentrim_Q12 \
  --pattern '*_clean.fastq' \
  --outdir /home/fenrir/results_dentrim_Q12/EMU_Q12 \
  --threads 8 \
  --keep-counts --keep-assignments \
  --emu-cmd /home/fenrir/anaconda3/envs/pipelinefenrir/bin/emu
```
# Diversidad
```

```
