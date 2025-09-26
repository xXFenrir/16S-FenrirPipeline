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

Desde [GutHub](https://github.com/epi2me-labs/pychopper):
```
conda install -c nanoporetech -c conda-forge -c bioconda "nanoporetech::pychopper"
```

Desde [Anaconda](https://anaconda.org/bioconda/pychopper):
```
conda install bioconda::pychopper 
```
Posteriormente, se debe revisar el comando base para usar `Chopper`, este puede verse con:
```
chopper -h
```
Mostrando así la composición básica del comando.
```
usage: pychopper [-h] [-b primers] [-g phmm_file] [-c config_file] [-k {PCS109,PCS110,PCS111,PCS114,LSK114,PCB111,PCB114}] [-q cutoff] [-Q min_qual] [-z min_len] [-r report_pdf]
                 [-u unclass_output] [-l len_fail_output] [-w rescue_output] [-S stats_output] [-K qc_fail_output] [-Y autotune_nr] [-L autotune_samples] [-A scores_output]
                 [-m method] [-x rescue] [-p] [-t threads] [-B batch_size] [-D read stats] [-y] [-U]
                 input_fastx [output_fastx]
```
En este caso se debe crear una carpeta y su respectiva ruta, para que allí sean descargados los resultados después de la limpieza de la muestra. Adicionalmente, que se tenga una archivo `.fasta` con los primers y otro `.txt` con la orientación de los mismos. De esta manera, ya es posible hacer uso del comando base de `Chopper`.
```
mkdir -p resultados

pychopper -m edlib \
  -b resources/primers_stesen.fasta \ #Identifica la ruta y el archivo donde está la secuencia de los primers
  -c resources/primer_stesen.txt \ #Identica la ruta y el archivo con la orientación de los mismos
  -Q 9 -z 1200 -t 4 \ #`-Q 9` indica que es QScore mínimo es 9, `-z 1200` que la longitud mínima es de 1200pb, `-t` indica la cantidad de núcleos de la CPU que usará `PyChopper` en simultáneo.
  -r results/sample1_report.pdf \ #Reporte gráfico de los resultados en formato PDF
  -u results/sample1_unclassified.fastq \ #fastq de las lecturas que se excluyeron
  -w results/sample1_rescued.fastq \ #fastq de las lecturas recuperadas
  "/home/fenrir/Documentos/Muestras 16S/sterile_sentinels/SRR26147154.fastq.gz" \ #Dirección de la muestra a tratar
  "results/sample1_oriented_trimmed.fastq" 
```
