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
