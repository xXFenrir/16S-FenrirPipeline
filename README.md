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
En este caso se debe crear una carpeta y su respectiva ruta, para que allí sean descargados los resultados después de la limpieza de la muestra. Adicionalmente, que se tenga una archivo `.fasta` con los primers y otro `.txt` con la orientación de los mismos. De esta manera, ya es posible hacer uso del comando base de `Pychopper`.
```
mkdir -p resultados

pychopper -m edlib \
  -b resources/primers_stesen.fasta \ #Identifica la ruta y el archivo donde está la secuencia de los primers
  -c resources/primer_stesen.txt \ #Identica la ruta y el archivo con la orientación de los mismos
  -Q 9 -z 1300 -t 4 \ #`-Q 9` indica que es QScore mínimo es 9, `-z 1200` que la longitud mínima es de 1200pb, `-t` indica la cantidad de núcleos de la CPU que usará `PyChopper` en simultáneo.
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

# Taxonomía
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
