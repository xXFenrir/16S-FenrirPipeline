# Anexo H. Fase piloto: construcción del pipeline con datos públicos

Antes de tener los datos de gulupa (corrida del 15 de julio de 2026), el pipeline se construyó y probó en octubre de 2025 con lecturas públicas de un estudio de suelos agrícolas secuenciado con MinION. Este anexo documenta ese conjunto de datos, cómo lo procesó el estudio original y qué se tomó de él.

## H.1 Conjunto de datos

| Dato | Valor |
|---|---|
| Estudio | Erlandson, S. R., Ewing, P. M., Osborne, S. L., y Lehman, R. M. (2024). Sterile sentinels and MinION sequencing capture active soil microbial communities that differentiate crop rotations. *Environmental Microbiome*. https://doi.org/10.1186/s40793-024-00571-8 |
| Tema | Comunidades microbianas activas del suelo bajo distintas rotaciones de cultivo, capturadas con bolsas de suelo estéril enterradas (*sterile sentinels*) |
| Secuenciación | Oxford Nanopore MinION, gen 16S rRNA completo (primers 27F/1492R) |
| Datos crudos | NCBI SRA, BioProject PRJNA1020132 (según la declaración de disponibilidad de datos del artículo) |
| Código del estudio | https://github.com/serlandson/sterile_sentinels |
| Corridas usadas en los comandos documentados | SRR26147154, SRR26147165 (y el resto de la carpeta `sterile_sentinels`) ⏳ confirmar que pertenecen a PRJNA1020132 |
| Carpeta local | `/home/fenrir/Documentos/Muestras_16S/sterile_sentinels` (antes `Muestras 16S`, con espacio) |

Se eligió porque comparte el problema de la tesis: suelo agrícola, efecto del manejo sobre la comunidad microbiana y 16S de longitud completa con Nanopore.

## H.2 Qué se probó en la fase piloto

Los FASTQ públicos ya venían con basecalling y demultiplexing, así que la fase piloto cubrió de las estadísticas en adelante:

| Etapa | Script de desarrollo | Parámetros probados |
|---|---|---|
| Estadísticas de los FASTQ originales | `fastq_estads.py` | Primers `primers_stesen.fasta`, búsqueda recursiva |
| Limpieza | Bucle Bash y luego `dentrim.py` (Filtlong → Pychopper) | Longitud 1300–1700 pb y Pychopper Q9; luego 1000–1700 pb y Q12; `-q 0.52`, `-Y 0` |
| Taxonomía | `EMU.py` | EMU con la base por defecto, nivel especie, 8 hilos |
| Diversidad | `metricasd.py` (sin R), luego `metricas.py` | Observed, Chao1, Shannon, Simpson; Bray-Curtis, Jaccard, Hellinger; PCoA |

Los scripts están en [`historico/`](historico/). Resultados de la fase piloto copiados por el script de recolección: `anexos/datos/piloto/` ⏳.

## H.3 Primers y procesamiento del estudio original

Fuente: archivo `bacteria minion read processing.txt` del repositorio del estudio.

Primers con cola de adaptador. Se separa la cola (22 nt) de la región que hibrida con el 16S:

| Nombre | Secuencia completa (5'→3') | Región 16S |
|---|---|---|
| Bac27_F | `TTTCTGTTGGTGCTGATATTGC` + `AGRGTTYGATYMTGGCTCAG` | 27F |
| Univ_1492 | `ACTTGCCTGTCGCTCTATCTTC` + `TACCTTGTTACGACTT` | 1492R |

Procesamiento del estudio original:

| Etapa | Herramienta y parámetros |
|---|---|
| Basecalling | Guppy, modelo `dna_r9.4.1_450bps_hac` (celdas R9.4.1) |
| Demultiplexing | `guppy_barcoder --barcode_kits EXP-PBC096 --require_barcodes_both_ends --trim_barcodes` |
| Filtro de calidad y longitud | `NanoFilt --quality 12 --length 1000 --maxlength 1700` |
| Recorte de primers | `cutadapt -g AGRGTTYGATYMTGGCTCAG -a TACCTTGTTACGACTT --rc -e 0.15 -m 1000` |
| Taxonomía | `emu abundance --type sr --db emu/default_16S --N 50 --keep-counts` |
| Análisis estadístico | Scripts en R |

## H.4 Comparación con el pipeline de la tesis

| Etapa | Estudio de referencia (2024) | Pipeline final de la tesis (2026) |
|---|---|---|
| Química / celda | R9.4.1 | Kit V14 (SQK-NBD114-96) |
| Basecalling | Guppy HAC | Dorado 2.1.0 HAC |
| Demultiplexing | guppy_barcoder, paso aparte | Dorado, durante el basecalling |
| Rango de longitud | 1000–1700 pb (NanoFilt) | 1000–1700 pb (Filtlong) |
| Calidad mínima | Q12 (NanoFilt) | Q8 (Dorado y Pychopper); Filtlong: ver [verificación 1](C_protocolo_comandos_parametros.md#verificaciones-pendientes-de-la-etapa-de-limpieza) |
| Primers | cutadapt, hasta 15 % de error | Pychopper con edlib, `-q 0.52`, y reorientación de lecturas |
| Taxonomía | EMU, base por defecto, `--type sr` | EMU, base por defecto, nivel especie |
| Estadística | R | Python (scikit-bio, SciPy) |

El rango de longitud final (1000–1700 pb) coincide con el del estudio de referencia.
