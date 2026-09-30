# Versiones de desarrollo del pipeline (octubre 2025)

Estos cuatro scripts son las **primeras versiones** del pipeline, escritas y probadas durante la fase piloto con datos públicos (ver [Anexo H](../H_fase_piloto_datos_publicos.md)). En septiembre de 2026 se retiraron de la raíz del repositorio porque las reemplazaron los scripts actuales (ver [`scripts/`](../../scripts/)). Se recuperaron sin modificaciones del historial de Git para documentar la evolución del pipeline ([Anexo I](../I_evolucion_y_decisiones.md)).

> **No usar para reproducir los resultados de la tesis.** Los resultados finales se generaron con los scripts de `scripts/`.

| Archivo | Nombre original en el repositorio | Última modificación | Commit | Reemplazado por |
|---|---|---|---|---|
| `fastq_estads.py` | `Reporte de estadísticas` | 2025-10-08 | [`1367f89`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/550577f0e764f53581d5a633ac041bf24b315e73/Reporte%20de%20estad%C3%ADsticas) | `stats_fastq.py` |
| `dentrim.py` | `Denoising y trimming` | 2025-10-08 | [`fca71f0`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/f71b85420f59eadb1ccd967b3398b3a232a871e2/Denoising%20y%20trimming) | `dentrim_bam.py` |
| `EMU.py` | `Taxonomía y diversidad` (antes `Taxonomía`) | 2025-10-21 | [`4e8ce08`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/132ffb5dd7720616b8733a0fe5b4fc35c4ea2ad3/Taxonom%C3%ADa%20y%20diversidad) | `EMU_propio.py` + `rebuild_counts.py` + `rarefaccion.py` |
| `metricas.py` | `Métricas alfa y beta` | 2025-10-22 | [`21065c6`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/24c54ec51da566f8c3c0be8fd786ee22b61af574/M%C3%A9tricas%20alfa%20y%20beta) | `diversidad_mod.py` |

Los enlaces de la columna *Commit* apuntan al archivo tal como estaba justo antes de borrarse. Solo funcionan para quien tenga acceso al repositorio mientras sea privado.

## Qué hacía cada versión y qué cambió

### `fastq_estads.py` → `stats_fastq.py`
- Calculaba por FASTQ: lecturas, bases, longitud promedio/mínima/máxima, N50, %GC, QScore promedio, **% de bases ≥Q20 y ≥Q30**, y % de lecturas con primer detectado.
- Búsqueda de primers: ventana de 120 nt en cada extremo, hasta 2 discrepancias (distancia de Hamming, las `N` del primer no penalizan), en ambas hebras, sobre las primeras 10 000 lecturas.
- Salida TXT alineado + XLSX (pandas → openpyxl → LibreOffice headless como respaldos) con columnas decimales en formato `0.00`.
- La versión actual reemplaza Q20%/Q30% por el **QScore mínimo y máximo por lectura** y agrega el filtro `--name-pattern`.

### `dentrim.py` → `dentrim_bam.py`
- Recibía FASTQ (no BAM). Ejecutaba `filtlong --min_length/--max_length | pychopper -m edlib -Q -z -Y 0 -q 0.52` (**Filtlong antes que Pychopper**).
- Guardaba solo `<muestra>_oriented_trimmed.fastq` y saltaba (y borraba) las muestras que fallaban.
- La versión actual parte del BAM de Dorado (`samtools fastq`), procesa por barcode, agrega `--fl-min-mean-q` y escribe `resumen_limpieza.tsv` con métricas antes/después.

### `EMU.py` → `EMU_propio.py` (+ `rebuild_counts.py`, `rarefaccion.py`)
- Envolvía `emu abundance` por muestra con **ejecución atómica**: todo se construía en un directorio `.__build__` y solo se movía a `--outdir` si todas las muestras terminaban bien.
- Agregaba `feature_table_counts.tsv`, `feature_table_relabund.tsv`, `taxonomy.tsv` (formato QIIME `k__;p__;…;s__`) y `manifest.tsv`.
- Calculaba además diversidad alfa (Observed, Chao1, Shannon, Simpson) y, con `--do-pcoa`, Bray-Curtis/Jaccard + PCoA clásico implementado con NumPy.
- En la versión actual, EMU y diversidad quedan en scripts separados, las tablas se renombran al español y la reconstrucción de conteos pasa a `rebuild_counts.py`.

### `metricas.py` → `diversidad_mod.py`
- Leía las tablas de `EMU.py` y calculaba alfa (Observed, Chao1, Shannon, Simpson), Bray-Curtis, Jaccard y PCoA **con implementaciones propias** (NumPy).
- Figuras estáticas (matplotlib) e interactivas (plotly): PCoA y barras apiladas de los N géneros más abundantes.
- No tenía pruebas estadísticas. La versión actual usa **scikit-bio** (`alpha_diversity`, `pcoa`, `permanova`) y **scipy** (`mannwhitneyu`), y agrega PERMANOVA, elipses de confianza al 95 % y comparación entre sistemas agrícolas.

## Otros fragmentos que solo existieron dentro del README

Durante octubre de 2025 el README también contenía código embebido que ya no está en ningún archivo:

- Un **bucle Bash** `filtlong | pychopper` que guardaba todos los reportes de Pychopper (`-r -S -A -K -l -u -w`). Commit [`fca71f0`](https://github.com/xXFenrir/16S-FenrirPipeline/blob/fca71f0807e393a35cd3aa29c4f47425061f7732/README.md).
- **`metricasd.py`** ("Informe Microbiota (sin R)"): alfa, beta (incluida la distancia de Hellinger) y un reporte HTML con gráficos SVG hechos a mano, escrito en Python puro por **conflictos con paquetes de R/Bioconductor**. Mismo commit.
- Un script de **grafos taxonómicos** con `networkx` y `graphviz` (jerarquía superkingdom → species). Mismo commit.

Para verlos: `git show fca71f0:README.md`.
