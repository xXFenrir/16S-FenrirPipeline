# Anexo E. Control de calidad de las lecturas

Métricas de calidad por muestra en cada punto del pipeline: después del basecalling, antes y después de la limpieza, y de las lecturas limpias que entran a la clasificación taxonómica. En el documento suele bastar con el resumen por sistema agrícola; aquí van las tablas completas.

Las tablas se completan con los archivos que copia [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh) en `anexos/datos/calidad/` y `anexos/datos/corrida/`.

## E.1 Definición de las métricas

| Métrica | Definición |
|---|---|
| Lecturas | Número de registros del FASTQ. |
| Bases | Suma de las longitudes de todas las lecturas. |
| Longitud promedio / mediana / mínima / máxima | Estadísticos de la distribución de longitudes de lectura (pb). |
| N50 | Longitud *L* tal que las lecturas de longitud ≥ *L* suman al menos el 50 % de las bases. |
| %GC | (G + C) / bases × 100. |
| QScore promedio | Calidad media en escala Phred. Ver la nota sobre cómo se promedia. |
| QScore mínimo / máximo por lectura | El menor y el mayor QScore medio entre las lecturas de la muestra (`stats_fastq.py`). |
| % primers | Porcentaje de lecturas con algún primer en los extremos. En la versión de desarrollo: ventana de 120 nt en cada extremo, hasta 2 discrepancias, ambas hebras, primeras 10 000 lecturas. |
| % removido | (valor antes − valor después) / valor antes × 100, para lecturas y bases (`resumen_limpieza.tsv`). |

**Nota sobre el QScore promedio.** Hay dos maneras de promediar calidades y dan números distintos:

1. **Promedio de probabilidades de error**, convertido de nuevo a Phred: $\bar{Q} = -10 \log_{10}\left(\frac{1}{n}\sum_i 10^{-Q_i/10}\right)$. Es la que usan Dorado (`--min-qscore`) y Pychopper (`-Q`) para filtrar.
2. **Promedio aritmético de los valores Phred**: $\bar{Q} = \frac{1}{n}\sum_i Q_i$. Es la que usaba la versión de desarrollo `fastq_estads.py` (suma de Phred / bases). Da valores más altos que la primera.

Conviene indicar en el documento cuál usan `stats_fastq.py` y `dentrim_bam.py` en su versión final, para que los QScore de las tablas se puedan comparar con los umbrales de filtrado (Q8).

## E.2 Lecturas por barcode tras el basecalling

Resumen de `sequencing_summary.txt` (Dorado, `--emit-summary`) hecho por el script de recolección: lecturas, bases, longitud media y QScore medio por barcode, incluidas las lecturas sin clasificar.

Archivo: `anexos/datos/corrida/lecturas_por_barcode_hac.tsv` ⏳

| Barcode | Lecturas | Bases | Longitud media | QScore medio |
|---|---|---|---|---|
| ⏳ | | | | |
| unclassified | ⏳ | | | |

## E.3 Antes y después de la limpieza

Fuente: `resumen_limpieza.tsv` de `dentrim_bam.py` → `anexos/datos/calidad/limpieza_hac/resumen_limpieza.tsv` ⏳

| Barcode | Lecturas antes | Lecturas después | % removido | Bases antes | Bases después | Longitud media antes → después | N50 antes → después | QScore antes → después |
|---|---|---|---|---|---|---|---|---|
| ⏳ | | | | | | | | |

Para el documento conviene resumir esta tabla por sistema agrícola (media ± desviación estándar de lecturas limpias y % removido).

## E.4 Lecturas limpias

Fuente: `estadisticas_hac.txt` / `.xlsx` de `stats_fastq.py` → `anexos/datos/calidad/` ⏳

| Muestra | Lecturas | Bases | Long. promedio | Long. mín. | Long. máx. | N50 | %GC | QScore promedio | QScore mín./máx. por lectura | % primers |
|---|---|---|---|---|---|---|---|---|---|---|
| ⏳ | | | | | | | | | | |

## E.5 Comparación HAC y SUP

Si en la tesis se compara el modelo HAC con SUP ([Anexo I](I_evolucion_y_decisiones.md#i4-modelos-de-basecalling-hac-y-sup)), el script copia también los resultados de la limpieza de la corrida SUP (`Limpieza/sup_8_trim_edlib/`) en `anexos/datos/calidad/limpieza_sup/`.

| Métrica (mediana por muestra) | HAC | SUP |
|---|---|---|
| Lecturas limpias | ⏳ | ⏳ |
| % removido | ⏳ | ⏳ |
| QScore promedio | ⏳ | ⏳ |
| Especies detectadas | ⏳ | ⏳ |

Los parámetros de limpieza de ambas corridas no fueron iguales (ver [Anexo I](I_evolucion_y_decisiones.md#i3-evolución-de-los-parámetros)). Si se comparan, hay que decirlo o repetir la limpieza con los mismos parámetros.
