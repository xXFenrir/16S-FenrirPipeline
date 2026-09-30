# Anexo I. Evolución del pipeline y decisiones metodológicas

Cómo cambió el pipeline entre octubre de 2025 y septiembre de 2026, qué alternativas se probaron y por qué se eligieron las herramientas y parámetros finales. Todo lo que aparece aquí está respaldado por el historial de Git del repositorio; lo que falta confirmar está marcado con ⏳.

## I.1 Cronología

| Fecha | Hito | Commit |
|---|---|---|
| 2025-08-13 | Fork de Dorado en GitHub (`xXFenrir/dorado`), sin cambios propios | — |
| 2025-10-08 | Primera versión del README: basecalling, demultiplexing, estadísticas y limpieza; scripts `Reporte de estadísticas` y `Denoising y trimming`. La limpieza pasa de un bucle Bash a `dentrim.py` | `fca71f0` |
| 2025-10-18 | Envoltorio de EMU con agregación de tablas y diversidad alfa; el mismo día se agregan Bray-Curtis, Jaccard y PCoA | `cfa98c0`, `ad8cb61` |
| 2025-10-21 | Ejecución atómica en el envoltorio de EMU; pruebas con lecturas limpias a Q12 | `4e8ce08`, `4675d5d` |
| 2025-10-21 → 22 | Script de métricas separado (`metricas.py`) con figuras estáticas e interactivas | `d9c90df`, `21065c6` |
| 2026-07-15 | Corrida de MinION de las muestras de gulupa | — |
| 2026-07-21 | Dorado 1.1.1 → 2.1.0; basecalling HAC sobre POD5 | `0bb6247`, `e37cb69` |
| 2026-07-23 | Basecalling con demultiplexing integrado (`--kit-name`) y `--min-qscore 8`; limpieza de la corrida SUP con `dentrim_bam.py` | `c40cbe8`, `f189e07` |
| 2026-07-27 | Curvas de rarefacción sobre la corrida SUP | `6d3597b` |
| 2026-09-21 | README con los scripts finales; se retiran los scripts de desarrollo | `c9e7bff` → `c520392` |

## I.2 Decisiones metodológicas

| Decisión | Alternativas consideradas | Motivo | Fuente |
|---|---|---|---|
| Basecalling y demultiplexing con Dorado | Guppy; Porechop para demultiplexing | Dorado es la herramienta recomendada por ONT y sucesora de Guppy | README |
| Demultiplexing durante el basecalling (`--kit-name`) | `dorado demux --no-classify` en un paso aparte, o en tubería con `\|` | Recomendación del fabricante; un solo comando | README (jul. y sep. 2026) |
| Recorte de primers con Pychopper | `dorado basecaller --trim all --primer-sequences` (probado el 2026-07-21) | ⏳ Según el README, Pychopper además orienta las lecturas y reporta las clasificadas, rescatadas y rechazadas | README |
| Filtlong además de Pychopper | Solo Pychopper | Pychopper solo fija una longitud mínima; Filtlong permite un rango | README |
| Backend `edlib` en Pychopper | `phmm` (perfiles HMM) | ⏳ | Nombre de la carpeta `hac_8_trim_edlib` |
| Clasificación directa con EMU, sin OTU ni ASV | Agrupar en OTU/ASV | Con lecturas del gen completo se puede clasificar contra una base de referencia sin agrupar | README |
| Reconstruir conteos (`rebuild_counts.py`) | Usar los conteos de EMU | Los conteos de EMU no siempre suman el total real de lecturas limpias | README |
| Mínimo de 500 lecturas por muestra | Sin mínimo | Evitar muestras demasiado pequeñas para ser representativas | README |
| Estadística en Python | R / Bioconductor (phyloseq) | Conflictos de paquetes de R y Bioconductor | README (oct. 2025) |
| scikit-bio para alfa, PCoA y PERMANOVA | Implementaciones propias con NumPy (versiones de desarrollo) | ⏳ La versión final agrega PERMANOVA y pruebas entre sistemas, que las de desarrollo no tenían | Comparación de versiones |
| Modelo HAC para los resultados finales | SUP | ⏳ | Ver I.4 |
| Filogenia | Árbol + UniFrac | Aparece como paso opcional en el README; no se implementó | README |

## I.3 Evolución de los parámetros

| Parámetro | Fase piloto (oct. 2025) | Prueba SUP (jul. 2026) | Final HAC (sep. 2026) |
|---|---|---|---|
| Dorado | 1.1.1 (FASTQ públicos, sin basecalling propio) | — (salida de MinKNOW) | 2.1.0, `hac` |
| Kit | SQK-16S114-24 (ejemplo en el README) | SQK-NBD114-96 | SQK-NBD114-96 |
| QScore mínimo en basecalling | — | ⏳ | 8 |
| Orden de la limpieza | Filtlong → Pychopper | ⏳ | Pychopper → Filtlong |
| Rango de longitud (pb) | 1300–1700, luego 1000–1700 | 1200–1600 | 1000–1700 |
| Pychopper `-Q` | 9, luego 12 | 8 | 8 |
| Pychopper `-q` | 0.52 | 0.3 | 0.52 |
| Filtlong `--min_mean_q` | — | — | 12 (ver [verificación 1](C_protocolo_comandos_parametros.md#verificaciones-pendientes-de-la-etapa-de-limpieza)) |
| Barcodes procesados | — | hasta 73 | hasta 96 |
| Mínimo de lecturas por muestra | — | 500 (rarefacción) | 500 (EMU) |
| Rarefacción | — | Paso fijo de 200 lecturas | 30 puntos, paso proporcional a la profundidad, 10 iteraciones |
| Diversidad alfa | Observed, Chao1, Shannon, Simpson | — | Simpson, riqueza observada |
| Diversidad beta | Bray-Curtis, Jaccard (+ Hellinger) | — | Bray-Curtis, Jaccard |
| Pruebas estadísticas | Ninguna | — | Mann-Whitney U, PERMANOVA (999 permutaciones) |

## I.4 Modelos de basecalling: HAC y SUP

La corrida se procesó con dos modelos de basecalling:

| | SUP | HAC |
|---|---|---|
| Origen | Carpeta de MinKNOW `sup_8/Microbioma_15072026/.../bam_pass` | Dorado 2.1.0 sobre los POD5 (`hac_8/`) |
| Limpieza | `Limpieza/sup_8_trim_edlib` (1200–1600 pb, `-q 0.3`) | `Limpieza/hac_8_trim_edlib` (1000–1700 pb, `-q 0.52`) |
| Taxonomía | `EMU_propio/EMUsup_results` | `EMU_propio/EMUhac_results` |
| Uso en la tesis | ⏳ | Resultados finales |

SUP (*super accuracy*) es más exacto por base que HAC pero mucho más lento. `rebuild_counts.py` y `rarefaccion.py` aceptan `--model hac` o `--model sup`. Motivo para usar HAC en los resultados finales: ⏳.

Si en el documento se comparan ambos modelos, hay que tener en cuenta que la limpieza no usó los mismos parámetros (tabla de arriba).

## I.5 Versiones de los scripts

| Versión de desarrollo (2025) | Versión final (2026) |
|---|---|
| `fastq_estads.py` | `stats_fastq.py` |
| `dentrim.py` | `dentrim_bam.py` |
| `EMU.py` | `EMU_propio.py`, `rebuild_counts.py`, `rarefaccion.py` |
| `metricas.py` (y `metricasd.py`, solo en el README) | `diversidad_mod.py` |

Qué cambió en cada uno: [`historico/README.md`](historico/README.md).
