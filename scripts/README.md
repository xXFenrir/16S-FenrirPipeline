# Scripts del pipeline (versión final)

Esta carpeta contiene los scripts con los que se generaron los resultados de la tesis. Se copian desde el PC con [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh).

| Orden | Script | Etapa | Entrada principal | Salida principal |
|---|---|---|---|---|
| 1 | `stats_fastq.py` | Estadísticas básicas de FASTQ | FASTQ/FASTQ.GZ (crudos o limpios) | Tabla TXT + XLSX por muestra |
| 2 | `dentrim_bam.py` | Denoising y trimming | Carpetas `barcodeXX/` con BAM de Dorado | `<barcode>_limpio.fastq` + `resumen_limpieza.tsv` |
| 3 | `EMU_propio.py` | Clasificación taxonómica (EMU) | `*_limpio.fastq` | `tabla_abundancia_relativa.tsv`, `taxonomia.tsv` |
| 4 | `rebuild_counts.py` | Reconstrucción de conteos enteros | Abundancias relativas + `resumen_limpieza.tsv` | `tabla_conteos.tsv` |
| 5 | `rarefaccion.py` | Curvas de rarefacción | `tabla_conteos.tsv` | Una curva PNG por muestra |
| 6 | `diversidad_mod.py` | Diversidad alfa/beta, PCoA, PERMANOVA | `tabla_abundancia_relativa.tsv` + metadatos | Tablas TSV, figuras PNG, resultado de PERMANOVA |

`stats_fastq.py` puede correrse antes y después de la limpieza para comparar.

Los comandos completos y la explicación de cada parámetro están en el [Anexo C](../anexos/C_protocolo_comandos_parametros.md). Las versiones de desarrollo de 2025 están en [`anexos/historico/`](../anexos/historico/).
