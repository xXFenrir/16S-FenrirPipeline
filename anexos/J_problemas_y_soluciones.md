# Anexo J. Problemas encontrados y soluciones

Registro de los errores y obstáculos técnicos que aparecieron al construir el pipeline y cómo se resolvieron. Sirve para justificar decisiones de diseño de los scripts y como guía para quien repita el análisis.

## J.1 Dorado ya no acepta FAST5

**Problema.** Las versiones actuales de Dorado solo leen POD5; los datos en FAST5 no se pueden procesar directamente.

**Solución.** Convertir antes con la herramienta de ONT:

```bash
pod5 convert fast5 /ruta/fast5/*.fast5 --output pod5_out/
```

La corrida de gulupa ya se entregó en POD5, así que no hizo falta convertirla.

## J.2 Pychopper necesita `edlib` en el mismo intérprete de Python

**Problema.** Con `-m edlib`, Pychopper importa el módulo `edlib` de Python; si el módulo está en otro entorno, falla.

**Solución.** Crear un entorno propio para la limpieza (`dentrim_env`) con samtools, Pychopper, edlib y Filtlong, y ejecutar `dentrim_bam.py` desde él.

## J.3 Pychopper no acepta la configuración de primers escrita en la línea de comandos

**Problema** (prueba del 2026-07-23 sobre `barcode11` de la corrida SUP). Se pasó la configuración directamente después de `-c`:

```text
(dentrim_env) fenrir@fenrir-Vector-16-HX-AI-A2XWHG:~/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/Limpieza$ pychopper -m edlib \
  -b /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/Limpieza/primers_gulupa/primers.fasta \
  -c "+:16sF,-16sR|+:16sR,-16sF" \
  -Q 8 -z 1000 -t 20 \
  /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/Limpieza/sup_8_limp/barcode11/barcode11_raw.fastq \
  /home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/Limpieza/sup_8_limp/barcode11/test_oriented.fastq
Traceback (most recent call last):
  File "/home/fenrir/miniconda3/envs/dentrim_env/bin/pychopper", line 10, in <module>
    sys.exit(main())
  File "/home/fenrir/miniconda3/envs/dentrim_env/lib/python3.9/site-packages/pychopper/scripts/pychopper.py", line 318, in main
    CONFIG = open(args.c, "r").readline().strip()
FileNotFoundError: [Errno 2] No such file or directory: '+:16sF,-16sR|+:16sR,-16sF'
```

**Causa.** `-c` espera la **ruta de un archivo**; Pychopper abre el archivo y lee su primera línea (`pychopper.py`, línea 318).

**Solución.** Guardar la configuración en `primers_config.txt` y pasar esa ruta (`--pconfig` en `dentrim_bam.py`).

**Pendiente.** En esa configuración las dos entradas empiezan por `+:`. El formato de Pychopper marca la segunda con `-:` (`+:16sF,-16sR|-:16sR,-16sF`) para que las lecturas reversas se reorienten. Revisar el contenido final de `primers_config.txt` ([verificación 2 del Anexo C](C_protocolo_comandos_parametros.md#verificaciones-pendientes-de-la-etapa-de-limpieza)).

## J.4 Rutas con espacios

**Problema.** La carpeta de la fase piloto se llamaba `Muestras 16S` (con espacio), lo que rompe comandos armados como texto.

**Solución.** Se renombró a `Muestras_16S` (commits `9cb15b1` y `0a92303`, 2025-10-20). Además, `dentrim.py` cita cada ruta con `shlex.quote()` antes de armar la tubería de Filtlong y Pychopper.

## J.5 Muestras que se quedan sin lecturas en la limpieza

**Problema.** Si una muestra no conserva lecturas tras el filtrado, Pychopper termina con error y puede quedar un FASTQ parcial.

**Solución.** Desde `dentrim.py` se ejecuta la tubería con `set -o pipefail`; si el código de salida no es 0, se borra la salida parcial, se marca la muestra como saltada y se continúa con la siguiente. Al final se listan las muestras saltadas.

## J.6 Resultados a medias cuando falla una muestra en EMU

**Problema.** Si EMU falla a mitad de un lote, las tablas globales pueden mezclar muestras completas con incompletas.

**Solución.** Ejecución atómica: cada muestra se procesa en un directorio temporal y solo se mueve a su ubicación final si terminó sin errores. En la versión de desarrollo todo el lote se construía en un directorio `.__build__` y se movía a la salida final únicamente si todo salía bien.

## J.7 Conteos de EMU que no suman el total de lecturas

**Problema.** Las abundancias de EMU son relativas y sus conteos estimados no siempre coinciden con el total de lecturas limpias de la muestra, lo que afecta la rarefacción.

**Solución.** `rebuild_counts.py` reconstruye los conteos multiplicando la abundancia relativa por el total real de lecturas limpias (de `resumen_limpieza.tsv`) y redondeando sin cambiar el total.

## J.8 Conflictos de paquetes de R y Bioconductor

**Problema.** Al preparar el análisis de diversidad aparecieron conflictos con paquetes de R y Bioconductor (octubre de 2025).

**Solución.** Hacer toda la estadística en Python: primero con implementaciones propias (`metricasd.py`, `metricas.py`) y en la versión final con scikit-bio y SciPy.

## J.9 Varios entornos conda y el ejecutable de EMU

**Problema.** Con varios entornos conda, `emu` podía no estar en el `PATH` o ser otra instalación.

**Solución.** El envoltorio de EMU de desarrollo aceptaba `--emu-cmd` con la ruta exacta del ejecutable (por ejemplo `/home/fenrir/anaconda3/envs/pipelinefenrir/bin/emu`) y verificaba `emu --version` antes de empezar.

## J.10 Generar el XLSX sin depender de una sola librería

**Problema.** No siempre estaban instalados `pandas` u `openpyxl` para escribir Excel.

**Solución.** El script de estadísticas intenta en orden `pandas`, `openpyxl` y LibreOffice en modo headless; si ninguno funciona, deja un TSV e indica cómo convertirlo.
