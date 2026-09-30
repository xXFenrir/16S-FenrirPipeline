# Entornos conda

Definiciones exportadas de los entornos conda usados en el pipeline. Se generan en el PC con [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh), que exporta dos archivos por entorno:

| Archivo | Contenido | Uso |
|---|---|---|
| `<entorno>.yml` | Todos los paquetes con versión exacta (`conda env export --no-builds`) | Reproducir el entorno lo más fielmente posible en Linux |
| `<entorno>_historial.yml` | Solo los paquetes que se instalaron explícitamente (`conda env export --from-history`) | Recrear el entorno en otro sistema operativo o cuando el primero falla |

Para recrear un entorno:

```bash
conda env create -f envs/dentrim_env.yml
```

Entornos que aparecen en la documentación del proyecto:

| Entorno | Periodo | Herramientas |
|---|---|---|
| `dentrim_env` | 2026 | samtools, Pychopper, edlib, Filtlong (Python 3.9) |
| `emu` | 2025–2026 | EMU (creado con `conda create -n emu -c conda-forge -c bioconda emu`) |
| `pipelinefenrir` | 2025 (fase piloto) | EMU y dependencias de los scripts de desarrollo |

El entorno con el que se ejecutaron `EMU_propio.py`, `rebuild_counts.py`, `rarefaccion.py` y `diversidad_mod.py` queda identificado en `anexos/datos/entorno/versiones_software.tsv` después de correr el script de recolección.

Dorado no se instala con conda: es un binario precompilado de Oxford Nanopore (ver el [Anexo B](../anexos/B_entorno_y_software.md)).
