# Archivos de configuración

Archivos pequeños de entrada que el pipeline necesita, además de los datos. Se copian desde el PC con [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh).

| Archivo | Lo usa | Contenido |
|---|---|---|
| `primers.fasta` | `dentrim_bam.py` (Pychopper `-b`) y `stats_fastq.py` (`--primers`) | Secuencias de los primers 16S, con los nombres `16sF` y `16sR` |
| `primers_config.txt` | `dentrim_bam.py` (Pychopper `-c`) | Una línea con las configuraciones válidas de primers por hebra |
| `piloto/primers_stesen.fasta`, `piloto/primers_stesen.txt` | Fase piloto (octubre 2025) | Primers 27F/1492R del conjunto de datos público *sterile sentinels* |

## Formato de `primers_config.txt`

Pychopper lee **solo la primera línea** del archivo (por eso no acepta la configuración escrita directamente después de `-c`; ver el [Anexo J](../anexos/J_problemas_y_soluciones.md)). El formato documentado por Pychopper es:

```
+:16sF,-16sR|-:16sR,-16sF
```

- `+:16sF,-16sR` → lectura en sentido directo: `16sF` al inicio y el reverso complementario de `16sR` al final.
- `-:16sR,-16sF` → lectura en sentido reverso: Pychopper la **reorienta** (reverso complementario) antes de escribirla.

> ⚠️ Si la segunda configuración empieza por `+:` (por ejemplo `+:16sF,-16sR|+:16sR,-16sF`), Pychopper marca también esas lecturas como hebra `+` y **no las reorienta**. Ver la verificación en el [Anexo C](../anexos/C_protocolo_comandos_parametros.md#verificaciones-pendientes-de-la-etapa-de-limpieza).
