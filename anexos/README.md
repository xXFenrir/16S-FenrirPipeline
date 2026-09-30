# Anexos de la tesis

Material complementario de la tesis sobre el rizobioma de cultivos de gulupa bajo tres sistemas de manejo agrícola, analizado con secuenciación del gen 16S rRNA completo en Oxford Nanopore. Aquí está lo que no cabe en el cuerpo del documento: comandos completos, parámetros, versiones de software, tablas por muestra, resultados complementarios y la historia del desarrollo del pipeline.

## Índice

| Anexo | Título | Contenido | Estado |
|---|---|---|---|
| [A](A_disponibilidad_y_reproducibilidad.md) | Disponibilidad del código y los datos | Estructura del repositorio, dónde están los datos, cómo reproducir el análisis | ✅ Redactado |
| [B](B_entorno_y_software.md) | Entorno computacional y versiones de software | Equipo, entornos conda, herramientas con versión y referencia, base de datos de EMU | 🟡 Faltan versiones exactas y datos del equipo |
| [C](C_protocolo_comandos_parametros.md) | Protocolo bioinformático detallado | Diagrama del flujo, comandos exactos y tabla de parámetros de cada etapa | ✅ Redactado · ⚠️ 2 verificaciones |
| [D](D_corrida_muestras_primers.md) | Corrida, muestras y primers | Datos de la corrida de MinION, correspondencia barcode–finca–sistema, primers | 🟡 Faltan metadatos, primers y reporte de MinKNOW |
| [E](E_control_de_calidad.md) | Control de calidad de las lecturas | Definición de métricas; tablas por muestra tras basecalling, antes/después de la limpieza | 🟡 Faltan tablas |
| [F](F_resultados_taxonomicos.md) | Resultados taxonómicos complementarios | Tablas completas de abundancia y taxonomía, curvas de rarefacción | 🟡 Faltan tablas y figuras |
| [G](G_diversidad_complementaria.md) | Diversidad: definiciones y resultados | Fórmulas de índices y pruebas; alfa por muestra, distancias, PCoA, PERMANOVA | ✅ Fórmulas · 🟡 Faltan tablas · ⚠️ 1 verificación |
| [H](H_fase_piloto_datos_publicos.md) | Fase piloto con datos públicos | Conjunto *sterile sentinels* (Erlandson et al., 2024) y comparación con su pipeline | ✅ Redactado |
| [I](I_evolucion_y_decisiones.md) | Evolución del pipeline y decisiones | Cronología, decisiones con alternativas y motivos, evolución de parámetros, HAC y SUP | ✅ Redactado · ⏳ 3 motivos por completar |
| [J](J_problemas_y_soluciones.md) | Problemas encontrados y soluciones | Errores reales del desarrollo y cómo se resolvieron | ✅ Redactado |
| [historico/](historico/README.md) | Versiones de desarrollo de los scripts | Los cuatro scripts de octubre de 2025, recuperados del historial de Git | ✅ |

Leyenda: ✅ listo con la información del repositorio · 🟡 texto listo, faltan archivos del equipo de análisis · ⏳ dato por completar · ⚠️ punto técnico por confirmar.

## Cómo completar los anexos desde el PC

Todo lo marcado con 🟡 y ⏳ está en el equipo donde se hizo el análisis. El script [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh) lo busca en las rutas del proyecto y lo copia al repositorio. Solo lee y copia: no borra ni modifica nada fuera del repositorio, no copia datos crudos (POD5, BAM, FASTQ) ni archivos de más de 20 MB, y no sube nada a GitHub.

```bash
# 1. Clonar el repositorio (solo la primera vez) y entrar a la rama con los anexos
git clone https://github.com/xXFenrir/16S-FenrirPipeline.git
cd 16S-FenrirPipeline
git checkout claude/keen-lamport-t4jomc

# 2. Ver qué encontraría, sin copiar nada
bash herramientas/recolectar_anexos.sh --simular

# 3. Copiar
bash herramientas/recolectar_anexos.sh

# 4. Revisar el resultado
less anexos/INVENTARIO_RECOLECCION.md
git status

# 5. Subir a GitHub
git add -A
git commit -m "Agregar scripts finales, entornos y datos de los anexos"
git push
```

Si la carpeta de la tesis no está en `~/Documentos/Tesis`, se indica con `--tesis /ruta/a/Tesis`.

El inventario (`anexos/INVENTARIO_RECOLECCION.md`) lista lo copiado, lo que no se encontró, lo omitido por tamaño, y el resultado de las tres verificaciones de abajo.

## Puntos por verificar antes de la entrega

1. **Filtro de calidad de Filtlong** ([Anexo C, verificación 1](C_protocolo_comandos_parametros.md#verificaciones-pendientes-de-la-etapa-de-limpieza)). Filtlong interpreta `--min_mean_q` como porcentaje de exactitud, no como Phred. Si `dentrim_bam.py` le pasa `12` tal cual, ese filtro no descarta nada y el filtro de calidad efectivo es Q8.
2. **Orientación de las lecturas** ([Anexo C, verificación 2](C_protocolo_comandos_parametros.md#verificaciones-pendientes-de-la-etapa-de-limpieza)). Confirmar que `primers_config.txt` marca la segunda configuración con `-:`; si no, las lecturas reversas no se reorientaron (la taxonomía no cambia).
3. **Distancia de Jaccard** ([Anexo G, verificación 3](G_diversidad_complementaria.md#verificación-3-jaccard-calculado-con-pdist-sobre-abundancias)). Con SciPy anterior a 1.15, `pdist(abundancias, "jaccard")` no calcula presencia/ausencia.
4. **Orden de Pychopper y Filtlong.** El README dice en la sección de Filtlong que se ejecuta *antes* de Pychopper, pero la lista de pasos (y los parámetros `--post-minlen`/`--post-maxlen`) indican que va *después*. Dejar una sola versión, la que hace el script.
5. **Nombres de los sistemas.** El texto usa "Campesino" y "Agroecológico"; los metadatos y el comando usan "Campesina" y "Agroecológica". Unificar en el documento.
6. **Comparaciones múltiples.** Las pruebas de Mann-Whitney se hacen para tres pares de sistemas. Indicar si se corrigieron los valores p (por ejemplo Holm o Benjamini-Hochberg) o por qué no.
7. **Dispersión en PERMANOVA** (opcional). PERMANOVA puede dar significativo si los grupos difieren en dispersión y no en centroide. Una prueba PERMDISP (`skbio.stats.distance.permdisp`) ayuda a interpretar el resultado.
8. **Motivos pendientes** ([Anexo I](I_evolucion_y_decisiones.md#i2-decisiones-metodológicas)): por qué HAC y no SUP en los resultados finales, por qué `edlib` y no `phmm`, y el motivo del cambio a Pychopper para recortar primers.
9. **Privacidad de los metadatos** antes de hacer público el repositorio ([Anexo D](D_corrida_muestras_primers.md#d2-muestras-y-diseño)).
10. **Publicación**: hacer público el repositorio, crear una versión con DOI y depositar las lecturas en NCBI SRA ([Anexo A](A_disponibilidad_y_reproducibilidad.md)).

## Cómo citar los anexos en el documento

Cada anexo tiene una dirección estable una vez que la rama se integra a `main`, por ejemplo:

> Anexo C. Protocolo bioinformático detallado: comandos y parámetros. Disponible en https://github.com/xXFenrir/16S-FenrirPipeline/blob/main/anexos/C_protocolo_comandos_parametros.md

Si la universidad pide los anexos dentro del PDF, los archivos `.md` se pueden convertir con Pandoc (`pandoc anexos/C_protocolo_comandos_parametros.md -o anexo_C.docx`) y ajustar el formato en Word.
