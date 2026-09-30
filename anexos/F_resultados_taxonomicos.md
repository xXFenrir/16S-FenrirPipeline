# Anexo F. Resultados taxonómicos complementarios

Tablas completas de la clasificación con EMU y curvas de rarefacción por muestra. En el documento se muestran los taxones más abundantes y un resumen de la profundidad de secuenciación; aquí queda el detalle.

Los archivos se copian con [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh) en `anexos/datos/taxonomia/`.

## F.1 Tablas de abundancia y taxonomía

| Archivo | Contenido | Formato |
|---|---|---|
| `tabla_abundancia_relativa.tsv` | Abundancia relativa de cada taxón en cada muestra (fracción entre 0 y 1) | Filas: taxones; columnas: muestras |
| `tabla_conteos.tsv` | Conteos enteros reconstruidos (`rebuild_counts.py`) | Filas: taxones; columnas: muestras |
| `taxonomia.tsv` | Linaje completo de cada taxón (superkingdom → species) | Una fila por ID de taxón |

Cada taxón se identifica como `rank|tax_id|nombre` (por ejemplo `species|562|Escherichia coli`), para que dos taxones con el mismo nombre y distinto `tax_id` de NCBI no se mezclen.

Tablas de la corrida HAC: `anexos/datos/taxonomia/hac/` ⏳. Si existe la corrida SUP: `anexos/datos/taxonomia/sup/` ⏳.

## F.2 Resumen de la clasificación

| Dato | HAC |
|---|---|
| Muestras clasificadas | ⏳ |
| Muestras descartadas por tener menos de 500 lecturas limpias | ⏳ (listado) |
| Especies detectadas en total | ⏳ |
| Especies por muestra (mediana y rango) | ⏳ |
| Filos, familias y géneros detectados | ⏳ |

## F.3 Curvas de rarefacción

Una curva por muestra (`rarefaccion.py`, 30 profundidades, 10 iteraciones). Figuras: `anexos/datos/taxonomia/hac/figuras/` ⏳.

Cómo leerlas: el eje X es el número de lecturas submuestreadas y el eje Y el número medio de taxones observados. Si la curva se aplana antes de llegar a la profundidad total de la muestra, secuenciar más no habría agregado muchos taxones nuevos. Si sigue subiendo, la riqueza de esa muestra está subestimada.

| Muestra | Lecturas | Taxones observados | ¿Alcanza meseta? |
|---|---|---|---|
| ⏳ | | | |

## F.4 Composición a otros niveles taxonómicos

EMU clasifica a nivel de especie. Para mostrar la composición por filo, familia o género se pueden agrupar las abundancias de `tabla_abundancia_relativa.tsv` con el linaje de `taxonomia.tsv`, o usar el comando de EMU:

```bash
emu collapse-taxonomy <archivo_rel-abundance.tsv> genus
```

Si en el documento se incluyen solo los 10–20 taxones más abundantes, la tabla completa de este anexo permite consultar el resto.
