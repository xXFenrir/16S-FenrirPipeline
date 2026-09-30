# 16S-FenrirPipeline — Algoritmo de interacciones bacteria-fago

Predicción *in silico* de interacciones bacteria-fago a partir de los perfiles taxonómicos 16S obtenidos con el pipeline de la rama `main` (Objetivos 2 y 3 del trabajo de grado *Análisis del microbioma de Passiflora edulis f. edulis mediante secuenciación dirigida al gen 16S usando MinION, y predicción de las posibles interacciones bacteria-fago*, Johann Sebastian Gallego Sierra, Universidad El Bosque).

## Ramas del repositorio

| Rama | Contenido |
|---|---|
| [`main`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/main) | Pipeline 16S ONT: basecalling, estadísticas, limpieza, taxonomía y diversidad (Objetivos 1 y 3) |
| [`algoritmo`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/algoritmo) | Predicción de interacciones bacteria-fago con DeepPBI-KG (Objetivos 2 y 3) |
| [`anexos`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/anexos) | Anexos del documento: tablas de resultados, figuras, matrices de predicción y grafos |

## Estructura de esta rama

```
01_preparacion_entradas/    Descarga de genomas NCBI, anotación Prokka y alineamiento BLAST
02_prediccion_DeepPBI-KG/   Script de predicción adaptado y gráficas de probabilidades
03_filtrado_y_redes/        Umbral 0.85, tablas para Gephi y estadísticas de red
```

## Flujo general

1. Tabla de abundancias relativas de EMU (rama `main`) → filtro por abundancia promedio ≥ 0.001.
2. Descarga de un genoma de referencia por taxón desde NCBI.
3. Anotación con Prokka y alineamiento con BLASTn contra los 3513 fagos de DeepPBI-KG.
4. Predicción con la red neuronal preentrenada de DeepPBI-KG para todos los pares fago-bacteria.
5. Filtro con umbral de 0.85 sobre el promedio de las salidas de genes clave y genoma completo.
6. Redes en Gephi (completa y de los cuatro taxones dominantes) y verificación bibliográfica.

## Requisitos

Además de los scripts de esta rama se necesita el repositorio original de DeepPBI-KG (https://github.com/Tongqing-Wei/DeepPBI-KG), del que se toman la carpeta `model/` (pesos y escaladores), el panel de fagos de referencia y su entorno (`requirements.txt`), junto con Prokka, BLAST+ y NCBI `datasets`.

## Descripción

Las abundancias relativas que entrega EMU se acoplan al modelo **DeepPBI-KG** (Wei et al., 2024; código original en https://github.com/Tongqing-Wei/DeepPBI-KG), que estima la probabilidad de interacción lítica entre cada par fago-bacteria a partir de genes clave (Random Forest + red neuronal profunda) y de características de genoma completo. Se usan los pesos, escaladores y el panel de 3513 fagos de referencia del repositorio original. Los scripts de esta rama son adaptaciones o scripts propios; el código base de DeepPBI-KG pertenece a sus autores.

- **Objetivo 2:** validación sobre el conjunto de referencia PRJNA1020132 (192 genomas bacterianos, 674.496 pares evaluados).
- **Objetivo 3:** aplicación sobre la rizósfera de gulupa con las abundancias del modelo HAC (79 genomas bacterianos).

## 1. Preparación de entradas (`01_preparacion_entradas/`)

- `descarga_genomas.py` filtra la tabla de abundancias relativas de EMU por abundancia promedio (≥ 0.001 por defecto) y descarga, con NCBI `datasets`, un genoma por taxón con búsqueda en cascada: TaxID exacto → nombre científico → genoma representante del género. Prioriza ensamblajes completos y, si no hay, cromosoma/scaffold.
- `bact_data.py` es la versión anterior del mismo procedimiento, usada en el Objetivo 2 (rutas fijas).
- `prokka_blast_mod.sh` es la versión modificada de `prokka_blast.sh` de DeepPBI-KG. Valida que Prokka y BLAST estén disponibles antes de empezar, anota cada genoma con Prokka y lo alinea con BLASTn contra la base de fagos/bacterias, además de registrar el avance en consola.
- `anotar_todo.py` ejecuta Prokka por lotes sobre carpetas de genomas de fagos y bacterias.

```
python3 descarga_genomas.py \
  -t EMUhac_results/feature_table_relabund_hac_results.tsv \
  -ob bacterias_fna \
  -min 0.001
```

**Resultado esperado:** un `.fna` por taxón (`taxid_<id>_<especie>.fna` o `GENERO_<género>_<taxid>.fna`), una carpeta Prokka por genoma (con el `.gbk` de CDS) y un reporte BLAST `.out` por genoma.

## 2. Predicción (`02_prediccion_DeepPBI-KG/`)

`DeepPBI-KG_propio.py` es el script de predicción de DeepPBI-KG adaptado para evaluar todos los pares fago-bacteria por lotes en CPU, con procesamiento en paralelo. Usa el mismo escalador y los mismos pesos preentrenados que el modelo original, y entrega para cada par `key_gene_output` (potencial de infección a nivel de genes clave) y `wgs_output` (coexistencia a nivel de genoma completo).

```
python3 DeepPBI-KG_propio.py \
  --phage_annotation fagos_annot --bacterium_annotation bacterias_annot \
  --phage_align_res fagos_align --bacterium_align_res bacterias_align \
  --phage_raw_data fagos_fna --bacterium_raw_data bacterias_fna \
  --model model --template template --output output/result.csv
```

`visualizacion_resultados.py` (Objetivo 2) y `visualizacion_resultados_mod.py` (Objetivo 3, `--dir_base`, `--umbral 0.85`) generan el diagrama de densidad por cuadrantes (WGS vs. genes clave) y el embudo de reducción de interacciones.

## 3. Filtrado y redes (`03_filtrado_y_redes/`)

- `filtro_redes.py` aplica el umbral (`--threshold 0.85`) y genera tres tablas: solo genes clave, solo WGS y la métrica combinada (promedio de ambas), que es la que se usa para las redes.
- `integrate_seq_mod.py` compara los tres conjuntos filtrados (diagrama de Venn y gráficas).
- `generar_nodos_gephi.py` / `generar_nodos_gephi_v2.py` convierten los identificadores de NCBI en nombres legibles y exportan las aristas (`Source`, `Target`, `Weight`) y los nodos para Gephi.
- `filtrar_top_taxones.py` / `filtrar_top_taxones_comprobacion.py` reducen la red a los cuatro taxones bacterianos más abundantes; `rank_top_taxones_en_interacciones.py` identifica esos taxones dentro de la tabla de interacciones.
- `stats_redes.py` calcula nodos, aristas, densidad, grado promedio y los nodos de mayor grado de cada red.

**Resultado esperado:** matrices de interacciones confiables (umbral 0.85), redes completas y filtradas a los cuatro taxones dominantes (visualizadas en Gephi con distribución dirigida por fuerzas y tamaño de nodo proporcional al grado), y un resumen de estadísticas de red.
