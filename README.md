# 16S-FenrirPipeline — Anexos del trabajo de grado

Trabajo de grado: *Análisis del microbioma de Passiflora edulis f. edulis mediante secuenciación dirigida al gen 16S usando MinION, y predicción de las posibles interacciones bacteria-fago*. Johann Sebastian Gallego Sierra, Universidad El Bosque.

## Ramas del repositorio

| Rama | Contenido |
|---|---|
| [`main`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/main) | Pipeline 16S ONT: basecalling, estadísticas, limpieza, taxonomía y diversidad (Objetivos 1 y 3) |
| [`algoritmo`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/algoritmo) | Predicción de interacciones bacteria-fago con DeepPBI-KG (Objetivos 2 y 3) |
| [`anexos`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/anexos) | Anexos del documento: tablas de resultados, figuras, matrices de predicción y grafos |

## Índice de anexos

Los anexos 1 a 9 corresponden a los ya citados en el documento. Los anexos 10 a 19 amplían los resultados que no pudieron incluirse en el cuerpo del texto.

| Anexo | Carpeta | Contenido principal | Sección de la tesis |
|---|---|---|---|
| 1 | [Anexo_01_Metadata_PRJNA1020132](Anexo_01_Metadata_PRJNA1020132) | Metadata del BioProject, SraRunTable y datos de muestra del artículo | 10.1.3, 10.1.5 |
| 2 | [Anexo_02_Estadisticas_pre_limpieza](Anexo_02_Estadisticas_pre_limpieza) | Estadísticas básicas de las 88 muestras crudas | 10.1.5 |
| 3 | [Anexo_03_Estadisticas_post_limpieza](Anexo_03_Estadisticas_post_limpieza) | Estadísticas tras Filtlong + Pychopper y retención por muestra | 10.1.6 |
| 4 | [Anexo_04_Abundancia_generos_referencia](Anexo_04_Abundancia_generos_referencia) | Abundancias relativas por género del artículo | 10.1.7 |
| 5 | [Anexo_05_Abundancia_especies_referencia](Anexo_05_Abundancia_especies_referencia) | Abundancias relativas por especie del artículo | 10.1.7 |
| 6 | [Anexo_06_Diferencias_por_genero](Anexo_06_Diferencias_por_genero) | Medias propias vs. artículo y correlaciones (género) | 10.1.7 |
| 7 | [Anexo_07_Diferencias_por_especie](Anexo_07_Diferencias_por_especie) | Medias propias vs. artículo y correlaciones (especie) | 10.1.7 |
| 8 | [Anexo_08_Comparacion_modelos_Dorado](Anexo_08_Comparacion_modelos_Dorado) | Comparación FAST/HAC/SUP global y por muestra | 10.3.1 |
| 9 | [Anexo_09_Reportes_pycoQC_HAC_SUP](Anexo_09_Reportes_pycoQC_HAC_SUP) | Reportes interactivos de calidad pycoQC (HAC y SUP) | 10.3.1 |
| 10 | [Anexo_10_Diversidad_PRJNA1020132](Anexo_10_Diversidad_PRJNA1020132) | Diversidad alfa/beta, PCoA y PERMANOVA del conjunto de referencia | 10.1.8 |
| 11 | [Anexo_11_Limpieza_gulupa_HAC_SUP](Anexo_11_Limpieza_gulupa_HAC_SUP) | Limpieza por muestra de gulupa y reportes NanoPlot | 10.3.2 |
| 12 | [Anexo_12_Taxonomia_EMU_gulupa](Anexo_12_Taxonomia_EMU_gulupa) | Tablas de EMU (abundancias, conteos, taxonomía, top 20) y rarefacción | 10.3.3 |
| 13 | [Anexo_13_Diversidad_gulupa_HAC_SUP](Anexo_13_Diversidad_gulupa_HAC_SUP) | Diversidad alfa/beta, PCoA y PERMANOVA de gulupa | 10.3.4 |
| 14 | [Anexo_14_Genomas_bacterianos_Obj2](Anexo_14_Genomas_bacterianos_Obj2) | 192 genomas descargados para el Objetivo 2 | 10.2.2 |
| 15 | [Anexo_15_Predicciones_DeepPBI-KG_Obj2](Anexo_15_Predicciones_DeepPBI-KG_Obj2) | Matriz completa de predicciones y tablas filtradas (umbral 0.85) | 10.2.4, 10.2.5 |
| 16 | [Anexo_16_Redes_Gephi_Obj2](Anexo_16_Redes_Gephi_Obj2) | Grafos de Gephi (red completa y top 4) | 10.2.4 |
| 17 | [Anexo_17_Genomas_bacterianos_gulupa](Anexo_17_Genomas_bacterianos_gulupa) | 79 genomas descargados para gulupa (HAC) | 9.3.6, 10.3.5 |
| 18 | [Anexo_18_Predicciones_DeepPBI-KG_gulupa](Anexo_18_Predicciones_DeepPBI-KG_gulupa) | Matriz completa de predicciones y tablas filtradas (umbral 0.85) | 10.3.5, 10.3.6 |
| 19 | [Anexo_19_Redes_Gephi_gulupa](Anexo_19_Redes_Gephi_gulupa) | Grafos de Gephi y estadísticas de todas las redes | 10.3.5 |


## Notas

- Las tablas *Top5_interacciones_por_taxon* (anexos 15 y 18) están ordenadas por puntaje y tienen columnas vacías de "Nivel de evidencia" y "Referencia" para completar con la verificación bibliográfica.
- Los archivos .zip contienen matrices CSV de gran tamaño; se descomprimen con cualquier gestor de archivos.
- Los archivos .gephi se abren con Gephi 0.10 o superior; los .html se abren con cualquier navegador.
