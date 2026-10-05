# 16S-FenrirPipeline — Anexos del trabajo de grado

Trabajo de grado: *Análisis del microbioma de Passiflora edulis f. edulis mediante secuenciación dirigida al gen 16S usando MinION, y predicción de las posibles interacciones bacteria-fago*. Johann Sebastian Gallego Sierra, Universidad El Bosque.

## Ramas del repositorio

| Rama | Contenido |
|---|---|
| [`pipeline`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/main) | Pipeline 16S ONT: basecalling, estadísticas, limpieza, taxonomía y diversidad (Objetivos 1 y 3) |
| [`algoritmo`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/algoritmo) | Predicción de interacciones bacteria-fago con DeepPBI-KG (Objetivos 2 y 3) |
| [`anexos`](https://github.com/xXFenrir/16S-FenrirPipeline/tree/anexos) | Anexos del documento: tablas de resultados, figuras, matrices de predicción y grafos |

## Índice de anexos

Los anexos 1 a 7 están descritos en el apartado de Anexos del documento, y los anexos 1, 2, 3, 8 y 9 se citan en el texto. Los anexos 10 a 21 amplían los resultados y la información que no pudieron incluirse en el cuerpo del texto. La última columna indica la sección de Resultados del documento con la que se relaciona cada anexo.

| Anexo | Carpeta | Contenido principal | Sección de Resultados donde se usa |
|---|---|---|---|
| 1 | [Anexo_01_Metadata_PRJNA1020132](Anexo_01_Metadata_PRJNA1020132) | Metadata del BioProject, SraRunTable y datos de muestra del artículo | Obj. 1 › Selección de repositorios; Diseño conceptual: Reporte de estadísticas básicas |
| 2 | [Anexo_02_Estadisticas_pre_limpieza](Anexo_02_Estadisticas_pre_limpieza) | Estadísticas básicas de las 88 muestras crudas | Obj. 1 › Diseño conceptual: Reporte de estadísticas básicas; Diseño conceptual: Denoising y trimming |
| 3 | [Anexo_03_Estadisticas_post_limpieza](Anexo_03_Estadisticas_post_limpieza) | Estadísticas tras Filtlong + Pychopper y retención por muestra | Obj. 1 › Diseño conceptual: Denoising y trimming |
| 4 | [Anexo_04_Abundancia_generos_referencia](Anexo_04_Abundancia_generos_referencia) | Abundancias relativas por género del artículo | Obj. 1 › Diseño conceptual: Taxonomía |
| 5 | [Anexo_05_Abundancia_especies_referencia](Anexo_05_Abundancia_especies_referencia) | Abundancias relativas por especie del artículo | Obj. 1 › Diseño conceptual: Taxonomía |
| 6 | [Anexo_06_Diferencias_por_genero](Anexo_06_Diferencias_por_genero) | Medias propias vs. artículo y correlación de Pearson (género) | Obj. 1 › Diseño conceptual: Taxonomía |
| 7 | [Anexo_07_Diferencias_por_especie](Anexo_07_Diferencias_por_especie) | Medias propias vs. artículo y correlación de Pearson (especie) | Obj. 1 › Diseño conceptual: Taxonomía |
| 8 | [Anexo_08_Comparacion_modelos_Dorado](Anexo_08_Comparacion_modelos_Dorado) | Comparación FAST/HAC/SUP global y por muestra | Obj. 3 › Basecalling y evaluación de calidad de secuencias |
| 9 | [Anexo_09_Reportes_pycoQC_HAC_SUP](Anexo_09_Reportes_pycoQC_HAC_SUP) | Reportes interactivos de calidad pycoQC (HAC y SUP) | Obj. 3 › Basecalling y evaluación de calidad de secuencias |
| 10 | [Anexo_10_Diversidad_PRJNA1020132](Anexo_10_Diversidad_PRJNA1020132) | Shannon y riqueza observada, Bray-Curtis/Jaccard, PCoA y PERMANOVA del conjunto de referencia | Obj. 1 › Diseño conceptual: Diversidad |
| 11 | [Anexo_11_Limpieza_gulupa_HAC_SUP](Anexo_11_Limpieza_gulupa_HAC_SUP) | Limpieza por muestra de gulupa y reportes NanoPlot | Obj. 3 › Filtrado de calidad y remoción de ruido (Denoising y Trimming) |
| 12 | [Anexo_12_Taxonomia_EMU_gulupa](Anexo_12_Taxonomia_EMU_gulupa) | Tablas de EMU (abundancias, conteos, taxonomía, top 20) y rarefacción | Obj. 3 › Clasificación taxonómica y análisis de rarefacción |
| 13 | [Anexo_13_Diversidad_gulupa_HAC_SUP](Anexo_13_Diversidad_gulupa_HAC_SUP) | Shannon y riqueza observada, Bray-Curtis/Jaccard, PCoA y PERMANOVA de gulupa | Obj. 3 › Diversidad Alfa y Beta |
| 14 | [Anexo_14_Genomas_bacterianos_Obj2](Anexo_14_Genomas_bacterianos_Obj2) | 192 genomas descargados para el Objetivo 2 | Obj. 2 › Preparación de datos de entrada |
| 15 | [Anexo_15_Predicciones_DeepPBI-KG_Obj2](Anexo_15_Predicciones_DeepPBI-KG_Obj2) | Matriz completa de predicciones y tablas filtradas (umbral 0.85) | Obj. 2 › Filtrado de probabilidades; Comprobación con literatura |
| 16 | [Anexo_16_Redes_Gephi_Obj2](Anexo_16_Redes_Gephi_Obj2) | Grafos de Gephi (red completa y top 4) | Obj. 2 › Filtrado de probabilidades |
| 17 | [Anexo_17_Genomas_bacterianos_gulupa](Anexo_17_Genomas_bacterianos_gulupa) | 79 genomas descargados para gulupa (HAC) | Obj. 3 › Interacciones bacteria-fago en muestras de gulupa |
| 18 | [Anexo_18_Predicciones_DeepPBI-KG_gulupa](Anexo_18_Predicciones_DeepPBI-KG_gulupa) | Matriz completa de predicciones y tablas filtradas (umbral 0.85) | Obj. 3 › Interacciones bacteria-fago en muestras de gulupa; Comprobación de interacciones en gulupa |
| 19 | [Anexo_19_Redes_Gephi_gulupa](Anexo_19_Redes_Gephi_gulupa) | Grafos de Gephi y estadísticas de todas las redes | Obj. 3 › Interacciones bacteria-fago en muestras de gulupa |
| 20 | [Anexo_20_Diferencias_por_familia](Anexo_20_Diferencias_por_familia) | Abundancias relativas por familia (propias y del artículo), medias y correlación de Pearson (familia) | Obj. 1 › Diseño conceptual: Taxonomía |
| 21 | [Anexo_21_Herramientas_pipelines_referencia](Anexo_21_Herramientas_pipelines_referencia) | Herramientas por etapa en seis pipelines 16S de referencia | Obj. 1 › Caracterización de herramientas bioinformáticas (Tabla 1) |

## Descripción de los anexos

**Anexo 1. Metadata del proyecto PRJNA1020132**

Contiene la información contextual y experimental del proyecto de referencia empleado para la comparación taxonómica, incluyendo códigos de acceso, descripción de muestras, condiciones de secuenciación y origen de los datos públicos utilizados como referencia (Erlandson et al., 2024).

**Anexo 2. Estadísticas básicas antes del proceso de limpieza**

Resume las métricas de calidad obtenidas en las lecturas crudas de las 88 muestras del proyecto PRJNA1020132, previas al filtrado, como número total de lecturas, longitud promedio, contenido GC, N50 y Q-scores. Permite evaluar la calidad inicial de las secuencias obtenidas mediante tecnología *Oxford Nanopore Technologies*.

**Anexo 3. Estadísticas básicas después del proceso de limpieza**

Presenta las mismas métricas del Anexo 2 calculadas sobre las lecturas que superaron el recorte de primers con Pychopper y el filtrado por longitud (1000–1700 pb) y calidad con Filtlong. Incluye, para cada muestra, el porcentaje de lecturas y bases retenidas y removidas, lo que permite valorar el efecto de la limpieza.

**Anexo 4. Abundancia relativa de géneros del conjunto de referencia**

Contiene las abundancias relativas por género reportadas por Erlandson et al. (2024) para cada muestra del proyecto PRJNA1020132. Sirve como punto de comparación para validar la clasificación taxonómica obtenida con el pipeline propio.

**Anexo 5. Abundancia relativa de especies del conjunto de referencia**

Contiene las abundancias relativas por especie reportadas por Erlandson et al. (2024) para cada muestra del proyecto PRJNA1020132. Complementa el Anexo 4 en la comparación a nivel de especie.

**Anexo 6. Diferencia de abundancias relativas por género**

Compara la abundancia relativa media de los 344 géneros presentes en ambos análisis, obtenida con el pipeline propio y reportada en el artículo, junto con la diferencia entre ellas. Incluye la correlación de Pearson entre ambas series y las abundancias por género de cada muestra procesada en este estudio.

**Anexo 7. Diferencia de abundancias relativas por especie**

Presenta la misma comparación del Anexo 6 a nivel de especie, para las 853 especies detectadas tanto en este estudio como en el artículo de referencia. Incluye la correlación de Pearson y las abundancias por especie de cada muestra.

**Anexo 8. Comparación de los modelos de basecalling de Dorado**

Resume el desempeño de los modelos FAST, HAC y SUP de Dorado en la corrida de las muestras de gulupa: lecturas totales, aprobadas (*pass*) y rechazadas (*fail*), bases, longitud media, N50 y QScore medio. Incluye el detalle por muestra y gráficas de lecturas aprobadas y de calidad, que sustentan la elección de los modelos HAC y SUP.

**Anexo 9. Reportes de calidad pycoQC de los modelos HAC y SUP**

Contiene los reportes interactivos de pycoQC generados a partir del resumen de secuenciación de Dorado para los modelos HAC y SUP. Muestran el resumen general de la corrida, la distribución de longitud y calidad de las lecturas, su evolución durante el tiempo de secuenciación y el número de lecturas por código de barras.

**Anexo 10. Diversidad microbiana del proyecto PRJNA1020132**

Reúne el análisis de diversidad del conjunto de referencia procesado con el pipeline propio: índices de Shannon y riqueza observada por muestra, prueba de Wilcoxon (Mann-Whitney U) entre las rotaciones de cultivo CS y CSSwP, matrices de disimilitud de Bray-Curtis y Jaccard, coordenadas del PCoA y resultado de la PERMANOVA. Incluye además las tablas de abundancia, conteos y taxonomía de EMU con las que se calcularon.

**Anexo 11. Limpieza de las lecturas de gulupa con los modelos HAC y SUP**

Detalla, para cada una de las 73 muestras de gulupa y para ambos modelos, el número de lecturas antes de la limpieza, después de Pychopper y después de Filtlong, junto con las bases, la longitud, el N50 y el QScore antes y después del proceso. Incluye los reportes de NanoPlot de las lecturas limpias y gráficas comparativas entre HAC y SUP.

**Anexo 12. Clasificación taxonómica con EMU de las muestras de gulupa**

Contiene los resultados de EMU para las muestras de gulupa con al menos 500 lecturas limpias (28 con HAC y 26 con SUP): abundancias relativas, conteos y taxonomía de cada taxón, y los 20 taxones más abundantes por especie, género y familia. Incluye las curvas de rarefacción y un resumen general por modelo.

**Anexo 13. Diversidad microbiana del rizobioma de gulupa**

Presenta el análisis de diversidad de las muestras de gulupa con los modelos HAC y SUP, agrupadas por sistema de manejo (Empresarial, Campesina y Agroecológica): índices de Shannon y riqueza observada por muestra y por sistema, valores p de la prueba de Wilcoxon (Mann-Whitney U), matrices de Bray-Curtis y Jaccard, PCoA y PERMANOVA.

**Anexo 14. Genomas bacterianos descargados para el Objetivo 2**

Lista los 192 genomas bacterianos descargados de NCBI para los taxones del proyecto PRJNA1020132 con abundancia relativa promedio de al menos 0,1 %, indicando su TaxID y si corresponden a la especie o a un genoma representante del género. Incluye el histograma de abundancias con el umbral de corte empleado.

**Anexo 15. Predicciones de DeepPBI-KG para el Objetivo 2**

Contiene la matriz completa de predicciones de DeepPBI-KG para los 674 496 pares fago-bacteria evaluados (3513 fagos y 192 bacterias), con las probabilidades de los modelos de genes clave y de genoma completo. Incluye las interacciones que superaron el umbral de 0,85, su reducción a los cuatro taxones dominantes, las cinco interacciones con mayor puntaje por taxón y las gráficas de densidad y reducción.

**Anexo 16. Redes de interacción bacteria-fago del Objetivo 2**

Contiene los grafos de Gephi de las interacciones predichas para el conjunto de referencia: la red completa con umbral de 0,85 y la red reducida a los cuatro taxones bacterianos dominantes, en formato editable y como imagen.

**Anexo 17. Genomas bacterianos descargados para gulupa**

Lista los 79 genomas bacterianos descargados de NCBI para los taxones de gulupa (modelo HAC) con abundancia relativa promedio de al menos 0,1 %, junto con el resultado de cada búsqueda. Incluye el registro completo de la descarga.

**Anexo 18. Predicciones de DeepPBI-KG para gulupa**

Contiene la matriz completa de predicciones para los 277 527 pares evaluados entre los 3513 fagos de referencia y las 79 bacterias de gulupa, con las probabilidades de genes clave y de genoma completo. Incluye las interacciones que superaron el umbral de 0,85, su reducción a los cuatro taxones dominantes, hasta cinco interacciones con mayor puntaje por taxón y las gráficas de densidad y reducción.

**Anexo 19. Redes de interacción bacteria-fago de gulupa**

Contiene los grafos de Gephi de las interacciones predichas para gulupa con el modelo HAC: la red completa y la red de los cuatro taxones dominantes, en formato editable, PDF e imagen. Incluye las estadísticas de todas las redes construidas (nodos, aristas, densidad y grado promedio).

**Anexo 20. Diferencia de abundancias relativas por familia**

Compara la abundancia relativa media de las 132 familias presentes en ambos análisis, obtenida con el pipeline propio y reportada en el artículo, junto con la correlación de Pearson entre ambas series. Incluye las abundancias por familia de cada muestra, tanto propias como del artículo.

**Anexo 21. Herramientas bioinformáticas reportadas por artículo en cada etapa del procesamiento del gen 16S**

Resume las herramientas empleadas en seis estudios que analizan el gen 16S rRNA secuenciado con *Oxford Nanopore Technologies*, organizadas por etapa: basecalling, demultiplexing, denoising y trimming, clustering, asignación taxonómica, filogenia y diversidad. Corresponde a la Tabla 1 del documento; el guion (-) indica que la etapa no se reportó o no aplica.

## Notas

- Las tablas *Top5_interacciones_por_taxon* (anexos 15 y 18) están ordenadas por puntaje y tienen columnas vacías de "Nivel de evidencia" y "Referencia" para completar con la verificación bibliográfica.
- Los archivos .zip contienen matrices CSV de gran tamaño; se descomprimen con cualquier gestor de archivos.
- Los archivos .gephi se abren con Gephi 0.10 o superior; los .html se abren con cualquier navegador.
