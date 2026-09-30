# Anexo G. Diversidad: definiciones y resultados complementarios

Fórmulas de los índices y pruebas usados en `diversidad_mod.py`, y tablas completas de diversidad alfa, distancias, PCoA y PERMANOVA. En el documento se muestran las figuras y los resultados principales; aquí van las definiciones y los valores por muestra.

Los archivos se copian con [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh) en `anexos/datos/diversidad/`.

## G.1 Definiciones y fórmulas

Notación: $p_i$ es la abundancia relativa del taxón $i$ en una muestra, $S$ el número de taxones y $x_{ij}$ la abundancia del taxón $i$ en la muestra $j$.

**Riqueza observada.** Número de taxones presentes en la muestra:

$$S_{obs} = \sum_{i} \mathbb{1}(x_i > 0)$$

**Índice de Simpson** (tal como lo calcula scikit-bio con la métrica `simpson`, también llamado índice de Gini-Simpson). Probabilidad de que dos lecturas tomadas al azar pertenezcan a taxones distintos; va de 0 (un solo taxón) a casi 1 (muchos taxones con abundancias parecidas):

$$D = 1 - \sum_{i=1}^{S} p_i^2$$

**Disimilitud de Bray-Curtis** entre las muestras $j$ y $k$. Cuantitativa: tiene en cuenta cuánto cambia la abundancia de cada taxón. Va de 0 (idénticas) a 1 (sin taxones en común):

$$BC_{jk} = \frac{\sum_i \lvert x_{ij} - x_{ik} \rvert}{\sum_i (x_{ij} + x_{ik})}$$

**Distancia de Jaccard.** Cualitativa: solo considera presencia o ausencia. $A$ y $B$ son los conjuntos de taxones presentes en cada muestra:

$$J_{jk} = 1 - \frac{\lvert A \cap B \rvert}{\lvert A \cup B \rvert}$$

**Análisis de Coordenadas Principales (PCoA).** A partir de la matriz de distancias $D$ entre $n$ muestras se construye $B = -\tfrac{1}{2} J D^{(2)} J$, donde $D^{(2)}$ son las distancias al cuadrado y $J = I - \tfrac{1}{n}\mathbf{1}\mathbf{1}^T$ es la matriz de centrado. Los vectores propios de $B$, escalados por la raíz de sus valores propios $\lambda_k$, son las coordenadas de las muestras. La proporción de variación que explica cada eje es $\lambda_k / \sum_k \lambda_k$; scikit-bio divide entre la suma de todos los valores propios (la traza de $B$). Bray-Curtis no es una métrica euclidiana, así que pueden aparecer valores propios negativos pequeños; scikit-bio avisa si son grandes.

**PERMANOVA** (Anderson, 2001). Prueba si los centroides de los grupos difieren en el espacio de la matriz de distancias. Con $N$ muestras en $g$ grupos:

$$F = \frac{SS_A / (g - 1)}{SS_W / (N - g)}, \qquad SS_T = \frac{1}{N} \sum_{i<j} d_{ij}^2, \qquad SS_W = \sum_{\text{grupos}} \frac{1}{n_g} \sum_{i<j \in g} d_{ij}^2, \qquad SS_A = SS_T - SS_W$$

El valor p se obtiene permutando las etiquetas de grupo 999 veces: $p = (\#\{F_{perm} \ge F_{obs}\} + 1) / (999 + 1)$. Con 999 permutaciones, el menor valor p posible es 0,001.

**Prueba U de Mann-Whitney.** Compara un índice alfa entre dos sistemas sin suponer normalidad. Con $n_1$ y $n_2$ muestras y $R_1$ la suma de rangos del grupo 1:

$$U_1 = R_1 - \frac{n_1 (n_1 + 1)}{2}$$

Se aplica a cada par de sistemas: Empresarial–Campesina, Empresarial–Agroecológica y Campesina–Agroecológica.

**Elipses de confianza al 95 %.** En la figura de PCoA se dibuja una elipse por sistema agrícola, calculada a partir de la dispersión de sus muestras en los dos primeros ejes.

## G.2 Diversidad alfa por muestra

Fuente: tabla de diversidad alfa de `diversidad_mod.py` ⏳

| Muestra | Sistema | Riqueza observada | Simpson |
|---|---|---|---|
| ⏳ | | | |

Resumen por sistema:

| Sistema | n | Riqueza (mediana [rango]) | Simpson (mediana [rango]) |
|---|---|---|---|
| Empresarial | ⏳ | ⏳ | ⏳ |
| Campesina | ⏳ | ⏳ | ⏳ |
| Agroecológica | ⏳ | ⏳ | ⏳ |

## G.3 Comparaciones de diversidad alfa entre sistemas

| Índice | Par de sistemas | U | p |
|---|---|---|---|
| Riqueza observada | Empresarial – Campesina | ⏳ | ⏳ |
| | Empresarial – Agroecológica | ⏳ | ⏳ |
| | Campesina – Agroecológica | ⏳ | ⏳ |
| Simpson | Empresarial – Campesina | ⏳ | ⏳ |
| | Empresarial – Agroecológica | ⏳ | ⏳ |
| | Campesina – Agroecológica | ⏳ | ⏳ |

## G.4 Matrices de distancia

Matrices completas de Bray-Curtis y Jaccard (muestra × muestra) y sus heatmaps: `anexos/datos/diversidad/` ⏳.

### Verificación 3. Jaccard calculado con `pdist` sobre abundancias

`scipy.spatial.distance.pdist(X, "jaccard")` solo calcula la distancia de presencia/ausencia descrita arriba si recibe datos booleanos (`X > 0`) **o** si la versión de SciPy es 1.15 o posterior. En versiones anteriores, con abundancias cuenta como diferente todo taxón cuya abundancia no sea exactamente igual en las dos muestras, así que casi todos los pares dan ≈1.

Prueba con dos muestras que tienen los mismos tres taxones con abundancias distintas (el valor correcto es 0):

| Versión de SciPy | `pdist(abundancias, "jaccard")` | `pdist(abundancias > 0, "jaccard")` |
|---|---|---|
| 1.11.4 | 1.0 | 0.0 |
| 1.13.1 | 1.0 | 0.0 |
| 1.14.1 | 1.0 | 0.0 |
| 1.15.3 | 0.0 | 0.0 |
| 1.17.1 | 0.0 | 0.0 |

Qué revisar: si `diversidad_mod.py` pasa la tabla de abundancias directamente a `pdist(..., "jaccard")`, confirmar que la versión de SciPy del entorno es ≥ 1.15. Si es menor, hay que convertir la tabla a presencia/ausencia (`X > 0`) y volver a calcular la matriz de Jaccard y su heatmap. El script de recolección muestra la línea del código, la versión de SciPy y la mediana de la matriz de Jaccard obtenida; si la mediana está muy cerca de 1, es una señal de este problema. Bray-Curtis no se ve afectada.

## G.5 PCoA

| Eje | Valor propio | % de variación explicada |
|---|---|---|
| PC1 | ⏳ | ⏳ |
| PC2 | ⏳ | ⏳ |
| PC3 | ⏳ | ⏳ |

Coordenadas por muestra: `anexos/datos/diversidad/` ⏳.

## G.6 PERMANOVA

| Distancia | Grupos | N | Permutaciones | pseudo-F | p |
|---|---|---|---|---|---|
| Bray-Curtis | Sistema (3) | ⏳ | 999 | ⏳ | ⏳ |

Salida completa de scikit-bio: `anexos/datos/diversidad/` ⏳.
