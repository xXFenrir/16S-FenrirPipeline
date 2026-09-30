# Anexo B. Entorno computacional y versiones de software

Este anexo reúne lo necesario para repetir el análisis con el mismo software: el equipo, los entornos conda, la versión de cada herramienta y la base de datos de referencia.

Los valores marcados con ⏳ se completan con los archivos que genera [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh) en `anexos/datos/entorno/`.

## B.1 Equipo y sistema operativo

| Elemento | Valor | Fuente |
|---|---|---|
| Equipo | Portátil MSI Vector 16 HX AI (A2XWHG) | Nombre del equipo en la terminal (`fenrir-Vector-16-HX-AI-A2XWHG`); confirmar en `sistema.txt` |
| Procesador | ⏳ | `sistema.txt` |
| Memoria RAM | ⏳ | `sistema.txt` |
| GPU (usada por Dorado) | ⏳ | `sistema.txt` |
| Sistema operativo | Linux, distribución ⏳ | `sistema.txt` |
| Gestor de entornos | Anaconda (`~/anaconda3`) en 2025; Miniconda (`~/miniconda3`) en 2026 | Rutas en los comandos documentados |
| Hilos usados | 20 en limpieza y taxonomía (8 en la fase piloto) | Comandos del [Anexo C](C_protocolo_comandos_parametros.md) |

## B.2 Entornos conda

| Entorno | Herramientas | Python | Archivo |
|---|---|---|---|
| `dentrim_env` | samtools, Pychopper, edlib, Filtlong | 3.9 | `envs/dentrim_env.yml` ⏳ |
| `emu` | EMU, minimap2 | ⏳ | `envs/emu.yml` ⏳ |
| ⏳ | pandas, NumPy, SciPy, scikit-bio, seaborn, matplotlib (scripts de taxonomía y diversidad) | ⏳ | ⏳ |
| `pipelinefenrir` (solo fase piloto) | EMU y scripts de desarrollo | ⏳ | `envs/pipelinefenrir.yml` ⏳ |

Instalación documentada durante el desarrollo:

```bash
# Limpieza (entorno dentrim_env)
conda install -c nanoporetech -c conda-forge -c bioconda "nanoporetech::pychopper"
conda install bioconda::filtlong

# Taxonomía
conda create -n emu -c conda-forge -c bioconda emu
```

## B.3 Software

| Herramienta | Versión | Uso en el pipeline | Instalación | Referencia |
|---|---|---|---|---|
| Dorado | **2.1.0** (1.1.1 en la fase piloto) | Basecalling (modelo HAC) y demultiplexing | Binario precompilado de ONT | Oxford Nanopore Technologies. Dorado. https://github.com/nanoporetech/dorado |
| pod5 | ⏳ | Conversión FAST5 → POD5 (solo si hace falta) | ⏳ | Oxford Nanopore Technologies. https://github.com/nanoporetech/pod5-file-format |
| samtools | ⏳ | BAM → FASTQ | bioconda | Danecek et al., 2021 |
| Pychopper | ⏳ | Detección de primers, orientación y recorte | canal `nanoporetech` | Oxford Nanopore Technologies / EPI2ME Labs. https://github.com/epi2me-labs/pychopper |
| edlib | ⏳ | Alineamiento de primers (backend de Pychopper) | conda | Šošić y Šikić, 2017 |
| Filtlong | ⏳ | Filtro por longitud y calidad | bioconda | Wick, R. Filtlong. https://github.com/rrwick/Filtlong |
| EMU | ⏳ | Clasificación taxonómica a nivel de especie | bioconda | Curry et al., 2022 |
| minimap2 | ⏳ | Alineamiento de lecturas (lo usa EMU) | bioconda | Li, 2018 |
| Python | 3.9 en `dentrim_env`; ⏳ en los demás | Scripts propios | conda | — |
| pandas | ⏳ | Tablas | conda/pip | McKinney, 2010 |
| NumPy | ⏳ | Cálculo numérico | conda/pip | Harris et al., 2020 |
| SciPy | ⏳ | Mann-Whitney U, distancias (`pdist`) | conda/pip | Virtanen et al., 2020 |
| scikit-bio | ⏳ | Diversidad alfa, PCoA, PERMANOVA | conda/pip | The scikit-bio development team. https://scikit.bio |
| matplotlib | ⏳ | Figuras | conda/pip | Hunter, 2007 |
| seaborn | ⏳ | Boxplots y heatmaps | conda/pip | Waskom, 2021 |
| openpyxl | ⏳ | Exportar XLSX (`stats_fastq.py`) | conda/pip | https://openpyxl.readthedocs.io |

El script de recolección genera `anexos/datos/entorno/versiones_software.tsv` con la versión exacta de cada paquete en cada entorno (a partir de `conda list`) y la salida de `dorado --version`.

## B.4 Base de datos de referencia

Se usó la **base de datos por defecto de EMU**: combina rrnDB v5.6 y NCBI 16S RefSeq del 17 de septiembre de 2020, con la taxonomía de NCBI de la misma fecha, y contiene 49 301 secuencias de 17 555 especies de bacterias y arqueas (documentación de EMU).

Descarga documentada en octubre de 2025:

```bash
export EMU_DATABASE_DIR="/home/fenrir/emu_db"
mkdir -p "$EMU_DATABASE_DIR"
cd "$EMU_DATABASE_DIR"
conda install -c conda-forge osfclient
osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
tar -xvf emu.tar
```

En los comandos finales la base está en `/home/fenrir/Documentos/Tesis/emu_db`. El script de recolección guarda la lista de archivos y la suma MD5 de `species_taxid.fasta` y `taxonomy.tsv` en `anexos/datos/entorno/emu_db_info.txt`, para dejar constancia de la versión exacta.

La base de datos de EMU cambió de formato en la versión 3.0 de EMU. La versión instalada debe corresponder con la base descargada ([documentación de EMU](https://github.com/treangenlab/emu#1-download-database)).

Valores por defecto de `emu abundance` según la documentación de EMU v3.6.2 (confirmar con la versión instalada): `--type map-ont`, `--min-abundance 0.0001`, `--N 50` alineamientos por lectura, `--K 500M`, `--max-align-len 2000`.

## B.5 Referencias

- Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32–46.
- Curry, K. D., Wang, Q., Nute, M. G., et al. (2022). Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data. *Nature Methods*, 19, 845–853. https://doi.org/10.1038/s41592-022-01520-4
- Danecek, P., Bonfield, J. K., Liddle, J., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008.
- Harris, C. R., Millman, K. J., van der Walt, S. J., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362.
- Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90–95.
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100.
- McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61.
- O'Leary, N. A., et al. (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Research*, 44(D1), D733–D745.
- Šošić, M., y Šikić, M. (2017). Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. *Bioinformatics*, 33(9), 1394–1395.
- Stoddard, S. F., Smith, B. J., Hein, R., Roller, B. R. K., y Schmidt, T. M. (2015). rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development. *Nucleic Acids Research*, 43(D1), D593–D598.
- Virtanen, P., Gommers, R., Oliphant, T. E., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.
- Waskom, M. L. (2021). seaborn: statistical data visualization. *Journal of Open Source Software*, 6(60), 3021.

Ajustar al estilo de citación de la tesis. Las herramientas sin artículo (Dorado, Pychopper, Filtlong, scikit-bio) se citan por su repositorio y versión.
