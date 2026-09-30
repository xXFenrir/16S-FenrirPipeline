# Anexo A. Disponibilidad del código y los datos, y reproducibilidad

## A.1 Repositorio

- Dirección: https://github.com/xXFenrir/16S-FenrirPipeline
- Visibilidad al 30 de septiembre de 2026: **privado**. Para que el jurado y los lectores puedan consultarlo hay que hacerlo público. Lo ideal es además publicar una versión congelada (un *Release* de GitHub) y archivarla en Zenodo para obtener un DOI; la integración de GitHub con Zenodo requiere que el repositorio sea público. En la tesis se cita esa versión.

Estructura:

```
16S-FenrirPipeline/
├── README.md              Descripción del pipeline etapa por etapa
├── scripts/               Scripts finales (stats_fastq.py, dentrim_bam.py, EMU_propio.py,
│                          rebuild_counts.py, rarefaccion.py, diversidad_mod.py)
├── config/                primers.fasta y primers_config.txt
├── envs/                  Entornos conda exportados (.yml)
├── anexos/
│   ├── README.md          Índice de anexos y lista de verificación
│   ├── A_ … J_*.md        Anexos
│   ├── datos/             Tablas, figuras y reportes copiados del equipo de análisis
│   └── historico/         Versiones de desarrollo de los scripts (2025)
└── herramientas/
    └── recolectar_anexos.sh   Copia al repositorio los archivos de los anexos
```

## A.2 Datos

El repositorio **no** incluye los datos de secuenciación: POD5, BAM y FASTQ ocupan varios GB y GitHub rechaza archivos de más de 100 MB. El archivo `.gitignore` los excluye para evitar subirlos por error.

| Datos | Ubicación en el equipo de análisis | Dónde deberían quedar publicados |
|---|---|---|
| Señal cruda (POD5) | `datos_gulupa/20260715_1807_MN30942_FAZ24575_f8347d8d/pod5/` | Opcional (muy pesada) |
| Lecturas por barcode tras basecalling (BAM) | `datos_gulupa/data_gulupa_qs8/hac_8/` | NCBI SRA, como FASTQ por muestra |
| Lecturas limpias (FASTQ) | `datos_gulupa/data_gulupa_qs8/Limpieza/hac_8_trim_edlib/` | NCBI SRA o material suplementario |
| Base de datos de EMU | `Documentos/Tesis/emu_db/` | No hace falta: se descarga de OSF ([Anexo B](B_entorno_y_software.md#b4-base-de-datos-de-referencia)) |
| Tablas y figuras de resultados | `datos_gulupa/data_gulupa_qs8/EMU_propio/` | Este repositorio (`anexos/datos/`) |

Las rutas de la tabla son relativas a `/home/fenrir/Documentos/Tesis/`.

Texto sugerido para la sección de disponibilidad de datos de la tesis (completar los identificadores):

> Las lecturas de secuenciación del gen 16S rRNA están disponibles en NCBI Sequence Read Archive bajo el BioProject PRJNA_______. Los scripts, los entornos de software y los anexos con los parámetros y resultados complementarios están disponibles en https://github.com/xXFenrir/16S-FenrirPipeline (versión ______, DOI ______).

## A.3 Cómo reproducir el análisis

1. Clonar el repositorio:
   ```bash
   git clone https://github.com/xXFenrir/16S-FenrirPipeline.git
   cd 16S-FenrirPipeline
   ```
2. Crear los entornos conda a partir de `envs/`:
   ```bash
   conda env create -f envs/dentrim_env.yml
   conda env create -f envs/emu.yml
   ```
   Si falla por diferencias de sistema operativo, usar los archivos `*_historial.yml`, que solo listan los paquetes instalados explícitamente.
3. Descargar Dorado 2.1.0 y la base de datos de EMU ([Anexo B](B_entorno_y_software.md)).
4. Ejecutar los comandos del [Anexo C](C_protocolo_comandos_parametros.md) en orden, cambiando el prefijo `/home/fenrir/Documentos/Tesis` por la carpeta local.
5. Comparar las salidas con las tablas de `anexos/datos/`.

Fuentes de variación esperables al repetir el análisis:
- El basecalling con GPU puede dar diferencias mínimas entre equipos o versiones de controladores; con otra versión de Dorado o del modelo HAC los resultados cambian.
- La rarefacción usa submuestreo aleatorio; las curvas solo son idénticas si el script fija la semilla.
- PERMANOVA usa permutaciones; si no se fija la semilla, el valor p puede variar ligeramente entre ejecuciones.
