# Anexo D. Corrida de secuenciación, muestras y primers

Datos de la corrida de MinION, correspondencia entre barcodes y muestras, y primers usados. Los valores marcados con ⏳ se completan con los archivos que copia [`herramientas/recolectar_anexos.sh`](../herramientas/recolectar_anexos.sh) en `anexos/datos/corrida/` y `anexos/datos/metadatos/`.

## D.1 Corrida de secuenciación

| Dato | Valor | Fuente |
|---|---|---|
| Fecha de la corrida | 15 de julio de 2026 | Nombre de la carpeta `20260715_1807_MN30942_FAZ24575_f8347d8d` |
| Nombre del experimento en MinKNOW | `Microbioma_15072026` | Estructura de carpetas de MinKNOW |
| Secuenciador | MinION, ID `MN30942` | Nombre de la carpeta |
| Celda de flujo | ID `FAZ24575`; tipo ⏳ (el kit V14 corresponde a celdas R10.4.1, FLO-MIN114) | Nombre de la carpeta; campo `protocol` de `final_summary_*.txt` |
| Kit de preparación y barcodes | SQK-NBD114-96 (Native Barcoding Kit 96 V14) | Comando de Dorado |
| Identificador de la corrida | `f8347d8d…` | Nombre de la carpeta |
| Duración, lecturas, bases y N50 de la corrida | ⏳ | Reporte de MinKNOW (`report_*.html` / `.json`) |
| Basecalling en vivo durante la corrida | ⏳ (modelo y umbral de calidad) | `final_summary_*.txt` |

La misma corrida (`f8347d8d`) aparece también como `20260715_2307_MN30942_FAZ24575_f8347d8d` dentro de `sup_8/Microbioma_15072026/`. La diferencia de 5 horas coincide con el desfase de Colombia respecto a UTC (UTC−5), así que probablemente es la misma hora de inicio expresada en UTC. Confirmar con los campos `started` de ambos `final_summary_*.txt`.

## D.2 Muestras y diseño

Las muestras de suelo (rizobioma) de cultivos de gulupa se agrupan según el sistema de manejo agrícola de la finca de origen. En los metadatos la columna se llama `Sistema` y sus valores son:

| Sistema | Número de muestras |
|---|---|
| Empresarial | ⏳ |
| Campesina | ⏳ |
| Agroecológica | ⏳ |
| **Total** | ⏳ |

`diversidad_mod.py` traduce cada barcode a su sistema en dos pasos: barcode → ID de finca con el CSV puente `Mapa Barcodes Microbioma.csv`, e ID de finca → `Sistema` con el Excel de metadatos.

Tabla de correspondencia (se completa desde `anexos/datos/metadatos/`):

| Barcode | ID de finca | Sistema | Lecturas tras basecalling | Lecturas limpias | ¿Entra al análisis? |
|---|---|---|---|---|---|
| barcode01 | ⏳ | ⏳ | ⏳ | ⏳ | ⏳ |
| … | | | | | |

Criterios de exclusión aplicados por el pipeline:
- Muestras con menos de 500 lecturas limpias: `EMU_propio.py --min-reads-input 500` no las clasifica. Listado ⏳.
- Barcodes procesados: hasta `barcode96` en la corrida HAC (`--max-barcode 96`) y hasta `barcode73` en la prueba con SUP (`--max-barcode 73`). Número real de barcodes con muestra ⏳.

> **Privacidad.** Antes de hacer público el repositorio, revisar que los metadatos no incluyan nombres de propietarios, teléfonos ni coordenadas exactas de las fincas. Si los incluyen, publicar solo los códigos de finca y el sistema de manejo.

## D.3 Primers

| Nombre | Secuencia (5'→3') | Uso | Referencia |
|---|---|---|---|
| `16sF` | ⏳ | Primer directo | ⏳ |
| `16sR` | ⏳ | Primer reverso | ⏳ |

Fuente: `config/primers.fasta` (copia de `Limpieza/primers_gulupa/primers.fasta`).

Configuración de orientación para Pychopper (`config/primers_config.txt`). El formato que documenta Pychopper es:

```
+:16sF,-16sR|-:16sR,-16sF
```

Ver el [formato explicado](../config/README.md#formato-de-primers_configtxt) y la [verificación 2 del Anexo C](C_protocolo_comandos_parametros.md#verificaciones-pendientes-de-la-etapa-de-limpieza) antes de describir la orientación de las lecturas en el documento.

Los primers de la fase piloto (27F/1492R, datos públicos) están en el [Anexo H](H_fase_piloto_datos_publicos.md#h3-primers-y-procesamiento-del-estudio-original).

## D.4 Archivos que acompañan este anexo

Se generan con el script de recolección:

| Carpeta | Contenido |
|---|---|
| `anexos/datos/corrida/` | `final_summary_*.txt`, reportes de MinKNOW (`report_*`), `barcode_alignment_*.tsv`, `throughput_*.csv`, `pore_activity_*.csv`, y `lecturas_por_barcode_*.tsv` (resumen por barcode de `sequencing_summary.txt`) |
| `anexos/datos/metadatos/` | `Mapa Barcodes Microbioma.csv` y el Excel de metadatos |
| `config/` | `primers.fasta`, `primers_config.txt` |
