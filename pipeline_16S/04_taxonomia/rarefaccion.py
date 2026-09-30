#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rarefaccion.py — Curva de rarefacción (HAC), con paso adaptativo por muestra
(cada curva usa un número de puntos proporcional a su propia profundidad,
para no simular mesetas falsas en muestras con pocas lecturas).
Líneas continuas, ID de finca al final de la curva y leyenda 'ID - Sistema'.
"""

import argparse
import sys
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def imprimir_error(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def cargar_mapeo_doble(ruta_puente: Path, ruta_meta: Path) -> dict:
    """
    Traduce Barcode -> {'id': ID de finca, 'sistema': Sistema Agrícola}.
    ruta_puente: CSV con columnas '#;ID;Barcode' (separador ';'), p. ej. Mapa_Barcodes_Microbioma.csv.
    ruta_meta:   Sistemas_Agrícolas_y_Muestras.xlsx (columnas 'ID', 'Sistema').
    """
    mapeo = {}
    try:
        # 1. Puente: Barcode -> ID (CSV, separador ; con respaldo a ,)
        df_puente = pd.read_csv(ruta_puente, sep=';')
        if not {"ID", "Barcode"}.issubset(df_puente.columns):
            df_puente = pd.read_csv(ruta_puente, sep=',')
        if not {"ID", "Barcode"}.issubset(df_puente.columns):
            raise SystemExit(
                f"ERROR: {ruta_puente} debe tener columnas 'ID' y 'Barcode'. "
                f"Columnas encontradas: {list(df_puente.columns)}"
            )
        barcode_a_id = {}
        for _, fila in df_puente.iterrows():
            try:
                num_bc = int(fila["Barcode"])
                barcode_a_id[num_bc] = str(fila["ID"]).strip()
            except (ValueError, TypeError):
                continue

        # 2. Metadatos: ID -> Sistema
        df_meta = pd.read_excel(ruta_meta, sheet_name=0, engine='openpyxl').dropna(subset=['ID'])
        id_a_sistema = dict(zip(df_meta['ID'].astype(str).str.strip(), df_meta['Sistema'].astype(str).str.strip()))

        # Construir diccionario: Barcode_Num -> {'id':..., 'sistema':...}
        for num_bc, finca_id in barcode_a_id.items():
            sistema = id_a_sistema.get(finca_id, "Desconocido")
            mapeo[num_bc] = {"id": finca_id, "sistema": sistema}

    except Exception as e:
        imprimir_error(f"[ERROR] Fallo leyendo puente/metadata: {e}")

    return mapeo

def obtener_info_muestra(nombre_columna: str, mapeo: dict):
    """
    Retorna:
    - etiqueta_final: el ID de la finca (p. ej. 'J95') para el texto al final de la curva.
    - etiqueta_leyenda: 'ID - Sistema' (p. ej. 'J95 - Empresarial') para la leyenda.
    Si no hay mapeo disponible para esa muestra, usa 'BCXX' como respaldo.
    """
    coincidencia = re.search(r'barcode0*(\d+)', nombre_columna.lower())
    num_bc = int(coincidencia.group(1)) if coincidencia else None
    str_bc = f"BC{num_bc:02d}" if num_bc is not None else nombre_columna

    info = mapeo.get(num_bc) if num_bc is not None else None
    if info:
        etiqueta_final = info["id"]
        etiqueta_leyenda = f"{info['id']} - {info['sistema']}"
    else:
        etiqueta_final = str_bc
        etiqueta_leyenda = str_bc

    return etiqueta_final, etiqueta_leyenda

def cargar_tabla_conteos(ruta_conteos: Path) -> pd.DataFrame:
    ext = ruta_conteos.suffix.lower()
    if ext == ".xlsx":
        df = pd.read_excel(ruta_conteos, index_col=0)
    else:
        df = pd.read_csv(ruta_conteos, sep="\t", index_col=0)
    if any("species" in str(x).lower() or "|" in str(x) for x in df.columns):
        df = df.T
    # Descartamos columnas de resumen agregadas por agrupar_counts_sistema.py
    # (Total_counts, Frecuencia_muestras, Frecuencia_<Sistema>...), que no son
    # muestras reales y arruinarían la rarefacción si se tratan como tal.
    columnas_resumen = [c for c in df.columns if str(c).startswith(("Total_counts", "Frecuencia_"))]
    if columnas_resumen:
        imprimir_error(f"[INFO] Ignorando columnas de resumen (no son muestras): {columnas_resumen}")
        df = df.drop(columns=columnas_resumen)
    return df

def rarefactar_muestra(vector_conteos: np.ndarray, n_puntos: int = 30, iteraciones: int = 10):
    """
    Calcula la curva de rarefacción con un número de puntos comparable
    entre muestras, sin importar su profundidad total (en vez de un
    tamaño de paso fijo, que le da muy pocos puntos a las muestras con pocas
    lecturas y hace que su curva parezca "meseta" solo por falta de
    resolución). 'iteraciones' es el número de submuestreos promediados
    en cada punto: más iteraciones = curva menos ruidosa y más confiable
    para juzgar si de verdad saturó o no.
    """
    lecturas_totales = int(vector_conteos.sum())
    if lecturas_totales == 0:
        return np.array([0]), np.array([0])

    lecturas_expandidas = np.repeat(np.arange(len(vector_conteos)), vector_conteos.astype(int))
    tamano_paso = max(1, lecturas_totales // n_puntos)
    profundidades = np.arange(1, lecturas_totales + 1, tamano_paso)
    if profundidades[-1] != lecturas_totales:
        profundidades = np.append(profundidades, lecturas_totales)

    riqueza_observada = []
    for d in profundidades:
        riqueza_sub = []
        for _ in range(iteraciones):
            muestreado = np.random.choice(lecturas_expandidas, size=d, replace=False)
            riqueza_sub.append(len(np.unique(muestreado)))
        riqueza_observada.append(np.mean(riqueza_sub))

    return profundidades, np.array(riqueza_observada)


def graficar_rarefaccion(df_conteos, mapeo, prefijo_salida, titulo, n_puntos=30, iteraciones=10, lecturas_minimas=1):
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(14, 8), dpi=300)

    # Estética de ejes (cuadro completo, como en la referencia)
    ax.tick_params(direction='in', length=6, width=1, labelsize=11)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    muestras = [col for col in df_conteos.columns if df_conteos[col].sum() >= lecturas_minimas]

    def extraer_num(x):
        m = re.search(r'barcode0*(\d+)', x.lower())
        return int(m.group(1)) if m else 999
    muestras.sort(key=extraer_num)

    if not muestras:
        return

    imprimir_error(f"[INFO] Trazando '{titulo}' para {len(muestras)} muestras...")

    # Forzamos una paleta de 20 colores distintos para máxima diferenciación
    cmap = cm.get_cmap('tab20')
    x_max = 0

    for i, col_muestra in enumerate(muestras):
        conteos = df_conteos[col_muestra].values
        conteos = conteos[conteos > 0]
        if len(conteos) == 0: continue

        profundidades, riqueza = rarefactar_muestra(conteos, n_puntos=n_puntos, iteraciones=iteraciones)
        etiqueta_final, etiqueta_leyenda = obtener_info_muestra(col_muestra, mapeo)
        x_max = max(x_max, profundidades[-1])

        color = cmap(i % 20)
        ax.plot(profundidades, riqueza, label=etiqueta_leyenda, color=color, linestyle='-', linewidth=1.6, alpha=0.85)

        # ID de la finca al final de la línea (ahora que las curvas terminan en
        # profundidades distintas, ya no se amontonan como cuando todas llegaban a 10000)
        ax.annotate(etiqueta_final,
                    xy=(profundidades[-1], riqueza[-1]),
                    xytext=(5, 0),
                    textcoords="offset points",
                    fontsize=8,
                    color=color,
                    fontweight='bold',
                    va='center')

    ax.set_xlim(left=0, right=x_max * 1.12)

    ax.set_xlabel("Tamaño de la Muestra (Número de Lecturas)", fontsize=13, labelpad=12)
    ax.set_ylabel("Riqueza Observada (Especies)", fontsize=13, labelpad=12)
    ax.set_title(titulo, fontsize=15, fontweight='bold', pad=15)

    # Leyenda a dos columnas fuera del gráfico (sin etiquetas flotantes sobre las curvas)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.72, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False, ncol=2, fontsize=9)

    png_salida = prefijo_salida.with_suffix(".png")
    plt.savefig(png_salida, bbox_inches="tight")
    plt.close()
    imprimir_error(f"[OK] Gráfica guardada: {png_salida.name}")

def main():
    base = Path("/home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/EMU_propio")

    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["hac", "sup"], default="hac",
                        help="Modelo de basecalling: define los valores por defecto de --counts, --outdir y el título.")
    parser.add_argument("--dataset", choices=["todo", "results"], default="todo",
                        help="'todo' -> EMU{model}_todo (todas las muestras); "
                             "'results' -> EMU{model}_results (solo muestras filtradas >500 lecturas).")
    parser.add_argument("--counts", type=Path, default=None,
                        help="Por defecto: EMU{model}_{dataset}/tabla_conteos_enriquecida.xlsx.")
    parser.add_argument("--outdir", type=Path, default=None,
                        help="Por defecto: EMU{model}_{dataset}/Rarefacción.")
    parser.add_argument("--meta", type=Path,
                        default=base / "Diversidad" / "Sistemas Agrícolas y Muestras.xlsx",
                        help="Maestro con columnas ID + Sistema.")
    parser.add_argument("--bridge", type=Path,
                        default=Path("/home/fenrir/Documentos/Tesis/datos_gulupa/Mapa Barcodes Microbioma.csv"),
                        help="CSV puente barcode->ID (columnas '#;ID;Barcode'). Es el mismo para HAC y SUP.")
    parser.add_argument("--n-points", type=int, default=30,
                        help="Número de puntos por curva (paso adaptativo: profundidad_total // n_points).")
    parser.add_argument("--iterations", type=int, default=10,
                        help="Submuestreos promediados por punto (más = curva más suave/confiable).")
    args = parser.parse_args()

    if args.counts is None:
        args.counts = base / f"EMU{args.model}_{args.dataset}" / "tabla_conteos_enriquecida.xlsx"
    if args.outdir is None:
        args.outdir = base / f"EMU{args.model}_{args.dataset}" / "Rarefacción"

    imprimir_error(f"[INFO] Modelo: {args.model} | dataset: {args.dataset} | counts: {args.counts} | outdir: {args.outdir}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    mapeo = cargar_mapeo_doble(args.bridge, args.meta)
    imprimir_error(f"[INFO] Muestras con ID/Sistema resueltos desde el puente: {len(mapeo)}")
    df_conteos = cargar_tabla_conteos(args.counts)

    sufijo = args.model.upper() if args.dataset == "todo" else f"{args.model.upper()}_{args.dataset}"
    # Una sola gráfica, con todas las muestras (paso adaptativo por muestra)
    graficar_rarefaccion(
        df_conteos, mapeo,
        args.outdir / f"Rarefaccion_{sufijo}",
        titulo=f"Curva de Rarefacción ({sufijo})",
        n_puntos=args.n_points, iteraciones=args.iterations, lecturas_minimas=1
    )

if __name__ == "__main__":
    main()
