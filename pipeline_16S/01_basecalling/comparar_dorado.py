#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
comparar_dorado.py

Compara resultados de basecalling FAST, HAC y SUP de Dorado
a partir de los archivos sequencing_summary.txt separando por proyectos:
- Gulupa (Barcodes 01-73)
- Eliana (Barcodes 74-87)
"""

import argparse
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


###############################################################################
# FUNCIONES AUXILIARES
###############################################################################

def detectar_separador(csv_file):
    with open(csv_file, "r", encoding="utf-8") as f:
        primera = f.readline()
    if primera.count(";") > primera.count(","):
        return ";"
    return ","


def normalizar_barcode(valor):
    if pd.isna(valor):
        return None
    valor = str(valor).strip().lower()
    m = re.search(r'(\d+)', valor)
    if m is None:
        return None
    numero = int(m.group(1))
    return f"barcode{numero:02d}"


def formato_etiqueta_barcode(bc, muestra=""):
    """Convierte 'barcode01' a 'B01' y le anexa la muestra si existe."""
    if bc == "unclassified":
        return "Unclassified"
    m = re.search(r'(\d+)', str(bc))
    if m:
        num = int(m.group(1))
        etiqueta_corta = f"B{num:02d}"
    else:
        etiqueta_corta = str(bc)

    if muestra and muestra not in ["No asignada", "nan", "None", ""]:
        return f"{etiqueta_corta} - {muestra}"
    return etiqueta_corta


def detectar_columnas(df):
    posibles_barcodes = ["barcode", "barcodes", "barcode_id", "barcodeid", "codigo", "code"]
    posibles_muestra = ["sample", "sample_id", "sampleid", "muestra", "id", "nombre", "sample name"]

    barcode_col = None
    muestra_col = None

    for c in df.columns:
        nombre = c.lower().strip()
        if barcode_col is None:
            for p in posibles_barcodes:
                if p in nombre:
                    barcode_col = c
                    break
        if muestra_col is None:
            for p in posibles_muestra:
                if p in nombre:
                    muestra_col = c
                    break

    if barcode_col is None:
        raise ValueError("No pude identificar automáticamente la columna del Barcode.")
    if muestra_col is None:
        raise ValueError("No pude identificar automáticamente la columna de la muestra.")

    return barcode_col, muestra_col


def leer_mapa(csv_file):
    sep = detectar_separador(csv_file)
    mapa = pd.read_csv(csv_file, sep=sep)
    barcode_col, muestra_col = detectar_columnas(mapa)
    mapa[barcode_col] = mapa[barcode_col].apply(normalizar_barcode)

    diccionario = {}
    for _, fila in mapa.iterrows():
        bc = fila[barcode_col]
        muestra = fila[muestra_col]
        if pd.notna(bc):
            diccionario[bc] = str(muestra)

    return diccionario


def calcular_N50(longitudes):
    if len(longitudes) == 0:
        return 0
    datos = np.sort(np.asarray(longitudes))[::-1]
    mitad = datos.sum() / 2
    acumulado = np.cumsum(datos)
    indice = np.where(acumulado >= mitad)[0][0]
    return int(datos[indice])


def porcentaje(a, b):
    if b == 0:
        return 0
    return round(100 * a / b, 2)


def comprobar_archivo(ruta):
    if not os.path.isfile(ruta):
        raise FileNotFoundError(f"\nNo existe el archivo:\n{ruta}")


def obtener_grupo_barcode(bc):
    """Clasifica el barcode según el rango asignado."""
    if bc == "unclassified":
        return "unclassified"
    m = re.search(r'(\d+)', str(bc))
    if m:
        num = int(m.group(1))
        if 1 <= num <= 73:
            return "gulupa"
        elif 74 <= num <= 87:
            return "eliana"
    return "otros"


###############################################################################
# PROCESAMIENTO DE SUMMARY
###############################################################################

def procesar_summary(summary_file, metodo, mapa_muestras, chunksize=250000):
    print(f"\nLeyendo {metodo}...")

    resumen_global = {
        "Metodo": metodo,
        "Lecturas_totales": 0,
        "Pass_reads": 0,
        "Fail_reads": 0,
        "Bases_totales": 0,
        "Pass_bases": 0,
        "Qscore_sum": 0.0,
        "Longitudes": [],
        "min_start_time": float('inf'),
        "max_end_time": 0.0
    }

    barcodes = defaultdict(lambda: {
        "reads": 0,
        "pass": 0,
        "fail": 0,
        "bases": 0,
        "pass_bases": 0,
        "qscore_sum": 0.0,
        "longitudes": []
    })

    columnas = [
        "passes_filtering",
        "sequence_length_template",
        "mean_qscore_template",
        "barcode_arrangement",
        "alias",
        "start_time",
        "duration"
    ]

    for chunk in pd.read_csv(
            summary_file,
            sep="\t",
            usecols=lambda x: x in columnas,
            chunksize=chunksize,
            low_memory=False):

        chunk["sequence_length_template"] = pd.to_numeric(chunk["sequence_length_template"], errors="coerce").fillna(0)
        chunk["mean_qscore_template"] = pd.to_numeric(chunk["mean_qscore_template"], errors="coerce").fillna(0)
        chunk["passes_filtering"] = chunk["passes_filtering"].astype(str).str.upper()

        if "barcode_arrangement" in chunk.columns:
            chunk["barcode"] = chunk["barcode_arrangement"]
        elif "alias" in chunk.columns:
            chunk["barcode"] = chunk["alias"]
        else:
            chunk["barcode"] = "unknown"

        chunk["barcode"] = chunk["barcode"].fillna("unknown").astype(str).str.lower()

        es_pass = chunk["passes_filtering"] == "TRUE"
        es_fail = chunk["passes_filtering"] == "FALSE"

        resumen_global["Lecturas_totales"] += len(chunk)
        resumen_global["Pass_reads"] += es_pass.sum()
        resumen_global["Fail_reads"] += es_fail.sum()

        resumen_global["Bases_totales"] += chunk["sequence_length_template"].sum()
        resumen_global["Pass_bases"] += chunk.loc[es_pass, "sequence_length_template"].sum()

        resumen_global["Qscore_sum"] += chunk["mean_qscore_template"].sum()
        resumen_global["Longitudes"].extend(chunk["sequence_length_template"].tolist())

        if "start_time" in chunk.columns:
            chunk_start = pd.to_numeric(chunk["start_time"], errors="coerce")
            if not chunk_start.empty and chunk_start.dropna().size > 0:
                resumen_global["min_start_time"] = min(resumen_global["min_start_time"], chunk_start.min())
                if "duration" in chunk.columns:
                    chunk_dur = pd.to_numeric(chunk["duration"], errors="coerce").fillna(0)
                    chunk_end = chunk_start + chunk_dur
                    resumen_global["max_end_time"] = max(resumen_global["max_end_time"], chunk_end.max())
                else:
                    resumen_global["max_end_time"] = max(resumen_global["max_end_time"], chunk_start.max())

        for bc, grupo in chunk.groupby("barcode"):
            nombre = normalizar_barcode(bc) or "unclassified"
            datos = barcodes[nombre]

            g_pass = grupo["passes_filtering"] == "TRUE"

            datos["reads"] += len(grupo)
            datos["pass"] += g_pass.sum()
            datos["fail"] += (grupo["passes_filtering"] == "FALSE").sum()
            datos["bases"] += grupo["sequence_length_template"].sum()
            datos["pass_bases"] += grupo.loc[g_pass, "sequence_length_template"].sum()
            datos["qscore_sum"] += grupo["mean_qscore_template"].sum()
            datos["longitudes"].extend(grupo["sequence_length_template"].tolist())

    def estructurar_resumen(res_dict):
        res_dict["% Pass"] = porcentaje(res_dict["Pass_reads"], res_dict["Lecturas_totales"])
        res_dict["QScore_medio"] = round(res_dict["Qscore_sum"] / res_dict["Lecturas_totales"], 2) if res_dict["Lecturas_totales"] > 0 else 0
        res_dict["Longitud_media"] = round(np.mean(res_dict["Longitudes"]), 2) if res_dict["Longitudes"] else 0
        res_dict["N50"] = calcular_N50(res_dict["Longitudes"])

        if res_dict["min_start_time"] != float('inf') and res_dict["max_end_time"] > 0:
            duracion_segundos = res_dict["max_end_time"] - res_dict["min_start_time"]
            horas = int(duracion_segundos // 3600)
            minutos = int((duracion_segundos % 3600) // 60)
            segundos = int(duracion_segundos % 60)
            res_dict["Tiempo_corrida"] = f"{horas:02d}h {minutos:02d}m {segundos:02d}s"
        else:
            res_dict["Tiempo_corrida"] = "N/A"
        return res_dict

    resumen_global = estructurar_resumen(resumen_global)

    def calcular_sub_resumen(grupo_objetivo):
        sub = {
            "Metodo": metodo, "Lecturas_totales": 0, "Pass_reads": 0, "Fail_reads": 0,
            "Bases_totales": 0, "Pass_bases": 0, "Qscore_sum": 0.0, "Longitudes": [],
            "min_start_time": resumen_global["min_start_time"], "max_end_time": resumen_global["max_end_time"]
        }
        for bc, datos in barcodes.items():
            grp = obtener_grupo_barcode(bc)
            if grp == grupo_objetivo or grp == "unclassified":
                sub["Lecturas_totales"] += datos["reads"]
                sub["Pass_reads"] += datos["pass"]
                sub["Fail_reads"] += datos["fail"]
                sub["Bases_totales"] += datos["bases"]
                sub["Pass_bases"] += datos["pass_bases"]
                sub["Qscore_sum"] += datos["qscore_sum"]
                sub["Longitudes"].extend(datos["longitudes"])
        return estructurar_resumen(sub)

    resumen_gulupa = calcular_sub_resumen("gulupa")
    resumen_eliana = calcular_sub_resumen("eliana")

    dict_barcodes = {}
    todos_barcodes = sorted(set(list(barcodes.keys()) + list(mapa_muestras.keys())))

    for bc in todos_barcodes:
        datos = barcodes[bc]
        muestra = mapa_muestras.get(bc, "No asignada" if bc == "unclassified" else "")

        dict_barcodes[bc] = {
            "Muestra": muestra,
            "Lecturas": datos["reads"],
            "PASS": datos["pass"],
            "FAIL": datos["fail"],
            "%PASS": porcentaje(datos["pass"], datos["reads"]),
            "Bases": datos["bases"],
            "Longitud_media": round(np.mean(datos["longitudes"]), 2) if datos["reads"] else 0,
            "Longitud_mediana": round(np.median(datos["longitudes"]), 2) if datos["reads"] else 0,
            "N50": calcular_N50(datos["longitudes"]),
            "Qscore": round(datos["qscore_sum"] / datos["reads"], 2) if datos["reads"] else 0
        }

    return resumen_global, resumen_gulupa, resumen_eliana, dict_barcodes


def construir_tabla_barcodes_filtrada(dict_fast, dict_hac, dict_sup, mapa_muestras, grupo_filtro, incluir_fast=True):
    modelos = [("HAC", dict_hac), ("SUP", dict_sup)]
    if incluir_fast and dict_fast:
        modelos.insert(0, ("FAST", dict_fast))

    todos_barcodes = sorted(set([bc for _, d in modelos for bc in d.keys()]))
    barcodes_filtrados = [bc for bc in todos_barcodes if obtener_grupo_barcode(bc) == grupo_filtro or bc == "unclassified"]

    metricas = [
        "Lecturas", "PASS", "FAIL", "%PASS", "Bases",
        "Longitud_media", "Longitud_mediana", "N50", "Qscore"
    ]

    columns = pd.MultiIndex.from_tuples(
        [("", "Barcode"), ("", "Muestra")] +
        [(nombre_m, metrica) for nombre_m, _ in modelos for metrica in metricas]
    )

    filas = []
    for bc in barcodes_filtrados:
        muestra = mapa_muestras.get(bc, "No asignada" if bc == "unclassified" else "")
        fila = [bc, muestra]

        for _, d in modelos:
            data_metodo = d.get(bc, {})
            for m in metricas:
                fila.append(data_metodo.get(m, 0))

        filas.append(fila)

    return pd.DataFrame(filas, columns=columns)


###############################################################################
# GENERACIÓN DE GRÁFICOS (SEABORN - VERTICAL)
###############################################################################

def generar_graficos_por_grupo(dict_fast, dict_hac, dict_sup, mapa_muestras, out_dir, grupo_nombre, titulo_grupo, incluir_fast=True):
    print(f"\nGenerando gráficos verticales para {titulo_grupo}...")

    registros = []
    metodos = [("HAC", dict_hac), ("SUP", dict_sup)]
    if incluir_fast and dict_fast:
        metodos.insert(0, ("FAST", dict_fast))

    for metodo_nombre, d in metodos:
        for bc, datos in d.items():
            if bc == "unclassified" or obtener_grupo_barcode(bc) != grupo_nombre:
                continue

            muestra = mapa_muestras.get(bc, "")
            etiqueta = formato_etiqueta_barcode(bc, muestra)

            registros.append({
                "Barcode": bc,
                "Etiqueta": etiqueta,
                "Metodo": metodo_nombre,
                "PASS": datos.get("PASS", 0),
                "Qscore": datos.get("Qscore", 0)
            })

    df_long = pd.DataFrame(registros)

    if df_long.empty:
        print(f"No hay datos suficientes para graficar {titulo_grupo}.")
        return

    # Ordenar barcodes
    df_long = df_long.sort_values(by="Barcode", ascending=True)

    sns.set_theme(style="whitegrid", font_scale=1.0)
    
    paleta_personalizada = {
        "FAST": "#E63946",  # Rojo
        "HAC": "#1D3557",   # Azul
        "SUP": "#2A9D8F"    # Verde
    }

    # Calcular altura dinámica según la cantidad de muestras para que se vea claro
    num_muestras = df_long["Etiqueta"].nunique()
    altura_fig = max(8, num_muestras * 0.35)

    suffix = "" if incluir_fast else "_sin_fast"

    # 1. Gráfico Lecturas PASS (Horizontal/Vertical)
    plt.figure(figsize=(12, altura_fig))
    ax1 = sns.barplot(
        data=df_long, y="Etiqueta", x="PASS", hue="Metodo",
        palette=paleta_personalizada, edgecolor="black", linewidth=0.5
    )
    plt.title(f"Lecturas PASS por Muestra - {titulo_grupo}", fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("Cantidad de Lecturas PASS", fontsize=12, labelpad=10)
    plt.ylabel("Barcode / Muestra", fontsize=12, labelpad=10)
    ax1.xaxis.set_major_formatter('{x:,.0f}')
    plt.legend(title="Método", frameon=True, facecolor="white", loc="lower right")
    plt.tight_layout()

    fig1_path = os.path.join(out_dir, f"{titulo_grupo}_lecturas_pass{suffix}.png")
    plt.savefig(fig1_path, dpi=300)
    plt.close()

    # 2. Gráfico QScore (Horizontal/Vertical)
    plt.figure(figsize=(12, altura_fig))
    ax2 = sns.barplot(
        data=df_long, y="Etiqueta", x="Qscore", hue="Metodo",
        palette=paleta_personalizada, edgecolor="black", linewidth=0.5
    )
    plt.title(f"Calidad Media (QScore) por Muestra - {titulo_grupo}", fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("QScore (Phred)", fontsize=12, labelpad=10)
    plt.ylabel("Barcode / Muestra", fontsize=12, labelpad=10)
    plt.axvline(10, color="gray", linestyle="--", alpha=0.7, label="Umbral Q10")

    q_min = max(0, df_long["Qscore"].min() - 1)
    q_max = df_long["Qscore"].max() + 1.5
    plt.xlim(q_min, q_max)

    plt.legend(title="Método", frameon=True, facecolor="white", loc="lower right")
    plt.tight_layout()

    fig2_path = os.path.join(out_dir, f"{titulo_grupo}_qscore{suffix}.png")
    plt.savefig(fig2_path, dpi=300)
    plt.close()

    print(f"Gráficos de {titulo_grupo} guardados:\n  - {fig1_path}\n  - {fig2_path}")


###############################################################################
# PROGRAMA PRINCIPAL
###############################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Comparación FAST vs HAC vs SUP (Dorado)")
    parser.add_argument("--fast", required=False, default=None, help="sequencing_summary FAST (opcional)")
    parser.add_argument("--hac", required=True, help="sequencing_summary HAC")
    parser.add_argument("--sup", required=True, help="sequencing_summary SUP")
    parser.add_argument("--mapa", required=True, help="CSV Barcode -> Muestra")
    parser.add_argument("--out", required=True, help="Directorio de salida")
    parser.add_argument("--sin-fast", action="store_true", help="Omitir el modelo FAST en tablas y gráficos")

    args = parser.parse_args()

    if not args.sin_fast and args.fast:
        comprobar_archivo(args.fast)
    
    comprobar_archivo(args.hac)
    comprobar_archivo(args.sup)
    comprobar_archivo(args.mapa)

    os.makedirs(args.out, exist_ok=True)

    MAPA_MUESTRAS = leer_mapa(args.mapa)
    print(f"Se encontraron {len(MAPA_MUESTRAS)} muestras en el mapa.")

    incluir_fast = (not args.sin_fast) and (args.fast is not None)

    # 1. Procesar summaries
    dict_fast, g_fast, gulupa_fast, eliana_fast = None, None, None, None
    if incluir_fast:
        g_fast, gulupa_fast, eliana_fast, dict_fast = procesar_summary(args.fast, "FAST", MAPA_MUESTRAS)
    
    g_hac,  gulupa_hac,  eliana_hac,  dict_hac  = procesar_summary(args.hac, "HAC", MAPA_MUESTRAS)
    g_sup,  gulupa_sup,  eliana_sup,  dict_sup  = procesar_summary(args.sup, "SUP", MAPA_MUESTRAS)

    # 2. Formatear DataFrames de Resumen
    columnas_resumen = {
        "Metodo": "Método", "Lecturas_totales": "Lecturas totales", "Pass_reads": "Pass reads",
        "Fail_reads": "Fail reads", "% Pass": "% Pass", "Bases_totales": "Bases totales",
        "Pass_bases": "Pass bases", "Longitud_media": "Longitud media", "N50": "N50",
        "QScore_medio": "QScore medio", "Tiempo_corrida": "Tiempo de corrida"
    }

    lista_global = [g_fast, g_hac, g_sup] if incluir_fast else [g_hac, g_sup]
    lista_gulupa = [gulupa_fast, gulupa_hac, gulupa_sup] if incluir_fast else [gulupa_hac, gulupa_sup]
    lista_eliana = [eliana_fast, eliana_hac, eliana_sup] if incluir_fast else [eliana_hac, eliana_sup]

    df_res_global = pd.DataFrame(lista_global)[list(columnas_resumen.keys())].rename(columns=columnas_resumen)
    df_res_gulupa = pd.DataFrame(lista_gulupa)[list(columnas_resumen.keys())].rename(columns=columnas_resumen)
    df_res_eliana = pd.DataFrame(lista_eliana)[list(columnas_resumen.keys())].rename(columns=columnas_resumen)

    # 3. Tablas detalladas filtradas
    df_gulupa_barcodes = construir_tabla_barcodes_filtrada(dict_fast, dict_hac, dict_sup, MAPA_MUESTRAS, "gulupa", incluir_fast=incluir_fast)
    df_eliana_barcodes = construir_tabla_barcodes_filtrada(dict_fast, dict_hac, dict_sup, MAPA_MUESTRAS, "eliana", incluir_fast=incluir_fast)

    # 4. Generar Gráficos independientes
    generar_graficos_por_grupo(dict_fast, dict_hac, dict_sup, MAPA_MUESTRAS, args.out, "gulupa", "Gulupa", incluir_fast=incluir_fast)
    generar_graficos_por_grupo(dict_fast, dict_hac, dict_sup, MAPA_MUESTRAS, args.out, "eliana", "Eliana", incluir_fast=incluir_fast)

    # 5. Exportar a Excel
    nombre_excel = "Comparacion_Dorado_sin_FAST.xlsx" if args.sin_fast else "Comparacion_Dorado.xlsx"
    archivo = os.path.join(args.out, nombre_excel)

    with pd.ExcelWriter(archivo, engine="openpyxl") as writer:
        df_res_global.to_excel(writer, sheet_name="Resumen_Global", index=False)
        df_res_gulupa.to_excel(writer, sheet_name="Resumen_Gulupa", index=False)
        df_res_eliana.to_excel(writer, sheet_name="Resumen_Eliana", index=False)
        df_gulupa_barcodes.to_excel(writer, sheet_name="Muestras_Gulupa")
        df_eliana_barcodes.to_excel(writer, sheet_name="Muestras_Eliana")

    print(f"\n¡Proceso completado con éxito! Archivos guardados en:\n{args.out}")
