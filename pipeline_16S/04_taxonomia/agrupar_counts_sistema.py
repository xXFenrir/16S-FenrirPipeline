#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
agrupar_counts_sistema.py — Enriquece la tabla de conteos (tabla_conteos)
con columnas de Total y Frecuencia por taxón, y genera además una tabla agrupada
por Sistema (sumando conteos de las muestras que pertenecen a cada sistema).

Entradas esperadas:
  - Tabla de conteos: filas = _feature_id (taxón), columnas = SampleID (muestras).
    Acepta .tsv o .xlsx.
  - Metadata, en uno de dos modos (igual que en rarefaccion.py / diversidad_mod.py):
      MODO DIRECTO: un solo archivo .xlsx con columnas SampleID y Sistema
      (p. ej. la hoja 'MAPEO_MUESTRAS' de METADATA_GULUPA.xlsx). Se usa pasando
      solo --metadata.
      MODO DOBLE (puente + maestro): --bridge es un .xlsx con columnas
      'barcode' e 'ID' (traduce barcodeXX -> ID de finca) y --metadata es el
      maestro con columnas 'ID' y 'Sistema' (p. ej. Sistemas Agrícolas y Muestras.xlsx).
      Se activa automáticamente si se pasa --bridge.

Salida (en --outdir):
  - tabla_conteos_enriquecida.tsv / .xlsx
      Tabla original (feature_id x muestras) + columnas nuevas:
        * Total_counts          -> suma de conteos del taxón en todas las muestras.
        * Frecuencia_muestras   -> en cuántas muestras aparece (conteo > 0), en total.
        * Frecuencia_<Sistema>  -> una columna por cada sistema (p. ej.
          Frecuencia_Agroecológica, Frecuencia_Campesina, Frecuencia_Empresarial):
          en cuántas muestras DE ESE sistema aparece el taxón.
      Las muestras sin Sistema asignado en la metadata se agrupan aparte en
      'Frecuencia_Sin_metadata', para que ninguna quede oculta.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd


def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def load_counts(path: Path) -> pd.DataFrame:
    """Carga la tabla de conteos (feature_id x muestras) desde .tsv o .xlsx."""
    ext = path.suffix.lower()
    if ext == ".xlsx":
        df = pd.read_excel(path, index_col=0)
    else:
        df = pd.read_csv(path, sep="\t", index_col=0)
    # Nos aseguramos de trabajar solo con columnas numéricas (por si ya
    # traía columnas de Total/Frecuencia de una corrida anterior).
    df = df.loc[:, df.apply(pd.api.types.is_numeric_dtype)]
    return df


def load_metadata_mapping(path: Path, sheet, sample_col: str, group_col: str) -> pd.Series:
    """Lee la metadata y devuelve un mapeo SampleID -> Grupo (Sistema)."""
    df_meta = pd.read_excel(path, sheet_name=sheet)
    if sample_col not in df_meta.columns:
        raise SystemExit(
            f"ERROR: la columna '{sample_col}' no existe en la hoja de metadata. "
            f"Columnas disponibles: {list(df_meta.columns)}"
        )
    if group_col not in df_meta.columns:
        raise SystemExit(
            f"ERROR: la columna '{group_col}' no existe en la hoja de metadata. "
            f"Columnas disponibles: {list(df_meta.columns)}"
        )
    df_meta = df_meta.dropna(subset=[sample_col, group_col])
    mapping = dict(
        zip(
            df_meta[sample_col].astype(str).str.strip(),
            df_meta[group_col].astype(str).str.strip(),
        )
    )
    return mapping


def load_double_mapping(bridge_path: Path, bridge_sheet, meta_path: Path, sample_columns) -> dict:
    """
    Construye SampleID -> Sistema usando dos archivos, igual que en
    rarefaccion.py / diversidad_mod.py:
      - bridge_path: .xlsx con columnas 'barcode' e 'ID' (traduce barcodeXX -> ID de finca).
      - meta_path:   .xlsx (maestro) con columnas 'ID' y 'Sistema'.
    """
    df_bridge = pd.read_excel(bridge_path, sheet_name=bridge_sheet)
    if "barcode" not in df_bridge.columns or "ID" not in df_bridge.columns:
        raise SystemExit(
            f"ERROR: el archivo puente debe tener columnas 'barcode' e 'ID'. "
            f"Columnas encontradas: {list(df_bridge.columns)}"
        )
    # Filtramos filas basura (resúmenes, totales, etc. que a veces quedan
    # mezclados en la misma hoja): solo nos quedamos con barcodes con formato válido.
    valid = df_bridge["barcode"].astype(str).str.match(r"^barcode\d+$", case=False, na=False)
    df_bridge = df_bridge.loc[valid]
    barcode_to_id = dict(
        zip(df_bridge["barcode"].astype(str).str.strip(), df_bridge["ID"].astype(str).str.strip())
    )

    df_meta = pd.read_excel(meta_path)
    if "ID" not in df_meta.columns or "Sistema" not in df_meta.columns:
        raise SystemExit(
            f"ERROR: el archivo maestro debe tener columnas 'ID' y 'Sistema'. "
            f"Columnas encontradas: {list(df_meta.columns)}"
        )
    df_meta = df_meta.dropna(subset=["ID"])
    id_to_sistema = dict(
        zip(
            df_meta["ID"].astype(str).str.strip(),
            df_meta["Sistema"].astype(str).str.strip(),
        )
    )

    mapping = {}
    sin_bridge = []
    for sample in sample_columns:
        m = re.search(r"barcode0*(\d+)", str(sample), re.IGNORECASE)
        if not m:
            continue
        bc_num = int(m.group(1))
        candidatos = [f"barcode{bc_num:02d}", f"barcode{bc_num}", f"BC{bc_num:02d}"]
        finca_id = next((barcode_to_id[c] for c in candidatos if c in barcode_to_id), None)
        if finca_id is None:
            sin_bridge.append(sample)
            continue
        sistema = id_to_sistema.get(finca_id)
        if sistema:
            mapping[sample] = sistema
        else:
            eprint(f"[WARN] '{sample}' -> ID '{finca_id}' no tiene Sistema en el maestro.")

    if sin_bridge:
        eprint(f"[WARN] {len(sin_bridge)} muestra(s) no encontradas en el archivo puente: {sin_bridge}")

    return mapping


def add_total_and_frequency(counts: pd.DataFrame) -> pd.DataFrame:
    """Agrega columnas Total_counts y Frecuencia_muestras a la tabla por muestra."""
    out = counts.copy()
    out["Total_counts"] = counts.sum(axis=1)
    out["Frecuencia_muestras"] = (counts > 0).sum(axis=1)
    out = out.sort_values("Total_counts", ascending=False)
    return out


def add_system_frequencies(counts: pd.DataFrame, mapping: dict) -> pd.DataFrame:
    """
    Para cada sistema presente en 'mapping', agrega una columna
    'Frecuencia_<Sistema>' = en cuántas muestras de ESE sistema aparece el
    taxón (conteo > 0). Las muestras sin sistema asignado no se pierden:
    se reportan en el log y se agrupan bajo 'Frecuencia_Sin_metadata'.
    Devuelve solo las columnas de frecuencia nuevas (para luego concatenarlas
    a la tabla principal).
    """
    sample_to_group = {s: mapping.get(s, "Sin_metadata") for s in counts.columns}
    sin_meta = [s for s, g in sample_to_group.items() if g == "Sin_metadata"]
    if sin_meta:
        eprint(
            f"[WARN] {len(sin_meta)} muestra(s) sin Sistema asignado en la metadata "
            f"(no cuentan en las frecuencias por sistema, se agrupan aparte): {sin_meta}"
        )

    presence = counts > 0
    systems = sorted(set(g for g in sample_to_group.values() if g != "Sin_metadata"))

    freq_cols = {}
    for sys_name in systems:
        cols_sys = [s for s, g in sample_to_group.items() if g == sys_name]
        freq_cols[f"Frecuencia_{sys_name}"] = presence[cols_sys].sum(axis=1)

    return pd.DataFrame(freq_cols, index=counts.index)


def write_tsv_and_xlsx(df: pd.DataFrame, outdir: Path, basename: str):
    tsv_path = outdir / f"{basename}.tsv"
    xlsx_path = outdir / f"{basename}.xlsx"
    df.to_csv(tsv_path, sep="\t")
    n_rows, n_cols = df.shape
    if n_rows > 1_048_575 or n_cols > 16_383:
        eprint(f"[WARN] {xlsx_path.name}: {n_rows}x{n_cols} excede el límite de Excel; se omite el .xlsx.")
    else:
        df.to_excel(xlsx_path, engine="openpyxl")
    eprint(f"[OK] Escrito: {tsv_path}")
    if xlsx_path.exists():
        eprint(f"[OK] Escrito: {xlsx_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Agrega Total/Frecuencia a la tabla de conteos y la agrupa por Sistema."
    )
    base = Path("/home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8/EMU_propio")

    ap.add_argument("--model", choices=["hac", "sup"], default="hac",
                     help="Modelo de basecalling: define los valores por defecto de --counts, --outdir y --bridge-sheet.")
    ap.add_argument("--dataset", choices=["todo", "results"], default="todo",
                     help="'todo' -> EMU{model}_todo (todas las muestras); "
                          "'results' -> EMU{model}_results (solo muestras filtradas >500 lecturas).")
    ap.add_argument("--counts", type=Path, default=None,
                     help="Tabla de conteos (.tsv o .xlsx). Por defecto: EMU{model}_{dataset}/tabla_conteos.tsv.")
    ap.add_argument("--metadata", type=Path,
                     default=base / "Diversidad" / "Sistemas Agrícolas y Muestras.xlsx",
                     help="Maestro con columnas ID + Sistema (Sistemas Agrícolas y Muestras.xlsx).")
    ap.add_argument("--bridge", type=Path,
                     default=base / "Diversidad" / "Muestras_raref.xlsx",
                     help="Archivo puente barcode->ID (Muestras_raref.xlsx). Pasa --no-bridge para usar --metadata en modo directo.")
    ap.add_argument("--bridge-sheet", default=None,
                     help="Hoja de --bridge a usar. Por defecto: 'HAC' o 'SUP' según --model.")
    ap.add_argument("--no-bridge", dest="use_bridge", action="store_false", default=True,
                     help="Desactiva el modo doble (puente+maestro) y usa --metadata en modo directo.")
    ap.add_argument("--sheet", default="MAPEO_MUESTRAS",
                     help="Hoja a usar en --metadata cuando se usa el modo directo (--no-bridge).")
    ap.add_argument("--sample-col", default="SampleID", help="Columna con el identificador de muestra (modo directo).")
    ap.add_argument("--group-col", default="Sistema", help="Columna de agrupación (modo directo, por defecto 'Sistema').")
    ap.add_argument("--outdir", type=Path, default=None,
                     help="Carpeta de salida. Por defecto: EMU_propio/EMU{model}_{dataset}.")
    args = ap.parse_args()

    # Defaults derivados del modelo/dataset (solo si no se pasaron explícitamente)
    if args.counts is None:
        args.counts = base / f"EMU{args.model}_{args.dataset}" / "tabla_conteos.tsv"
    if args.outdir is None:
        args.outdir = base / f"EMU{args.model}_{args.dataset}"
    if args.bridge_sheet is None:
        args.bridge_sheet = args.model.upper()

    eprint(f"[INFO] Modelo: {args.model} | dataset: {args.dataset} | counts: {args.counts} | outdir: {args.outdir} | bridge-sheet: {args.bridge_sheet}")

    ensure_dir(args.outdir)

    # Permite pasar el índice de la hoja como número si se prefiere (p. ej. --sheet 1)
    sheet = args.sheet
    if isinstance(sheet, str) and sheet.isdigit():
        sheet = int(sheet)

    eprint(f"[INFO] Cargando tabla de conteos: {args.counts}")
    counts = load_counts(args.counts)
    eprint(f"[INFO] Tabla de conteos: {counts.shape[0]} taxones x {counts.shape[1]} muestras")

    if args.use_bridge:
        eprint(f"[INFO] Modo doble: puente={args.bridge} (hoja {args.bridge_sheet}) + maestro={args.metadata}")
        mapping = load_double_mapping(args.bridge, args.bridge_sheet, args.metadata, counts.columns)
    else:
        eprint(f"[INFO] Modo directo: metadata={args.metadata} (hoja: {sheet})")
        mapping = load_metadata_mapping(args.metadata, sheet, args.sample_col, args.group_col)
    eprint(f"[INFO] Muestras con Sistema asignado: {len(mapping)} / {counts.shape[1]}")

    # Tabla única: muestras + Total_counts + Frecuencia_muestras + Frecuencia_<Sistema>
    tabla = add_total_and_frequency(counts)
    freq_sistemas = add_system_frequencies(counts, mapping)
    tabla = pd.concat([tabla, freq_sistemas.loc[tabla.index]], axis=1)

    write_tsv_and_xlsx(tabla, args.outdir, "tabla_conteos_enriquecida")

    eprint("\n[FIN] Listo. Resumen:")
    eprint(f" - tabla_conteos_enriquecida.(tsv/xlsx): {tabla.shape[0]} taxones, "
           f"{counts.shape[1]} muestras + Total_counts + Frecuencia_muestras + "
           f"{[c for c in freq_sistemas.columns]}")


if __name__ == "__main__":
    main()
