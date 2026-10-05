#!/usr/bin/env python3

# Rehace tabla_conteos.tsv: abundancia relativa de EMU x total de lecturas limpias de cada muestra,
# redondeando por mayor residuo para que cada muestra sume exactamente su total.

import argparse, sys, re, gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd

COLUMNAS_TAXONOMICAS = ["superkingdom","phylum","class","order","family","genus","species"]

def imprimir_error(*a, **k): print(*a, file=sys.stderr, **k)
def asegurar_directorio(p: Path): p.mkdir(parents=True, exist_ok=True)

def _escribir_xlsx(df: pd.DataFrame, path: Path, index: bool = True):
    # excel no aguanta más de ~1M filas
    n_filas, n_cols = df.shape
    if n_filas > 1_048_575 or n_cols > 16_383:
        imprimir_error(f"[WARN] {path.name}: {n_filas} filas x {n_cols} columnas excede el límite de Excel; se omite el .xlsx.")
        return None
    df.to_excel(path, index=index, engine="openpyxl")
    return path

def _leer_tsv(path: Path, **kw):
    return pd.read_csv(path, sep="\t", **kw)

def _detectar_columna_relativa(df: pd.DataFrame) -> Optional[str]:
    candidatas = [c for c in df.columns if re.search(r"abund", c, re.I)]
    for c in candidatas:
        try:
            s = pd.to_numeric(df[c], errors="coerce").fillna(0)
            suma = float(s.sum())
            if 0.9 <= suma <= 1.1:
                return c
        except Exception:
            pass
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]):
            suma = float(pd.to_numeric(df[c], errors="coerce").fillna(0).sum())
            if 0.9 <= suma <= 1.1:
                return c
    return None

def _construir_id_taxon(row: pd.Series, nivel: str) -> str:
    tax_id = str(row.get("tax_id", "NA")).replace("|","_")
    nombre = str(row.get(nivel, row.get(nivel.lower(),"Unassigned")) or "Unassigned").replace("|","_")
    return f"{nivel.lower()}|{tax_id}|{nombre}"

def _cadena_taxonomica_qiime(row: pd.Series) -> str:
    etiquetas = ["k__","p__","c__","o__","f__","g__","s__"]
    out = []
    for lab, key in zip(etiquetas, COLUMNAS_TAXONOMICAS):
        val = str(row.get(key, "Unassigned") or "Unassigned")
        out.append(f"{lab}{val}")
    return "; ".join(out)

def _contar_lecturas_fastq(fq: Path) -> int:
    abrir = gzip.open if fq.suffix == ".gz" else open
    n = 0
    with abrir(fq, "rt", errors="ignore") as fh:
        for _ in fh: n += 1
    return n // 4

def _cargar_totales_lecturas(path: Optional[Path], col: str) -> Dict[str, int]:
    # devuelve {'barcode01': lecturas, ...} desde resumen_limpieza de dentrim_bam.py
    # (por defecto 'lecturas after filtlong', que son las que entraron a EMU)
    if path is None:
        return {}
    ext = path.suffix.lower()
    if ext == ".xlsx":
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep="\t")
    if "barcode" not in df.columns:
        raise SystemExit(f"ERROR: {path} debe tener una columna 'barcode'. Columnas: {list(df.columns)}")
    if col not in df.columns:
        raise SystemExit(f"ERROR: la columna '{col}' no existe en {path}. Columnas disponibles: {list(df.columns)}")
    totales = {}
    for bc, val in zip(df["barcode"].astype(str).str.strip(), pd.to_numeric(df[col], errors="coerce")):
        if pd.notna(val):
            totales[bc] = int(val)
    return totales

def _rutas_por_muestra(emu_outdir: Path) -> Dict[str, Dict[str, Path]]:
    muestras = {}
    manifiesto = emu_outdir / "manifiesto.tsv"
    man = None
    if manifiesto.exists():
        try:
            man = _leer_tsv(manifiesto)
        except Exception:
            man = None

    for sdir in (emu_outdir / "samples").glob("*"):
        if not sdir.is_dir():
            continue
        muestra = sdir.name
        rel = None; asignaciones = None
        for f in sdir.glob("*.tsv"):
            n = f.name.lower()
            if "rel" in n and "abundance" in n:
                rel = f
            elif "assign" in n:
                asignaciones = f
        fq = None
        if man is not None and "sample" in man.columns and "fastq" in man.columns:
            fila = man.loc[man["sample"].astype(str)==muestra]
            if len(fila) == 1:
                fq = Path(str(fila["fastq"].iloc[0]))
        muestras[muestra] = {"rel": rel, "assign": asignaciones, "fastq": fq}
    return muestras

def _convertir_a_enteros(p: np.ndarray, total: int) -> np.ndarray:
    if total <= 0 or p.size == 0:
        return np.zeros_like(p, dtype=int)
    x = np.clip(p, 0, None).astype(float)
    if x.sum() == 0:
        return np.zeros_like(x, dtype=int)
    x = x / x.sum()
    bruto = x * total
    base = np.floor(bruto).astype(int)
    resto = total - base.sum()
    if resto > 0:
        frac = bruto - base
        idx = np.argsort(-frac)[:resto]
        base[idx] += 1
    return base

def reconstruir_conteos(emu_outdir: Path, nivel: str, profundidad_pseudo: int, totales_lecturas: Dict[str, int]) -> Tuple[Path, Path]:
    muestras = _rutas_por_muestra(emu_outdir)
    if not muestras:
        raise SystemExit("No encontré subcarpetas en outdir/samples.")

    tablas_conteos = []
    mapa_taxonomia = {}

    for muestra, rutas in muestras.items():
        tsv_relativa = rutas["rel"]
        if tsv_relativa is None or not tsv_relativa.exists():
            imprimir_error(f"[WARN] Salto {muestra}: no encontré *_rel-abundance.tsv")
            continue
        df = _leer_tsv(tsv_relativa)

        col_min = {c.lower(): c for c in df.columns}
        for std in COLUMNAS_TAXONOMICAS:
            src = col_min.get(std)
            if src and src != std:
                df.rename(columns={src: std}, inplace=True)

        col_relativa = _detectar_columna_relativa(df)
        if col_relativa is None:
            imprimir_error(f"[WARN] Salto {muestra}: no detecté columna de relativas en {tsv_relativa.name}")
            continue

        total = None
        origen = None

        # el total se busca primero en el resumen de la limpieza
        if totales_lecturas:
            m = re.search(r"barcode0*(\d+)", muestra, re.IGNORECASE)
            if m:
                clave_bc = f"barcode{int(m.group(1)):02d}"
                if clave_bc in totales_lecturas:
                    total = totales_lecturas[clave_bc]
                    origen = f"resumen dentrim ({clave_bc})"

        # si no, en los assignments de EMU
        if total is None and rutas["assign"] and rutas["assign"].exists():
            try:
                nfilas = sum(1 for _ in open(rutas["assign"], "r")) - 1
                total = max(0, nfilas)
                origen = "assignments"
            except Exception:
                total = None
        # o contando el FASTQ
        if total is None and rutas["fastq"] and Path(rutas["fastq"]).exists():
            try:
                total = _contar_lecturas_fastq(Path(rutas["fastq"]))
                origen = "fastq"
            except Exception:
                total = None
        # si nada de eso existe se usa una profundidad fija, que ya no es un conteo real
        if total is None:
            total = int(profundidad_pseudo)
            origen = "pseudo-depth (¡respaldo, no es un total real!)"
            imprimir_error(f"[WARN] {muestra}: no encontré total real (ni resumen dentrim, ni assignments, ni FASTQ); "
                   f"uso profundidad fija {total}.")
        else:
            imprimir_error(f"[INFO] {muestra}: total de lecturas = {total} (origen: {origen})")

        ids_taxon = []
        cadenas_tax = []
        for _, row in df.iterrows():
            id_taxon = _construir_id_taxon(row, nivel)
            ids_taxon.append(id_taxon)
            cadenas_tax.append(_cadena_taxonomica_qiime(row))
        df["_id_taxon"] = ids_taxon
        df["_cadena_tax"] = cadenas_tax

        for id_taxon, t in zip(ids_taxon, cadenas_tax):
            mapa_taxonomia[id_taxon] = t

        p = pd.to_numeric(df[col_relativa], errors="coerce").fillna(0).to_numpy()
        c = _convertir_a_enteros(p, total)
        tablas_conteos.append(pd.DataFrame({"_id_taxon": df["_id_taxon"], muestra: c}).set_index("_id_taxon"))

    if not tablas_conteos:
        raise SystemExit("No pude reconstruir ninguna muestra.")

    CNT = pd.concat(tablas_conteos, axis=1).fillna(0).astype(int)

    ruta_conteos = emu_outdir / "tabla_conteos.tsv"
    CNT.to_csv(ruta_conteos, sep="\t")

    ruta_conteos_xlsx = emu_outdir / "tabla_conteos.xlsx"
    _escribir_xlsx(CNT, ruta_conteos_xlsx, index=True)

    df_taxonomia = pd.DataFrame({"id_taxon": list(mapa_taxonomia.keys()),
                           "taxonomia": list(mapa_taxonomia.values())}).sort_values("id_taxon")
    ruta_taxonomia = emu_outdir / "taxonomia.tsv"
    if not ruta_taxonomia.exists():
        df_taxonomia.to_csv(ruta_taxonomia, sep="\t", index=False)

    # de paso la tabla de relativas en xlsx para revisarla más fácil
    ruta_relativa_xlsx = None
    tsv_relativa_global = emu_outdir / "tabla_abundancia_relativa.tsv"
    if tsv_relativa_global.exists():
        try:
            df_relativa = _leer_tsv(tsv_relativa_global, index_col=0)
            ruta_relativa_xlsx = emu_outdir / "tabla_abundancia_relativa.xlsx"
            _escribir_xlsx(df_relativa, ruta_relativa_xlsx, index=True)
        except Exception as e:
            imprimir_error(f"[WARN] No pude convertir {tsv_relativa_global.name} a xlsx: {e}")
    else:
        imprimir_error(f"[INFO] No encontré {tsv_relativa_global.name} en {emu_outdir}; se omite su conversión a xlsx.")

    return ruta_conteos, ruta_taxonomia

def main():
    base = Path("/home/fenrir/Documentos/Tesis/datos_gulupa/data_gulupa_qs8")
    ap = argparse.ArgumentParser(description="Reconstruye tabla_conteos.tsv desde rel-abundance + totales.")
    ap.add_argument("--model", choices=["hac", "sup"], default="hac",
                     help="Modelo de basecalling: define los valores por defecto de --emu-outdir y --read-totals.")
    ap.add_argument("--dataset", choices=["todo", "results"], default="todo",
                     help="'todo' -> EMU{model}_todo (todas las muestras); "
                          "'results' -> EMU{model}_results (solo muestras filtradas >500 lecturas).")
    ap.add_argument("--emu-outdir", type=Path, default=None,
                     help="Directorio de salida de EMU (contiene samples/, manifiesto.tsv). "
                          "Por defecto: EMU_propio/EMU{model}_{dataset}.")
    ap.add_argument("--rank", default="species", choices=COLUMNAS_TAXONOMICAS, help="Nivel taxonómico para construir id_taxon.")
    ap.add_argument("--pseudo-depth", type=int, default=10000, help="Profundidad fija SOLO si no hay resumen dentrim, ni assignments, ni FASTQ.")
    ap.add_argument("--read-totals", type=Path, default=None,
                     help="Resumen de dentrim_bam.py (barcode + lecturas reales). Tiene prioridad sobre assignments/FASTQ/pseudo-depth. "
                          "Por defecto: Limpieza/{model}_8_trim_edlib/resumen_limpieza_{model}_edlib.xlsx.")
    ap.add_argument("--read-totals-col", default="lecturas after filtlong",
                     help="Columna de --read-totals a usar como total real de lecturas por muestra.")
    ap.add_argument("--no-read-totals", dest="use_read_totals", action="store_false", default=True,
                     help="Ignora --read-totals aunque exista (usa assignments/FASTQ/pseudo-depth).")
    args = ap.parse_args()

    # rutas por defecto según modelo y dataset
    if args.emu_outdir is None:
        args.emu_outdir = base / "EMU_propio" / f"EMU{args.model}_{args.dataset}"
    if args.read_totals is None:
        args.read_totals = base / "Limpieza" / f"{args.model}_8_trim_edlib" / f"resumen_limpieza_{args.model}_edlib.xlsx"

    imprimir_error(f"[INFO] Modelo: {args.model} | dataset: {args.dataset} | emu-outdir: {args.emu_outdir} | read-totals: {args.read_totals}")

    totales_lecturas = {}
    if args.use_read_totals and args.read_totals and args.read_totals.exists():
        totales_lecturas = _cargar_totales_lecturas(args.read_totals, args.read_totals_col)
        imprimir_error(f"[INFO] Totales reales cargados desde {args.read_totals}: {len(totales_lecturas)} barcodes.")
    elif args.use_read_totals:
        imprimir_error(f"[WARN] No encontré --read-totals en {args.read_totals}; se usará assignments/FASTQ/pseudo-depth.")

    ruta_conteos, ruta_taxonomia = reconstruir_conteos(args.emu_outdir.resolve(), args.rank, args.pseudo_depth, totales_lecturas)
    ruta_conteos_xlsx = ruta_conteos.with_suffix(".xlsx")
    ruta_relativa_xlsx = args.emu_outdir.resolve() / "tabla_abundancia_relativa.xlsx"

    imprimir_error(f"[OK] Escribí:\n - {ruta_conteos}")
    if ruta_conteos_xlsx.exists():
        imprimir_error(f" - {ruta_conteos_xlsx}")
    imprimir_error(f" - {ruta_taxonomia if ruta_taxonomia.exists() else '(taxonomia ya existía)'}")
    if ruta_relativa_xlsx.exists():
        imprimir_error(f" - {ruta_relativa_xlsx}")

if __name__ == "__main__":
    main()
