#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EMU_propio.py — Ejecuta EMU sobre FASTQ 16S y agrega resultados (tablas de conteos/relativas + taxonomía),
con ejecución atómica, limpieza en fallos y filtro previo de lecturas mínimas por muestra.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import shlex
import shutil
import subprocess as sp
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ---------- utilidades básicas ----------

def imprimir_error(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def asegurar_directorio(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def buscar_comando(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def escribir_df_atomico(df: pd.DataFrame, path: Path, sep: str = "\t", index: bool = True):
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep=sep, index=index)
    os.replace(tmp, path)

def mover_arbol_atomico(src: Path, dst: Path):
    try:
        os.replace(src, dst)
    except OSError:
        if dst.exists():
            shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(src, dst, dirs_exist_ok=False)
        shutil.rmtree(src, ignore_errors=True)

def nombre_base_seguro(p: Path) -> str:
    stem = p.name
    m = re.search(r"(SRR|ERR|DRR)\d{5,}", stem)
    if m:
        return m.group(0)
    base = re.sub(r"(\.fastq(\.gz)?|\.fq(\.gz)?)$", "", stem, flags=re.IGNORECASE)
    parent = p.parent.name
    return f"{parent}_{base}"

def leer_tabla_opcional(path: Path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as e:
        imprimir_error(f"[WARN] No se pudo leer {path}: {e}")
        return None

def contar_lecturas_fastq(ruta_fastq: Path) -> int:
    """Cuenta el número total de lecturas en un archivo FASTQ (.fastq o .fastq.gz)."""
    lineas = 0
    try:
        if ruta_fastq.suffix == ".gz":
            with gzip.open(ruta_fastq, "rb") as f:
                for _ in f:
                    lineas += 1
        else:
            with open(ruta_fastq, "rb") as f:
                for _ in f:
                    lineas += 1
        return lineas // 4
    except Exception as e:
        imprimir_error(f"[WARN] No se pudo contar lecturas en {ruta_fastq}: {e}")
        return 0

# ---------- taxonomía y métricas ----------

COLUMNAS_TAXONOMICAS = ["superkingdom","phylum","class","order","family","genus","species"]

def _normalizar_columnas(cols: List[str]) -> Dict[str, str]:
    return {c.lower(): c for c in cols}

def a_cadena_taxonomica_qiime(row: pd.Series) -> str:
    etiquetas = ["k__","p__","c__","o__","f__","g__","s__"]
    out = []
    for lab, key in zip(etiquetas, COLUMNAS_TAXONOMICAS):
        val = str(row.get(key, "Unassigned") or "Unassigned")
        out.append(f"{lab}{val}")
    return "; ".join(out)

# ---------- EMU por muestra ----------

def ejecutar_emu_por_muestra(
    comando_emu: str,
    db: Path,
    fastq: Path,
    carpeta_temp_muestra: Path,
    muestra: str,
    hilos: int = 1,
    conservar_conteos: bool = False,
    conservar_asignaciones: bool = False,
) -> Tuple[Path, Optional[Path], Optional[Path]]:
    asegurar_directorio(carpeta_temp_muestra)
    cmd = [
        comando_emu, "abundance",
        "--db", str(db),
        "--threads", str(hilos),
        "--output-dir", str(carpeta_temp_muestra),
        "--output-basename", muestra,
    ]
    if conservar_conteos:
        cmd.append("--keep-counts")
    if conservar_asignaciones:
        cmd.append("--keep-read-assignments")
    cmd.append(str(fastq))

    imprimir_error("[EMU]", " ".join(shlex.quote(c) for c in cmd))
    sp.run(cmd, check=True)

    tsv_relativa = tsv_conteos = tsv_asignaciones = None
    for f in carpeta_temp_muestra.glob("*.tsv"):
        nombre = f.name.lower()
        if "rel" in nombre and "abundance" in nombre:
            tsv_relativa = f
        elif "count" in nombre:
            tsv_conteos = f
        elif "assign" in nombre or "read-assign" in nombre:
            tsv_asignaciones = f

    if tsv_relativa is None:
        candidato = carpeta_temp_muestra / f"{muestra}_rel-abundance.tsv"
        if candidato.exists():
            tsv_relativa = candidato

    if tsv_relativa is None:
        raise RuntimeError(f"No se encontró la tabla de abundancia relativa para {muestra}.")

    return tsv_relativa, tsv_conteos, tsv_asignaciones

# ---------- agregación y validaciones ----------

def filtrar_por_nivel(df: pd.DataFrame, nivel: str) -> pd.DataFrame:
    col_min = _normalizar_columnas(df.columns)
    col = col_min.get(nivel.lower())
    if col is None:
        return df
    mantener = df[col].fillna("").astype(str) != ""
    return df.loc[mantener].copy()

def construir_id_taxon(row: pd.Series, nivel: str) -> str:
    r = nivel.lower()
    tax_id = str(row.get("tax_id", "NA"))
    nombre = None
    for k in COLUMNAS_TAXONOMICAS:
        if k == r:
            nombre = str(row.get(k, "Unassigned") or "Unassigned")
            break
    if nombre is None:
        nombre = str(row.get(r, "Unassigned") or "Unassigned")
    seguro = nombre.replace("|", "_")
    tid = tax_id.replace("|", "_")
    return f"{r}|{tid}|{seguro}"

def _detectar_columna_relativa(df: pd.DataFrame) -> Optional[str]:
    candidatas = [c for c in df.columns if re.search(r"abund", c, re.I)]
    if candidatas:
        for c in candidatas:
            s = df[c]
            try:
                vals = s.astype(float)
                suma = float(vals.sum())
                if 0.9 <= suma <= 1.1:
                    return c
            except Exception:
                pass
        return candidatas[0]
    num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    for c in num:
        suma = float(df[c].sum())
        if 0.9 <= suma <= 1.1:
            return c
    return None

def _detectar_columna_conteo(df: pd.DataFrame) -> Optional[str]:
    candidatas = [c for c in df.columns if re.fullmatch(r"count[s]?", c, re.I)]
    for c in candidatas:
        if pd.api.types.is_numeric_dtype(df[c]):
            return c
    num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    for c in num:
        suma = float(df[c].sum())
        if suma > 1.1:
            return c
    return None

def agregar_tablas(
    relativas_por_muestra: Dict[str, Path],
    conteos_por_muestra: Dict[str, Path],
    outdir: Path,
    nivel: str,
    abundancia_minima: float = 0.0,
    lecturas_minimas_por_muestra: int = 0,
) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
    asegurar_directorio(outdir)
    tablas_relativas, tablas_conteos = [], []
    filas_taxonomia: Dict[str, str] = {}

    # relativas
    for muestra, tsv in relativas_por_muestra.items():
        df = leer_tabla_opcional(tsv)
        if df is None:
            continue
        mapa_col = _normalizar_columnas(df.columns)
        for std in COLUMNAS_TAXONOMICAS:
            if std not in df.columns:
                src = mapa_col.get(std)
                if src:
                    df.rename(columns={src: std}, inplace=True)
        df = filtrar_por_nivel(df, nivel)
        col_relativa = _detectar_columna_relativa(df)
        if col_relativa is None:
            raise RuntimeError(f"No se detectó columna de abundancia relativa en {tsv}")
        if abundancia_minima > 0:
            df = df[df[col_relativa] >= float(abundancia_minima)].copy()

        ids_taxon, cadenas_tax = [], []
        for _, row in df.iterrows():
            id_taxon = construir_id_taxon(row, nivel)
            ids_taxon.append(id_taxon)
            cadenas_tax.append(a_cadena_taxonomica_qiime(row))
        df["_id_taxon"] = ids_taxon
        df["_cadena_tax"] = cadenas_tax
        for id_taxon, cadena in zip(df["_id_taxon"], df["_cadena_tax"]):
            filas_taxonomia[id_taxon] = cadena
        tablas_relativas.append(df[["_id_taxon", col_relativa]].rename(columns={col_relativa: muestra}).set_index("_id_taxon"))

    # conteos
    for muestra, tsv in conteos_por_muestra.items():
        df = leer_tabla_opcional(tsv)
        if df is None:
            continue
        mapa_col = _normalizar_columnas(df.columns)
        for std in COLUMNAS_TAXONOMICAS:
            if std not in df.columns:
                src = mapa_col.get(std)
                if src:
                    df.rename(columns={src: std}, inplace=True)
        df = filtrar_por_nivel(df, nivel)
        col_conteo = _detectar_columna_conteo(df)
        if col_conteo is None:
            imprimir_error(f"[WARN] No se detectó columna de conteos en {tsv}")
            continue

        total = float(df[col_conteo].sum())
        if lecturas_minimas_por_muestra > 0 and total < lecturas_minimas_por_muestra:
            imprimir_error(f"[WARN] Muestra {muestra} descartada por pocas lecturas: {total} < {lecturas_minimas_por_muestra}")
            continue

        ids_taxon, cadenas_tax = [], []
        for _, row in df.iterrows():
            id_taxon = construir_id_taxon(row, nivel)
            ids_taxon.append(id_taxon)
            cadenas_tax.append(a_cadena_taxonomica_qiime(row))
        df["_id_taxon"] = ids_taxon
        df["_cadena_tax"] = cadenas_tax
        for id_taxon, cadena in zip(df["_id_taxon"], df["_cadena_tax"]):
            filas_taxonomia[id_taxon] = cadena
        tablas_conteos.append(df[["_id_taxon", col_conteo]].rename(columns={col_conteo: muestra}).set_index("_id_taxon"))

    ruta_relativas = ruta_conteos = ruta_taxonomia = None

    if tablas_relativas:
        df_relativas = pd.concat(tablas_relativas, axis=1).fillna(0.0)
        if (df_relativas.sum().sum() <= 0) or (df_relativas.shape[1] == 0):
            raise RuntimeError("Tabla de abundancias relativas vacía o con solo ceros.")
        ruta_relativas = outdir / "tabla_abundancia_relativa.tsv"
        escribir_df_atomico(df_relativas, ruta_relativas, sep="\t", index=True)

    if tablas_conteos:
        df_conteos = pd.concat(tablas_conteos, axis=1).fillna(0.0)
        if (df_conteos.sum().sum() <= 0) or (df_conteos.shape[1] == 0):
            raise RuntimeError("Tabla de conteos vacía o con solo ceros.")
        ruta_conteos = outdir / "tabla_conteos.tsv"
        escribir_df_atomico(df_conteos, ruta_conteos, sep="\t", index=True)

    if filas_taxonomia:
        df_taxonomia = pd.DataFrame({"id_taxon": list(filas_taxonomia.keys()),
                               "taxonomia": list(filas_taxonomia.values())}).sort_values("id_taxon")
        ruta_taxonomia = outdir / "taxonomia.tsv"
        escribir_df_atomico(df_taxonomia, ruta_taxonomia, sep="\t", index=False)

    return (ruta_conteos, ruta_relativas, ruta_taxonomia)

# ---------- descubrimiento inputs ----------

def descubrir_fastqs(carpeta_entrada: Optional[Path], patron: Optional[str], glob_entrada: Optional[str]) -> List[Path]:
    rutas: List[Path] = []
    if glob_entrada:
        for p in Path(".").glob(glob_entrada):
            if p.is_file():
                rutas.append(p.resolve())
    if carpeta_entrada and patron:
        for p in carpeta_entrada.rglob(patron):
            if p.is_file():
                rutas.append(p.resolve())
    return sorted(set(rutas))

# ---------- pre-chequeos/atomicidad ----------

def verificaciones_previas(args) -> Tuple[str, Path, List[Path]]:
    comando_emu = args.emu_cmd or buscar_comando("emu")
    if not comando_emu:
        raise RuntimeError("No se encontró 'emu' en PATH. Especifica --emu-cmd o ajusta tu entorno.")
    try:
        sp.run([comando_emu, "--version"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    except Exception:
        imprimir_error("[WARN] No se pudo verificar 'emu --version' (continuo de todas formas).")

    if not args.db.exists() or not args.db.is_dir():
        raise RuntimeError(f"La BD de EMU no existe o no es carpeta: {args.db}")
    if not any(args.db.iterdir()):
        raise RuntimeError(f"La BD de EMU está vacía: {args.db}")

    outdir = args.outdir.resolve()
    carpeta_build = outdir.parent / (outdir.name + ".__build__")
    if carpeta_build.exists():
        shutil.rmtree(carpeta_build, ignore_errors=True)
    asegurar_directorio(carpeta_build)

    fastqs_crudos = descubrir_fastqs(args.input_dir, args.pattern, args.input_glob)
    if not fastqs_crudos:
        raise RuntimeError("No se encontraron FASTQ con los criterios dados.")

    # Filtro previo de lecturas mínimas por FASTQ
    fastqs = []
    lecturas_minimas = args.min_reads_input
    imprimir_error(f"[INFO] Evaluando lecturas mínimas en los archivos FASTQ (Umbral: >= {lecturas_minimas})...")
    for fq in fastqs_crudos:
        n_lecturas = contar_lecturas_fastq(fq)
        muestra = nombre_base_seguro(fq)
        if n_lecturas >= lecturas_minimas:
            fastqs.append(fq)
            imprimir_error(f"  -> [Aceptado] Muestra '{muestra}': {n_lecturas:,} lecturas.")
        else:
            imprimir_error(f"  -> [Omitido]  Muestra '{muestra}': {n_lecturas:,} lecturas (menor a {lecturas_minimas}).")

    if not fastqs:
        raise RuntimeError(f"Ningún archivo FASTQ superó el umbral de {lecturas_minimas} lecturas.")

    return comando_emu, carpeta_build, fastqs

def finalizar_exito(carpeta_build: Path, outdir: Path, forzar: bool):
    if outdir.exists():
        if not forzar:
            raise RuntimeError(f"--outdir ya existe: {outdir}. Use --force para reemplazar.")
        shutil.rmtree(outdir, ignore_errors=True)
    mover_arbol_atomico(carpeta_build, outdir)

def abortar_y_limpiar(carpeta_build: Path, msg: str, code: int = 1):
    try:
        if carpeta_build and carpeta_build.exists():
            shutil.rmtree(carpeta_build, ignore_errors=True)
    finally:
        imprimir_error(f"[ERROR] {msg}")
        sys.exit(code)

# ---------- CLI ----------

def construir_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pipeline EMU: clasificación taxonómica y métricas de diversidad (16S) con filtro previo de lecturas mínimas.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--db", required=True, type=Path, help="Ruta a la base de datos de EMU.")
    p.add_argument("--emu-cmd", default=None, help="Ruta al ejecutable 'emu' (si no está en PATH).")

    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument("--input-dir", type=Path, help="Carpeta raíz con FASTQ.")
    inp.add_argument("--input-glob", type=str, help="Glob alternativo para FASTQ.")

    p.add_argument("--pattern", type=str, default="*_limpio.fastq", help="Patrón (rglob) para FASTQ bajo --input-dir.")
    p.add_argument("--outdir", required=True, type=Path, help="Directorio de salida final.")
    p.add_argument("--threads", type=int, default=8, help="Hilos para EMU.")
    p.add_argument("--rank", type=str, default="species", choices=[*COLUMNAS_TAXONOMICAS], help="Nivel taxonómico para agregación.")
    p.add_argument("--min-reads-input", type=int, default=500, help="Número mínimo de lecturas en el FASTQ de entrada para ser procesado.")
    p.add_argument("--min-abundance", type=float, default=0.0, help="Umbral de abundancia relativa mínima.")
    p.add_argument("--keep-counts", action="store_true", help="Pedir a EMU que guarde conteos por taxón.")
    p.add_argument("--keep-assignments", action="store_true", help="Pedir a EMU que guarde asignaciones por lectura.")
    p.add_argument("--stop-on-error", dest="stop_on_error", action="store_true", default=True, help="Abortar todo si una muestra falla.")
    p.add_argument("--no-stop-on-error", dest="stop_on_error", action="store_false", help="No abortar todo; continuar saltando muestras fallidas.")
    p.add_argument("--force", action="store_true", help="Si --outdir existe, reemplazarlo al finalizar exitosamente.")
    return p

# ---------- main ----------

def main() -> int:
    args = construir_parser().parse_args()

    try:
        comando_emu, carpeta_build, fastqs = verificaciones_previas(args)
    except Exception as e:
        abortar_y_limpiar(carpeta_build=None, msg=str(e), code=2)

    imprimir_error("=== EMU pipeline (ejecución atómica) ===")
    imprimir_error(f"- DB: {args.db}")
    imprimir_error(f"- EMU: {comando_emu}")
    imprimir_error(f"- FASTQ a procesar: {len(fastqs)}")
    imprimir_error(f"- Rank: {args.rank} | min-abundance: {args.min_abundance}")
    imprimir_error(f"- Mínimo de lecturas por muestra: {args.min_reads_input}")
    imprimir_error(f"- Build dir: {carpeta_build}")

    carpeta_muestras = carpeta_build / "samples"
    asegurar_directorio(carpeta_muestras)

    relativas_por_muestra: Dict[str, Path] = {}
    conteos_por_muestra: Dict[str, Path] = {}

    for fq in fastqs:
        fq = Path(fq)
        muestra = nombre_base_seguro(fq)
        carpeta_muestra_temp = carpeta_muestras / (muestra + ".__tmp__")
        carpeta_muestra_final = carpeta_muestras / muestra

        try:
            asegurar_directorio(carpeta_muestra_temp)
            tsv_relativa, tsv_conteos, tsv_asignaciones = ejecutar_emu_por_muestra(
                comando_emu=comando_emu, db=args.db, fastq=fq,
                carpeta_temp_muestra=carpeta_muestra_temp, muestra=muestra,
                hilos=args.threads,
                conservar_conteos=args.keep_counts,
                conservar_asignaciones=args.keep_assignments,
            )
            mover_arbol_atomico(carpeta_muestra_temp, carpeta_muestra_final)

            if tsv_relativa:
                relativas_por_muestra[muestra] = carpeta_muestra_final / Path(tsv_relativa).name
            if tsv_conteos:
                conteos_por_muestra[muestra] = carpeta_muestra_final / Path(tsv_conteos).name

        except Exception as e:
            shutil.rmtree(carpeta_muestra_temp, ignore_errors=True)
            msg = f"Fallo procesando la muestra '{muestra}': {e}"
            if args.stop_on_error:
                abortar_y_limpiar(carpeta_build=carpeta_build, msg=msg, code=3)
            else:
                imprimir_error("[WARN]", msg)
                continue

    if not relativas_por_muestra:
        abortar_y_limpiar(carpeta_build=carpeta_build, msg="No hay muestras válidas con tabla de abundancias.", code=4)

    # Agregados globales
    try:
        ruta_conteos, ruta_relativas, ruta_taxonomia = agregar_tablas(
            relativas_por_muestra=relativas_por_muestra,
            conteos_por_muestra=conteos_por_muestra,
            outdir=carpeta_build,
            nivel=args.rank,
            abundancia_minima=args.min_abundance,
        )
    except Exception as e:
        abortar_y_limpiar(carpeta_build=carpeta_build, msg=f"Error agregando tablas: {e}", code=6)

    # Finalizar: mover build -> outdir
    try:
        finalizar_exito(carpeta_build, args.outdir.resolve(), forzar=args.force)
    except Exception as e:
        abortar_y_limpiar(carpeta_build=carpeta_build, msg=f"No se pudo mover salidas a {args.outdir}: {e}", code=8)

    imprimir_error("\n=== OK. Salidas en:", args.outdir)
    return 0

if __name__ == "__main__":
    sys.exit(main())
