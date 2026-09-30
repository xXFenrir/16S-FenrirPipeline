#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EMU.py — Ejecuta EMU sobre FASTQ 16S, agrega resultados y calcula diversidad alfa/beta + PCoA,
con ejecución atómica y limpieza en fallos (no deja salidas parciales).

Salidas (solo si TODO OK, en --outdir):
  - Por muestra (subcarpeta):
      sample_rel-abundance.tsv, sample_counts.tsv (si --keep-counts),
      sample_read-assignments.tsv (si --keep-assignments)
  - Globales:
      manifest.tsv
      feature_table_counts.tsv
      feature_table_relabund.tsv
      taxonomy.tsv
      alpha_diversity.tsv
      braycurtis_dm.tsv, pcoa_braycurtis_coords.tsv, pcoa_braycurtis_eigvals.tsv, pcoa_braycurtis_variance.tsv
      jaccard_dm.tsv,    pcoa_jaccard_coords.tsv,    pcoa_jaccard_eigvals.tsv,    pcoa_jaccard_variance.tsv
"""

from __future__ import annotations

import argparse
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

try:
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except Exception:
    _HAS_MPL = False

# ---------- utilidades básicas ----------

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def atomic_write_df(df: pd.DataFrame, path: Path, sep: str = "\t", index: bool = True):
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep=sep, index=index)
    os.replace(tmp, path)

def atomic_move_tree(src: Path, dst: Path):
    try:
        os.replace(src, dst)  # rename atómico (misma partición)
    except OSError:
        # fallback a copytree + rmtree si diferente FS
        if dst.exists():
            shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(src, dst, dirs_exist_ok=False)
        shutil.rmtree(src, ignore_errors=True)

def safe_basename_noext(p: Path) -> str:
    stem = p.name
    m = re.search(r"(SRR|ERR|DRR)\d{5,}", stem)
    if m:
        return m.group(0)
    base = re.sub(r"(\.fastq(\.gz)?|\.fq(\.gz)?)$", "", stem, flags=re.IGNORECASE)
    parent = p.parent.name
    return f"{parent}_{base}"

def read_table_maybe(path: Path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as e:
        eprint(f"[WARN] No se pudo leer {path}: {e}")
        return None

# ---------- taxonomía y métricas ----------

TAX_COLS = ["superkingdom","phylum","class","order","family","genus","species"]

def _norm_cols(cols: List[str]) -> Dict[str, str]:
    return {c.lower(): c for c in cols}

def to_qiime_tax_string(row: pd.Series) -> str:
    labels = ["k__","p__","c__","o__","f__","g__","s__"]
    out = []
    for lab, key in zip(labels, TAX_COLS):
        val = str(row.get(key, "Unassigned") or "Unassigned")
        out.append(f"{lab}{val}")
    return "; ".join(out)

def chao1(counts: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=float)
    Sobs = np.sum(counts > 0)
    f1 = np.sum(counts == 1)
    f2 = np.sum(counts == 2)
    if f2 > 0:
        return float(Sobs + (f1 * f1) / (2.0 * f2))
    return float(Sobs + (f1 * (f1 - 1.0)) / 2.0)

def shannon_entropy(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    return float(-np.sum(p * np.log(p)))

def simpson_index(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    return float(1.0 - np.sum(p * p))

# ---------- EMU por muestra ----------

def run_emu_for_sample(
    emu_cmd: str,
    db: Path,
    fastq: Path,
    outdir_sample_tmp: Path,
    sample: str,
    threads: int = 1,
    keep_counts: bool = False,
    keep_assignments: bool = False,
) -> Tuple[Path, Optional[Path], Optional[Path]]:
    ensure_dir(outdir_sample_tmp)
    cmd = [
        emu_cmd, "abundance",
        "--db", str(db),
        "--threads", str(threads),
        "--output-dir", str(outdir_sample_tmp),
        "--output-basename", sample,
    ]
    if keep_counts:
        cmd.append("--keep-counts")
    if keep_assignments:
        cmd.append("--keep-read-assignments")
    cmd.append(str(fastq))

    eprint("[EMU]", " ".join(shlex.quote(c) for c in cmd))
    sp.run(cmd, check=True)

    rel_tsv = counts_tsv = assigns_tsv = None
    for f in outdir_sample_tmp.glob("*.tsv"):
        name = f.name.lower()
        if "rel" in name and "abundance" in name:
            rel_tsv = f
        elif "count" in name:
            counts_tsv = f
        elif "assign" in name or "read-assign" in name:
            assigns_tsv = f

    if rel_tsv is None:
        # fallback
        cand = outdir_sample_tmp / f"{sample}_rel-abundance.tsv"
        if cand.exists():
            rel_tsv = cand

    if rel_tsv is None:
        raise RuntimeError(f"No se encontró la tabla de abundancia relativa para {sample}.")

    return rel_tsv, counts_tsv, assigns_tsv

# ---------- agregación y validaciones ----------

def select_rank(df: pd.DataFrame, rank: str) -> pd.DataFrame:
    col_low = _norm_cols(df.columns)
    at = col_low.get(rank.lower())
    if at is None:
        return df
    keep = df[at].fillna("").astype(str) != ""
    return df.loc[keep].copy()

def build_feature_id(row: pd.Series, rank: str) -> str:
    r = rank.lower()
    tax_id = str(row.get("tax_id", "NA"))
    name = None
    for k in TAX_COLS:
        if k == r:
            name = str(row.get(k, "Unassigned") or "Unassigned")
            break
    if name is None:
        name = str(row.get(r, "Unassigned") or "Unassigned")
    safe = name.replace("|", "_")
    tid = tax_id.replace("|", "_")
    return f"{r}|{tid}|{safe}"

def _detect_rel_col(df: pd.DataFrame) -> Optional[str]:
    candidates = [c for c in df.columns if re.search(r"abund", c, re.I)]
    if candidates:
        for c in candidates:
            s = df[c]
            try:
                vals = s.astype(float)
                sm = float(vals.sum())
                if 0.9 <= sm <= 1.1:
                    return c
            except Exception:
                pass
        return candidates[0]
    num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    for c in num:
        sm = float(df[c].sum())
        if 0.9 <= sm <= 1.1:
            return c
    return None

def _detect_count_col(df: pd.DataFrame) -> Optional[str]:
    candidates = [c for c in df.columns if re.fullmatch(r"count[s]?", c, re.I)]
    for c in candidates:
        if pd.api.types.is_numeric_dtype(df[c]):
            return c
    num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    for c in num:
        sm = float(df[c].sum())
        if sm > 1.1:
            return c
    return None

def aggregate_tables(
    per_sample_rel: Dict[str, Path],
    per_sample_counts: Dict[str, Path],
    outdir: Path,
    rank: str,
    min_abundance: float = 0.0,
    min_reads_per_sample: int = 0,
) -> Tuple[Optional[Path], Optional[Path], Optional[Path], Optional[Path]]:
    ensure_dir(outdir)
    rel_tabs, cnt_tabs = [], []
    tax_rows: Dict[str, str] = {}
    kept_samples: List[str] = []

    # relativas
    for sample, tsv in per_sample_rel.items():
        df = read_table_maybe(tsv)
        if df is None:
            continue
        colmap = _norm_cols(df.columns)
        for std in TAX_COLS:
            if std not in df.columns:
                src = colmap.get(std)
                if src:
                    df.rename(columns={src: std}, inplace=True)
        df = select_rank(df, rank)
        rel_col = _detect_rel_col(df)
        if rel_col is None:
            raise RuntimeError(f"No se detectó columna de abundancia relativa en {tsv}")
        if min_abundance > 0:
            df = df[df[rel_col] >= float(min_abundance)].copy()

        feats, tax_strings = [], []
        for _, row in df.iterrows():
            fid = build_feature_id(row, rank)
            feats.append(fid)
            tax_strings.append(to_qiime_tax_string(row))
        df["_feature_id"] = feats
        df["_tax_str"] = tax_strings
        for fid, tstr in zip(df["_feature_id"], df["_tax_str"]):
            tax_rows[fid] = tstr
        rel_tabs.append(df[["_feature_id", rel_col]].rename(columns={rel_col: sample}).set_index("_feature_id"))
        kept_samples.append(sample)

    # conteos
    for sample, tsv in per_sample_counts.items():
        df = read_table_maybe(tsv)
        if df is None:
            continue
        colmap = _norm_cols(df.columns)
        for std in TAX_COLS:
            if std not in df.columns:
                src = colmap.get(std)
                if src:
                    df.rename(columns={src: std}, inplace=True)
        df = select_rank(df, rank)
        count_col = _detect_count_col(df)
        if count_col is None:
            eprint(f"[WARN] No se detectó columna de conteos en {tsv}")
            continue
        # filtro por lecturas mínimas por muestra (si aplica)
        tot = float(df[count_col].sum())
        if min_reads_per_sample > 0 and tot < min_reads_per_sample:
            eprint(f"[WARN] Muestra {sample} descartada por pocas lecturas: {tot} < {min_reads_per_sample}")
            continue

        feats, tax_strings = [], []
        for _, row in df.iterrows():
            fid = build_feature_id(row, rank)
            feats.append(fid)
            tax_strings.append(to_qiime_tax_string(row))
        df["_feature_id"] = feats
        df["_tax_str"] = tax_strings
        for fid, tstr in zip(df["_feature_id"], df["_tax_str"]):
            tax_rows[fid] = tstr
        cnt_tabs.append(df[["_feature_id", count_col]].rename(columns={count_col: sample}).set_index("_feature_id"))

    rel_path = cnt_path = tax_path = alpha_path = None

    if rel_tabs:
        rel_df = pd.concat(rel_tabs, axis=1).fillna(0.0)
        if (rel_df.sum().sum() <= 0) or (rel_df.shape[1] == 0):
            raise RuntimeError("Tabla de abundancias relativas vacía o con solo ceros.")
        rel_path = outdir / "feature_table_relabund.tsv"
        atomic_write_df(rel_df, rel_path, sep="\t", index=True)

    if cnt_tabs:
        cnt_df = pd.concat(cnt_tabs, axis=1).fillna(0.0)
        if (cnt_df.sum().sum() <= 0) or (cnt_df.shape[1] == 0):
            raise RuntimeError("Tabla de conteos vacía o con solo ceros.")
        cnt_path = outdir / "feature_table_counts.tsv"
        atomic_write_df(cnt_df, cnt_path, sep="\t", index=True)

    if tax_rows:
        tax_df = pd.DataFrame({"feature_id": list(tax_rows.keys()),
                               "taxonomy": list(tax_rows.values())}).sort_values("feature_id")
        tax_path = outdir / "taxonomy.tsv"
        atomic_write_df(tax_df, tax_path, sep="\t", index=False)

    # alpha-diversidad
    if cnt_tabs or rel_tabs:
        table = None
        if cnt_tabs:
            table = pd.concat(cnt_tabs, axis=1).fillna(0.0)
        elif rel_tabs:
            table = pd.concat(rel_tabs, axis=1).fillna(0.0)

        rows = []
        for sample in table.columns:
            x = table[sample].values
            if cnt_tabs:
                total = x.sum()
                p = x / total if total > 0 else x
                obs = int(np.sum(x > 0))
                c1 = chao1(x)
            else:
                p = x
                obs = int(np.sum(p > 0))
                c1 = float(obs)
            H = shannon_entropy(p)
            S = simpson_index(p)
            rows.append({"sample": sample, "observed": obs, "chao1": c1, "shannon": H, "simpson": S})
        alpha_df = pd.DataFrame(rows).set_index("sample")
        if alpha_df.isna().any().any():
            raise RuntimeError("Alpha-diversidad contiene NaNs inesperados.")
        alpha_path = outdir / "alpha_diversity.tsv"
        atomic_write_df(alpha_df, alpha_path, sep="\t", index=True)

    return (cnt_path, rel_path, tax_path, alpha_path)

# ---------- Beta + PCoA ----------

def _bray_curtis_dm(X: np.ndarray) -> np.ndarray:
    n = X.shape[0]
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        xi = X[i]
        for j in range(i + 1, n):
            xj = X[j]
            num = np.abs(xi - xj).sum()
            den = (xi + xj).sum()
            D[i, j] = D[j, i] = 0.0 if den == 0 else num / den
    return D

def _jaccard_dm_bin(Xbin: np.ndarray) -> np.ndarray:
    n = Xbin.shape[0]
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        ai = Xbin[i] > 0
        for j in range(i + 1, n):
            aj = Xbin[j] > 0
            inter = np.logical_and(ai, aj).sum()
            uni = np.logical_or(ai, aj).sum()
            D[i, j] = D[j, i] = 1.0 if uni == 0 else 1.0 - (inter / uni)
    return D

def _classical_pcoa(D: np.ndarray):
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pos = eigvals > 1e-12
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]
    coords = eigvecs_pos * np.sqrt(eigvals_pos)
    var = eigvals_pos / eigvals_pos.sum() if eigvals_pos.size else np.array([])
    return coords, eigvals_pos, var

def _save_dm_and_pcoa(D: np.ndarray, sample_ids: List[str], outdir: Path, prefix: str) -> Dict[str, Optional[Path]]:
    ensure_dir(outdir)
    dm_path = outdir / f"{prefix}_dm.tsv"
    atomic_write_df(pd.DataFrame(D, index=sample_ids, columns=sample_ids), dm_path, sep="\t", index=True)

    coords, eigvals, var = _classical_pcoa(D)
    k = min(3, coords.shape[1]) if coords.size else 0

    coords_path = outdir / f"pcoa_{prefix}_coords.tsv"
    eig_path    = outdir / f"pcoa_{prefix}_eigvals.tsv"
    var_path    = outdir / f"pcoa_{prefix}_variance.tsv"

    atomic_write_df(pd.DataFrame(coords[:, :k], index=sample_ids,
                                 columns=[f"PC{i+1}" for i in range(k)]), coords_path, sep="\t", index=True)
    atomic_write_df(pd.DataFrame({"eigenvalue": eigvals}), eig_path, sep="\t", index=False)
    atomic_write_df(pd.DataFrame({"proportion_explained": var}), var_path, sep="\t", index=False)

    fig_path = None
    if _HAS_MPL and k >= 2:
        plt.figure()
        plt.scatter(coords[:, 0], coords[:, 1])
        for i, sid in enumerate(sample_ids):
            plt.text(coords[i, 0], coords[i, 1], str(sid))
        pc1 = f"{(var[0]*100):.1f}%" if var.size > 0 else ""
        pc2 = f"{(var[1]*100):.1f}%" if var.size > 1 else ""
        plt.xlabel(f"PC1 ({pc1})"); plt.ylabel(f"PC2 ({pc2})")
        plt.title(f"PCoA — {prefix}")
        plt.tight_layout()
        fig_path = outdir / f"pcoa_{prefix}.png"
        plt.savefig(fig_path, dpi=150)
        plt.close()

    return {"dm": dm_path, "coords": coords_path, "eigvals": eig_path, "variance": var_path, "figure": fig_path}

def write_beta_and_pcoa_from_table(feature_table_tsv: Path, outdir: Path,
                                   mode: str, normalize_to_relative: bool = False) -> Dict[str, Optional[Path]]:
    df = pd.read_csv(feature_table_tsv, sep="\t", index_col=0)
    if any(pd.Series(df.index.astype(str)).str.match(r"^(SRR|ERR|DRR)", na=False)):
        df = df.T
    X = df.values.T  # samples x features
    sample_ids = list(df.columns)

    if mode.lower() == "braycurtis":
        if normalize_to_relative:
            row_sums = X.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1.0
            X = X / row_sums
        D = _bray_curtis_dm(X)
        return _save_dm_and_pcoa(D, sample_ids, outdir, "braycurtis")
    elif mode.lower() == "jaccard":
        Xbin = (X > 0).astype(int)
        D = _jaccard_dm_bin(Xbin)
        return _save_dm_and_pcoa(D, sample_ids, outdir, "jaccard")
    else:
        raise ValueError("mode debe ser 'braycurtis' o 'jaccard'")

# ---------- descubrimiento inputs ----------

def discover_fastqs(input_dir: Optional[Path], pattern: Optional[str], input_glob: Optional[str]) -> List[Path]:
    paths: List[Path] = []
    if input_glob:
        for p in Path(".").glob(input_glob):
            if p.is_file():
                paths.append(p.resolve())
    if input_dir and pattern:
        for p in input_dir.rglob(pattern):
            if p.is_file():
                paths.append(p.resolve())
    return sorted(set(paths))

def merge_assignments(assign_paths: Dict[str, Path], out_csv: Path) -> Optional[Path]:
    rows = []
    for sample, p in assign_paths.items():
        df = read_table_maybe(p)
        if df is None:
            continue
        df.insert(0, "sample", sample)
        rows.append(df)
    if not rows:
        return None
    big = pd.concat(rows, axis=0, ignore_index=True)
    atomic_write_df(big, out_csv, sep=",", index=False)
    return out_csv

# ---------- pre-chequeos/atomicidad ----------

def preflight_checks(args) -> Tuple[str, Path, List[Path]]:
    # EMU
    emu_cmd = args.emu_cmd or which("emu")
    if not emu_cmd:
        raise RuntimeError("No se encontró 'emu' en PATH. Especifica --emu-cmd o ajusta tu entorno.")
    try:
        sp.run([emu_cmd, "--version"], stdout=sp.PIPE, stderr=sp.PIPE, check=True)
    except Exception:
        eprint("[WARN] No se pudo verificar 'emu --version' (continuo de todas formas).")

    # BD
    if not args.db.exists() or not args.db.is_dir():
        raise RuntimeError(f"La BD de EMU no existe o no es carpeta: {args.db}")
    if not any(args.db.iterdir()):
        raise RuntimeError(f"La BD de EMU está vacía: {args.db}")

    # OUTDIR build
    outdir = args.outdir.resolve()
    build = outdir.parent / (outdir.name + ".__build__")
    if build.exists():
        shutil.rmtree(build, ignore_errors=True)
    ensure_dir(build)

    # FASTQ
    fastqs = discover_fastqs(args.input_dir, args.pattern, args.input_glob)
    if not fastqs:
        raise RuntimeError("No se encontraron FASTQ con los criterios dados.")
    return emu_cmd, build, fastqs

def finalize_success(build: Path, outdir: Path, force: bool):
    if outdir.exists():
        if not force:
            raise RuntimeError(f"--outdir ya existe: {outdir}. Use --force para reemplazar.")
        shutil.rmtree(outdir, ignore_errors=True)
    atomic_move_tree(build, outdir)

def abort_and_cleanup(build: Path, msg: str, code: int = 1):
    try:
        if build and build.exists():
            shutil.rmtree(build, ignore_errors=True)
    finally:
        eprint(f"[ERROR] {msg}")
        sys.exit(code)

# ---------- CLI ----------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pipeline EMU: clasificación taxonómica y métricas de diversidad (16S) con ejecución atómica.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--db", required=True, type=Path, help="Ruta a la base de datos de EMU.")
    p.add_argument("--emu-cmd", default=None, help="Ruta al ejecutable 'emu' (si no está en PATH).")

    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument("--input-dir", type=Path, help="Carpeta raíz con FASTQ.")
    inp.add_argument("--input-glob", type=str, help="Glob alternativo para FASTQ (p.ej., 'results/*/*_final.fastq').")

    p.add_argument("--pattern", type=str, default="*_final.fastq", help="Patrón (rglob) para FASTQ bajo --input-dir.")
    p.add_argument("--outdir", required=True, type=Path, help="Directorio de salida final.")
    p.add_argument("--threads", type=int, default=8, help="Hilos para EMU.")
    p.add_argument("--rank", type=str, default="species", choices=[*TAX_COLS], help="Nivel taxonómico para agregación.")
    p.add_argument("--min-abundance", type=float, default=0.0, help="Umbral de abundancia relativa mínima.")
    p.add_argument("--keep-counts", action="store_true", help="Pedir a EMU que guarde conteos por taxón.")
    p.add_argument("--keep-assignments", action="store_true", help="Pedir a EMU que guarde asignaciones por lectura.")
    p.add_argument("--do-pcoa", action="store_true", help="Calcular Bray–Curtis/Jaccard y PCoA.")
    p.add_argument("--pcoa-relative", action="store_true",
                   help="Para Bray–Curtis, usar relativas (si existen) o normalizar counts a relativas.")
    p.add_argument("--min-reads-per-sample", type=int, default=0,
                   help="Descarta muestras con menos de N lecturas en counts antes de métricas.")
    p.add_argument("--stop-on-error", dest="stop_on_error", action="store_true", default=True,
                   help="Abortar todo si una muestra falla (por defecto ON).")
    p.add_argument("--no-stop-on-error", dest="stop_on_error", action="store_false",
                   help="No abortar todo; continuar saltando muestras fallidas.")
    p.add_argument("--force", action="store_true",
                   help="Si --outdir existe, reemplazarlo al finalizar exitosamente.")
    return p

# ---------- main ----------

def main() -> int:
    args = build_parser().parse_args()

    try:
        emu_cmd, build_dir, fastqs = preflight_checks(args)
    except Exception as e:
        abort_and_cleanup(build=None, msg=str(e), code=2)  # no build aún

    eprint("=== EMU pipeline (ejecución atómica) ===")
    eprint(f"- DB: {args.db}")
    eprint(f"- EMU: {emu_cmd}")
    eprint(f"- FASTQ: {len(fastqs)}")
    eprint(f"- Rank: {args.rank} | min-abundance: {args.min_abundance}")
    eprint(f"- Stop on error: {args.stop_on_error}")
    eprint(f"- Build dir: {build_dir}")

    # subcarpeta para salidas por muestra
    samples_root = build_dir / "samples"
    ensure_dir(samples_root)

    per_sample_rel: Dict[str, Path] = {}
    per_sample_cnt: Dict[str, Path] = {}
    per_sample_assign: Dict[str, Path] = {}
    manifest_rows = []

    # Ejecutar EMU por muestra dentro del build
    for fq in fastqs:
        fq = Path(fq)
        sample = safe_basename_noext(fq)
        sample_tmp = samples_root / (sample + ".__tmp__")
        sample_final = samples_root / sample

        try:
            ensure_dir(sample_tmp)
            rel_tsv, cnt_tsv, asg_tsv = run_emu_for_sample(
                emu_cmd=emu_cmd, db=args.db, fastq=fq,
                outdir_sample_tmp=sample_tmp, sample=sample,
                threads=args.threads,
                keep_counts=args.keep_counts,
                keep_assignments=args.keep_assignments,
            )
            # mover carpeta tmp -> final (atómico)
            atomic_move_tree(sample_tmp, sample_final)

            manifest_rows.append({
                "sample": sample,
                "fastq": str(fq),
                "outdir_sample": str(sample_final),
                "rel_abundance_tsv": str(sample_final / Path(rel_tsv).name),
                "counts_tsv": str(sample_final / Path(cnt_tsv).name) if cnt_tsv else "",
                "assignments_tsv": str(sample_final / Path(asg_tsv).name) if asg_tsv else "",
            })

            if rel_tsv:
                per_sample_rel[sample] = sample_final / Path(rel_tsv).name
            if cnt_tsv:
                per_sample_cnt[sample] = sample_final / Path(cnt_tsv).name
            if asg_tsv:
                per_sample_assign[sample] = sample_final / Path(asg_tsv).name

        except Exception as e:
            # limpiar tmp de esa muestra y decidir continuar o abortar
            shutil.rmtree(sample_tmp, ignore_errors=True)
            msg = f"Fallo procesando la muestra '{sample}': {e}"
            if args.stop_on_error:
                abort_and_cleanup(build=build_dir, msg=msg, code=3)
            else:
                eprint("[WARN]", msg)
                continue

    # Verificación: ¿quedó al menos 1 muestra válida?
    if not per_sample_rel:
        abort_and_cleanup(build=build_dir, msg="No hay muestras válidas con tabla de abundancias.", code=4)

    # manifest
    try:
        manifest_df = pd.DataFrame(manifest_rows).sort_values("sample")
        atomic_write_df(manifest_df, build_dir / "manifest.tsv", sep="\t", index=False)
    except Exception as e:
        abort_and_cleanup(build=build_dir, msg=f"Error escribiendo manifest: {e}", code=5)

    # Agregados globales
    try:
        cnt_path, rel_path, tax_path, alpha_path = aggregate_tables(
            per_sample_rel=per_sample_rel,
            per_sample_counts=per_sample_cnt,
            outdir=build_dir,
            rank=args.rank,
            min_abundance=args.min_abundance,
            min_reads_per_sample=args.min_reads_per_sample,
        )
        if not (cnt_path or rel_path):
            raise RuntimeError("No se pudo construir ninguna feature table.")
    except Exception as e:
        abort_and_cleanup(build=build_dir, msg=f"Error agregando tablas: {e}", code=6)

    # Beta + PCoA
    if args.do_pcoa:
        try:
            if args.pcoa_relative and rel_path and rel_path.exists():
                write_beta_and_pcoa_from_table(rel_path, build_dir, mode="braycurtis", normalize_to_relative=False)
            else:
                if cnt_path and cnt_path.exists():
                    write_beta_and_pcoa_from_table(cnt_path, build_dir, mode="braycurtis", normalize_to_relative=True)
                elif rel_path and rel_path.exists():
                    write_beta_and_pcoa_from_table(rel_path, build_dir, mode="braycurtis", normalize_to_relative=False)

            if cnt_path and cnt_path.exists():
                write_beta_and_pcoa_from_table(cnt_path, build_dir, mode="jaccard")
            elif rel_path and rel_path.exists():
                write_beta_and_pcoa_from_table(rel_path, build_dir, mode="jaccard")
        except Exception as e:
            abort_and_cleanup(build=build_dir, msg=f"Error en beta-diversidad/PCoA: {e}", code=7)

    # Finalizar: mover build -> outdir
    try:
        finalize_success(build_dir, args.outdir.resolve(), force=args.force)
    except Exception as e:
        abort_and_cleanup(build=build_dir, msg=f"No se pudo mover salidas a {args.outdir}: {e}", code=8)

    eprint("\n=== OK. Salidas en:", args.outdir)
    return 0

if __name__ == "__main__":
    sys.exit(main())

