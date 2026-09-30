#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys, re, gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd

TAX_COLS = ["superkingdom","phylum","class","order","family","genus","species"]

def eprint(*a, **k): print(*a, file=sys.stderr, **k)
def ensure_dir(p: Path): p.mkdir(parents=True, exist_ok=True)

def _read_tsv(path: Path, **kw):
    return pd.read_csv(path, sep="\t", **kw)

def _detect_rel_col(df: pd.DataFrame) -> Optional[str]:
    # prioriza nombres con "abund"
    cands = [c for c in df.columns if re.search(r"abund", c, re.I)]
    for c in cands:
        try:
            s = pd.to_numeric(df[c], errors="coerce").fillna(0)
            sm = float(s.sum())
            if 0.9 <= sm <= 1.1:
                return c
        except Exception:
            pass
    # si no, cualquier numérica que sume ~1
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]):
            sm = float(pd.to_numeric(df[c], errors="coerce").fillna(0).sum())
            if 0.9 <= sm <= 1.1:
                return c
    return None

def _build_feature_id(row: pd.Series, rank: str) -> str:
    tax_id = str(row.get("tax_id", "NA")).replace("|","_")
    name = str(row.get(rank, row.get(rank.lower(),"Unassigned")) or "Unassigned").replace("|","_")
    return f"{rank.lower()}|{tax_id}|{name}"

def _qiime_tax_string(row: pd.Series) -> str:
    labels = ["k__","p__","c__","o__","f__","g__","s__"]
    out = []
    for lab, key in zip(labels, TAX_COLS):
        val = str(row.get(key, "Unassigned") or "Unassigned")
        out.append(f"{lab}{val}")
    return "; ".join(out)

def _count_fastq_reads(fq: Path) -> int:
    openf = gzip.open if fq.suffix == ".gz" else open
    n = 0
    with openf(fq, "rt", errors="ignore") as fh:
        for _ in fh: n += 1
    return n // 4

def _sample_paths(emu_outdir: Path) -> Dict[str, Dict[str, Path]]:
    """Devuelve por muestra: rel_tsv, assign_tsv (si existe), fastq (si está en manifest)."""
    samples = {}
    # manifest.tsv para ubicar FASTQ
    manifest = emu_outdir / "manifest.tsv"
    man = None
    if manifest.exists():
        try:
            man = _read_tsv(manifest)
        except Exception:
            man = None

    for sdir in (emu_outdir / "samples").glob("*"):
        if not sdir.is_dir():
            continue
        sample = sdir.name
        rel = None; assign = None
        for f in sdir.glob("*.tsv"):
            n = f.name.lower()
            if "rel" in n and "abundance" in n:
                rel = f
            elif "assign" in n:
                assign = f
        fq = None
        if man is not None and "sample" in man.columns and "fastq" in man.columns:
            row = man.loc[man["sample"].astype(str)==sample]
            if len(row) == 1:
                fq = Path(str(row["fastq"].iloc[0]))
        samples[sample] = {"rel": rel, "assign": assign, "fastq": fq}
    return samples

def _integerize(p: np.ndarray, total: int) -> np.ndarray:
    """Convierte proporciones p a enteros que suman total (método largest remainder)."""
    if total <= 0 or p.size == 0:
        return np.zeros_like(p, dtype=int)
    x = np.clip(p, 0, None).astype(float)
    if x.sum() == 0:
        return np.zeros_like(x, dtype=int)
    x = x / x.sum()
    raw = x * total
    base = np.floor(raw).astype(int)
    rem = total - base.sum()
    if rem > 0:
        frac = raw - base
        idx = np.argsort(-frac)[:rem]
        base[idx] += 1
    return base

def rebuild_counts(emu_outdir: Path, rank: str, pseudo_depth: int) -> Tuple[Path, Path]:
    samples = _sample_paths(emu_outdir)
    if not samples:
        raise SystemExit("No encontré subcarpetas en outdir/samples.")

    counts_tabs = []
    taxonomy_map = {}

    for sample, paths in samples.items():
        rel_tsv = paths["rel"]
        if rel_tsv is None or not rel_tsv.exists():
            eprint(f"[WARN] Salto {sample}: no encontré *_rel-abundance.tsv")
            continue
        df = _read_tsv(rel_tsv)

        # normaliza nombres taxonómicos
        low = {c.lower(): c for c in df.columns}
        for std in TAX_COLS:
            src = low.get(std)
            if src and src != std:
                df.rename(columns={src: std}, inplace=True)

        rel_col = _detect_rel_col(df)
        if rel_col is None:
            eprint(f"[WARN] Salto {sample}: no detecté columna de relativas en {rel_tsv.name}")
            continue

        # total lecturas
        total = None
        if paths["assign"] and paths["assign"].exists():
            try:
                # contar filas menos el header
                nrows = sum(1 for _ in open(paths["assign"], "r")) - 1
                total = max(0, nrows)
            except Exception:
                total = None
        if total is None and paths["fastq"] and Path(paths["fastq"]).exists():
            try:
                total = _count_fastq_reads(Path(paths["fastq"]))
            except Exception:
                total = None
        if total is None:
            total = int(pseudo_depth)
            eprint(f"[INFO] {sample}: uso profundidad fija {total} (no hallé assignments/FASTQ).")

        # construir feature_id + taxonomy string
        feats = []
        taxstr = []
        for _, row in df.iterrows():
            fid = _build_feature_id(row, rank)
            feats.append(fid)
            taxstr.append(_qiime_tax_string(row))
        df["_feature_id"] = feats
        df["_tax_str"] = taxstr

        # guardar taxonomy global (último gana, es el mismo fid)
        for fid, t in zip(feats, taxstr):
            taxonomy_map[fid] = t

        # convertir relativas -> counts
        p = pd.to_numeric(df[rel_col], errors="coerce").fillna(0).to_numpy()
        c = _integerize(p, total)
        counts_tabs.append(pd.DataFrame({"_feature_id": df["_feature_id"], sample: c}).set_index("_feature_id"))

    if not counts_tabs:
        raise SystemExit("No pude reconstruir ninguna muestra.")

    CNT = pd.concat(counts_tabs, axis=1).fillna(0).astype(int)

    # escribir salidas
    out_counts = emu_outdir / "feature_table_counts.tsv"
    CNT.to_csv(out_counts, sep="\t")

    tax_df = pd.DataFrame({"feature_id": list(taxonomy_map.keys()),
                           "taxonomy": list(taxonomy_map.values())}).sort_values("feature_id")
    out_tax = emu_outdir / "taxonomy.tsv"
    if not out_tax.exists():  # no pisar si ya existe
        tax_df.to_csv(out_tax, sep="\t", index=False)

    return out_counts, out_tax

def main():
    ap = argparse.ArgumentParser(description="Reconstruye feature_table_counts.tsv desde rel-abundance + totales.")
    ap.add_argument("--emu-outdir", required=True, type=Path, help="Directorio de salida de EMU (contiene samples/, manifest.tsv).")
    ap.add_argument("--rank", default="species",
                    choices=TAX_COLS, help="Nivel taxonómico para construir feature_id.")
    ap.add_argument("--pseudo-depth", type=int, default=10000,
                    help="Profundidad fija si no hay assignments ni FASTQ.")
    args = ap.parse_args()

    out_counts, out_tax = rebuild_counts(args.emu_outdir.resolve(), args.rank, args.pseudo_depth)
    eprint(f"[OK] Escribí:\n - {out_counts}\n - {out_tax if out_tax.exists() else '(taxonomy ya existía)'}")

if __name__ == "__main__":
    main()

