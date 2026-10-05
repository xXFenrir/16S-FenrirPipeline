#!/usr/bin/env python3

# Limpieza de los FASTQ del repositorio (Objetivo 1), igual que en el artículo:
# primero pychopper (orienta y recorta entre primers) y luego filtlong (calidad y 1000-1700 pb).
# Deja solo <muestra>_clean.fastq por muestra y dentrim_summary.csv en la carpeta de salida.

import argparse, csv, glob, gzip, os, shlex, shutil, subprocess, sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def _must_exist(path: str, kind: str) -> None:
    if kind == "dir" and not os.path.isdir(path):
        sys.exit(f"ERROR: No existe directorio: {path}")
    if kind == "file" and not os.path.isfile(path):
        sys.exit(f"ERROR: No existe archivo: {path}")

def _which_or_die(cmd: str) -> None:
    if shutil.which(cmd) is None:
        sys.exit(f"ERROR: No se encontró '{cmd}' en PATH")

def _pychopper_python() -> Optional[str]:
    pc = shutil.which("pychopper")
    if not pc: return None
    try:
        with open(pc, "rb") as fh:
            first = fh.readline().decode("utf-8", errors="ignore").strip()
        if first.startswith("#!"):
            return first[2:].strip().split()[0]
    except Exception:
        pass
    return None

def _assert_edlib_available_for_pychopper() -> None:
    interp = _pychopper_python()
    if not interp: return
    test = subprocess.run([interp, "-c", "import edlib"], capture_output=True)
    if test.returncode != 0:
        sys.exit(
            "ERROR: El 'pychopper' requiere el módulo 'edlib' en su intérprete.\n"
            f"       Intérprete: {interp}\n"
            "       Instálalo (e.g. `python -m pip install edlib`) o usa un env conda con pychopper+edlib."
        )

def _sample_name(f: str) -> str:
    name = Path(f).name
    for suf in (".fastq.gz",".fq.gz",".fastq",".fq"):
        if name.endswith(suf): return name[:-len(suf)]
    return Path(f).stem

def _is_gz(p: str) -> bool: return p.endswith(".gz")

def _fastq_count(path: str) -> int:
    n = 0
    if _is_gz(path):
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
            for _ in fh: n += 1
    else:
        with open(path, "rt", encoding="utf-8", errors="ignore") as fh:
            for _ in fh: n += 1
    return n // 4

def _safe_unlink(path: str) -> None:
    try: os.remove(path)
    except FileNotFoundError: pass

def _run_bash(cmd: str, verbose: bool = True) -> int:
    if verbose: print(f"$ {cmd}")
    return subprocess.run(cmd, shell=True, executable="/bin/bash").returncode

def _write_metrics_tsv(path: str, metrics: Dict[str, int]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["bucket","count_reads"])
        for k,v in metrics.items(): w.writerow([k,v])


def process_sample_paperlike(
    fq: str,
    sdir: str,
    primers: str,
    pconfig: str,
    qscore_pc: int,
    qcut_pc: float,
    mapper: str,
    threads: int,
    post_minlen: int,
    post_maxlen: int,
    fl_min_mean_q: int,
    save_reports: bool,
    keep_intermediates: bool,
    verbose: bool
) -> Tuple[str, Dict[str,int]]:
    sample = _sample_name(fq)
    os.makedirs(sdir, exist_ok=True)
    oriented_tmp = os.path.join(sdir, f"{sample}__oriented_tmp.fastq")
    clean_path   = os.path.join(sdir, f"{sample}_clean.fastq")

    # 1) pychopper
    cmd1 = (
        f"pychopper -m {mapper} -b {shlex.quote(primers)} -c {shlex.quote(pconfig)} "
        f"-Q {qscore_pc} -z {post_minlen} -t {threads} -Y 0 -q {qcut_pc} "
    )
    if save_reports and keep_intermediates:
        report_pdf = os.path.join(sdir, f"{sample}_report.pdf")
        stats_tsv  = os.path.join(sdir, f"{sample}_stats.tsv")
        scores_tsv = os.path.join(sdir, f"{sample}_scores.tsv")
        cmd1 += (
            f"-r {shlex.quote(report_pdf)} -S {shlex.quote(stats_tsv)} -A {shlex.quote(scores_tsv)} "
            f"-K {shlex.quote(os.path.join(sdir, f'{sample}_qc_fail.fastq'))} "
            f"-l {shlex.quote(os.path.join(sdir, f'{sample}_len_fail.fastq'))} "
            f"-u {shlex.quote(os.path.join(sdir, f'{sample}_unclassified.fastq'))} "
            f"-w {shlex.quote(os.path.join(sdir, f'{sample}_rescued.fastq'))} "
        )
    cmd1 += f"{shlex.quote(fq)} {shlex.quote(oriented_tmp)}"
    if _run_bash(cmd1, verbose) != 0:
        raise RuntimeError(f"pychopper falló en {sample}")

    # 2) filtlong (el artículo usó 1000-1700 pb)
    cmd2 = (
        f"filtlong --min_mean_q {fl_min_mean_q} --min_length {post_minlen} --max_length {post_maxlen} "
        f"{shlex.quote(oriented_tmp)} > {shlex.quote(clean_path)}"
    )
    if _run_bash(cmd2, verbose) != 0:
        raise RuntimeError(f"filtlong (post) falló en {sample}")

    # conteos para el resumen
    metrics = {
        "in_raw": _fastq_count(fq),
        "after_pychopper": _fastq_count(oriented_tmp) if os.path.exists(oriented_tmp) else 0,
        "final": _fastq_count(clean_path) if os.path.exists(clean_path) else 0,
    }

    if not keep_intermediates:
        _safe_unlink(oriented_tmp)

    return clean_path, metrics


def run_dentrim(
    input_dir: str,
    outdir: str,
    primers: str,
    pconfig: str,
    post_minlen: int,
    post_maxlen: int,
    qscore_pc: int,
    qcut_pc: float,
    fl_min_mean_q: int,
    mapper: str,
    threads: int,
    save_reports: bool,
    keep_intermediates: bool,
    verbose: bool
) -> Dict[str, List[str]]:
    _must_exist(input_dir, "dir")
    _must_exist(primers, "file")
    _must_exist(pconfig, "file")
    if post_minlen > post_maxlen:
        sys.exit("ERROR: post_minlen > post_maxlen")
    _which_or_die("filtlong"); _which_or_die("pychopper"); _which_or_die("bash")
    _assert_edlib_available_for_pychopper()
    os.makedirs(outdir, exist_ok=True)

    patterns = ["*.fastq.gz","*.fq.gz","*.fastq","*.fq"]
    inputs: List[str] = []
    for pat in patterns: inputs.extend(glob.glob(os.path.join(input_dir, pat)))
    if not inputs: sys.exit(f"No se encontraron FASTQ/FQ en: {input_dir}")

    if verbose:
        print("== Parámetros ==")
        print(f"Input:         {input_dir}")
        print(f"Output:        {outdir}")
        print(f"Post min/max:  {post_minlen}/{post_maxlen}")
        print(f"pychopper -Q:  {qscore_pc}   (umbral de calidad por lectura para clasificación)")
        print(f"pychopper -q:  {qcut_pc}     (cutoff del clasificador)")
        print(f"filtlong Q:    {fl_min_mean_q} (min_mean_q)")
        print(f"Mapper:        {mapper}")
        print(f"Hilos:         {threads}")
        print(f"Reports:       {save_reports and keep_intermediates}\n")

    processed, failed = [], []
    summary_rows: List[Dict[str,str]] = []

    for fq in sorted(inputs):
        sample = _sample_name(fq)
        print(f"\n== Procesando: {sample} ==")
        sdir = os.path.join(outdir, sample); os.makedirs(sdir, exist_ok=True)
        try:
            final_path, metrics = process_sample_paperlike(
                fq=fq, sdir=sdir, primers=primers, pconfig=pconfig,
                qscore_pc=qscore_pc, qcut_pc=qcut_pc, mapper=mapper, threads=threads,
                post_minlen=post_minlen, post_maxlen=post_maxlen,
                fl_min_mean_q=fl_min_mean_q, save_reports=save_reports,
                keep_intermediates=keep_intermediates, verbose=verbose
            )

            if keep_intermediates:
                _write_metrics_tsv(os.path.join(sdir, f"{sample}_metrics.tsv"), metrics)

            in_raw = metrics.get("in_raw",0)
            after_pc = metrics.get("after_pychopper",0)
            final_ct = metrics.get("final",0)

            pct_after_pc = (after_pc/in_raw*100.0) if in_raw else 0.0
            pct_final    = (final_ct/in_raw*100.0)  if in_raw else 0.0

            summary_rows.append({
                "sample": sample,
                "in_raw": str(in_raw),
                "after_pychopper": str(after_pc),
                "final": str(final_ct),
                "pct_after_pychopper_vs_raw": f"{pct_after_pc:.2f}",
                "pct_final_vs_raw": f"{pct_final:.2f}",
                "final_path": final_path
            })
            processed.append(sample)
            print(f"✅ {sample} listo. %final_vs_raw={pct_final:.2f}%")
        except Exception as e:
            failed.append(sample)
            print(f"⚠️  {sample} falló: {e}", file=sys.stderr)

    summary_csv = os.path.join(outdir, "dentrim_summary.csv")
    fieldnames = ["sample","in_raw","after_pychopper","final",
                  "pct_after_pychopper_vs_raw","pct_final_vs_raw","final_path"]
    with open(summary_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames); w.writeheader()
        for r in summary_rows: w.writerow(r)

    print(f"\n📄 Resumen global: {summary_csv}")
    if failed: print("⚠️  Muestras con error:", ", ".join(failed))
    return {"processed": processed, "failed": failed, "outdir": outdir}


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="pychopper → filtlong (paper-like). Solo *_clean.fastq por muestra + dentrim_summary.csv."
    )
    p.add_argument("-i","--input", required=True, help="Carpeta con FASTQ/FQ (.gz opcional)")
    p.add_argument("-o","--outdir", required=True, help="Carpeta de salida")
    p.add_argument("--primers", required=True, help="FASTA de primers (pychopper -b)")
    p.add_argument("--pconfig", required=True, help="Config de primers (TXT de pychopper -c)")

    # defaults como en el artículo
    p.add_argument("--post-minlen", type=int, default=1000, help="MinLen final (filtlong)")
    p.add_argument("--post-maxlen", type=int, default=1700, help="MaxLen final (filtlong)")
    p.add_argument("-Q","--qscore", type=int, default=12, help="QScore mínimo de pychopper (-Q)")
    p.add_argument("-q","--qcut", type=float, default=0.52, help="Cutoff del clasificador de pychopper (-q)")
    p.add_argument("--fl-min-mean-q", type=int, default=12, help="Calidad media mínima para filtlong (min_mean_q)")
    p.add_argument("-m","--mapper", default="edlib", choices=["edlib","hmmer"], help="Motor de mapeo de primers")
    p.add_argument("-t","--threads", type=int, default=8, help="Hilos para pychopper")

    p.add_argument("--no-reports", action="store_true", help="(Solo si usas --keep-intermediates) No guardar reportes de pychopper")
    p.add_argument("--keep-intermediates", action="store_true", help="Conservar orientado/buckets/reportes y métricas por muestra")
    p.add_argument("-v","--verbose", action="store_true", help="Imprimir comandos")
    return p

def main() -> None:
    args = _build_parser().parse_args()
    save_reports = not args.no_reports
    run_dentrim(
        input_dir=args.input,
        outdir=args.outdir,
        primers=args.primers,
        pconfig=args.pconfig,
        post_minlen=args.post_minlen,
        post_maxlen=args.post_maxlen,
        qscore_pc=args.qscore,
        qcut_pc=args.qcut,
        fl_min_mean_q=args.fl_min_mean_q,
        mapper=args.mapper,
        threads=args.threads,
        save_reports=save_reports,
        keep_intermediates=args.keep_intermediates,
        verbose=args.verbose
    )

if __name__ == "__main__":
    main()

