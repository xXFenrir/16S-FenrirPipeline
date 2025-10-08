#!/usr/bin/env python3
"""
dentrim.py — Filtra y orienta lecturas 16S (filtlong + pychopper) y produce SOLO el FASTQ final.

Provee:
  - run_dentrim(input_dir, outdir, minlen, maxlen, qscore, threads, primers, pconfig, verbose=True)
  - CLI si se ejecuta como script.

Requisitos en PATH: filtlong, pychopper, bash
"""

import argparse
import glob
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List


def _must_exist(path: str, kind: str) -> None:
    if kind == "dir" and not os.path.isdir(path):
        sys.exit(f"ERROR: No existe directorio: {path}")
    if kind == "file" and not os.path.isfile(path):
        sys.exit(f"ERROR: No existe archivo: {path}")


def _which_or_die(cmd: str) -> None:
    if shutil.which(cmd) is None:
        sys.exit(f"ERROR: No se encontró '{cmd}' en PATH")


def _sample_name(f: str) -> str:
    name = Path(f).name
    for suf in (".gz", ".fastq", ".fq"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    return name


def run_dentrim(
    input_dir: str,
    outdir: str,
    minlen: int,
    maxlen: int,
    qscore: int,
    threads: int,
    primers: str,
    pconfig: str,
    verbose: bool = True,
) -> Dict[str, List[str]]:
    """
    Ejecuta el pipeline sobre todos los FASTQ/FQ del directorio.

    Devuelve:
      {
        "processed": [lista de muestras procesadas],
        "skipped":   [lista de muestras saltadas],
        "outdir":    ruta de salida (str)
      }
    """
    # Validaciones de entrada
    _must_exist(input_dir, "dir")
    _must_exist(primers, "file")
    _must_exist(pconfig, "file")
    if minlen > maxlen:
        sys.exit(f"ERROR: minlen ({minlen}) > maxlen ({maxlen})")

    # Dependencias
    _which_or_die("filtlong")
    _which_or_die("pychopper")
    _which_or_die("bash")

    os.makedirs(outdir, exist_ok=True)

    # Recolectar entradas
    patterns = ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]
    inputs: List[str] = []
    for pat in patterns:
        inputs.extend(glob.glob(os.path.join(input_dir, pat)))
    if not inputs:
        sys.exit(f"No se encontraron FASTQ en: {input_dir}")

    if verbose:
        print("== Parámetros ==")
        print(f"Input:   {input_dir}")
        print(f"Output:  {outdir}")
        print(f"MinLen:  {minlen}")
        print(f"MaxLen:  {maxlen}")
        print(f"QScore:  {qscore}")
        print(f"Threads: {threads}")
        print(f"Primers: {primers}")
        print(f"PConfig: {pconfig}\n")

    processed: List[str] = []
    skipped: List[str] = []

    for fq in inputs:
        sample = _sample_name(fq)
        sdir = os.path.join(outdir, sample)
        os.makedirs(sdir, exist_ok=True)
        out_fastq = os.path.join(sdir, f"{sample}_oriented_trimmed.fastq")
        if verbose:
            print(f"Procesando {sample} …")

        # Quoted paths (maneja espacios)
        fq_q = shlex.quote(fq)
        primers_q = shlex.quote(primers)
        pconfig_q = shlex.quote(pconfig)
        out_q = shlex.quote(out_fastq)

        # filtlong lee el archivo directamente; pychopper desde stdin.
        cmd = (
            "set -o pipefail; "
            f"filtlong --min_length {minlen} --max_length {maxlen} {fq_q}"
            f" | pychopper -m edlib -b {primers_q} -c {pconfig_q}"
            f"     -Q {qscore} -z {minlen} -t {threads} -Y 0 -q 0.52"
            f"     - {out_q}"
        )
        proc = subprocess.run(cmd, shell=True, executable="/bin/bash")

        if proc.returncode != 0:
            # Pudo quedar en 0 lecturas o error interno en pychopper; limpiamos y marcamos como saltado
            try:
                os.remove(out_fastq)
            except FileNotFoundError:
                pass
            skipped.append(sample)
            if verbose:
                print(f"⚠️  {sample}: sin lecturas tras filtrado u otro error. Saltando.", file=sys.stderr)
            continue

        processed.append(sample)

    if verbose:
        print(f"\n✅ Listo. Resultados en: {outdir}")
        if skipped:
            print("⚠️  Muestras saltadas:", ", ".join(skipped))

    return {"processed": processed, "skipped": skipped, "outdir": outdir}


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Filtra y orienta lecturas 16S (filtlong+pychopper) y deja SOLO el FASTQ final."
    )
    p.add_argument("-i", "--input", required=True, help="Carpeta con FASTQ/FQ (.gz opcional)")
    p.add_argument("-o", "--outdir", required=True, help="Carpeta de salida")
    p.add_argument("--minlen", required=True, type=int, help="Longitud mínima (filtlong, pychopper -z)")
    p.add_argument("--maxlen", required=True, type=int, help="Longitud máxima (filtlong)")
    p.add_argument("-Q", "--qscore", required=True, type=int, help="QScore mínimo (pychopper -Q)")
    p.add_argument("-t", "--threads", required=True, type=int, help="Núcleos/hilos (pychopper -t)")
    p.add_argument("--primers", required=True, help="FASTA de primers")
    p.add_argument("--pconfig", required=True, help="Config de primers (TXT de pychopper)")
    return p


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    run_dentrim(
        input_dir=args.input,
        outdir=args.outdir,
        minlen=args.minlen,
        maxlen=args.maxlen,
        qscore=args.qscore,
        threads=args.threads,
        primers=args.primers,
        pconfig=args.pconfig,
        verbose=True,
    )


if __name__ == "__main__":
    main()
