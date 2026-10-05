#!/usr/bin/env python3

# Limpieza de los BAM de Dorado (una carpeta por barcode): samtools fastq -> pychopper -> filtlong

import argparse, csv, glob, math, os, shlex, shutil, subprocess, sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any


def _calcular_estadisticas_fastq(ruta_fastq: str) -> Dict[str, Any]:
    if not os.path.exists(ruta_fastq) or os.path.getsize(ruta_fastq) == 0:
        return {
            "reads": 0, "bases": 0, "mean_len": 0,
            "median_len": 0, "n50": 0, "mean_q": 0.0
        }

    longitudes = []
    bases_totales = 0
    prob_error_total = 0.0

    with open(ruta_fastq, "rt", encoding="utf-8", errors="ignore") as fh:
        while True:
            encabezado = fh.readline()
            if not encabezado:
                break
            secuencia = fh.readline().strip()
            mas = fh.readline()
            calidad = fh.readline().strip()

            if not calidad:
                break

            l = len(secuencia)
            longitudes.append(l)
            bases_totales += l

            # el QScore se saca de la prob. de error de cada base, no promediando los Phred
            for caracter in calidad:
                q = ord(caracter) - 33
                prob_error_total += 10.0 ** (-q / 10.0)

    num_lecturas = len(longitudes)
    if num_lecturas == 0 or bases_totales == 0:
        return {
            "reads": 0, "bases": 0, "mean_len": 0,
            "median_len": 0, "n50": 0, "mean_q": 0.0
        }

    longitud_media = int(round(bases_totales / num_lecturas))

    longitudes.sort()
    if num_lecturas % 2 == 1:
        longitud_mediana = longitudes[num_lecturas // 2]
    else:
        longitud_mediana = int(round((longitudes[num_lecturas // 2 - 1] + longitudes[num_lecturas // 2]) / 2.0))

    mitad_bases = bases_totales / 2.0
    bases_acumuladas = 0
    n50 = 0
    for l in longitudes:
        bases_acumuladas += l
        if bases_acumuladas >= mitad_bases:
            n50 = l
            break

    prob_error_media = prob_error_total / bases_totales
    if prob_error_media > 0:
        qscore_medio = round(-10.0 * math.log10(prob_error_media), 2)
    else:
        qscore_medio = 60.0

    return {
        "reads": num_lecturas,
        "bases": bases_totales,
        "mean_len": longitud_media,
        "median_len": longitud_mediana,
        "n50": n50,
        "mean_q": qscore_medio
    }


def _debe_existir(ruta: str, tipo: str) -> None:
    if tipo == "dir" and not os.path.isdir(ruta):
        sys.exit(f"ERROR: No existe directorio: {ruta}")
    if tipo == "file" and not os.path.isfile(ruta):
        sys.exit(f"ERROR: No existe archivo: {ruta}")

def _requerir_comando(cmd: str) -> None:
    if shutil.which(cmd) is None:
        sys.exit(f"ERROR: No se encontró '{cmd}' en PATH")

def _verificar_edlib_para_pychopper() -> None:
    prueba = subprocess.run([sys.executable, "-c", "import edlib"], capture_output=True)
    if prueba.returncode != 0:
        sys.exit(
            "ERROR: El 'pychopper' requiere el módulo 'edlib' en su intérprete.\n"
            f"       Intérprete actual: {sys.executable}\n"
            "       Instálalo corriendo: python -m pip install edlib"
        )

def _eliminar_si_existe(ruta: str) -> None:
    try: os.remove(ruta)
    except FileNotFoundError: pass

def _ejecutar_bash(cmd: str, detallado: bool = True) -> int:
    if detallado: print(f"$ {cmd}")
    return subprocess.run(cmd, shell=True, executable="/bin/bash").returncode


def procesar_muestra(
    archivo_bam: str,
    id_barcode: str,
    carpeta_muestra: str,
    cebadores: str,
    config_cebadores: str,
    qscore_pychopper: int,
    qcut_pychopper: float,
    alineador: str,
    hilos: int,
    longitud_min: int,
    longitud_max: int,
    calidad_media_min_filtlong: int,
    guardar_reportes: bool,
    conservar_intermedios: bool,
    detallado: bool
) -> Tuple[str, Dict[str, Any], Dict[str, Any], Dict[str, Any]]:

    os.makedirs(carpeta_muestra, exist_ok=True)
    fastq_crudo    = os.path.join(carpeta_muestra, f"{id_barcode}_crudo.fastq")
    orientado_temp = os.path.join(carpeta_muestra, f"{id_barcode}__orientado_tmp.fastq")
    ruta_limpio    = os.path.join(carpeta_muestra, f"{id_barcode}_limpio.fastq")

    # 0) BAM -> FASTQ
    comando_bam_a_fastq = f"samtools fastq {shlex.quote(archivo_bam)} > {shlex.quote(fastq_crudo)} 2>/dev/null || true"
    _ejecutar_bash(comando_bam_a_fastq, detallado)

    estadisticas_crudas = _calcular_estadisticas_fastq(fastq_crudo)
    if estadisticas_crudas["reads"] == 0:
        raise RuntimeError(f"No se pudieron extraer lecturas FASTQ de {id_barcode}")

    # 1) pychopper (-g con phmm, -b con edlib)
    if alineador == "phmm" or alineador == "hmmer":
        bandera_primer = f"-g {shlex.quote(cebadores)}"
        alineador_real = "phmm"
    else:
        bandera_primer = f"-b {shlex.quote(cebadores)}"
        alineador_real = "edlib"

    comando_pychopper = (
        f"pychopper -m {alineador_real} {bandera_primer} -c {shlex.quote(config_cebadores)} "
        f"-Q {qscore_pychopper} -z {longitud_min} -t {hilos} -Y 0 -q {qcut_pychopper} "
    )
    if guardar_reportes and conservar_intermedios:
        reporte_pdf = os.path.join(carpeta_muestra, f"{id_barcode}_reporte.pdf")
        estadisticas_tsv = os.path.join(carpeta_muestra, f"{id_barcode}_estadisticas.tsv")
        puntajes_tsv = os.path.join(carpeta_muestra, f"{id_barcode}_puntajes.tsv")
        comando_pychopper += (
            f"-r {shlex.quote(reporte_pdf)} -S {shlex.quote(estadisticas_tsv)} -A {shlex.quote(puntajes_tsv)} "
            f"-K {shlex.quote(os.path.join(carpeta_muestra, f'{id_barcode}_fallo_calidad.fastq'))} "
            f"-l {shlex.quote(os.path.join(carpeta_muestra, f'{id_barcode}_fallo_longitud.fastq'))} "
            f"-u {shlex.quote(os.path.join(carpeta_muestra, f'{id_barcode}_sin_clasificar.fastq'))} "
            f"-w {shlex.quote(os.path.join(carpeta_muestra, f'{id_barcode}_rescatadas.fastq'))} "
        )
    comando_pychopper += f"{shlex.quote(fastq_crudo)} {shlex.quote(orientado_temp)}"

    if _ejecutar_bash(comando_pychopper, detallado) != 0:
        _eliminar_si_existe(fastq_crudo)
        raise RuntimeError(f"pychopper falló en {id_barcode}")

    estadisticas_pychopper = _calcular_estadisticas_fastq(orientado_temp)

    # 2) filtlong
    comando_filtlong = (
        f"filtlong --min_mean_q {calidad_media_min_filtlong} --min_length {longitud_min} --max_length {longitud_max} "
        f"{shlex.quote(orientado_temp)} > {shlex.quote(ruta_limpio)}"
    )
    if _ejecutar_bash(comando_filtlong, detallado) != 0:
        _eliminar_si_existe(fastq_crudo)
        raise RuntimeError(f"filtlong falló en {id_barcode}")

    estadisticas_finales = _calcular_estadisticas_fastq(ruta_limpio)

    # borrar temporales
    _eliminar_si_existe(fastq_crudo)
    if not conservar_intermedios:
        _eliminar_si_existe(orientado_temp)

    return ruta_limpio, estadisticas_crudas, estadisticas_pychopper, estadisticas_finales


def ejecutar_limpieza(
    carpeta_entrada: str,
    carpeta_salida: str,
    cebadores: str,
    config_cebadores: str,
    max_barcodes: int,
    longitud_min: int,
    longitud_max: int,
    qscore_pychopper: int,
    qcut_pychopper: float,
    calidad_media_min_filtlong: int,
    alineador: str,
    hilos: int,
    guardar_reportes: bool,
    conservar_intermedios: bool,
    detallado: bool
) -> Dict[str, List[str]]:
    _debe_existir(carpeta_entrada, "dir")
    _debe_existir(cebadores, "file")
    _debe_existir(config_cebadores, "file")
    if longitud_min > longitud_max:
        sys.exit("ERROR: post_minlen > post_maxlen")
    _requerir_comando("samtools"); _requerir_comando("filtlong"); _requerir_comando("pychopper"); _requerir_comando("bash")
    _verificar_edlib_para_pychopper()
    os.makedirs(carpeta_salida, exist_ok=True)

    muestras_objetivo: List[Tuple[str, str]] = []
    for i in range(1, max_barcodes + 1):
        nombre_barcode = f"barcode{i:02d}"
        carpeta_barcode = os.path.join(carpeta_entrada, nombre_barcode)

        if os.path.isdir(carpeta_barcode):
            archivos_bam = glob.glob(os.path.join(carpeta_barcode, "*.bam"))
            if archivos_bam:
                muestras_objetivo.append((nombre_barcode, archivos_bam[0]))

    if not muestras_objetivo:
        sys.exit(f"No se encontraron carpetas con archivos BAM en: {carpeta_entrada}")

    if detallado:
        print("== Parámetros ==")
        print(f"Input Dir:     {carpeta_entrada}")
        print(f"Output Dir:    {carpeta_salida}")
        print(f"Muestras:      {len(muestras_objetivo)} barcodes detectados")
        print(f"Post min/max:  {longitud_min}/{longitud_max}")
        print(f"pychopper -Q:  {qscore_pychopper}")
        print(f"filtlong Q:    {calidad_media_min_filtlong}")
        print(f"Hilos:         {hilos}\n")

    procesados, fallidos = [], []
    filas_resumen: List[Dict[str, Any]] = []

    for id_barcode, ruta_bam in muestras_objetivo:
        nombre_muestra = os.path.basename(ruta_bam).replace(".bam", "")
        print(f"\n== Procesando: {id_barcode} ({nombre_muestra}.bam) ==")
        carpeta_muestra = os.path.join(carpeta_salida, id_barcode)
        try:
            ruta_final, estadisticas_crudas, estadisticas_pychopper, estadisticas_finales = procesar_muestra(
                archivo_bam=ruta_bam, id_barcode=id_barcode, carpeta_muestra=carpeta_muestra,
                cebadores=cebadores, config_cebadores=config_cebadores,
                qscore_pychopper=qscore_pychopper, qcut_pychopper=qcut_pychopper, alineador=alineador, hilos=hilos,
                longitud_min=longitud_min, longitud_max=longitud_max,
                calidad_media_min_filtlong=calidad_media_min_filtlong, guardar_reportes=guardar_reportes,
                conservar_intermedios=conservar_intermedios, detallado=detallado
            )

            # % de remoción
            lecturas_crudas   = estadisticas_crudas["reads"]
            lecturas_pychopper  = estadisticas_pychopper["reads"]
            lecturas_finales  = estadisticas_finales["reads"]

            bases_crudas = estadisticas_crudas["bases"]
            bases_finales  = estadisticas_finales["bases"]

            pct_remocion_pychopper = round(((lecturas_crudas - lecturas_pychopper) / lecturas_crudas * 100.0), 2) if lecturas_crudas else 0.0
            pct_remocion_total     = round(((lecturas_crudas - lecturas_finales) / lecturas_crudas * 100.0), 2) if lecturas_crudas else 0.0
            pct_remocion_bases     = round(((bases_crudas - bases_finales) / bases_crudas * 100.0), 2) if bases_crudas else 0.0

            filas_resumen.append({
                "barcode": id_barcode,
                "muestra": nombre_muestra,
                "lecturas raw": lecturas_crudas,
                "lecturas after pychopper": lecturas_pychopper,
                "% remocion reads p": pct_remocion_pychopper,
                "lecturas after filtlong": lecturas_finales,
                "% remocion reads p+f": pct_remocion_total,
                "bases raw": bases_crudas,
                "bases after": bases_finales,
                "% remocion bases": pct_remocion_bases,
                "longitud media raw": estadisticas_crudas["mean_len"],
                "longitud media after": estadisticas_finales["mean_len"],
                "longitud mediana raw": estadisticas_crudas["median_len"],
                "longitud mediana after": estadisticas_finales["median_len"],
                "N50 raw": estadisticas_crudas["n50"],
                "N50 after": estadisticas_finales["n50"],
                "QScore promedio raw": estadisticas_crudas["mean_q"],
                "QScore promedio after": estadisticas_finales["mean_q"]
            })
            procesados.append(id_barcode)
            print(f"✅ {id_barcode} listo. Lecturas finales: {lecturas_finales} (Remoción total: {pct_remocion_total}%)")
        except Exception as e:
            fallidos.append(id_barcode)
            print(f"⚠️  {id_barcode} falló: {e}", file=sys.stderr)

    resumen_tsv = os.path.join(carpeta_salida, "resumen_limpieza.tsv")
    columnas = [
        "barcode", "muestra", "lecturas raw", "lecturas after pychopper",
        "% remocion reads p", "lecturas after filtlong", "% remocion reads p+f",
        "bases raw", "bases after", "% remocion bases",
        "longitud media raw", "longitud media after",
        "longitud mediana raw", "longitud mediana after",
        "N50 raw", "N50 after",
        "QScore promedio raw", "QScore promedio after"
    ]
    with open(resumen_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columnas, delimiter="\t")
        w.writeheader()
        for r in filas_resumen:
            w.writerow(r)

    print(f"\n📄 Resumen global generado en: {resumen_tsv}")
    if fallidos: print("⚠️  Barcodes con error:", ", ".join(fallidos))
    return {"processed": procesados, "failed": fallidos, "outdir": carpeta_salida}


def _construir_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Pipeline BAM -> FASTQ -> pychopper -> filtlong por carpeta de barcode."
    )
    p.add_argument("-i","--input", required=True, help="Ruta hasta la carpeta que contiene las carpetas 'barcodeXX'")
    p.add_argument("-o","--outdir", required=True, help="Carpeta de salida")
    p.add_argument("--primers", required=True, help="FASTA de primers (pychopper -b)")
    p.add_argument("--pconfig", required=True, help="Config de primers (TXT de pychopper -c)")
    p.add_argument("--max-barcode", type=int, default=73, help="Número máximo de barcode a procesar (default: 73)")

    p.add_argument("--post-minlen", type=int, default=1000, help="MinLen final (filtlong)")
    p.add_argument("--post-maxlen", type=int, default=1700, help="MaxLen final (filtlong)")
    p.add_argument("-Q","--qscore", type=int, default=8, help="QScore mínimo de pychopper (-Q)")
    p.add_argument("-q","--qcut", type=float, default=0.52, help="Cutoff del clasificador de pychopper (-q)")
    p.add_argument("--fl-min-mean-q", type=int, default=12, help="Calidad media mínima para filtlong (min_mean_q)")
    p.add_argument("-m","--mapper", default="edlib", choices=["edlib","hmmer"], help="Motor de mapeo de primers")
    p.add_argument("-t","--threads", type=int, default=20, help="Hilos para pychopper")

    p.add_argument("--no-reports", action="store_true", help="No guardar reportes de pychopper")
    p.add_argument("--keep-intermediates", action="store_true", help="Conservar orientado/buckets/reportes")
    p.add_argument("-v","--verbose", action="store_true", help="Imprimir comandos")
    return p

def main() -> None:
    args = _construir_parser().parse_args()
    guardar_reportes = not args.no_reports
    ejecutar_limpieza(
        carpeta_entrada=args.input,
        carpeta_salida=args.outdir,
        cebadores=args.primers,
        config_cebadores=args.pconfig,
        max_barcodes=args.max_barcode,
        longitud_min=args.post_minlen,
        longitud_max=args.post_maxlen,
        qscore_pychopper=args.qscore,
        qcut_pychopper=args.qcut,
        calidad_media_min_filtlong=args.fl_min_mean_q,
        alineador=args.mapper,
        hilos=args.threads,
        guardar_reportes=guardar_reportes,
        conservar_intermedios=args.keep_intermediates,
        detallado=args.verbose
    )

if __name__ == "__main__":
    main()
