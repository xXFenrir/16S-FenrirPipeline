#!/usr/bin/env bash
# recolectar_anexos.sh — Copia desde el equipo de análisis hacia este repositorio los archivos
# que completan los anexos de la tesis (scripts finales, entornos, tablas, figuras y reportes).
#
# Uso (desde la carpeta del repositorio clonado, en el PC donde se hizo el análisis):
#   bash herramientas/recolectar_anexos.sh [opciones]
#
# Opciones:
#   --tesis RUTA   Carpeta de la tesis (por defecto: ~/Documentos/Tesis)
#   --max-mb N     Tamaño máximo por archivo, en MB (por defecto: 20)
#   --simular      Solo muestra lo que haría; no crea ni copia nada
#   -h, --help     Muestra esta ayuda
#
# Garantías:
#   - Fuera del repositorio solo LEE. No borra, mueve ni modifica nada del equipo.
#   - Nunca copia datos crudos (POD5, FAST5, BAM, SAM, FASTQ) ni archivos mayores que --max-mb.
#   - No hace commit ni push: al terminar se revisa con `git status`.
#   - Escribe el resumen en anexos/INVENTARIO_RECOLECCION.md.

set -uo pipefail
export LC_ALL=C   # salida de lscpu, free, sort, etc. en formato estable

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESIS="$HOME/Documentos/Tesis"
MAX_MB=20
SIMULAR=0

ayuda() { sed -n '2,19p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while (( $# )); do
  case "$1" in
    --tesis)   TESIS="${2:?Falta la ruta después de --tesis}"; shift 2 ;;
    --max-mb)  MAX_MB="${2:?Falta el número después de --max-mb}"; shift 2 ;;
    --simular) SIMULAR=1; shift ;;
    -h|--help) ayuda; exit 0 ;;
    *) echo "Opción desconocida: $1" >&2; ayuda >&2; exit 2 ;;
  esac
done

if [[ ! -d "$REPO/anexos" || ! -f "$REPO/README.md" ]]; then
  echo "No encuentro el repositorio. Ejecuta el script desde el repositorio clonado: bash herramientas/recolectar_anexos.sh" >&2
  exit 1
fi
if [[ ! -d "$TESIS" ]]; then
  echo "No existe la carpeta de la tesis: $TESIS  (usa --tesis /ruta/a/Tesis)" >&2
  exit 1
fi
if ! [[ "$MAX_MB" =~ ^[0-9]+$ ]]; then
  echo "--max-mb debe ser un número entero" >&2
  exit 2
fi

GULUPA="$TESIS/datos_gulupa"
QS8="$GULUPA/data_gulupa_qs8"
DATOS="$REPO/anexos/datos"
INVENTARIO="$REPO/anexos/INVENTARIO_RECOLECCION.md"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
touch "$TMP"/{copiados,faltantes,omitidos,avisos,scripts,verificaciones,otros_scripts,corrida}

SCRIPTS_FINALES=(stats_fastq.py dentrim_bam.py EMU_propio.py rebuild_counts.py rarefaccion.py diversidad_mod.py)
declare -A SCRIPT_SRC=()

# ---------------------------------------------------------------- utilidades

paso()      { printf '\n== %s\n' "$*" >&2; }
registrar() { printf '%s\n' "$2" >> "$TMP/$1"; }
corta()     { local p="$1"; printf '%s' "${p/#$HOME/\~}"; }          # ruta legible para el inventario
en_repo()   { printf '%s' "${1#"$REPO"/}"; }
tam_mb()    { local b; b=$(stat -c %s "$1" 2>/dev/null || echo 0); echo $(( (b + 1048575) / 1048576 )); }
crear_dir() { (( SIMULAR )) || mkdir -p "$1"; }

es_dato_crudo() {
  local n="${1,,}"
  case "$n" in
    *.pod5|*.fast5|*.bam|*.bai|*.sam|*.cram|*.fastq|*.fq|*.fastq.gz|*.fq.gz) return 0 ;;
  esac
  return 1
}

# copiar ORIGEN DESTINO — copia un archivo respetando las reglas de tamaño y tipo
copiar() {
  local src="$1" dst="$2" mb
  if [[ ! -f "$src" ]]; then registrar faltantes "- \`$(corta "$src")\`"; return 1; fi
  if es_dato_crudo "$src"; then registrar omitidos "- \`$(corta "$src")\` (datos de secuenciación)"; return 1; fi
  mb=$(tam_mb "$src")
  if (( mb > MAX_MB )); then registrar omitidos "- \`$(corta "$src")\` (${mb} MB, límite ${MAX_MB} MB)"; return 1; fi
  if (( SIMULAR )); then registrar copiados "- \`$(corta "$src")\` → \`$(en_repo "$dst")\`"; return 0; fi
  mkdir -p "$(dirname "$dst")"
  if cp -p "$src" "$dst"; then
    registrar copiados "- \`$(corta "$src")\` → \`$(en_repo "$dst")\`"
  else
    registrar avisos "- No se pudo copiar \`$(corta "$src")\`"
    return 1
  fi
}

# copiar_arbol ORIGEN DESTINO PROFUNDIDAD EXCLUIR PATRON... — copia conservando rutas relativas
copiar_arbol() {
  local src="$1" dst="$2" prof="$3" excluir="$4"; shift 4
  if [[ ! -d "$src" ]]; then registrar faltantes "- Carpeta \`$(corta "$src")\`"; return 0; fi
  local -a filtro=() ; local p n=0 f
  for p in "$@"; do (( ${#filtro[@]} )) && filtro+=(-o); filtro+=(-iname "$p"); done
  local -a no=(-not -path "$REPO/*")
  [[ -n "$excluir" ]] && no+=(-not -path "$excluir")
  while IFS= read -r -d '' f; do
    copiar "$f" "$dst/${f#"$src"/}" && n=$((n + 1))
  done < <(find "$src" -maxdepth "$prof" -type f \( "${filtro[@]}" \) "${no[@]}" -print0 2>/dev/null | sort -z)
  (( n == 0 )) && registrar faltantes "- Ningún archivo \`$*\` en \`$(corta "$src")\`"
  return 0
}

# ---------------------------------------------------------------- 1. scripts finales

recolectar_scripts() {
  paso "Scripts finales"
  local -a expr=() ; local s
  for s in "${SCRIPTS_FINALES[@]}"; do (( ${#expr[@]} )) && expr+=(-o); expr+=(-name "$s"); done
  local -a raices=("$TESIS")
  [[ -d "$HOME/scriptsbioinf" ]] && raices+=("$HOME/scriptsbioinf")
  find "${raices[@]}" -type f \( "${expr[@]}" \) -not -path "$REPO/*" -not -path '*/.git/*' \
       -printf '%T@\t%p\n' 2>/dev/null | sort -t$'\t' -k1,1nr > "$TMP/scripts_encontrados"

  for s in "${SCRIPTS_FINALES[@]}"; do
    local -a rutas=() ; local ruta
    while IFS=$'\t' read -r _ ruta; do
      [[ "$(basename "$ruta")" == "$s" ]] && rutas+=("$ruta")
    done < "$TMP/scripts_encontrados"
    if (( ${#rutas[@]} == 0 )); then
      registrar scripts "| \`$s\` | ❌ no encontrado | |"
      registrar faltantes "- Script \`$s\`"
      continue
    fi
    SCRIPT_SRC[$s]="${rutas[0]}"                       # el más reciente
    copiar "${rutas[0]}" "$REPO/scripts/$s"
    local nota="" distintos
    distintos=$(for ruta in "${rutas[@]}"; do md5sum "$ruta" | cut -c1-32; done | sort -u | wc -l)
    if (( ${#rutas[@]} > 1 )); then
      nota="${#rutas[@]} copias, ${distintos} versión(es) distinta(s)"
      if (( distintos > 1 )); then
        registrar avisos "- Hay ${distintos} versiones distintas de \`$s\`; se copió la más reciente. Confirmar que es la que produjo los resultados:"
        for ruta in "${rutas[@]}"; do
          registrar avisos "  - \`$(corta "$ruta")\` — $(date -r "$ruta" '+%Y-%m-%d %H:%M'), md5 $(md5sum "$ruta" | cut -c1-8)"
        done
      fi
    fi
    registrar scripts "| \`$s\` | \`$(corta "${rutas[0]}")\` ($(date -r "${rutas[0]}" '+%Y-%m-%d')) | ${nota} |"
  done

  # Otros scripts de análisis que podrían ser relevantes (solo se listan)
  find "$TESIS" -maxdepth 6 -type f \( -name '*.py' -o -name '*.sh' -o -name '*.R' -o -name '*.Rmd' -o -name '*.ipynb' \) \
       -not -path "$REPO/*" -not -path '*/.git/*' -not -path '*/site-packages/*' -not -path '*/envs/*' -print0 2>/dev/null \
    | sort -z | while IFS= read -r -d '' f; do
        b="$(basename "$f")"
        [[ " ${SCRIPTS_FINALES[*]} " == *" $b "* ]] && continue
        printf -- '- `%s`\n' "$(corta "$f")"
      done | head -n 100 >> "$TMP/otros_scripts"
}

# ---------------------------------------------------------------- 2. configuración y metadatos

recolectar_config() {
  paso "Primers y configuración"
  copiar "$QS8/Limpieza/primers_gulupa/primers.fasta"        "$REPO/config/primers.fasta"
  copiar "$QS8/Limpieza/primers_gulupa/primers_config.txt"   "$REPO/config/primers_config.txt"
  [[ -f "$GULUPA/data_gulupa/primers_seq/primers_seq_gulupa.fasta" ]] &&
    copiar "$GULUPA/data_gulupa/primers_seq/primers_seq_gulupa.fasta" "$REPO/config/primers_seq_gulupa.fasta"
  local f
  for f in "$HOME/resources/primers_stesen.fasta" "$HOME/resources/primers_stesen.txt"; do
    [[ -f "$f" ]] && copiar "$f" "$REPO/config/piloto/$(basename "$f")"
  done

  paso "Metadatos"
  copiar "$GULUPA/Mapa Barcodes Microbioma.csv" "$DATOS/metadatos/Mapa Barcodes Microbioma.csv"
  local n=0
  while IFS= read -r -d '' f; do
    copiar "$f" "$DATOS/metadatos/$(basename "$f")" && n=$((n + 1))
  done < <(find "$GULUPA" "$TESIS" -maxdepth 3 -type f \( -iname '*metadat*' -o -iname '*finca*' \) \
             \( -iname '*.xlsx' -o -iname '*.xls' -o -iname '*.csv' -o -iname '*.tsv' \) -not -path "$REPO/*" -print0 2>/dev/null | sort -zu)
  (( n == 0 )) && registrar faltantes "- Excel de metadatos (ningún archivo con 'metadat' o 'finca' en el nombre en \`$(corta "$GULUPA")\`)"
  # Otras hojas de cálculo que podrían ser los metadatos: solo se listan
  find "$GULUPA" -maxdepth 2 -type f \( -iname '*.xlsx' -o -iname '*.xls' -o -iname '*.ods' \) -not -iname '*metadat*' -not -iname '*finca*' \
       -not -path "$REPO/*" -print0 2>/dev/null | sort -z | while IFS= read -r -d '' f; do
    registrar avisos "- Hoja de cálculo no copiada (¿son los metadatos?): \`$(corta "$f")\`"
  done
  registrar avisos "- **Privacidad:** revisar \`anexos/datos/metadatos/\` antes de hacer público el repositorio (nombres de propietarios, teléfonos, coordenadas exactas)."
}

# ---------------------------------------------------------------- 3. corrida (MinKNOW + Dorado)

resumir_sequencing_summary() {  # ARCHIVO SALIDA
  local cuerpo
  cuerpo=$(awk -F'\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
      L = h["sequence_length_template"]; Q = h["mean_qscore_template"]; P = h["passes_filtering"]
      B = ("barcode_arrangement" in h) ? h["barcode_arrangement"] : (("barcode" in h) ? h["barcode"] : (("alias" in h) ? h["alias"] : 0))
      if (!L || !Q) exit 3
      next
    }
    {
      b = B ? $B : "sin_barcode"; n[b]++; bases[b] += $L; q[b] += $Q
      if (P && $P ~ /^(TRUE|True|true|1)$/) pass[b]++
    }
    END {
      for (b in n) printf "%s\t%d\t%s\t%d\t%.1f\t%.2f\n", b, n[b], (P ? pass[b] + 0 : "NA"), bases[b], bases[b] / n[b], q[b] / n[b]
    }
  ' "$1" | sort -V) || return 1
  [[ -n "$cuerpo" ]] || return 1
  if (( SIMULAR )); then return 0; fi
  mkdir -p "$(dirname "$2")"
  { printf 'barcode\tlecturas\tlecturas_pass\tbases\tlongitud_media\tqscore_medio_lecturas\n'; printf '%s\n' "$cuerpo"; } > "$2"
}

recolectar_corrida() {
  paso "Reportes de la corrida"
  local f dir etiqueta salida n=0
  while IFS= read -r -d '' f; do
    dir="$(basename "$(dirname "$f")")"
    copiar "$f" "$DATOS/corrida/$dir/$(basename "$f")" && n=$((n + 1))
    if [[ "$(basename "$f")" == final_summary_* ]]; then
      {
        printf -- '- `%s`\n' "$(corta "$f")"
        grep -E '^(instrument|flow_cell_id|flow_cell_product_code|protocol|sample_id|protocol_group_id|started|acquisition_stopped|processing_stopped|basecalling_enabled|pod5_files_in_final_dest|bam_files_in_final_dest)=' "$f" \
          | sed 's/^/  - `/; s/$/`/'
      } >> "$TMP/corrida"
    fi
  done < <(find "$GULUPA" -maxdepth 6 -type f \( -name 'final_summary_*' -o -name 'report_*.html' -o -name 'report_*.json' \
             -o -name 'report_*.pdf' -o -name 'barcode_alignment_*' -o -name 'throughput_*' -o -name 'pore_activity_*' -o -name 'sample_sheet_*' \) \
             -not -path "$REPO/*" -print0 2>/dev/null | sort -z)
  (( n == 0 )) && registrar faltantes "- Reportes de MinKNOW (\`final_summary_*\`, \`report_*\`) en \`$(corta "$GULUPA")\`"

  paso "Resumen de sequencing_summary por barcode"
  n=0
  while IFS= read -r -d '' f; do
    etiqueta="$(dirname "${f#"$GULUPA"/}")"; etiqueta="${etiqueta//\//_}"
    salida="$DATOS/corrida/lecturas_por_barcode_${etiqueta}.tsv"
    if resumir_sequencing_summary "$f" "$salida"; then
      registrar copiados "- Resumen por barcode de \`$(corta "$f")\` → \`$(en_repo "$salida")\`"
      n=$((n + 1))
    else
      registrar avisos "- No se pudo resumir \`$(corta "$f")\` (columnas inesperadas)"
    fi
  done < <(find "$GULUPA" -maxdepth 7 -type f -name 'sequencing_summary*.txt' -not -path "$REPO/*" -print0 2>/dev/null | sort -z)
  (( n == 0 )) && registrar faltantes "- \`sequencing_summary*.txt\` en \`$(corta "$GULUPA")\`"
}

# ---------------------------------------------------------------- 4. resultados

recolectar_resultados() {
  paso "Limpieza y control de calidad"
  local run d
  for run in hac sup; do
    d="$QS8/Limpieza/${run}_8_trim_edlib"
    [[ -d "$d" ]] || { [[ "$run" == hac ]] && registrar faltantes "- Carpeta \`$(corta "$d")\`"; continue; }
    copiar_arbol "$d" "$DATOS/calidad/limpieza_${run}" 1 "" '*.tsv' '*.csv' '*.txt' '*.xlsx'
  done
  local n=0 f
  while IFS= read -r -d '' f; do
    copiar "$f" "$DATOS/calidad/$(basename "$f")" && n=$((n + 1))
  done < <(find "$TESIS" -type f \( -name 'estadisticas_*.txt' -o -name 'estadisticas_*.xlsx' -o -name 'estadisticas_*.tsv' \) \
             -not -path "$REPO/*" -print0 2>/dev/null | sort -z)
  (( n == 0 )) && registrar faltantes "- Tablas de \`stats_fastq.py\` (\`estadisticas_*.txt/.xlsx\`) en \`$(corta "$TESIS")\`"

  paso "Taxonomía"
  for run in hac sup; do
    d="$QS8/EMU_propio/EMU${run}_results"
    [[ -d "$d" ]] || { [[ "$run" == hac ]] && registrar faltantes "- Carpeta \`$(corta "$d")\`"; continue; }
    copiar_arbol "$d" "$DATOS/taxonomia/$run" 1 "" '*.tsv' '*.csv' '*.txt' '*.xlsx'
  done
  copiar_arbol "$QS8/EMU_propio" "$DATOS/taxonomia/figuras" 5 "$QS8/EMU_propio/Diversidad/*" '*.png' '*.pdf' '*.svg'

  paso "Diversidad"
  copiar_arbol "$QS8/EMU_propio/Diversidad" "$DATOS/diversidad" 4 "" \
    '*.tsv' '*.csv' '*.txt' '*.xlsx' '*.png' '*.pdf' '*.svg' '*.html'

  paso "Fase piloto"
  [[ -d "$HOME/og_stats" ]] && copiar_arbol "$HOME/og_stats" "$DATOS/piloto/og_stats" 1 "" '*.txt' '*.xlsx' '*.tsv'
  [[ -d "$HOME/results_dentrim_Q12/EMU_Q12" ]] &&
    copiar_arbol "$HOME/results_dentrim_Q12/EMU_Q12" "$DATOS/piloto/EMU_Q12" 2 "" '*.tsv' '*.png' '*.html' '*.txt'
  [[ -f "$HOME/metricas/metadata_Q12.tsv" ]] && copiar "$HOME/metricas/metadata_Q12.tsv" "$DATOS/piloto/metadata_Q12.tsv"
  return 0
}

# ---------------------------------------------------------------- 5. entorno

encontrar_conda() {
  local c
  for c in "${CONDA_EXE:-}" "$(command -v conda 2>/dev/null)" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
           "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda"; do
    [[ -n "$c" && -x "$c" ]] && { printf '%s' "$c"; return 0; }
  done
  return 1
}

recolectar_entorno() {
  paso "Entorno computacional"
  local out="$DATOS/entorno"
  crear_dir "$out"
  local sistema
  sistema=$(
    echo "Fecha de recolección: $(date '+%Y-%m-%d %H:%M')"
    ( [[ -r /etc/os-release ]] && . /etc/os-release && echo "Sistema operativo: ${PRETTY_NAME:-desconocido}" )
    echo "Kernel: $(uname -sr)"
    echo "Equipo: $(cat /sys/class/dmi/id/sys_vendor 2>/dev/null) $(cat /sys/class/dmi/id/product_name 2>/dev/null)"
    echo "Procesador: $(lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -n 1)"
    echo "Hilos de CPU: $(nproc 2>/dev/null)"
    echo "Memoria RAM: $(free -h 2>/dev/null | awk '/^Mem:/ {print $2}')"
    if command -v nvidia-smi >/dev/null 2>&1; then
      echo "GPU: $(nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null | paste -sd ';' -)"
      echo "CUDA (nvidia-smi): $(nvidia-smi 2>/dev/null | sed -n 's/.*CUDA Version: *\([0-9.]*\).*/\1/p' | head -n 1)"
    else
      echo "GPU: nvidia-smi no disponible"
    fi
  )
  if (( SIMULAR )); then registrar copiados "- Datos del equipo → \`$(en_repo "$out/sistema.txt")\`"
  else printf '%s\n' "$sistema" > "$out/sistema.txt"; registrar copiados "- Datos del equipo → \`$(en_repo "$out/sistema.txt")\`"; fi

  # Versiones: Dorado (binario) + paquetes clave de cada entorno conda
  local tsv="$TMP/versiones.tsv" bin v
  printf 'entorno\tpaquete\tversion\tcanal\n' > "$tsv"
  while IFS= read -r -d '' bin; do
    v=$(timeout 30 "$bin" --version 2>&1 | head -n 1)
    printf '(binario %s)\tdorado\t%s\tONT\n' "$(corta "$bin")" "$v" >> "$tsv"
  done < <( { command -v dorado 2>/dev/null | tr '\n' '\0'
              find "$HOME" -maxdepth 3 -type f -path '*dorado*/bin/dorado' -print0 2>/dev/null; } | sort -zu )

  local conda prefijo nombre
  local -a prefijos=()
  if conda=$(encontrar_conda); then
    while IFS= read -r prefijo; do [[ -d "$prefijo/conda-meta" ]] && prefijos+=("$prefijo"); done \
      < <("$conda" env list 2>/dev/null | awk '!/^#/ && NF {print $NF}')
  else
    registrar avisos "- No se encontró \`conda\`: no se exportaron entornos ni versiones de paquetes."
  fi
  local claves='^(python|pychopper|filtlong|samtools|edlib|parasail-python|emu|minimap2|pandas|numpy|scipy|scikit-bio|seaborn|matplotlib|openpyxl|pod5|osfclient|plotly|networkx|biom-format|pysam)$'
  for prefijo in "${prefijos[@]}"; do
    nombre="$(basename "$prefijo")"
    "$conda" list -p "$prefijo" 2>/dev/null \
      | awk -v env="$nombre" -v re="$claves" '!/^#/ && NF >= 2 && $1 ~ re { printf "%s\t%s\t%s\t%s\n", env, $1, $2, (NF >= 4 ? $4 : "defaults") }' >> "$tsv"
    if (( ! SIMULAR )); then
      mkdir -p "$out/conda_list" "$REPO/envs"
      "$conda" list -p "$prefijo" > "$out/conda_list/$nombre.txt" 2>/dev/null
      "$conda" env export -p "$prefijo" --no-builds 2>/dev/null | grep -v '^prefix:' > "$REPO/envs/$nombre.yml"
      "$conda" env export -p "$prefijo" --from-history 2>/dev/null | grep -v '^prefix:' > "$REPO/envs/${nombre}_historial.yml"
    fi
    registrar copiados "- Entorno conda \`$nombre\` → \`envs/$nombre.yml\`, \`envs/${nombre}_historial.yml\`"
  done
  if (( SIMULAR )); then registrar copiados "- Versiones de software → \`$(en_repo "$out/versiones_software.tsv")\`"
  else cp "$tsv" "$out/versiones_software.tsv"; registrar copiados "- Versiones de software → \`$(en_repo "$out/versiones_software.tsv")\`"; fi
  cp "$tsv" "$TMP/versiones_final.tsv"

  # Ayuda (--help) de cada script: documenta todos sus parámetros y valores por defecto
  local s src py ok
  for s in "${SCRIPTS_FINALES[@]}"; do
    src="${SCRIPT_SRC[$s]:-}"; [[ -n "$src" ]] || continue
    ok=""
    for py in python3 "${prefijos[@]/%//bin/python}"; do
      command -v "$py" >/dev/null 2>&1 || [[ -x "$py" ]] || continue
      if timeout 60 "$py" "$src" --help > "$TMP/ayuda.txt" 2>&1 && grep -qi 'usage' "$TMP/ayuda.txt"; then ok="$py"; break; fi
    done
    if [[ -n "$ok" ]]; then
      if (( ! SIMULAR )); then
        mkdir -p "$out/ayuda_scripts"
        { echo "# $s --help  (intérprete: $(corta "$(command -v "$ok" 2>/dev/null || echo "$ok")"))"; cat "$TMP/ayuda.txt"; } > "$out/ayuda_scripts/$s.txt"
      fi
      registrar copiados "- Ayuda de \`$s\` → \`anexos/datos/entorno/ayuda_scripts/$s.txt\` (funciona con \`$(corta "$ok")\`)"
    else
      registrar avisos "- No se pudo obtener \`--help\` de \`$s\` con ningún intérprete."
    fi
  done

  # Base de datos de EMU: lista de archivos y sumas MD5
  local db info="" f
  for db in "$TESIS/emu_db" "$HOME/emu_db"; do
    [[ -d "$db" ]] || continue
    info+="## $(corta "$db")"$'\n'"$(ls -la "$db" 2>/dev/null)"$'\n'
    for f in species_taxid.fasta taxonomy.tsv; do
      [[ -f "$db/$f" ]] && info+="md5 $(md5sum "$db/$f" | cut -c1-32)  $f"$'\n'
    done
    info+=$'\n'
  done
  if [[ -n "$info" ]]; then
    (( SIMULAR )) || printf '%s' "$info" > "$out/emu_db_info.txt"
    registrar copiados "- Información de la base de EMU → \`anexos/datos/entorno/emu_db_info.txt\`"
  else
    registrar faltantes "- Base de datos de EMU (\`$(corta "$TESIS/emu_db")\` o \`~/emu_db\`)"
  fi
}

# ---------------------------------------------------------------- 6. verificaciones

verificar() {
  paso "Verificaciones"
  local v="$TMP/verificaciones" f

  # 1. Filtlong --min_mean_q
  {
    echo "### 1. Escala de \`--fl-min-mean-q\` en Filtlong (Anexo C, verificación 1)"
    echo
    f="${SCRIPT_SRC[dentrim_bam.py]:-}"
    if [[ -n "$f" ]]; then
      echo "Líneas de \`dentrim_bam.py\` que tratan la calidad media de Filtlong:"
      echo
      echo '```'
      grep -n -i -E 'min_mean_q|min-mean-q|fl_min|mean_q' "$f" || echo "(ninguna coincidencia)"
      echo '```'
      echo
      echo "Si el valor (12) llega sin convertir a \`--min_mean_q\`, Filtlong lo interpreta como 12 % de exactitud y no filtra nada. Q12 equivale a \`--min_mean_q 93.69\`."
    else
      echo "No se encontró \`dentrim_bam.py\`."
    fi
    echo
  } >> "$v"

  # 2. Orientación: primers_config.txt + anotación strand= de Pychopper en las lecturas limpias
  {
    echo "### 2. Orientación de las lecturas (Anexo C, verificación 2)"
    echo
    f="$QS8/Limpieza/primers_gulupa/primers_config.txt"
    if [[ -f "$f" ]]; then
      local cfg; cfg="$(head -n 1 "$f" | tr -d '\r')"
      echo "Primera línea de \`primers_config.txt\`: \`$cfg\`"
      echo
      if [[ "|$cfg" == *"|-:"* ]]; then
        echo "✅ Hay una configuración marcada con \`-:\`: Pychopper reorienta las lecturas reversas."
      else
        echo "⚠️ Ninguna configuración empieza por \`-:\`: las lecturas reversas **no** se reorientan. Formato esperado: \`+:16sF,-16sR|-:16sR,-16sF\`."
      fi
    else
      echo "No se encontró \`primers_config.txt\`."
    fi
    echo
    echo "Anotación \`strand=\` de Pychopper en las primeras 50 000 lecturas de cada FASTQ limpio:"
    echo
    local total_menos=0 total_mas=0 fila mas menos tot lector filas=""
    while IFS= read -r -d '' f; do
      lector="cat"; [[ "$f" == *.gz ]] && lector="gzip -dc"
      fila=$($lector "$f" 2>/dev/null | head -n 200000 | awk 'NR % 4 == 1 { t++; if ($0 ~ /strand=\+|TS:A:\+/) p++; else if ($0 ~ /strand=-|TS:A:-/) m++ } END { printf "%d %d %d", p + 0, m + 0, t + 0 }')
      read -r mas menos tot <<< "$fila"
      total_mas=$((total_mas + mas)); total_menos=$((total_menos + menos))
      filas+="| \`$(basename "$f")\` | $mas | $menos | $tot |"$'\n'
    done < <(find "$QS8/Limpieza/hac_8_trim_edlib" -maxdepth 2 -type f \( -name '*_limpio.fastq' -o -name '*_limpio.fastq.gz' \) -print0 2>/dev/null | sort -zV)
    if [[ -n "$filas" ]]; then
      echo "| FASTQ | strand=+ | strand=- | lecturas revisadas |"
      echo "|---|---|---|---|"
      printf '%s' "$filas"
      echo
    fi
    if (( total_mas + total_menos == 0 )); then
      echo "Sin anotaciones \`strand=\` (o sin FASTQ limpios): no se puede verificar por esta vía."
    elif (( total_menos == 0 )); then
      echo "⚠️ Ninguna lectura con \`strand=-\`: las lecturas reversas no se reorientaron."
    else
      echo "✅ Hay lecturas con \`strand=-\` (${total_menos} de $((total_mas + total_menos))): la reorientación funciona."
    fi
    echo
  } >> "$v"

  # 3. Jaccard: línea del código, versión de SciPy y mediana de la matriz
  {
    echo "### 3. Distancia de Jaccard (Anexo G, verificación 3)"
    echo
    f="${SCRIPT_SRC[diversidad_mod.py]:-}"
    if [[ -n "$f" ]]; then
      echo "Líneas de \`diversidad_mod.py\` con \`jaccard\`:"
      echo
      echo '```'
      grep -n -i 'jaccard' "$f" || echo "(ninguna coincidencia)"
      echo '```'
    else
      echo "No se encontró \`diversidad_mod.py\`."
    fi
    echo
    echo "Versiones de SciPy instaladas (el cálculo es de presencia/ausencia con abundancias solo desde SciPy 1.15):"
    echo
    local scipys
    scipys=$(awk -F'\t' '$2 == "scipy" { printf "- entorno `%s`: SciPy %s\n", $1, $3 }' "$TMP/versiones_final.tsv" 2>/dev/null)
    printf '%s\n' "${scipys:-- (no se encontró SciPy en los entornos conda)}"
    echo
    echo "Resumen de las matrices de Jaccard encontradas (valores fuera de la diagonal):"
    echo
    local matrices=""
    while IFS= read -r -d '' f; do
      fila=$(python3 - "$f" <<'PY' 2>/dev/null
import csv, statistics, sys, os
ruta = sys.argv[1]
sep = "," if ruta.endswith(".csv") else "\t"
with open(ruta, newline="") as fh:
    filas = [r for r in csv.reader(fh, delimiter=sep) if r]
cuerpo = filas[1:]
if not cuerpo or any(len(r) != len(filas[0]) for r in cuerpo) or len(cuerpo) != len(filas[0]) - 1:
    sys.exit(0)  # no es una matriz cuadrada con encabezado e índice
vals = []
for i, fila in enumerate(cuerpo):
    for j, x in enumerate(fila[1:]):
        if i != j:
            try:
                vals.append(float(x))
            except ValueError:
                pass
if vals:
    print(f"| `{os.path.basename(ruta)}` | {len(vals)} | {min(vals):.3f} | {statistics.median(vals):.3f} | {max(vals):.3f} |")
PY
)
      [[ -n "$fila" ]] && matrices+="$fila"$'\n'
    done < <(find "$QS8/EMU_propio/Diversidad" -maxdepth 4 -type f -iname '*jaccard*' \( -iname '*.tsv' -o -iname '*.csv' \) -print0 2>/dev/null | sort -z)
    if [[ -n "$matrices" ]]; then
      echo "| Archivo | pares | mínimo | mediana | máximo |"
      echo "|---|---|---|---|---|"
      printf '%s' "$matrices"
      echo
      echo "Una mediana muy cercana a 1 junto con SciPy < 1.15 indica que se calculó sobre abundancias y no sobre presencia/ausencia."
    else
      echo "No se encontraron matrices de Jaccard en \`$(corta "$QS8/EMU_propio/Diversidad")\`."
    fi
    echo
  } >> "$v"
}

# ---------------------------------------------------------------- inventario

escribir_inventario() {
  local c_cop c_fal c_omi c_avi
  c_cop=$(grep -c '^- ' "$TMP/copiados"); c_fal=$(grep -c '^- ' "$TMP/faltantes")
  c_omi=$(grep -c '^- ' "$TMP/omitidos"); c_avi=$(grep -c '^- ' "$TMP/avisos")
  {
    echo "# Inventario de la recolección"
    echo
    echo "Generado por \`herramientas/recolectar_anexos.sh\` el $(date '+%Y-%m-%d %H:%M')."
    echo "Carpeta de la tesis: \`$(corta "$TESIS")\` · Límite por archivo: ${MAX_MB} MB · Modo: $( (( SIMULAR )) && echo 'simulación (no se copió nada)' || echo 'copia')"
    echo
    echo "| | Cantidad |"
    echo "|---|---|"
    echo "| Copiados o generados | $c_cop |"
    echo "| No encontrados | $c_fal |"
    echo "| Omitidos (tamaño o datos crudos) | $c_omi |"
    echo "| Avisos | $c_avi |"
    echo
    echo "## Verificaciones"
    echo
    cat "$TMP/verificaciones"
    echo "## Scripts finales"
    echo
    echo "| Script | Copiado desde | Nota |"
    echo "|---|---|---|"
    cat "$TMP/scripts"
    echo
    if [[ -s "$TMP/corrida" ]]; then
      echo "## Datos de la corrida (final_summary)"
      echo
      cat "$TMP/corrida"
      echo
    fi
    echo "## Avisos"
    echo
    if [[ -s "$TMP/avisos" ]]; then cat "$TMP/avisos"; else echo "Ninguno."; fi
    echo
    echo "## No encontrados"
    echo
    echo "Si alguno existe en otra ruta, cópialo a mano a la carpeta correspondiente."
    echo
    if [[ -s "$TMP/faltantes" ]]; then cat "$TMP/faltantes"; else echo "Ninguno."; fi
    echo
    echo "## Omitidos"
    echo
    if [[ -s "$TMP/omitidos" ]]; then cat "$TMP/omitidos"; else echo "Ninguno."; fi
    echo
    echo "## Copiados o generados"
    echo
    if [[ -s "$TMP/copiados" ]]; then cat "$TMP/copiados"; else echo "Ninguno."; fi
    echo
    if [[ -s "$TMP/otros_scripts" ]]; then
      echo "## Otros scripts en la carpeta de la tesis (no copiados)"
      echo
      echo "Revisar si alguno se usó para resultados de la tesis; si es así, copiarlo a \`scripts/\`."
      echo
      cat "$TMP/otros_scripts"
      echo
    fi
    echo "## Siguientes pasos"
    echo
    echo "1. Revisar los metadatos en \`anexos/datos/metadatos/\` (privacidad)."
    echo "2. Completar los ⏳ de los anexos con estos archivos."
    echo "3. Subir: \`git status\`, \`git add -A\`, \`git commit -m \"Agregar scripts finales, entornos y datos de los anexos\"\`, \`git push\`."
  }
}

# ---------------------------------------------------------------- principal

(( SIMULAR )) && echo "MODO SIMULACIÓN: no se va a copiar ni escribir nada." >&2
echo "Repositorio: $REPO" >&2
echo "Carpeta de la tesis: $TESIS" >&2

recolectar_scripts
recolectar_config
recolectar_corrida
recolectar_resultados
recolectar_entorno
verificar

if (( SIMULAR )); then
  escribir_inventario
else
  escribir_inventario > "$INVENTARIO"
  echo >&2
  echo "Listo. Inventario: $(en_repo "$INVENTARIO")" >&2
  echo "Revisa los cambios con: git status" >&2
fi
