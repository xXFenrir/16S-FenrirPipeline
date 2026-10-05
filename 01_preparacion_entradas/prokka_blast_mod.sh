#!/bin/bash

# use getopts Parse command line arguments
while getopts ":i:o:p:d:b:O:" opt; do
  case $opt in
    i)
      input_dir="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;
    p)
      prokka="$OPTARG"
      ;;
    d)
      db="$OPTARG"
      ;;
    b)
      blast="$OPTARG"
      ;;
    O)
      out="$OPTARG"
      ;;
    \?)
      echo "invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "option -$OPTARG need a parameter." >&2
      exit 1
      ;;
  esac
done

# prokka batch annotation and blast
# con "$input_dir"/* ya no falla por los espacios en el nombre de la carpeta
for file in "$input_dir"/*; do
    # acepta .fasta y .fna
    if [[ "$file" == *.fasta || "$file" == *.fna ]]; then
        val=$(echo "${file##*/}" | cut -d '.' -f 1)
        
        # --cpus 0 usa todos los hilos
        "$prokka" "$file" --outdir "$output_dir/prokka_$val" --prefix "$val" --force --cpus 0
        
        "$blast" -query "$file" -db "$db" -outfmt 6 -max_target_seqs 1 -out "$out/$val.out" -num_threads 16
    fi
done
