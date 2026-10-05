# Deja solo los taxones elegidos (top 4) en un CSV Source/Target/Weight que ya tiene nombres legibles

import pandas as pd
import argparse

# minúsculas y corrige el typo weisella -> weissella
def normalizar(nombre):
    n = nombre.strip().lower()
    n = n.replace("weisella", "weissella")
    return n

def main():
    parser = argparse.ArgumentParser(description="Filtra un CSV Source/Target/Weight de Gephi a solo ciertos taxones bacterianos.")
    parser.add_argument('--input', type=str, required=True, help="CSV de entrada (Source, Target, Weight) con nombres legibles.")
    parser.add_argument('--output', type=str, required=True, help="CSV de salida filtrado.")
    parser.add_argument('--taxones', type=str, required=True,
                         help="Lista de taxones separados por coma, ej: 'Lactococcus lactis,Weissella soli,Weissella oryzae,Lactobacillus coryniformis'")
    args = parser.parse_args()

    taxones_normalizados = set(normalizar(t) for t in args.taxones.split(','))

    df = pd.read_csv(args.input)
    df_filtrado = df[df['Target'].apply(normalizar).isin(taxones_normalizados)].copy()

    df_filtrado.to_csv(args.output, index=False)

    print(f"Filtrado {args.input} -> {args.output}")
    print(f"  {len(df)} aristas totales -> {len(df_filtrado)} aristas retenidas")
    print(f"  Taxones encontrados en los datos: {sorted(df_filtrado['Target'].unique())}")

    faltantes = taxones_normalizados - set(df_filtrado['Target'].apply(normalizar).unique())
    if faltantes:
        print(f"  [!] Taxones pedidos que NO aparecieron en esta combinación: {sorted(faltantes)}")

if __name__ == "__main__":
    main()
