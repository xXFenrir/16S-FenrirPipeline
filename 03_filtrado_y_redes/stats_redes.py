import pandas as pd
import os

BASE = "/home/fenrir/Documentos/Tesis/Algoritmo"
CMD = os.path.join(BASE, "con mis datos")

COMBOS = ["HAC/todo", "HAC/results", "SUP/todo", "SUP/results"]
REDES = {
    "key_gene": "interacciones_key_gene.csv",
    "wgs": "interacciones_wgs.csv",
    "comp": "interacciones_comp.csv",
}


def analizar_red(edges_path, universo_size):
    df = pd.read_csv(edges_path)
    n_aristas = len(df)
    fagos = df['Source'].value_counts()
    bacterias = df['Target'].value_counts()
    n_fagos = fagos.shape[0]
    n_bacterias = bacterias.shape[0]
    n_nodos = n_fagos + n_bacterias
    densidad_pct = (n_aristas / universo_size * 100) if universo_size else 0.0
    grado_prom = (2 * n_aristas / n_nodos) if n_nodos else 0.0

    if n_fagos:
        fago_max_id, fago_max_grado = fagos.idxmax(), fagos.max()
    else:
        fago_max_id, fago_max_grado = None, 0

    if n_bacterias:
        bact_max_id, bact_max_grado = bacterias.idxmax(), bacterias.max()
    else:
        bact_max_id, bact_max_grado = None, 0

    return {
        'n_fagos': n_fagos,
        'n_bacterias': n_bacterias,
        'n_nodos': n_nodos,
        'n_aristas': n_aristas,
        'densidad_pct': round(densidad_pct, 4),
        'grado_prom': round(grado_prom, 2),
        'fago_max_grado_id': fago_max_id,
        'fago_max_grado': fago_max_grado,
        'bact_max_grado_id': bact_max_id,
        'bact_max_grado': bact_max_grado,
    }


def main():
    filas = []
    for combo in COMBOS:
        combo_dir = os.path.join(CMD, combo, "output")
        result_csv = os.path.join(combo_dir, "result.csv")
        try:
            universo = pd.read_csv(result_csv, usecols=['phage'])
            universo_size = len(universo)
        except FileNotFoundError:
            print(f"[AVISO] No se encontro {result_csv}, se omite {combo}")
            continue

        red_dir = os.path.join(combo_dir, "red_gephi")
        for red_nombre, archivo in REDES.items():
            edges_path = os.path.join(red_dir, archivo)
            if not os.path.exists(edges_path):
                print(f"[AVISO] No se encontro {edges_path}, se omite")
                continue

            stats = analizar_red(edges_path, universo_size)
            stats['combo'] = combo
            stats['red'] = red_nombre
            filas.append(stats)

            print(f"=== {combo} / {red_nombre} ===")
            print(f"  n_fagos={stats['n_fagos']}  n_bacterias={stats['n_bacterias']}  n_nodos={stats['n_nodos']}")
            print(f"  n_aristas={stats['n_aristas']}  densidad_pct={stats['densidad_pct']}%  grado_prom={stats['grado_prom']}")
            print(f"  fago_max_grado={stats['fago_max_grado_id']} ({stats['fago_max_grado']} interacciones)")
            print(f"  bact_max_grado={stats['bact_max_grado_id']} ({stats['bact_max_grado']} interacciones)")
            print()

    if filas:
        resumen = pd.DataFrame(filas)
        cols = ['combo', 'red', 'n_fagos', 'n_bacterias', 'n_nodos', 'n_aristas',
                'densidad_pct', 'grado_prom', 'fago_max_grado_id', 'fago_max_grado',
                'bact_max_grado_id', 'bact_max_grado']
        resumen = resumen[cols]
        out_path = os.path.join(CMD, "resumen_estadisticas_redes.csv")
        resumen.to_csv(out_path, index=False)
        print(f"Resumen guardado en: {out_path}")
    else:
        print("No se generaron estadisticas (no se encontraron archivos de red).")


if __name__ == '__main__':
    main()
