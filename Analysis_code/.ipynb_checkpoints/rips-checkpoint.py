#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genera complejos de Rips y diagramas de persistencia para todos los CSV
que contienen centroides en una carpeta dada.

Ejecución
---------
$ python calcular_rips.py /ruta/a/centroides [--radio 1000] [--workers 4]

Argumentos
----------
ruta/a/centroides : carpeta con CSV (debe contener 'X_centroid', 'Y_centroid').
--radio            : max_edge_length del complejo de Rips (float, default 1000).
--workers          : núcleos a usar (int, default = todos los disponibles).

Salidas
-------
Para cada <nombre>.csv:
    <nombre>_complejo_rips.png
    <nombre>_diagrama_persistencia.png
    <nombre>.csv                (tabla birth–death)

Todos los ficheros se guardan en:
    <ruta_centroides>/resultados/rips_<radio>/
"""

import os
import sys
import time
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from multiprocessing import cpu_count

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gudhi as gd
from tqdm import tqdm


# --------------------------------------------------------------------------- #
#  FUNCIÓN QUE PROCESA UN ÚNICO CSV (se ejecuta en cada proceso)
# --------------------------------------------------------------------------- #
def _procesar_csv(nombre_csv: str, ruta_in: str, ruta_out: str, radio: float):
    """Lee un CSV, calcula Rips + persistencia y guarda resultados."""
    ruta_completa = os.path.join(ruta_in, nombre_csv)

    # --- Leer centroides -----------------------------------------------------
    df = pd.read_csv(ruta_completa)
    centroides_x = df["X_centroid"].to_numpy()
    centroides_y = df["Y_centroid"].to_numpy()
    puntos = np.column_stack((centroides_x, centroides_y))

    # --- Complejo de Rips ----------------------------------------------------
    rips_complex = gd.RipsComplex(points=puntos, max_edge_length=radio)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)

    nombre_base = os.path.splitext(nombre_csv)[0]

    # --- Imagen del complejo de Rips ----------------------------------------
    plt.figure(figsize=(24, 10))
    plt.scatter(centroides_x, centroides_y, color="black",
                label="Centroides", s=5)
    for simplex in simplex_tree.get_skeleton(1):
        if len(simplex[0]) == 2:
            i, j = simplex[0]
            plt.plot([centroides_x[i], centroides_x[j]],
                     [centroides_y[i], centroides_y[j]],
                     color="gray", linewidth=0.5)
    plt.title(f"Complejo de Rips (r={radio}) · {nombre_csv}")
    plt.xlabel("X"); plt.ylabel("Y"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(ruta_out,
                             f"{nombre_base}_complejo_rips.png"))
    plt.close()

    # --- Diagrama de persistencia -------------------------------------------
    diag = simplex_tree.persistence()
    plt.figure(figsize=(6, 6))
    gd.plot_persistence_diagram(diag)
    plt.title(f"Diagrama de Persistencia (r={radio}) · {nombre_csv}")
    plt.xlabel("Birth"); plt.ylabel("Death"); plt.tight_layout()
    plt.savefig(os.path.join(ruta_out,
                             f"{nombre_base}_diagrama_persistencia.png"))
    plt.close()

    # --- CSV con pares birth-death ------------------------------------------
    diagram_df = pd.DataFrame(
        [[dim, b, d] for dim, (b, d) in diag if dim <= 2],
        columns=["dimension", "birth", "death"]
    )
    diagram_df.to_csv(os.path.join(ruta_out, f"{nombre_base}.csv"),
                      index=False)
    return nombre_csv  # para saber cuál terminó


# --------------------------------------------------------------------------- #
#  FUNCIÓN PRINCIPAL
# --------------------------------------------------------------------------- #
def calcular_rips_y_persistencia(ruta_centroides: str,
                                 radio: float,
                                 n_workers: int) -> str:
    """Prepara carpetas, lanza procesos y muestra progreso."""
    ruta_resultados = os.path.join(ruta_centroides, "resultados")
    ruta_rips = os.path.join(ruta_resultados, f"rips_{radio}")
    os.makedirs(ruta_rips, exist_ok=True)

    # Archivos a procesar
    archivos_csv = sorted(
        f for f in os.listdir(ruta_centroides)
        if f.lower().endswith(".csv")
    )
    if not archivos_csv:
        print("  No se encontraron CSV en la ruta indicada.")
        return ruta_rips

    # --------------------------------------------------------------------- #
    #  Paralelismo con ProcessPoolExecutor
    # --------------------------------------------------------------------- #
    inicio = time.time()
    tarea = partial(_procesar_csv,
                    ruta_in=ruta_centroides,
                    ruta_out=ruta_rips,
                    radio=radio)

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(tarea, csv): csv for csv in archivos_csv}
        for _ in tqdm(as_completed(futures), total=len(futures),
                      desc=f"Procesando ({n_workers} núcleos)"):
            pass  # tqdm se va actualizando

    print(f"\n Resultados guardados en: {ruta_rips}")
    print(f"  Tiempo total: {time.time() - inicio:.2f} s")
    return ruta_rips


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def parse_args():
    parser = argparse.ArgumentParser(description="Generar complejos de Rips y diagramas de persistencia de archivos CSV")
    parser.add_argument("ruta_centroides", type=str, help="Ruta a la carpeta con archivos CSV")
    parser.add_argument("--radio", type=int, default=1000, help="Valor máximo de radio para el complejo de Rips (default=1000)")
    parser.add_argument("--workers", type=int, default=2, help="Número de núcleos para procesamiento paralelo (default=2)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    if not os.path.isdir(args.ruta_centroides):
        print(f" La ruta '{args.ruta_centroides}' no existe o no es directorio.")
        sys.exit(1)

    calcular_rips_y_persistencia(args.ruta_centroides,
                                 radio=args.radio,
                                 n_workers=args.workers)
