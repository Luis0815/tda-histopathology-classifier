#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calcula distancias de Wasserstein (y opcionalmente Bottleneck) entre diagramas
de persistencia almacenados como CSV (col. 'dimension', 'birth', 'death').

Uso
----
$ python calcular_distancias.py /ruta/a/rips_1000 [--workers 4]

- /ruta/a/rips_1000  : carpeta con muchos <nombre>.csv (diagramas)
- --workers N        : núcleos a usar (default=4)
"""

import os
import sys
import time
import argparse
from itertools import combinations_with_replacement
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import gudhi as gd
import gudhi.wasserstein as gw
from tqdm import tqdm


# --------------------------------------------------------------------------- #
#  FUNCIÓN QUE CARGA UN CSV Y DEVUELVE DOS ARRAYS (dim 0 y dim 1)
# --------------------------------------------------------------------------- #
def _cargar_csv(ruta_csv):
    df = pd.read_csv(ruta_csv)
    if {'birth', 'death', 'dimension'}.issubset(df.columns):
        diag0 = df[df["dimension"] == 0][["birth", "death"]].to_numpy(dtype=float)
        diag1 = df[df["dimension"] == 1][["birth", "death"]].to_numpy(dtype=float)
        return diag0, diag1
    else:
        raise ValueError(f"{os.path.basename(ruta_csv)} no tiene columnas "
                         "'dimension', 'birth', 'death'")


# --------------------------------------------------------------------------- #
#  FUNCIÓN QUE CALCULA DISTANCIAS ENTRE UN PAR 
# --------------------------------------------------------------------------- #
def _distancia_par(args):
    # Desempaquetar
    i, j, arch_i, arch_j, diag_i_0, diag_i_1, diag_j_0, diag_j_1 = args

    # --- Wasserstein --------------------------------------------------------
    wass0 = gw.wasserstein_distance(diag_i_0, diag_j_0, order=1)
    wass1 = gw.wasserstein_distance(diag_i_1, diag_j_1, order=1)

    # --- Bottleneck -----------------
    # bott0 = gd.bottleneck_distance(diag_i_0, diag_j_0)
    # bott1 = gd.bottleneck_distance(diag_i_1, diag_j_1)

    return i, j, wass0, wass1  # , bott0, bott1


# --------------------------------------------------------------------------- #
#  FUNCIÓN PRINCIPAL
# --------------------------------------------------------------------------- #
def calcular_distancias(ruta_directorio: str, workers: int = 2) -> str:
    t0 = time.time()

    # Carpeta de salida
    carpeta_salida = os.path.join(os.path.dirname(ruta_directorio),
                                  "distancias_wasserstein")
    os.makedirs(carpeta_salida, exist_ok=True)

    # ------------------- Cargar todos los diagramas -------------------------
    archivos = sorted(f for f in os.listdir(ruta_directorio)
                      if f.lower().endswith(".csv"))

    if not archivos:
        print(" No se encontraron archivos CSV en la ruta indicada.")
        return carpeta_salida

    diag0_list, diag1_list = [], []
    for nombre in archivos:
        try:
            d0, d1 = _cargar_csv(os.path.join(ruta_directorio, nombre))
            diag0_list.append(d0)
            diag1_list.append(d1)
        except Exception as e:
            print(f" {nombre}: {e}")
            return carpeta_salida

    n = len(archivos)
    dist_wass0 = np.zeros((n, n))
    dist_wass1 = np.zeros((n, n))
    # dist_bott0 = np.zeros((n, n))   
    # dist_bott1 = np.zeros((n, n))

    # ------------------- Preparar tareas para paralelo ----------------------
    tareas = []
    for (i, j) in combinations_with_replacement(range(n), 2):
        tareas.append((i, j,
                       archivos[i], archivos[j],
                       diag0_list[i], diag1_list[i],
                       diag0_list[j], diag1_list[j]))

    # ------------------- Ejecutar en paralelo -------------------------------
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_distancia_par, t) for t in tareas]

        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc=f"Calculando distancias ({workers} núcleos)"):
            i, j, w0, w1 = fut.result()               # , b0, b1
            dist_wass0[i, j] = dist_wass0[j, i] = w0
            dist_wass1[i, j] = dist_wass1[j, i] = w1
            # dist_bott0[i, j] = dist_bott0[j, i] = b0
            # dist_bott1[i, j] = dist_bott1[j, i] = b1

    # ------------------- Guardar resultados ---------------------------------
    pd.DataFrame(dist_wass0, index=archivos,
                 columns=archivos).to_csv(
        os.path.join(carpeta_salida, "wasserstein_dim0.csv"))
    pd.DataFrame(dist_wass1, index=archivos,
                 columns=archivos).to_csv(
        os.path.join(carpeta_salida, "wasserstein_dim1.csv"))

    # Si activas Bottleneck, descomenta estas líneas
    # pd.DataFrame(dist_bott0, index=archivos,
    #              columns=archivos).to_csv(
    #     os.path.join(carpeta_salida, "bottleneck_dim0.csv"))
    # pd.DataFrame(dist_bott1, index=archivos,
    #              columns=archivos).to_csv(
    #     os.path.join(carpeta_salida, "bottleneck_dim1.csv"))

    print(f"\n Distancias guardadas en: {carpeta_salida}")
    print(f"  Tiempo total: {time.time() - t0:.2f} s")
    return carpeta_salida


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(
        description="Calcular distancias de Wasserstein entre diagramas CSV.")
    p.add_argument("ruta", help="Carpeta con diagramas (CSV).")
    p.add_argument("--workers", type=int, default=4,
                   help="Núcleos a usar (default=4).")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not os.path.isdir(args.ruta):
        print(f" La ruta '{args.ruta}' no existe o no es un directorio.")
        sys.exit(1)

    calcular_distancias(args.ruta, workers=args.workers)
