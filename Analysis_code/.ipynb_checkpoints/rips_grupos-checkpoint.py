#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genera complejos de Rips y diagramas de persistencia para subconjuntos celulares
(tumorales, linfoides, mieloides y no tumorales) a partir de archivos CSV con centroides.

Ejecución
---------
$ python rips_grupos.py /ruta/a/csvs [--radio 1000] [--workers 2]

Argumentos
----------
ruta/a/csvs   : carpeta con archivos .csv con columnas 'X_centroid', 'Y_centroid', 'phenotype_key'
--radio       : radio máximo (max_edge_length) para el complejo de Rips (float, default 1000)
--workers     : núcleos para procesamiento paralelo (int, default = 2)

Salidas
-------
Por cada archivo y grupo:
    <nombre_archivo>_<grupo>_complejo_rips.png
    <nombre_archivo>_<grupo>_diagrama_persistencia.png
    <nombre_archivo>_<grupo>.csv  (tabla birth–death)

Todos los resultados se guardan en:
    <ruta_csvs>/resultados/rips_grupos_<radio>/
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import gudhi as gd
import matplotlib.pyplot as plt
from tqdm import tqdm
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed

# Grupos celulares
GRUPOS = {
    'tumorales': ['tumor cells', 'Ki67+ tumor cells'],
    'linfoides': ['NK', 'B cells', 'effector CD8+ T cells', 'memory CD8+ T cells',
                  'CD4+ T cells', 'regulatory T cells', 'memory CD4+ T cells'],
    'mieloides': ['neutrophils', 'other APCs', 'dendritic cells',
                  'M1/M0 macrophages', 'M2 macrophages'],
    'no_tumorales': ['endothelial cells', 'stromal cells']
}

# --------------------------------------------------------------------------- #
# FUNCIÓN PARA PROCESAR UN SOLO ARCHIVO
# --------------------------------------------------------------------------- #
def procesar_archivo(nombre_csv, ruta_in, ruta_out, radio):
    """Procesa un archivo CSV generando Rips y persistencia por grupo celular"""
    ruta_csv = os.path.join(ruta_in, nombre_csv)
    df = pd.read_csv(ruta_csv)

    for grupo, tipos in GRUPOS.items():
        df_grupo = df[df["phenotype"].isin(tipos)]
        if df_grupo.empty:
            continue

        puntos = df_grupo[["X_centroid", "Y_centroid"]].to_numpy()
        if len(puntos) < 2:
            continue

        rips_complex = gd.RipsComplex(points=puntos, max_edge_length=radio)
        simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)

        base = os.path.splitext(nombre_csv)[0]
        nombre_out = f"{base}_{grupo}"

        # Imagen del complejo
        plt.figure(figsize=(24, 10))
        plt.scatter(puntos[:, 0], puntos[:, 1], color="black", s=5, label="Centroides")
        for simplex in simplex_tree.get_skeleton(1):
            if len(simplex[0]) == 2:
                i, j = simplex[0]
                plt.plot([puntos[i][0], puntos[j][0]],
                         [puntos[i][1], puntos[j][1]],
                         color="gray", linewidth=0.5)
        plt.title(f"Rips ({grupo}) · {nombre_csv}")
        plt.xlabel("X"); plt.ylabel("Y"); plt.legend(); plt.tight_layout()
        plt.savefig(os.path.join(ruta_out, f"{nombre_out}_complejo_rips.png"))
        plt.close()

        # Diagrama de persistencia
        diag = simplex_tree.persistence()
        diagram_df = pd.DataFrame(
            [[dim, b, d] for dim, (b, d) in diag if dim <= 2],
            columns=["dimension", "birth", "death"]
        )
        diagram_df.to_csv(os.path.join(ruta_out, f"{nombre_out}.csv"),
                          index=False)

        plt.figure(figsize=(6, 6))
        gd.plot_persistence_diagram(diag)
        plt.title(f"Persistencia ({grupo}) · {nombre_csv}")
        plt.xlabel("Birth"); plt.ylabel("Death"); plt.tight_layout()
        plt.savefig(os.path.join(ruta_out, f"{nombre_out}_diagrama_persistencia.png"))
        plt.close()

    return nombre_csv

# --------------------------------------------------------------------------- #
# FUNCIÓN PRINCIPAL
# --------------------------------------------------------------------------- #
def calcular_rips_grupos(ruta_csvs: str, radio: float, n_workers: int):
    """Ejecuta el procesamiento paralelo de todos los CSV"""
    ruta_resultados = os.path.join(ruta_csvs, "resultados")
    ruta_out = os.path.join(ruta_resultados, f"rips_grupos_{radio}")
    os.makedirs(ruta_out, exist_ok=True)

    archivos = sorted(f for f in os.listdir(ruta_csvs) if f.endswith(".csv"))
    if not archivos:
        print("  No se encontraron archivos .csv en la ruta.")
        return ruta_out

    print(f"Procesando {len(archivos)} archivos con {n_workers} núcleos...")
    tarea = partial(procesar_archivo, ruta_in=ruta_csvs, ruta_out=ruta_out, radio=radio)

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(tarea, csv): csv for csv in archivos}
        for _ in tqdm(as_completed(futures), total=len(futures),
                      desc="Procesando archivos", unit="archivo"):
            pass

    print(f"\nResultados guardados en: {ruta_out}")
    return ruta_out

# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def parse_args():
    parser = argparse.ArgumentParser(description="Genera Rips y persistencia por grupo celular desde CSV")
    parser.add_argument("ruta_csvs", type=str, help="Ruta a carpeta con archivos CSV")
    parser.add_argument("--radio", type=int, default=1000, help="Radio máximo para Rips (default=1000)")
    parser.add_argument("--workers", type=int, default=2, help="Núcleos para procesamiento paralelo (default=2)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    if not os.path.isdir(args.ruta_csvs):
        print(f"La ruta '{args.ruta_csvs}' no existe o no es un directorio.")
        sys.exit(1)

    calcular_rips_grupos(args.ruta_csvs, radio=args.radio, n_workers=args.workers)
