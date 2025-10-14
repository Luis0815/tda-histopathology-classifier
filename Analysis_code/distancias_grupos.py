#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calcula distancias de Wasserstein (y opcionalmente Bottleneck) **por grupo
celular** a partir de diagramas de persistencia almacenados en CSV
(col. 'dimension', 'birth', 'death').

Ejecución
---------
$ python distancias_grupos.py /ruta/a/diagramas   \
                                 [--workers 4]       \
                                 [--bottleneck]

Argumentos
----------
/ruta/a/diagramas : carpeta con muchos <archivo>.csv generados previamente
                    por rips_grupos.py (nombre termina en _<grupo>.csv).
--workers N       : núcleos que se usarán (int, default = 4).
--bottleneck      : si se indica, también se calculan distancias Bottleneck.

Salidas
-------
Para cada grupo ('tumorales', 'linfoides', 'mieloides', 'no_tumorales'):

    distancias_wasserstein_dim0_<grupo>.csv
    distancias_wasserstein_dim1_<grupo>.csv
    (y opcionalmente)
    distancias_bottleneck_dim0_<grupo>.csv
    distancias_bottleneck_dim1_<grupo>.csv

Todos los ficheros se guardan en:
    <ruta_diagramas>/resultados/distancias_grupos/

Compatibilidad
--------------
Probado en Python ≥3.7.  Para Python <3.9 se usan anotaciones de typing del
módulo estándar (List, Tuple).
"""

import os
import re
import sys
import time
import argparse
from itertools import combinations_with_replacement
from typing import List, Tuple, Dict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import gudhi as gd
import gudhi.wasserstein as gw
from tqdm import tqdm

# --------------------------------------------------------------------------- #
#  CONSTANTES
# --------------------------------------------------------------------------- #
GRUPOS_PERMITIDOS: List[str] = [
    'tumorales', 'linfoides', 'mieloides', 'no_tumorales'
]

PATRON_GRUPO = re.compile(r'_(' + '|'.join(GRUPOS_PERMITIDOS) + r')\.csv$', re.I)


# --------------------------------------------------------------------------- #
#  UTILIDADES
# --------------------------------------------------------------------------- #
def extraer_grupo(nombre_archivo: str) -> str:
    """Devuelve el nombre del grupo según el sufijo del archivo o 'desconocido'."""
    m = PATRON_GRUPO.search(nombre_archivo)
    return m.group(1).lower() if m else "desconocido"


def cargar_diagramas(ruta_csv: str) -> Tuple[np.ndarray, np.ndarray]:
    """Lee un CSV y devuelve (diag_dim0, diag_dim1)."""
    df = pd.read_csv(ruta_csv, usecols=['dimension', 'birth', 'death'])
    diag0 = df[df["dimension"] == 0][['birth', 'death']].to_numpy(dtype=float)
    diag1 = df[df["dimension"] == 1][['birth', 'death']].to_numpy(dtype=float)
    return diag0, diag1


def distancia_par(args):
    """Función auxiliar para cálculo paralelo."""
    i, j, diag_i0, diag_i1, diag_j0, diag_j1, calc_bott = args
    w0 = gw.wasserstein_distance(diag_i0, diag_j0, order=1)
    w1 = gw.wasserstein_distance(diag_i1, diag_j1, order=1)
    if calc_bott:
        b0 = gd.bottleneck_distance(diag_i0, diag_j0)
        b1 = gd.bottleneck_distance(diag_i1, diag_j1)
        return i, j, w0, w1, b0, b1
    return i, j, w0, w1, None, None


# --------------------------------------------------------------------------- #
#  DISTANCIAS POR GRUPO
# --------------------------------------------------------------------------- #
def procesar_grupo(diags0: Dict[str, np.ndarray],
                   diags1: Dict[str, np.ndarray],
                   carpeta_out: str,
                   workers: int,
                   calc_bottleneck: bool):
    archivos = list(diags0.keys())
    n = len(archivos)

    # Inicializar matrices
    m_w0 = np.zeros((n, n))
    m_w1 = np.zeros((n, n))
    m_b0 = np.zeros((n, n))
    m_b1 = np.zeros((n, n))

    # Crear lista de tareas
    tareas = [
        (i, j,
         diags0[archivos[i]], diags1[archivos[i]],
         diags0[archivos[j]], diags1[archivos[j]],
         calc_bottleneck)
        for i, j in combinations_with_replacement(range(n), 2)
    ]

    with ProcessPoolExecutor(max_workers=workers) as pool:
        for res in tqdm(pool.map(distancia_par, tareas),
                        total=len(tareas), desc="  pares"):
            i, j, w0, w1, b0, b1 = res
            m_w0[i, j] = m_w0[j, i] = w0
            m_w1[i, j] = m_w1[j, i] = w1
            if calc_bottleneck:
                m_b0[i, j] = m_b0[j, i] = b0
                m_b1[i, j] = m_b1[j, i] = b1

    # Guardar matrices
    pd.DataFrame(m_w0, index=archivos, columns=archivos).to_csv(
        os.path.join(carpeta_out, "distancias_wasserstein_dim0.csv"))
    pd.DataFrame(m_w1, index=archivos, columns=archivos).to_csv(
        os.path.join(carpeta_out, "distancias_wasserstein_dim1.csv"))

    if calc_bottleneck:
        pd.DataFrame(m_b0, index=archivos, columns=archivos).to_csv(
            os.path.join(carpeta_out, "distancias_bottleneck_dim0.csv"))
        pd.DataFrame(m_b1, index=archivos, columns=archivos).to_csv(
            os.path.join(carpeta_out, "distancias_bottleneck_dim1.csv"))


# --------------------------------------------------------------------------- #
#  FUNCIÓN PRINCIPAL
# --------------------------------------------------------------------------- #
def distancias_por_grupo(ruta_dir: str,
                         workers: int = 2,
                         calc_bottleneck: bool = False):
    t0 = time.time()

    ruta_dir     = os.path.abspath(ruta_dir)              
    parent_dir   = os.path.dirname(ruta_dir)               
    carpeta_out_base = os.path.join(parent_dir, "distancias_grupos")
    os.makedirs(carpeta_out_base, exist_ok=True)

    # Diccionarios para almacenar diagramas por grupo
    datos0 = {g: {} for g in GRUPOS_PERMITIDOS}
    datos1 = {g: {} for g in GRUPOS_PERMITIDOS}

    # Recorrer CSV
    for f in os.listdir(ruta_dir):
        if not f.lower().endswith(".csv"):
            continue
        grupo = extraer_grupo(f)
        if grupo not in GRUPOS_PERMITIDOS:
            continue
        diag0, diag1 = cargar_diagramas(os.path.join(ruta_dir, f))
        datos0[grupo][f] = diag0
        datos1[grupo][f] = diag1

    # Procesar cada grupo
    for grupo in GRUPOS_PERMITIDOS:
        if not datos0[grupo]:
            print(f"  Grupo {grupo} sin archivos, se omite.")
            continue
        print(f" Calculando distancias para grupo '{grupo}' "
              f"({len(datos0[grupo])} archivos)…")
        carpeta_grupo = os.path.join(carpeta_out_base, grupo)
        os.makedirs(carpeta_grupo, exist_ok=True)
        procesar_grupo(datos0[grupo], datos1[grupo],
                       carpeta_grupo, workers, calc_bottleneck)

    print(f"\n Resultados guardados en: {carpeta_out_base}")
    print(f"  Tiempo total: {time.time() - t0:.2f} s")


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def get_args():
    p = argparse.ArgumentParser(
        description="Calcula distancias Wasserstein (y opc. Bottleneck) "
                    "por grupo celular.")
    p.add_argument("ruta", help="Carpeta con diagramas CSV")
    p.add_argument("--workers", type=int, default=4,
                   help="Núcleos a usar (default=4)")
    p.add_argument("--bottleneck", action="store_true",
                   help="Incluir distancias Bottleneck")
    return p.parse_args()


# --------------------------------------------------------------------------- #
#  MAIN
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    args = get_args()

    if not os.path.isdir(args.ruta):
        print(f" La ruta '{args.ruta}' no existe o no es un directorio.")
        sys.exit(1)

    distancias_por_grupo(args.ruta,
                         workers=args.workers,
                         calc_bottleneck=args.bottleneck)
