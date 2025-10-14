#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genera clustermaps a partir de archivos CSV que contienen matrices
de distancias entre centroides celulares, usando agrupamiento jerárquico.

Ejecución
---------
$ python generar_clustermaps.py /ruta/a/matrices [--metodo ward]

Argumentos
----------
/ruta/a/matrices : Carpeta con archivos CSV. Cada archivo debe contener
                   una matriz cuadrada de distancias (float) con nombres
                   de archivo como índices y columnas.

--metodo         : Método de linkage para clustermap. Opciones comunes:
                   - 'ward'     (default, requiere matriz euclidiana)
                   - 'single'   (mínima distancia)
                   - 'complete' (máxima distancia)
                   - 'average'  (promedio entre grupos)
                   - 'centroid' (centroide)
                   - 'median'   (mediana entre grupos)

Formato esperado
----------------
Cada archivo .csv debe contener:
- Índices y columnas: nombres de muestras (coinciden con archivos originales).
- Contenido: valores de distancia numéricos (float).
- El nombre del archivo puede contener pistas del tipo de muestra:
    'carcinoma', 'dysplasia', 'stroma_ad_dysplasia', 'stroma_ad_carcinoma'.

Colores asignados a tipos de muestra:
    carcinoma             → verde oscuro  (#6B990F)
    dysplasia             → azul oscuro   (#260F99)
    stroma-ad-carcinoma   → verde claro   (#E5FFB2)
    stroma-ad-dysplasia   → lila claro    (#BFB2FF)

Salidas
-------
Para cada archivo <nombre>.csv:
    visualizacion/clustermap_<nombre>.png

Todos los resultados se guardan en:
    <ruta>/visualizacion/
"""

import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------- #
# FUNCIONES AUXILIARES
# ---------------------------------------------------------------------------- #

def clean_filename(filename):
    """Quita prefijo y extensión del nombre de archivo."""
    filename = os.path.splitext(filename)[0]
    filename = filename.replace('filtrado_', '')
    return filename

def get_sample_type(filename):
    """Detecta tipo de muestra a partir del nombre del archivo."""
    filename = filename.lower()
    if 'stroma_ad_dysplasia' in filename:
        return 'stroma-ad-dysplasia'
    elif 'stroma_ad_carcinoma' in filename:
        return 'stroma-ad-carcinoma'
    elif 'dysplasia' in filename:
        return 'dysplasia'
    elif 'carcinoma' in filename:
        return 'carcinoma'
    else:
        return 'other'

# ---------------------------------------------------------------------------- #
# PLOTEO DEL CLUSTERMAP
# ---------------------------------------------------------------------------- #

def plot_clustermap(matrix, filenames, title, output_dir, metodo='ward'):
    """Genera y guarda el clustermap con anotaciones de tipo de muestra."""
    cleaned_filenames = [clean_filename(f) for f in filenames]
    sample_types = [get_sample_type(f) for f in cleaned_filenames]

    type_colors = {
        'dysplasia': '#260F99',
        'carcinoma': '#6B990F',
        'stroma-ad-dysplasia': '#BFB2FF',
        'stroma-ad-carcinoma': '#E5FFB2',
        'other': '#7f7f7f'
    }
    row_colors = [type_colors[tipo] for tipo in sample_types]

    g = sns.clustermap(matrix,
                       xticklabels=cleaned_filenames,
                       yticklabels=cleaned_filenames,
                       cmap='viridis',
                       figsize=(15, 15),
                       dendrogram_ratio=(.1, .2),
                       cbar_pos=(0.02, 0.8, 0.05, 0.18),
                       cbar_kws={'label': 'Distancia'},
                       row_colors=row_colors,
                       col_colors=row_colors,
                       method=metodo)

    # Etiquetas coloreadas y formateadas
    new_order = g.dendrogram_row.reordered_ind
    reordered_colors = [row_colors[i] for i in new_order]

    for label, color in zip(g.ax_heatmap.get_xticklabels(), reordered_colors):
        label.set_backgroundcolor(color)
        label.set_color('black')
        label.set_fontsize(6)
        label.set_rotation(90)
        label.set_ha('center')

    for label, color in zip(g.ax_heatmap.get_yticklabels(), reordered_colors):
        label.set_backgroundcolor(color)
        label.set_color('black')
        label.set_fontsize(6)
        label.set_rotation(0)
        label.set_va('center')

    # Título y leyenda
    g.fig.suptitle(f"Clustermap de {title}", fontsize=16, y=1.05)
    g.fig.subplots_adjust(right=0.8)
    legend_elements = [Patch(facecolor=color, label=label)
                       for label, color in type_colors.items()
                       if label in sample_types]

    g.ax_heatmap.legend(handles=legend_elements,
                        loc='upper left',
                        bbox_to_anchor=(1.01, 1.0),
                        title="Tipo de muestra",
                        fontsize=8,
                        title_fontsize=10)

    # Guardar
    plt.savefig(os.path.join(output_dir, f'clustermap_{title}.png'),
                bbox_inches='tight')
    plt.close()


# ---------------------------------------------------------------------------- #
# FUNCIÓN PRINCIPAL
# ---------------------------------------------------------------------------- #

def crear_visualizaciones(ruta_directorio, metodo='ward'):
    """Crea clustermaps para todos los CSV en una carpeta dada."""
    carpeta_visualizacion = os.path.join(ruta_directorio, "visualizacion")
    os.makedirs(carpeta_visualizacion, exist_ok=True)

    archivos_csv = [f for f in os.listdir(ruta_directorio)
                    if f.lower().endswith('.csv')]

    if not archivos_csv:
        print(" No se encontraron archivos CSV en la ruta indicada.")
        return

    for archivo in archivos_csv:
        nombre_base = clean_filename(os.path.splitext(archivo)[0])
        distancias = pd.read_csv(os.path.join(ruta_directorio, archivo), index_col=0)
        filenames = list(distancias.columns)
        plot_clustermap(distancias.values, filenames, nombre_base,
                        carpeta_visualizacion, metodo=metodo)

    print(f"Visualizaciones guardadas en: {carpeta_visualizacion}")


# ---------------------------------------------------------------------------- #
# CLI
# ---------------------------------------------------------------------------- #

def parse_args():
    parser = argparse.ArgumentParser(description="Genera clustermaps jerárquicos a partir de matrices de distancia.")
    parser.add_argument("ruta", type=str, help="Ruta a la carpeta con matrices CSV.")
    parser.add_argument("--metodo", type=str, default="ward",
                        help="Método de linkage: ward (default), single, complete, average, centroid, median")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not os.path.isdir(args.ruta):
        print(f"La ruta '{args.ruta}' no existe o no es un directorio válido.")
    else:
        crear_visualizaciones(args.ruta, metodo=args.metodo)
