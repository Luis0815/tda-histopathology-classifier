#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genera clustermaps jerárquicos a partir de matrices de distancia (CSV), con agrupación jerárquica
y anotaciones por tipo de muestra, estado Fanconi y origen anatómico.

Ejecución
---------
$ python clustermap_multiple.py /ruta/a/matrices [--metodo ward]

Argumentos
----------
/ruta/a/matrices : Ruta a la carpeta con archivos CSV. Cada archivo debe contener una matriz
                   cuadrada de distancias con nombres de muestra como índices y columnas.

--metodo         : Método de linkage para la agrupación jerárquica. Opciones válidas:
                   - 'ward'     (por defecto, requiere distancias euclidianas)
                   - 'single'   (mínima distancia entre grupos)
                   - 'complete' (máxima distancia entre grupos)
                   - 'average'  (promedio de distancias)
                   - 'centroid' (distancia entre centroides)
                   - 'median'   (mediana de distancias)

Salida
------
Para cada matriz, se generan tres clustermaps en subcarpetas:
    - visualizacion/por_tipo/
    - visualizacion/por_fanconi/
    - visualizacion/por_origen/
"""

import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# === Funciones de limpieza y clasificación ===

def clean_filename(filename):
    filename = os.path.splitext(filename)[0]
    filename = filename.replace('filtrado_', '')
    return filename

def get_sample_type(filename):
    filename = filename.lower()
    if 'stroma_ad_dysplasia' in filename or 'dysplasia-stroma' in filename:
        return 'stroma-ad-dysplasia'
    elif 'stroma_ad_carcinoma' in filename:
        return 'stroma-ad-carcinoma'
    elif 'dysplasia' in filename and 'lg' not in filename and 'hg' not in filename:
        return 'dysplasia'
    elif 'dysplasia-hg' in filename or 'dysplasia-lg' in filename:
        return 'dysplasia'
    elif 'carcinoma' in filename:
        return 'carcinoma'
    else:
        return 'other'

def get_fanconi_status(filename):
    return 'Fanconi' if 'FA' in filename else 'No Fanconi'

def get_origin(filename):
    if 'HNSCC' in filename:
        return 'Head and Neck'
    elif 'AGSCC' in filename:
        return 'Anogenital'
    else:
        return 'Otro'

# === Función para graficar clustermaps ===

def plot_clustermap(matrix, filenames, title, output_dir, labels, label_name, color_dict, metodo='ward'):
    cleaned_filenames = [clean_filename(f) for f in filenames]
    row_colors = [color_dict[labels[f]] for f in cleaned_filenames]

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

    new_row_order = g.dendrogram_row.reordered_ind
    reordered_colors = [row_colors[i] for i in new_row_order]

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

    g.fig.suptitle(f"{label_name}: {title}", fontsize=16, y=1.05)
    g.fig.subplots_adjust(right=0.8)

    legend_elements = [Patch(facecolor=color, label=label) 
                       for label, color in color_dict.items() 
                       if label in labels.values()]
    g.ax_heatmap.legend(
        handles=legend_elements,
        loc='upper left',
        bbox_to_anchor=(1.01, 1.0),
        title=label_name,
        fontsize=8,
        title_fontsize=10
    )

    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f'{label_name.lower().replace(" ", "_")}_clustermap_{title}.png'),
                bbox_inches='tight')
    plt.close()

# === Función principal ===

def crear_visualizaciones(ruta_directorio, metodo='ward'):
    carpeta_visualizacion = os.path.join(ruta_directorio, "visualizacion")
    os.makedirs(carpeta_visualizacion, exist_ok=True)

    archivos_csv = [f for f in os.listdir(ruta_directorio) if f.endswith('.csv')]

    for archivo in archivos_csv:
        nombre_base = clean_filename(os.path.splitext(archivo)[0])
        distancias = pd.read_csv(os.path.join(ruta_directorio, archivo), index_col=0)
        filenames = list(distancias.columns)
        cleaned_filenames = [clean_filename(f) for f in filenames]

        sample_types = {f: get_sample_type(f) for f in cleaned_filenames}
        fanconi_status = {f: get_fanconi_status(f) for f in cleaned_filenames}
        origins = {f: get_origin(f) for f in cleaned_filenames}

        type_colors = {
            'dysplasia': '#260F99',
            'carcinoma': '#6B990F',
            'stroma-ad-dysplasia': '#BFB2FF',
            'stroma-ad-carcinoma': '#E5FFB2',
            'other': '#7f7f7f'
        }

        fanconi_colors = {
            'Fanconi': '#d73027',
            'No Fanconi': '#4575b4'
        }

        origin_colors = {
            'Head and Neck': '#66c2a5',
            'Anogenital': '#fc8d62',
            'Otro': '#bdbdbd'
        }

        plot_clustermap(distancias.values, filenames, nombre_base,
                        os.path.join(carpeta_visualizacion, "por_tipo"),
                        sample_types, "Tipo de muestra", type_colors, metodo=metodo)

        plot_clustermap(distancias.values, filenames, nombre_base,
                        os.path.join(carpeta_visualizacion, "por_fanconi"),
                        fanconi_status, "Estado Fanconi", fanconi_colors, metodo=metodo)

        plot_clustermap(distancias.values, filenames, nombre_base,
                        os.path.join(carpeta_visualizacion, "por_origen"),
                        origins, "Origen anatómico", origin_colors, metodo=metodo)

    print("✔ Visualizaciones generadas en subcarpetas: por_tipo, por_fanconi, por_origen.")

# === CLI ===

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
