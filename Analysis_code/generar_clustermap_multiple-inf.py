#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genera clustermaps jerárquicos a partir de matrices de distancia (CSV), con agrupación jerárquica
y anotaciones por tipo de muestra y estado Fanconi.

Ejecución
---------
$ python clustermap_multiple.py /ruta/a/matrices [--metodo average]
"""

import os
import argparse
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# === Funciones de limpieza y clasificación ===

def clean_filename(filename):
    """Limpia el nombre del archivo (quita 'filtrado_' y la extensión)."""
    filename = os.path.splitext(filename)[0]
    filename = filename.replace('filtrado_', '')
    return filename


def get_sample_type(filename):
    """Clasifica el tipo de muestra con base en el nombre."""
    filename = filename.lower()

    if '_and_stroma' in filename or '-and-stroma' in filename:
        return 'and-stroma'
    elif 'stroma_ad' in filename and 'dysplasia' in filename:
        return 'stroma-ad-dysplasia'
    elif 'stroma_ad' in filename and 'carcinoma' in filename:
        return 'stroma-ad-carcinoma'
    elif 'dysplasia' in filename:
        return 'dysplasia'
    elif 'carcinoma' in filename:
        return 'carcinoma'
    else:
        return 'other'


def get_fanconi_status(filename):
    """
    Determina el estado Fanconi: basta con que el nombre contenga una 'F' (mayúscula).
    Ejemplo:
        Carcinoma_FAHNSCC_14_1.csv → Fanconi
        stroma_ad_HG_dysplasia_F79P1_1.csv → Fanconi
        carcinoma_invasive_14_1.csv → No Fanconi
    """
    name = os.path.splitext(os.path.basename(filename))[0]
    return 'Fanconi' if 'F' in name else 'No Fanconi'


# === (Comentado) Origen anatómico ===
# def get_origin(filename):
#     """Clasifica el origen anatómico."""
#     if 'HNSCC' in filename:
#         return 'Head and Neck'
#     elif 'AGSCC' in filename:
#         return 'Anogenital'
#     else:
#         return 'Otro'


# === Función para graficar clustermap combinado ===

def plot_clustermap_combinado(matrix_df, title, output_dir,
                              sample_types, fanconi_status,
                              type_colors, fanconi_colors,
                              # origin_colors, origins,   # ← comentado
                              metodo='ward'):
    """Genera el clustermap con anotaciones de tipo y Fanconi."""

    muestras = matrix_df.index.tolist()

    # DataFrame de colores laterales
    row_colors_df = pd.DataFrame({
        "Tipo": [type_colors.get(sample_types[m], '#999999') for m in muestras],
        "Fanconi": [fanconi_colors.get(fanconi_status[m], '#999999') for m in muestras],
        # "Origen": [origin_colors.get(origins[m], '#999999') for m in muestras]  # ← comentado
    }, index=muestras)

    # Crear clustermap
    g = sns.clustermap(
        matrix_df,
        cmap='viridis',
        figsize=(15, 18),
        row_colors=row_colors_df,
        col_colors=row_colors_df,
        method=metodo,
        cbar_kws={'label': 'Distancia'},
        dendrogram_ratio=(.1, .2),
        cbar_pos=(0.02, 0.8, 0.05, 0.18),
        xticklabels=True,
        yticklabels=False
    )

    plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=6, rotation=90)
    g.fig.suptitle(f"Clustermap combinado: {title}", fontsize=16, y=1.05)

    # Leyendas
    legend_tipo = [Patch(color=c, label=l) for l, c in type_colors.items() if l in sample_types.values()]
    legend_fanconi = [Patch(color=c, label=l) for l, c in fanconi_colors.items() if l in fanconi_status.values()]
    # legend_origen = [Patch(color=c, label=l) for l, c in origin_colors.items() if l in origins.values()]  # ← comentado

    legend_titles = [
        Patch(facecolor='white', edgecolor='white', label='Tipo de muestra'),
        *legend_tipo,
        Patch(facecolor='white', edgecolor='white', label='Estado Fanconi'),
        *legend_fanconi,
        # Patch(facecolor='white', edgecolor='white', label='Origen anatómico'),
        # *legend_origen,  # ← comentado
    ]

    g.fig.legend(
        handles=legend_titles,
        loc='center left',
        bbox_to_anchor=(0.9, 0.5),
        fontsize=9,
        title_fontsize=10,
        borderaxespad=2.0,
        frameon=False
    )

    g.fig.subplots_adjust(right=0.8)

    os.makedirs(output_dir, exist_ok=True)
    g.savefig(os.path.join(output_dir, f'clustermap_combinado_{title}.svg'), bbox_inches='tight')
    plt.close()


# === Función principal ===

def crear_visualizaciones(ruta_directorio, metodo='ward'):
    carpeta_visualizacion = os.path.join(ruta_directorio, "visualizacion", "combinado")
    os.makedirs(carpeta_visualizacion, exist_ok=True)

    archivos_csv = [f for f in os.listdir(ruta_directorio) if f.endswith('.csv')]

    for archivo in archivos_csv:
        path = os.path.join(ruta_directorio, archivo)
        nombre_base = clean_filename(os.path.splitext(archivo)[0])

        distancias = pd.read_csv(path, index_col=0)

        if distancias.shape[0] != distancias.shape[1]:
            print(f"⚠️  Saltando '{archivo}': matriz no cuadrada ({distancias.shape}).")
            continue
        if not (distancias.index.tolist() == distancias.columns.tolist()):
            print(f"⚠️  Saltando '{archivo}': nombres de filas y columnas no coinciden.")
            continue

        cleaned_names = [clean_filename(name) for name in distancias.index]
        distancias.index = cleaned_names
        distancias.columns = cleaned_names

        sample_types = {f: get_sample_type(f) for f in cleaned_names}
        fanconi_status = {f: get_fanconi_status(f) for f in cleaned_names}
        # origins = {f: get_origin(f) for f in cleaned_names}  # ← comentado

        # Colores
        type_colors = {
            'dysplasia': '#260F99',
            'carcinoma': '#6B990F',
            'stroma-ad-dysplasia': '#BFB2FF',
            'stroma-ad-carcinoma': '#E5FFB2',
            'and-stroma': '#000000',
            'other': '#7f7f7f'
        }

        fanconi_colors = {
            'Fanconi': '#d73027',
            'No Fanconi': '#4575b4'
        }

        # origin_colors = {
        #     'Head and Neck': '#66c2a5',
        #     'Anogenital': '#fc8d62',
        #     'Otro': '#bdbdbd'
        # }

        plot_clustermap_combinado(distancias, nombre_base, carpeta_visualizacion,
                                  sample_types, fanconi_status,
                                  type_colors, fanconi_colors,
                                  # origin_colors, origins,  # ← comentado
                                  metodo=metodo)

    print("✔ Visualizaciones generadas en:", carpeta_visualizacion)


# === CLI ===

def parse_args():
    parser = argparse.ArgumentParser(description="Genera un clustermap con anotaciones de tipo y Fanconi.")
    parser.add_argument("ruta", type=str, help="Ruta a la carpeta con matrices CSV.")
    parser.add_argument("--metodo", type=str, default="average",
                        help="Método de linkage: ward (default), single, complete, average, centroid, median")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not os.path.isdir(args.ruta):
        print(f"La ruta '{args.ruta}' no existe o no es un directorio válido.")
    else:
        crear_visualizaciones(args.ruta, metodo=args.metodo)
