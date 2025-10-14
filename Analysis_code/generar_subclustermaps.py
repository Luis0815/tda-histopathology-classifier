#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genera clustermaps jerárquicos con anotaciones por tipo de muestra,
estado Fanconi y origen anatómico, incluyendo:

1. Un mapa combinado completo.
2. Subgrupos: 
   - Carcinoma FA vs No FA
   - Stroma-ad FA vs No FA (tanto en carcinoma como displasia)
   - Displasia FA vs No FA
   - Carcinoma + Displasia FA vs No FA

Salida:
- visualizacion/combinado/
"""

import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# === Clasificadores ===
def clean_filename(filename):
    filename = os.path.splitext(filename)[0].replace('filtrado_', '')
    return filename

def get_sample_type(filename):
    fname = filename.lower()
    if 'stroma_ad_dysplasia' in fname or 'dysplasia-stroma' in fname:
        return 'stroma-ad-dysplasia'
    elif 'stroma_ad_carcinoma' in fname:
        return 'stroma-ad-carcinoma'
    elif 'dysplasia' in fname:
        return 'dysplasia'
    elif 'carcinoma' in fname:
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

# === Gráfico de clustermap ===
def plot_clustermap_combinado(matrix_df, title, output_dir,
                              sample_types, fanconi_status, origins,
                              type_colors, fanconi_colors, origin_colors,
                              metodo='ward'):

    muestras = matrix_df.index.tolist()
    row_colors_df = pd.DataFrame({
        "Tipo": [type_colors[sample_types[m]] for m in muestras],
        "Fanconi": [fanconi_colors[fanconi_status[m]] for m in muestras],
        "Origen": [origin_colors[origins[m]] for m in muestras]
    }, index=muestras)

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

    legend_tipo = [Patch(color=c, label=l) for l, c in type_colors.items() if l in sample_types.values()]
    legend_fanconi = [Patch(color=c, label=l) for l, c in fanconi_colors.items() if l in fanconi_status.values()]
    legend_origen = [Patch(color=c, label=l) for l, c in origin_colors.items() if l in origins.values()]

    g.fig.legend(
        handles=[
            Patch(facecolor='white', edgecolor='white', label='Tipo de muestra'), *legend_tipo,
            Patch(facecolor='white', edgecolor='white', label='Estado Fanconi'), *legend_fanconi,
            Patch(facecolor='white', edgecolor='white', label='Origen anatómico'), *legend_origen
        ],
        loc='center left', bbox_to_anchor=(0.9, 0.5), fontsize=9, frameon=False
    )

    g.fig.subplots_adjust(right=0.8)
    os.makedirs(output_dir, exist_ok=True)
    g.savefig(os.path.join(output_dir, f'clustermap_combinado_{title}.png'), bbox_inches='tight')
    plt.close()

# === Visualización general y subgrupos ===
def crear_visualizaciones(ruta_directorio, metodo='ward'):
    carpeta_salida = os.path.join(ruta_directorio, "visualizacion", "combinado")
    os.makedirs(carpeta_salida, exist_ok=True)

    archivos_csv = [f for f in os.listdir(ruta_directorio) if f.endswith('.csv')]

    for archivo in archivos_csv:
        path = os.path.join(ruta_directorio, archivo)
        nombre_base = clean_filename(os.path.splitext(archivo)[0])

        distancias = pd.read_csv(path, index_col=0)

        if distancias.shape[0] != distancias.shape[1]:
            print(f" Saltando '{archivo}': matriz no cuadrada ({distancias.shape}).")
            continue

        if not (distancias.index.tolist() == distancias.columns.tolist()):
            print(f" Saltando '{archivo}': nombres de filas y columnas no coinciden.")
            continue

        cleaned_names = [clean_filename(name) for name in distancias.index]
        distancias.index = cleaned_names
        distancias.columns = cleaned_names

        tipos = {f: get_sample_type(f) for f in cleaned_names}
        fanconi = {f: get_fanconi_status(f) for f in cleaned_names}
        origen = {f: get_origin(f) for f in cleaned_names}

        colores_tipo = {
            'dysplasia': '#260F99',
            'carcinoma': '#6B990F',
            'stroma-ad-dysplasia': '#BFB2FF',
            'stroma-ad-carcinoma': '#E5FFB2',
            'other': '#7f7f7f'
        }
        colores_fanconi = {'Fanconi': '#d73027', 'No Fanconi': '#4575b4'}
        colores_origen = {'Head and Neck': '#66c2a5', 'Anogenital': '#fc8d62', 'Otro': '#bdbdbd'}

        subconjuntos = {
            "completo": cleaned_names,
            "carcinoma": [s for s in cleaned_names if tipos[s] == 'carcinoma'],
            "dysplasia": [s for s in cleaned_names if tipos[s] == 'dysplasia'],
            "stroma-ad": [s for s in cleaned_names if tipos[s].startswith('stroma-ad')],
            "carcinoma_dysplasia": [s for s in cleaned_names if tipos[s] in ('carcinoma', 'dysplasia')]
        }

        for nombre_sub, muestras in subconjuntos.items():
            if len(muestras) < 3:
                continue
            submatriz = distancias.loc[muestras, muestras]
            plot_clustermap_combinado(submatriz, f"{nombre_base}_{nombre_sub}", carpeta_salida,
                                      {k: tipos[k] for k in muestras},
                                      {k: fanconi[k] for k in muestras},
                                      {k: origen[k] for k in muestras},
                                      colores_tipo, colores_fanconi, colores_origen,
                                      metodo=metodo)

    print("✔ Visualizaciones generadas en:", carpeta_salida)

# === CLI ===
def parse_args():
    parser = argparse.ArgumentParser(description="Genera clustermaps por subgrupo y anotaciones laterales.")
    parser.add_argument("ruta", type=str, help="Ruta a la carpeta con matrices CSV.")
    parser.add_argument("--metodo", type=str, default="ward",
                        help="Método de linkage: ward, single, complete, average, centroid, median")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if not os.path.isdir(args.ruta):
        print(f"La ruta '{args.ruta}' no es un directorio válido.")
    else:
        crear_visualizaciones(args.ruta, metodo=args.metodo)
