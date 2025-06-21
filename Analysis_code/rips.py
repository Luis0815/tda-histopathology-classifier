#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm

def calcular_rips_y_persistencia(ruta_centroides):
    ruta_persistencia = os.path.join(ruta_centroides, "imagenes")
    if not os.path.exists(ruta_persistencia):
        os.makedirs(ruta_persistencia)

    archivos_csv = [archivo for archivo in os.listdir(ruta_centroides) if archivo.endswith('.csv')]

    for archivo_csv in tqdm(archivos_csv, desc="Procesando archivos"):
        ruta_completa = os.path.join(ruta_centroides, archivo_csv)

        df = pd.read_csv(ruta_completa)
        centroides_x = df['X_centroid'].tolist()
        centroides_y = df['Y_centroid'].tolist()
        
        puntos = np.array(list(zip(centroides_x, centroides_y)))

        rips_complex = gd.RipsComplex(points=puntos, max_edge_length=80)
        simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)

        plt.figure(figsize=(24, 10))
        plt.scatter(centroides_x, centroides_y, color='black', label='Centroides', s=5)
        for simplex in simplex_tree.get_skeleton(1):
            if len(simplex[0]) == 2:  
                arista = simplex[0]
                x = [centroides_x[i] for i in arista]
                y = [centroides_y[i] for i in arista]
                plt.plot(x, y, color='gray', linestyle='-', linewidth=0.5)
        plt.xlabel('Coordenada X')
        plt.ylabel('Coordenada Y')
        plt.title(f'Complejo de Rips ({archivo_csv})')
        plt.legend()
        
        nombre_imagen_rips = f"{os.path.splitext(archivo_csv)[0]}_complejo_rips.png"
        ruta_imagen_rips = os.path.join(ruta_persistencia, nombre_imagen_rips)
        plt.tight_layout()
        plt.savefig(ruta_imagen_rips)
        plt.close()

        diag = simplex_tree.persistence()

    
        persistencia_aplanada = []
        for d in diag:
            dimension, (birth, death) = d
            if dimension <= 2:
                persistencia_aplanada.append([dimension, birth, death])

        diagram_df = pd.DataFrame(persistencia_aplanada, columns=['dimension', 'birth', 'death'])

        nombre_diagrama_csv = f"{os.path.splitext(archivo_csv)[0]}.csv"
        ruta_diagrama_csv = os.path.join(ruta_persistencia, nombre_diagrama_csv)
        diagram_df.to_csv(ruta_diagrama_csv, index=False)

        # Gráfico de Diagrama de Persistencia
        plt.figure(figsize=(6, 6))
        gd.plot_persistence_diagram(diag)
        plt.title(f'Diagrama de Persistencia ({archivo_csv})')
        plt.xlabel('Birth')
        plt.ylabel('Death')

        nombre_imagen_persistencia = f"{os.path.splitext(archivo_csv)[0]}_diagrama_persistencia.png"
        ruta_imagen_persistencia = os.path.join(ruta_persistencia, nombre_imagen_persistencia)
        plt.tight_layout()
        plt.savefig(ruta_imagen_persistencia)
        plt.close()

    print(f'Todos los resultados se guardaron en: {ruta_persistencia}')
    return ruta_persistencia

if __name__ == "__main__":
    if len(sys.argv) > 1:
        ruta_centroides = sys.argv[1]
        calcular_rips_y_persistencia(ruta_centroides)
    else:
        print("Por favor, especifica la ruta de los centroides.")
