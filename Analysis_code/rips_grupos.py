import os
import sys
import numpy as np
import pandas as pd
import gudhi as gd
import matplotlib.pyplot as plt
from tqdm import tqdm

# Definir los grupos celulares
grupos = {
    'tumorales': ['tumor cells', 'Ki67+ tumor cells'],  
    'linfoides': ['NK', 'B cells', 'effector CD8+ T cells', 'memory CD8+ T cells', 'CD4+ T cells', 'regulatory T cells', 'memory CD4+ T cells'],
    'mieloides': ['neutrophils', 'other APCs', 'dendritic cells', 'M1/M0 macrophages', 'M2 macrophages'],
    'no_tumorales': ['endothelial cells', 'stromal cells']
}

def calcular_rips_y_persistencia(ruta_centroides):
    # Crear carpeta para guardar resultados
    ruta_persistencia = os.path.join(ruta_centroides, "persistencia_grupos_80")
    os.makedirs(ruta_persistencia, exist_ok=True)

    # Obtener lista de archivos CSV
    archivos_csv = [f for f in os.listdir(ruta_centroides) if f.endswith(".csv")]

    # Iterar sobre cada archivo de centroides con tqdm
    for archivo_csv in tqdm(archivos_csv, desc="Procesando archivos", unit="archivo"):
        # Construir la ruta completa del archivo
        ruta_completa = os.path.join(ruta_centroides, archivo_csv)
        
        # Leer el archivo CSV
        df = pd.read_csv(ruta_completa)

        # Iterar sobre cada grupo con tqdm
        for grupo, tipos_celulares in grupos.items():
            # Filtrar las filas que corresponden a las células del grupo actual
            df_grupo = df[df['phenotype_key'].isin(tipos_celulares)]

            # Si el grupo tiene datos
            if df_grupo.empty:
                continue
            
            # Extraer las coordenadas de los centroides para el grupo
            centroides_x = df_grupo['X_centroid'].tolist()
            centroides_y = df_grupo['Y_centroid'].tolist()
            
            # Convertir las coordenadas de los centroides a un formato adecuado para Gudhi
            puntos = np.array(list(zip(centroides_x, centroides_y)))

            # Calcular el complejo de Rips con un radio específico
            rips_complex = gd.RipsComplex(points=puntos, max_edge_length=80)
            simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)  # Ahora con dimensión superior

            # Visualizar el complejo de Rips y el diagrama de persistencia
            plt.figure(figsize=(36, 15))
            plt.scatter(centroides_x, centroides_y, color='black', label='Centroides')
            for simplex in simplex_tree.get_skeleton(1):  # Obtener las aristas
                if len(simplex[0]) == 2:  # Solo aristas (simplexes de dimensión 1)
                    arista = simplex[0]
                    x = [centroides_x[i] for i in arista]
                    y = [centroides_y[i] for i in arista]
                    plt.plot(x, y, color='gray', linestyle='-', linewidth=1)
            plt.xlabel('Coordenada X')
            plt.ylabel('Coordenada Y')
            plt.title(f'Complejo de Rips ({archivo_csv}, {grupo})')
            plt.legend()

            # Guardar la imagen del complejo de Rips
            nombre_imagen_rips = f"{os.path.splitext(archivo_csv)[0]}_{grupo}_complejo_rips.png"
            ruta_imagen_rips = os.path.join(ruta_persistencia, nombre_imagen_rips)
            plt.tight_layout()
            plt.savefig(ruta_imagen_rips)
            plt.close()

            # Calcular el diagrama de persistencia
            diag = simplex_tree.persistence()

            # Guardar los datos de persistencia en un DataFrame
            persistencia_aplanada = []
            for d in diag:
                dimension, (birth, death) = d
                persistencia_aplanada.append([dimension, birth, death])

            diagram_df = pd.DataFrame(persistencia_aplanada, columns=['dimension', 'birth', 'death'])

            # Guardar el diagrama de persistencia en CSV
            nombre_diagrama_csv = f"{os.path.splitext(archivo_csv)[0]}_{grupo}.csv"
            ruta_diagrama_csv = os.path.join(ruta_persistencia, nombre_diagrama_csv)
            diagram_df.to_csv(ruta_diagrama_csv, index=False)

            # Subplot 2: Diagrama de Persistencia
            plt.figure(figsize=(6, 6))
            gd.plot_persistence_diagram(diag)
            plt.title(f'Diagrama de Persistencia ({archivo_csv}, {grupo})')
            plt.xlabel('Birth')
            plt.ylabel('Death')

            # Guardar la imagen del diagrama de persistencia
            nombre_imagen_persistencia = f"{os.path.splitext(archivo_csv)[0]}_{grupo}_diagrama_persistencia.png"
            ruta_imagen_persistencia = os.path.join(ruta_persistencia, nombre_imagen_persistencia)
            plt.tight_layout()
            plt.savefig(ruta_imagen_persistencia)
            plt.close()

            print(f'Guardado: {nombre_imagen_persistencia} en {ruta_persistencia}')

    return ruta_persistencia

if __name__ == "__main__":
    if len(sys.argv) > 1:
        ruta_centroides = sys.argv[1]
        calcular_rips_y_persistencia(ruta_centroides)
    else:
        print("Por favor, especifica la ruta de los centroides.")
