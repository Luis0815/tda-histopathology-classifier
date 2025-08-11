#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import gudhi as gd
import gudhi.wasserstein as gw
import sys
import re

# Grupos permitidos para clasificación
grupos_permitidos = ['tumorales', 'linfoides', 'mieloides', 'no_tumorales']

def obtener_grupo(archivo):
    """Extrae el grupo desde el nombre del archivo."""
    match = re.search(r'_(' + '|'.join(grupos_permitidos) + r')\.csv$', archivo)
    return match.group(1) if match else "Desconocido"

def distancias_por_grupo(ruta_directorio):
    """Calcula distancias por grupo y guarda resultados."""
    carpeta_salida = os.path.join(ruta_directorio, "distancias_por_grupo")
    os.makedirs(carpeta_salida, exist_ok=True)

    # Inicializar diccionarios para almacenar diagramas por grupo
    grupos = {g: {'diag_1': {}, 'diag_0': {}} for g in grupos_permitidos}

    # Listar todos los archivos CSV en el directorio
    archivos_csv = [f for f in os.listdir(ruta_directorio) if f.endswith('.csv')]

    # Leer y clasificar archivos CSV por grupo
    for archivo_csv in archivos_csv:
        try:
            ruta_completa = os.path.join(ruta_directorio, archivo_csv)
            df = pd.read_csv(ruta_completa, header=0)

            # Verificar columnas requeridas
            if 'birth' not in df.columns or 'death' not in df.columns:
                print(f"El archivo {archivo_csv} no contiene las columnas necesarias.")
                continue

            # Convertir columnas a float
            df['birth'] = df['birth'].astype(float)
            df['death'] = df['death'].astype(float)

            # Extraer diagramas para dimensiones 0 y 1
            diag_1 = df[df['dimension'] == 1][['birth', 'death']].to_numpy()
            diag_0 = df[df['dimension'] == 0][['birth', 'death']].to_numpy()

            # Asignar el archivo al grupo correspondiente
            grupo = obtener_grupo(archivo_csv)
            if grupo in grupos:
                grupos[grupo]['diag_1'][archivo_csv] = diag_1
                grupos[grupo]['diag_0'][archivo_csv] = diag_0
            else:
                print(f"Archivo {archivo_csv} no pertenece a un grupo reconocido.")

        except Exception as e:
            print(f"Error al procesar {archivo_csv}: {e}")

    # Función para calcular y guardar distancias por grupo
    def calcular_distancias(diagramas, diagramas_0, nombre_grupo):
        archivos = list(diagramas.keys())
        n = len(archivos)
        tolerancia = 1e-10

        # Inicializar matrices de distancias
        distancias_bottleneck_dim1 = np.zeros((n, n))
        distancias_bottleneck_dim0 = np.zeros((n, n))
        distancias_wasserstein_dim1 = np.zeros((n, n))
        distancias_wasserstein_dim0 = np.zeros((n, n))

        # Calcular distancias entre todos los pares de diagramas
        for i in range(n):
            for j in range(i, n):
                archivo_i = archivos[i]
                archivo_j = archivos[j]
                diag_i = diagramas[archivo_i]
                diag_j = diagramas[archivo_j]
                diag_i_0 = diagramas_0[archivo_i]
                diag_j_0 = diagramas_0[archivo_j]

                # Calcular distancias de Bottleneck
                distancia_bottleneck = gd.bottleneck_distance(diag_i, diag_j)
                distancia_bottleneck_0 = gd.bottleneck_distance(diag_i_0, diag_j_0)

                # Calcular distancias de Wasserstein
                distancia_wasserstein = gw.wasserstein_distance(diag_i, diag_j, order=1)
                distancia_wasserstein_0 = gw.wasserstein_distance(diag_i_0, diag_j_0, order=1)

                # Almacenar distancias en las matrices
                distancias_bottleneck_dim1[i, j] = distancias_bottleneck_dim1[j, i] = distancia_bottleneck
                distancias_bottleneck_dim0[i, j] = distancias_bottleneck_dim0[j, i] = distancia_bottleneck_0
                distancias_wasserstein_dim1[i, j] = distancias_wasserstein_dim1[j, i] = distancia_wasserstein
                distancias_wasserstein_dim0[i, j] = distancias_wasserstein_dim0[j, i] = distancia_wasserstein_0

        # Crear carpeta para resultados del grupo
        carpeta_grupo = os.path.join(carpeta_salida, nombre_grupo)
        os.makedirs(carpeta_grupo, exist_ok=True)

        # Guardar matrices de distancias por grupo
        pd.DataFrame(distancias_bottleneck_dim1, index=archivos, columns=archivos).to_csv(os.path.join(carpeta_grupo, f'distancias_bottleneck_dim1_{nombre_grupo}.csv'))
        pd.DataFrame(distancias_bottleneck_dim0, index=archivos, columns=archivos).to_csv(os.path.join(carpeta_grupo, f'distancias_bottleneck_dim0_{nombre_grupo}.csv'))
        pd.DataFrame(distancias_wasserstein_dim1, index=archivos, columns=archivos).to_csv(os.path.join(carpeta_grupo, f'distancias_wasserstein_dim1_{nombre_grupo}.csv'))
        pd.DataFrame(distancias_wasserstein_dim0, index=archivos, columns=archivos).to_csv(os.path.join(carpeta_grupo, f'distancias_wasserstein_dim0_{nombre_grupo}.csv'))

    # Calcular distancias para cada grupo
    for grupo in grupos_permitidos:
        if grupos[grupo]['diag_1']:
            print(f"Calculando distancias para {grupo}...")
            calcular_distancias(grupos[grupo]['diag_1'], grupos[grupo]['diag_0'], grupo)

    print(f"Distancias guardadas en {carpeta_salida}.")

# Ejecución principal
if __name__ == "__main__":
    if len(sys.argv) > 1:
        ruta_directorio = sys.argv[1]
        distancias_por_grupo(ruta_directorio)
    else:
        print("Por favor, especifica la ruta del directorio que contiene los diagramas de persistencia.")

