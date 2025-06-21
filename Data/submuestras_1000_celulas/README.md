#  Submuestras de Tejido Celular (1000 células por muestra)

Esta carpeta contiene los scripts para aplicar Análisis Topológico de Datos (TDA) a muestras histológicas, utilizando coordenadas celulares extraídas de imágenes. El análisis se realiza a partir de archivos `.csv` que contienen las posiciones de las células y sus respectivos fenotipos.

---

## Estructura de los archivos

Cada archivo `.csv` contiene las coordenadas espaciales de las células, junto con información adicional como tipo celular 

| X_centroid | Y_centroid | phenotype           |
|------------|------------|---------------------|
| 123.4      | 567.8      |     tumor cells     | 
| 234.5      | 678.9      | 	M2 macrophages  |
| ...        | ...        | ...                 |

---

## Notas

- Todas las submuestras tienen un tamaño estandarizado de ~3000 células para permitir comparaciones consistentes en los análisis de TDA.