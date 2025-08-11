# Análisis Topológico de Distribución Celular

Este documento presenta un flujo metodológico integral para analizar la organización espacial y topológica de células en muestras biológicas, utilizando herramientas de Análisis Topológico de Datos (TDA). El enfoque se basa en:

- La construcción de complejos de Rips a partir de datos de centroides celulares.
- El cálculo de distancias topológicas (Wasserstein y Bottleneck) entre diagramas de persistencia.
- La visualización mediante clustermaps jerárquicos con anotaciones biológicas relevantes.

Se realizan tres análisis principales con distintos niveles de detalle y enfoque:

1. Análisis general considerando todos los núcleos celulares.
2. Análisis segmentado por grupos celulares específicos (tumorales, linfoides, mieloides y no tumorales).
3. Análisis por combinaciones de estos grupos para explorar interacciones topológicas entre subpoblaciones celulares.

Cada análisis sigue un flujo en tres etapas: generación de diagramas, cálculo de distancias y visualización.

---

A continuación, se describen en detalle los pasos y scripts correspondientes a cada uno de los análisis realizados.

<details>
<summary><strong>Análisis 1 – Núcleos</strong></summary>

Este análisis permite obtener complejos de Rips, calcular distancias de Wasserstein y generar clustermaps jerárquicos con anotaciones para interpretar relaciones entre muestras.

#### Flujo general

1. **Generar Rips y diagramas de persistencia** (`rips.py`)
2. **Calcular distancias de Wasserstein** (`distancias_wasserstein.py`)
3. **Visualizar resultados en clustermaps** (`clustermap_multiple.py`)

---

#### 1.1 Generar Rips y Diagramas de Persistencia

**Script:** `rips.py`  
**Entrada:** Carpeta con archivos `.csv` que contengan las columnas `X_centroid` y `Y_centroid`.  
**Salida:**
- Imagen del complejo de Rips por muestra.
- Imagen del diagrama de persistencia por muestra.
- Archivo CSV con las coordenadas `birth` y `death`.

**Comando de ejecución:**
```bash
python rips.py /ruta/a/centroides --radio 1000 --workers 4
```

**Los resultados se guardarán en:**
```swift
/ruta/a/centroides/resultados/rips_<radio>/
```

#### 1.2. Calcular Distancias de Wasserstein  

**Script:** `distancias_wasserstein.py`  
**Entrada:** Carpeta `rips_<radio>` generada en el paso anterior.  
**Salida:**  
- `wasserstein_dim0.csv`  
- `wasserstein_dim1.csv`  
*(Matrices de distancias para dimensión 0 y 1.)*  

**Comando de ejecución:**
```bash
python distancias_wasserstein.py /ruta/a/centroides/resultados/rips_1000 --workers 4
```

**Los resultados se guardarán en:**
```bash
/ruta/a/centroides/resultados/distancias_wasserstein/
```

#### 1.3. Generar Clustermaps con Anotaciones  

**Script:** `clustermap_multiple.py`  
**Entrada:** Matrices de distancia (`wasserstein_dim0.csv` o `wasserstein_dim1.csv`).  
**Salida:** Clustermaps agrupados por:  
- Tipo de muestra.  
- Estado Fanconi.  
- Origen anatómico.  

**Comando de ejecución:**
```bash
python clustermap_multiple.py /ruta/a/centroides/resultados/distancias_wasserstein --metodo ward
```

**Las imágenes se guardarán en subcarpetas:**
```bash
visualizacion/por_tipo/
visualizacion/por_fanconi/
visualizacion/por_origen/
```
---

</details>

<details>
<summary><strong>Análisis 2 – Grupos Celulares</strong></summary>

Este análisis genera complejos de Rips y diagramas de persistencia para subconjuntos celulares específicos (tumorales, linfoides, mieloides y no tumorales), calcula distancias de Wasserstein por grupo, y genera clustermaps jerárquicos con anotaciones para interpretar relaciones entre muestras por grupo celular.

### Flujo general

1. **Generar Rips y diagramas de persistencia por grupo celular** (`rips_grupos.py`)  
2. **Calcular distancias de Wasserstein (y opcionalmente Bottleneck) por grupo celular** (`distancias_grupos.py`)  
3. **Visualizar resultados en clustermaps** (`clustermap_multiple.py`)

---

#### 2.1 Generar Rips y Diagramas de Persistencia por Grupo Celular

**Script:** `rips_grupos.py`  
**Entrada:** Carpeta con archivos `.csv` que contengan las columnas `X_centroid`, `Y_centroid` y `phenotype`.  
**Salida:**  
Por cada archivo y grupo celular (tumorales, linfoides, mieloides, no tumorales):
- Imagen del complejo de Rips (`<archivo>_<grupo>_complejo_rips.png`)
- Imagen del diagrama de persistencia (`<archivo>_<grupo>_diagrama_persistencia.png`)
- Archivo CSV con las coordenadas `birth` y `death` (`<archivo>_<grupo>.csv`)

**Comando de ejecución:**
```bash
python rips_grupos.py /ruta/a/csvs --radio 1000 --workers 4
```

**Los resultados se guardarán en:**
```bash
/ruta/a/centroides/resultados/distancias_wasserstein/
```

#### 2.2 Calcular Distancias de Wasserstein y Bottleneck por Grupo Celular

**Script:** `distancias_grupos.py`  
**Entrada:** Carpeta con archivos CSV generados en el paso anterior (`*_grupo.csv`).  
**Salida:** Para cada grupo celular:

- `distancias_wasserstein_dim0_<grupo>.csv`
- `distancias_wasserstein_dim1_<grupo>.csv`
- *(Opcional)* `distancias_bottleneck_dim0_<grupo>.csv`
- *(Opcional)* `distancias_bottleneck_dim1_<grupo>.csv`

**Comando de ejecución:**
```bash
python distancias_grupos.py /ruta/a/diagramas --workers 4 [--bottleneck]
```
**Los resultados se guardarán en:**
```bash
/ruta/a/diagramas/resultados/distancias_grupos/
```

#### 2.3 Generar Clustermaps con Anotaciones

**Script:** `clustermap_multiple.py`  
**Entrada:** Matrices de distancia (`distancias_wasserstein_dim0_<grupo>.csv`, etc.).  
**Salida:** Clustermaps agrupados por:  
- Tipo de muestra  
- Estado Fanconi  
- Origen anatómico  

**Comando de ejecución:**
```bash
python clustermap_multiple.py /ruta/a/diagramas/resultados/distancias_grupos --metodo ward
```
---

</details>

<details>
<summary><strong>Análisis 3: Comparación Topológica por Combinaciones de Grupos Celularess</strong></summary>
  


En este análisis se estudia la estructura topológica de la distribución espacial de células considerando **combinaciones específicas de grupos celulares** (tumorales, linfoides, mieloides y no tumorales).  
A diferencia de los análisis anteriores que trabajaban con **todas las células** o con **un único grupo**, aquí se generan subconjuntos de datos para cada combinación de grupos, y se repite el flujo de análisis topológico para cada caso.

El objetivo es identificar patrones diferenciales de conectividad y persistencia entre las combinaciones de grupos, evaluando si ciertas asociaciones celulares producen cambios estructurales relevantes en la topología.

---

### Flujo General del Análisis
1. **Generar Rips y diagramas de persistencia por cada combinación de grupos celulares** (`rips_grupos_combinaciones.py`)  
2. **Calcular distancias de Wasserstein (y opcionalmente Bottleneck) por combinaciones de grupos celulares ** (`distancias_grupos.py`)  
3. **Visualizar resultados en clustermaps** (`clustermap_multiple.py`)

#### 3.1 Generar Complejos de Rips y Diagramas de Persistencia para Combinaciones de Grupos

**Script:** `rips_combinaciones.py`  
**Descripción:** Genera complejos de Rips y diagramas de persistencia para combinaciones de grupos celulares (tumorales, linfoides, mieloides y no tumorales) a partir de archivos CSV con centroides.

**Entrada:** Carpeta con archivos `.csv` con columnas:  
- `X_centroid`  
- `Y_centroid`  
- `phenotype`  

**Salida:** Para cada archivo y combinación de grupos, se crea una carpeta con el nombre de la combinación, dentro de la cual se generan los archivos: 
- `<archivo>_<grupo_comb>_complejo_rips.png`  
- `<archivo>_<grupo_comb>_diagrama_persistencia.png`  
- `<archivo>_<grupo_comb>.csv` (tabla birth–death)  

**Comando de ejecución:**
```bash
python rips_grupos_combinaciones.py /ruta/a/csvs --radio 1000 --workers 2
```
**Los resultados se guardarán en:**
```bash
/ruta/a/csvs/resultados/rips_combinaciones_<radio>/
```


#### 3.2 Calcular Distancias de Wasserstein y Bottleneck por Combinación de Grupos

**Script:** `distancias_wasserstein.py`  
**Entrada:** Cada subcarpeta generada por `rips_combinaciones.py` (contiene archivos `.csv` con tablas birth–death).  

**Salida:**  
- `distancias_wasserstein_dim0.csv`  
- `distancias_wasserstein_dim1.csv`  
- (Opcional) `distancias_bottleneck_dim0.csv`  
- (Opcional) `distancias_bottleneck_dim1.csv`  

**Comando de ejecución (ejemplo para una carpeta):**
```bash
python distancias_wasserstein.py /ruta/a/csvs/resultados/rips_combinaciones_<radio>/<subcarpeta> --workers 4
```

#### 3.3 Generar Clustermaps por Combinación de Grupos

**Script:** `clustermap_multiple.py`  
**Entrada:** Matrices de distancia (`distancias_wasserstein_dim0.csv`, etc.) de cada subcarpeta procesada.  

**Salida:** Clustermaps agrupados por:  
- Tipo de muestra  
- Estado Fanconi  
- Origen anatómico  

**Comando de ejecución (ejemplo para una carpeta):**
```bash
python clustermap_multiple.py /ruta/a/csvs/resultados/rips_combinaciones_<radio>/<subcarpeta> --metodo ward
```

---
</details>
