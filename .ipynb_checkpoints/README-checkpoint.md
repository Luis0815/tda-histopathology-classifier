# Topological-Analysis-of-Multicellular-Patterns
Código y datos para analizar patrones multicelulares mediante TDA. Se estudia la formación de estructuras por adhesión celular y, por separado, se aplica el análisis a datos histológicos de distintos estadios del cáncer para explorar su organización espacial.

### Objetivos Principales:
1. **Clasificar la disposición espacial** de distintos tipos celulares usando características topológicas extraídas a través del TDA.
2. **Aplicar TDA a diversos conjuntos de datos**, incluyendo etapas de progresión del cáncer, para investigar cómo la disposición celular en estas etapas se correlaciona con la progresión de la enfermedad.

### Metodología:
El pipeline incluye los siguientes pasos:
- Extracción de los centróides de las ubicaciones celulares a partir de imágenes de muestras de tejido.
- Cálculo de **complejos de Rips** y **homología persistente** a partir de estos datos espaciales.
- Análisis de diagramas de persistencia utilizando las distancias de ***Wasserstein** para comparar diferentes patrones.
- Visualización de los resultados para interpretar las relaciones espaciales entre tipos celulares en diferentes etapas del desarrollo tumoral.

### Conjuntos de Datos:
- Datos extraídos de imágenes de cáncer incluyendo carcinoma, displasia de alto grado y sus respectivas regiones de estroma adyacente.

---

### Estructura del Repositorio y Uso

El repositorio se organiza en carpetas que agrupan los scripts de análisis y los conjuntos de datos. A continuación se describen las carpetas principales:

#### 📁 `analysis_code/`  
Contiene scripts para el análisis de datos histológicos. Estos scripts permiten:
- Procesar coordenadas celulares contenidas en archivos `.csv`.  
- Aplicar TDA a muestras con hasta 18 tipos celulares distintos.  
- Clasificar la organización espacial de distintas condiciones (carcinoma, displasia y regiones adyacentes).  

Cada subcarpeta también contiene un `README.md` con instrucciones específicas de uso.

#### 📁 Carpeta de datos
Contiene carpetas con distintos archivos .csv con coordenadas celulares y metadatos por muestra.

Cada subcarpeta también contiene un `README.md` con información detallada de cada submuestra.

#### 📁 Carpeta de resultados
Contiene los resultados de cada proceso realizado a cada submuestra.

---


## Bibliografía

1. **Topological Data Analysis of Spatial Patterning in Heterogeneous Cell Populations: Clustering and Sorting with Varying Cell-Cell Adhesion**  
   [Nature Scientific Reports](https://www.nature.com/articles/s41540-023-00302-8)
   
2. **Persistent Homology Based Characterization of the Breast Cancer Immune Microenvironment: A Feasibility Study**  
   [Dagstuhl Reports](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2020.11)
   
3. **An Introduction to Topological Data Analysis: Fundamental and Practical Aspects for Data Scientists**  
   [Semantic Scholar](https://www.semanticscholar.org/reader/aff16209e232d38fc94a5b0c72067b88d106453f)
   
4. **Comparison of Persistence Diagrams**  
   [arXiv](https://ar5iv.labs.arxiv.org/html/2003.01352)
