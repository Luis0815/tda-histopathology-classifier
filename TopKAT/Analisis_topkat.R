library(TopKAT)
library(patchwork)
library(survival)
library(survminer)
library(TDAstats)
library(ggtda)
library(TopKAT)
library(ggplot2)
library(readr)
library(reshape2)
library(viridis)
library(dplyr)
library(stringr)
ruta <- "/home/jupyter-luisraul/R-studio_data_and_results/datos_para_R_carcinoma_HYN_nivel.csv"
# Cargar el archivo
data1.df <- read.csv(ruta)

# Verificar la estructura del dataframe
head(data1.df)  # Mostrar las primeras filas
str(data1.df)   # Ver estructura (variables y tipos de datos)


# Plotting some images from the first 50
p1 <- data1.df %>% dplyr::filter(PID == 1) %>% 
  ggplot(aes(x = x, y = y, colour = type)) + geom_point() +
  theme_bw() + ggtitle("Sample 1") +
  theme(legend.position = )

p2 <- data1.df %>% dplyr::filter(PID == 2) %>% 
  ggplot(aes(x = x, y = y, colour = type)) + geom_point() +
  theme_bw() + ggtitle("Sample 2")

p3 <- data1.df %>% dplyr::filter(PID == 12) %>% 
  ggplot(aes(x = x, y = y, colour = type)) + geom_point() +
  theme_bw() + ggtitle("Sample 51")

p4 <- data1.df %>% dplyr::filter(PID == 20) %>% 
  ggplot(aes(x = x, y = y, colour = type)) + geom_point() +
  theme_bw() + ggtitle("Sample 52")

# Arrange the plots
(p1 + p2 + p3 + p4 & theme(legend.position = "bottom")) + plot_layout(guides = "collect")


# Compute the similarity matrix
simmat <- rips_similarity_matrix(data1.df, max.threshold = 1000, print.progress = TRUE)



# Plotting the persistence diagrams
pd1 <- plot_persistence(simmat$rips.list[[1]], title = " Carcinoma_AGSCC_2_1.csv")

pd2 <- plot_persistence(simmat$rips.list[[8]], title = " Carcinoma_AGSCC_2_16.csv")

pd3 <- plot_persistence(simmat$rips.list[[33]], title = " Carcinoma_FAHNSCC_14_3.csv
")

pd4 <- plot_persistence(simmat$rips.list[[17]], title = "Carcinoma_FAHNSCC_15_2.csv")

# Arrange the plots
(pd1 + pd2 + pd3 + pd4 & theme(legend.position = "bottom")) + plot_layout(guides = "collect")


# Visualize the kernel matrices for the connected components
K.cc <- simmat$K.list[[1]] %>%
  reshape2::melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  theme_bw() +
  xlab("Image 1") + ylab("Image 2") +
  ggtitle("Kernel Matrix for Connected Components") +
  theme(legend.text = element_text(angle = 45, hjust = 1))


# Visualize the kernel matrices for the loops
K.loop <- simmat$K.list[[2]] %>%
  reshape2::melt() %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  theme_bw() +
  xlab("Image 1") + ylab("Image 2") +
  ggtitle("Kernel Matrix for Loops") +
  theme(legend.text = element_text(angle = 45, hjust = 1))

# Arrange the plots
K.cc + K.loop & theme(legend.position = "bottom")

#crear y
# Leer el archivo generado en Python
df <- read.csv("/home/jupyter-luisraul/R-studio_data_and_results/datos_para_R_carcinoma_HYN_nivel.csv")

# Obtener nombres únicos de archivo
filenames <- unique(df$archivo)  # asegúrate de usar el nombre correcto del dataframe

# Crear el vector binario según las condiciones del nombre
y <- ifelse(grepl("F", filenames, ignore.case = TRUE), 1, 0)

# Confirmar distribución
cat("📊 Distribución de grupos (0 = carcinoma, 1 = carcinoma+stroma):\n")
print(table(y, useNA = "ifany"))

# Verificar longitud
length(y)


length(y)       # debe ser igual a
#length(cens)    # y al número de filas/columnas de simmat$K.list[[1]]
dim(simmat$K.list[[1]])

kernel_ids <- rownames(simmat$K.list[[1]])
print(kernel_ids)



# Applying TopKAT to the simulated data solo binario enfermo y no 
res <- TopKAT(y = y, X = NULL, K.list = simmat$K.list,
              omega.list = c(0, 0.5, 1), outcome.type = "binary")


# Output the p-value
res$overall.pval

res$p.vals


# Plot image and persistence diagram side-by-side
p3 + pd3


# View the birth/death distances for the loops
tail(simmat$rips.list[[17]])


# Save the corresponding persistence diagram
pd <- as.data.frame(simmat$rips.list[[17]])

# Threshold at t = 5
pd3.1 <- plot_persistence(pd %>% dplyr::filter(birth <= 10 & death <= 10), 
                          title = expression(paste("ε = 10")),
                          dims = c(15, 15))

# Threshold at t=10
pd3.2 <- plot_persistence(pd %>% dplyr::filter(birth <= 15 & death <= 15), 
                          title =  expression(paste("ε = 15")),
                          dims = c(15, 15))


# Threshold at t=15
pd3.3 <- plot_persistence(pd %>% dplyr::filter(birth <= 20 & death <= 20), 
                          title =  expression(paste("ε = 20")),
                          dims = c(15, 15))


# Threshold at t=20
pd3.4 <- plot_persistence(pd %>% dplyr::filter(birth <= 25 & death <= 25), 
                          title =  expression(paste("ε = 25")),
                          dims = c(15, 15))

# Arrange
(pd3.1 + pd3.2 + pd3.3 + pd3.4 & theme(legend.position = "bottom")) + plot_layout(guides = "collect")


res_scale_import <- scale_importance(pd.list = simmat$rips.list,
                                     y = y, 
                                     omega.list = c(0, 0.5, 1),
                                     threshold = 100, 
                                     PIDs = 1:84,
                                     outcome.type = "binary")

# Create a data.frame
res_scale_import.df <- data.frame(
  thresh = res_scale_import$threshold.seq,
  pval = res_scale_import$pvals
)

# Plot
res_scale_import.df %>% 
  ggplot(aes(x = thresh, y = pval)) +
  geom_point() +
  theme_bw() +
  xlab(expression(epsilon)) +
  ylab(expression(p ~ "-valor")) +
  ggtitle("") +
  geom_vline(xintercept = res_scale_import$min.thresh, linetype = "dashed") +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10), # top, right, bottom, left
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )


# Plot the simplicial complex at r=res_scale_import$min.thresh
sc1 <- plot_cells_with_scale(
  image = data1.df %>% dplyr::filter(PID == 1),
  threshold = res_scale_import$min.thresh,
  title = "Simplicial Complex for \n Sample 1"
)

sc2 <- plot_cells_with_scale(
  image = data1.df %>% dplyr::filter(PID == 2),
  threshold = res_scale_import$min.thresh,
  title = "Simplicial Complex for \n Sample 2"
)

sc3 <- plot_cells_with_scale(
  image = data1.df %>% dplyr::filter(PID == 16),
  threshold = res_scale_import$min.thresh,
  title = "Simplicial Complex for \n Sample 16"
)

sc4 <- plot_cells_with_scale(
  image = data1.df %>% dplyr::filter(PID == 24),
  threshold = res_scale_import$min.thresh,
  title = "Simplicial Complex for \n Sample 20"
)


# Arrange the plots

sc1 + sc2 + sc3 + sc4
sc4

# Connectivity matrices
c1 <- plot_cell_connections(
  image = data1.df %>% dplyr::filter(PID == 1),
  threshold = res_scale_import$min.thresh,
  title = "Cell Connectivity for Patient 1",
  type.column = "type",
  unique.types = unique(data1.df$type)
)

c2 <- plot_cell_connections(
  image = data1.df %>% dplyr::filter(PID == 2),
  threshold = res_scale_import$min.thresh,
  title = "Cell Connectivity for Patient 2",
  type.column = "type",
  unique.types = unique(data1.df$type)
)

c3 <- plot_cell_connections(
  image = data1.df %>% dplyr::filter(PID == 16),
  threshold = res_scale_import$min.thresh,
  title = "Cell Connectivity for Patient 16",
  type.column = "type",
  unique.types = unique(data1.df$type)
)

c4 <- plot_cell_connections(
  image = data1.df %>% dplyr::filter(PID == 20),
  threshold = res_scale_import$min.thresh,
  title = "Cell Connectivity for Patient 20",
  type.column = "type",
  unique.types = unique(data1.df$type)
)

# Arrange the plots
c1 + c2 + c3 + c4
c1
c2
c3
c4

res_scale_import$min.thresh

NUMBER <- 60

# Obtener el nombre del archivo asociado al paciente con PID 16
nombre_archivo <- data1.df %>%
  dplyr::filter(PID == NUMBER) %>%
  dplyr::pull(archivo) %>%
  unique()
c7 <- plot_cell_connections(
  image = data1.df %>% dplyr::filter(PID == NUMBER ),
  threshold = res_scale_import$min.thresh,
  title = paste0("Cell Connectivity (", nombre_archivo, ")"),
  type.column = "type",
  unique.types = unique(data1.df$type)
)
c7




# === Carpeta de salida ===
output_dir <- file.path(getwd(), "resultados_topkat_carcinomaFA_VS_NOFA")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# === Parámetros / checks ===
if (!exists("data1.df")) stop("data1.df no existe. Carga tu dataframe primero.")
if (!exists("res_scale_import")) stop("res_scale_import no existe. Ejecuta scale_importance antes.")

cell.types <- as.character(unique(data1.df$type))
threshold <- res_scale_import$min.thresh

# === Asociar a cada PID su nombre de archivo (tomando el primer archivo si hay >1) ===
pid_list <- unique(data1.df$PID)
pid_archivo <- sapply(pid_list, function(pid) {
  archivos <- unique(as.character(data1.df$archivo[data1.df$PID == pid]))
  if (length(archivos) == 0) return(NA_character_)
  return(archivos[1])
}, USE.NAMES = FALSE)

pid_df <- data.frame(PID = pid_list, archivo = pid_archivo, stringsAsFactors = FALSE)

# === Asignar grupo: 0 = carcinoma_FA, 1 = stroma_ad_carcinoma_FA, NA = otros/no-FA ===
# === Asignar grupo: 1 = carcinoma FA, 0 = carcinoma noFA, NA = otros ===
pid_df$grupo <- NA_integer_

for (i in seq_len(nrow(pid_df))) {
  name <- tolower(pid_df$archivo[i])
  has_carcinoma <- str_detect(name, "carcinoma")
  is_f <- str_detect(name, "f")  # si contiene 'F' -> FA
  
  if (has_carcinoma && is_f) {
    pid_df$grupo[i] <- 1   # carcinoma FA
  } else if (has_carcinoma && !is_f) {
    pid_df$grupo[i] <- 0   # carcinoma noFA
  } else {
    pid_df$grupo[i] <- NA  # excluir stroma u otros
  }
}

# Filtrar solo PIDs válidos
pid_df <- pid_df %>% filter(!is.na(grupo))

cat("Conteo por grupo (0=noFA, 1=FA):\n")
print(table(pid_df$grupo))


# Filtrar solo PIDs FA relevantes
pid_df <- pid_df %>% filter(!is.na(grupo))
if (nrow(pid_df) == 0) stop("No se encontraron muestras FA 'carcinoma' o 'stroma-ad carcinoma' con el patrón esperado.")

cat("PIDs usados:\n")
print(table(pid_df$grupo))

# === Inicializar matrices de suma ===
connect_carcinomaFA <- matrix(0, nrow = length(cell.types), ncol = length(cell.types),
                              dimnames = list(cell.types, cell.types))
connect_stromaFA <- connect_carcinomaFA

# === Iterar PIDs y sumar matrices ===
for (i in seq_len(nrow(pid_df))) {
  pid <- pid_df$PID[i]
  grupo <- pid_df$grupo[i]
  patient <- data1.df %>% filter(PID == pid)
  
  if (nrow(patient) == 0) {
    warning(sprintf("PID %s sin filas en data1.df — se salta.", pid))
    next
  }
  
  connect_i <- tryCatch({
    generate_connectivity(images.df = patient,
                          threshold = threshold,
                          type.column = "type",
                          unique.types = cell.types)
  }, error = function(e) {
    warning(sprintf("Error generate_connectivity PID %s: %s", pid, e$message))
    return(NULL)
  })
  if (is.null(connect_i)) next
  
  # asegurar orden/llenado en caso de tipos faltantes
  tmp <- matrix(0, nrow = length(cell.types), ncol = length(cell.types), dimnames = list(cell.types, cell.types))
  rmatch <- match(rownames(connect_i), cell.types)
  cmatch <- match(colnames(connect_i), cell.types)
  tmp[rmatch, cmatch] <- connect_i
  
  if (grupo == 0) {
    connect_carcinomaFA <- connect_carcinomaFA + tmp
  } else if (grupo == 1) {
    connect_stromaFA <- connect_stromaFA + tmp
  }
}

# === Promediar ===
n_carcinomaFA <- sum(pid_df$grupo == 0)
n_stromaFA <- sum(pid_df$grupo == 1)

if (n_carcinomaFA == 0) stop("No hay muestras carcinoma FA.")
if (n_stromaFA == 0) stop("No hay muestras stroma-ad carcinoma FA.")

connect_carcinomaFA <- connect_carcinomaFA / n_carcinomaFA
connect_stromaFA <- connect_stromaFA / n_stromaFA

# === Visualización (heatmaps) ===
plot_connectivity_matrix <- function(connect, title) {
  ggplot(reshape2::melt(connect), aes(Var1, Var2, fill = value)) +
    geom_tile(colour = "white") +
    viridis::scale_fill_viridis(option = "turbo") +
    labs(x = "Tipo celular 1", y = "Tipo celular 2", fill = "Conexiones") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 12)) +
    ggtitle(title)
}

p_carcinomaFA <- plot_connectivity_matrix(connect_carcinomaFA, paste0("Average Connectivity (Carcinoma FA) — n=", n_carcinomaFA))
p_stromaFA <- plot_connectivity_matrix(connect_stromaFA, paste0("Average Connectivity (Stroma-ad Carcinoma FA) — n=", n_stromaFA))

# === Guardar SVGs ===
ggsave(file.path(output_dir, "connectivity_carcinomaFA.svg"), p_carcinomaFA, device = "svg", width = 8, height = 6)
ggsave(file.path(output_dir, "connectivity_stromaFA.svg"), p_stromaFA, device = "svg", width = 8, height = 6)
ggsave(file.path(output_dir, "connectivity_combined_carcinomaFA_vs_stromaFA.svg"), p_carcinomaFA + p_stromaFA, device = "svg", width = 12, height = 6)

cat("✅ SVGs guardados en:", output_dir, "\n")

# === PDF con una página por PID (en la misma carpeta) ===
pdf(file.path(output_dir, "Reporte_Conectividad_carcinomaFA_vs_NOFA.pdf"), width = 8, height = 6)

# Portada/resumen
plot.new()
title(main = "Reporte: Carcinoma FA vs Carcinoma no-FA", cex.main = 1.2)
text(0, 0.9, paste("Total PIDs usados:", nrow(pid_df)), adj = 0)
text(0, 0.85, paste("Carcinoma FA:", n_carcinomaFA), adj = 0)
text(0, 0.80, paste("Stroma-ad FA:", n_stromaFA), adj = 0)
if (exists("res")) {
  text(0, 0.72, paste("TopKAT p-valor global:", signif(res$overall.pval, 4)), adj = 0)
}
if (exists("res_scale_import")) {
  text(0, 0.68, paste("Threshold óptimo:", res_scale_import$min.thresh), adj = 0)
}

# Página con mapas promedio
print(p_carcinomaFA + p_stromaFA)

# Páginas por PID (solo PIDs incluidos en pid_df)
for (pid in pid_df$PID) {
  nombre_archivo <- unique(data1.df$archivo[data1.df$PID == pid])
  titulo <- paste0("PID ", pid, " — ", nombre_archivo)
  p <- plot_cell_connections(
    image = data1.df %>% filter(PID == pid),
    threshold = threshold,
    title = titulo,
    type.column = "type",
    unique.types = cell.types
  )
  print(p)
}


dev.off()
cat("✅ PDF guardado en:", file.path(output_dir, "Reporte_Conectividad_carcinomaFA_vs_stromaFA.pdf"), "\n")
