
library(dplyr)
library(tidyr)
library(stringr)
library(future)
library(future.apply)

SSN_bySample.dir <- "/STORAGE/csbig/anakamura/TALL_Ranked_NoNormal"

# Crear el directorio si no existe
if (!dir.exists(SSN_bySample.dir)) {
  dir.create(SSN_bySample.dir)
}

SSN_file <- "/STORAGE/csbig/anakamura/TALL_SSN_NoNormal"
SSNs <- list.files(SSN_file)

# Función para procesar cada muestra
Rank_and_Bind <- function(i, EdgeNum = 10000) {
  
  # Inicializar un dataframe vacío para almacenar el top por muestra
  Top_sample <- data.frame(Edge = character(), Z_score_ = numeric(), stringsAsFactors = FALSE)
  
  for (doc in SSNs) {
    
    #print(which(SSNs == doc))
    
    file_path <- file.path(SSN_file, doc)
    
    # Leer los nombres de las columnas de este archivo específico
    column_names <- strsplit(readLines(file_path, n = 1), ",")[[1]]  # Leer solo la primera línea para los nombres de las columnas
    
    # Leer solo la fila correspondiente a la muestra de interés (i-ésima fila)
    sample_line <- readLines(file_path, n = i + 1)[i + 1]  # Leer la línea que contiene la muestra de interés (i + 1 porque la primera fila son los colnames)
    
    # Convertir la línea leída en un vector de valores numéricos
    sample_values <- as.numeric(unlist(strsplit(sample_line, ",")))  # Convertir la línea en un vector numérico
    
    # Crear un dataframe con los nombres de las columnas
    SSNn <- as.data.frame(t(sample_values), stringsAsFactors = FALSE)
    colnames(SSNn) <- column_names  # Asignar los nombres de las columnas
    
    # Convertir de ancho a largo y seleccionar los enlaces más grandes
    SSNn <- SSNn %>%
      pivot_longer(cols = everything(), names_to = "Edge", values_to = "Z_score_") %>%
      arrange(desc(abs(Z_score_))) %>%
      head(EdgeNum)  # Tomar los 10,000 enlaces más grandes por archivo
    
    # Obtener el valor mínimo actual en Top_sample
    if (nrow(Top_sample) > 0) {
      min_z_score <- min(abs(Top_sample$Z_score_))
    } else {
      min_z_score <- 0  # Si Top_sample está vacío, tomamos min_z_score como 0
    }
    
    # Filtrar los enlaces de SSNn que son mayores que el mínimo en Top_sample
    better_links <- SSNn %>% filter(abs(Z_score_) > min_z_score)
    
    if (nrow(better_links) > 0) {
      # Hacer bind de Top_sample con los better_links
      Top_sample <- bind_rows(Top_sample, better_links)
      
      # Ordenar por valor absoluto y seleccionar los top 10,000 enlaces
      Top_sample <- Top_sample %>%
        arrange(desc(abs(Z_score_))) %>%
        head(EdgeNum)  # Mantener solo los 10,000 enlaces más grandes
    }
  }
  
  # Mutar y guardar el archivo final
  Top_sample <- Top_sample %>%
    mutate(Source = sub("_.*", "", Edge),  # Palabra antes del guion bajo
           Target = sub(".*_", "", Edge))
  
  # Guardar el archivo con el top 10,000 enlaces para la muestra i
  output_file <- file.path(SSN_bySample.dir, paste0(i, "_SSN.csv"))
  write.csv(Top_sample, output_file, row.names = FALSE)
}

# Planificación paralela con 40 núcleos
future::plan(multisession, workers = 40)

# Procesar las muestras (puedes ajustar el rango según tus datos)
results <- future_lapply(1:243, FUN = Rank_and_Bind, future.seed = TRUE)

# Volver a la ejecución secuencial
plan(sequential)
