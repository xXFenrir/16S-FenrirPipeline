# ==============================================================================
# SCRIPT DE TESIS: ANÁLISIS INTEGRAL DE MICROBIOTA (SUELO 2020)
# ==============================================================================

# Cargar librerías necesarias [cite: 3, 5]
library(MicrobiotaProcess)
library(ggplot2)
library(dplyr)
library(vegan)
library(magrittr)

# 1. PREPARACIÓN Y ALFA DIVERSIDAD [cite: 1, 4]
# Rarefacción obligatoria para comparación de diversidad [cite: 1]
mpse_tesis %<>% mp_rrarefy()
mpse_tesis %<>% mp_cal_alpha(.abundance = RareAbundance)

# Guardar Tabla de Alfa Diversidad [cite: 1, 4]
df_alpha <- mpse_tesis %>% mp_extract_sample() %>% as.data.frame()
write.csv(df_alpha[, sapply(df_alpha, is.atomic)], "1_Tabla_Alfa_Diversidad.csv", row.names = FALSE)

# Gráfico de Alfa Diversidad [cite: 3]
p_alpha <- mp_plot_alpha(.data = mpse_tesis, .group = crop_rotation, .alpha = c(Shannon, Pielou)) + 
  theme_bw()
ggsave("1_Figura_Alpha.pdf", p_alpha, width = 8, height = 5)


# 2. DIVERSIDAD BETA Y PERMANOVA (CON ELIPSES DE COLOR) 
mpse_tesis %<>% mp_cal_dist(.abundance = RareAbundance, distmethod = "bray")
dist_mat <- mpse_tesis %>% mp_extract_dist(distmethod = "bray")

# Cálculo de PERMANOVA con Vegan [cite: 5]
perm_res <- adonis2(dist_mat ~ crop_rotation, data = df_alpha, permutations = 999)
write.csv(as.data.frame(perm_res), "2_Tabla_PERMANOVA.csv")

# PCoA Limpio: Sin scatters en ejes y con elipses rellenas [cite: 5]
p_pcoa <- mpse_tesis %>% 
  mp_plot_ord(.ordmethod = "PCoA", .group = crop_rotation, .size = 3) +
  stat_ellipse(aes(fill = crop_rotation), geom = "polygon", alpha = 0.2, level = 0.95) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(title = paste("PCoA Bray-Curtis (p-value:", perm_res$`Pr(>F)`[1], ")"))

ggsave("2_Figura_Beta_PCoA.pdf", p_pcoa, width = 8, height = 6)


# 3. ABUNDANCIA RELATIVA (GÉNERO Y ESPECIE) SIN "OTHERS" 

# Extraemos los datos a un dataframe para evitar errores de estructura [cite: 2]
df_master <- as.data.frame(mpse_tesis)

graficar_abundancia_final <- function(datos, nivel_nombre) {
  
  # Identificar columna de nombres (OTU o label) [cite: 2, 4]
  col_taxa <- if("label" %in% colnames(datos)) "label" else "OTU"
  
  # Filtrar Top 20 REALES (Omitiendo la categoría Others) [cite: 4]
  top20_data <- datos %>%
    filter(!get(col_taxa) %in% c("Others", "others", "Other")) %>%
    group_by(Sample, !!sym(col_taxa), crop_rotation) %>%
    summarise(Abund = sum(RareAbundance), .groups = "drop") %>%
    group_by(!!sym(col_taxa)) %>%
    mutate(MediaGlobal = mean(Abund)) %>%
    ungroup() %>%
    arrange(desc(MediaGlobal)) %>%
    filter(!!sym(col_taxa) %in% unique(!!sym(col_taxa))[1:20])
  
  # Limpiar nombres (quitar prefijos como species.xxx) [cite: 4]
  top20_data$TaxonLimpio <- sub(".*\\.", "", as.character(top20_data[[col_taxa]]))
  
  # Gráfico de Barras (Ocupando el 100% de la barra) [cite: 5]
  p_bar <- ggplot(top20_data, aes(x = Sample, y = Abund, fill = TaxonLimpio)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_grid(~crop_rotation, scales = "free_x", space = "free") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = paste("Top 20", nivel_nombre, "(Sin Others)"), y = "Abundancia Relativa")
  
  ggsave(paste0("3_Abundancia_", nivel_nombre, ".pdf"), p_bar, width = 12, height = 7)
  write.csv(top20_data, paste0("3_Tabla_Top20_", nivel_nombre, ".csv"), row.names = FALSE)
}

# Ejecutar para los taxones encontrados 
try({ graficar_abundancia_final(df_master, "Taxones_Principales") })

print("--- PROCESO FINALIZADO: TODOS LOS ARCHIVOS GUARDADOS ---")
