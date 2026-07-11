# ── 0. LIBRERÍAS ──────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(sf)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(gridExtra)
library(grid)
library(writexl)
library(spatstat)
library(spatstat.geom)
library(stpp)
library(plotly)
library(GET)
library(fields)
library(RColorBrewer)
library(jsonlite)
library(truncnorm)

setwd("C:/Users/ivanf/Desktop/Octavo/DTIC/Archivos Mabel")

#BASES DE DATOS
data_combinada_catalogo <- read_csv("data_combinada_catalogo (1).csv", col_types = cols(),      
                                    locale = locale(encoding = "UTF-8"))
data_combinada_final <- read_csv("data_combinada_final.csv",col_types = cols(),      
                                 locale = locale(encoding = "UTF-8"))
presencias_agrupadas<- read_csv("presencias_agrupadas_1km.csv",col_types = cols(),      
                                locale = locale(encoding = "UTF-8"))
data_pseudoausencia_final<- read_csv("data_pseudoausencias_final.csv",col_types = cols(),      
                                     locale = locale(encoding = "UTF-8"))
predicciones<- read_csv("predicciones_multiespecie_095.csv",col_types = cols(),      
                        locale = locale(encoding = "UTF-8"))
Ecopal<- read_csv("ECOPAL.csv",col_types = cols(),      
                  locale = locale(encoding = "UTF-8"))
# Mapa oficial base de Ecuador
ecuador <- ne_countries(
  scale = "medium",
  country = "Ecuador",
  returnclass = "sf"
) %>%
  st_transform(4326)

filtrar_dentro_ecuador <- function(data, lon_col, lat_col) {
  
  puntos_sf <- data %>%
    mutate(
      lon_tmp = as.numeric(.data[[lon_col]]),
      lat_tmp = as.numeric(.data[[lat_col]])
    ) %>%
    filter(!is.na(lon_tmp), !is.na(lat_tmp)) %>%
    st_as_sf(
      coords = c("lon_tmp", "lat_tmp"),
      crs = 4326,
      remove = FALSE
    )
  
  dentro <- st_within(puntos_sf, ecuador)
  
  puntos_sf %>%
    mutate(en_ecuador = lengths(dentro) > 0) %>%
    filter(en_ecuador) %>%
    st_drop_geometry() %>%
    select(-lon_tmp, -lat_tmp, -en_ecuador)
}

# Para presencias agrupadas
presencias_agrupadas<- filtrar_dentro_ecuador(presencias_agrupadas,lon_col = "longitud",lat_col = "latitud")

# Para data_combinada_final
data_combinada_final <- filtrar_dentro_ecuador(data_combinada_final,lon_col = "longitud",lat_col = "latitud")

# Para ECOPAL
Ecopal<- filtrar_dentro_ecuador(Ecopal,lon_col = "x",lat_col = "y")

coords_RB <- data_pseudoausencia_final %>%
  mutate(
    geo_parsed = lapply(.geo, fromJSON),
    longitud = sapply(geo_parsed, function(g) g$coordinates[1]),
    latitud  = sapply(geo_parsed, function(g) g$coordinates[2])
  ) %>%
  select(longitud, latitud)

coords_RB <- filtrar_dentro_ecuador(coords_RB, lon_col = "longitud", lat_col = "latitud")

presencias_estandarizada <- presencias_agrupadas %>%
  mutate(
    numero_col = paste0("AGR_", grupo_geografico),   # ID artificial
    ESPECIE = especie,
    AÑO = as.numeric(anio),
    mes = as.numeric(mes),
    `COORD DECIMALES  LATITUD`  = latitud,
    `COORD DECIMALES  LONGITUD` = longitud,
    `ALTITUD msnm` = NA_real_  
  ) %>%
  select(
    numero_col, ESPECIE, AÑO, mes,
    `COORD DECIMALES  LATITUD`,
    `COORD DECIMALES  LONGITUD`,
    `ALTITUD msnm`
  )

base_final <- dplyr::bind_rows(presencias_estandarizada) %>%
  dplyr::mutate(`COORD DECIMALES  LONGITUD` = base::ifelse(`COORD DECIMALES  LONGITUD` > 0, 
                                                           -`COORD DECIMALES  LONGITUD`, 
                                                           `COORD DECIMALES  LONGITUD`)) %>%
  dplyr::filter(`COORD DECIMALES  LATITUD` >= -5 & `COORD DECIMALES  LATITUD` <= 2,
                `COORD DECIMALES  LONGITUD` >= -81 & `COORD DECIMALES  LONGITUD` <= -75)

# Creamos los intervalos para la distribución de los datos por decada y cada 20 años
base_final <- base_final %>%
  mutate(
    # Intervalo cada 20 años   
    periodo_20 = cut(AÑO, 
                     breaks = c(1950, 1969, 1989, 2009, 2019),  # Definir correctamente los cortes de 20 años
                     labels = c("1950-1969", "1970-1989", "1990-2009", "2010-2019"),
                     right = FALSE),  # Para incluir el año inicial y excluir el final
    
    # Intervalo cada 10 años
    periodo_decadas = cut(AÑO, 
                          breaks = c(1915, 1925, 1935, 1945, 1955,1965,1975,1985,1995,2005,2015,2025),  # Cortes para décadas
                          labels = c("1915-1924", "1925-1934", "1935-1944", "1945-1954", "1955-1964",
                                     "1965-1974", "1975-1984", "1985-1994", "1995-2004", "2005-2014", "2015-2025"),
                          right = FALSE)  # Incluir el año inicial, excluir el final
  )

summary(base_final$periodo_20)

summary(base_final$periodo_decadas)

#ESPECIE A ANALIZAR
sort(table(base_final$ESPECIE), decreasing = TRUE)

especie_tratada <- "sylvatica"#"2. anthonyi"#3. "bilinguis"#1. "sylvatica"
Base_especie <- base_final %>% filter(ESPECIE == especie_tratada)
table(Base_especie$periodo_decadas, Base_especie$ESPECIE)
TG_pseudoausencias <- base_final %>% filter(ESPECIE != especie_tratada) %>% filter(AÑO >=min(Base_especie$AÑO))

coords_RB <- data_pseudoausencia_final %>%
  mutate(
    geo_parsed = lapply(.geo, fromJSON),
    longitud = sapply(geo_parsed, function(g) g$coordinates[1]),
    latitud  = sapply(geo_parsed, function(g) g$coordinates[2])
  ) %>%
  select(longitud, latitud) %>%
  filter(
    longitud >= -81.5,
    longitud <= -75,
    latitud  >= -5.5,
    latitud  <= 2
  )

periodos <- levels(Base_especie$periodo_decadas)
set.seed(123)

RB_pseudoausencias <- bind_rows(
  lapply(periodos, function(p) {
    coords_RB[sample(1:nrow(coords_RB), floor(0.05*(table(TG_pseudoausencias$periodo_decadas)[p])/0.95)),] %>%
      mutate(periodo_decadas = p)
  })
)

presencias_df <- Base_especie %>%
  transmute(
    longitud = `COORD DECIMALES  LONGITUD`,
    latitud  = `COORD DECIMALES  LATITUD`,
    periodo_decadas,
    tipo = "presencia",
    presencia = 1,
    Year = AÑO
  )

TG_df <- TG_pseudoausencias %>%
  transmute(
    longitud = `COORD DECIMALES  LONGITUD`,
    latitud  = `COORD DECIMALES  LATITUD`,
    periodo_decadas,
    tipo = "TGB",
    presencia = 0,
    Year = AÑO
  )

RB_df <- RB_pseudoausencias %>%
  transmute(
    longitud,
    latitud,
    periodo_decadas,
    tipo = "RB",
    presencia = 0
  )

Final_df <- bind_rows(
  presencias_df,
  TG_df,
  RB_df
)


tabla_periodo_decadas <- Final_df %>%
  filter(presencia == 1) %>%          # solo presencias
  group_by(periodo_decadas) %>%
  summarise(Conteo = n()) %>%
  mutate(Tipo = "Décadas")

tabla_periodo_decadas
# Filtrar el df completo (con pseudoausencias) desde la década umbral
decada_minima <- tabla_periodo_decadas$periodo_decadas[
  min(which(tabla_periodo_decadas$Conteo >= 15))
]

Final_df <- Final_df %>%
  filter(periodo_decadas >= decada_minima)

table(Final_df$periodo_decadas, Final_df$tipo)

########################################################################################
###################### ZONAS PROTEGIDAS Y VIAS DEL ECUADOR##################################
areas_protegidas <- st_read("C:/Users/ivanf/Desktop/Octavo/DTIC/Archivos Mabel/SNAP_AreasProtegidas_Ecuador.shp")
areas_4326 <- st_transform(areas_protegidas, 4326)

# Ruta de la carpeta donde están todos los archivos del shapefile
ruta_capa <- "C:/Users/ivanf/Desktop/Octavo/DTIC/Archivos Mabel/via_l/via_l.shp"
# Leer la capa
vias <- st_read(ruta_capa)

vias_principales_secundarias <- vias %>%
  filter(
    hct %in% c(14, 15),     # primaria y secundaria
    typ %in% c(1, 17, 41),  # vía, circunvalación, autopista
    rst %in% c(1, 2),       # pavimento o lastrada
    wtc == 1                # todo el año
  )%>%
  st_transform(4326)

##########################################################################################
############## COMPROBACION COORDENADAS VARIABLES AMABIANTALES
Ecopal
datos_filtrados <- Final_df %>%
  filter(!is.na(longitud) & !is.na(latitud)) %>%
  mutate(presencia = ifelse(is.na(presencia), 0, presencia))
table(datos_filtrados$periodo_decadas, datos_filtrados$tipo)


vars_ecopal <- c(
  "bio01", "bio02", "bio03", "bio04", "bio07",
  "bio12", "bio13", "bio14", "bio15", "bio18", "bio19",
  "slope", "aspect", "hillshade", "tri", "watdist", "landcover"
)

# ============================================================
# 1) Match ECOPAL por coordenadas redondeadas a 1 decimal
# ============================================================

datos_modelo <- datos_filtrados %>%
  mutate(
    longitud = as.numeric(longitud),
    latitud  = as.numeric(latitud),
    lon_r = round(longitud, 1),
    lat_r = round(latitud, 1)
  ) %>%
  left_join(
    Ecopal %>%
      mutate(
        lon_r = round(x, 1),
        lat_r = round(y, 1)
      ) %>%
      group_by(lon_r, lat_r) %>%
      summarise(
        across(all_of(vars_ecopal), ~ mean(.x, na.rm = TRUE)),
        n_ecopal_celda = n(),
        .groups = "drop"
      ),
    by = c("lon_r", "lat_r")
  ) %>%
  mutate(
    metodo_ecopal = ifelse(is.na(bio01), NA_character_, "redondeo_1_decimal"),
    dist_ecopal_m = NA_real_
  )

# ============================================================
# 2) Para los que no hicieron match, usar vecino ECOPAL más cercano
# ============================================================

idx_faltantes <- which(is.na(datos_modelo$bio01))

if (length(idx_faltantes) > 0) {
  
  puntos_faltantes_sf <- st_as_sf(
    datos_modelo[idx_faltantes, ],
    coords = c("longitud", "latitud"),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(32717)
  
  ecopal_sf <- st_as_sf(
    Ecopal %>%
      mutate(
        x = as.numeric(x),
        y = as.numeric(y)
      ),
    coords = c("x", "y"),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(32717)
  
  idx_cercano <- st_nearest_feature(puntos_faltantes_sf, ecopal_sf)
  
  dist_m <- st_distance(
    puntos_faltantes_sf,
    ecopal_sf[idx_cercano, ],
    by_element = TRUE
  )
  
  datos_modelo[idx_faltantes, vars_ecopal] <- st_drop_geometry(ecopal_sf)[idx_cercano, vars_ecopal]
  datos_modelo$metodo_ecopal[idx_faltantes] <- "vecino_mas_cercano"
  datos_modelo$dist_ecopal_m[idx_faltantes] <- as.numeric(dist_m)
}


###################### USANDO SIEMPRE  EL VECINO MAS CERCANO
# 1. Convertir datos_filtrados a sf
puntos_sf <- st_as_sf(
  datos_filtrados %>%
    mutate(
      longitud = as.numeric(longitud),
      latitud = as.numeric(latitud)
    ),
  coords = c("longitud", "latitud"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(32717)

# 2. Convertir ECOPAL a sf
ecopal_sf <- st_as_sf(
  Ecopal %>%
    mutate(
      x = as.numeric(x),
      y = as.numeric(y)
    ),
  coords = c("x", "y"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(32717)

# 3. Buscar vecino más cercano
idx_cercano <- st_nearest_feature(puntos_sf, ecopal_sf)

dist_m <- st_distance(
  puntos_sf,
  ecopal_sf[idx_cercano, ],
  by_element = TRUE
)

# 4. Crear base final con ECOPAL asignado por vecino más cercano
datos_modelo_nn <- datos_filtrados %>%
  mutate(
    distancia_ecopal_m = as.numeric(dist_m),
    metodo_ecopal = "vecino_mas_cercano"
  ) %>%
  bind_cols(
    Ecopal[idx_cercano, ] %>%
      select(all_of(vars_ecopal))
  )

# =========================================================
# 1. Preparar la red vial en CRS proyectado para distancias
# =========================================================
vias_dist <- vias_principales_secundarias %>%
  st_transform(32717)

# =========================================================
# 2. Función para agregar distancia a la vía más cercana
# =========================================================
agregar_distancia_via <- function(df, lon = "longitud", lat = "latitud") {
  
  # Convertir observaciones a objeto sf
  puntos_sf <- st_as_sf(
    df,
    coords = c(lon, lat),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(st_crs(vias_dist))
  
  # Índice de la vía más cercana para cada punto
  idx_via_cercana <- st_nearest_feature(puntos_sf, vias_dist)
  
  # Distancia en metros a la vía más cercana
  dist_via_m <- st_distance(
    puntos_sf,
    vias_dist[idx_via_cercana, ],
    by_element = TRUE
  )
  
  # Agregar variables al data frame
  df %>%
    mutate(
      distancia_via_m = as.numeric(dist_via_m),
      distancia_via_km = as.numeric(dist_via_m) / 1000,
      log_distancia_via = log1p(as.numeric(dist_via_m)),
      cerca_via_1km = ifelse(as.numeric(dist_via_m) <= 1000, 1, 0),
      cerca_via_5km = ifelse(as.numeric(dist_via_m) <= 5000, 1, 0)
    )
}

# =========================================================
# 3. Aplicar a las dos bases del modelo
# =========================================================
datos_modelo <- agregar_distancia_via(datos_modelo)

datos_modelo_nn <- agregar_distancia_via(datos_modelo_nn)

datos_modelo %>%
  summarise(
    n = n(),
    min = min(distancia_via_m, na.rm = TRUE),
    q1 = quantile(distancia_via_m, 0.25, na.rm = TRUE),
    mediana = median(distancia_via_m, na.rm = TRUE),
    media = mean(distancia_via_m, na.rm = TRUE),
    q3 = quantile(distancia_via_m, 0.75, na.rm = TRUE),
    max = max(distancia_via_m, na.rm = TRUE),
    hasta_1km = sum(distancia_via_m <= 1000, na.rm = TRUE),
    hasta_5km = sum(distancia_via_m <= 5000, na.rm = TRUE),
    mayor_10km = sum(distancia_via_m > 10000, na.rm = TRUE)
  )

datos_modelo_nn %>%
  summarise(
    n = n(),
    min = min(distancia_via_m, na.rm = TRUE),
    q1 = quantile(distancia_via_m, 0.25, na.rm = TRUE),
    mediana = median(distancia_via_m, na.rm = TRUE),
    media = mean(distancia_via_m, na.rm = TRUE),
    q3 = quantile(distancia_via_m, 0.75, na.rm = TRUE),
    max = max(distancia_via_m, na.rm = TRUE),
    hasta_1km = sum(distancia_via_m <= 1000, na.rm = TRUE),
    hasta_5km = sum(distancia_via_m <= 5000, na.rm = TRUE),
    mayor_10km = sum(distancia_via_m > 10000, na.rm = TRUE)
  )

datos_modelo %>%
  group_by(tipo) %>%
  summarise(
    n = n(),
    mediana = median(distancia_via_m, na.rm = TRUE),
    media = mean(distancia_via_m, na.rm = TRUE),
    max = max(distancia_via_m, na.rm = TRUE),
    cerca_1km = sum(cerca_via_1km == 1, na.rm = TRUE),
    cerca_5km = sum(cerca_via_5km == 1, na.rm = TRUE),
    .groups = "drop"
  )

datos_modelo_nn %>%
  group_by(tipo) %>%
  summarise(
    n = n(),
    mediana = median(distancia_via_m, na.rm = TRUE),
    media = mean(distancia_via_m, na.rm = TRUE),
    max = max(distancia_via_m, na.rm = TRUE),
    cerca_1km = sum(cerca_via_1km == 1, na.rm = TRUE),
    cerca_5km = sum(cerca_via_5km == 1, na.rm = TRUE),
    .groups = "drop"
  )


# =========================================================
# 1. Preparar zonas protegidas para cálculo de distancias
# =========================================================

areas_protegidas_dist <- areas_protegidas %>%
  st_make_valid() %>%
  st_transform(4326) %>%
  # Opcional: recortar visualmente/analíticamente al Ecuador continental
  # Si no quieres recortar, comenta st_crop().
  st_crop(
    xmin = -82,
    xmax = -74.5,
    ymin = -5.5,
    ymax = 2
  ) %>%
  st_transform(32717)

# =========================================================
# 2. Función para agregar distancia a zona protegida más cercana
# =========================================================

agregar_distancia_zp <- function(df, lon = "longitud", lat = "latitud") {
  
  puntos_sf <- st_as_sf(
    df,
    coords = c(lon, lat),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(st_crs(areas_protegidas_dist))
  
  # Zona protegida más cercana
  idx_zp_cercana <- st_nearest_feature(
    puntos_sf,
    areas_protegidas_dist
  )
  
  # Distancia en metros
  dist_zp_m <- st_distance(
    puntos_sf,
    areas_protegidas_dist[idx_zp_cercana, ],
    by_element = TRUE
  )
  
  # Atributos de la zona protegida más cercana
  areas_atrib <- areas_protegidas_dist %>%
    st_drop_geometry()
  
  extraer_col <- function(base, columna, idx) {
    if (columna %in% names(base)) {
      as.character(base[[columna]][idx])
    } else {
      rep(NA_character_, length(idx))
    }
  }
  
  df %>%
    mutate(
      dist_res_m = as.numeric(dist_zp_m),
      dist_res = dist_res_m / 1000,
      log_dist_res = log1p(dist_res_m),
      
      cerca_zp_1km = ifelse(dist_res_m <= 1000, 1, 0),
      cerca_zp_5km = ifelse(dist_res_m <= 5000, 1, 0),
      
      nombre_zp_cercana = extraer_col(areas_atrib, "Nombre", idx_zp_cercana),
      categoria_zp_cercana = extraer_col(areas_atrib, "Categoría", idx_zp_cercana),
      id_zp_cercana = extraer_col(areas_atrib, "ID_Área", idx_zp_cercana)
    )
}

# =========================================================
# 3. Aplicar a las dos bases
# =========================================================

datos_modelo <- agregar_distancia_zp(datos_modelo)

datos_modelo_nn <- agregar_distancia_zp(datos_modelo_nn)

datos_modelo %>%
  summarise(
    n = n(),
    min_dist_zp_m = min(dist_res_m, na.rm = TRUE),
    q1_dist_zp_m = quantile(dist_res_m, 0.25, na.rm = TRUE),
    mediana_dist_zp_m = median(dist_res_m, na.rm = TRUE),
    media_dist_zp_m = mean(dist_res_m, na.rm = TRUE),
    q3_dist_zp_m = quantile(dist_res_m, 0.75, na.rm = TRUE),
    max_dist_zp_m = max(dist_res_m, na.rm = TRUE),
    dentro_o_toca_zp = sum(dist_res_m == 0, na.rm = TRUE),
    hasta_1km = sum(dist_res_m <= 1000, na.rm = TRUE),
    hasta_5km = sum(dist_res_m <= 5000, na.rm = TRUE),
    mayor_10km = sum(dist_res_m > 10000, na.rm = TRUE)
  )

datos_modelo %>%
  group_by(tipo) %>%
  summarise(
    n = n(),
    mediana_dist_zp_km = median(dist_res, na.rm = TRUE),
    media_dist_zp_km = mean(dist_res, na.rm = TRUE),
    dentro_o_toca_zp = sum(dist_res_m == 0, na.rm = TRUE),
    hasta_5km = sum(dist_res_m <= 5000, na.rm = TRUE),
    .groups = "drop"
  )

#-------------------------------Agregar w_sample: Predicción de Danny para TGB, 1 para los demas

# Verificar que exista la especie
if (!especie_tratada %in% names(predicciones)) {
  stop(paste("La especie", especie_tratada, "no existe en predicciones. Revisa names(predicciones)."))
}

pred_sdm <- predicciones %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    rho_sdm = as.numeric(.data[[especie_tratada]])
  ) %>%
  select(x, y, rho_sdm) %>%
  filter(
    !is.na(x),
    !is.na(y),
    !is.na(rho_sdm)
  )

# Revisión rápida
summary(pred_sdm$rho_sdm)
range(pred_sdm$x, na.rm = TRUE)
range(pred_sdm$y, na.rm = TRUE)

# =========================================================
# 2A. Agregar rho_sdm y w_sample a datos_modelo
#     Método: redondeo a 1 decimal + vecino más cercano
# =========================================================

agregar_w_sample_mix <- function(df, pred_sdm, lon = "longitud", lat = "latitud") {
  
  # 1) Match por celda redondeada
  pred_sdm_celda <- pred_sdm %>%
    mutate(
      lon_r_sdm = round(x, 1),
      lat_r_sdm = round(y, 1)
    ) %>%
    group_by(lon_r_sdm, lat_r_sdm) %>%
    summarise(
      rho_sdm = mean(rho_sdm, na.rm = TRUE),
      n_sdm_celda = n(),
      .groups = "drop"
    )
  
  df_out <- df %>%
    mutate(
      lon_r_sdm = round(.data[[lon]], 1),
      lat_r_sdm = round(.data[[lat]], 1)
    ) %>%
    left_join(
      pred_sdm_celda,
      by = c("lon_r_sdm", "lat_r_sdm")
    ) %>%
    mutate(
      metodo_sdm = ifelse(is.na(rho_sdm), NA_character_, "redondeo_1_decimal"),
      dist_sdm_m = NA_real_
    )
  
  # 2) Para faltantes, usar vecino más cercano
  idx_faltantes <- which(is.na(df_out$rho_sdm))
  
  if (length(idx_faltantes) > 0) {
    
    puntos_faltantes_sf <- st_as_sf(
      df_out[idx_faltantes, ],
      coords = c(lon, lat),
      crs = 4326,
      remove = FALSE
    ) %>%
      st_transform(32717)
    
    pred_sdm_sf <- st_as_sf(
      pred_sdm,
      coords = c("x", "y"),
      crs = 4326,
      remove = FALSE
    ) %>%
      st_transform(32717)
    
    idx_cercano <- st_nearest_feature(
      puntos_faltantes_sf,
      pred_sdm_sf
    )
    
    dist_m <- st_distance(
      puntos_faltantes_sf,
      pred_sdm_sf[idx_cercano, ],
      by_element = TRUE
    )
    
    pred_sdm_tabla <- pred_sdm_sf %>%
      st_drop_geometry()
    
    df_out$rho_sdm[idx_faltantes] <- pred_sdm_tabla$rho_sdm[idx_cercano]
    df_out$metodo_sdm[idx_faltantes] <- "vecino_mas_cercano"
    df_out$dist_sdm_m[idx_faltantes] <- as.numeric(dist_m)
  }
  
  # 3) Crear peso muestral
  df_out <- df_out %>%
    mutate(
      w_sample = case_when(
        tipo == "TGB" ~ rho_sdm,
        tipo == "RB" ~ 1,
        tipo == "Presencia" ~ 1,
        TRUE ~ 1
      )
    )
  
  return(df_out)
}

datos_modelo <- agregar_w_sample_mix(
  df = datos_modelo,
  pred_sdm = pred_sdm
)

# =========================================================
# 2B. Agregar rho_sdm y w_sample a datos_modelo_nn
#     Método: siempre vecino más cercano
# =========================================================

agregar_w_sample_nn <- function(df, pred_sdm, lon = "longitud", lat = "latitud") {
  
  puntos_sf <- st_as_sf(
    df,
    coords = c(lon, lat),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(32717)
  
  pred_sdm_sf <- st_as_sf(
    pred_sdm,
    coords = c("x", "y"),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(32717)
  
  idx_cercano <- st_nearest_feature(
    puntos_sf,
    pred_sdm_sf
  )
  
  dist_m <- st_distance(
    puntos_sf,
    pred_sdm_sf[idx_cercano, ],
    by_element = TRUE
  )
  
  pred_sdm_tabla <- pred_sdm_sf %>%
    st_drop_geometry()
  
  df %>%
    mutate(
      rho_sdm = pred_sdm_tabla$rho_sdm[idx_cercano],
      metodo_sdm = "vecino_mas_cercano",
      dist_sdm_m = as.numeric(dist_m),
      
      w_sample = case_when(
        tipo == "TGB" ~ rho_sdm,
        tipo == "RB" ~ 1,
        tipo == "Presencia" ~ 1,
        TRUE ~ 1
      )
    )
}

datos_modelo_nn <- agregar_w_sample_nn(
  df = datos_modelo_nn,
  pred_sdm = pred_sdm
)


datos_modelo %>%
  summarise(
    n = n(),
    sin_rho_sdm = sum(is.na(rho_sdm)),
    sin_w_sample = sum(is.na(w_sample)),
    min_w = min(w_sample, na.rm = TRUE),
    mediana_w = median(w_sample, na.rm = TRUE),
    media_w = mean(w_sample, na.rm = TRUE),
    max_w = max(w_sample, na.rm = TRUE)
  )

datos_modelo_nn %>%
  summarise(
    n = n(),
    sin_rho_sdm = sum(is.na(rho_sdm)),
    sin_w_sample = sum(is.na(w_sample)),
    min_w = min(w_sample, na.rm = TRUE),
    mediana_w = median(w_sample, na.rm = TRUE),
    media_w = mean(w_sample, na.rm = TRUE),
    max_w = max(w_sample, na.rm = TRUE)
  )


datos_modelo %>%
  group_by(tipo) %>%
  summarise(
    n = n(),
    min_w = min(w_sample, na.rm = TRUE),
    q1_w = quantile(w_sample, 0.25, na.rm = TRUE),
    mediana_w = median(w_sample, na.rm = TRUE),
    media_w = mean(w_sample, na.rm = TRUE),
    q3_w = quantile(w_sample, 0.75, na.rm = TRUE),
    max_w = max(w_sample, na.rm = TRUE),
    .groups = "drop"
  )

datos_modelo_nn %>%
  group_by(tipo) %>%
  summarise(
    n = n(),
    min_w = min(w_sample, na.rm = TRUE),
    q1_w = quantile(w_sample, 0.25, na.rm = TRUE),
    mediana_w = median(w_sample, na.rm = TRUE),
    media_w = mean(w_sample, na.rm = TRUE),
    q3_w = quantile(w_sample, 0.75, na.rm = TRUE),
    max_w = max(w_sample, na.rm = TRUE),
    .groups = "drop"
  )



datos_modelo %>%
  count(tipo, metodo_sdm)

datos_modelo_nn %>%
  count(tipo, metodo_sdm)


datos_modelo %>%
  filter(metodo_sdm == "vecino_mas_cercano") %>%
  summarise(
    n = n(),
    min = min(dist_sdm_m, na.rm = TRUE),
    mediana = median(dist_sdm_m, na.rm = TRUE),
    media = mean(dist_sdm_m, na.rm = TRUE),
    max = max(dist_sdm_m, na.rm = TRUE),
    mayor_5km = sum(dist_sdm_m > 5000, na.rm = TRUE),
    mayor_10km = sum(dist_sdm_m > 10000, na.rm = TRUE)
  )

datos_modelo_nn %>%
  summarise(
    n = n(),
    min = min(dist_sdm_m, na.rm = TRUE),
    mediana = median(dist_sdm_m, na.rm = TRUE),
    media = mean(dist_sdm_m, na.rm = TRUE),
    max = max(dist_sdm_m, na.rm = TRUE),
    mayor_5km = sum(dist_sdm_m > 5000, na.rm = TRUE),
    mayor_10km = sum(dist_sdm_m > 10000, na.rm = TRUE)
  )

#########################################################################################
#---------------Preparar datos para los siguientes cripts--------------------------------
#########################################################################################

asignar_year_rb <- function(base_datos, seed = 2025) {
  
  set.seed(seed)
  
  base_datos <- as.data.frame(base_datos)
  
  # Normalizar nombre Year/year
  if ("year" %in% names(base_datos) && !"Year" %in% names(base_datos)) {
    base_datos$Year <- base_datos$year
  }
  if ("Year" %in% names(base_datos) && !"year" %in% names(base_datos)) {
    base_datos$year <- base_datos$Year
  }
  
  # Identificar RB sin año
  idx_rb_na <- which(
    base_datos$tipo %in% c("RB", "Pseudoausencia_RB") &
      is.na(base_datos$Year)
  )
  
  if (length(idx_rb_na) == 0) {
    message("No hay RB con Year faltante.")
    return(base_datos)
  }
  
  # Caso 1: si existe periodo_decadas tipo "1985-1994"
  if ("periodo_decadas" %in% names(base_datos)) {
    
    periodo <- as.character(base_datos$periodo_decadas[idx_rb_na])
    
    año_ini <- as.numeric(sub("^(\\d{4}).*$", "\\1", periodo))
    año_fin <- as.numeric(sub("^\\d{4}\\D+(\\d{4}).*$", "\\1", periodo))
    
    # Si no logra leer el año final, usar inicio + 9
    año_fin[is.na(año_fin)] <- año_ini[is.na(año_fin)] + 9
    
    # Asignar año aleatorio dentro de la década
    base_datos$Year[idx_rb_na] <- mapply(
      function(a, b) {
        if (is.na(a) || is.na(b)) return(NA_real_)
        sample(seq(a, b), size = 1)
      },
      año_ini,
      año_fin
    )
    
  } else {
    
    # Caso 2: si no existe periodo_decadas, muestrear años observados
    años_obs <- base_datos$Year[!is.na(base_datos$Year)]
    
    base_datos$Year[idx_rb_na] <- sample(
      años_obs,
      size = length(idx_rb_na),
      replace = TRUE
    )
  }
  
  # Mantener year sincronizado
  base_datos$year <- base_datos$Year
  
  return(base_datos)
}


datos_modelo <- asignar_year_rb(datos_modelo)
datos_modelo_nn <- asignar_year_rb(datos_modelo_nn)


# 1. Para sylvatica
# save(list = c("datos_modelo"), file = "datos_modelo_sylvatica.RData", envir = .GlobalEnv)
# save(list = c("datos_modelo_nn"), file = "datos_modelo_nn_sylvatica.RData", envir = .GlobalEnv)

# 2. Para anthonyi
# save(list = c("datos_modelo"), file = "datos_modelo_anthonyi.RData", envir = .GlobalEnv)
# save(list = c("datos_modelo_nn"), file = "datos_modelo_nn_anthonyi.RData", envir = .GlobalEnv)

# 3. Para bilinguis
# save(list = c("datos_modelo"), file = "datos_modelo_bilinguis.RData", envir = .GlobalEnv)
# save(list = c("datos_modelo_nn"), file = "datos_modelo_nn_bilinguis.RData", envir = .GlobalEnv)




#----------------------------------------GRAFICOS ---------------------------------------

#grafico general
ggplot() +
  geom_sf(data = ecuador, fill = "gray95", color = "gray60") +
  geom_point(
    data = datos_filtrados,
    aes(x = longitud, y = latitud, color = tipo),
    size = 2,
    alpha = 0.75
  ) +
  coord_sf(
    xlim = c(-92, -75),
    ylim = c(-5.5, 2.5),
    expand = FALSE
  ) +
  labs(
    title = "Datos filtrados en Ecuador",
    subtitle = "Presencias y pseudoausencias",
    x = "Longitud",
    y = "Latitud",
    color = "Tipo"
  ) +
  theme_minimal()

nombre_leyenda <- paste0("Presencia de la especie ", especie_tratada)
grafico <- ggplot() +
  # ── Fondo: Ecuador ────────────────────────────────────────────────────────
  geom_sf(data = ecuador,
          fill = "#F0EDE4", color = "#AAAAAA", linewidth = 0.4) +
  
  # ── Áreas protegidas ──────────────────────────────────────────────────────
  geom_sf(data = areas_4326,
          aes(fill = "Áreas protegidas"),
          color = "#1A6B2E", alpha = 0.55, linewidth = 0.2) +
  
  # ── Vías principales/secundarias ──────────────────────────────────────────
  geom_sf(data = vias_principales_secundarias,
          aes(color = "Vías principales"),
          linewidth = 0.45, alpha = 0.75) +
  
  # ── Presencias de la especie ──────────────────────────────────────────────
  geom_point(
    data = datos_filtrados %>% filter(presencia == 1),
    #data = datos_modelo %>% filter(distancia_via_m >=10000),
    #data = datos_modelo %>% filter(dist_res >= 10),
    aes(x = longitud, y = latitud, shape = nombre_leyenda),
    color = "#E07B00", size = 1.4, alpha = 0.85
  ) +
  
  # ── Escalas manuales para la leyenda ─────────────────────────────────────
  scale_fill_manual(
    name   = NULL,
    values = c("Áreas protegidas" = "#2E8B4A")
  ) +
  scale_color_manual(
    name   = NULL,
    values = c("Vías principales" = "#2255AA")
  ) +
  scale_shape_manual(
    name   = NULL,
    values =setNames(16, nombre_leyenda)
  ) +
  
  # ── Guías de leyenda unificadas ───────────────────────────────────────────
  guides(
    fill  = guide_legend(override.aes = list(alpha = 0.6, size = 5)),
    color = guide_legend(override.aes = list(linewidth = 1.2)),
    shape = guide_legend(override.aes = list(color = "#E07B00", size = 3))
  ) +
  
  # ── Extensión espacial ────────────────────────────────────────────────────
  coord_sf(xlim = c(-81.5, -74.5), ylim = c(-5.5, 2), expand = FALSE) +
  
  # ── Etiquetas ─────────────────────────────────────────────────────────────
  labs(
    title    = "Carreteras y Zonas Protegidas del Ecuador",
    subtitle = paste0("Distribución de presencias de la especie ",especie_tratada," respecto a la infraestructura vial y áreas de conservación"),
    x = "Longitud", y = "Latitud"
  ) +
  
  # ── Tema ──────────────────────────────────────────────────────────────────
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 9, color = "gray35", hjust = 0,
                                    margin = margin(b = 8)),
    plot.caption     = element_text(size = 7.5, color = "gray50", hjust = 1),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_line(color = "#DCDCDC", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#D6EAF8", color = NA),  # fondo mar
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(10, 10, 10, 10)
  )

#ggsave("mapa_ecuador.png", grafico, width = 14, height = 10, dpi = 300)
ggsave(paste0("mapa_ecuador_",especie_tratada,".png"), grafico,
       width = 10, height = 10,  # más cuadrado, igual al aspecto del mapa
       dpi = 300)
grafico
