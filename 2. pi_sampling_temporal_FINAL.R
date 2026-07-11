# CODIGO PI SAMPLING

#Cargar datos 

# 1. sylvatica
# especie_tratada <- "sylvatica"
# load("datos_modelo_sylvatica.RData")
# load("datos_modelo_nn_sylvatica.RData")

# 2. anthonyi
# especie_tratada <- "anthonyi"
# load("datos_modelo_anthonyi.RData")    #Base hecha con redondeo + Vecino cercano(NN)
# load("datos_modelo_nn_anthonyi.RData") #Base solo con vecino cercano (NN)

# 3. bilinguis
especie_tratada <- "bilinguis"
load("datos_modelo_bilinguis.RData")
load("datos_modelo_nn_bilinguis.RData")

# INSTALACIÓN Y CARGA DE PAQUETES
pkgs_req <- c("readxl", "dplyr", "lubridate", "sf", "spatstat", "spatstat.geom", "stpp", "GET", "fields", "jsonlite", "spdep")
pkgs_opt <- c("ggplot2", "viridis", "rnaturalearth", "rnaturalearthdata", "gridExtra", "writexl", "plotly", "RColorBrewer")

for (p in pkgs_req) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

use_ggplot <- all(sapply(pkgs_opt, requireNamespace, quietly = TRUE))
if (use_ggplot) {
  for (p in pkgs_opt)
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  cat("  Gráficos ggplot2 disponibles.\n")
} else {
  cat("  ggplot2 no disponible — usando gráficos base como fallback.\n")
}

set.seed(2025)

ecuador <- ne_countries(
  scale = "medium",
  country = "Ecuador",
  returnclass = "sf"
)

# FUNCIÓN PARA PREPARAR CADA BASE PARA EL PIPELINE
preparar_df_pi <- function(base_datos, quitar_na = TRUE) {
  
  # Normalizar nombre de columna año: acepta 'year' o 'Year'
  base_datos <- as.data.frame(base_datos)
  if ("year" %in% names(base_datos) && !"Year" %in% names(base_datos))
    base_datos$Year <- base_datos$year
  if ("Year" %in% names(base_datos) && !"year" %in% names(base_datos))
    base_datos$year <- base_datos$Year
  
  df <- base_datos %>%
    as.data.frame() %>%
    dplyr::mutate(
      # Coordenadas
      lon = as.numeric(longitud),
      lat = as.numeric(latitud),
      
      # Variable respuesta del modelo de muestreo
      M = as.numeric(presencia),
      
      # Año normalizado como covariable de detección temporal
      # year_norm = 0 en el año mínimo del dataset, 1 en el máximo
      # Captura el sesgo MNAR temporal: ciencia ciudadana post-2014
      # aumenta la probabilidad de detección en años recientes
      year_min  = min(as.numeric(Year), na.rm = TRUE),
      year_max  = max(as.numeric(Year), na.rm = TRUE),
      year_norm = (as.numeric(Year) - year_min) / 
        pmax(year_max - year_min, 1),
      
      # Año
      Year = as.numeric(Year),
      
      # Variables de accesibilidad
      dist_road = as.numeric(distancia_via_km),
      dist_res = as.numeric(dist_res),
      
      # Peso muestral para estimar pi
      # Presencias y RB deben pesar 1; TGB usa rho_sdm.
      w_sample = as.numeric(w_sample),
      
      # Landcover como factor, no como variable continua
      landcover = as.factor(landcover),
      
      # Aspect es circular: no conviene meterlo crudo como variable lineal
      aspect_rad = as.numeric(aspect) * pi / 180,
      aspect_sin = sin(aspect_rad),
      aspect_cos = cos(aspect_rad)
    )
  
  # Variables numéricas que sí se pueden estandarizar
  vars_z <- c(
    # accesibilidad
    "dist_road", "dist_res",
    
    # bioclimáticas
    "bio01", "bio02", "bio03", "bio04", "bio07",
    "bio12", "bio13", "bio14", "bio15", "bio18", "bio19",
    
    # topográficas
    "slope", "hillshade", "tri", "watdist",
    
    # aspect transformado
    "aspect_sin", "aspect_cos",
    
    # temporal — year_norm ya está en [0,1], no necesita estandarización adicional
    # pero se incluye en vars_z para que tenga su versión _z también
    "year_norm"
  )
  
  vars_z <- vars_z[vars_z %in% names(df)]
  
  for (v in vars_z) {
    df[[paste0(v, "_z")]] <- as.numeric(scale(as.numeric(df[[v]])))
  }
  
  # Selección final para los scripts
  columnas_finales <- c(
    # Identificación básica
    "lon", "lat", "M", "Year", "tipo",
    
    # Temporal normalizado — covariable de detección
    "year_norm",
    
    # Pesos de pseudoausencias para estimar pi
    "w_sample",
    
    # Variables de accesibilidad crudas
    "dist_road", "dist_res",
    
    # Variables ambientales crudas
    "bio01", "bio02", "bio03", "bio04", "bio07",
    "bio12", "bio13", "bio14", "bio15", "bio18", "bio19",
    "slope", "aspect", "aspect_sin", "aspect_cos",
    "hillshade", "tri", "watdist", "landcover",
    
    # Variables estandarizadas
    paste0(vars_z, "_z"),
    
    # Variables útiles de control / trazabilidad si existen
    "rho_sdm", "metodo_sdm", "dist_sdm_m",
    "metodo_ecopal", "dist_ecopal_m",
    "dist_res_m", "distancia_via_m", "distancia_via_km"
  )
  
  columnas_finales <- columnas_finales[columnas_finales %in% names(df)]
  
  df <- df %>%
    dplyr::select(dplyr::all_of(columnas_finales))
  
  if (quitar_na) {
    df <- df %>%
      dplyr::filter(
        stats::complete.cases(
          lon, lat, M,
          Year,          # necesario para year_norm — evita desajuste de dimensiones
          dist_road, dist_res,
          w_sample,
          slope, tri, hillshade, watdist,
          aspect_sin, aspect_cos,
          landcover
        )
      )
  }
  
  df <- df %>% as.data.frame()
  
  return(df)
}

# Funcion para SAR-Logit
sar_logit_pl <- function(df, listw, formula_base,
                         lambda_max = 1.0,
                         tol = 1e-5, max_iter = 200,
                         verbose = TRUE) {
  # Inicializar con logístico estándar
  pi_hat <- predict(glm(formula_base, data = df, family = binomial()),
                    type = "response")
  
  rho_a     <- 0
  converged <- FALSE
  hist_rho  <- hist_diff <- hist_L <- numeric(max_iter)
  
  for (iter in seq_len(max_iter)) {
    
    # Paso A: lag espacial determinístico (SIN ε)
    df$Wpi <- as.numeric(lag.listw(listw, pi_hat))
    
    # Paso B: logístico con lag como covariable adicional
    # Pesos de muestra: presencias=1, RB=1, TGB=rho_SDM
    # Esto evita que TGB en zonas de alta probabilidad ecologica contaminen pi
    mod_i  <- glm(update(formula_base, . ~ . + Wpi),
                  data = df, family = binomial(), weights =w_sample)
    
    pi_new <- predict(mod_i, type = "response")
    rho_a  <- coef(mod_i)["Wpi"]
    L_i    <- abs(rho_a) * lambda_max / 4
    
    # Paso C: norma sup
    diff_i           <- max(abs(pi_new - pi_hat))
    hist_rho[iter]   <- rho_a
    hist_diff[iter]  <- diff_i
    hist_L[iter]     <- L_i
    
    # Advertir si la condición de contracción se viola
    if (L_i >= 1 && iter <= 3)
      cat(sprintf("  AVISO iter %d: L = %.4f >= 1 (rho_a = %.4f)\n",
                  iter, L_i, rho_a))
    
    if (verbose && (iter == 1 || iter %% 25 == 0))
      cat(sprintf("  iter %3d | rho_a = %+.5f | Delta_pi = %.2e | L = %.4f\n",
                  iter, rho_a, diff_i, L_i))
    
    # Actualizar
    pi_hat <- pi_new
    
    if (diff_i < tol) {
      converged <- TRUE
      if (verbose)
        cat(sprintf("  Convergido: iter = %d  | Delta_pi = %.2e\n",
                    iter, diff_i))
      break
    }
  }
  
  if (!converged) warning("SAR-logit no convergió en ", max_iter, " iter.")
  
  list(model    = mod_i,
       pi_hat   = pi_hat,
       rho_a    = rho_a,
       L        = abs(rho_a) * lambda_max / 4,
       n_iter   = iter,
       converged= converged,
       hist_rho = hist_rho[seq_len(iter)],
       hist_diff= hist_diff[seq_len(iter)],
       hist_L   = hist_L[seq_len(iter)])
}

#Funcion del Pipeline completo pi_samplinng
correr_pipeline_pi <- function(df,
                               formula_std,
                               nombre_modelo = "modelo_pi",
                               k = 5,
                               umbral = 0.05,
                               generar_graficos = TRUE,
                               exportar = FALSE,
                               carpeta_salida = ".",
                               ecuador = NULL) {
  
  DIR_OUT     <- "."
  # =============================================================================
  # 1.  ESTRUCTURA DE DATOS REALES
  # =============================================================================
  
  cat("\n=== SECCIÓN 1: Datos reales =============================================\n")
  
  df <- as.data.frame(df)
  
  n        <- nrow(df)
  n_pres   <- sum(df$M == 1, na.rm = TRUE)
  n_pseudo <- sum(df$M == 0, na.rm = TRUE)
  n_RB     <- sum(df$tipo %in% c("RB", "Pseudoausencia_RB"), na.rm = TRUE)
  n_TGB    <- sum(df$tipo %in% c("TGB", "Pseudoausencia_TG", "Pseudoausencia_TGB"), na.rm = TRUE)
  
  cat(sprintf("  Registros totales : %d\n", n))
  cat(sprintf("  Presencias  (M=1) : %d\n", n_pres))
  cat(sprintf("  Pseudoaus. RB     : %d  (M=0, q proporcional 1)\n", n_RB))
  cat(sprintf("  Pseudoaus. TGB    : %d  (M=0, q proporcional p_sp)\n", n_TGB))
  cat(sprintf("  Proporcion TGB    : %.0f%%\n", 100 * n_TGB / max(n_pseudo, 1)))
  cat(sprintf("  Prevalencia obs.  : %.3f\n", mean(df$M)))
  
  # =============================================================================
  # 2.  MATRIZ DE PESOS ESPACIALES W  (k-NN k=5, row-estandarizada)
  # =============================================================================
  
  cat("\n=== SECCIÓN 2: Matriz de pesos W ========================================\n")
  
  coords_mat      <- cbind(df$lon, df$lat)
  nb_knn5         <- knn2nb(knearneigh(coords_mat, k = k))
  lw              <- nb2listw(nb_knn5, style = "W")
  lambda_max_W    <- 1.0
  contraccion_lim <- 4.0
  
  cat(sprintf("  k-NN = %d  |  λ_max(W) ≤ %.1f  |  condición contracción: |ρ_a| < %.1f\n",
              k, lambda_max_W, contraccion_lim))
  
  # =============================================================================
  # 2b. PESOS DE MUESTRA DEL SDM PARA LAS PSEUDOAUSENCIAS TGB
  # =============================================================================
  
  cat("\n=== SECCIÓN 2b: Pesos de muestra SDM para TGB ==========================\n")
  
  df$w_sample <- as.numeric(df$w_sample)
  
  if (any(is.na(df$w_sample))) {
    stop("Existen NA en w_sample. Revisa la construcción de pesos para TGB/RB/presencias.")
  }
  
  cat(sprintf("  Peso medio presencias  (w_j=1)      : %.3f\n",
              mean(df$w_sample[df$M == 1], na.rm = TRUE)))
  cat(sprintf("  Peso medio RB          (w_j=1)      : %.3f\n",
              mean(df$w_sample[df$tipo %in% c("RB", "Pseudoausencia_RB")], na.rm = TRUE)))
  cat(sprintf("  Peso medio TGB         (w_j=rho_SDM): %.3f\n",
              mean(df$w_sample[df$tipo %in% c("TGB", "Pseudoausencia_TG", "Pseudoausencia_TGB")], na.rm = TRUE)))
  cat("  TGB con w_sample < 0.1 (ruido minimo):",
      sum(df$w_sample[df$tipo %in% c("TGB", "Pseudoausencia_TG", "Pseudoausencia_TGB")] < 0.1, na.rm = TRUE), "\n")
  cat("  TGB con w_sample > 0.5 (ruido alto, peso reducido):",
      sum(df$w_sample[df$tipo %in% c("TGB", "Pseudoausencia_TG", "Pseudoausencia_TGB")] > 0.5, na.rm = TRUE), "\n")
  
  # =============================================================================
  # 3.  ESCALÓN 1: LOGÍSTICO ESTÁNDAR
  # =============================================================================
  
  cat("\n=== ESCALÓN 1: Logístico estándar =======================================\n")
  
  mod_logit        <- glm(formula_std, data = df, family = binomial())
  df$pi_logit      <- predict(mod_logit, type = "response")
  df$resid_pearson <- residuals(mod_logit, type = "pearson")
  
  aic_logit  <- AIC(mod_logit)
  pR2_logit  <- 1 - mod_logit$deviance / mod_logit$null.deviance
  
  cat("\n  Coeficientes:\n")
  print(round(summary(mod_logit)$coefficients, 4))
  cat(sprintf("\n  AIC = %.2f  |  Pseudo-R² (McFadden) = %.4f\n",
              aic_logit, pR2_logit))
  cat(sprintf("  π̂ rango: [%.4f, %.4f]  media = %.4f\n",
              min(df$pi_logit), max(df$pi_logit), mean(df$pi_logit)))
  
  # =============================================================================
  # 4.  ESCALÓN 1 — DIAGNÓSTICO: TEST DE MORAN I
  # =============================================================================
  
  cat("\n=== DIAGNÓSTICO 1: Test de Moran I (residuos logístico) =================\n")
  
  moran_logit <- moran.test(df$resid_pearson, listw = lw, alternative = "greater")
  MI_logit    <- moran_logit$estimate["Moran I statistic"]
  pval_logit  <- moran_logit$p.value
  
  cat(sprintf("  Moran I = %.4f  |  E[I] = %.4f  |  p-valor = %.4e\n",
              MI_logit,
              moran_logit$estimate["Expectation"],
              pval_logit))
  
  UMBRAL         <- umbral
  spatial_needed <- pval_logit < UMBRAL
  
  if (spatial_needed) {
    cat("  RESULTADO: autocorrelacion espacial significativa (p < 0.05)\n")
    cat("  -> Subiendo al ESCALON 2: SAR-logit\n")
  } else {
    cat("  RESULTADO: sin autocorrelacion significativa (p >= 0.05)\n")
    cat("  -> Logistico estandar es suficiente para pi\n")
  }
  
  # =============================================================================
  # 5.  ESCALÓN 2: SAR-LOGIT DETERMINÍSTICO SIN TÉRMINO DE ERROR
  # =============================================================================
  
  cat("\n=== ESCALÓN 2: SAR-logit (PL iterativo) =================================\n")
  
  sar <- sar_logit_pl(df = df, lw, formula_std,
                      lambda_max = lambda_max_W,
                      tol = 1e-5, max_iter = 200, verbose = TRUE)
  
  df$pi_sar    <- sar$pi_hat
  df$resid_sar <- residuals(sar$model, type = "pearson")
  aic_sar      <- AIC(sar$model)
  
  cat(sprintf("\n  rho_a = %+.5f  |  L = %.5f  (%s)\n",
              sar$rho_a, sar$L,
              ifelse(sar$L < 1, "contraccion OK", "CONTRACCION VIOLADA")))
  cat(sprintf("  Iteraciones = %d  |  Convergencia = %s\n",
              sar$n_iter, ifelse(sar$converged, "SI", "NO")))
  cat(sprintf("  AIC: logístico = %.2f  |  SAR = %.2f  |  delta = %+.2f\n",
              aic_logit, aic_sar, aic_sar - aic_logit))
  
  cat("\n  Coeficientes SAR-logit (última iteración):\n")
  print(round(summary(sar$model)$coefficients, 4))
  
  # =============================================================================
  # 6.  DIAGNÓSTICO 2: MORAN I POST SAR-LOGIT
  # =============================================================================
  
  cat("\n=== DIAGNÓSTICO 2: Moran I (residuos SAR-logit) =========================\n")
  
  moran_sar  <- moran.test(df$resid_sar, listw = lw, alternative = "greater")
  MI_sar     <- moran_sar$estimate["Moran I statistic"]
  pval_sar   <- moran_sar$p.value
  reduccion  <- 100 * (1 - abs(MI_sar) / abs(MI_logit))
  
  cat(sprintf("  Moran I = %.4f  |  p-valor = %.4e\n", MI_sar, pval_sar))
  cat(sprintf("  Reduccion de autocorrelacion: %.1f%%  (%.4f -> %.4f)\n",
              reduccion, MI_logit, MI_sar))
  
  if (pval_sar < UMBRAL) {
    cat("  RESULTADO: autocorrelacion residual persiste (p < 0.05)\n")
    cat("  -> Considerar LGCP con R-INLA (ver Sección 12)\n")
  } else {
    cat("  RESULTADO: autocorrelacion residual eliminada (p >= 0.05)\n")
    cat("  -> SAR-logit es suficiente\n")
  }
  
  # =============================================================================
  # 7.  ESCALÓN 3: IPW ESTABILIZADO
  # =============================================================================
  
  cat("\n=== ESCALÓN 3: IPW estabilizado =========================================\n")
  
  pi_final     <- if (spatial_needed && sar$converged) df$pi_sar else df$pi_logit
  modelo_final <- if (spatial_needed && sar$converged) "SAR-logit" else "Logístico"
  df$pi_final  <- pi_final
  cat(sprintf("  Modelo pi final: %s\n\n", modelo_final))
  
  idx_pres   <- which(df$M == 1)
  pi_pres    <- df$pi_final[idx_pres]
  n_pres_eff <- length(idx_pres)
  
  w_raw  <- 1 / pi_pres
  w99    <- quantile(w_raw, 0.99)
  w_tr   <- pmin(w_raw, w99)
  w_stab <- w_tr * n_pres_eff / sum(w_tr)
  
  df$w_ipw           <- NA
  df$w_ipw[idx_pres] <- w_stab
  n_trunc            <- sum(w_raw > w99)
  
  # Effective Sample Size
  ESS     <- sum(w_stab)^2 / sum(w_stab^2)
  ESS_pct <- 100 * ESS / n_pres_eff
  
  cat(sprintf("  Pesos brutos    : min=%.3f  mediana=%.3f  max=%.3f\n",
              min(w_raw), median(w_raw), max(w_raw)))
  cat(sprintf("  w99 (umbral)    : %.4f\n", w99))
  cat(sprintf("  Truncados       : %d  (%.1f%%)\n",
              n_trunc, 100 * n_trunc / n_pres_eff))
  cat(sprintf("  Pesos estables  : min=%.4f  mediana=%.4f  max=%.4f  suma=%.0f\n",
              min(w_stab), median(w_stab), max(w_stab), sum(w_stab)))
  cat(sprintf("  ESS             : %.1f  (%.1f%%)\n", ESS, ESS_pct))
  
  # =============================================================================
  # 8.  CALIBRACIÓN Y VALIDACIÓN DE π̂
  # =============================================================================
  
  cat("\n=== SECCIÓN 8: Calibración de pi ========================================\n")
  
  q_breaks   <- quantile(df$pi_final, probs = seq(0, 1, 0.2))
  df$qpi     <- cut(df$pi_final, breaks = q_breaks,
                    include.lowest = TRUE, labels = paste0("Q", 1:5))
  cal        <- aggregate(M ~ qpi, data = df, FUN = mean)
  cal$pi_med <- tapply(df$pi_final, df$qpi, mean)
  cal$delta  <- abs(cal$M - cal$pi_med)
  
  cat("  Calibración por quintiles de pi:\n")
  cat(sprintf("  %-6s  %-10s  %-10s  %-8s  %s\n",
              "Quint.", "pi medio", "Obs M=1", "|delta|", "OK"))
  for (i in seq_len(nrow(cal))) {
    flag <- ifelse(cal$delta[i] < 0.10, "OK", "revisar")
    cat(sprintf("  %-6s  %-10.4f  %-10.4f  %-8.4f  %s\n",
                as.character(cal$qpi[i]), cal$pi_med[i], cal$M[i],
                cal$delta[i], flag))
  }
  cat(sprintf("\n  ESS = %.1f (%.0f%%)  %s\n", ESS, ESS_pct,
              ifelse(ESS_pct > 50, "-> ESS aceptable", "-> ESS bajo, revisar")))
  
  # =============================================================================
  # 9.  GRÁFICOS
  # =============================================================================
  
  cat("\n=== SECCIÓN 9: Gráficos =================================================\n")
  
  plots <- NULL
  
  if (generar_graficos) {
    
    if (use_ggplot) {
      
      if (!is.null(ecuador)) {
        p1 <- ggplot() +
          geom_sf(
            data = ecuador,
            fill = "gray95",
            color = "gray50",
            linewidth = 0.3
          ) +
          geom_point(
            data = df,
            aes(x = lon, y = lat, colour = pi_final),
            size = 1.4,
            alpha = 0.75
          ) +
          scale_colour_viridis_c(
            option = "plasma",
            direction = -1,
            name = expression(hat(pi))
          ) +
          coord_sf(
            xlim = c(-81.5, -74.5),
            ylim = c(-5.5, 2),
            expand = FALSE
          ) +
          labs(
            title = TeX(paste("$\\pi$ final —", modelo_final)),
            subtitle = TeX(sprintf("$\\rho_a$ = %.4f  |  L = %.4f", sar$rho_a, sar$L)),
            x = "Longitud",
            y = "Latitud"
          ) +
          theme_bw(base_size = 11)
      } else {
        p1 <- ggplot(df, aes(lon, lat, colour = pi_final)) +
          geom_point(size = 1.4, alpha = 0.75) +
          scale_colour_viridis_c(option = "plasma", direction = -1,
                                 name = expression(hat(pi))) +
          labs(title    = TeX(paste("$\\pi$ final —", modelo_final)),
               subtitle = TeX(sprintf("$\\rho_a$ = %.4f  |  L = %.4f", sar$rho_a, sar$L)),
               x = "Longitud", y = "Latitud") +
          theme_bw(base_size = 11)
      }
      
      p2 <- ggplot(df, aes(pi_logit, pi_sar, colour = factor(M))) +
        geom_point(alpha = 0.3, size = 1.2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        scale_colour_manual(values = c("0" = "#999999", "1" = "#006400"),
                            labels = c("Pseudoausencia", "Presencia")) +
        labs(title = "Logístico vs. SAR-logit",
             x = expression(hat(pi)[logit]),
             y = expression(hat(pi)[SAR]), colour = NULL) +
        theme_bw(base_size = 11) + theme(legend.position = "bottom")
      
      p3 <- ggplot(data.frame(w = w_stab), aes(x = w)) +
        geom_histogram(bins = 40, fill = "#006400", alpha = 0.7, colour = "white") +
        geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
        labs(title    = "Pesos IPW estabilizados",
             subtitle = TeX(sprintf("ESS = %.0f (%.0f%%)  |  $w_{99}$ = %.2f",
                                ESS, ESS_pct, w99)),
             x = expression(tilde(w)[i]), y = "Conteo") +
        theme_bw(base_size = 11)
      
      conv_df <- data.frame(iter = seq_along(sar$hist_diff),
                            diff = sar$hist_diff)
      p4 <- ggplot(conv_df, aes(iter, diff)) +
        geom_line(colour = "#006400", linewidth = 0.9) +
        geom_hline(yintercept = 1e-5, linetype = "dashed", colour = "red") +
        scale_y_log10() +
        labs(title = TeX("Convergencia SAR-logit  $(||pi_{nuevo} - pi_{viejo}||_{\\inf})$"),
             x = "Iteración", y = TeX("Delta $\\pi_{max}$  (log10)")) +
        theme_bw(base_size = 11)
      
      resid_long <- data.frame(
        lon    = rep(df$lon, 2), lat = rep(df$lat, 2),
        resid  = c(df$resid_pearson, df$resid_sar),
        modelo = rep(c("Logístico estándar", "SAR-logit"), each = n)
      )
      p5 <- ggplot(resid_long, aes(lon, lat, colour = resid)) +
        geom_point(size = 0.9, alpha = 0.6) +
        scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                               midpoint = 0, name = "Residuo\nPearson") +
        facet_wrap(~modelo) +
        labs(title = "Residuos de Pearson: antes vs. después SAR-logit",
             x = "Longitud", y = "Latitud") +
        theme_bw(base_size = 10)
      
      if (requireNamespace("gridExtra", quietly=TRUE)) {
        print(gridExtra::grid.arrange(p1, p2, ncol=2))
        print(gridExtra::grid.arrange(p3, p4, ncol=2))
      } else if (requireNamespace("patchwork", quietly=TRUE)) {
        print(p1 + p2)
        print(p3 + p4)
      } else {
        print(p1); print(p2); print(p3); print(p4)
      }
      print(p5)
      
      plots <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5)
      ggsave(file.path(DIR_OUT,paste0("pi_final_",especie_tratada,".png")), p1, width=10, height=5, dpi=150)
      ggsave(file.path(DIR_OUT,paste0("logvsSar_",especie_tratada,".png")), p2, width=10, height=5, dpi=150)
      ggsave(file.path(DIR_OUT,paste0("histograma_ipw_",especie_tratada,".png")), p3, width=10, height=5, dpi=150)
      ggsave(file.path(DIR_OUT,paste0("convergencia_",especie_tratada,".png")), p4, width=10, height=5, dpi=150)
      ggsave(file.path(DIR_OUT,paste0("residuosP_",especie_tratada,".png")), p5, width=10, height=5, dpi=150)
      
      
    } else {
      
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
      
      col_pi <- colorRampPalette(c("yellow","orange","red","darkred"))(20)
      cuts   <- cut(df$pi_final, 20)
      plot(df$lon, df$lat, col = col_pi[as.integer(cuts)],
           pch = 16, cex = 0.5,
           main = paste("pi final:", modelo_final),
           xlab = "Longitud", ylab = "Latitud")
      
      plot(df$pi_logit, df$pi_sar, pch = 16, cex = 0.4,
           col = ifelse(df$M == 1, "#006400", "#999999"),
           main = "Logístico vs. SAR-logit",
           xlab = "pi logístico", ylab = "pi SAR")
      abline(0, 1, lty = 2, col = "black")
      
      hist(w_stab, breaks = 30, col = "#006400", border = "white",
           main = sprintf("IPW estabilizados  (ESS = %.0f, %.0f%%)",
                          ESS, ESS_pct),
           xlab = "w_i", ylab = "Conteo")
      abline(v = 1, col = "red", lty = 2)
      
      plot(sar$hist_diff, type = "l", col = "#006400", lwd = 2, log = "y",
           main = "Convergencia SAR-logit",
           xlab = "Iteración", ylab = "Delta pi_max (log10)")
      abline(h = 1e-5, col = "red", lty = 2)
      
      par(mfrow = c(1, 1))
    }
    
    cat("  Graficos generados.\n")
    
  } else {
    cat("  Gráficos omitidos por generar_graficos = FALSE.\n")
  }
  
  # =============================================================================
  # 10.  TABLA RESUMEN DEL PIPELINE
  # =============================================================================
  
  cat("\n=== SECCIÓN 10: Tabla resumen ===========================================\n")
  
  resumen <- data.frame(
    modelo = nombre_modelo,
    formula = paste(deparse(formula_std), collapse = " "),
    Paso = c(
      "Datos","Datos","Datos",
      "Escalon 1","Escalon 1","Escalon 1",
      "Diagnostico 1","Diagnostico 1",
      "Escalon 2","Escalon 2","Escalon 2","Escalon 2","Escalon 2",
      "Diagnostico 2","Diagnostico 2",
      "Escalon 3","Escalon 3","Escalon 3",
      "Final"
    ),
    Estadistico = c(
      "n total","n presencias","Prevalencia",
      "AIC logístico","Pseudo-R2 McFadden","pi media logit",
      "Moran I (logístico)","p-valor Moran",
      "rho_a","L Lipschitz","|rho_a|<4 (v6) OK","Convergencia","AIC SAR",
      "Moran I (SAR)","p-valor Moran SAR",
      "w99 (p99)","% truncados","ESS %",
      "Modelo pi final"
    ),
    Valor = c(
      n, sum(df$M==1), round(mean(df$M),3),
      round(aic_logit,2), round(pR2_logit,4), round(mean(df$pi_logit),4),
      round(MI_logit,4), formatC(pval_logit, format="e", digits=3),
      round(sar$rho_a,5), round(sar$L,5),
      ifelse(sar$L<1,"SI","NO"), ifelse(sar$converged,"SI","NO"),
      round(aic_sar,2),
      round(MI_sar,4), formatC(pval_sar, format="e", digits=3),
      round(w99,4), paste0(round(100*n_trunc/n_pres_eff,1),"%"),
      paste0(round(ESS_pct,1),"%"),
      modelo_final
    ), stringsAsFactors = FALSE
  )
  
  print(resumen, row.names = FALSE)
  
  # Resumen ancho para comparar modelos
  resumen_comparacion <- data.frame(
    modelo = nombre_modelo,
    formula = paste(deparse(formula_std), collapse = " "),
    n_total = n,
    n_presencias = sum(df$M == 1),
    prevalencia = mean(df$M),
    AIC_logit = aic_logit,
    pseudoR2_logit = pR2_logit,
    pi_media_logit = mean(df$pi_logit),
    Moran_logit_I = as.numeric(MI_logit),
    Moran_logit_p = pval_logit,
    rho_a = as.numeric(sar$rho_a),
    L = as.numeric(sar$L),
    contraccion_OK = sar$L < 1,
    convergencia_SAR = sar$converged,
    n_iter_SAR = sar$n_iter,
    AIC_SAR = aic_sar,
    Moran_SAR_I = as.numeric(MI_sar),
    Moran_SAR_p = pval_sar,
    reduccion_Moran_pct = as.numeric(reduccion),
    w99 = as.numeric(w99),
    pct_truncados = 100 * n_trunc / n_pres_eff,
    ESS = ESS,
    ESS_pct = ESS_pct,
    modelo_final = modelo_final
  )
  
  # =============================================================================
  # 11.  EXPORTAR RESULTADOS
  # =============================================================================
  
  if (exportar) {
    cat("\n=== SECCIÓN 11: Exportar ================================================\n")
    
    cols_export <- c(
      "lon","lat","Year","tipo","M",
      "pi_logit","pi_sar","pi_final","w_ipw",
      "resid_pearson","resid_sar"
      #"w_sample"
    )
    cols_export <- cols_export[cols_export %in% names(df)]
    
    df_out <- df[, cols_export]
    
    archivo_csv <- file.path(carpeta_salida, paste0("pi_estimado_", nombre_modelo,"_",especie_tratada, ".csv"))
    archivo_rds <- file.path(carpeta_salida, paste0("params_pi_model_", nombre_modelo,"_",especie_tratada, ".rds"))
    
    write.csv(df_out, archivo_csv, row.names = FALSE)
    
    params <- list(
      modelo_final=modelo_final, rho_a=sar$rho_a, L=sar$L,
      lambda_max_W=lambda_max_W, convergido=sar$converged,
      n_iter=sar$n_iter, aic_logit=aic_logit, aic_sar=aic_sar,
      moran_logit_I=MI_logit, moran_logit_p=pval_logit,
      moran_sar_I=MI_sar, moran_sar_p=pval_sar,
      w99=w99, pct_trunc=100*n_trunc/n_pres_eff, ESS=ESS, ESS_pct=ESS_pct,
      formula_std=formula_std,
      resumen_comparacion=resumen_comparacion
    )
    saveRDS(params, archivo_rds)
    
    cat("  Archivo exportado:", archivo_csv, "\n")
    cat("  Parámetros exportados:", archivo_rds, "\n")
  }
  
  cat("=========================================================================\n")
  cat(sprintf("  PIPELINE COMPLETADO — modelo: %s\n", modelo_final))
  cat("=========================================================================\n")
  
  return(list(
    df = df,
    formula = formula_std,
    logit = mod_logit,
    sar = sar,
    lw = lw,
    moran_logit = moran_logit,
    moran_sar = moran_sar,
    calibracion = cal,
    resumen = resumen,
    resumen_comparacion = resumen_comparacion,
    plots = plots
  ))
}


# EJECUCION CON VARIAS FORMULAS PROBADAS
formulas_pi <- list(
  # Modelo espacial puro (referencia — modelo actual)
  f5_sin_watdist = M ~ dist_road_z + dist_res_z + slope_z + tri_z + 
    aspect_sin + aspect_cos + hillshade_z + landcover,
  
  # Modelo espacio-temporal: añade year_norm como covariable de detección
  # Hipótesis: la probabilidad de detección aumenta con el año
  # (ciencia ciudadana post-2014 aumenta muestreo)
  f5_temporal = M ~ dist_road_z + dist_res_z + slope_z + tri_z + 
    aspect_sin + aspect_cos + hillshade_z + landcover + year_norm,
  
  # Modelo temporal mínimo: solo accesibilidad + año
  # Más parsimonioso — útil si las bioclimáticas no aportan a la detección
  f_temporal_min = M ~ dist_road_z + dist_res_z + year_norm,
  
  f7 = M ~ dist_road_z + dist_res_z + tri_z +
    aspect_sin + aspect_cos + hillshade_z
)

bases_pi <- list(
  MIX = preparar_df_pi(datos_modelo, quitar_na = TRUE),
  NN  = preparar_df_pi(datos_modelo_nn, quitar_na = TRUE)
)

resultados_pi <- list()

for (b in names(bases_pi)) {
  for (f in names(formulas_pi)) {
    
    nombre <- paste(b, f, sep = "_")
    
    resultados_pi[[nombre]] <- correr_pipeline_pi(
      df = bases_pi[[b]],
      formula_std = formulas_pi[[f]],
      nombre_modelo = nombre,
      k = 5,
      umbral = 0.05,
      generar_graficos = FALSE,
      exportar = FALSE,
      ecuador = NULL
    )
  }
}

tabla_comparacion_pi <- dplyr::bind_rows(
  lapply(resultados_pi, function(x) x$resumen_comparacion)
)

tabla_comparacion_pi %>%
  dplyr::arrange(desc(ESS_pct), AIC_SAR)


#CORRER SOLO CON EL MEJOR MODELO ESPACIAL (referencia)
mejor_espacial <- correr_pipeline_pi(
  df = bases_pi$NN,
  formula_std = formulas_pi$f5_sin_watdist,
  nombre_modelo = "NN_f5_sin_watdist",
  k = 5,
  umbral = 0.05,
  generar_graficos = TRUE,
  exportar = TRUE,
  carpeta_salida = ".",
  ecuador = ecuador
)

# CORRER CON MODELO ESPACIO-TEMPORAL
# Incluye year_norm para capturar el sesgo MNAR temporal (ciencia ciudadana)
mejor_temporal <- correr_pipeline_pi(
  df = bases_pi$NN,
  formula_std = formulas_pi$f5_temporal,
  nombre_modelo = "NN_f5_temporal",
  k = 5,
  umbral = 0.05,
  generar_graficos = TRUE,
  exportar = TRUE,
  carpeta_salida = ".",
  ecuador = ecuador
)

# Comparación diagnóstica: ¿el modelo temporal produce mayor variación en pi_hat?
cat("\n=== COMPARACION MODELOS ESPACIAL vs ESPACIO-TEMPORAL ==================\n")
pi_esp  <- mejor_espacial$df$pi_final[mejor_espacial$df$M == 1]
pi_temp <- mejor_temporal$df$pi_final[mejor_temporal$df$M == 1]
w_esp   <- mejor_espacial$df$w_ipw[mejor_espacial$df$M == 1]
w_temp  <- mejor_temporal$df$w_ipw[mejor_temporal$df$M == 1]
cat(sprintf("  pi_hat espacial  : rango [%.3f, %.3f]  CV_w=%.3f\n",
            min(pi_esp), max(pi_esp), sd(w_esp)/mean(w_esp)))
cat(sprintf("  pi_hat temporal  : rango [%.3f, %.3f]  CV_w=%.3f\n",
            min(pi_temp), max(pi_temp), sd(w_temp)/mean(w_temp)))
# Test: ¿pi_hat difiere entre pre y post 2014?
yr_esp  <- mejor_espacial$df$Year[mejor_espacial$df$M == 1]
yr_temp <- mejor_temporal$df$Year[mejor_temporal$df$M == 1]
t_esp  <- t.test(pi_esp  ~ (yr_esp  >= 2014))
t_temp <- t.test(pi_temp ~ (yr_temp >= 2014))
cat(sprintf("  t-test pre/post 2014 — espacial : p=%.4f  (%.3f vs %.3f)\n",
            t_esp$p.value, t_esp$estimate[1], t_esp$estimate[2]))
cat(sprintf("  t-test pre/post 2014 — temporal : p=%.4f  (%.3f vs %.3f)\n",
            t_temp$p.value, t_temp$estimate[1], t_temp$estimate[2]))
cat("  Si p < 0.05 en temporal pero no en espacial: el modelo temporal\n")
cat("  captura el sesgo MNAR que el espacial omitia.\n")
cat("=========================================================================\n")

# Usar el modelo temporal para el test GET
# Sustituir 'mejor' por el modelo que tenga mayor CV_w
mejor <- if (sd(w_temp)/mean(w_temp) > sd(w_esp)/mean(w_esp)) {
  cat("  -> Usando modelo TEMPORAL para el test GET (mayor variacion en pi_hat)\n")
  mejor_temporal
} else {
  cat("  -> Usando modelo ESPACIAL para el test GET (no hay ganancia con year_norm)\n")
  mejor_espacial
}
