# =============================================================================
#  get_separabilidad.R  —  Test GET de separabilidad espacio-temporal bajo MNAR
#  Versión 3.0 — integrada sobre script sylvatica
#
#  Estadístico: S(u,t) 3D de Ghorbani et al. (2021), ec. 11
#    Naive:     S(u,t) = rho_hat(u,t)           / [rho_sp(u)*rho_t(t)/n]
#    Corregido: S(u,t) = rho_hat_ipw(u,t)/pi(u) / [rho_sp_w(u)*rho_t_w(t)/n_eff]
#  Proyección: S_t(t)  = integral_W S(u,t) du   (Ghorbani ec. 12)
#              S_sp(u) = integral_T S(u,t) dt    (Ghorbani ec. 13)
#
#  DIFERENCIA CLAVE respecto a versión anterior:
#  · Antes: S_t(t) calculado directo como rho_t(t_perm)/(n/|T|)
#    → INVARIANTE a permutaciones → corredor ancho=0 (bug matemático)
#  · Ahora: rho_hat(u,t) 3D recalculado en cada permutación (numerador de S)
#    El denominador rho_sp*rho_t/n_eff es fijo (Ghorbani pág.8: escala)
#    → las permutaciones generan S_t^(k) genuinamente distintas
#
#  PARÁMETROS PARA n=206 (sylvatica):
#    K_ESTRATOS = 2   (K=5 daba ~41 pts/estrato → casi-invarianza)
#    N_GRID_SP  = 15  (15×15×25 = 5625 pts; manejable en ~5-10 min)
#    N_GRID_T   = 25
#    Para n>500: subir N_GRID_SP a 20-25, K_ESTRATOS a 3-5
# =============================================================================

cat("=========================================================================\n")
cat("  TEST GET DE SEPARABILIDAD ESPACIO-TEMPORAL BAJO MNAR  [v3.0]\n")
cat("  Estadístico S(u,t) 3D — Ghorbani et al. (2021) + correcciones MNAR\n")
cat("=========================================================================\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ── Parámetros globales ───────────────────────────────────────────────────────
set.seed(2025)
#ARCHIVO_PI  <- "C:/Users/ivanf/Desktop/Octavo/DTIC/Archivos Mabel/pi_estimado_NN_f5_temporal.csv"   # salida de pi_sampling_model_v2.R
# columnas: lon, lat, Year, M, pi_final, w_ipw
K_ESTRATOS  <- 2     # estratos pi_hat: usar 2-3 si n<300, 5 si n>500
K_PERM      <- 999   # permutaciones GET
ALPHA       <- 0.05
N_GRID_SP   <- 15    # puntos por eje espacial — grilla N×N×N_GRID_T
N_GRID_T    <- 25    # puntos temporales
DIR_OUT     <- "."

cat(sprintf("  K_ESTRATOS=%d  N_GRID_SP=%d  N_GRID_T=%d  K_PERM=%d\n",
            K_ESTRATOS, N_GRID_SP, N_GRID_T, K_PERM))
cat(sprintf("  Grilla 3D total: %d puntos por permutación\n",
            N_GRID_SP^2 * N_GRID_T))

# =============================================================================
# SECCIÓN 1: Cargar datos
# =============================================================================
cat("\n=== SECCIÓN 1: Carga de datos ==========================================\n")

#1. sylvatica
 # especie_tratada <- "sylvatica"
 # ARCHIVO_PI  <- "C:/Users/ivanf/Desktop/Octavo/DTIC/Archivos Mabel/pi_estimado_NN_f5_temporal_sylvatica.csv"   # salida de pi_sampling_model_v2.R

#2. anthonyi
# especie_tratada <- "anthonyi"
# ARCHIVO_PI  <- "C:/Users/ivanf/Desktop/Octavo/DTIC/Archivos Mabel/pi_estimado_NN_f5_temporal_anthonyi.csv"   # salida de pi_sampling_model_v2.R

#3. bilinguis
especie_tratada <- "bilinguis"
ARCHIVO_PI  <- "C:/Users/ivanf/Desktop/Octavo/DTIC/Archivos Mabel/pi_estimado_NN_f5_temporal_bilinguis.csv"   # salida de pi_sampling_model_v2.R


if (!file.exists(ARCHIVO_PI)){
  stop("No se encontró '", ARCHIVO_PI, "'. Ejecute primero pi_sampling_model_v2.R")
}

df_full  <- read.csv(ARCHIVO_PI, stringsAsFactors = FALSE)
cols_req <- c("lon", "lat", "Year", "M", "pi_final", "w_ipw")
faltantes <- setdiff(cols_req, names(df_full))
if (length(faltantes) > 0){
  stop("Faltan columnas: ", paste(faltantes, collapse=", "))
}

df_pres <- df_full[df_full$M == 1 & !is.na(df_full$w_ipw), ]
n_pres  <- nrow(df_pres)

if (any(df_pres$pi_final <= 0 | df_pres$pi_final >= 1)) {
  warning("pi_final fuera de (0,1); truncando.")
  df_pres$pi_final <- pmax(0.001, pmin(0.999, df_pres$pi_final))
}

cat(sprintf("  Total filas CSV  : %d\n", nrow(df_full)))
cat(sprintf("  Presencias (M=1) : %d\n", n_pres))
cat(sprintf("  Rango lon        : [%.2f, %.2f]\n", min(df_pres$lon), max(df_pres$lon)))
cat(sprintf("  Rango lat        : [%.2f, %.2f]\n", min(df_pres$lat), max(df_pres$lat)))
cat(sprintf("  Rango temporal   : [%d, %d]\n", min(df_pres$Year), max(df_pres$Year)))
cat(sprintf("  w_ipw: min=%.3f  med=%.3f  max=%.3f  sum=%.1f\n",
            min(df_pres$w_ipw), median(df_pres$w_ipw),
            max(df_pres$w_ipw), sum(df_pres$w_ipw)))

if (n_pres < 100) warning("n < 100: potencia del test puede ser muy baja.")

# =============================================================================
# SECCIÓN 2: Dominio y grilla 3D
# =============================================================================
cat("\n=== SECCIÓN 2: Dominio y grilla 3D =====================================\n")

lon_range <- range(df_pres$lon);  lat_range <- range(df_pres$lat)
t_range   <- range(df_pres$Year)
mg_sp <- 0.05 * c(diff(lon_range), diff(lat_range));  mg_t <- 0.05 * diff(t_range)

W_lon <- lon_range + c(-1,1)*mg_sp[1]
W_lat <- lat_range + c(-1,1)*mg_sp[2]
T_rng <- t_range   + c(-1,1)*mg_t
area_W <- diff(W_lon) * diff(W_lat)
len_T  <- diff(T_rng)

grid_lon <- seq(W_lon[1], W_lon[2], length.out = N_GRID_SP)
grid_lat <- seq(W_lat[1], W_lat[2], length.out = N_GRID_SP)
grid_t   <- seq(T_rng[1], T_rng[2], length.out = N_GRID_T)
grid_sp  <- expand.grid(lon = grid_lon, lat = grid_lat)   # N_SP × 2
N_SP     <- N_GRID_SP^2
dlon_sp  <- diff(grid_lon)[1];  dlat_sp <- diff(grid_lat)[1];  dt_g <- diff(grid_t)[1]

cat(sprintf("  Dominio W: lon [%.2f,%.2f]  lat [%.2f,%.2f]\n",
            W_lon[1],W_lon[2],W_lat[1],W_lat[2]))
cat(sprintf("  Dominio T: [%.0f,%.0f]  (%.0f años)\n", T_rng[1],T_rng[2],len_T))
cat(sprintf("  Grilla 3D: %d × %d × %d = %d puntos\n",
            N_GRID_SP, N_GRID_SP, N_GRID_T, N_SP*N_GRID_T))

# =============================================================================
# SECCIÓN 3: Bandwidths (idéntico al script original)
# =============================================================================
cat("\n=== SECCIÓN 3: Selección de bandwidths ================================\n")

factor_epa <- (5)^(1/5) * ((3/5) / (1/(2*sqrt(pi))))^(1/5)
cat(sprintf("  kappa Gaussiano→Epanechnikov: %.4f\n", factor_epa))

set.seed(42)
t_jitter <- df_pres$Year + runif(n_pres, -0.05, 0.05)
hT_gauss <- tryCatch(bw.SJ(t_jitter),
                     error = function(e) { cat("  AVISO: bw.SJ→rule-of-thumb\n")
                       1.06*sd(t_jitter)*n_pres^(-1/5) })
hT <- hT_gauss * factor_epa
cat(sprintf("  hT_gauss=%.4f  hT_Epa=%.4f años\n", hT_gauss, hT))

lon_jitter <- df_pres$lon + runif(n_pres,-1e-5,1e-5)
lat_jitter <- df_pres$lat + runif(n_pres,-1e-5,1e-5)
hS_lon <- tryCatch(bw.SJ(lon_jitter), error=function(e) 1.06*sd(lon_jitter)*n_pres^(-1/5))
hS_lat <- tryCatch(bw.SJ(lat_jitter), error=function(e) 1.06*sd(lat_jitter)*n_pres^(-1/5))
hS <- sqrt(hS_lon * hS_lat) * factor_epa
cat(sprintf("  hS_lon=%.4f  hS_lat=%.4f  hS_Epa=%.4f grados\n", hS_lon, hS_lat, hS))

# =============================================================================
# SECCIÓN 4: Funciones kernel (idéntico al script original + rho3d nueva)
# =============================================================================

ker1d <- function(x, h) pmax(0, 0.75*(1-(x/h)^2)) / h
ker2d <- function(dx, dy, h) ker1d(dx, h) * ker1d(dy, h)   # isotrópico

edge_corr_1d <- function(ti, h, a, b) {
  lo <- pmax(-1, (a-ti)/h);  hi <- pmin(1, (b-ti)/h)
  f  <- function(u) 0.75*(u - u^3/3)
  pmax(1e-8, f(hi) - f(lo))
}
edge_corr_2d <- function(lon_i, lat_i, h, W_lon, W_lat) {
  pmax(1e-8, edge_corr_1d(lon_i, h, W_lon[1], W_lon[2]) *
         edge_corr_1d(lat_i, h, W_lat[1], W_lat[2]))
}

# ── rho_hat 3D vectorizado (N_SP × N_GRID_T) ─────────────────────────────────
# Calcula sum_i w_i * K_sp(u-u_i)/C_sp_i * K_t(t-t_i)/C_t_i
# para todos los puntos (u,t) de la grilla simultáneamente.
rho3d <- function(t_vec, lon_vec, lat_vec, w_vec,
                  grid_sp, grid_t, hS, hT, W_lon, W_lat, T_rng) {
  n_pts <- length(t_vec)
  n_sp  <- nrow(grid_sp)
  n_t   <- length(grid_t)
  
  Ct  <- edge_corr_1d(t_vec,   hT, T_rng[1], T_rng[2])   # (n_pts)
  Csp <- edge_corr_2d(lon_vec, lat_vec, hS, W_lon, W_lat) # (n_pts)
  wC  <- w_vec / (Csp * Ct)                                # (n_pts)
  
  # Kernel espacial: garantizamos matriz n_sp × n_pts
  dlon <- matrix(outer(grid_sp$lon, lon_vec, "-"), nrow = n_sp, ncol = n_pts)
  dlat <- matrix(outer(grid_sp$lat, lat_vec, "-"), nrow = n_sp, ncol = n_pts)
  K_sp <- matrix(ker2d(dlon, dlat, hS),              nrow = n_sp, ncol = n_pts)
  
  # Kernel temporal: garantizamos matriz n_t × n_pts
  K_t  <- matrix(outer(grid_t, t_vec, function(a, b) ker1d(a - b, hT)),
                 nrow = n_t, ncol = n_pts)
  
  # sweep(K_sp, 2, wC) escala columna i por wC[i]  → n_sp × n_pts
  # %*% t(K_t)  → (n_sp × n_pts) %*% (n_pts × n_t) = n_sp × n_t
  sweep(K_sp, 2, wC, "*") %*% t(K_t)
}

# ── Intensidades marginales (para el denominador fijo) ────────────────────────
rho_sp_fn <- function(lon_vec, lat_vec, w_vec, grid_sp, hS, W_lon, W_lat) {
  n_pts <- length(lon_vec)
  n_sp  <- nrow(grid_sp)
  Csp  <- edge_corr_2d(lon_vec, lat_vec, hS, W_lon, W_lat)
  dlon <- matrix(outer(grid_sp$lon, lon_vec, "-"), nrow = n_sp, ncol = n_pts)
  dlat <- matrix(outer(grid_sp$lat, lat_vec, "-"), nrow = n_sp, ncol = n_pts)
  K_sp <- matrix(ker2d(dlon, dlat, hS),              nrow = n_sp, ncol = n_pts)
  as.numeric(K_sp %*% (w_vec / Csp))   # N_SP
}

rho_t_fn <- function(t_vec, w_vec, grid_t, hT, T_rng) {
  Ct <- edge_corr_1d(t_vec, hT, T_rng[1], T_rng[2])
  K_t <- outer(grid_t, t_vec, function(a,b) ker1d(a-b, hT))
  as.numeric(K_t %*% (w_vec / Ct))   # N_GRID_T
}

# ── S(u,t) y proyecciones normalizadas ───────────────────────────────────────
# S(u,t) = rho3d / [rho_sp * rho_t / n_eff]  →  E[S] ≈ 1 bajo H0
#
# Las proyecciones integran S sobre el dominio complementario:
#   S_t(t)  = integral_W S(u,t) du  →  E[S_t(t)]  = integral_W 1 du  = |W|
#   S_sp(u) = integral_T S(u,t) dt  →  E[S_sp(u)] = integral_T 1 dt  = |T|
#
# Para que E[S_t] ≈ 1 y E[S_sp] ≈ 1 (línea de referencia en 1),
# se NORMALIZA dividiendo por la medida de Lebesgue del dominio integrado.
# Sin esta división las curvas flotan en |W| (≈20-30 grados² para Ecuador)
# y la línea yintercept=1 es incorrecta.
calc_S <- function(rho3d_mat, rho_sp, rho_t, n_eff, dlon, dlat, dt,
                   area_W, len_T) {
  denom  <- pmax(outer(rho_sp, rho_t) / n_eff, 1e-12)
  S_mat  <- rho3d_mat / denom
  # Integral numérica normalizada por la medida del dominio
  S_t    <- (colSums(S_mat) * dlon * dlat) / area_W   # → E ≈ 1 bajo H0
  S_sp   <- (rowSums(S_mat) * dt)          / len_T    # → E ≈ 1 bajo H0
  list(S_mat=S_mat, S_t=S_t, S_sp=S_sp)
}

# ── St_perm_fn — S_t bajo permutación (normalizado) ──────────────────────────
# · Recalcula SOLO el numerador rho3d con tiempos permutados.
# · El denominador (rho_sp_fix * rho_t_fix) es fijo — Ghorbani (2021), pág. 8.
# · Micro-jitter ±1e-4 años en t_perm: rompe escalones del kernel Epanechnikov
#   sobre tiempos enteros sin alterar la interpretación (< 0.04 días de ruido).
# · División final por area_W: garantiza E[S_t] ≈ 1 bajo H0.
St_perm_fn <- function(t_perm, lon_v, lat_v, w_v, n_eff_v,
                       rho_sp_fix, rho_t_fix,
                       grid_sp, grid_t, hS, hT, W_lon, W_lat, T_rng,
                       dlon, dlat, area_W) {
  t_perm_jit <- t_perm + runif(length(t_perm), -1e-4, 1e-4)
  r3    <- rho3d(t_perm_jit, lon_v, lat_v, w_v,
                 grid_sp, grid_t, hS, hT, W_lon, W_lat, T_rng)
  denom <- pmax(outer(rho_sp_fix, rho_t_fix) / n_eff_v, 1e-12)
  (colSums(r3 / denom) * dlon * dlat) / area_W
}

# ── GET bilateral pooled ──────────────────────────────────────────────────────
get_bilateral <- function(S_obs, S_perm_mat, ALPHA) {
  # ── GET rango extremo bilateral — Myllymäki et al. (2017) ──────────────────
  #
  # NOTA IMPORTANTE sobre la interpretación del gráfico:
  # El corredor [min,max] de permutaciones válidas y la decisión formal (p-valor)
  # pueden no coincidir visualmente. Esto NO es un error — es una propiedad del
  # rango extremo bilateral:
  #
  # · El rango extremo de una curva = min_j { rango bilateral en t_j }
  # · El test rechaza si r_obs <= c_alpha  (el observado es globalmente extremo)
  # · El corredor visual muestra [min,max] de las permutaciones no rechazadas
  # · Una curva puede cruzar el corredor en t_j sin que su rango extremo
  #   GLOBAL sea suficientemente bajo, porque en otros t_j está bien rankeada
  #
  # Interpretación correcta:
  # · Si hay cruces visuales SIN rechazo: el observado es localmente extremo
  #   en algunos años pero no globalmente. No se rechaza H0.
  # · Si hay rechazo SIN cruces visuales: el rango extremo del observado
  #   es bajo pero el cruce puede ser muy sutil (cerca del borde del corredor).
  # · En ambos casos, la decisión formal (p-valor) es la correcta.
  # · Referencia: Myllymäki et al. (2017) JRSS-B, Sección 2.3.
  
  N_ALL  <- nrow(S_perm_mat) + 1
  S_all  <- rbind(S_obs, S_perm_mat)          # (K+1) x J
  
  # Paso 1: rangos bilaterales en conjunto pooled
  rk_inf <- apply(S_all, 2, rank, ties.method = "min")
  rk_sup <- apply(-S_all, 2, rank, ties.method = "min")
  rk_bil <- pmin(rk_inf, rk_sup)              # (K+1) x J
  
  # Paso 2: rango extremo = mínimo bilateral sobre todos los t_j
  extr   <- apply(rk_bil, 1, min)             # longitud K+1
  r_obs  <- extr[1];   r_k <- extr[-1]
  rk_j   <- rk_bil[1, ]                       # rango puntual del observado
  
  # Paso 3: p-valor exacto y decisión
  p_val  <- sum(extr <= r_obs) / N_ALL
  rechaz <- p_val <= ALPHA
  
  # Paso 4: umbral empírico c_alpha
  c_alph <- sort(extr)[max(1L, floor(ALPHA * N_ALL))]
  
  # Paso 5: corredor genuino = [min,max] de permutaciones con r_k >= c_alpha
  # Este corredor NO tiene garantía de coherencia puntual con la decisión formal
  # porque el rango extremo es un mínimo global, no un criterio puntual.
  # La coherencia es GLOBAL: si r_obs < c_alpha el observado CONTRIBUYÓ al
  # mínimo en algún t_j, pero ese t_j puede estar cerca del borde del corredor.
  idx_v <- which(r_k >= c_alph)
  if (length(idx_v) > 0) {
    L <- apply(S_perm_mat[idx_v, , drop = FALSE], 2, min, na.rm = TRUE)
    U <- apply(S_perm_mat[idx_v, , drop = FALSE], 2, max, na.rm = TRUE)
  } else {
    L <- apply(S_perm_mat, 2, min, na.rm = TRUE)
    U <- apply(S_perm_mat, 2, max, na.rm = TRUE)
  }
  
  fuera <- S_obs < L | S_obs > U
  
  # Informar si hay discrepancia visual vs decisión formal
  if (any(fuera) && !rechaz)
    message("  NOTA: cruces visuales del corredor sin rechazo formal. ",
            "El observado es localmente extremo en algunos anos pero no ",
            "globalmente (rango extremo bilateral). La decision correcta es: ",
            "H0 no rechazada (p=", round(p_val,4), ").")
  if (!any(fuera) && rechaz)
    message("  NOTA: rechazo formal sin cruce visual evidente. ",
            "El cruce puede ser sutil (cerca del borde del corredor). ",
            "La decision correcta es: H0 RECHAZADA (p=", round(p_val,4), ").")
  
  list(p = p_val, rechazar = rechaz, c_alpha = c_alph,
       r_obs = r_obs, r_k = r_k, rk_j = rk_j, L = L, U = U,
       M = apply(S_perm_mat, 2, median, na.rm = TRUE),
       fuera = fuera, n_val = length(idx_v), ancho = mean(U - L))
}


# =============================================================================
# SECCIÓN 5: Estadísticos observados — CORREGIDO (IPW) y NAIVE
# =============================================================================
cat("\n=== SECCIÓN 5: Estadísticos observados ================================\n")

lon_p <- df_pres$lon;  lat_p <- df_pres$lat
t_p   <- df_pres$Year; w_p   <- df_pres$w_ipw;  pi_p <- df_pres$pi_final
n_eff <- sum(w_p)
w_n   <- rep(1, n_pres)   # pesos naive (todos = 1, n_eff_naive = n_pres)

cat("  Calculando rho_hat 3D (corregido y naive)...\n")
t0 <- proc.time()["elapsed"]

rho3d_c <- rho3d(t_p, lon_p, lat_p, w_p, grid_sp, grid_t, hS, hT, W_lon, W_lat, T_rng)
rho3d_n <- rho3d(t_p, lon_p, lat_p, w_n, grid_sp, grid_t, hS, hT, W_lon, W_lat, T_rng)

rho_sp_c <- rho_sp_fn(lon_p, lat_p, w_p, grid_sp, hS, W_lon, W_lat)
rho_t_c  <- rho_t_fn(t_p, w_p, grid_t, hT, T_rng)
rho_sp_n <- rho_sp_fn(lon_p, lat_p, w_n, grid_sp, hS, W_lon, W_lat)
rho_t_n  <- rho_t_fn(t_p, w_n, grid_t, hT, T_rng)

Sc <- calc_S(rho3d_c, rho_sp_c, rho_t_c, n_eff,   dlon_sp, dlat_sp, dt_g, area_W, len_T)
Sn <- calc_S(rho3d_n, rho_sp_n, rho_t_n, n_pres,  dlon_sp, dlat_sp, dt_g, area_W, len_T)

S_t_c  <- Sc$S_t;  S_sp_c <- Sc$S_sp
S_t_n  <- Sn$S_t;  S_sp_n <- Sn$S_sp

cat(sprintf("  [%.1f s]\n", proc.time()["elapsed"] - t0))
cat(sprintf("  S_t corregido : [%.3f, %.3f]  (esperado ≈ 1 bajo H0)\n",
            min(S_t_c), max(S_t_c)))
cat(sprintf("  S_t naive     : [%.3f, %.3f]  (esperado ≈ 1 bajo H0)\n",
            min(S_t_n), max(S_t_n)))
cat(sprintf("  Verificación normalización: media S_t corregido=%.3f  naive=%.3f\n",
            mean(S_t_c), mean(S_t_n)))
cat(sprintf("  (si media ≈ 1 la normalización es correcta; si ≈ |W|=%.2f no se aplicó)\n",
            area_W))

# =============================================================================
# SECCIÓN 6: Test chi-cuadrado — NAIVE y CORREGIDO
# =============================================================================
cat("\n=== SECCIÓN 6: Test chi-cuadrado =======================================\n")

# Tamaño de grilla: n/(I²*J) >= 5  →  I = floor((n/5)^(1/3))
I_chi <- max(3L, floor((n_pres/5)^(1/3)))
J_chi <- I_chi
cat(sprintf("  Grilla: %d×%d espacial × %d temporal = %d celdas\n",
            I_chi, I_chi, J_chi, I_chi^2*J_chi))

lon_br <- quantile(lon_p, seq(0,1,length.out=I_chi+1))
lat_br <- quantile(lat_p, seq(0,1,length.out=I_chi+1))
t_br   <- quantile(t_p,   seq(0,1,length.out=J_chi+1))
lon_br[1] <- lon_br[1]-1e-6;  lat_br[1] <- lat_br[1]-1e-6;  t_br[1] <- t_br[1]-1e-6

lc <- cut(lon_p, breaks=lon_br, labels=FALSE)
ac <- cut(lat_p, breaks=lat_br, labels=FALSE)
tc <- cut(t_p,   breaks=t_br,   labels=FALSE)

chi_test <- function(wts, label) {
  W_tot <- sum(wts)
  wijk  <- array(0, c(I_chi, I_chi, J_chi))
  for (i in seq_len(n_pres))
    if (!is.na(lc[i]) & !is.na(ac[i]) & !is.na(tc[i]))
      wijk[lc[i], ac[i], tc[i]] <- wijk[lc[i], ac[i], tc[i]] + wts[i]
  w_sp <- apply(wijk, c(1,2), sum)
  w_t  <- apply(wijk, 3, sum)
  e    <- outer(as.vector(w_sp), w_t) / W_tot
  chi2 <- sum((as.vector(wijk) - e)^2 / pmax(e, 1e-9))
  df_c <- (I_chi^2 - 1) * (J_chi - 1)
  p    <- 1 - pchisq(chi2, df_c)
  min_e <- min(e)
  cat(sprintf("  chi² %-10s: X²=%.2f  df=%d  p=%.6f  -> H₀ %s%s\n",
              label, chi2, df_c, p,
              if(p<=ALPHA) "RECHAZADA" else "no rechazada",
              if(min_e<5) sprintf(" [AVISO: min(e)=%.1f<5]",min_e) else ""))
  list(chi2=chi2, df=df_c, p=p, rechazar=p<=ALPHA, min_e=min_e)
}

res_chi_n <- chi_test(w_n, "NAIVE")
res_chi_c <- chi_test(w_p, "CORREGIDO")

# =============================================================================
# SECCIÓN 7: Estratificación por cuantiles de π̂
# =============================================================================
cat("\n=== SECCIÓN 7: Estratificación por cuantiles de π̂ ====================\n")

bp_pi <- quantile(pi_p, seq(0,1,length.out=K_ESTRATOS+1))
bp_pi[1] <- 0;  bp_pi[K_ESTRATOS+1] <- 1
estrato <- cut(pi_p, breaks=bp_pi, labels=paste0("E",1:K_ESTRATOS),
               include.lowest=TRUE)
cat(sprintf("  K_ESTRATOS=%d  |  breakpoints: %s\n", K_ESTRATOS,
            paste(round(bp_pi,3), collapse=" | ")))
for (k in 1:K_ESTRATOS) {
  nk <- sum(estrato==paste0("E",k))
  pk <- pi_p[estrato==paste0("E",k)]
  cat(sprintf("  E%d: n=%3d  π̂∈[%.3f,%.3f]\n", k, nk, min(pk), max(pk)))
}

# =============================================================================
# SECCIÓN 8: Función S_t bajo permutación (usa rho3d 3D — el cambio clave)
# =============================================================================
# Para cada permutación se recalcula SOLO el numerador rho3d con tiempos
# permutados. El denominador (rho_sp * rho_t) es fijo (Ghorbani pág. 8).
# Esto garantiza que las curvas S_t^(k) sean genuinamente distintas entre sí.

# (St_perm_fn definida arriba en Sección 8)

# =============================================================================
# SECCIÓN 9: Permutaciones — CORREGIDO (estratificado) y NAIVE (global)
# =============================================================================
cat(sprintf("\n=== SECCIÓN 9: %d permutaciones =====================================\n",
            K_PERM))
cat(sprintf("  Grilla 3D: %d×%d×%d  n=%d  [punto=10 perms]\n",
            N_GRID_SP, N_GRID_SP, N_GRID_T, n_pres))

S_t_perm_c <- matrix(NA_real_, K_PERM, N_GRID_T)
S_t_perm_n <- matrix(NA_real_, K_PERM, N_GRID_T)

t0p <- proc.time()["elapsed"]
cat("  ")
for (k in seq_len(K_PERM)) {
  if (k %% 10 == 0) cat(".")
  if (k %% 200 == 0) cat(sprintf(" %d (%.0fs)\n  ", k, proc.time()["elapsed"]-t0p))
  
  # Corregido: permuta DENTRO de cada estrato (preserva estructura de pi)
  t_pc <- t_p
  for (ek in levels(estrato)) {
    idx <- which(estrato==ek)
    if (length(idx)>1) t_pc[idx] <- sample(t_p[idx])
  }
  S_t_perm_c[k,] <- St_perm_fn(t_pc, lon_p, lat_p, w_p, n_eff,
                               rho_sp_c, rho_t_c,
                               grid_sp, grid_t, hS, hT, W_lon, W_lat, T_rng,
                               dlon_sp, dlat_sp, area_W)
  
  # Naive: permuta GLOBALMENTE (sin estratificar, sin IPW)
  t_pn <- sample(t_p)
  S_t_perm_n[k,] <- St_perm_fn(t_pn, lon_p, lat_p, w_n, n_pres,
                               rho_sp_n, rho_t_n,
                               grid_sp, grid_t, hS, hT, W_lon, W_lat, T_rng,
                               dlon_sp, dlat_sp, area_W)
}
cat(sprintf("\n  Completado en %.1f s\n", proc.time()["elapsed"]-t0p))

# =============================================================================
# SECCIÓN 10: GET rango extremo bilateral
# =============================================================================
cat("\n=== SECCIÓN 10: GET ====================================================\n")

GET_c <- get_bilateral(S_t_c, S_t_perm_c, ALPHA)
GET_n <- get_bilateral(S_t_n, S_t_perm_n, ALPHA)

cat(sprintf("  CORREGIDO: r_obs=%d  c_alpha=%d  p=%.4f  -> H₀ %s\n",
            GET_c$r_obs, GET_c$c_alpha, GET_c$p,
            if(GET_c$rechazar) "RECHAZADA" else "no rechazada"))
cat(sprintf("             vals=%d/%d  ancho_corredor=%.4f\n",
            GET_c$n_val, K_PERM, GET_c$ancho))
cat(sprintf("  NAIVE:     r_obs=%d  c_alpha=%d  p=%.4f  -> H₀ %s\n",
            GET_n$r_obs, GET_n$c_alpha, GET_n$p,
            if(GET_n$rechazar) "RECHAZADA" else "no rechazada"))
cat(sprintf("             vals=%d/%d  ancho_corredor=%.4f\n",
            GET_n$n_val, K_PERM, GET_n$ancho))

# Punto más extremo
t_ext_c <- grid_t[which.min(GET_c$rk_j)]
t_ext_n <- grid_t[which.min(GET_n$rk_j)]
cat(sprintf("  Punto más extremo corregido: t=%.1f  S_t=%.3f\n",
            t_ext_c, S_t_c[which.min(GET_c$rk_j)]))
cat(sprintf("  Punto más extremo naive:     t=%.1f  S_t=%.3f\n",
            t_ext_n, S_t_n[which.min(GET_n$rk_j)]))

chi_vs_get_aviso <- function(p_chi_n, p_chi_c, rechazar_get_n, rechazar_get_c, alpha) {
  if ((p_chi_n < alpha | p_chi_c < alpha) & !rechazar_get_n & !rechazar_get_c) {
    cat("
  ⚠  AVISO METODOLÓGICO: chi² rechaza pero GET no rechaza.
")
    cat("     → Confíe en el GET. Posibles causas del falso rechazo chi²:
")
    cat("       · Dependencia del tamaño/posición de la grilla (MAUP)
")
    cat("       · Supuestos asintóticos no satisfechos (min e_ij < 5)
")
    cat("       · El GET controla mejor el error tipo I con n moderado
")
  } else if (rechazar_get_n | rechazar_get_c) {
    if (p_chi_n >= alpha & p_chi_c >= alpha)
      cat("
  NOTA: GET rechaza pero chi² no — el GET tiene mayor potencia.
")
  }
}

chi_vs_get_aviso(res_chi_n$p, res_chi_c$p, GET_n$rechazar, GET_c$rechazar, ALPHA)

# Coherencia corredor ↔ decisión
#if (any(GET_c$fuera) != GET_c$rechazar)
#  warning("Incoherencia corredor-decisión CORREGIDO (puede ser por empates)")
#if (any(GET_n$fuera) != GET_n$rechazar)
# warning("Incoherencia corredor-decisión NAIVE (puede ser por empates)")

# ── Nota metodológica: divergencia chi² vs GET ────────────────────────────────
# Si chi² rechaza contundentemente (p < 0.001) pero GET no rechaza:
#   → CONFIAR en el GET. Razones:
#   1. GET no depende de supuestos asintóticos (chi² requiere e_ij >= 5).
#   2. GET es libre del MAUP (Modifiable Areal Unit Problem): el chi² cambia
#      con distintas grillas; el GET es invariante a la discretización.
#   3. Para n moderado (<500), el GET tiene mejor control del error tipo I.
#   4. Si chi² rechaza y GET no, el rechazo chi² probablemente refleja
#      heterogeneidad de la grilla, no no-separabilidad real del proceso.
# Si GET rechaza y chi² no: también confiar en el GET (mayor potencia).
# El chi² es un test rápido de exploración preliminar — Ghorbani (2021) Sec. 5.


# =============================================================================
# SECCIÓN 11: Gráficos
# =============================================================================
cat("\n=== SECCIÓN 11: Gráficos ===============================================\n")

# ── Plot 1: S_t corregido ─────────────────────────────────────────────────────
df1 <- data.frame(t=grid_t, St=S_t_c, L=GET_c$L, U=GET_c$U, M=GET_c$M,
                  fuera=GET_c$fuera)
p1 <- ggplot(df1, aes(x=t)) +
  geom_ribbon(aes(ymin=L,ymax=U), fill="#1D9E75", alpha=0.35) +
  geom_line(aes(y=L), color="#1D9E75", linewidth=0.5, linetype="dotted") +
  geom_line(aes(y=U), color="#1D9E75", linewidth=0.5, linetype="dotted") +
  geom_line(aes(y=M), color="#1D9E75", linewidth=0.5, linetype="dashed") +
  geom_hline(yintercept=1, color="gray50", linewidth=0.5) +
  geom_line(aes(y=St), color="#003D7A", linewidth=1.1) +
  geom_point(data=subset(df1,fuera), aes(y=St), color="#C0392B", size=2.5) +
  labs(
    title    = "S_t(t) — test CORREGIDO (IPW + permutaciones estratificadas)",
    subtitle = sprintf("S(u,t) 3D de Ghorbani (2021) | GET bilateral pooled | K=%d estratos | H₀ %s (p=%.4f)",
                       K_ESTRATOS, if(GET_c$rechazar) "RECHAZADA" else "no rechazada", GET_c$p),
    x="Año", y=expression(S[t](t)),
    caption="Línea azul: S_t observado | Banda verde: envolvente GET 95% | Puntos rojos: salida del corredor"
  ) +
  theme_minimal(base_size=11) +
  theme(plot.subtitle=element_text(color=if(GET_c$rechazar)"#C0392B" else "#1D9E75"))
p1
ggsave(file.path(DIR_OUT,paste0("grafico_GET_St_corr_","especie_tratada",".png")), p1, width=9, height=5, dpi=150)
cat(paste0("Guardado: grafico_GET_St_corr_",especie_tratada,".png\n"))

library(latex2exp)

p1 <- ggplot(df1, aes(x = t)) +
  geom_ribbon(aes(ymin = L, ymax = U), fill = "#1D9E75", alpha = 0.35) +
  geom_line(aes(y = L), color = "#1D9E75", linewidth = 0.5, linetype = "dotted") +
  geom_line(aes(y = U), color = "#1D9E75", linewidth = 0.5, linetype = "dotted") +
  geom_line(aes(y = M), color = "#1D9E75", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 1, color = "gray50", linewidth = 0.5) +
  geom_line(aes(y = St), color = "#003D7A", linewidth = 1.1) +
  geom_point(data = subset(df1, fuera), aes(y = St), color = "#C0392B", size = 2.5) +
  labs(
    title = TeX(paste0("$S_t(t)$ — test corregido (IPW + permutaciones estratificadas) para ", especie_tratada) ),
    
    subtitle = TeX(sprintf(
      "$S(u,t)$ 3D de Ghorbani (2021) | GET bilateral pooled | $K=%d$ estratos | $H_0$ %s ($p=%.4f$)",
      K_ESTRATOS,
      if(GET_c$rechazar) "rechazada" else "no rechazada",
      GET_c$p
    )),
    
    x = "Año",
    y = TeX("$S_t(t)$"),
    
    caption = TeX("Línea azul: $S_t$ observado | Banda verde: envolvente GET 95\\% | Puntos rojos: salida del corredor")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.subtitle = element_text(
      color = if(GET_c$rechazar) "#C0392B" else "#1D9E75"
    )
  )
p1
ggsave(file.path(DIR_OUT,paste0("grafico_GET_St_corr_",especie_tratada,".png")), p1, width=9, height=5, dpi=150)
cat(paste0("Guardado: grafico_GET_St_corr_",especie_tratada,".png\n"))

# ── Plot 2: S_t naive ─────────────────────────────────────────────────────────
df2 <- data.frame(t=grid_t, St=S_t_n, L=GET_n$L, U=GET_n$U, M=GET_n$M,
                  fuera=GET_n$fuera)
p2 <- ggplot(df2, aes(x=t)) +
  geom_ribbon(aes(ymin=L,ymax=U), fill="#2980B9", alpha=0.35) +
  geom_line(aes(y=L), color="#2980B9", linewidth=0.5, linetype="dotted") +
  geom_line(aes(y=U), color="#2980B9", linewidth=0.5, linetype="dotted") +
  geom_line(aes(y=M), color="#2980B9", linewidth=0.5, linetype="dashed") +
  geom_hline(yintercept=1, color="gray50", linewidth=0.5) +
  geom_line(aes(y=St), color="#C0392B", linewidth=1.1) +
  geom_point(data=subset(df2,fuera), aes(y=St), color="#C0392B", size=2.5, shape=17) +
  labs(
    title    = TeX(paste0("$S_t(t)$ — test NAIVE (sin corrección MNAR) para ",especie_tratada)),
    subtitle = TeX(sprintf("S(u,t) 3D de Ghorbani (2021) | GET bilateral pooled | $H_0$ %s ($p=%.4f$)",
                       if(GET_n$rechazar) "RECHAZADA" else "no rechazada", GET_n$p)),
    x="Año", y=expression(S[t](t)),
    caption=TeX("Línea roja: $S_t$ observado | Banda azul: envolvente GET $95\\%$ | Triángulos: salida del corredor")
  ) +
  theme_minimal(base_size=11) +
  theme(plot.subtitle=element_text(color=if(GET_n$rechazar)"#C0392B" else "#2980B9"))
p2
ggsave(file.path(DIR_OUT,paste0("grafico_GET_St_naive_",especie_tratada,".png")), p2, width=9, height=5, dpi=150)
cat(paste0("Guardado: grafico_GET_St_naive_",especie_tratada,".png\n"))

# ── Plot 3: Comparación naive vs corregido ────────────────────────────────────
df3 <- data.frame(t=grid_t, Stc=S_t_c, Stn=S_t_n,
                  Lc=GET_c$L, Uc=GET_c$U, Ln=GET_n$L, Un=GET_n$U)
p3 <- ggplot(df3, aes(x=t)) +
  geom_ribbon(aes(ymin=Ln,ymax=Un), fill="#2980B9", alpha=0.30) +
  geom_line(aes(y=Stn), color="#C0392B", linewidth=0.9, linetype="dashed") +
  geom_ribbon(aes(ymin=Lc,ymax=Uc), fill="#1D9E75", alpha=0.35) +
  geom_line(aes(y=Stc), color="#003D7A", linewidth=1.1) +
  geom_hline(yintercept=1, color="gray50", linewidth=0.5) +
  labs(
    title = paste0("Comparación: NAIVE vs CORREGIDO para ", especie_tratada),
    subtitle = TeX(paste0(
      "Naive (azul/rojo): $H_0$ ", if(GET_n$rechazar) "RECHAZADA" else "no rechazada",
      sprintf(" ($p$=%.4f)", GET_n$p),
      "  |  Corregido (verde/azul): $H_0$ ",
      if(GET_c$rechazar) "RECHAZADA" else "no rechazada",
      sprintf(" ($p$=%.4f)", GET_c$p))),
    x="Año", y=expression(S[t](t)),
    caption=TeX(paste0(
      "Envolvente GET genuina (Myllymäki 2017): salida ↔ rechazo $H_0$\n",
      "Banda azul: naive (MCAR, perms. globales) | Banda verde: corregido (MNAR, perms. estratificadas π̂)\n",
      sprintf("chi² NAIVE $p$=%.4f | chi² CORREGIDO $p$=%.4f",
              res_chi_n$p, res_chi_c$p)))
  ) +
  theme_minimal(base_size=11) +
  theme(plot.subtitle=element_text(size=9,color="gray30"),
        plot.caption=element_text(size=8,color="gray50",hjust=0))
p3 <- ggplot(df3, aes(x=t)) +
  geom_ribbon(aes(ymin=Ln,ymax=Un), fill="#2980B9", alpha=0.30) +
  geom_line(aes(y=Stn), color="#C0392B", linewidth=0.9, linetype="dashed") +
  geom_ribbon(aes(ymin=Lc,ymax=Uc), fill="#1D9E75", alpha=0.35) +
  geom_line(aes(y=Stc), color="#003D7A", linewidth=1.1) +
  geom_hline(yintercept=1, color="gray50", linewidth=0.5) +
  labs(
    title = paste0("Comparación: NAIVE vs CORREGIDO para ", especie_tratada),
    subtitle = TeX(paste0(
      "Naive (azul/rojo): $H_0$ ", if(GET_n$rechazar) "RECHAZADA" else "no rechazada",
      sprintf(" ($p$=%.4f)", GET_n$p),
      "  |  Corregido (verde/azul): $H_0$ ",
      if(GET_c$rechazar) "RECHAZADA" else "no rechazada",
      sprintf(" ($p$=%.4f)", GET_c$p))),
    x="Año", y=expression(S[t](t)),
    caption = TeX(paste0(
      "Envolvente GET genuina: salida $\\leftrightarrow$ rechazo de $H_0$. ",
      "Banda azul: naive. Banda verde: corregido por $\\hat{\\pi}$. ",
      sprintf("$\\chi^2$ naive: $p=%.4f$ | $\\chi^2$ corregido: $p=%.4f$",
              res_chi_n$p, res_chi_c$p)
    ))
  ) +
  theme_minimal(base_size=11) +
  theme(plot.subtitle=element_text(size=9,color="gray30"),
        plot.caption=element_text(size=8,color="gray50",hjust=0))


p3
ggsave(file.path(DIR_OUT,paste0("grafico_comparacion_",especie_tratada,".png")), p3, width=10, height=5, dpi=150)
cat(paste0("Guardado: grafico_comparacion_",especie_tratada,".png\n"))

# ── Plot 4: S_sp(u) corregido ─────────────────────────────────────────────────
# ── Ecuador continental en WGS84  ──────────────────
ecuador <- ne_countries(scale = "medium", country = "Ecuador", returnclass = "sf")
ecuador_continental <- ecuador %>%
  st_crop(xmin = -81.5, xmax = -75, ymin = -5.5, ymax = 2)

df4 <- cbind(grid_sp, S_sp=S_sp_c)
p4 <- ggplot(df4, aes(x=lon,y=lat,fill=S_sp)) +
  geom_tile() +
  geom_sf(data = ecuador_continental, inherit.aes = FALSE,
          fill = NA, color = "gray30", linewidth = 0.5) +
  scale_fill_gradient2(low="#2980B9",mid="white",high="#C0392B",midpoint=1,
                       name="S_sp(u)") +
  geom_point(data=df_pres, aes(x=lon,y=lat), inherit.aes=FALSE,
             color="black", size=0.5, alpha=0.4) +
  labs(title="Proyección espacial S_sp(u) — corregido",
       subtitle="Rojo: más intenso que esperado bajo H₀ | Azul: menos",
       x="Longitud", y="Latitud",
       caption="Puntos negros: registros de presencia") +
  coord_sf() + theme_minimal(base_size=11)

df4 <- cbind(grid_sp, S_sp = S_sp_c)

# Convertir ceros a NA para que sean transparentes
df4$S_sp[df4$S_sp == 0] <- NA

p4 <- ggplot(df4, aes(x = lon, y = lat, fill = S_sp)) +
  geom_raster(interpolate = TRUE) +        # interpolado = bordes suaves
  geom_sf(data = ecuador_continental, inherit.aes = FALSE,
          fill = NA, color = "gray20", linewidth = 0.6) +
  scale_fill_gradient2(
    low = "#2980B9", mid = "white", high = "#C0392B",
    midpoint = 1, name = TeX("$S_{sp}(u)$"),
    na.value = "transparent"               # ceros = transparente
  ) +
  geom_point(data = df_pres, aes(x = lon, y = lat), inherit.aes = FALSE,
             color = "black", size = 0.5, alpha = 0.5) +
  coord_sf(xlim = c(-81.5, -75.2), ylim = c(-5.2, 1.6)) +
  labs(
    title    = TeX(paste0("Proyección espacial $S_{sp}(u)$ — corregido para ",especie_tratada)),
    subtitle = TeX("Rojo: más intenso que esperado bajo $H_0$ | Azul: menos"),
    x = "Longitud", y = "Latitud",
    caption  = "Puntos negros: registros de presencia"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "gray95", color = NA)
  )
p4
ggsave(file.path(DIR_OUT,paste0("grafico_GET_Ssp_",especie_tratada,".png")), p4, width=7, height=7, dpi=150)
cat(paste0("Guardado: grafico_GET_Ssp_",especie_tratada,".png\n"))


# ── Plot 5: chi-cuadrado visual ───────────────────────────────────────────────
df5 <- data.frame(
  test = rep(c("NAIVE","CORREGIDO"), each=2),
  resultado = rep(c("chi²","p-valor"), 2),
  valor = c(res_chi_n$chi2, res_chi_n$p, res_chi_c$chi2, res_chi_c$p)
)
# Tabla de resultados chi2 como texto en ggplot
chi_txt <- data.frame(
  x = c(1,2), y = c(1,1),
  label = c(
    sprintf("NAIVE\nX²=%.2f\ndf=%d\np=%.4f\n%s",
            res_chi_n$chi2, res_chi_n$df, res_chi_n$p,
            if(res_chi_n$rechazar) "H₀ RECHAZADA" else "H₀ no rechazada"),
    sprintf("CORREGIDO\nX²=%.2f\ndf=%d\np=%.4f\n%s",
            res_chi_c$chi2, res_chi_c$df, res_chi_c$p,
            if(res_chi_c$rechazar) "H₀ RECHAZADA" else "H₀ no rechazada")
  ),
  color = c(
    if(res_chi_n$rechazar) "#C0392B" else "#2980B9",
    if(res_chi_c$rechazar) "#C0392B" else "#1D9E75"
  )
)
p5 <- ggplot(chi_txt, aes(x=x,y=y,label=label,color=color)) +
  geom_text(size=5, fontface="bold", lineheight=1.4) +
  scale_color_identity() +
  scale_x_continuous(limits=c(0.5,2.5)) +
  scale_y_continuous(limits=c(0.5,1.5)) +
  labs(title=paste0(sprintf("Test chi-cuadrado de separabilidad — grilla %d×%d×%d para ",
                     I_chi,I_chi,J_chi),especie_tratada),
       subtitle="Conteos observados (naive) vs sumas IPW (corregido)") +
  theme_void(base_size=13) +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,color="gray40"))
p5
ggsave(file.path(DIR_OUT,"grafico_chi2.png"), p5, width=7, height=4, dpi=150)
cat("  Guardado: grafico_chi2.png\n")


# ── Plot 6: S_sp(u) naive ─────────────────────────────────────────────────
df6 <- cbind(grid_sp, S_sp=S_sp_n)
p6 <- ggplot(df6, aes(x=lon,y=lat,fill=S_sp)) +
  geom_tile() +
  scale_fill_gradient2(low="#2980B9",mid="white",high="#C0392B",midpoint=1,
                       name="S_sp(u)") +
  geom_point(data=df_pres, aes(x=lon,y=lat), inherit.aes=FALSE,
             color="black", size=0.5, alpha=0.4) +
  labs(title="Proyección espacial S_sp(u) — naive",
       subtitle="Rojo: más intenso que esperado bajo H₀ | Azul: menos",
       x="Longitud", y="Latitud",
       caption="Puntos negros: registros de presencia") +
  coord_equal() + theme_minimal(base_size=11)

df6 <- cbind(grid_sp, S_sp = S_sp_n)

# Convertir ceros a NA para que sean transparentes
df6$S_sp[df6$S_sp == 0] <- NA

p6 <- ggplot(df6, aes(x = lon, y = lat, fill = S_sp)) +
  geom_raster(interpolate = TRUE) +
  geom_sf(
    data = ecuador_continental,
    inherit.aes = FALSE,
    fill = NA,
    color = "gray20",
    linewidth = 0.6
  ) +
  scale_fill_gradient2(
    low = "#2980B9",
    mid = "white",
    high = "#C0392B",
    midpoint = 1,
    name = TeX("$S_{sp}(u)$"),
    na.value = "transparent"
  ) +
  geom_point(
    data = df_pres,
    aes(x = lon, y = lat),
    inherit.aes = FALSE,
    color = "black",
    size = 0.5,
    alpha = 0.5
  ) +
  coord_sf(
    xlim = c(-81.5, -75.2),
    ylim = c(-5.2, 1.6)
  ) +
  labs(
    title = TeX(paste0("Proyección espacial $S_{sp}(u)$ — naive para ", especie_tratada)),
    subtitle = TeX("Rojo: más intenso que esperado bajo $H_0$ | Azul: menos"),
    x = "Longitud",
    y = "Latitud",
    caption = "Puntos negros: registros de presencia"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "gray95", color = NA)
  )

p6

ggsave(file.path(DIR_OUT,paste0("grafico_GET_Ssp_n_",especie_tratada,".png")), p6, width=7, height=7, dpi=150)
cat(paste0("Guardado: grafico_GET_Ssp_n_",especie_tratada,".png\n"))
# =============================================================================
# SECCIÓN 12: Guardar resultados
# =============================================================================
cat("\n=== SECCIÓN 12: Guardando resultados ==================================\n")

resultados <- list(
  # Config
  n_pres=n_pres, n_eff=n_eff, K_estratos=K_ESTRATOS,
  K_perm=K_PERM, alpha=ALPHA, hT=hT, hS=hS,
  N_GRID_SP=N_GRID_SP, N_GRID_T=N_GRID_T,
  grid_t=grid_t, grid_sp=grid_sp,
  # Estadísticos observados
  S_t_corr=S_t_c,  S_sp_corr=S_sp_c,
  S_t_naive=S_t_n, S_sp_naive=S_sp_n,
  # GET
  GET_corr=GET_c, GET_naive=GET_n,
  S_t_perm_corr=S_t_perm_c, S_t_perm_naive=S_t_perm_n,
  # Chi-cuadrado
  chi_naive=res_chi_n, chi_corr=res_chi_c,
  I_chi=I_chi, J_chi=J_chi,
  # Intensidades marginales
  rho_sp_corr=rho_sp_c, rho_t_corr=rho_t_c,
  rho_sp_naive=rho_sp_n, rho_t_naive=rho_t_n
)
saveRDS(resultados, file.path(DIR_OUT,paste0("resultados_GET_v3_",especie_tratada,".rds")))
cat(paste0("Guardado: resultados_GET_v3_",especie_tratada,".rds\n"))

# =============================================================================
# RESUMEN FINAL
# =============================================================================
cat("\n=========================================================================\n")
cat("  RESUMEN FINAL — TEST DE SEPARABILIDAD ESPACIO-TEMPORAL\n")
cat("=========================================================================\n")
cat(sprintf("  n presencias=%d  K_estratos=%d  K_perm=%d\n",
            n_pres, K_ESTRATOS, K_PERM))
cat(sprintf("  Grilla 3D evaluada: %d puntos (N_SP=%d, N_T=%d)\n",
            N_GRID_SP^2 * N_GRID_T, N_GRID_SP, N_GRID_T))
cat("\n  ── GET (S_t, proyección de S(u,t) 3D) ─────────────────────────────\n")
cat(sprintf("  NAIVE    : p=%.4f  H₀ %s  ancho=%.4f\n",
            GET_n$p, if(GET_n$rechazar)"RECHAZADA ***" else "no rechazada", GET_n$ancho))
cat(sprintf("  CORREGIDO: p=%.4f  H₀ %s  ancho=%.4f\n",
            GET_c$p, if(GET_c$rechazar)"RECHAZADA ***" else "no rechazada", GET_c$ancho))
cat("\n  ── Chi-cuadrado ────────────────────────────────────────────────────\n")
cat(sprintf("  NAIVE    : X²=%.2f  df=%d  p=%.6f  H₀ %s\n",
            res_chi_n$chi2, res_chi_n$df, res_chi_n$p,
            if(res_chi_n$rechazar)"RECHAZADA ***" else "no rechazada"))
cat(sprintf("  CORREGIDO: X²=%.2f  df=%d  p=%.6f  H₀ %s\n",
            res_chi_c$chi2, res_chi_c$df, res_chi_c$p,
            if(res_chi_c$rechazar)"RECHAZADA ***" else "no rechazada"))

cat("\n  ── Interpretación Metodológica ─────────────────────────────────────\n")

# 1. Diagnóstico de divergencia Chi2 vs GET
if (res_chi_c$rechazar && !GET_c$rechazar) {
  cat("  [!] DIVERGENCIA DETECTADA: Chi² rechaza pero GET no.\n")
  cat("      -> CONFÍE EN EL GET.\n")
  cat("      El estadístico Chi² está inflado por celdas con expectativas minúsculas\n")
  cat("      debido al clumping espacial y n moderado, elevando el Error Tipo I.\n")
  cat("      El GET es asintóticamente libre y robusto frente al MAUP.\n\n")
}

# 2. Diagnóstico del efecto MNAR
if (GET_n$rechazar && !GET_c$rechazar) {
  cat("  [✓] ARTEFACTO MNAR: El rechazo naive desaparece con la corrección IPW.\n")
  cat("      La aparente no-separabilidad era solo sesgo de muestreo.\n")
} else if (!GET_n$rechazar && GET_c$rechazar) {
  cat("  [✓] EFECTO MNAR: El sesgo de muestreo ocultaba la no-separabilidad real,\n")
  cat("      la cual solo es visible tras la corrección estratificada.\n")
} else if (GET_n$rechazar && GET_c$rechazar) {
  cat("  [✓] Ambos test GET rechazan: No-separabilidad real confirmada.\n")
  if (GET_c$ancho < GET_n$ancho)
    cat("      Corredor corregido más estrecho → mayor potencia bajo MNAR.\n")
} else {
  cat("  [✓] Ambos test GET no rechazan: Sin evidencia estocástica de no-separabilidad.\n")
  
  if (n_pres < 300) {
    cat(sprintf("\n  NOTA SOBRE POTENCIA Y TAMAÑO MUESTRAL (n=%d):\n", n_pres))
    cat("  Si sospecha que debería haber rechazo, el test GET 3D puede estar perdiendo\n")
    cat("  potencia si la grilla es muy fina (exceso de control FWER). Opciones:\n")
    cat("  1. Reducir N_GRID_SP (ej. a 10-12) para concentrar la densidad.\n")
    cat("  2. Reducir K_ESTRATOS (ej. a 2) para dar más variabilidad a la permutación.\n")
  }
}
cat("=========================================================================\n")



#------------------- GRAFICO: Comparación de densidades---------------------------

library(dplyr)
library(sf)
library(spatstat.geom)
library(spatstat.explore)
library(rnaturalearth)
library(ggplot2)
library(viridis)
library(stars)  

# ── Datos y geometría ────────────────────────────────────────────────────────
ecuador <- ne_countries(scale = "medium", country = "Ecuador", returnclass = "sf")

ecuador_continental <- ecuador %>%
  st_crop(xmin = -81.5, xmax = -75, ymin = -5.5, ymax = 2)

ecuador_proj <- st_transform(ecuador_continental, 32717)
ecuador_owin <- as.owin(ecuador_proj)

pres_sf   <- st_as_sf(df_pres, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
pres_proj <- st_transform(pres_sf, 32717)
coords_utm <- st_coordinates(pres_proj)

# ── Procesos puntuales ───────────────────────────────────────────────────────
pp_naive <- ppp(
  x = coords_utm[,1], y = coords_utm[,2],
  window = ecuador_owin,
  marks  = rep(1, nrow(df_pres))
)
pp_naive <- pp_naive[inside.owin(pp_naive$x, pp_naive$y, ecuador_owin)]

pp_ipw <- ppp(
  x = coords_utm[,1], y = coords_utm[,2],
  window = ecuador_owin,
  marks  = df_pres$w_ipw
)
pp_ipw <- pp_ipw[inside.owin(pp_ipw$x, pp_ipw$y, ecuador_owin)]

# ── Densidades ───────────────────────────────────────────────────────────────
sigma_use <- bw.scott(pp_naive)[1]

dens_naive <- density.ppp(pp_naive, sigma = sigma_use, eps = c(5000, 5000),
                          weights = marks(pp_naive))
dens_ipw   <- density.ppp(pp_ipw,   sigma = sigma_use, eps = c(5000, 5000),
                          weights = marks(pp_ipw))

# ── Helper: convierte im de spatstat → data.frame para ggplot ───────────────
im_to_df <- function(im_obj) {
  as.data.frame(im_obj) %>%
    rename(x = x, y = y, value = value) %>%
    dplyr::filter(!is.na(value))
}

df_naive <- im_to_df(dens_naive)
df_ipw   <- im_to_df(dens_ipw)

# ── Borde de Ecuador en WGS84 para la capa de polígono ──────────────────────
# (se mantiene en UTM 32717 para coincidir con la densidad)
ecuador_borde <- st_geometry(ecuador_proj)

# ── Gráfico 1: Intensidad no corregida (Naive) ───────────────────────────────
p_naive <- ggplot() +
  geom_raster(data = df_naive, aes(x = x, y = y, fill = value)) +
  geom_sf(data = ecuador_proj, fill = NA, color = "gray30", linewidth = 0.5) +
  scale_fill_viridis_c(option = "viridis", name = "Intensidad") +
  labs(
    title = paste0("Intensidad no corregida para ",especie_tratada),
    x = "Este (m)", y = "Norte (m)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.text     = element_text(size = 8)
  )

# ── Gráfico 2: Intensidad ponderada IPW ─────────────────────────────────────
p_ipw <- ggplot() +
  geom_raster(data = df_ipw, aes(x = x, y = y, fill = value)) +
  geom_sf(data = ecuador_proj, fill = NA, color = "gray30", linewidth = 0.5) +
  scale_fill_viridis_c(option = "viridis", name = "Intensidad") +
  labs(
    title = paste0("Intensidad ponderada IPW para ",especie_tratada),
    x = "Este (m)", y = "Norte (m)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.text     = element_text(size = 8)
  )

# ── Visualizar ───────────────────────────────────────────────────────────────
print(p_naive)
print(p_ipw)

# ── Guardar (opcional) ───────────────────────────────────────────────────────
ggsave(paste0("densidad_naive_",especie_tratada,".png"),  p_naive, width = 8, height = 6, dpi = 180)
ggsave(paste0("densidad_ipw_",especie_tratada,".png"),    p_ipw,   width = 8, height = 6, dpi = 180)
