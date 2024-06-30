packages = c('readxl', 'mapSpain', 'tsModel', 'grid', 'tsibble', 'plot3D', 'splines', 'MASS', 'tidyverse',  'splines', 'gnm', 'dlnm', 'mixmeta', 'furrr', 'lubridate', 'ggplot2', 'gridExtra', 'mgcv', 'data.table', 'INLA', 'sf', 'sp', 'spdep', 'ggpubr', 'hydroGOF', 'geofacet', 'ggthemes')
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos='http://cran.rediris.es')
    library(x, character.only = TRUE)
  }
})

source("./data/funciones.R")

################################################################
# PREPARACIÓN DE LOS DATOS

# Preparo los datos para poder ajustar el modelo en INLA
################################################################

# Exporto los datos de la cartografía de Asturias

Asturias <- esp_get_munic_siane(region = "Asturias") %>%
  mutate(Provincia = esp_dict_translate(ine.prov.name, "es"))

# Guardo los datos de los municipios para crear la estructura de vecindad

shp_AS <- Asturias %>%
  mutate(AS = as.factor(name)) %>%
  group_by(name, cmun) %>%  summarize(geometry = st_union(geometry))

shp_AS <- shp_AS %>%
  mutate(cmun = str_remove(cmun, "^0+"))

shp_AS <- shp_AS %>%
  arrange(as.numeric(cmun))

p <- st_as_sf(shp_AS)

nb.map <- poly2nb(as_Spatial(p$geometry))

# Creo el grafo con los vecinos que hará falta a la hora de ajustar
# el modelo de INLA

nb2INLA("Asturias.graph", nb.map)

H <- inla.read.graph(filename = "Asturias.graph")

# Represento matriz con los vecinos

image(inla.graph2matrix(H), xlab = "", ylab = "") #cada cuadradito es un municipio y los de la misma fila son sus vecinos

# Represento el mapa de Asturias con el grafo

plot(st_geometry(Asturias),
     axes = T,
     col = "darkgreen")

plot(
  nb.map,
  Asturias$geometry,
  add = T,
  lwd = 2,
  cex = 2
)

# Represento el mapa de Asturias con sus municipios
# puedo seleccionar un municipio y R me dice
# quiénes son sus vecinos

plot_map_neig <- function(neig, p) {
  plot(st_geometry(p), col = "white")
  plot(st_geometry(p[neig, ]), col = "red", add = TRUE)
  neighbors <- st_geometry(p[nb.map[[neig]], ])
  plot(neighbors, col = "pink", add = TRUE)
  cat("Has seleccionado", p$name[neig], "y sus vecinos son:", "\n")
  cat(p$name[nb.map[[neig]]], "\n")
}

# Ejemplo:

plot_map_neig(neig = 44, p)

# Importo la serie según el municipio

load("./data/serie.municipios.RData")

# Para cada serie defino la variable de periodo vacacional

ld <- split(datos.serie.mun, datos.serie.mun$MUNIC_DE)
names(ld) <- unique(datos.serie.mun$MUNIC_DE)


for(i in 1:78){
  ld[[i]] <- ld[[i]] %>%
    mutate(vacaciones = ifelse(mm == 12 & dd == 25, "Otras vacaciones",
                               ifelse(mm == 1 & dd == 1, "Otras vacaciones",
                                      ifelse(mm == 1 & dd == 6|mm == 5 & dd == 1|mm == 9 & dd == 8|mm == 10 & dd == 12|mm == 11 & dd == 1|mm == 12 & dd == 6|mm == 12 & dd == 8, "Otras vacaciones", ifelse(mm == 7 & dd >= 15|mm == 8 & dd <= 15, "Verano", "No vacaciones")))))
  
  ld[[i]][c(106:107, 456:457, 813:814, 1198:1199, 1548:1549, 1933:1934, 2290:2291, 2647:2648, 3025:3026, 3382:3383, 3739:3740, 4117:4118, 4474:4475, 4859:4860, 5216:5217, 5566:5567, 5951:5952, 6308:6309, 6658:6659, 7043:7044, 7400:7401, 7750:7751, 8135:8136, 8492:8493, 8877:8878, 9227:9228, 9584:9585, 9969:9970, 10319:10320, 10676:10677, 11061:11062, 11411:11412, 11796:11797, 12153:12154, 12510:12511), "vacaciones"] <- "Otras vacaciones"
  
  ld[[i]]$vacaciones <- as.factor(ld[[i]]$vacaciones)
  ld[[i]]$vacaciones <- relevel(ld[[i]]$vacaciones, ref = "No vacaciones")
}

# Uno las 78 series, al final tengo una base de datos con una serie temporal por cada municipio

datos.serie.mun <- list_rbind(ld)

# número de retardos que quiero considerar para temperatura

nlag = 4

# crossbasis para temperatura media

lag_tmed <- tsModel::Lag(datos.serie.mun$tmed,
                         group = datos.serie.mun$MUNIC_DE,
                         k = 0:nlag) #para cada municipio obtengo la temperatura en los diferentes retardos

# nodos para la relación lag-respuesta

lagknot = logknots(4, 2) #2 nodos en 4 retardos equiespaciados en escala logaritmica

var <- lag_tmed

#creo la crossbasis, para la relación dosis-respuesta uso un polinomio de grado 2, para la lag-respuesta un natural spline cúbico con 2 nodos en 4 retardos

basis_tmed <- crossbasis(
  var,
  argvar = list(fun = "poly", degree = 2, int = F),
  arglag = list(fun = "ns", knots = lagknot, int = T)
)

colnames(basis_tmed) = paste0("basis_tmed.", colnames(basis_tmed))

# para cada municipio miro los retardos 1 y 2 días antes y después de un determinado día vacacional.

lag_vacas <- tsModel::Lag(datos.serie.mun$vacaciones,
                          group = datos.serie.mun$MUNIC_DE,
                          k = -2:2)
var.vacas <- lag_vacas

#creo la crossbasis para las vacaciones, para la relación dosis-respuesta uso una variable indicadora para cada tipo de vacaciones (menos la referencia), mientras que para la lag-respuesta creo intervalos de retardos y asumo un efecto constante dentro de cada intervalo, en este caso un efecto constante 1 y 2 días antes del día considerado, uno en el mismo día y uno 1 y 2 días después del día considerado.

basis_vacas <- crossbasis(
  var.vacas,
  argvar = list(fun = "integer", int = F),
  arglag = list(fun = "strata", breaks = c(0, 1))
)

colnames(basis_vacas) = paste0("basis_vacas.", colnames(basis_vacas))

# Creo una variable índice para cada año, me hará falta después para controlar la tendencia en el modelo

datos.serie.mun$year_index <- datos.serie.mun$yy - 1986

# Creo una dummy para cada día de la semana (menos el día de referencia)

datos.serie.mun <- datos.serie.mun %>%
  mutate(
    lunes = ifelse(dow == "lunes", 1, 0),
    martes = ifelse(dow == "martes", 1, 0),
    miercoles = ifelse(dow == "miércoles", 1, 0),
    viernes = ifelse(dow == "viernes", 1, 0),
    sabado = ifelse(dow == "sábado", 1, 0),
    domingo = ifelse(dow == "domingo", 1, 0)
  )

################################################################
# DEFINICIÓN DEL MODELO 

# Defino el modelo que me permite estimar los efectos
# globales ajustando por los efectos espaciales de los municipios
################################################################

S <- U <- rep(1:78, each = 12784) # variables que representan el efecto espacial y heterogeneo en el modelo BYM

# Formula

formula <- Y ~ 1 + f(
  S,
  model       = "besag",
  graph       = H,
  scale.model = TRUE,
  replicate = T2,
  hyper       =
    list(prec = list(
      prior = "loggamma", param = c(0.01, 0.01)
    ))
) +
  f(U,
    model       = "iid",
    replicate = T2,
    hyper       =
      list(prec = list(
        prior = "loggamma", param = c(0.01, 0.01)
      ))) + basis_tmed + basis_vacas + lunes + martes + miercoles + viernes + sabado + domingo + T2

# Ajuste del modelo

# mod.suicides.nb <- inla(
#   formula,
#   family          = "nbinomial",
#   data            = datos.serie.mun,
#   control.compute = list(
#     dic = TRUE,
#     waic = TRUE,
#     cpo = TRUE,
#     config = T
#   ),
#   control.predictor = list(
#     link = 1,
#     compute = TRUE,
#     cdf = c(log(1))
#   ),
#   control.inla = list(strategy = 'adaptive'),
#   control.fixed = list(
#     correlation.matrix = TRUE,
#     prec.intercept = 1,
#     prec = 1
#   ),
#   verbose = FALSE
# )
# 
# save(mod.suicides.nb, file = "./data/ModeloINLAConDOWNB.RData")
load("./data/ModeloINLAConDOWNB.RData")

################################################################
# EXPLORACIÓN DEL MODELO

# Estudio el efecto de las variables consideradas obtenido 
# mediante el modelo espacio-temporal ajustado
################################################################ 

# Extraigo los coeficientes y las matrices de varianza-covarianza

coef <- mod.suicides.nb$summary.fixed$mean
vcov <- mod.suicides.nb$misc$lincomb.derived.covariance.matrix

## DÍA DE LA SEMANA

summary(mod.suicides.nb)

diasum <- mod.suicides.nb$summary.fixed[16:21, c(1, 3, 5)]

yall <- diasum[, 1] # Guardo los coeficientes correspondientes

eyall <- exp(yall) # Convierto los coeficientes en Riesgos Relativos
esall <- exp(diasum[, 2:3]) # Obtengo el intervalo de credibilidad de los RR

# Represento el efecto del día de la semana

p.dia <- data.frame(RR = eyall,
                    LowRR = esall[, 1],
                    HighRR = esall[, 2])

p.jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

p.dia <- rbind(p.dia[1:3, ], p.jueves, p.dia[4:6, ])

p.dia <- p.dia %>%
  mutate(dia = c(
    "lunes",
    "martes",
    "miércoles",
    "jueves",
    "viernes",
    "sábado",
    "domingo"
  ))

ggplot(p.dia, aes(
  x = factor(dia, levels = c("lunes", "martes", "miércoles", "jueves", "viernes", "sábado", "domingo")),
  y = RR,
  group = 1
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    colour = "purple3",
    size = 0.8
  ) +
  geom_point(colour = "purple4", size = 3) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "\nDía de la semana", y = "RR\n", title = "") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## PERIODO VACACIONAL

# Obtengo el efecto general de las vacaciones

indt.vacas <- grep("basis_vacas", mod.suicides.nb$names.fixed)

vpredcen <- crosspred(
  basis_vacas,
  coef = coef[indt.vacas],
  vcov = vcov[indt.vacas, indt.vacas],
  at = 1:3,
  cen = 1,
  model.link = "log",
  lag = c(-2, 2)
)

vpredcen$matRRlow
vpredcen$matRRfit
vpredcen$matRRhigh

# Represento el efecto del periodo vacacional

# Verano

p <- data.frame(
  RR = vpredcen$matRRfit[3, ],
  LowRR = vpredcen$matRRlow[3, ],
  HighRR = vpredcen$matRRhigh[3, ]
)

p <- p %>%
  mutate(lag = -2:2)

ggplot(p, aes(x = lag, y = RR)) +
  geom_errorbar(aes(ymin = LowRR, ymax = HighRR, width = 0.2),
                colour = "purple",
                size = 0.8) +
  geom_line(colour = "purple", size = 0.8) +
  geom_point(colour = "purple4", size = 3) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "\nLag", y = "RR\n", title = "Verano") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(ylim = c(0, 3))

# Otras vacaciones

p <- data.frame(
  RR = vpredcen$matRRfit[2, ],
  LowRR = vpredcen$matRRlow[2, ],
  HighRR = vpredcen$matRRhigh[2, ]
)

p <- p %>%
  mutate(lag = -2:2)

ggplot(p, aes(x = lag, y = RR)) +
  geom_errorbar(aes(ymin = LowRR, ymax = HighRR, width = 0.2),
                colour = "purple",
                size = 0.8) +
  geom_line(colour = "purple", size = 0.8) +
  geom_point(colour = "purple4", size = 3) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "\nLag", y = "RR\n", title = "Otras vacaciones") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA))) +
  coord_cartesian(ylim = c(0, 3))

## TEMPERATURA

# Me quedo sólo con los que tienen que ver con temperatura

indt <- grep("basis_tmed", mod.suicides.nb$names.fixed)

# Malla para las predicciones

tpred <- seq(0, 26, by = 1)

predt <- crosspred(
  basis_tmed,
  coef = coef[indt],
  vcov = vcov[indt, indt],
  model.link = "log",
  at = tpred
) 

# Estimo la temperatura de mínima mortalidad

mmt <- predt$predvar[which.min(predt$allRRfit)]

# Estimo el intrevalo de credibilidad para la TMM

p <- crossreduce(
  basis_tmed,
  coef = coef[indt],
  vcov = vcov[indt, indt],
  model.link = "log",
  at = tpred
)

xvar <- seq(0, 26, by = 1)
bvar <- do.call("onebasis", c(list(x = xvar), attr(basis_tmed, "argvar")))

pci <- quantile(
  findmin(
    basis = bvar,
    coef = p$coefficients,
    vcov = p$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

# Predicción de efectos respecto de TMM

predt <- crosspred(
  basis_tmed,
  coef = coef[indt],
  vcov = vcov[indt, indt],
  model.link = "log",
  at = tpred,
  cen = mmt
)

predt$allRRlow
predt$allRRhigh
predt$matRRlow

# Represento el efecto global de la temperatura media

RR_overall <- data.frame(
  RR = predt$allRRfit,
  LowRR = predt$allRRlow,
  HighRR = predt$allRRhigh
)

RR_overall <- rownames_to_column(RR_overall, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overall$tmean[RR_overall$RR == 1]

xlab <- pretty(RR_overall$tmean, n = 10)
ylab <- pretty(c(RR_overall$LowRR, RR_overall$HighRR), n = 10)

RR_overall %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "purple",
              alpha = 0.5) +
  geom_line(aes(tmean, RR), colour = "purple4", size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "purple4",
    size = 3,
    colour = "purple4",
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "\nTemperatura media (ºC)", y = "RR\n")
