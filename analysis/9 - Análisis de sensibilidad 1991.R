packages = c('readxl', 'mapSpain', 'tsModel', 'grid', 'tsibble', 'plot3D', 'splines', 'MASS', 'tidyverse',  'splines', 'gnm', 'dlnm', 'mixmeta', 'furrr', 'lubridate', 'ggplot2', 'gridExtra', 'mgcv')
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos='http://cran.rediris.es')
    library(x, character.only = TRUE)
  }
})


source("./data/funciones.R")
load("./data/serie.global.RData")

# Asigno un total de 0 suicidios a todos los días de 1991

datos.serie <- datos.serie %>%
  mutate(Total = ifelse(yy == 1991, 0, Total))

# Defino la variable de periodo vacacional

datos.serie <- datos.serie %>%
  mutate(vacaciones = ifelse(mm == 12 & dd == 25, "Navidad",
                             ifelse(mm == 1 & dd == 1, "Año nuevo",
                                    ifelse(mm == 1 & dd == 6|mm == 5 & dd == 1|mm == 9 & dd == 8|mm == 10 & dd == 12|mm == 11 & dd == 1|mm == 12 & dd == 6|mm == 12 & dd == 8, "Otras vacaciones", ifelse(mm == 7 & dd >= 15|mm == 8 & dd <= 15, "Verano", "No vacaciones"))))) 

datos.serie[c(106:107, 456:457, 813:814, 1198:1199, 1548:1549, 1933:1934, 2290:2291, 2647:2648, 3025:3026, 3382:3383, 3739:3740, 4117:4118, 4474:4475, 4859:4860, 5216:5217, 5566:5567, 5951:5952, 6308:6309, 6658:6659, 7043:7044, 7400:7401, 7750:7751, 8135:8136, 8492:8493, 8877:8878, 9227:9228, 9584:9585, 9969:9970, 10319:10320, 10676:10677, 11061:11062, 11411:11412, 11796:11797, 12153:12154, 12510:12511), "vacaciones"] <- "Jueves/viernes santo"

datos.serie$vacaciones <- as.factor(datos.serie$vacaciones)
datos.serie$vacaciones <- relevel(datos.serie$vacaciones, ref = "No vacaciones")

################################################################
# PREPARACIÓN PARA LA SELECCIÓN DEL MODELO

# Selecciono el modelo que mejor se ajusta a los datos según
# el criterio QAIC. Empiezo por una lista de modelos con diferentes
# características, se evalúa el QAIC de cada uno y al final nos
# quedamos con el modelo que presenta el QAIC más bajo.
################################################################

## Defino todas las posibles características del modelo 

# Número de retardos para la temperatura

nlag_vec <- 2:6

# Posibles tipos de relación dosis-respuesta para la temperatura

x_funlist <- list("poly", "ns")

# Posibles tipos de relación lag-respuesta para las vacaciones

x_funlistv <- list("integer", "strata") 

# Los posibles grados de los polinomios o ubicación de los nodos de los splines para la relación dosis-respuesta de la temperatura

x_dflist <- list(
  degree = 1,
  degree = 2,
  degree = 3,
  knots = c(10, 90),
  knots = c(10, 50),
  knots = c(50, 90),
  knots = c(25, 50, 90)
)

# Los posibles intervalos de retardos a considerar en el caso de una relación lag-respuesta de tipo "strata" para las vacaciones

x_breakslist <- list(
  breaks = c(0, 1),
  breaks = c(-1, 0, 2),
  breaks = c(-1, 2),
  breaks = c(-1, 0, 1),
  breaks = c(0, 2)
)

# Posible número de nodos a considerar para la relación lag-respuesta de la temperatura

lagnk_vec <- 1:2 

# Todas las posibles combinaciones

nproof <- expand.grid(
  funx = 1:length(x_funlist),
  dfx = 1:length(x_dflist),
  nlag = 1:length(nlag_vec),
  lagnk = 1:length(lagnk_vec),
  funxv = 1:length(x_funlistv),
  breaksx = 1:length(x_breakslist)
)

dim(nproof)[1] # número inicial de modelos a probar: 1400


#Elimino de nproof las opciones que no me interesa o no se pueden evaluar: un polinomio de grado mayor de 3, splines que no son cúbicos y los casos con demasiados nodos para poca ventana, ya que por ejemplo en menos de 4 días es poco probable que el efecto cambie más de 1 vez o que en menos de 7 días el efecto cambie más de 2 veces.

idno <- ((nproof$funx == 1 &
            nproof$dfx > 3) |
           (nproof$funx == 2 & nproof$dfx <= 3) |
           (nproof$nlag < 3 & nproof$lagnk > 1)
)
nproof <- nproof[idno == F, ]
dim(nproof)[1] # número final de modelos a probar: 560

# Genero una lista con las opciones válidas, en la que guardaré los resultados

proofs <- vector("list", nrow(nproof))
names(proofs) <- paste("opcion", 1:nrow(nproof), sep = "")
for (i in 1:nrow(nproof)) {
  proofs[[i]] <- vector("list", 6)
  names(proofs[[i]]) <- c("funx", "dfx", "nlag", "lagnk", "funxv", "breaksx")
  proofs[[i]][[1]] <- x_funlist[nproof[i, 1]]
  proofs[[i]][[2]] <- x_dflist[nproof[i, 2]]
  proofs[[i]][[3]] <- nlag_vec[nproof[i, 3]]
  proofs[[i]][[4]] <- lagnk_vec[nproof[i, 4]]
  proofs[[i]][[5]] <- x_funlistv[nproof[i, 5]]
  proofs[[i]][[6]] <- x_breakslist[nproof[i, 6]]
}

names(proofs) <- paste("proofs", 1:length(proofs), sep = "")

################################################################
# SELECCIÓN DEL MODELO CONDICIONAL

# Selecciono el modelo que mejor se ajusta a los datos 
# usando la función gnm y empleando el uso de crossbasis,
# por lo tanto un modelo condicional de retardos distribuidos.
################################################################ 

# Mediante una función ajusto cada uno de los posibles modelos y guardo el QAIC correspondiente. El modelo final será el que tiene el QAIC más bajo y lo guardo para no tener que realizar este proceso todas las veces.

datos.serie$stratum <- as.factor(as.factor(datos.serie$yy))
ind <- tapply(datos.serie$Total, datos.serie$stratum, sum)[datos.serie$stratum]

# 
# QAICsel <- function(x) {
#   dfj <- x$dfx
#   mi.argvarj <- if (x$funx == "ns") {
#     knotsx <- quantile(datos.serie$temp.pob.sum,
#                        probs = dfj[[1]] / 100,
#                        na.rm = T)
#     boundx <- quantile(datos.serie$temp.pob.sum,
#                        probs = c(0, 100) / 100,
#                        na.rm = T) #mínimo y max como nodos de contorno
#     list(
#       fun = "ns",
#       knots = knotsx,
#       Boundary.knots = boundx,
#       int = F
#     )
#   } else{
#     list(fun = "poly",
#          degree = dfj[[1]],
#          int = F)
#   }
#   
#   nlagj <- x$nlag
#   klagj <- logknots(nlagj, x$lagnk)
#   mi.arglagj <- list(fun = "ns", knots = klagj, int = T)
#   
#   cbTj <- crossbasis(
#     datos.serie$temp.pob.sum,
#     lag = nlagj,
#     argvar = mi.argvarj,
#     arglag = mi.arglagj
#   )
#   vacan <- as.factor(as.numeric(datos.serie$vacaciones))
#   
#   breaksj <- x$breaksx
#   mi.arglagjv <- if (x$funxv == "integer") {
#     list(fun = "integer", int = T)
#   } else{
#     list(fun = "strata", breaks = breaksj[[1]])
#   }
#   
#   cbvacaj <- crossbasis(
#     vacan,
#     lag = c(-2, 2),
#     argvar = list(fun = "integer", int = F),
#     arglag = mi.arglagjv
#   )
#   
#   modelgnmj <- try(gnm(
#     Total ~ cbTj + cbvacaj + dow,
#     data = datos.serie,
#     family = quasipoisson,
#     eliminate = factor(stratum),
#     na.action = "na.exclude",
#     subset = ind > 0
#   ))
#   cat(length(modelgnmj$coefficients))
#   QAICj <- try(QAICM(modelgnmj, "dev"))
#   lres <- list(
#     fit = modelgnmj,
#     QAIC = QAICj,
#     cbT = cbTj,
#     cbvaca = cbvacaj
#   )
#   return(lres)
# }
# 
# # identificación  de la opción con minimo QAIC.
# 
# fitlist <- map(proofs, QAICsel)
# QAICs <- map_df(fitlist, function(x) {
#   data.frame(QAIC = x$QAIC)
# })
# 
# proofs[which.min(QAICs$QAIC)]
# 
# 
# best <- fitlist[which.min(QAICs$QAIC)]
# 
# # rescate del modelo
# 
# model.gnm2 <- map_depth(best, 1, function(x) {
#   return(x$fit)
# }) # o similar (esto no lo he probado)
# 
# model.gnm2 <- model.gnm2$proofs93
# cbvacaj <- best$proofs93$cbvaca
# cbTj <- best$proofs93$cbT
# save(model.gnm2, cbvacaj, cbTj, file = "./data/modelo.sensibilidad91.gnm.RData")

################################################################
# EXPLORACIÓN DEL MODELO FINALMENTE SELECCIONADO

# Estudio el efecto de las variables consideradas obtenido 
# mediante el modelo finalmente seleccionado
################################################################ 

load("./data/modelo.sensibilidad91.gnm.RData")

# El modelo con menor QAIC tiene las siguientes características:
# temperatura: dosis-respuesta polinomio de grado 2, lag-respuesta natural cubic splines con 2 nodos en 4 retardos
# vacaciones: lag-respuesta con strata e intervalos -2/-1, 0 y 1/2

## DÍA DE LA SEMANA

summary(model.gnm2)

# Represento el efecto del día de la semana

yall <- coef(model.gnm2)[24:29]

# Sall.gnm <- confint(model.gnm2)
# save(Sall.gnm, file = "./data/ICDiaGnmSensibilidad91.RData")
load("./data/ICDiaGnmSensibilidad91.RData")
Sall2 <- Sall.gnm[24:29, 1:2]

eyall <- exp(yall)
esall <- exp(Sall2)

p.dia <- data.frame(RR = eyall,
                    LowRR = esall[, 1],
                    HighRR = esall[, 2])

p.jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

p.dia <- rbind(p.dia[2:4, ], p.jueves, p.dia[6, ], p.dia[5, ], p.dia[1, ])

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

vpredcen <- crosspred(cbvacaj, model.gnm2, at = 1:6, cen = 1)

vpredcen$matRRlow #si mayor que 1 es significativo (riesgo)
vpredcen$matRRfit
vpredcen$matRRhigh #si menor que 1 es significativo (protector)

# Año Nuevo

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
  labs(x = "\nLag", y = "RR\n", title = "Año Nuevo") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

# Jueves/Viernes Santo

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
  labs(x = "\nLag", y = "RR\n", title = "Jueves/Viernes Santo") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

# Navidad

p <- data.frame(
  RR = vpredcen$matRRfit[4, ],
  LowRR = vpredcen$matRRlow[4, ],
  HighRR = vpredcen$matRRhigh[4, ]
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
  labs(x = "\nLag", y = "RR\n", title = "Navidad") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

# Verano

p <- data.frame(
  RR = vpredcen$matRRfit[6, ],
  LowRR = vpredcen$matRRlow[6, ],
  HighRR = vpredcen$matRRhigh[6, ]
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
  guides(colour = guide_legend(override.aes = list(shape = NA)))

# Otras vacaciones

p <- data.frame(
  RR = vpredcen$matRRfit[5, ],
  LowRR = vpredcen$matRRlow[5, ],
  HighRR = vpredcen$matRRhigh[5, ]
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
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## TEMPERATURA

# malla para la predicción

tpred <- seq(0, 26, by = 1)

pred <- crosspred(cbTj, model.gnm2, at = tpred)

# Estimo la temperatura de mínima mortalidad y su IC95%

mmt <- pred$predvar[which.min(pred$allRRfit)]

pci2 <- quantile(findmin(
  cbTj,
  model.gnm2,
  from = 0,
  to = 26,
  by = 1,
  sim = T,
  nsim = 10000
),
c(2.5, 97.5) / 100)

# Predicción de efectos respecto de TMM

predcen <- crosspred(cbTj, model.gnm2, at = tpred, cen = mmt)

# Efecto global

predcen$allRRlow
predcen$allRRfit
predcen$allRRhigh

# Efectos retardados

predcen$matRRlow
predcen$matRRfit
predcen$matRRhigh

# Represento el efecto de la temperatura

RR_overall <- data.frame(
  RR = predcen$allRRfit,
  LowRR = predcen$allRRlow,
  HighRR = predcen$allRRhigh
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

# Miro el efecto en el punto de calor extremo (percentil 99.5)

calor <- quantile(datos.serie$temp.pob.sum,
                  probs = c(0.995),
                  na.rm = T)

predcenP95 <- crosspred(cbTj, model.gnm2, at = calor, cen = mmt)

plot(predcenP95,
     ptype = "slices",
     ci = "bars",
     var = calor)

c(predcenP95$allRRlow,
  predcenP95$allRRfit,
  predcenP95$allRRhigh)

predcenP95$matRRlow
predcenP95$matRRfit
predcenP95$matRRhigh

################################################################
# VALIDACIÓN DEL MODELO FINALMENTE SELECCIONADO

# Validación del modelo mediante el estudio de los residuos
# y de la multicolinealidad
################################################################

# Estudio de la correlación de los residuos mediante función de autocorrelación parcial

res <- residuals(model.gnm2)
pacf(as.numeric(na.omit(res)), lag.max = 4)
