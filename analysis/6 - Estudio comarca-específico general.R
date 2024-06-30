packages = c('tsModel', 'grid', 'tsibble', 'plot3D', 'splines', 'MASS',  'tidyverse', 'gnm', 'dlnm', 'mixmeta', 'furrr', 'lubridate', 'ggplot2',  'gridExtra', 'mgcv')
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos='http://cran.rediris.es')
    library(x, character.only = TRUE)
  }
})

options(scipen = 999)

# Cargo la serie temporal según la comarca

source("./data/funciones.R")
load("./data/serie.comarcas.RData")

# Creo una lista donde cada elemento es la serie de cada comarca

ld <- split(datos.serie.com, datos.serie.com$comarca_DE)
names(ld) <- unique(datos.serie.com$comarca_DE)

comarca <- c("aviles",
             "caudal",
             "eonavia",
             "gijon",
             "nalon",
             "narcea",
             "oriente",
             "oviedo")

## Efecto de la temperatura media

################################################################
# PRIMERA ETAPA

# Realizo la primera etapa del estudio comarca-específico.
# Consiste en ajustar un modelo para cada comarca y guardar
# los resultados relativos a la temperatura media
################################################################

# Creo el objeto donde guardaré los coeficientes del efecto de temperatura de cada comarca

yall <- matrix(NA, length(ld), 2, dimnames = list(comarca, paste("b", seq(2), sep = "")))

# Creo el objeto donde guardaré las matrices de varianza-covarianza para cada comarca

Sall <- vector("list", length(ld))
names(Sall) <- comarca

# Creo el objeto donde guardaré los coeficientes de la relación dosis-respuesta para la representación de los efectos de cada comarca

bvarj <- vector("list", 8)
names(bvarj) <- comarca

# Ajusto el modelo para cada comarca y guardo los resultados

system.time({
  for(i in seq(ld)) {
    sub <- ld[[i]]
    
    # Creo la variable de periodo vacacional
    suppressWarnings({
      sub <- sub %>%
        mutate(vacaciones = ifelse(mm == 12 & dd == 25, "Otras vacaciones",ifelse(mm == 1 & dd == 1, "Otras vacaciones", ifelse(mm == 1 & dd == 6|mm == 5 & dd == 1|mm == 9 & dd == 8|mm == 10 & dd == 12|mm == 11 & dd == 1|mm == 12 & dd == 6|mm == 12 & dd == 8, "Otras vacaciones", ifelse(mm == 7 & dd >= 15|mm == 8 & dd <= 15, "Verano", "No vacaciones")))))
      
      sub[c(106:107, 456:457, 813:814, 1198:1199, 1548:1549, 1933:1934, 2290:2291, 2647:2648, 3025:3026, 3382:3383, 3739:3740, 4117:4118, 4474:4475, 4859:4860, 5216:5217, 5566:5567, 5951:5952, 6308:6309, 6658:6659, 7043:7044, 7400:7401, 7750:7751, 8135:8136, 8492:8493, 8877:8878, 9227:9228, 9584:9585, 9969:9970, 10319:10320, 10676:10677, 11061:11062, 11411:11412, 11796:11797, 12153:12154, 12510:12511), "vacaciones"] <- "Otras vacaciones"
      
      sub$vacaciones <- as.factor(sub$vacaciones)
      sub$vacaciones <- relevel(sub$vacaciones, ref = "No vacaciones")
      
      # Creo el estrato para el modelo condicional
      sub$stratum <- as.factor(as.factor(sub$yy))
      ind <- tapply(sub$Total, sub$stratum, sum)[sub$stratum]
      
      # Defino la relación dosis-respuesta para temperatura
      mi.argvarj <- list(fun = "poly", degree = 2, int = F)
      
      # Defino el número de retardos
      nlagj <- 4
      
      # Defino la relación lag-respuesta para temperatura
      lagnk <- 2; klagj <- logknots(nlagj, lagnk)
      mi.arglagj <- list(fun = "ns", knots = klagj, int = T)
      
      # Defino la crossbasis
      cbTj <- crossbasis(sub$tmed, lag = nlagj, argvar = mi.argvarj, arglag = mi.arglagj)
      
      # Defino la relación lag-respuesta para vacaciones
      vacan <- as.factor(as.numeric(sub$vacaciones))
      mi.arglagjv <- list(fun = "strata", breaks = c(0, 1))
      
      # Defino la crossbasis para vacaciones
      cbvacaj <- crossbasis(vacan, lag = c(-2,2), argvar = list(fun= "integer", int = F), arglag = mi.arglagjv)
      
      dfseas <- 1
      ny <- length(unique(sub$yy))
    })
    
    # Ajusto el modelo
    mfirst <- gnm(Total ~ cbTj + cbvacaj + dow, data = sub, family = quasipoisson, eliminate = factor(stratum), na.action = "na.exclude", subset = ind > 0)
    
    # Malla para las predicciones
    tpred <- seq(0, 26, by = 1)
    
    bvarj[[i]] <- do.call("onebasis", c(list(x = tpred), attr(cbTj, "argvar")))
    
    pred <- crosspred(cbTj, mfirst, at = tpred)                 
    
    # Calculo el punto de mínima mortalidad
    mmt <- pred$predvar[which.min(pred$allRRfit)]

    # Guardo los coeficientes que resumen el efecto global de temperatura y que me harán falta para combinar los resultados
    suppressWarnings({
      crall <- crossreduce(cbTj, mfirst, cen = mmt, at = tpred)
    })
    
    yall[i,] <- coef(crall)  
    Sall[[i]] <- vcov(crall) 
    
  }
})

################################################################
# SEGUNDA ETAPA

# Realizo la segunda etapa del estudio comarca-específico.
# Consiste en combinar los resultados de cada comarca mediante
# meta-análisis y obtener un único efecto global.
################################################################

# Método de estimación de los coeficientes: reml = Restricted maximum likelihood (REML) estimator

method <- "reml"

# Ajusto el modelo meta-analítico

mvall <- mixmeta(yall ~ 1, Sall, method = method) 
summary(mvall)

# Creo la base para las predicciones

xvar <- seq(0, 26, by = 1)
bvar <- do.call("onebasis", c(list(x = xvar), attr(cbTj, "argvar")))

cpall <- crosspred(
  bvar,
  coef = coef(mvall),
  vcov = vcov(mvall),
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

# Calculo el punto de mínima mortalidad y su IC del 95%

mmt <- cpall$predvar[which.min(cpall$allRRfit)] 

pci <- quantile(
  findmin(
    basis = bvar,
    coef = coef(mvall),
    vcov = vcov(mvall),
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

# Obtengo el efecto global de temperatura media usando como referencia el tmm

cpall <- crosspred(
  bvar,
  coef = coef(mvall),
  vcov = vcov(mvall),
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# Represento el efecto global obtenido combinando los resultados de las comarcas

RR_overall <- data.frame(
  RR = cpall$allRRfit,
  LowRR = cpall$allRRlow,
  HighRR = cpall$allRRhigh
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
    size = 2,
    colour = "purple4",
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "\nTemperatura media (ºC)", y = "RR\n") +
  ggtitle("General")

# Límite inferior del 95%IC del RR, RR y límite superior del 95%IC

cpall$allRRlow
cpall$allRRfit
cpall$allRRhigh

################################################################
# BLUP

# Estimo las BLUPs de cada comarca en base al efecto estimado
# en cada comarca y al efecto global obtenido combinando los
# resultados
################################################################

# Estimación de las BLUPs

bl <- blup(mvall, vcov = T)

## Comarca de Avilés

# Sin corregir por el sesgo

cpalla1 <- crosspred(
  bvarj[[1]],
  coef = yall[1, ],
  vcov = Sall$aviles,
  model.link = "log",
  from = 0,
  to = 26,
  by = 1
)

mmt <- cpalla1$predvar[which.min(cpalla1$allRRfit)]

cpalla1 <- crosspred(
  bvarj[[1]],
  coef = yall[1, ],
  vcov = Sall$aviles,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# BLUP

cpalla2 <- crosspred(
  bvarj[[1]],
  coef = bl$aviles$blup,
  vcov = bl$aviles$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpalla2$predvar[which.min(cpalla2$allRRfit)]

pci <- quantile(
  findmin(
    basis = bvarj[[1]],
    coef = bl$aviles$blup,
    vcov = bl$aviles$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

cpalla2 <- crosspred(
  bvarj[[1]],
  coef = bl$aviles$blup,
  vcov = bl$aviles$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

cpalla2$allRRlow
cpalla2$allRRfit
cpalla2$allRRhigh

# Representación de la curva sin corregir y de la BLUP

RR_overalla <- data.frame(
  RR = cpalla1$allRRfit,
  LowRR = cpalla1$allRRlow,
  HighRR = cpalla1$allRRhigh
)

RR_overalla <- rownames_to_column(RR_overalla, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overalla$tmean[RR_overalla$RR == 1]

xlab <- pretty(RR_overalla$tmean, n = 10)
ylab <- pretty(c(RR_overalla$LowRR, RR_overalla$HighRR), n = 10)

RR_overalla2 <- data.frame(
  RR2 = cpalla2$allRRfit,
  LowRR2 = cpalla2$allRRlow,
  HighRR2 = cpalla2$allRRhigh
)

RR_overalla2 <- rownames_to_column(RR_overalla2, "tmean2") %>% mutate_at(1, as.numeric)

mmt2 <- RR_overalla2$tmean2[RR_overalla2$RR2 == 1]

xlab2 <- pretty(RR_overalla2$tmean2, n = 10)
ylab2 <- pretty(c(RR_overalla2$LowRR, RR_overalla2$HighRR), n = 10)

colors <- c("Comarca" = "olivedrab4",
            "BLUP" = "royalblue4")

RR_overalla %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "olivedrab3",
              alpha = 0.5) +
  geom_line(aes(tmean, RR, colour = "Comarca"), size = 1) +
  geom_ribbon(
    aes(ymin = RR_overalla2$LowRR2, ymax = RR_overalla2$HighRR2),
    fill = "royalblue3",
    alpha = 0.5
  ) +
  geom_line(aes(RR_overalla2$tmean2, RR_overalla2$RR2, colour = "BLUP"),
            size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "olivedrab4",
    size = 2,
    colour = "olivedrab4",
    show.legend = FALSE
  ) +
  geom_point(
    aes(mmt2, 1),
    shape = 21,
    fill = "royalblue4",
    size = 2,
    colour = "royalblue4",
    show.legend = T
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nTemperatura media (ºC)",
    y = "RR\n",
    title = "Avilés",
    color = ""
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca del Caudal

# Sin corregir por el sesgo

cpallc1 <- crosspred(
  bvarj[[2]],
  coef = yall[2, ],
  vcov = Sall$caudal,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallc1$predvar[which.min(cpallc1$allRRfit)]

cpallc1 <- crosspred(
  bvarj[[2]],
  coef = yall[2, ],
  vcov = Sall$caudal,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# BLUP

cpallc2 <- crosspred(
  bvarj[[2]],
  coef = bl$caudal$blup,
  vcov = bl$caudal$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallc2$predvar[which.min(cpallc2$allRRfit)]

pci <- quantile(
  findmin(
    basis = bvarj[[2]],
    coef = bl$caudal$blup,
    vcov = bl$caudal$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

cpallc2 <- crosspred(
  bvarj[[2]],
  coef = bl$caudal$blup,
  vcov = bl$caudal$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

cpallc2$allRRlow
cpallc2$allRRfit
cpallc2$allRRhigh

# Representación de la curva sin corregir y de la BLUP

RR_overallc <- data.frame(
  RR = cpallc1$allRRfit,
  LowRR = cpallc1$allRRlow,
  HighRR = cpallc1$allRRhigh
)

RR_overallc <- rownames_to_column(RR_overallc, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overallc$tmean[RR_overallc$RR == 1]

xlab <- pretty(RR_overallc$tmean, n = 10)

RR_overallc2 <- data.frame(
  RR2 = cpallc2$allRRfit,
  LowRR2 = cpallc2$allRRlow,
  HighRR2 = cpallc2$allRRhigh
)

RR_overallc2 <- rownames_to_column(RR_overallc2, "tmean2") %>% mutate_at(1, as.numeric)

mmt2 <- RR_overallc2$tmean2[RR_overallc2$RR2 == 1]

xlab2 <- pretty(RR_overallc2$tmean2, n = 10)
ylab2 <- pretty(c(RR_overallc2$LowRR, RR_overallc2$HighRR), n = 10)

colors <- c("Comarca" = "olivedrab4",
            "BLUP" = "royalblue4")

RR_overallc %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "olivedrab3",
              alpha = 0.5) +
  geom_line(aes(tmean, RR, colour = "Comarca"), size = 1) +
  geom_ribbon(
    aes(ymin = RR_overallc2$LowRR2, ymax = RR_overallc2$HighRR2),
    fill = "royalblue3",
    alpha = 0.5
  ) +
  geom_line(aes(RR_overallc2$tmean2, RR_overallc2$RR2, colour = "BLUP"),
            size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "olivedrab4",
    size = 2,
    colour = "olivedrab4",
    show.legend = FALSE
  ) +
  geom_point(
    aes(mmt2, 1),
    shape = 21,
    fill = "royalblue4",
    size = 2,
    colour = "royalblue4",
    show.legend = T
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nTemperatura media (ºC)",
    y = "RR\n",
    title = "Caudal",
    color = ""
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca de Eo-Navia

# Sin corregir por el sesgo

cpalle1 <- crosspred(
  bvarj[[3]],
  coef = yall[3, ],
  vcov = Sall$eonavia,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpalle1$predvar[which.min(cpalle1$allRRfit)]

cpalle1 <- crosspred(
  bvarj[[3]],
  coef = yall[3, ],
  vcov = Sall$eonavia,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# BLUP

cpalle2 <- crosspred(
  bvarj[[3]],
  coef = bl$eonavia$blup,
  vcov = bl$eonavia$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpalle2$predvar[which.min(cpalle2$allRRfit)]

pci <- quantile(
  findmin(
    basis = bvarj[[3]],
    coef = bl$eonavia$blup,
    vcov = bl$eonavia$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

cpalle2 <- crosspred(
  bvarj[[3]],
  coef = bl$eonavia$blup,
  vcov = bl$eonavia$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

cpalle2$allRRlow
cpalle2$allRRfit
cpalle2$allRRhigh

# Representación de la curva sin corregir y de la BLUP

RR_overalle <- data.frame(
  RR = cpalle1$allRRfit,
  LowRR = cpalle1$allRRlow,
  HighRR = cpalle1$allRRhigh
)

RR_overalle <- rownames_to_column(RR_overalle, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overalle$tmean[RR_overalle$RR == 1]

xlab <- pretty(RR_overalle$tmean, n = 10)

RR_overalle2 <- data.frame(
  RR2 = cpalle2$allRRfit,
  LowRR2 = cpalle2$allRRlow,
  HighRR2 = cpalle2$allRRhigh
)

RR_overalle2 <- rownames_to_column(RR_overalle2, "tmean2") %>% mutate_at(1, as.numeric)

mmt2 <- RR_overalle2$tmean2[RR_overalle2$RR2 == 1]

xlab2 <- pretty(RR_overalle2$tmean2, n = 10)
ylab2 <- pretty(c(RR_overalle2$LowRR, RR_overalle2$HighRR), n = 10)

colors <- c("Comarca" = "olivedrab4",
            "BLUP" = "royalblue4")

RR_overalle %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "olivedrab3",
              alpha = 0.5) +
  geom_line(aes(tmean, RR, colour = "Comarca"), size = 1) +
  geom_ribbon(
    aes(ymin = RR_overalle2$LowRR2, ymax = RR_overalle2$HighRR2),
    fill = "royalblue3",
    alpha = 0.5
  ) +
  geom_line(aes(RR_overalle2$tmean2, RR_overalle2$RR2, colour = "BLUP"),
            size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "olivedrab4",
    size = 2,
    colour = "olivedrab4",
    show.legend = FALSE
  ) +
  geom_point(
    aes(mmt2, 1),
    shape = 21,
    fill = "royalblue4",
    size = 2,
    colour = "royalblue4",
    show.legend = T
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nTemperatura media (ºC)",
    y = "RR\n",
    title = "Eo-Navia",
    color = ""
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca de Gijón

# Sin corregir por el sesgo

cpallg1 <- crosspred(
  bvarj[[4]],
  coef = yall[4, ],
  vcov = Sall$gijon,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallg1$predvar[which.min(cpallg1$allRRfit)]

cpallg1 <- crosspred(
  bvarj[[4]],
  coef = yall[4, ],
  vcov = Sall$gijon,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# BLUP

cpallg2 <- crosspred(
  bvarj[[4]],
  coef = bl$gijon$blup,
  vcov = bl$gijon$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallg2$predvar[which.min(cpallg2$allRRfit)]

pci <- quantile(
  findmin(
    basis = bvarj[[4]],
    coef = bl$gijon$blup,
    vcov = bl$gijon$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

cpallg2 <- crosspred(
  bvarj[[4]],
  coef = bl$gijon$blup,
  vcov = bl$gijon$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

cpallg2$allRRlow
cpallg2$allRRfit
cpallg2$allRRhigh

# Representación de las dos curvas

RR_overallg <- data.frame(
  RR = cpallg1$allRRfit,
  LowRR = cpallg1$allRRlow,
  HighRR = cpallg1$allRRhigh
)

RR_overallg <- rownames_to_column(RR_overallg, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overallg$tmean[RR_overallg$RR == 1]

xlab <- pretty(RR_overallg$tmean, n = 10)

RR_overallg2 <- data.frame(
  RR2 = cpallg2$allRRfit,
  LowRR2 = cpallg2$allRRlow,
  HighRR2 = cpallg2$allRRhigh
)

RR_overallg2 <- rownames_to_column(RR_overallg2, "tmean2") %>% mutate_at(1, as.numeric)

mmt2 <- RR_overallg2$tmean2[RR_overallg2$RR2 == 1]

xlab2 <- pretty(RR_overallg2$tmean2, n = 10)
ylab2 <- pretty(c(RR_overallg2$LowRR, RR_overallg2$HighRR), n = 10)

colors <- c("Comarca" = "olivedrab4",
            "BLUP" = "royalblue4")

RR_overallg %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "olivedrab3",
              alpha = 0.5) +
  geom_line(aes(tmean, RR, colour = "Comarca"), size = 1) +
  geom_ribbon(
    aes(ymin = RR_overallg2$LowRR2, ymax = RR_overallg2$HighRR2),
    fill = "royalblue3",
    alpha = 0.5
  ) +
  geom_line(aes(RR_overallg2$tmean2, RR_overallg2$RR2, colour = "BLUP"),
            size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "olivedrab4",
    size = 2,
    colour = "olivedrab4",
    show.legend = FALSE
  ) +
  geom_point(
    aes(mmt2, 1),
    shape = 21,
    fill = "royalblue4",
    size = 2,
    colour = "royalblue4",
    show.legend = T
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nTemperatura media (ºC)",
    y = "RR\n",
    title = "Gijón",
    color = ""
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca del Nalón

# Sin corregir

cpalln1 <- crosspred(
  bvarj[[5]],
  coef = yall[5, ],
  vcov = Sall$nalon,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpalln1$predvar[which.min(cpalln1$allRRfit)]

cpalln1 <- crosspred(
  bvarj[[5]],
  coef = yall[5, ],
  vcov = Sall$nalon,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# BLUP

cpalln2 <- crosspred(
  bvarj[[5]],
  coef = bl$nalon$blup,
  vcov = bl$nalon$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpalln2$predvar[which.min(cpalln2$allRRfit)]

pci <- quantile(
  findmin(
    basis = bvarj[[5]],
    coef = bl$nalon$blup,
    vcov = bl$nalon$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

cpalln2 <- crosspred(
  bvarj[[5]],
  coef = bl$nalon$blup,
  vcov = bl$nalon$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

cpalln2$allRRlow
cpalln2$allRRfit
cpalln2$allRRhigh

# Representación de las dos curvas

RR_overalln <- data.frame(
  RR = cpalln1$allRRfit,
  LowRR = cpalln1$allRRlow,
  HighRR = cpalln1$allRRhigh
)

RR_overalln <- rownames_to_column(RR_overalln, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overalln$tmean[RR_overalln$RR == 1]

xlab <- pretty(RR_overalln$tmean, n = 10)

RR_overalln2 <- data.frame(
  RR2 = cpalln2$allRRfit,
  LowRR2 = cpalln2$allRRlow,
  HighRR2 = cpalln2$allRRhigh
)

RR_overalln2 <- rownames_to_column(RR_overalln2, "tmean2") %>% mutate_at(1, as.numeric)

mmt2 <- RR_overalln2$tmean2[RR_overalln2$RR2 == 1]

xlab2 <- pretty(RR_overalln2$tmean2, n = 10)
ylab2 <- pretty(c(RR_overalln2$LowRR, RR_overalln2$HighRR), n = 10)

colors <- c("Comarca" = "olivedrab4",
            "BLUP" = "royalblue4")

RR_overalln %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "olivedrab3",
              alpha = 0.5) +
  geom_line(aes(tmean, RR, colour = "Comarca"), size = 1) +
  geom_ribbon(
    aes(ymin = RR_overalln2$LowRR2, ymax = RR_overalln2$HighRR2),
    fill = "royalblue3",
    alpha = 0.5
  ) +
  geom_line(aes(RR_overalln2$tmean2, RR_overalln2$RR2, colour = "BLUP"),
            size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "olivedrab4",
    size = 2,
    colour = "olivedrab4",
    show.legend = FALSE
  ) +
  geom_point(
    aes(mmt2, 1),
    shape = 21,
    fill = "royalblue4",
    size = 2,
    colour = "royalblue4",
    show.legend = T
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nTemperatura media (ºC)",
    y = "RR\n",
    title = "Nalón",
    color = ""
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca de Narcea

# Sin corregir

cpallnr1 <- crosspred(
  bvarj[[6]],
  coef = yall[6, ],
  vcov = Sall$narcea,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallnr1$predvar[which.min(cpallnr1$allRRfit)]

cpallnr1 <- crosspred(
  bvarj[[6]],
  coef = yall[6, ],
  vcov = Sall$narcea,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# BLUP

cpallnr2 <- crosspred(
  bvarj[[6]],
  coef = bl$narcea$blup,
  vcov = bl$narcea$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallnr2$predvar[which.min(cpallnr2$allRRfit)]

pci <- quantile(
  findmin(
    basis = bvarj[[6]],
    coef = bl$narcea$blup,
    vcov = bl$narcea$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

cpallnr2 <- crosspred(
  bvarj[[6]],
  coef = bl$narcea$blup,
  vcov = bl$narcea$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

cpallnr2$allRRlow
cpallnr2$allRRfit
cpallnr2$allRRhigh

# Representación de las dos curvas

RR_overallnr <- data.frame(
  RR = cpallnr1$allRRfit,
  LowRR = cpallnr1$allRRlow,
  HighRR = cpallnr1$allRRhigh
)

RR_overallnr <- rownames_to_column(RR_overallnr, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overallnr$tmean[RR_overallnr$RR == 1]

xlab <- pretty(RR_overallnr$tmean, n = 10)

RR_overallnr2 <- data.frame(
  RR2 = cpallnr2$allRRfit,
  LowRR2 = cpallnr2$allRRlow,
  HighRR2 = cpallnr2$allRRhigh
)

RR_overallnr2 <- rownames_to_column(RR_overallnr2, "tmean2") %>% mutate_at(1, as.numeric)

mmt2 <- RR_overallnr2$tmean2[RR_overallnr2$RR2 == 1]

xlab2 <- pretty(RR_overallnr2$tmean2, n = 10)
ylab2 <- pretty(c(RR_overallnr2$LowRR, RR_overallnr2$HighRR), n = 10)

colors <- c("Comarca" = "olivedrab4",
            "BLUP" = "royalblue4")

RR_overallnr %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "olivedrab3",
              alpha = 0.5) +
  geom_line(aes(tmean, RR, colour = "Comarca"), size = 1) +
  geom_ribbon(
    aes(ymin = RR_overallnr2$LowRR2, ymax = RR_overallnr2$HighRR2),
    fill = "royalblue3",
    alpha = 0.5
  ) +
  geom_line(aes(RR_overallnr2$tmean2, RR_overallnr2$RR2, colour = "BLUP"),
            size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "olivedrab4",
    size = 2,
    colour = "olivedrab4",
    show.legend = FALSE
  ) +
  geom_point(
    aes(mmt2, 1),
    shape = 21,
    fill = "royalblue4",
    size = 2,
    colour = "royalblue4",
    show.legend = T
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nTemperatura media (ºC)",
    y = "RR\n",
    title = "Narcea",
    color = ""
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca del Oriente

# Sin corregir

cpallo1 <- crosspred(
  bvarj[[7]],
  coef = yall[7, ],
  vcov = Sall$oriente,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallo1$predvar[which.min(cpallo1$allRRfit)]

cpallo1 <- crosspred(
  bvarj[[7]],
  coef = yall[7, ],
  vcov = Sall$oriente,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# BLUP

cpallo2 <- crosspred(
  bvarj[[7]],
  coef = bl$oriente$blup,
  vcov = bl$oriente$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallo2$predvar[which.min(cpallo2$allRRfit)]

pci <- quantile(
  findmin(
    basis = bvarj[[7]],
    coef = bl$oriente$blup,
    vcov = bl$oriente$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

cpallo2 <- crosspred(
  bvarj[[7]],
  coef = bl$oriente$blup,
  vcov = bl$oriente$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

cpallo2$allRRlow
cpallo2$allRRfit
cpallo2$allRRhigh

# Representación de las dos curvas

RR_overallo <- data.frame(
  RR = cpallo1$allRRfit,
  LowRR = cpallo1$allRRlow,
  HighRR = cpallo1$allRRhigh
)

RR_overallo <- rownames_to_column(RR_overallo, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overallo$tmean[RR_overallo$RR == 1]

xlab <- pretty(RR_overallo$tmean, n = 10)

RR_overallo2 <- data.frame(
  RR2 = cpallo2$allRRfit,
  LowRR2 = cpallo2$allRRlow,
  HighRR2 = cpallo2$allRRhigh
)

RR_overallo2 <- rownames_to_column(RR_overallo2, "tmean2") %>% mutate_at(1, as.numeric)

mmt2 <- RR_overallo2$tmean2[RR_overallo2$RR2 == 1]

xlab2 <- pretty(RR_overallo2$tmean2, n = 10)
ylab2 <- pretty(c(RR_overallo2$LowRR, RR_overallo2$HighRR), n = 10)

colors <- c("Comarca" = "olivedrab4",
            "BLUP" = "royalblue4")

RR_overallo %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "olivedrab3",
              alpha = 0.5) +
  geom_line(aes(tmean, RR, colour = "Comarca"), size = 1) +
  geom_ribbon(
    aes(ymin = RR_overallo2$LowRR2, ymax = RR_overallo2$HighRR2),
    fill = "royalblue3",
    alpha = 0.5
  ) +
  geom_line(aes(RR_overallo2$tmean2, RR_overallo2$RR2, colour = "BLUP"),
            size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "olivedrab4",
    size = 2,
    colour = "olivedrab4",
    show.legend = FALSE
  ) +
  geom_point(
    aes(mmt2, 1),
    shape = 21,
    fill = "royalblue4",
    size = 2,
    colour = "royalblue4",
    show.legend = T
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nTemperatura media (ºC)",
    y = "RR\n",
    title = "Oriente",
    color = ""
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca del Oviedo

# Sin corregir

cpallov1 <- crosspred(
  bvarj[[8]],
  coef = yall[8, ],
  vcov = Sall$oviedo,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallov1$predvar[which.min(cpallov1$allRRfit)]

cpallov1 <- crosspred(
  bvarj[[8]],
  coef = yall[8, ],
  vcov = Sall$oviedo,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

# BLUP

cpallov2 <- crosspred(
  bvarj[[8]],
  coef = bl$oviedo$blup,
  vcov = bl$oviedo$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26
)

mmt <- cpallov2$predvar[which.min(cpallov2$allRRfit)]

pci <- quantile(
  findmin(
    basis = bvarj[[8]],
    coef = bl$oviedo$blup,
    vcov = bl$oviedo$vcov,
    from = 0,
    to = 26,
    by = 1,
    sim = T,
    nsim = 10000
  ),
  c(2.5, 97.5) / 100
)

cpallov2 <- crosspred(
  bvarj[[8]],
  coef = bl$oviedo$blup,
  vcov = bl$oviedo$vcov,
  model.link = "log",
  by = 1,
  from = 0,
  to = 26,
  cen = mmt
)

cpallov2$allRRlow
cpallov2$allRRfit
cpallov2$allRRhigh

# Representación de las dos curvas

RR_overallov <- data.frame(
  RR = cpallov1$allRRfit,
  LowRR = cpallov1$allRRlow,
  HighRR = cpallov1$allRRhigh
)

RR_overallov <- rownames_to_column(RR_overallov, "tmean") %>% mutate_at(1, as.numeric)

mmt <- RR_overallov$tmean[RR_overallov$RR == 1]

xlab <- pretty(RR_overallov$tmean, n = 10)

RR_overallov2 <- data.frame(
  RR2 = cpallov2$allRRfit,
  LowRR2 = cpallov2$allRRlow,
  HighRR2 = cpallov2$allRRhigh
)

RR_overallov2 <- rownames_to_column(RR_overallov2, "tmean2") %>% mutate_at(1, as.numeric)

mmt2 <- RR_overallov2$tmean2[RR_overallov2$RR2 == 1]

xlab2 <- pretty(RR_overallov2$tmean2, n = 10)
ylab2 <- pretty(c(RR_overallov2$LowRR, RR_overallov2$HighRR), n = 10)

colors <- c("Comarca" = "olivedrab4",
            "BLUP" = "royalblue4")

RR_overallov %>%
  ggplot(aes(tmean, RR)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  geom_ribbon(aes(ymin = LowRR, ymax = HighRR),
              fill = "olivedrab3",
              alpha = 0.5) +
  geom_line(aes(tmean, RR, colour = "Comarca"), size = 1) +
  geom_ribbon(
    aes(ymin = RR_overallov2$LowRR2, ymax = RR_overallov2$HighRR2),
    fill = "royalblue3",
    alpha = 0.5
  ) +
  geom_line(aes(RR_overallov2$tmean2, RR_overallov2$RR2, colour = "BLUP"),
            size = 1) +
  geom_point(
    aes(mmt, 1),
    shape = 21,
    fill = "olivedrab4",
    size = 2,
    colour = "olivedrab4",
    show.legend = FALSE
  ) +
  geom_point(
    aes(mmt2, 1),
    shape = 21,
    fill = "royalblue4",
    size = 2,
    colour = "royalblue4",
    show.legend = T
  ) +
  scale_x_continuous(breaks = xlab) +
  scale_y_continuous(breaks = ylab) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nTemperatura media (ºC)",
    y = "RR\n",
    title = "Oviedo",
    color = ""
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Efecto de las vacaciones de verano

################################################################
# PRIMERA ETAPA

# Realizo la primera etapa del estudio comarca-específico.
# Consiste en ajustar un modelo para cada comarca y guardar
# los resultados relativos a las vacaciones de verano
################################################################

yall <- matrix(NA, length(ld), 3, dimnames = list(comarca, paste("b", seq(3), sep = "")))

Sall <- vector("list", length(ld))
names(Sall) <- comarca

blagj <- vector("list", 8)
names(blagj) <- comarca

system.time({
  for(i in seq(ld)) {
    sub <- ld[[i]]
    
    suppressWarnings({
      sub <- sub %>%
        mutate(vacaciones = ifelse(mm == 12 & dd == 25, "Otras vacaciones", ifelse(mm == 1 & dd == 1, "Otras vacaciones", ifelse(mm == 1 & dd == 6|mm == 5 & dd == 1|mm == 9 & dd == 8|mm == 10 & dd == 12|mm == 11 & dd == 1|mm == 12 & dd == 6|mm == 12 & dd == 8, "Otras vacaciones", ifelse(mm == 7 & dd >= 15|mm == 8 & dd <= 15, "Verano", "No vacaciones")))))
      
      sub[c(106:107, 456:457, 813:814, 1198:1199, 1548:1549, 1933:1934, 2290:2291, 2647:2648, 3025:3026, 3382:3383, 3739:3740, 4117:4118, 4474:4475, 4859:4860, 5216:5217, 5566:5567, 5951:5952, 6308:6309, 6658:6659, 7043:7044, 7400:7401, 7750:7751, 8135:8136, 8492:8493, 8877:8878, 9227:9228, 9584:9585, 9969:9970, 10319:10320, 10676:10677, 11061:11062, 11411:11412, 11796:11797, 12153:12154, 12510:12511), "vacaciones"] <- "Otras vacaciones"
      
      sub$vacaciones <- as.factor(sub$vacaciones)
      sub$vacaciones <- relevel(sub$vacaciones, ref = "No vacaciones")
      
      sub$stratum <- as.factor(as.factor(sub$yy))
      ind <- tapply(sub$Total, sub$stratum, sum)[sub$stratum]
      
      mi.argvarj <- list(fun = "poly", degree = 2, int = F)
      
      nlagj <- 4
      lagnk <- 2; klagj <- logknots(nlagj, lagnk)
      mi.arglagj <- list(fun = "ns", knots = klagj, int = T)
      cbTj <- crossbasis(sub$tmed, lag = nlagj, argvar = mi.argvarj, arglag = mi.arglagj)
      
      vacan <- as.factor(as.numeric(sub$vacaciones))
      mi.arglagjv <- list(fun = "strata", breaks = c(0, 1))
      cbvacaj <- crossbasis(vacan, lag = c(-2, 2), argvar = list(fun = "integer", int = F), arglag = mi.arglagjv)
      
      dfseas <- 1
      ny <- length(unique(sub$yy))
    })
    
    mfirst <- gnm(Total ~ cbTj + cbvacaj + dow, data = sub, family = quasipoisson, eliminate = factor(stratum), na.action = "na.exclude", subset = ind > 0)
    tpred <- seq(0, 26, by = 1)
    
    xlag <- -2:2
    blagj[[i]] <- do.call("onebasis", c(list(x = xlag), attr(cbvacaj, "arglag")))

    suppressWarnings({
      crall <- crossreduce(cbvacaj, mfirst, type = "var", value = 3, cen = 1)
    })

    yall[i,] <- coef(crall) 
    Sall[[i]] <- vcov(crall) 
    
  }
})

################################################################
# SEGUNDA ETAPA

# Realizo la segunda etapa del estudio comarca-específico.
# Consiste en combinar los resultados de cada comarca mediante
# meta-análisis y obtener un único efecto general.
################################################################

# Ajuste del modelo meta-analítico

method <- "reml"

mvall <- mixmeta(yall ~ 1, Sall, method = method)
summary(mvall)

xlag <- -2:2
blag <- do.call("onebasis", c(list(x = xlag), attr(cbvacaj, "arglag")))

cpall <- crosspred(
  blag,
  coef = coef(mvall),
  vcov = vcov(mvall),
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpall$allRRlow
cpall$allRRfit
cpall$allRRhigh

# Representación efecto general combinado

p <- data.frame(
  RR = cpall$allRRfit,
  LowRR = cpall$allRRlow,
  HighRR = cpall$allRRhigh
)

p <- rownames_to_column(p, "lag") %>% mutate_at(1, as.numeric)

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
  labs(x = "\nLag", y = "RR\n", title = "Vacaciones de verano") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

################################################################
# BLUP

# Estimo las BLUPs de cada comarca en base al efecto estimado
# en cada comarca y al efecto general obtenido combinando los
# resultados
################################################################

bl <- blup(mvall, vcov = T)

## Comarca de Avilés

# Sin corregir

cpalla1 <- crosspred(
  blagj[[1]],
  coef = yall[1, ],
  vcov = Sall$aviles,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpalla2 <- crosspred(
  blagj[[1]],
  coef = bl$aviles$blup,
  vcov = bl$aviles$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpalla2$allRRlow
cpalla2$allRRfit
cpalla2$allRRhigh

# Representación efectos

pa1 <- data.frame(
  RR = cpalla1$allRRfit,
  LowRR = cpalla1$allRRlow,
  HighRR = cpalla1$allRRhigh
)

pa1 <- rownames_to_column(pa1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pa1)))

pa2 <- data.frame(
  RR = cpalla2$allRRfit,
  LowRR = cpalla2$allRRlow,
  HighRR = cpalla2$allRRhigh
)

pa2 <- rownames_to_column(pa2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pa2)))

pa <- rbind(pa1, pa2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pa, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Avilés"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapae = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca del Caudal

# Sin corregir

cpallc1 <- crosspred(
  blagj[[2]],
  coef = yall[2, ],
  vcov = Sall$caudal,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallc2 <- crosspred(
  blagj[[2]],
  coef = bl$caudal$blup,
  vcov = bl$caudal$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallc2$allRRlow
cpallc2$allRRfit
cpallc2$allRRhigh

# Representación efectos

pc1 <- data.frame(
  RR = cpallc1$allRRfit,
  LowRR = cpallc1$allRRlow,
  HighRR = cpallc1$allRRhigh
)

pc1 <- rownames_to_column(pc1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo =  rep("Comarca", nrow(pc1)))

pc2 <- data.frame(
  RR = cpallc2$allRRfit,
  LowRR = cpallc2$allRRlow,
  HighRR = cpallc2$allRRhigh
)

pc2 <- rownames_to_column(pc2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pc2)))

pc <- rbind(pc1, pc2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pc, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Caudal"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapce = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca de Eo-Navia

# Sin corregir

cpalle1 <- crosspred(
  blagj[[3]],
  coef = yall[3, ],
  vcov = Sall$eonavia,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpalle2 <- crosspred(
  blagj[[3]],
  coef = bl$eonavia$blup,
  vcov = bl$eonavia$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpalle2$allRRlow
cpalle2$allRRfit
cpalle2$allRRhigh

# Representación efectos

pe1 <- data.frame(
  RR = cpalle1$allRRfit,
  LowRR = cpalle1$allRRlow,
  HighRR = cpalle1$allRRhigh
)

pe1 <- rownames_to_column(pe1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pe1)))

pe2 <- data.frame(
  RR = cpalle2$allRRfit,
  LowRR = cpalle2$allRRlow,
  HighRR = cpalle2$allRRhigh
)

pe2 <- rownames_to_column(pe2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pe2)))

pe <- rbind(pe1, pe2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pe, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Eo-Navia"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapee = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca de Gijón

# Sin corregir

cpallg1 <- crosspred(
  blagj[[4]],
  coef = yall[4, ],
  vcov = Sall$gijon,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallg2 <- crosspred(
  blagj[[4]],
  coef = bl$gijon$blup,
  vcov = bl$gijon$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallg2$allRRlow
cpallg2$allRRfit
cpallg2$allRRhigh

# Representación efectos

pg1 <- data.frame(
  RR = cpallg1$allRRfit,
  LowRR = cpallg1$allRRlow,
  HighRR = cpallg1$allRRhigh
)

pg1 <- rownames_to_column(pg1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pg1)))

pg2 <- data.frame(
  RR = cpallg2$allRRfit,
  LowRR = cpallg2$allRRlow,
  HighRR = cpallg2$allRRhigh
)

pg2 <- rownames_to_column(pg2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pg2)))

pg <- rbind(pg1, pg2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pg, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Gijón"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapge = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca del Nalón

# Sin corregir

cpalln1 <- crosspred(
  blagj[[5]],
  coef = yall[5, ],
  vcov = Sall$nalon,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpalln2 <- crosspred(
  blagj[[5]],
  coef = bl$nalon$blup,
  vcov = bl$nalon$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpalln2$allRRlow
cpalln2$allRRfit
cpalln2$allRRhigh

# Representación efectos

pn1 <- data.frame(
  RR = cpalln1$allRRfit,
  LowRR = cpalln1$allRRlow,
  HighRR = cpalln1$allRRhigh
)

pn1 <- rownames_to_column(pn1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo =  rep("Comarca", nrow(pn1)))

pn2 <- data.frame(
  RR = cpalln2$allRRfit,
  LowRR = cpalln2$allRRlow,
  HighRR = cpalln2$allRRhigh
)

pn2 <- rownames_to_column(pn2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo =  rep("BLUP", nrow(pn2)))

pn <- rbind(pn1, pn2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pn, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Nalón"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapne = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca de Narcea

# Sin corregir

cpallnr1 <- crosspred(
  blagj[[6]],
  coef = yall[6, ],
  vcov = Sall$narcea,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallnr2 <- crosspred(
  blagj[[6]],
  coef = bl$narcea$blup,
  vcov = bl$narcea$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallnr2$allRRlow
cpallnr2$allRRfit
cpallnr2$allRRhigh

# Representación efectos

pnr1 <- data.frame(
  RR = cpallnr1$allRRfit,
  LowRR = cpallnr1$allRRlow,
  HighRR = cpallnr1$allRRhigh
)

pnr1 <- rownames_to_column(pnr1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pnr1)))

pnr2 <- data.frame(
  RR = cpallnr2$allRRfit,
  LowRR = cpallnr2$allRRlow,
  HighRR = cpallnr2$allRRhigh
)

pnr2 <- rownames_to_column(pnr2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pnr2)))

pnr <- rbind(pnr1, pnr2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pnr, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Narcea"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapnre = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca del Oriente

# Sin corregir

cpallo1 <- crosspred(
  blagj[[7]],
  coef = yall[7, ],
  vcov = Sall$oriente,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallo2 <- crosspred(
  blagj[[7]],
  coef = bl$oriente$blup,
  vcov = bl$oriente$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallo2$allRRlow
cpallo2$allRRfit
cpallo2$allRRhigh

# Representación efectos

po1 <- data.frame(
  RR = cpallo1$allRRfit,
  LowRR = cpallo1$allRRlow,
  HighRR = cpallo1$allRRhigh
)

po1 <- rownames_to_column(po1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(po1)))

po2 <- data.frame(
  RR = cpallo2$allRRfit,
  LowRR = cpallo2$allRRlow,
  HighRR = cpallo2$allRRhigh
)

po2 <- rownames_to_column(po2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(po2)))

po <- rbind(po1, po2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(po, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Oriente"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapoe = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca del Oviedo

# Sin corregir

cpallov1 <- crosspred(
  blagj[[8]],
  coef = yall[8, ],
  vcov = Sall$oviedo,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallov2 <- crosspred(
  blagj[[8]],
  coef = bl$oviedo$blup,
  vcov = bl$oviedo$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallov2$allRRlow
cpallov2$allRRfit
cpallov2$allRRhigh

# Representación efectos

pov1 <- data.frame(
  RR = cpallov1$allRRfit,
  LowRR = cpallov1$allRRlow,
  HighRR = cpallov1$allRRhigh
)

pov1 <- rownames_to_column(pov1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pov1)))

pov2 <- data.frame(
  RR = cpallov2$allRRfit,
  LowRR = cpallov2$allRRlow,
  HighRR = cpallov2$allRRhigh
)

pov2 <- rownames_to_column(pov2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pov2)))

pov <- rbind(pov1, pov2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pov, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Oviedo"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapove = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Efecto de las otras vacaciones

################################################################
# PRIMERA ETAPA

# Realizo la primera etapa del estudio comarca-específico.
# Consiste en ajustar un modelo para cada comarca y guardar
# los resultados relativos a las otras vacaciones
################################################################

yall <- matrix(NA, length(ld), 3, dimnames = list(comarca, paste("b", seq(3), sep = "")))

Sall <- vector("list", length(ld))
names(Sall) <- comarca

blagj <- vector("list", 8)
names(blagj) <- comarca

system.time({
  for(i in seq(ld)) {
    sub <- ld[[i]]
    
    suppressWarnings({
      sub <- sub %>%
        mutate(vacaciones = ifelse(mm == 12 & dd == 25, "Otras vacaciones", ifelse(mm == 1 & dd == 1, "Otras vacaciones", ifelse(mm == 1 & dd == 6|mm == 5 & dd == 1|mm == 9 & dd == 8|mm == 10 & dd == 12|mm == 11 & dd == 1|mm == 12 & dd == 6|mm == 12 & dd == 8, "Otras vacaciones", ifelse(mm == 7 & dd >= 15|mm == 8 & dd <= 15, "Verano", "No vacaciones")))))
      
      sub[c(106:107, 456:457, 813:814, 1198:1199, 1548:1549, 1933:1934, 2290:2291, 2647:2648, 3025:3026, 3382:3383, 3739:3740, 4117:4118, 4474:4475, 4859:4860, 5216:5217, 5566:5567, 5951:5952, 6308:6309, 6658:6659, 7043:7044, 7400:7401, 7750:7751, 8135:8136, 8492:8493, 8877:8878, 9227:9228, 9584:9585, 9969:9970, 10319:10320, 10676:10677, 11061:11062, 11411:11412, 11796:11797, 12153:12154, 12510:12511), "vacaciones"] <- "Otras vacaciones"
      
      sub$vacaciones <- as.factor(sub$vacaciones)
      sub$vacaciones <- relevel(sub$vacaciones, ref = "No vacaciones")
      
      sub$stratum <- as.factor(as.factor(sub$yy))
      ind <- tapply(sub$Total, sub$stratum, sum)[sub$stratum]
      
      mi.argvarj <- list(fun = "poly", degree = 2, int = F)
      
      nlagj <- 4
      lagnk <- 2; klagj <- logknots(nlagj,lagnk)
      mi.arglagj <- list(fun = "ns", knots = klagj, int = T)
      cbTj <- crossbasis(sub$tmed, lag = nlagj, argvar = mi.argvarj, arglag = mi.arglagj)
      
      vacan <- as.factor(as.numeric(sub$vacaciones))
      mi.arglagjv <- list(fun = "strata", breaks=c(0, 1))
      cbvacaj <- crossbasis(vacan, lag = c(-2, 2), argvar = list(fun = "integer", int = F), arglag = mi.arglagjv)
      
      dfseas <- 1
      ny <- length(unique(sub$yy))
    })
    
    mfirst <- gnm(Total ~ cbTj + cbvacaj + dow, data = sub, family = quasipoisson, eliminate = factor(stratum), na.action = "na.exclude", subset = ind > 0)
    tpred <- seq(0, 26, by = 1)
    
    xlag <- -2:2
    blagj[[i]] <- do.call("onebasis", c(list(x = xlag), attr(cbvacaj, "arglag")))
    
    suppressWarnings({
      crall <- crossreduce(cbvacaj, mfirst, type = "var", value = 2, cen = 1)
    })

    yall[i,] <- coef(crall) 
    Sall[[i]] <- vcov(crall) 
    
  }
})

################################################################
# SEGUNDA ETAPA

# Realizo la segunda etapa del estudio comarca-específico.
# Consiste en combinar los resultados de cada comarca mediante
# meta-análisis y obtener un único efecto general.
################################################################

# Ajuste del modelo meta-analítico

method <- "reml"
mvall <- mixmeta(yall ~ 1, Sall, method = method)
summary(mvall)

xlag <- -2:2
blag <- do.call("onebasis", c(list(x = xlag), attr(cbvacaj, "arglag")))

cpall <- crosspred(
  blag,
  coef = coef(mvall),
  vcov = vcov(mvall),
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpall$allRRlow
cpall$allRRfit
cpall$allRRhigh

# Representación del efecto general combinado

p <- data.frame(
  RR = cpall$allRRfit,
  LowRR = cpall$allRRlow,
  HighRR = cpall$allRRhigh
)

p <- rownames_to_column(p, "lag") %>% mutate_at(1, as.numeric)

ggplot(p, aes(x = lag, y = RR)) +
  geom_errorbar(aes(ymin = LowRR, ymax = HighRR, width = 0.1),
                colour = "purple",
                size = 0.8) +
  geom_line(colour = "purple", size = 0.8) +
  geom_point(colour = "purple4", size = 3) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "\nRetardo", y = "RR\n", title = "Otras vacaciones") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))


################################################################
# BLUP

# Estimo las BLUPs de cada comarca en base al efecto estimado
# en cada comarca y al efecto general obtenido combinando los
# resultados
################################################################

bl <- blup(mvall, vcov = T)

## Comarca de Avilés

# Sin corregir

cpalla1 <- crosspred(
  blagj[[1]],
  coef = yall[1, ],
  vcov = Sall$aviles,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpalla2 <- crosspred(
  blagj[[1]],
  coef = bl$aviles$blup,
  vcov = bl$aviles$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpalla2$allRRlow
cpalla2$allRRfit
cpalla2$allRRhigh

# Representación de los dos efectos

pa1 <- data.frame(
  RR = cpalla1$allRRfit,
  LowRR = cpalla1$allRRlow,
  HighRR = cpalla1$allRRhigh
)

pa1 <- rownames_to_column(pa1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo =  rep("Comarca", nrow(pa1)))

pa2 <- data.frame(
  RR = cpalla2$allRRfit,
  LowRR = cpalla2$allRRlow,
  HighRR = cpalla2$allRRhigh
)
pa2 <- rownames_to_column(pa2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pa2)))

pa <- rbind(pa1, pa2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pa, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.1),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Avilés"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapae = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca del Caudal

# Sin corregir

cpallc1 <- crosspred(
  blagj[[2]],
  coef = yall[2, ],
  vcov = Sall$caudal,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallc2 <- crosspred(
  blagj[[2]],
  coef = bl$caudal$blup,
  vcov = bl$caudal$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallc2$allRRlow
cpallc2$allRRfit
cpallc2$allRRhigh

# Representación de los dos efectos

pc1 <- data.frame(
  RR = cpallc1$allRRfit,
  LowRR = cpallc1$allRRlow,
  HighRR = cpallc1$allRRhigh
)

pc1 <- rownames_to_column(pc1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pc1)))

pc2 <- data.frame(
  RR = cpallc2$allRRfit,
  LowRR = cpallc2$allRRlow,
  HighRR = cpallc2$allRRhigh
)

pc2 <- rownames_to_column(pc2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pc2)))

pc <- rbind(pc1, pc2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pc, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Caudal"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapce = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca de Eo-Navia

# Sin corregir

cpalle1 <- crosspred(
  blagj[[3]],
  coef = yall[3, ],
  vcov = Sall$eonavia,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpalle2 <- crosspred(
  blagj[[3]],
  coef = bl$eonavia$blup,
  vcov = bl$eonavia$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpalle2$allRRlow
cpalle2$allRRfit
cpalle2$allRRhigh

# Representación de los dos efectos

pe1 <- data.frame(
  RR = cpalle1$allRRfit,
  LowRR = cpalle1$allRRlow,
  HighRR = cpalle1$allRRhigh
)

pe1 <- rownames_to_column(pe1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pe1)))

pe2 <- data.frame(
  RR = cpalle2$allRRfit,
  LowRR = cpalle2$allRRlow,
  HighRR = cpalle2$allRRhigh
)

pe2 <- rownames_to_column(pe2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pe2)))

pe <- rbind(pe1, pe2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pe, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Eo-Navia"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapee = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca de Gijón

# Sin corregir

cpallg1 <- crosspred(
  blagj[[4]],
  coef = yall[4, ],
  vcov = Sall$gijon,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallg2 <- crosspred(
  blagj[[4]],
  coef = bl$gijon$blup,
  vcov = bl$gijon$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallg2$allRRlow
cpallg2$allRRfit
cpallg2$allRRhigh

# Representación de los dos efectos

pg1 <- data.frame(
  RR = cpallg1$allRRfit,
  LowRR = cpallg1$allRRlow,
  HighRR = cpallg1$allRRhigh
)

pg1 <- rownames_to_column(pg1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pg1)))

pg2 <- data.frame(
  RR = cpallg2$allRRfit,
  LowRR = cpallg2$allRRlow,
  HighRR = cpallg2$allRRhigh
)

pg2 <- rownames_to_column(pg2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pg2)))

pg <- rbind(pg1, pg2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pg, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Gijón"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapge = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca del Nalón

# Sin corregir

cpalln1 <- crosspred(
  blagj[[5]],
  coef = yall[5, ],
  vcov = Sall$nalon,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpalln2 <- crosspred(
  blagj[[5]],
  coef = bl$nalon$blup,
  vcov = bl$nalon$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpalln2$allRRlow
cpalln2$allRRfit
cpalln2$allRRhigh

# Representación de los dos efectos

pn1 <- data.frame(
  RR = cpalln1$allRRfit,
  LowRR = cpalln1$allRRlow,
  HighRR = cpalln1$allRRhigh
)

pn1 <- rownames_to_column(pn1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pn1)))

pn2 <- data.frame(
  RR = cpalln2$allRRfit,
  LowRR = cpalln2$allRRlow,
  HighRR = cpalln2$allRRhigh
)

pn2 <- rownames_to_column(pn2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pn2)))

pn <- rbind(pn1, pn2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pn, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Nalón"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapne = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca del Narcea

# Sin corregir

cpallnr1 <- crosspred(
  blagj[[6]],
  coef = yall[6, ],
  vcov = Sall$narcea,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallnr2 <- crosspred(
  blagj[[6]],
  coef = bl$narcea$blup,
  vcov = bl$narcea$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallnr2$allRRlow
cpallnr2$allRRfit
cpallnr2$allRRhigh

# Representación de los dos efectos

pnr1 <- data.frame(
  RR = cpallnr1$allRRfit,
  LowRR = cpallnr1$allRRlow,
  HighRR = cpallnr1$allRRhigh
)

pnr1 <- rownames_to_column(pnr1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pnr1)))

pnr2 <- data.frame(
  RR = cpallnr2$allRRfit,
  LowRR = cpallnr2$allRRlow,
  HighRR = cpallnr2$allRRhigh
)

pnr2 <- rownames_to_column(pnr2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pnr2)))

pnr <- rbind(pnr1, pnr2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pnr, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Narcea"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapnre = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca del Oriente

# Sin corregir

cpallo1 <- crosspred(
  blagj[[7]],
  coef = yall[7, ],
  vcov = Sall$oriente,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallo2 <- crosspred(
  blagj[[7]],
  coef = bl$oriente$blup,
  vcov = bl$oriente$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)
cpallo2$allRRlow
cpallo2$allRRfit
cpallo2$allRRhigh

# Representación de los dos efectos

po1 <- data.frame(
  RR = cpallo1$allRRfit,
  LowRR = cpallo1$allRRlow,
  HighRR = cpallo1$allRRhigh
)

po1 <- rownames_to_column(po1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(po1)))

po2 <- data.frame(
  RR = cpallo2$allRRfit,
  LowRR = cpallo2$allRRlow,
  HighRR = cpallo2$allRRhigh
)

po2 <- rownames_to_column(po2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(po2)))

po <- rbind(po1, po2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(po, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Oriente"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapoe = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Comarca de Oviedo

# Sin corregir

cpallov1 <- crosspred(
  blagj[[8]],
  coef = yall[8, ],
  vcov = Sall$oviedo,
  model.link = "log",
  at = -2:2,
  cen = 1
)

# BLUP

cpallov2 <- crosspred(
  blagj[[8]],
  coef = bl$oviedo$blup,
  vcov = bl$oviedo$vcov,
  model.link = "log",
  at = -2:2,
  cen = 1
)

cpallov2$allRRlow
cpallov2$allRRfit
cpallov2$allRRhigh

# Representación de los dos efectos

pov1 <- data.frame(
  RR = cpallov1$allRRfit,
  LowRR = cpallov1$allRRlow,
  HighRR = cpallov1$allRRhigh
)

pov1 <- rownames_to_column(pov1, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("Comarca", nrow(pov1)))

pov2 <- data.frame(
  RR = cpallov2$allRRfit,
  LowRR = cpallov2$allRRlow,
  HighRR = cpallov2$allRRhigh
)

pov2 <- rownames_to_column(pov2, "lag") %>% mutate_at(1, as.numeric) %>% mutate(grupo = rep("BLUP", nrow(pov2)))

pov <- rbind(pov1, pov2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pov, aes(x = lag, y = RR, color = grupo)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_line(size = 0.8, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = lag, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nLag",
    y = "RR\n",
    colour = "",
    title = "Oviedo"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shapove = NA))) +
  coord_cartesian(ylim = c(0, 2))

## Efecto del día de la semana

################################################################
# PRIMERA ETAPA

# Realizo la primera etapa del estudio comarca-específico.
# Consiste en ajustar un modelo para cada comarca y guardar
# los resultados relativos al día de la semana
################################################################

yall <- matrix(NA, length(ld), 6, dimnames = list(comarca, paste("b", seq(6), sep = "")))

Sall <- vector("list", length(ld))
names(Sall) <- comarca

# Objetos donde guardaré los IC95%

loweric <- vector("list", length(ld))
names(loweric) <- comarca
upperic <- vector("list", length(ld))
names(upperic) <- comarca

system.time({
  for(i in seq(ld)) {
    datos <- ld[[i]]
    
    suppressWarnings({
      datos <- datos %>%
        mutate(vacaciones = ifelse(mm == 12 & dd == 25, "Otras vacaciones", ifelse(mm == 1 & dd == 1, "Otras vacaciones", ifelse(mm == 1 & dd == 6|mm == 5 & dd == 1|mm == 9 & dd == 8|mm == 10 & dd == 12|mm == 11 & dd == 1|mm == 12 & dd == 6|mm == 12 & dd == 8, "Otras vacaciones", ifelse(mm == 7 & dd >= 15|mm == 8 & dd <= 15, "Verano", "No vacaciones")))))
      
      datos[c(106:107, 456:457, 813:814, 1198:1199, 1548:1549, 1933:1934, 2290:2291, 2647:2648, 3025:3026, 3382:3383, 3739:3740, 4117:4118, 4474:4475, 4859:4860, 5216:5217, 5566:5567, 5951:5952, 6308:6309, 6658:6659, 7043:7044, 7400:7401, 7750:7751, 8135:8136, 8492:8493, 8877:8878, 9227:9228, 9584:9585, 9969:9970, 10319:10320, 10676:10677, 11061:11062, 11411:11412, 11796:11797, 12153:12154, 12510:12511), "vacaciones"] <- "Otras vacaciones"
      
      datos$vacaciones <- as.factor(datos$vacaciones)
      datos$vacaciones <- relevel(datos$vacaciones, ref = "No vacaciones")
      
      datos$stratum <- as.factor(as.factor(datos$yy))
      ind <- tapply(datos$Total, datos$stratum, sum)[datos$stratum]
      
      mi.argvarj <- list(fun = "poly", degree = 2, int = F)
      
      nlagj <- 4
      lagnk <- 2; klagj <- logknots(nlagj,lagnk)
      mi.arglagj <- list(fun = "ns", knots = klagj, int = T)
      cbTj <- crossbasis(datos$tmed, lag = nlagj, argvar = mi.argvarj, arglag = mi.arglagj)
      
      vacan <- as.factor(as.numeric(datos$vacaciones))
      mi.arglagjv <- list(fun = "strata", breaks = c(0, 1))
      cbvacaj <- crossbasis(vacan, lag = c(-2, 2), argvar = list(fun = "integer", int = F), arglag = mi.arglagjv)
      
      dfseas <- 1
      ny <- length(unique(datos$yy))
    })
    
    datos <- as.data.frame(datos)
    mfirst <- gnm(Total ~ cbTj + cbvacaj + dow, data = datos, family = quasipoisson, eliminate = factor(stratum), na.action = "na.exclude", datosset = ind > 0)
    
    yall[i,] <- coef(mfirst)[15:20] 
    Sall[[i]] <- vcov(mfirst)[15:20, 15:20]
    # loweric[[i]] <- confint.gnm.diasemana(mfirst)[,1]
    # upperic[[i]] <- confint.gnm.diasemana(mfirst)[,2]
    
  }
})

# Guardo los IC95% para no tener que calcularlos todas las veces

# save(loweric, upperic, file = "./data/ICDiaComarcaGnm.RData")
load("./data/ICDiaComarcaGnm.RData")

################################################################
# SEGUNDA ETAPA

# Realizo la segunda etapa del estudio comarca-específico.
# Consiste en combinar los resultados de cada comarca mediante
# meta-análisis y obtener un único efecto general.
################################################################

# Ajuste del modelo meta-analítico

method <- "reml"
mvall <- mixmeta(yall ~ 1, Sall, method = method)
summary(mvall)

# Representación del efecto general combinado

spall <- predict(mvall, ci = T)

RR <- exp(spall$fit[1, ])
LowRR <- exp(spall$ci.lb[1, ])
HighRR <- exp(spall$ci.ub[1, ])

p.dia <- data.frame(RR = RR, LowRR = LowRR, HighRR = HighRR)
p.jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)
p.dia <- rbind(p.dia[2:4, ], p.jueves, p.dia[6, ], p.dia[5, ], p.dia[1, ])
p.dia <- p.dia %>%
  mutate(dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"))

ggplot(p.dia, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  group = 1
)) +
  geom_errorbar(aes(ymin = LowRR, ymax = HighRR, width = 0.2),
                colour = "purple",
                size = 0.8) +
  geom_point(colour = "purple4", size = 3) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "\nDía de la semana", y = "RR\n", title = "") +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

p.dia

################################################################
# BLUP

# Estimo las BLUPs de cada comarca en base al efecto estimado
# en cada comarca y al efecto general obtenido combinando los
# resultados
################################################################

bl <- blup(mvall, vcov = T, pi = T)

## Comarca de Avilés

# Sin corregir

RRaviles1 <- exp(yall)[1, ]
LowRRaviles1 <- exp(loweric$aviles)
HighRRaviles1 <- exp(upperic$aviles)

# BLUP

RRaviles2 <- exp(bl$aviles$blup)
LowRRaviles2 <- exp(bl$aviles$pi.lb)
HighRRaviles2 <- exp(bl$aviles$pi.ub)

# Representación de los dos efectos

paviles1 <- data.frame(RR = RRaviles1, LowRR = LowRRaviles1, HighRR = HighRRaviles1)

p.aviles1jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

paviles1 <- rbind(paviles1[2:4, ], p.aviles1jueves, paviles1[6, ], paviles1[5, ], paviles1[1, ])

paviles1 <- paviles1 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("Comarca", nrow(paviles1))
  )

paviles2 <- data.frame(RR = RRaviles2, LowRR = LowRRaviles2, HighRR = HighRRaviles2)

p.aviles2jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

paviles2 <- rbind(paviles2[2:4, ], p.aviles2jueves, paviles2[6, ], paviles2[5, ], paviles2[1, ])

paviles2 <- paviles2 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("BLUP", nrow(paviles2))
  )

paviles <- rbind(paviles1, paviles2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(paviles, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  color = grupo
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_point(aes(x = dia, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nDía de la semana",
    y = "RR\n",
    colour = "",
    title = "Avilés"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca del Caudal

# Sin corregir

RRcaudal1 <- exp(yall)[2, ]
LowRRcaudal1 <- exp(loweric$caudal)
HighRRcaudal1 <- exp(upperic$caudal)

# BLUP

RRcaudal2 <- exp(bl$caudal$blup)
LowRRcaudal2 <- exp(bl$caudal$pi.lb)
HighRRcaudal2 <- exp(bl$caudal$pi.ub)

# Representación de los dos efectos

pcaudal1 <- data.frame(RR = RRcaudal1, LowRR = LowRRcaudal1, HighRR = HighRRcaudal1)

p.caudal1jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

pcaudal1 <- rbind(pcaudal1[2:4, ], p.caudal1jueves, pcaudal1[6, ], pcaudal1[5, ], pcaudal1[1, ])

pcaudal1 <- pcaudal1 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("Comarca", nrow(pcaudal1))
  )

pcaudal2 <- data.frame(RR = RRcaudal2, LowRR = LowRRcaudal2, HighRR = HighRRcaudal2)

p.caudal2jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

pcaudal2 <- rbind(pcaudal2[2:4, ], p.caudal2jueves, pcaudal2[6, ], pcaudal2[5, ], pcaudal2[1, ])

pcaudal2 <- pcaudal2 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("BLUP", nrow(pcaudal2))
  )

pcaudal <- rbind(pcaudal1, pcaudal2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pcaudal, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  color = grupo
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_point(aes(x = dia, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nDía de la semana",
    y = "RR\n",
    colour = "",
    title = "Caudal"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca de Eo-Navia

# Sin corregir

RReonavia1 <- exp(yall)[3, ]
LowRReonavia1 <- exp(loweric$eonavia)
HighRReonavia1 <- exp(upperic$eonavia)

# BLUP

RReonavia2 <- exp(bl$eonavia$blup)
LowRReonavia2 <- exp(bl$eonavia$pi.lb)
HighRReonavia2 <- exp(bl$eonavia$pi.ub)

# Representación de los dos efectos

peonavia1 <- data.frame(RR = RReonavia1, LowRR = LowRReonavia1, HighRR =
                          HighRReonavia1)

p.eonavia1jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

peonavia1 <- rbind(peonavia1[2:4, ], p.eonavia1jueves, peonavia1[6, ], peonavia1[5, ], peonavia1[1, ])

peonavia1 <- peonavia1 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("Comarca", nrow(peonavia1))
  )

peonavia2 <- data.frame(RR = RReonavia2, LowRR = LowRReonavia2, HighRR =
                          HighRReonavia2)

p.eonavia2jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

peonavia2 <- rbind(peonavia2[2:4, ], p.eonavia2jueves, peonavia2[6, ], peonavia2[5, ], peonavia2[1, ])

peonavia2 <- peonavia2 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("BLUP", nrow(peonavia2))
  )

peonavia <- rbind(peonavia1, peonavia2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(peonavia, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  color = grupo
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_point(aes(x = dia, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nDía de la semana",
    y = "RR\n",
    colour = "",
    title = "Eo-Navia"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca de Gijón

# Sin corregir

RRgijon1 <- exp(yall)[4, ]
LowRRgijon1 <- exp(loweric$gijon)
HighRRgijon1 <- exp(upperic$gijon)

# BLUP

RRgijon2 <- exp(bl$gijon$blup)
LowRRgijon2 <- exp(bl$gijon$pi.lb)
HighRRgijon2 <- exp(bl$gijon$pi.ub)

# Representación de los dos efectos

pgijon1 <- data.frame(RR = RRgijon1, LowRR = LowRRgijon1, HighRR = HighRRgijon1)

p.gijon1jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

pgijon1 <- rbind(pgijon1[2:4, ], p.gijon1jueves, pgijon1[6, ], pgijon1[5, ], pgijon1[1, ])

pgijon1 <- pgijon1 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("Comarca", nrow(pgijon1))
  )

pgijon2 <- data.frame(RR = RRgijon2, LowRR = LowRRgijon2, HighRR = HighRRgijon2)

p.gijon2jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

pgijon2 <- rbind(pgijon2[2:4, ], p.gijon2jueves, pgijon2[6, ], pgijon2[5, ], pgijon2[1, ])

pgijon2 <- pgijon2 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("BLUP", nrow(pgijon2))
  )

pgijon <- rbind(pgijon1, pgijon2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pgijon, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  color = grupo
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_point(aes(x = dia, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nDía de la semana",
    y = "RR\n",
    colour = "",
    title = "Gijón"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca del Nalón

# Sin corregir

RRnalon1 <- exp(yall)[5, ]
LowRRnalon1 <- exp(loweric$nalon)
HighRRnalon1 <- exp(upperic$nalon)

# BLUP

RRnalon2 <- exp(bl$nalon$blup)
LowRRnalon2 <- exp(bl$nalon$pi.lb)
HighRRnalon2 <- exp(bl$nalon$pi.ub)

# Representación de los dos efectos

pnalon1 <- data.frame(RR = RRnalon1, LowRR = LowRRnalon1, HighRR = HighRRnalon1)

p.nalon1jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

pnalon1 <- rbind(pnalon1[2:4, ], p.nalon1jueves, pnalon1[6, ], pnalon1[5, ], pnalon1[1, ])

pnalon1 <- pnalon1 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("Comarca", nrow(pnalon1))
  )

pnalon2 <- data.frame(RR = RRnalon2, LowRR = LowRRnalon2, HighRR = HighRRnalon2)

p.nalon2jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

pnalon2 <- rbind(pnalon2[2:4, ], p.nalon2jueves, pnalon2[6, ], pnalon2[5, ], pnalon2[1, ])

pnalon2 <- pnalon2 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("BLUP", nrow(pnalon2))
  )

pnalon <- rbind(pnalon1, pnalon2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pnalon, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  color = grupo
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_point(aes(x = dia, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nDía de la semana",
    y = "RR\n",
    colour = "",
    title = "Nalón"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca de Narcea

# Sin corregir

RRnarcea1 <- exp(yall)[6, ]
LowRRnarcea1 <- exp(loweric$narcea)
HighRRnarcea1 <- exp(upperic$narcea)

# BLUP

RRnarcea2 <- exp(bl$narcea$blup)
LowRRnarcea2 <- exp(bl$narcea$pi.lb)
HighRRnarcea2 <- exp(bl$narcea$pi.ub)

# Representación de los dos efectos

pnarcea1 <- data.frame(RR = RRnarcea1, LowRR = LowRRnarcea1, HighRR = HighRRnarcea1)

p.narcea1jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

pnarcea1 <- rbind(pnarcea1[2:4, ], p.narcea1jueves, pnarcea1[6, ], pnarcea1[5, ], pnarcea1[1, ])

pnarcea1 <- pnarcea1 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("Comarca", nrow(pnarcea1))
  )

pnarcea2 <- data.frame(RR = RRnarcea2, LowRR = LowRRnarcea2, HighRR = HighRRnarcea2)

p.narcea2jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

pnarcea2 <- rbind(pnarcea2[2:4, ], p.narcea2jueves, pnarcea2[6, ], pnarcea2[5, ], pnarcea2[1, ])

pnarcea2 <- pnarcea2 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("BLUP", nrow(pnarcea2))
  )
pnarcea <- rbind(pnarcea1, pnarcea2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(pnarcea, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  color = grupo
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_point(aes(x = dia, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nDía de la semana",
    y = "RR\n",
    colour = "",
    title = "Narcea"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca del Oriente

# Sin corregir

RRoriente1 <- exp(yall)[7, ]
LowRRoriente1 <- exp(loweric$oriente)
HighRRoriente1 <- exp(upperic$oriente)

# BLUP

RRoriente2 <- exp(bl$oriente$blup)
LowRRoriente2 <- exp(bl$oriente$pi.lb)
HighRRoriente2 <- exp(bl$oriente$pi.ub)

# Representación de los dos efectos

poriente1 <- data.frame(RR = RRoriente1, LowRR = LowRRoriente1, HighRR =
                          HighRRoriente1)

p.oriente1jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

poriente1 <- rbind(poriente1[2:4, ], p.oriente1jueves, poriente1[6, ], poriente1[5, ], poriente1[1, ])

poriente1 <- poriente1 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("Comarca", nrow(poriente1))
  )

poriente2 <- data.frame(RR = RRoriente2, LowRR = LowRRoriente2, HighRR =
                          HighRRoriente2)

p.oriente2jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

poriente2 <- rbind(poriente2[2:4, ], p.oriente2jueves, poriente2[6, ], poriente2[5, ], poriente2[1, ])

poriente2 <- poriente2 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("BLUP", nrow(poriente2))
  )

poriente <- rbind(poriente1, poriente2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(poriente, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  color = grupo
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_point(aes(x = dia, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nDía de la semana",
    y = "RR\n",
    colour = "",
    title = "Oriente"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))

## Comarca de Oviedo

# Sin corregir

RRoviedo1 <- exp(yall)[8, ]
LowRRoviedo1 <- exp(loweric$oviedo)
HighRRoviedo1 <- exp(upperic$oviedo)

# BLUP

RRoviedo2 <- exp(bl$oviedo$blup)
LowRRoviedo2 <- exp(bl$oviedo$pi.lb)
HighRRoviedo2 <- exp(bl$oviedo$pi.ub)

# Representación de los dos efectos

poviedo1 <- data.frame(RR = RRoviedo1, LowRR = LowRRoviedo1, HighRR = HighRRoviedo1)

p.oviedo1jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

poviedo1 <- rbind(poviedo1[2:4, ], p.oviedo1jueves, poviedo1[6, ], poviedo1[5, ], poviedo1[1, ])

poviedo1 <- poviedo1 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("Comarca", nrow(poviedo1))
  )

poviedo2 <- data.frame(RR = RRoviedo2, LowRR = LowRRoviedo2, HighRR = HighRRoviedo2)

p.oviedo2jueves <- data.frame(RR = 1, LowRR = NA, HighRR = NA)

poviedo2 <- rbind(poviedo2[2:4, ], p.oviedo2jueves, poviedo2[6, ], poviedo2[5, ], poviedo2[1, ])

poviedo2 <- poviedo2 %>%
  mutate(
    dia = c("LU", "MA", "MI", "JU", "VI", "SA", "DO"),
    grupo = rep("BLUP", nrow(poviedo2))
  )

poviedo <- rbind(poviedo1, poviedo2)

colors <- c("Comarca" = "olivedrab3",
            "BLUP" = "royalblue3")

ggplot(poviedo, aes(
  x = factor(dia, levels = c("LU", "MA", "MI", "JU", "VI", "SA", "DO")),
  y = RR,
  color = grupo
)) +
  geom_errorbar(
    aes(ymin = LowRR, ymax = HighRR, width = 0.2),
    size = 0.8,
    position = position_dodge(width = 0.2)
  ) +
  geom_point(aes(x = dia, y = RR),
             size = 3,
             position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 1,
             size = 0.5,
             linetype = "dashed") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    x = "\nDía de la semana",
    y = "RR\n",
    colour = "",
    title = "Oviedo"
  ) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(shape = NA)))
