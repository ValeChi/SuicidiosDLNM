packages = c('readxl', 'mapSpain', 'tidyverse', 'lubridate', 'ggplot2',  'gridExtra', 'ggrepel', "stringr", "climaemet", "giscoR", "mice", "naniar", "maps", "sf", "units", "rnaturalearth", "forcats")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos='http://cran.rediris.es')
    library(x, character.only = TRUE)
  }
})

################################################################
# EXPORTACIÓN DE DATOS DE TEMPERATURA

# Exporto los datos meteorológicos de la AEMET y los guardo 
# en un archivo RData para no tener que repetir el proceso.
################################################################


## Obtener api key de la AEMET
#
#browseURL("https://opendata.aemet.es/centrodedescargas/obtencionAPIKey")
# 
## Acceder a los datos con la api key obtenida
#
# aemet_api_key("MiApiKey",overwrite=T,install=T)
# 
## Guardo los datos de todas las estaciones meteorológicas
#
# stations <- aemet_stations() 
# 
## Selecciono las estaciones de Asturias
#
# estaciones <- stations %>% filter(provincia == "ASTURIAS")
#
## Creo una lista donde cada elemento es una de las estaciones de Asturias y, para cada una,
## exporto los datos en el periodo de tiempo estudiado, es decir de enero de 1987 a diciembre de 2021.
#
# StatCod <- unique(estaciones$indicativo)
# 
# dlist <- vector("list", length(StatCod))
#
# names(dlist) <- StatCod
#
# for(i in 1:length(StatCod)){
# dlist[[i]] <- aemet_daily_clim(station = StatCod[i], start = "1987-01-01", end = "2021-12-31")
# }
#
## Guardo la lista con los datos de cada una de las estaciones para el periodo estudiado 
## y guardo también los datos informativos de cada una de las estaciones (latitud, longitud, nombre, etc.)
#
# save(dlist, file = "./data/TemperaturaAsturias.RData")
# save(estaciones, file = "./data/EstacionesAsturias.RData")


################################################################
# ESTUDIO DE DATOS DE TEMPERATURA

# Miro la cantidad de datos faltantes para la variable temperatura
# y represento las estaciones en el mapa de Asturias
################################################################


# Importo los datos que me bajé desde la AEMET y que me guardé

load("./data/TemperaturaAsturias.RData")
load("./data/EstacionesAsturias.RData")

# Miro el porcentaje de NA de temperatura por cada estación, el criterio
# de inclusión es presentar menos del 60% de datos faltantes

cat("El porcentaje de NA de para cada estación es:", "\n")
for (i in 1:13) {
  if (nrow(dlist[[i]]) > 0) {
    db <- as.data.frame(dlist[[i]])
    ts <- seq.POSIXt(as.POSIXct("1987-01-02"), as.POSIXct("2022-01-01"), by = "day")
    ts <- data.frame(fecha = ts)
    ts$fecha <- as.Date(ts$fecha)
    db <- ts %>%
      left_join(db, by = "fecha")
    if (dim(prop.table(table(is.na(db$tmed)))) > 1) {
      cat(dlist[[i]]$indicativo[1],
          dlist[[i]]$nombre[1],
          ":",
          round(prop.table(table(
            is.na(db$tmed)
          ))[[2]] * 100, 2),
          "\n")
    } else {
      cat(
        "La estación",
        dlist[[i]]$indicativo[1],
        dlist[[i]]$nombre[1],
        "presenta todos los datos",
        "\n"
      )
    }
  } else {
    cat("La estación número", i, "no presenta ningún dato", "\n")
  }
}

# Miro la correlación con la temperatura registrada en Oviedo por la estación
# 1249I, la única que presenta todos los datos

cat("El coeficiente de Pearson y su IC del 95% para cada estación es;",
    "\n")
for (i in 1:13) {
  if (nrow(dlist[[i]]) > 0) {
    db <- as.data.frame(dlist[[i]])
    ts <- seq.POSIXt(as.POSIXct("1987-01-02"), as.POSIXct("2022-01-01"), by = "day")
    ts <- data.frame(fecha = ts)
    ts$fecha <- as.Date(ts$fecha)
    db <- ts %>%
      left_join(db, by = "fecha")
    if (i != 10) {
      cat(
        dlist[[i]]$indicativo[1],
        dlist[[i]]$nombre[1],
        ":",
        round(cor.test(dlist[[10]]$tmed, db$tmed)$estimate, 3),
        "(",
        round(cor.test(dlist[[10]]$tmed, db$tmed)$conf.int[1], 3),
        "-",
        round(cor.test(dlist[[10]]$tmed, db$tmed)$conf.int[2], 3),
        ")\n"
      )
    }
    else {
      NULL
    }
  }
}


# Ubico todas las estaciones en el mapa de Asturias

## Creo los nombres para las etiquetas de la figura

estaciones$estacion <- paste(estaciones$indicativo, estaciones$nombre, sep = ": ")

## Obtengo la cartografía de Asturias a partir del paquete mapSpain

Asturias <- esp_get_munic_siane(region = "Asturias") %>%
  mutate(Provincia = esp_dict_translate(ine.prov.name, "es"))

## Represento el mapa de Asturias con la ubicación de todas las estaciones

ggplot(Asturias) +
  geom_sf(aes(fill = Provincia), color = "grey70") +
  geom_point(
    data = estaciones,
    aes(x = longitud, y = latitud),
    color = 'black',
    fill = "red",
    size = 5,
    shape = 21
  ) +
  labs(title = "Estaciones meteorológicas de Asturias") +
  scale_fill_discrete(type =
                        hcl.colors(4, "Greens")) +
  theme(legend.position = 'none') +
  xlab("Longitud") +
  ylab("Latitud") +
  geom_label_repel(
    data = estaciones,
    aes(x = longitud, y = latitud, label = estacion),
    box.padding   = 0.40,
    point.padding = 0.5,
    segment.color = 'black',
    size = 1.5,
    min.segment.length = 0.1,
    nudge_y = 0.10
  )

################################################################
# IMPUTACIÓN DE DATOS DE TEMPERATURA

# Imputo los datos faltantes de las estaciones seleccionadas
# para poder trabajar con datos completos.
################################################################

## Para cada estación selecciono las variables que quiero usar como predictores para imputar los datos faltantes
# 
## ESTACIÓN 1212E, Aeropuerto de Asturias
# 
# e.1212e <- as.data.frame(dlist$"1212E")
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1212e <- fecha %>%
#   left_join(e.1212e)
# e.1212e <- e.1212e %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia, sol)
# colnames(e.1212e)[-1] <- paste(colnames(e.1212e)[-1], "Aeropuerto", sep = "_")
# 
## ESTACIÓN 1283U, Cabo Busto
# 
# e.1283u <- as.data.frame(dlist$"1283U") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1283u <- fecha %>%
#   left_join(e.1283u)
# e.1283u <- e.1283u %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia, sol)
# colnames(e.1283u)[-1] <- paste(colnames(e.1283u)[-1], "CaboBusto", sep = "_")
# 
## ESTACIÓN 1210X, Cabo Peñas
# 
# e.1210x <- as.data.frame(dlist$"1210X") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1210x <- fecha %>%
#   left_join(e.1210x)
# e.1210x <- e.1210x %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia)
# colnames(e.1210x)[-1] <- paste(colnames(e.1210x)[-1], "CaboPeñas", sep = "_")
# 
## ESTACIÓN 1207U, Gijón Campus
# 
# e.1207u <- as.data.frame(dlist$"1207U") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1207u <- fecha %>%
#   left_join(e.1207u)
# e.1207u <- e.1207u %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia, sol)
# colnames(e.1207u)[-1] <- paste(colnames(e.1207u)[-1], "GijonCampus", sep = "_")
# 
## ESTACIÓN 1208H, Gijón Puerto
# 
# e.1208h <- as.data.frame(dlist$"1208H") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1208h <- fecha %>%
#   left_join(e.1208h)
# e.1208h <- e.1208h %>%
#   dplyr::select(fecha, tmed, tmin, tmax, sol) #la estación de Gijón puerto no registra la velocidad media del viento
# colnames(e.1208h)[-1] <- paste(colnames(e.1208h)[-1], "GijonPuerto", sep = "_")
# 
## ESTACIÓN 1208, Gijón 
# 
# e.1208 <- as.data.frame(dlist$"1208") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1208 <- fecha %>%
#   left_join(e.1208)
# e.1208 <- e.1208 %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia, sol)
# colnames(e.1208)[-1] <- paste(colnames(e.1208)[-1], "Gijon", sep = "_")
# 
## ESTACIÓN 1183X, Llanes
# 
# e.1183x <- as.data.frame(dlist$"1183X") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1183x <- fecha %>%
#   left_join(e.1183x)
# e.1183x <- e.1183x %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia, sol)
# colnames(e.1183x)[-1] <- paste(colnames(e.1183x)[-1], "Llanes", sep = "_")
# 
## ESTACIÓN 1249I, Oviedo
# 
# e.1249i <- as.data.frame(dlist$"1249I") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1249i <- fecha %>%
#   left_join(e.1249i)
# e.1249i <- e.1249i %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia, sol)
# colnames(e.1249i)[-1] <- paste(colnames(e.1249i)[-1], "Oviedo", sep = "_")
# 
## ESTACIÓN 1221D, Pajares-Valgrande
# 
# e.1221d <- as.data.frame(dlist$"1221D") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1221d <- fecha %>%
#   left_join(e.1221d)
# e.1221d <- e.1221d %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia)
# colnames(e.1221d)[-1] <- paste(colnames(e.1221d)[-1], "PajaresValgrande", sep = "_")
# 
## ESTACIÓN 1542, Puerto de Leitariegos
# 
# e.1542 <- as.data.frame(dlist$"1542") 
# fecha <- data.frame(fecha = ts)
# fecha$fecha <- as.Date(fecha$fecha)
# e.1542 <- fecha %>%
#   left_join(e.1542)
# e.1542 <- e.1542 %>%
#   dplyr::select(fecha, tmed, tmin, tmax, velmedia, sol)
# colnames(e.1542)[-1] <- paste(colnames(e.1542)[-1], "PuertoDeLeitariegos", sep = "_")
# 
## Junto las variables de todas las estaciones y a partir de esa base de datos voy a imputar los datos faltantes
# 
# datos.imputacion <- e.1212e %>%
#   left_join(e.1283u) %>%
#   left_join(e.1210x) %>%
#   left_join(e.1207u) %>%
#   left_join(e.1208h) %>%
#   left_join(e.1208) %>%
#   left_join(e.1183x) %>%
#   left_join(e.1249i) %>%
#   left_join(e.1221d) %>%
#   left_join(e.1542)
# 
## Extraigo el año y el mes de la fecha para poder tener en cuenta estacionalidad y tendencia de la serie
## a la hora de realizar la imputación múltiple
# 
# datos.imputacion$yy <- year(datos.imputacion$fecha)
# datos.imputacion$mm <- month(datos.imputacion$fecha)
# 
## Gráfico de valores faltantes de temperatura
# 
# gg_miss_var(datos.imputacion[, c(2, 7, 12, 16, 21, 25, 30, 35, 40, 44)]) + theme_bw()
# 
## Imputo los datos usando el método de la regresión lineal con 30 iteraciones y creando 5 bases de datos completas
# 
# set.seed(123)
# datos.imputacion2 <- mice(datos.imputacion[, -c(1, 3, 4, 8, 9, 13, 14, 17, 18, 22, 23, 26, 27, 31, 32, 36, 37, 41, 42, 45, 46)], method = "norm", maxit = 30, m = 5)
# 
# datos.imputacion2$predictorMatrix
# #plot(datos.imputacion2)
# 
## Guardo en objetos separados las 5 bases de datos finales obtenidas
#
# for(i in 1:5) {
# imputed <- paste("imputed", i, sep = ".")
# assign(imputed, as.data.frame(complete(datos.imputacion2, action = i)[ , c(1, 4, 7, 9, 12, 14, 17, 20, 23, 25)]))
# }
# 
## Selecciono las 5 mediciones de cada estación y obtengo una final calculando la media
#
# for(i in 1:10) {
#   imputed.1m <- imputed.1 %>%
#     select(i)
#   imputed.2m <- imputed.2 %>%
#     select(i)
#   imputed.3m <- imputed.3 %>%
#     select(i)
#   imputed.4m <- imputed.4 %>%
#     select(i)
#   imputed.5m <- imputed.5 %>%
#     select(i)
#   db <- cbind(imputed.1m, imputed.2m, imputed.3m, imputed.4m, imputed.5m)
#   nombre <- paste("tmed.imputada.estacion", i, sep = "")
#   assign(nombre, rowMeans(db))
# 
# }
# 
## Creo una base de datos con la temperatura observada y la imputada para cada estación y guardo la base
## para no tener que repetir el proceso todas las veces
#
# datos.imputacion.final <- data.frame(fecha = datos.imputacion$fecha, tmed_Aeropuerto = datos.imputacion$tmed_Aeropuerto, tmed_Aeropuerto_i = tmed.imputada.estacion1, tmed_CaboBusto = datos.imputacion$tmed_CaboBusto, tmed_CaboBusto_i = tmed.imputada.estacion2, tmed_CaboPeñas = datos.imputacion$tmed_CaboPeñas, tmed_CaboPeñas_i = tmed.imputada.estacion3, tmed_GijonCampus = datos.imputacion$tmed_GijonCampus, tmed_GijonCampus_i = tmed.imputada.estacion4, tmed_GijonPuerto = datos.imputacion$tmed_GijonPuerto, tmed_GijonPuerto_i = tmed.imputada.estacion5, tmed_Gijon = datos.imputacion$tmed_Gijon, tmed_Gijon_i = tmed.imputada.estacion6, tmed_Llanes = datos.imputacion$tmed_Llanes, tmed_Llanes_i = tmed.imputada.estacion7, tmed_Oviedo = datos.imputacion$tmed_Oviedo, tmed_Oviedo_i = tmed.imputada.estacion8, tmed_PajaresValgrande = datos.imputacion$tmed_PajaresValgrande, tmed_PajaresValgrande_i = tmed.imputada.estacion9, tmed_PuertoDeLeitariegos = datos.imputacion$tmed_PuertoDeLeitariegos, tmed_PuertoDeLeitariegos_i = tmed.imputada.estacion10)
# 
# save(datos.imputacion.final, file="./data/datos.imputacion.final.RData")

# Importo los datos imputados de cada estación

load("./data/datos.imputacion.final.RData")

# Represento la serie de datos observados y la de datos imputados para ver si la imputación se puede considerar buena

## Cambio el nombre de las variables para poder cambiar de formato wide a long y hacer la gráfica

datos.imputacion.final <- datos.imputacion.final %>%
  rename(tmed.Aeropuerto = tmed_Aeropuerto,
         tmed_imputada.Aeropuerto = tmed_Aeropuerto_i,
         tmed.CaboBusto = tmed_CaboBusto,
         tmed_imputada.CaboBusto = tmed_CaboBusto_i,
         tmed.CaboPeñas = tmed_CaboPeñas,
         tmed_imputada.CaboPeñas = tmed_CaboPeñas_i,
         tmed.GijonCampus = tmed_GijonCampus,
         tmed_imputada.GijonCampus = tmed_GijonCampus_i,
         tmed.GijonPuerto = tmed_GijonPuerto,
         tmed_imputada.GijonPuerto = tmed_GijonPuerto_i,
         tmed.Gijon = tmed_Gijon,
         tmed_imputada.Gijon = tmed_Gijon_i,
         tmed.Llanes = tmed_Llanes,
         tmed_imputada.Llanes = tmed_Llanes_i,
         tmed.Oviedo = tmed_Oviedo,
         tmed_imputada.Oviedo = tmed_Oviedo_i,
         tmed.PajaresValgrande = tmed_PajaresValgrande,
         tmed_imputada.PajaresValgrande = tmed_PajaresValgrande_i,
         tmed.PuertoDeLeitariegos = tmed_PuertoDeLeitariegos,
         tmed_imputada.PuertoDeLeitariegos = tmed_PuertoDeLeitariegos_i)

## Paso al formato long

long.imputados <- reshape(datos.imputacion.final,
                          varying = c(2:21),
                          direction = "long")

## Represento las series de valores observados e imputados para cada estación

ggplot(long.imputados) +
  geom_line(aes(fecha, tmed), col = "red", alpha = 0.6) +
  geom_line(aes(fecha, tmed_imputada),
            col = "blue",
            alpha = 0.4) +
  ylim(-15, 30) +
  facet_wrap(time ~ ., scales = "free") +
  ggtitle("Series temporales de temperatura: valores observados (rojo) y valores imputados (azul)") +
  theme_minimal() +
  xlab("Fecha") +
  ylab("Temperatura media (ºC)") #La serie de datos imputados se parece a la de datos observados, la variabilidad está bien y podemos decir que los datos parecen bien imputados.

# Represento las funciones de densidad de los datos observados y de los imputados. Si la imputación es buena, las densidades deberían tener una forma parecida.

long.imputados <- long.imputados %>%
  rename(
    tmed.Observados = tmed,
    tmed.Imputados = tmed_imputada,
    estación = time
  ) %>%
  dplyr::select(-id)

long.imputados2 <- reshape(long.imputados,
                           varying = c(3, 4),
                           direction = "long")

dplot <- long.imputados2 %>% mutate(Datos = factor(time))

ggplot(dplot) +
  geom_density(aes(x = tmed, fill = Datos, col = Datos),
               alpha = 0.2,
               size = 0.5) +
  facet_wrap(estación ~ ., scales = "free", nrow = 5) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Figura 2: gráficos de densidad para la temperatura media") +
  labs(x = "Temperatura media", y = "Densidad") +
  theme_minimal() #Las densidades se parecen, la imputación es buena.

################################################################
# CÁLCULO DE LAS DISTANCIAS

# Calculo las distancias entre las estaciones meteorológicas y 
# cada uno de los municipios de Asturias
################################################################

# Importo los datos de la población de cada municipio que me harán falta después

load("./data/pob_media.RData")

pob <- as.data.frame(pob)

# Junto los datos de población con los demás datos de los municipios

Asturias <- Asturias %>%
  mutate(cmun = str_remove(cmun, "^0+"))
Asturias <- Asturias %>%
  left_join(pob) 

# Calculo la distancia entre cada municipio y cada estación 

estaciones.sf <- st_as_sf(estaciones,
                          coords = c("longitud", "latitud"),
                          crs = 4258)

dist.Aeropuerto <- st_distance(slice(estaciones.sf, 1), Asturias)
dist.CaboBusto <- st_distance(slice(estaciones.sf, 2), Asturias)
dist.CaboPeñas <- st_distance(slice(estaciones.sf, 3), Asturias)
dist.GijonCampus <- st_distance(slice(estaciones.sf, 5), Asturias)
dist.GijonPuerto <- st_distance(slice(estaciones.sf, 7), Asturias)
dist.Gijon <- st_distance(slice(estaciones.sf, 8), Asturias)
dist.Llanes <- st_distance(slice(estaciones.sf, 9), Asturias)
dist.Oviedo <- st_distance(slice(estaciones.sf, 10), Asturias)
dist.PajaresValgrande <- st_distance(slice(estaciones.sf, 12), Asturias)
dist.PuertoDeLeitariegos <- st_distance(slice(estaciones.sf, 13), Asturias)

# Convierto la distancia en km

Asturias <- Asturias %>%
  mutate(
    dist.Aeropuerto = as.vector(set_units(dist.Aeropuerto, "km")),
    dist.CaboBusto = as.vector(set_units(dist.CaboBusto, "km")),
    dist.CaboPeñas = as.vector(set_units(dist.CaboPeñas, "km")),
    dist.GijonCampus = as.vector(set_units(dist.GijonCampus, "km")),
    dist.GijonPuerto = as.vector(set_units(dist.GijonPuerto, "km")),
    dist.Gijon = as.vector(set_units(dist.Gijon, "km")),
    dist.Llanes = as.vector(set_units(dist.Llanes, "km")),
    dist.Oviedo = as.vector(set_units(dist.Oviedo, "km")),
    dist.PajaresValgrande = as.vector(set_units(dist.PajaresValgrande, "km")),
    dist.PuertoDeLeitariegos = as.vector(set_units(dist.PuertoDeLeitariegos, "km"))
  )

# Para los municipios en los que se encuentran las estaciones, R clcula una distancia de 0 km y eso llevaría a tener unos valores infinitos en los pasos siguientes. Por esa razón a las distancias iguales a 0 asigno un valor muy muy pequeño (0.0000001)

Asturias <- Asturias %>%
  mutate(
    dist.Aeropuerto = ifelse(dist.Aeropuerto == 0, 0.0000001, dist.Aeropuerto),
    dist.GijonCampus = ifelse(dist.GijonCampus == 0, 0.0000001, dist.GijonCampus),
    dist.Gijon = ifelse(dist.Gijon == 0, 0.0000001, dist.Gijon),
    dist.CaboPeñas = ifelse(dist.CaboPeñas == 0, 0.0000001, dist.CaboPeñas),
    dist.CaboBusto = ifelse(dist.CaboBusto == 0, 0.0000001, dist.CaboBusto),
    dist.PajaresValgrande = ifelse(dist.PajaresValgrande == 0, 0.0000001, dist.PajaresValgrande),
    dist.Oviedo = ifelse(dist.Oviedo == 0, 0.0000001, dist.Oviedo),
  )

# Calculo el radio, es decir la máxima distancia que consideraré para asignar las estaciones de referencia a los municipios

# Extraigo la distancia mínima asociada a cada municipio

Asturias <- Asturias %>%
  mutate(dist.min = apply(as.data.frame(Asturias)[, 11:20], 1, FUN = min))

# El máximo de las distancias mínimas es mi radio (53.44 km)

radio <- max(Asturias$dist.min)

# Asigno a cada municipio las estaciones que se encuentran dentro del radio de distancia

Asturias <- Asturias %>%
  mutate(
    Aeropuerto = ifelse(dist.Aeropuerto <= radio, 1, 0),
    CaboBusto = ifelse(dist.CaboBusto <= radio, 1, 0),
    CaboPeñas = ifelse(dist.CaboPeñas <= radio, 1, 0),
    GijonCampus = ifelse(dist.GijonCampus <= radio, 1, 0),
    GijonPuerto = ifelse(dist.GijonPuerto <= radio, 1, 0),
    Gijon = ifelse(dist.Gijon <= radio, 1, 0),
    Llanes = ifelse(dist.Llanes <= radio, 1, 0),
    Oviedo = ifelse(dist.Oviedo <= radio, 1, 0),
    PajaresValgrande = ifelse(dist.PajaresValgrande <= radio, 1, 0),
    PuertoDeLeitariegos = ifelse(dist.PuertoDeLeitariegos <= radio, 1, 0)
  )

################################################################
# CÁLCULO DE LA TEMPERATURA MEDIA POR MUNICIPIO

# Estimo la temperatura media de cada municipio según el
# método de la distancia inversa ponderada
################################################################

# Para cada municipio, cojo la temperatura imputada de cada estación y la divido por la distancia observada entre ese municipio y esa estación. Obtengo 78 bases de datos con la temperatura de las 10 estaciones dividida por la distancia entre el municipio y la estación.

dlist <- vector("list", length(Asturias$cmun))
names(dlist) <- Asturias$cmun
for (i in 1:length(Asturias$cmun)) {
  dlist[[i]] <- datos.imputacion.final[, c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21)]
  for (j in 1:10) {
    dlist[[i]][j] <- dlist[[i]][j] / as.data.frame(Asturias)[i, 10 + j]
  }
}

# Para cada municipio, elimino las columnas de temperatura de las estaciones que no le han sido asignadas. Al final me quedo con 78 bases de datos, cada una corresponde a un municipio y contiene la temperatura imputada registrada por las estaciones que le han sido asignadas, dividida por la distancia observada entre el municipio y esa estación.

dlist2 <- vector("list", length(Asturias$cmun))
names(dlist2) <- Asturias$cmun
for (i in 1:78) {
  for (j in 1:10) {
    if (as.data.frame(Asturias)[i, 21 + j] == 1) {
      dlist2[[i]][j] <- dlist[[i]][j]
    } else {
      dlist2[[i]][j] <- NA
    }
  }
}

dlist3 <- lapply(dlist2, function(x)
  x[!is.na(x)])

# Para cada municipio, sumo los datos de las estaciones de referencia (el numerador de la ecuación final)

dlist4 <- vector("list", length(Asturias$cmun))
names(dlist4) <- Asturias$cmun
for (i in 1:78) {
  dlist4[[i]] <- Reduce("+", dlist3[[i]])
  
}

temp.municipios <- as.data.frame(dlist4)

# Para cada municipio, calculo la inversa de la distancia entre el municipio y cada estación de referencia y sumo los valores (denominador de la ecuación)

dlist5 <- vector("list", length(Asturias$cmun))
names(dlist5) <- Asturias$cmun
for (i in 1:78) {
  for (j in 1:10) {
    if (as.data.frame(Asturias)[i, 21 + j] == 1) {
      dlist5[[i]][j] <- 1 / as.data.frame(Asturias)[i, 10 + j]
    } else {
      dlist5[[i]][j] <- NA
    }
  }
}

dlist5 <- lapply(dlist5, function(x)
  x[!is.na(x)])

dlist6 <- vector("list", length(Asturias$cmun))

names(dlist6) <- Asturias$cmun

for (i in 1:78) {
  dlist6[[i]] <- Reduce("+", dlist5[[i]])
  
}

dist.municipios <- unlist(dlist6)

# Obtengo la temperatura media de cada municipio

temp.municipios.final <- sweep(temp.municipios, 2, dist.municipios, FUN = '/')

# Represento la serie de temperatura por municipio

temp.municipios.final2 <- temp.municipios.final
temp.municipios.final2$fecha <- datos.imputacion.final$fecha

names(temp.municipios.final2) <- c(
  "tmed.1Allande",
  "tmed.2Aller",
  "tmed.3Amieva",
  "tmed.4Avilés",
  "tmed.5BelmonteDeMiranda",
  "tmed.6Bimenes",
  "tmed.7Boal",
  "tmed.8Cabrales",
  "tmed.9Cabranes",
  "tmed.10Candamo",
  "tmed.11CangasDelNarcea",
  "tmed.12CangasDeOnís",
  "tmed.13Caravia",
  "tmed.14Carreño",
  "tmed.15Caso",
  "tmed.16Castrillón",
  "tmed.17Castropol",
  "tmed.18Cuaña",
  "tmed.19Colunga",
  "tmed.20CorveraDeAsturias",
  "tmed.21Cudillero",
  "tmed.22Degaña",
  "tmed.23ElFranco",
  "tmed.24Gijón",
  "tmed.25Gozón",
  "tmed.26Grado",
  "tmed.27GrandasDeSalime",
  "tmed.28Ibias",
  "tmed.29Eilao",
  "tmed.30Illas",
  "tmed.31Langreo",
  "tmed.32Laviana",
  "tmed.33Lena",
  "tmed.34Valdés",
  "tmed.35Llanera",
  "tmed.36Llanes",
  "tmed.37Mieres",
  "tmed.38Morcín",
  "tmed.39MurosDeNalón",
  "tmed.40Nava",
  "tmed.41Navia",
  "tmed.42Noreña",
  "tmed.43Onís",
  "tmed.44Oviedo",
  "tmed.45Parres",
  "tmed.46PeñamelleraAlta",
  "tmed.47PeñamelleraBaja",
  "tmed.48Pesoz",
  "tmed.49Piloña",
  "tmed.50Ponga",
  "tmed.51Pravia",
  "tmed.52Proaza",
  "tmed.53Quirós",
  "tmed.54LesRegueres",
  "tmed.55Ribadedeva",
  "tmed.56Ribadesella",
  "tmed.57RiberaDeArriba",
  "tmed.58Riosa",
  "tmed.59Salas",
  "tmed.60SanMartínDelReyAurelio",
  "tmed.61SanMartínDeOscos",
  "tmed.62SantaEulaliaDeOscos",
  "tmed.63SanTirsoDeAbres",
  "tmed.64SantoAdriano",
  "tmed.65Sariego",
  "tmed.66Siero",
  "tmed.67Sobrescobio",
  "tmed.68Somiedo",
  "tmed.69SotoDelBarco",
  "tmed.70TapiaDeCasariego",
  "tmed.71Taramundi",
  "tmed.72Teverga",
  "tmed.73Tineo",
  "tmed.74Vegadeo",
  "tmed.75VillanuevaDeOscos",
  "tmed.76Villaviciosa",
  "tmed.77Villayón",
  "tmed.78YernesYTameza",
  "fecha"
)

long.tmunic <- reshape(temp.municipios.final2,
                       varying = c(1:78),
                       direction = "long")

long.tmunic$time <- as.factor(long.tmunic$time)

long.tmunic$time <- fct_relevel(long.tmunic$time, gtools::mixedsort(levels(long.tmunic$time)))

ggplot(long.tmunic) +
  geom_line(aes(fecha, tmed), col = "red", alpha = 0.6) +
  ylim(-15, 30) +
  ggforce::facet_wrap_paginate(
    ~ time,
    scales = "free",
    nrow = 6,
    ncol = 7,
    page = 2
  ) + # distribuyo las gráficas en 2 páginas (son 78 series)
  ggtitle("Serie temporal de temperatura media por municipios") +
  theme_minimal() +
  xlab("Fecha") +
  ylab("Temperatura media (ºC)")

rm(list = ls()[!ls() %in% c(
  "Asturias",
  "datos.imputacion.final",
  "estaciones_sf",
  "poblacion",
  "temp.municipios.final",
  "radio"
)])

# Guardo la base con la temperatura media registrada en cada municipio

temp.municipios.final$fecha <- datos.imputacion.final$fecha
temp.municipios.final <- temp.municipios.final %>%
  dplyr::select(fecha, everything())

#save(temp.municipios.final, file = "data/TmedMunicipios.RData")

################################################################
# CÁLCULO DE LA TEMPERATURA GLOBAL

# Estimo la temperatura media global a partir de la temperatura
# media estimada para cada municipio, usando el método de la 
# ponderación por población
################################################################

# Multiplico la temperatura media de cada municipio por su población y sumo los valores de los 78 municipios (numerador de la ecuación)

temp.municipios.final <- temp.municipios.final[, -c(1, 80)]
temp.pob <- sweep(temp.municipios.final, 2, as.vector(as.data.frame(Asturias)[, 9]), FUN = '*')

temp.pob.sum <- rowSums(temp.pob)

# Sumo los valores de población de los 78 municipios (denominador de la ecuación)

pob <- as.vector(as.data.frame(Asturias)[, 9])
pob.suma <- sum(pob)
temp.pob.sum <- as.data.frame(temp.pob.sum)

# Calculo el numerador dividido por el denominador obteniendo la temperatura media global de Asturias

temp.global.final <- sweep(temp.pob.sum, 2, pob.suma, FUN = '/')

rm(list = ls()[!ls() %in% c(
  "Asturias",
  "datos.imputacion.final",
  "estaciones_sf",
  "poblacion",
  "temp.municipios.final",
  "radio",
  "temp.global.final"
)])

# Guardo la base de datos con la temperatura media global

temp.global.final$fecha <- datos.imputacion.final$fecha
temp.global.final <- temp.global.final %>%
  dplyr::select(fecha, everything())

#save(temp.global.final, file = "./data/TmedGlobal.RData")

################################################################
# CÁLCULO DE LA TEMPERATURA POR COMARCAS

# Estimo la temperatura media por comarcas a partir de la temperatura
# media estimada para cada municipio, usando el método de la 
# ponderación por población
################################################################

# Creo una base de datos para cada comarca, donde cada columna es la temperatura media de los municipios que forman parte de esa comarca

load("./data/TmedMunicipios.RData")

temp.aviles <- temp.municipios.final %>%
  dplyr::select(1, X4, X10, X16, X20, X21, X25, X30, X39, X51, X69)

temp.caudal <- temp.municipios.final %>%
  dplyr::select(1, X2, X33, X37)

temp.eonavia <- temp.municipios.final %>%
  dplyr::select(1, X7, X17, X18, X23, X27, X29, X41, X48, X61, X63, X62, X70, X71, X34, X74, X75, X77)

temp.gijon <- temp.municipios.final %>%
  dplyr::select(1, X14, X24, X76)

temp.nalon <- temp.municipios.final %>%
  dplyr::select(1, X15, X31, X32, X60, X67)

temp.narcea <- temp.municipios.final %>%
  dplyr::select(1, X1, X11, X22, X28, X73)

temp.oriente <- temp.municipios.final %>%
  dplyr::select(1, X3, X8, X12, X13, X19, X36, X43, X45, X46, X47, X49, X50, X55, X56)

temp.oviedo <- temp.municipios.final %>%
  dplyr::select(1, X5, X6, X9, X26, X35, X38, X40, X42, X44, X52, X53, X54, X57, X58, X59, X64, X65, X66, X68, X72, X78)

# Selecciono la información de los municipios que forman cada comarca

asturias.aviles <- Asturias %>%
  filter(cmun == "4"|cmun == "10"|cmun == "16"|cmun == "20"|cmun == "21"|cmun == "25"|cmun == "30"|cmun == "39"|cmun == "51"|cmun == "69")

asturias.caudal <- Asturias %>%
  filter(cmun == "2"|cmun == "33"|cmun == "37")

asturias.eonavia <- Asturias %>%
  filter(cmun == "7"|cmun == "17"|cmun == "18"|cmun == "23"|cmun == "27"|cmun == "29"|cmun == "41"|cmun == "48"|cmun == "61"|cmun == "63"|cmun == "62"|cmun == "70"|cmun == "71"|cmun == "34"|cmun == "74"|cmun == "75"|cmun == "77")

asturias.gijon <- Asturias %>%
  filter(cmun == "14"|cmun == "24"|cmun == "76")

asturias.nalon <- Asturias %>%
  filter(cmun == "15"|cmun == "31"|cmun == "32"|cmun == "60"|cmun == "67")

asturias.narcea <- Asturias %>%
  filter(cmun == "1"|cmun == "11"|cmun == "22"|cmun == "28"|cmun == "73")

asturias.oriente <- Asturias %>%
  filter(cmun == "3"|cmun == "8"|cmun == "12"|cmun == "13"|cmun == "19"|cmun == "36"|cmun == "43"|cmun == "45"|cmun == "46"|cmun == "47"|cmun == "49"|cmun == "50"|cmun == "55"|cmun == "56")

asturias.oviedo <- Asturias %>%
  filter(cmun == "5"|cmun == "6"|cmun == "9"|cmun == "26"|cmun == "35"|cmun == "38"|cmun == "40"|cmun == "42"|cmun == "44"|cmun == "52"|cmun == "53"|cmun == "54"|cmun == "57"|cmun == "58"|cmun == "59"|cmun == "64"|cmun == "65"|cmun == "66"|cmun == "68"|cmun == "72"|cmun == "78")

# Estimación de la temperatura media de la comarca de AVILÉS, usando el método de la ponderación por población usado anteriormente

temp.aviles <- temp.aviles[, -c(1)]
temp.aviles2 <- sweep(temp.aviles, 2, as.vector(as.data.frame(asturias.aviles)[, 9]), FUN = '*')

temp.aviles2 <- rowSums(temp.aviles2)

pob.aviles <- as.vector(as.data.frame(asturias.aviles)[, 9])
pob.aviles <- sum(pob.aviles)
temp.aviles2 <- as.data.frame(temp.aviles2)

temp.aviles.final <- sweep(temp.aviles2, 2, pob.aviles, FUN = '/')

# Estimación de la temperatura media de la comarca de CAUDAL, usando el método de la ponderación por población usado anteriormente

temp.caudal <- temp.caudal[, -c(1)]
temp.caudal2 <- sweep(temp.caudal, 2, as.vector(as.data.frame(asturias.caudal)[, 9]), FUN = '*')

temp.caudal2 <- rowSums(temp.caudal2)

pob.caudal <- as.vector(as.data.frame(asturias.caudal)[, 9])
pob.caudal <- sum(pob.caudal)
temp.caudal2 <- as.data.frame(temp.caudal2)

temp.caudal.final <- sweep(temp.caudal2, 2, pob.caudal, FUN = '/')

# Estimación de la temperatura media de la comarca de EO-NAVIA, usando el método de la ponderación por población usado anteriormente

temp.eonavia <- temp.eonavia[, -c(1)]
temp.eonavia <- temp.eonavia %>%
  dplyr::select(1:X29, X34, X41:X61, X62, X63, everything(.))
temp.eonavia2 <- sweep(temp.eonavia, 2, as.vector(as.data.frame(asturias.eonavia)[, 9]), FUN = '*')

temp.eonavia2 <- rowSums(temp.eonavia2)

pob.eonavia <- as.vector(as.data.frame(asturias.eonavia)[, 9])
pob.eonavia <- sum(pob.eonavia)
temp.eonavia2 <- as.data.frame(temp.eonavia2)

temp.eonavia.final <- sweep(temp.eonavia2, 2, pob.eonavia, FUN = '/')

# Estimación de la temperatura media de la comarca de GIJÓN, usando el método de la ponderación por población usado anteriormente

temp.gijon <- temp.gijon[, -c(1)]
temp.gijon2 <- sweep(temp.gijon, 2, as.vector(as.data.frame(asturias.gijon)[, 9]), FUN = '*')

temp.gijon2 <- rowSums(temp.gijon2)

pob.gijon <- as.vector(as.data.frame(asturias.gijon)[, 9])
pob.gijon <- sum(pob.gijon)
temp.gijon2 <- as.data.frame(temp.gijon2)

temp.gijon.final <- sweep(temp.gijon2, 2, pob.gijon, FUN = '/')

# Estimación de la temperatura media de la comarca de NALÓN, usando el método de la ponderación por población usado anteriormente

temp.nalon <- temp.nalon[, -c(1)]
temp.nalon2 <- sweep(temp.nalon, 2, as.vector(as.data.frame(asturias.nalon)[, 9]), FUN = '*')

temp.nalon2 <- rowSums(temp.nalon2)

pob.nalon <- as.vector(as.data.frame(asturias.nalon)[, 9])
pob.nalon <- sum(pob.nalon)
temp.nalon2 <- as.data.frame(temp.nalon2)

temp.nalon.final <- sweep(temp.nalon2, 2, pob.nalon, FUN = '/')

# Estimación de la temperatura media de la comarca de NARCEA, usando el método de la ponderación por población usado anteriormente

temp.narcea <- temp.narcea[, -c(1)]
temp.narcea2 <- sweep(temp.narcea, 2, as.vector(as.data.frame(asturias.narcea)[, 9]), FUN = '*')

temp.narcea2 <- rowSums(temp.narcea2)

pob.narcea <- as.vector(as.data.frame(asturias.narcea)[, 9])
pob.narcea <- sum(pob.narcea)
temp.narcea2 <- as.data.frame(temp.narcea2)

temp.narcea.final <- sweep(temp.narcea2, 2, pob.narcea, FUN = '/')

# Estimación de la temperatura media de la comarca de ORIENTE, usando el método de la ponderación por población usado anteriormente

temp.oriente <- temp.oriente[, -c(1)]
temp.oriente2 <- sweep(temp.oriente, 2, as.vector(as.data.frame(asturias.oriente)[, 9]), FUN = '*')

temp.oriente2 <- rowSums(temp.oriente2)

pob.oriente <- as.vector(as.data.frame(asturias.oriente)[, 9])
pob.oriente <- sum(pob.oriente)
temp.oriente2 <- as.data.frame(temp.oriente2)

temp.oriente.final <- sweep(temp.oriente2, 2, pob.oriente, FUN = '/')

# Estimación de la temperatura media de la comarca de OVIEDO, usando el método de la ponderación por población usado anteriormente

temp.oviedo <- temp.oviedo[, -c(1)]
temp.oviedo2 <- sweep(temp.oviedo, 2, as.vector(as.data.frame(asturias.oviedo)[, 9]), FUN = '*')

temp.oviedo2 <- rowSums(temp.oviedo2)

pob.oviedo <- as.vector(as.data.frame(asturias.oviedo)[, 9])
pob.oviedo <- sum(pob.oviedo)
temp.oviedo2 <- as.data.frame(temp.oviedo2)

temp.oviedo.final <- sweep(temp.oviedo2, 2, pob.oviedo, FUN = '/')

fecha <- temp.global.final$fecha

# Construyo la base final con la temperatura media de cada comarca

temp.comarcas.final <- cbind(
  fecha,
  temp.aviles.final,
  temp.caudal.final,
  temp.eonavia.final,
  temp.gijon.final,
  temp.nalon.final,
  temp.narcea.final,
  temp.oriente.final,
  temp.oviedo.final
)

#save(temp.comarcas.final, file = "./data/TmedComarcas.RData")

