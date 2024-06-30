# Nota: este script no va a poder reproducirse debido a la privacidad de la base de datos original, protegida mediante contraseña. No puedo dejar visibles los datos originales sin antes pedir permiso a quién me los compartió, sin embargo el archivo final con los datos agregados está disponible para poder reproducir el análisis.

packages = c('readxl', 'mapSpain', 'tidyverse', 'lubridate', 'ggplot2',  'gridExtra', 'excel.link', "stringr", "forecast")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos='http://cran.rediris.es')
    library(x, character.only = TRUE)
  }
})

################################################################
# SERIE GLOBAL

# Creo la serie temporal diaria del número de suicidios y
# uno los datos con los de temperatura
################################################################

# Importación de datos de suicidios

datos <- xl.read.file("./data/Suicidios Asturias desde 87-valen.xlsx", password = "MiPassword")

# Convierto las fechas en formato fecha

datos$FECH_NAC <- as.Date(datos$FECH_NAC)
datos$FECH_DEF <- as.Date(datos$FECH_DEF)

# Calculo la edad al momento del suicidio en años cumplidos a partir de la fecha de nacimiento y de la de defunción

datos <- datos %>%
  mutate(edad = (FECH_NAC %--% FECH_DEF) / dyears(1))
datos$edad <- floor(datos$edad)

# Le pongo etiquetas al sexo

datos <- datos %>%
  mutate(SEXO = ifelse(SEXO == 1, "Hombres", "Mujeres"))

# Creo la variable de grupo de edad con 3 grupos de edad: < 16 años (jovenes), entre 16 y 65 años (edad laboral) y mayores de 65 años (edad de jubilación).

datos <- datos %>%
  mutate(grupo_edad = ifelse(
    edad >= 0 &
      edad <= 15,
    "< 16",
    ifelse(edad >= 16 & edad <= 65, "16-65", "> 65")
  ))

## Comprobación municipios de residencia y defunción: los que se suicidaron fuera de Asturias se excluirán 
## debido a que es menos probable que hayan estado expuestos a la temperatura de la comunidad

# Importo la cartografía de Asturias a partir del paquete mapSpain

Asturias <- esp_get_munic_siane(region = "Asturias") %>%
  mutate(Provincia = esp_dict_translate(ine.prov.name, "es"))

# Estudio los códigos del municipio de defunción: si algún código no corresponde a ningún municipio de Asturias quiere decir que el suicidio ocurrió en otra comunidad autónoma. Extraigo los códigos con 2 dígitos (así como aparecen en la base de datos)

Asturias2 <- Asturias %>%
  mutate(cmun = str_remove(cmun, "^0+"))

# Estudio del municipio de residencia

table(datos$MUNIC_RE %in% c(Asturias2$cmun)) # Todos los individuos residían en Asturias

# Caso1: municipio de residencia y defunción coinciden

table(datos$MUNIC_RE == datos$MUNIC_DE) # 3955 de 4810 (82.2%)

# Caso2: el municipio de residencia y el de defunción no coinciden pero ambos se encuentran en Asturias

solo.asturias <- datos %>%
  filter(MUNIC_RE %in% c(Asturias2$cmun) &
           MUNIC_DE %in% c(Asturias2$cmun))

table(solo.asturias$MUNIC_RE != solo.asturias$MUNIC_DE) # 674 de 4810 (14%)

# Caso3: el municipio de residencia y el de defunción no coinciden y el de defunción no pertenece a Asturias

no.asturias <- datos %>%
  filter(!datos$MUNIC_DE %in% c(Asturias2$cmun))

table(is.na(no.asturias$MUNIC_DE)) # 20 no se suicidaron en Asturias (0.4%) y de 161 (3.4%) no sabemos donde se suicidaron (corresponden a todos los suicidios de 1991).

# Creo una variable para identificar los 4 casos: coinciden, diferentes en Asturias, diferentes en otra CA, dato faltante

datos <- datos %>%
  mutate(municipios = ifelse(
    MUNIC_RE == MUNIC_DE,
    "coinciden",
    ifelse(
      MUNIC_RE != MUNIC_DE &
        MUNIC_DE %in% c(Asturias2$cmun),
      "diferentes Asturias",
      "diferentes otra"
    )
  ))

datos <- datos %>%
  mutate(municipios = ifelse(is.na(MUNIC_DE), "DEF faltante", municipios))

table(datos$municipios)

# Elimino los que no se murieron en Asturias ya que es muy probable que no hayan estado expuestos a la temperatura de Asturias.

datos <- datos %>%
  filter(municipios != "diferentes otra")

# Asumo que los que tienen dato faltante en el municipio de defunción se suicidaron en el mismo municipio de residencia

datos <- datos %>%
  mutate(MUNIC_DE = ifelse(is.na(MUNIC_DE), MUNIC_RE, MUNIC_DE))

# Creo la variable comarca de residencia agrupando los municipios

datos <- datos %>%
  dplyr::mutate(comarca_RE = ifelse(MUNIC_RE == "4"|MUNIC_RE == "10"|MUNIC_RE == "16"|MUNIC_RE == "20"|MUNIC_RE == "21"|MUNIC_RE == "25"|MUNIC_RE == "30"|MUNIC_RE == "39"|MUNIC_RE == "51"|MUNIC_RE == "69", "aviles", ifelse(MUNIC_RE == "2"|MUNIC_RE == "33"|MUNIC_RE == "37", "caudal", ifelse(MUNIC_RE == "7"|MUNIC_RE == "17"|MUNIC_RE == "18"|MUNIC_RE == "23"|MUNIC_RE == "27"|MUNIC_RE == "29"|MUNIC_RE == "41"|MUNIC_RE == "48"|MUNIC_RE == "61"|MUNIC_RE == "63"|MUNIC_RE == "62"|MUNIC_RE == "70"|MUNIC_RE == "71"|MUNIC_RE == "34"|MUNIC_RE == "74"|MUNIC_RE == "75"|MUNIC_RE == "77", "eonavia", ifelse(MUNIC_RE == "14"|MUNIC_RE == "24"|MUNIC_RE == "76", "gijon", ifelse(MUNIC_RE == "15"|MUNIC_RE == "31"|MUNIC_RE == "32"|MUNIC_RE == "60"|MUNIC_RE == "67", "nalon", ifelse(MUNIC_RE == "1"|MUNIC_RE == "11"|MUNIC_RE == "22"|MUNIC_RE == "28"|MUNIC_RE == "73", "narcea", ifelse(MUNIC_RE == "3"|MUNIC_RE == "8"|MUNIC_RE == "12"|MUNIC_RE == "13"|MUNIC_RE == "19"|MUNIC_RE == "36"|MUNIC_RE == "43"|MUNIC_RE == "45"|MUNIC_RE == "46"|MUNIC_RE == "47"|MUNIC_RE == "49"|MUNIC_RE == "50"|MUNIC_RE == "55"|MUNIC_RE == "56", "oriente", ifelse(MUNIC_RE == "5"|MUNIC_RE == "6"|MUNIC_RE == "9"|MUNIC_RE == "26"|MUNIC_RE == "35"|MUNIC_RE == "38"|MUNIC_RE == "40"|MUNIC_RE == "42"|MUNIC_RE == "44"|MUNIC_RE == "52"|MUNIC_RE == "53"|MUNIC_RE == "54"|MUNIC_RE == "57"|MUNIC_RE == "58"|MUNIC_RE == "59"|MUNIC_RE == "64"|MUNIC_RE == "65"|MUNIC_RE == "66"|MUNIC_RE == "68"|MUNIC_RE == "72"|MUNIC_RE == "78", "oviedo", NA)))))))))

# Creo la variable comarca de defunción agrupando los municipios

datos <- datos %>%
  dplyr::mutate(comarca_DE = ifelse(MUNIC_DE == "4"|MUNIC_DE == "10"|MUNIC_DE == "16"|MUNIC_DE == "20"|MUNIC_DE == "21"|MUNIC_DE == "25"|MUNIC_DE == "30"|MUNIC_DE == "39"|MUNIC_DE == "51"|MUNIC_DE == "69", "aviles", ifelse(MUNIC_DE == "2"|MUNIC_DE == "33"|MUNIC_DE == "37", "caudal", ifelse(MUNIC_DE == "7"|MUNIC_DE == "17"|MUNIC_DE == "18"|MUNIC_DE == "23"|MUNIC_DE == "27"|MUNIC_DE == "29"|MUNIC_DE == "41"|MUNIC_DE == "48"|MUNIC_DE == "61"|MUNIC_DE == "63"|MUNIC_DE == "62"|MUNIC_DE == "70"|MUNIC_DE == "71"|MUNIC_DE == "34"|MUNIC_DE == "74"|MUNIC_DE == "75"|MUNIC_DE == "77", "eonavia", ifelse(MUNIC_DE == "14"|MUNIC_DE == "24"|MUNIC_DE == "76", "gijon", ifelse(MUNIC_DE == "15"|MUNIC_DE == "31"|MUNIC_DE == "32"|MUNIC_DE == "60"|MUNIC_DE == "67", "nalon", ifelse(MUNIC_DE == "1"|MUNIC_DE == "11"|MUNIC_DE == "22"|MUNIC_DE == "28"|MUNIC_DE == "73", "narcea", ifelse(MUNIC_DE == "3"|MUNIC_DE == "8"|MUNIC_DE == "12"|MUNIC_DE == "13"|MUNIC_DE == "19"|MUNIC_DE == "36"|MUNIC_DE == "43"|MUNIC_DE == "45"|MUNIC_DE == "46"|MUNIC_DE == "47"|MUNIC_DE == "49"|MUNIC_DE == "50"|MUNIC_DE == "55"|MUNIC_DE == "56", "oriente", ifelse(MUNIC_DE == "5"|MUNIC_DE == "6"|MUNIC_DE == "9"|MUNIC_DE == "26"|MUNIC_DE == "35"|MUNIC_DE == "38"|MUNIC_DE == "40"|MUNIC_DE == "42"|MUNIC_DE == "44"|MUNIC_DE == "52"|MUNIC_DE == "53"|MUNIC_DE == "54"|MUNIC_DE == "57"|MUNIC_DE == "58"|MUNIC_DE == "59"|MUNIC_DE == "64"|MUNIC_DE == "65"|MUNIC_DE == "66"|MUNIC_DE == "68"|MUNIC_DE == "72"|MUNIC_DE == "78", "oviedo", NA)))))))))

## Creación de la serie temporal diaria

# Calculo el número total de suicidios según la fecha de defunción

datos.serie <- datos %>%
  group_by(FECH_DEF) %>%
  summarise(Total = n()) %>%
  ungroup() 

# Añado a la serie los días que van desde el 1 de enero de 1987 hasta el 31 de diciembre de 2021 en los que no se observó ningún suicidio, para tener la serie temporal completa

ts <- seq.POSIXt(as.POSIXct("1987-01-02"), as.POSIXct("2022-01-01"), by = "day")
ts <- data.frame(FECH_DEF = ts)
ts$FECH_DEF <- as.Date(ts$FECH_DEF)
datos.serie <- ts %>%
  left_join(datos.serie)

# A las fechas en las que no se observó ningún suicidio les asigno un valor de 0 suicidios totales

datos.serie <- datos.serie %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

datos.serie <- datos.serie %>%
  rename(fecha = FECH_DEF)

## Importación de datos de temperatura

load("./data/TmedGlobal.RData")

# Uno los datos de suicidio con los de temperatura

datos.serie <- datos.serie %>%
  left_join(temp.global.final)

# Creo las variables de año, mes, día y día de la semana

datos.serie$yy <- year(datos.serie$fecha)
datos.serie$mm <- month(datos.serie$fecha)
datos.serie$dd <- day(datos.serie$fecha)
datos.serie$dow <- weekdays(datos.serie$fecha)
datos.serie$dow <- as.factor(datos.serie$dow)
datos.serie$dow <- relevel(datos.serie$dow, ref = "jueves") # Establezco el jueves como día de referencia

################################################################
# SERIES SEGÚN EL MUNICIPIO DE DEFUNCIÓN

# Creo la serie temporal diaria del número de suicidios según el 
# municipio de defunción y uno los datos con los de temperatura
################################################################

# En este caso, elimino los que no se murieron en el mismo municipio de residencia ya que es muy probable que no hayan estado expuestos a la temperatura de ese municipio.

datos.mun <- datos %>%
  filter(municipios == "coinciden" | municipios == "DEF faltante")

# Quito la variable de municipio de defunción, me quedo sólo con el de residencia ya que los dos coinciden y, en los casos de 1991, asumo que el municipio de defunción es igual al de residencia

datos.mun <- datos.mun %>%
  dplyr::select(-MUNIC_DE) %>%
  rename(MUNIC_DE = MUNIC_RE)

# Calculo el número total de suicidios según la fecha de defunción y el municipio

datos.serie.mun <- datos.mun %>%
  group_by(FECH_DEF, MUNIC_DE) %>%
  summarise(Total = n()) %>%
  ungroup()

# Creo la variable fecha que va del 1 de enero del 1987 al 31 de diciembre de 2021 para cada uno de los 78 municipios

ts.mun <- rep(ts$FECH_DEF, 78)

# Creo un vector con los códigos de los municipios y lo ordeno

mun <- c(rep(1:78, 12784))
mun <- sort(mun)

# Creo un dataframe con la fecha y el identificador del municipio

fecha.mun <- data.frame(FECH_DEF = ts.mun, MUNIC_DE = mun)

# Obtengo una serie por cada municipio

datos.serie.mun <- fecha.mun %>%
  left_join(datos.serie.mun)

# A las fechas en las que no se observó ningún suicidio les asigno el valor de 0 suicidios

datos.serie.mun <- datos.serie.mun %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

# Importo los datos de temperatura media por municipios

load("./data/TmedMunicipios.RData")

# Cambio el nombre de las columnas añadiendo el código del municipio

names(temp.municipios.final)[-1] <- c(paste("tmed", 1:78, sep = "."))

# Cambio el formato de la base de datos de wide a long

long.tmun <- reshape(temp.municipios.final,
                     varying = c(2:79),
                     direction = "long")

long.tmun <- long.tmun %>%
  rename(FECH_DEF = fecha, MUNIC_DE = time)

long.tmun <- long.tmun %>%
  dplyr::select(-id)

# Uno los datos de suicidio con los de temperatura

datos.serie.mun <- long.tmun %>%
  left_join(datos.serie.mun)

# Creo las variables de año, mes, día y día de la semana

datos.serie.mun$FECH_DEF <- as.Date(datos.serie.mun$FECH_DEF)

datos.serie.mun$MUNIC_DE <- as.factor(datos.serie.mun$MUNIC_DE)

datos.serie.mun$yy <- year(datos.serie.mun$FECH_DEF)
datos.serie.mun$mm <- month(datos.serie.mun$FECH_DEF)
datos.serie.mun$dd <- day(datos.serie.mun$FECH_DEF)
datos.serie.mun$dow <- weekdays(datos.serie.mun$FECH_DEF)
datos.serie.mun$dow <- as.factor(datos.serie.mun$dow)
datos.serie.mun$dow <- relevel(datos.serie.mun$dow, ref = "jueves") # Establezco el jueves como día de referencia

datos.serie.mun <- datos.serie.mun %>%
  rename(fecha = FECH_DEF)

################################################################
# SERIES SEGÚN LA COMARCA DE DEFUNCIÓN

# Creo la serie temporal diaria del número de suicidios según la 
# comarca de defunción y uno los datos con los de temperatura
################################################################

# Excluyo las personas que se suicidaron en una comarca diferente a la de residencia ya que es menos probable que hayan estado expuestas a la temperatura de esa comarca

datos <- datos %>%
  filter(comarca_RE == comarca_DE)

# Calculo el número total de suicidios según la fecha de defunción y la comarca

datos.serie.com <- datos %>%
  group_by(FECH_DEF, comarca_DE) %>%
  summarise(Total = n()) %>%
  ungroup() 

# Creo la variable fecha que va del 1 de enero del 1987 al 31 de diciembre de 2021 para cada una de las 8 comarcas

ts.com <- rep(ts$FECH_DEF, 8)

# Creo un vector con los nombres de las comarcas y lo ordeno

com <- c(rep(
  c("aviles",
    "caudal",
    "eonavia",
    "gijon",
    "nalon",
    "narcea",
    "oriente",
    "oviedo"),
  12784
))

com <- sort(com)

# Creo un dataframe con la fecha y el nombre de la comarca

fecha.com <- data.frame(FECH_DEF = ts.com, comarca_DE = com)

# Obtengo una serie por cada comarca

datos.serie.com <- fecha.com %>%
  left_join(datos.serie.com)

# A las fechas en las que no se observó ningún suicidio les asigno el valor de 0 suicidios

datos.serie.com <- datos.serie.com %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

# Importo los datos de temperatura media por comarcas

load("./data/TmedComarcas.RData")

# Cambio el nombre de las columnas y añado el nombre de la comarca

names(temp.comarcas.final)[-1] <- c(paste(
  "tmed",
  c("aviles",
    "caudal",
    "eonavia",
    "gijon",
    "nalon",
    "narcea",
    "oriente",
    "oviedo"),
  sep = "."
))

# Cambio el formato de la base de datos de wide a long

long.tcom <- reshape(temp.comarcas.final,
                     varying = c(2:9),
                     direction = "long")

long.tcom <- long.tcom %>%
  rename(FECH_DEF = fecha, comarca_DE = time)

long.tcom <- long.tcom %>%
  dplyr::select(-id)

# Uno los datos de suicidio con los de temperatura

datos.serie.com <- long.tcom %>%
  left_join(datos.serie.com)

# Creo las variables de año, mes, día y día de la semana

datos.serie.com$FECH_DEF <- as.Date(datos.serie.com$FECH_DEF)

datos.serie.com$comarca_DE <- as.factor(datos.serie.com$comarca_DE)

datos.serie.com$yy <- year(datos.serie.com$FECH_DEF)
datos.serie.com$mm <- month(datos.serie.com$FECH_DEF)
datos.serie.com$dd <- day(datos.serie.com$FECH_DEF)
datos.serie.com$dow <- weekdays(datos.serie.com$FECH_DEF)
datos.serie.com$dow <- as.factor(datos.serie.com$dow)
datos.serie.com$dow <- relevel(datos.serie.com$dow, ref = "jueves") # Establezco el jueves como día de referencia

datos.serie.com <- datos.serie.com %>%
  rename(fecha = FECH_DEF)

################################################################
# SERIES SEGÚN EL SEXO

# Creo la serie temporal diaria del número de suicidios según el 
# sexo y la comarca y uno los datos con los de temperatura
################################################################

## HOMBRES

# Selecciono sólo los suicidios cometidos por hombres

datos.h <- datos %>%
  filter(SEXO == "Hombres")

# Calculo el número total de suicidios según la fecha de defunción y la comarca

datos.serie.h <- datos.h %>%
  group_by(FECH_DEF, comarca_DE) %>%
  summarise(Total = n()) %>%
  ungroup() 

# Obtengo una serie por cada comarca

datos.serie.h <- fecha.com %>%
  left_join(datos.serie.h)

# A las fechas en las que no se observó ningún suicidio les asigno el valor de 0 suicidios

datos.serie.h <- datos.serie.h %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

# Uno los datos de suicidio con los de temperatura

datos.serie.h <- long.tcom %>%
  left_join(datos.serie.h)

# Creo las variables de año, mes, día y día de la semana

datos.serie.h$FECH_DEF <- as.Date(datos.serie.h$FECH_DEF)

datos.serie.h$comarca_DE <- as.factor(datos.serie.h$comarca_DE)

datos.serie.h$yy <- year(datos.serie.h$FECH_DEF)
datos.serie.h$mm <- month(datos.serie.h$FECH_DEF)
datos.serie.h$dd <- day(datos.serie.h$FECH_DEF)
datos.serie.h$dow <- weekdays(datos.serie.h$FECH_DEF)
datos.serie.h$dow <- as.factor(datos.serie.h$dow)
datos.serie.h$dow <- relevel(datos.serie.h$dow, ref = "jueves") # Establezco el jueves como día de referencia

datos.serie.h <- datos.serie.h %>%
  rename(fecha = FECH_DEF)

## MUJERES

# Selecciono sólo los suicidios cometidos por mujeres

datos.m <- datos %>%
  filter(SEXO == "Mujeres")

# Calculo el número total de suicidios según la fecha de defunción y la comarca

datos.serie.m <- datos.m %>%
  group_by(FECH_DEF, comarca_DE) %>%
  summarise(Total = n()) %>%
  ungroup()

# Obtengo una serie por cada comarca

datos.serie.m <- fecha.com %>%
  left_join(datos.serie.m)

# A las fechas en las que no se observó ningún suicidio les asigno el valor de 0 suicidios

datos.serie.m <- datos.serie.m %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

# Uno los datos de suicidio con los de temperatura

datos.serie.m <- long.tcom %>%
  left_join(datos.serie.m)

# Creo las variables de año, mes, día y día de la semana

datos.serie.m$FECH_DEF <- as.Date(datos.serie.m$FECH_DEF)

datos.serie.m$comarca_DE <- as.factor(datos.serie.m$comarca_DE)

datos.serie.m$yy <- year(datos.serie.m$FECH_DEF)
datos.serie.m$mm <- month(datos.serie.m$FECH_DEF)
datos.serie.m$dd <- day(datos.serie.m$FECH_DEF)
datos.serie.m$dow <- weekdays(datos.serie.m$FECH_DEF)
datos.serie.m$dow <- as.factor(datos.serie.m$dow)
datos.serie.m$dow <- relevel(datos.serie.m$dow, ref = "jueves") # Establezco el jueves como día de referencia

datos.serie.m <- datos.serie.m %>%
  rename(fecha = FECH_DEF)

################################################################
# SERIES SEGÚN EL GRUPO DE EDAD

# Creo la serie temporal diaria del número de suicidios según el 
# grupo de edad y la comarca y uno los datos con los de temperatura
################################################################

## < 16 AÑOS

# Selecciono sólo los suicidios cometidos por gente con menos de 16 años

datos.015 <- datos %>%
  filter(grupo_edad == "< 16")

# Calculo el número total de suicidios según la fecha de defunción y la comarca

datos.serie.015 <- datos.015 %>%
  group_by(FECH_DEF, comarca_DE) %>%
  summarise(Total = n()) %>%
  ungroup() 

# Obtengo una serie por cada comarca

datos.serie.015 <- fecha.com %>%
  left_join(datos.serie.015)

# A las fechas en las que no se observó ningún suicidio les asigno el valor de 0 suicidios

datos.serie.015 <- datos.serie.015 %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

# Uno los datos de suicidio con los de temperatura

datos.serie.015 <- long.tcom %>%
  left_join(datos.serie.015)

# Creo las variables de año, mes, día y día de la semana

datos.serie.015$FECH_DEF <- as.Date(datos.serie.015$FECH_DEF)

datos.serie.015$comarca_DE <- as.factor(datos.serie.015$comarca_DE)

datos.serie.015$yy <- year(datos.serie.015$FECH_DEF)
datos.serie.015$mm <- month(datos.serie.015$FECH_DEF)
datos.serie.015$dd <- day(datos.serie.015$FECH_DEF)
datos.serie.015$dow <- weekdays(datos.serie.015$FECH_DEF)
datos.serie.015$dow <- as.factor(datos.serie.015$dow)
datos.serie.015$dow <- relevel(datos.serie.015$dow, ref = "jueves") # Establezco el jueves como día de referencia

datos.serie.015 <- datos.serie.015 %>%
  rename(fecha = FECH_DEF)

## 16-65 AÑOS

# Selecciono sólo los suicidios cometidos por gente de 16-65 años

datos.1665 <- datos %>%
  filter(grupo_edad == "16-65")

# Calculo el número total de suicidios según la fecha de defunción y la comarca

datos.serie.1665 <- datos.1665 %>%
  group_by(FECH_DEF, comarca_DE) %>%
  summarise(Total = n()) %>%
  ungroup() 

# Obtengo una serie por cada comarca

datos.serie.1665 <- fecha.com %>%
  left_join(datos.serie.1665)

# A las fechas en las que no se observó ningún suicidio les asigno el valor de 0 suicidios

datos.serie.1665 <- datos.serie.1665 %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

# Uno los datos de suicidio con los de temperatura

datos.serie.1665 <- long.tcom %>%
  left_join(datos.serie.1665)

# Creo las variables de año, mes, día y día de la semana

datos.serie.1665$FECH_DEF <- as.Date(datos.serie.1665$FECH_DEF)

datos.serie.1665$comarca_DE <- as.factor(datos.serie.1665$comarca_DE)

datos.serie.1665$yy <- year(datos.serie.1665$FECH_DEF)
datos.serie.1665$mm <- month(datos.serie.1665$FECH_DEF)
datos.serie.1665$dd <- day(datos.serie.1665$FECH_DEF)
datos.serie.1665$dow <- weekdays(datos.serie.1665$FECH_DEF)
datos.serie.1665$dow <- as.factor(datos.serie.1665$dow)
datos.serie.1665$dow <- relevel(datos.serie.1665$dow, ref = "jueves") # Establezco el jueves como día de referencia

datos.serie.1665 <- datos.serie.1665 %>%
  rename(fecha = FECH_DEF)

## > 65 AÑOS

# Selecciono sólo los suicidios cometidos por gente mayor de 65 años

datos.65mas <- datos %>%
  filter(grupo_edad == "> 65")

# Calculo el número total de suicidios según la fecha de defunción y la comarca

datos.serie.65mas <- datos.65mas %>%
  group_by(FECH_DEF, comarca_DE) %>%
  summarise(Total = n()) %>%
  ungroup()

# Obtengo una serie por cada comarca

datos.serie.65mas <- fecha.com %>%
  left_join(datos.serie.65mas)

# A las fechas en las que no se observó ningún suicidio les asigno el valor de 0 suicidios

datos.serie.65mas <- datos.serie.65mas %>%
  mutate(Total = ifelse(is.na(Total), 0, Total))

# Uno los datos de suicidio con los de temperatura

datos.serie.65mas <- long.tcom %>%
  left_join(datos.serie.65mas)

# Creo las variables de año, mes, día y día de la semana

datos.serie.65mas$FECH_DEF <- as.Date(datos.serie.65mas$FECH_DEF)

datos.serie.65mas$comarca_DE <- as.factor(datos.serie.65mas$comarca_DE)

datos.serie.65mas$yy <- year(datos.serie.65mas$FECH_DEF)
datos.serie.65mas$mm <- month(datos.serie.65mas$FECH_DEF)
datos.serie.65mas$dd <- day(datos.serie.65mas$FECH_DEF)
datos.serie.65mas$dow <- weekdays(datos.serie.65mas$FECH_DEF)
datos.serie.65mas$dow <- as.factor(datos.serie.65mas$dow)
datos.serie.65mas$dow <- relevel(datos.serie.65mas$dow, ref = "jueves") # Establezco el jueves como día de referencia

datos.serie.65mas <- datos.serie.65mas %>%
  rename(fecha = FECH_DEF)

# Guardo todas las series obtenidas

# save(datos.serie, file = "./data/serie.global.RData")
# save(datos.serie.h, datos.serie.m, datos.serie.015, datos.serie.1665, datos.serie.65mas, file="./data/serie.grupos.RData")
# save(datos.serie.mun, file = "./data/serie.municipios.RData")
# save(datos.serie.com, file = "./data/serie.comarcas.RData")