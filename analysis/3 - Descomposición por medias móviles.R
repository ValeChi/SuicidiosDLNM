packages = c('readxl', 'mapSpain', 'tidyverse', 'lubridate', 'ggplot2',  'gridExtra', 'excel.link', "stringr", "forecast")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos='http://cran.rediris.es')
    library(x, character.only = TRUE)
  }
})

load("./data/serie.global.RData")

################################################################
# ESTUDIO COMPONENTES DE LA SERIE DE SUICIDIOS

# Estudio de tendencia, estacionalidad y residuo de la serie
# de suicidios mediante descomposición por medias móviles
################################################################

## Serie diaria de suicidios

span <- 270 / nrow(datos.serie)

ggplot(datos.serie, aes(fecha, Total)) + geom_line() +
  stat_smooth(
    method = "loess",
    span = span,
    n = 1000,
    col = "#F44336",
    alpha = 0.5
  ) +
  scale_y_continuous("Nº de suicidios") +
  scale_x_date("Día", date_breaks = "1 year", date_labels = "%Y") + labs(title = "Nº diario de suicidios") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

## Serie mensual de suicidios para poder estudiar mejor sus componentes

# Obtengo el número de suicidios mensual

datos <- datos.serie %>%
  group_by(mm, yy) %>%
  summarise(Total = sum(Total)) %>%
  ungroup()

datos <- datos %>%
  arrange(yy)

# Convierto la serie en un objeto de serie temporal con estacionalidad de orden 12

serie <- ts(datos$Total, start = 1987, frequency = 12)

# Represento la serie mensual de suicidios

autoplot(serie,
         xlab = "Mes",
         ylab = "Nº de suicidios",
         main = "Nº mensual de suicidios") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(breaks = seq(from = 1986, to =  2023, by = 1))

# Estudio de la tendencia

serie.desc <- decompose(serie, type = "additive")

autoplot(
  as.ts(serie.desc$trend),
  xlab = "",
  ylab = "Nº de suicidios",
  main = "(A) Tendencia de la serie de nº de suicidios"
) +
  scale_x_continuous(breaks = seq(from = 1987, to =  2022, by = 2)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Estudio de la estacionalidad

ggseasonplot(
  serie,
  year.labels = TRUE,
  xlab = "",
  ylab = "Nº de suicidios",
  main = "Gráfico de estacionalidad de la serie\nde nº de suicidios",
  season.labels = c(
    "Ene",
    "Feb",
    "Mar",
    "Abr",
    "May",
    "Jun",
    "Jul",
    "Ago",
    "Sep",
    "Oct",
    "Nov",
    "Dic"
  )
) +
  theme_bw()

# Estudio del residuo

error.serie <- remainder(serie.desc)

sderror.serie <- sd(error.serie, na.rm = TRUE)

autoplot(
  error.serie,
  xlab = "",
  ylab = "Error",
  main = "Residuo de la serie de nº de suicidios",
  colour = "black"
) +
  geom_hline(
    yintercept = c(3, 2, -2, -3) * sderror.serie,
    colour = c("red", "green", "green", "red"),
    lty = 2
  ) +
  scale_x_continuous(breaks = seq(from = 1987, to = 2022, by = 2)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Miro en correspondencia de qué meses y años hay valores anómalos de la serie

abs(error.serie) > 3 * sderror.serie
