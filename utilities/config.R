library(reshape2)
library(tidyverse)
library(kableExtra)
library(gridExtra)
library(patchwork)

# PAQUETES BAYESIANOS
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)

# Temas para graficas en ggplot
sin_lineas <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
sin_leyenda <- theme(legend.position = "none")
sin_ejes <- theme(axis.ticks = element_blank(), 
        axis.text = element_blank())

# Tablas en formato Latex a partir de un tibble
fancy_table <- function(X, caption, format='html'){
    "
    Transformación de tibbles a formato latex
    Params:
    X (tibble) : Tabla de datos
    caption (string) : Caption de la tabla
    "
   X |> kbl(caption= caption, booktabs = T, align = "l", format = format) |>
        kable_styling(latex_options = c("stripped","hold_position"))
}

# capacidad predictiva. Tomada de las notas "Comparacion de modelos." del curso
calcula_metricas <- function(log_lik){
    r_eff <- relative_eff(exp(log_lik), cores = 2) 
    within(list(), {
      #loo  <- loo(log_lik, r_eff = r_eff)
      waic <- waic(log_lik, r_eff = r_eff)
    })
  }