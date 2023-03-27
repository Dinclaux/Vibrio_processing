# Reinitialize the session
rm(list=ls(all=TRUE))

#########################################################
###       Installing and loading required packages    ###
#########################################################
library(tidyverse)
library(readxl)
library(plotly)

#########################################################
###              working directory                    ###
#########################################################
# Changer le répertoire de travail
setwd("~/Script/TWB/Vibrio_processing")

# Lire les données brutes
raw_data <- read_excel("2022-09-30_Combinatoire VIBI_CONSTITUTIF.xlsx", sheet = "Bilan Constitutif") %>%
  drop_na(ratio_meanScarletOD)

# Suppression des colones inutiles
raw_data <- raw_data[,-1:-2]

# Modifs des colones en facteur
raw_data$Promoteur <- factor(raw_data$Promoteur)
raw_data$RBS <- factor(raw_data$RBS)
raw_data$Terminateur <- factor(raw_data$Terminateur)

donnees <- data.frame(Souches = raw_data$Souche, Valeurs = raw_data$ratio_meanScarletOD)

# Calculer les moyennes et les écarts-types
donnees_summary <- donnees %>%
  group_by(Souches) %>%
  summarize(mean_y = mean(Valeurs), sd_y = sd(Valeurs)) %>%
  arrange(Souches) # trier les données

# Mise à zéro des données négatives
donnees_summary$mean_y[donnees_summary$mean_y < 0] <- 0

# Définir une plage de valeurs à afficher
plage_min <- 0
plage_max <- 10000

# Filtrer les données selon la plage minimale et maximale
donnees_filtrees <- donnees_summary %>%
  filter(mean_y >= plage_min & mean_y <= plage_max)

# Identifier les données avec un écart-type supérieur à X% de la moyenne
donnees_filtrees$high_sd <- donnees_filtrees$sd_y > 0.2 * donnees_filtrees$mean_y

# Créer le graphique avec plotly
plot_ly(data = donnees_filtrees, x = ~reorder(Souches, mean_y), y = ~mean_y, type = 'bar', 
        error_y = list(type = "data", array = ~sd_y,color = "grey"), 
        marker = list(color = ~ifelse(high_sd, "red", "blue"))) %>%
  layout(xaxis = list(title = "Construits", tickangle = 90, showticklabels = TRUE), 
         yaxis = list(title = "Ratio Scarlet/OD"), 
         title = "Moyennes et écart-types des ratios Scarlet/OD des construits",
         barmode = "group")

