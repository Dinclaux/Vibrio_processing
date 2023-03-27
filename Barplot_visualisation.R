# Reinitialize the session
rm(list=ls(all=TRUE))

#########################################################
###       Installing and loading required packages    ###
#########################################################
library(tidyverse)
library(readxl)
library(tidyr)
library(openxlsx)

#########################################################
###              working directory                    ###
#########################################################
# Changer le répertoire de travail
setwd("~/Script/TWB/Vibrio_processing/Sans conta blanc/")

# Lire les données brutes
raw_data <- read_excel("2023-03-14_export_total_OD_0.25.xlsx", sheet = "export") %>%
  drop_na(ratio_meanScarletOD)

# Suppression des colones inutiles
raw_data <- raw_data[,-1]

# Utiliser la fonction separate() pour séparer la colonne en trois colonnes distinctes
raw_data <- separate(raw_data, "Souches", into = c("Promoteur", "RBS", "Terminateur"), sep = "_",fill = "right",remove=FALSE)

# Modifs des colones en facteur
raw_data$Promoteur <- factor(raw_data$Promoteur)
raw_data$RBS <- factor(raw_data$RBS)
raw_data$Terminateur <- factor(raw_data$Terminateur)

donnees <- data.frame(Souches = raw_data$Souches, Valeurs = raw_data$ratio_meanScarletOD)

# Calculer les moyennes et les écarts-types
donnees_summary <- donnees %>%
  group_by(Souches) %>%
  summarize(mean_y = mean(Valeurs), sd_y = sd(Valeurs), n = n()) %>%
  arrange(Souches) # trier les données

# Mise à zéro des données négatives
donnees_summary$mean_y[donnees_summary$mean_y < 0] <- 0

# Définir une plage de valeurs à afficher
plage_min <- 100
plage_max <- 5000

# Filtrer les données selon la plage minimale et maximale
donnees_filtrees <- donnees_summary %>%
  filter(mean_y >= plage_min & mean_y <= plage_max)

# Identifier les données avec un écart-type supérieur à X% de la moyenne
donnees_filtrees$high_sd <- donnees_filtrees$sd_y > 0.2 * donnees_filtrees$mean_y

# Créer le graphique
ggplot(data = donnees_filtrees, aes(x = reorder(Souches, mean_y), y = mean_y, fill = high_sd)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = c("blue", "red"), guide = guide_legend(title = "Construits"),labels = c("Valide", "Sd > 20%", "Simplica")) +
  geom_errorbar(aes(ymin = mean_y - sd_y, ymax = mean_y + sd_y),
                width = 0.4) +
  labs(x = "Construits", y = "Ratio Scarlet/OD") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


write.xlsx(donnees_summary, "donnees_summary.xlsx")