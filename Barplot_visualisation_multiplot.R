# Reinitialize the session
rm(list=ls(all=TRUE))

#########################################################
###       Installing and loading required packages    ###
#########################################################
library(tidyverse)
library(readxl)
library(tidyr)
library(openxlsx)
library(gridExtra)

#########################################################
###              working directory                    ###
#########################################################
# Changer le répertoire de travail
setwd("~/Script/TWB/Vibrio_processing/Sans conta blanc/Proc_barplot/")

# Récupérer la liste des fichiers Excel dans le dossier
file_names <- list.files(pattern = "*.xlsx")

# Initialiser la liste de données
data_list <- list()

# Boucle pour lire chaque fichier et stocker les données dans la liste
for (i in seq_along(file_names)) {
  # Charger les données
  raw_data <- read_excel(file_names[i], sheet = "export") %>%
    drop_na(ratio_meanScarletOD)
  
  # Supprimer les colonnes inutiles
  raw_data <- raw_data[,-1]
 

  # Séparer la colonne en trois colonnes distinctes
  raw_data <- separate(raw_data, "Souches", into = c("Promoteur", "RBS", "Terminateur"), sep = "_",fill = "right",remove=FALSE)
  
  # Modifier les colonnes en facteur
  raw_data$Promoteur <- factor(raw_data$Promoteur)
  raw_data$RBS <- factor(raw_data$RBS)
  raw_data$Terminateur <- factor(raw_data$Terminateur)
 
  extracton <- data.frame(Souches = raw_data$Souches, Valeurs = raw_data$ratio_meanScarletOD)
  # Calculer les moyennes et les écarts-types
  donnees_summary <- extracton %>%
    group_by(Souches) %>%
    summarize(mean_y = mean(Valeurs), sd_y = sd(Valeurs), n = n()) %>%
    arrange(Souches)
  
  # Mise à zéro des données négatives
  donnees_summary$mean_y[donnees_summary$mean_y < 0] <- 0
  # Stocker les données dans la liste
  donnees_summary <- cbind(donnees_summary,file_names[i])
  data_list[[i]] <- donnees_summary
}

# Fusionner toutes les données en un seul jeu de données
all_data <- dplyr::bind_rows(data_list)

# Définir une plage de valeurs à afficher
plage_min <- 100
plage_max <- 5000

# Filtrer les données selon la plage minimale et maximale
all_data <- all_data %>%
  filter(mean_y >= plage_min & mean_y <= plage_max)

# Identifier les données avec un écart-type supérieur à X% de la moyenne
all_data$high_sd <- all_data$sd_y > 0.2 * all_data$mean_y
all_data <- all_data %>% rename(Experiences = "file_names[i]")
all_data$Experiences <- (substr(all_data$Experiences, start = 24+1, stop = nchar(all_data$Experiences) - 5))
all_data$Experiences <- as.factor(all_data$Experiences)
#all_data$Experiences <-  factor(all_data$Experiences, levels = c("OD_0.2","OD_0.55",
#                                                                 "OD_0.25","OD_0.6",
#                                                                 "OD_0.3","OD_0.65",
#                                                                 "OD_0.35","OD_0.7",
#                                                                 "OD_0.4","OD_0.75",
#                                                                 "OD_0.45","OD_0.8",
#                                                                 "OD_0.5"))

p<-ggplot(data = all_data, aes(x = Souches, y = mean_y, fill = high_sd)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 2500))+
  scale_fill_manual(values = c("blue", "red"), guide = guide_legend(title = "Construits"),labels = c("Valide", "Sd > 20%", "Simplica")) +
  geom_errorbar(aes(ymin = mean_y - sd_y, ymax = mean_y + sd_y),
                width = 0.4) +
  labs(x = "Construits", y = "Ratio Scarlet/OD") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5,vjust = 0.2))

p + facet_wrap(~ Experiences, ncol = 1, scales = "fixed")

