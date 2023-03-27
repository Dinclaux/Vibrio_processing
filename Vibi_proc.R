#2022-09-28 mickael.dinclaux@inrae.fr

#v1.0 intial script 28/09/2022
#v1.1 add batch possibility and progressbar

#script for extract µ and Scarlet rate from bioteck's raw data

# Reinitialize the session
rm(list=ls(all=TRUE))


#########################################################
###       Installing and loading required packages    ###
#########################################################


if (!require("xlsx")){
  install.packages("xlsx", dependencies = TRUE)
  library(xlsx)
}
if (!require("dplyr")){
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("stringr")){
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}

if (!require("tidyverse")){
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}
if (!require("ggplot2")){
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
#########################################################
###              working directory                    ###
#########################################################


setwd("~/Script/TWB/Vibrio_processing/Sans conta blanc")

batch = c("220719-Data cinétique Constitutif 1 alpha.xlsx",
          "220720-Data cinétique Constitutif 1 gamma.xlsx",
          "220721-Data cinétique Constitutif 1 Delta.xlsx",
          "220721-Data cinétique Constitutif 2 beta.xlsx",
          "220726-Data cinétique Constitutif 3 reste Sauf ParaBAD et Ptac plaque alpha.xlsx",
          "220801-Data cinetique Constitutif 1 gamma bis.xlsx",
          "220802-Data cinetique Constitutif 1 beta bis.xlsx",
          "220823-Data cinetique Constitutif 2 alpha bis.xlsx",
          "230213_Constitutif1 alpha.xlsx",
          "230214_Constitutif1 gamma.xlsx",
          "230215_Constitutif1 delta.xlsx",
          "230216_Constitutif2 alpha.xlsx",
          "230216_Constitutif2 beta.xlsx",
          "230221_Constitutif1 beta bis.xlsx")
 
ODseq =c()
ODseq = seq(0.200,0.800, by = 0.05)

#########################################################
###              Do not modify                        ###
#########################################################

# Initializes the progress bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = (length(batch)*length(ODseq)+1), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
progressbar <- 0

for(c in 1:length(ODseq)){
  OD <- ODseq[c]
  if (exists("export_mixdata") == TRUE) {
    rm(export_mixdata)
  }
  
for(b in 1:length(batch)){
  progressbar<-  progressbar + 1
  
setTxtProgressBar(pb, progressbar)  # Sets the progress bar to the current state
filename = batch[b]


data_OD <- read.xlsx(filename, header = TRUE, sheetName = "Do600", check.names = FALSE)
data_scarlet <- read.xlsx(filename, header = TRUE, sheetName = "mScarlet", check.names = FALSE)

# Supprimer les lignes avec des valeurs manquantes
data_OD <- na.omit(data_OD)
data_scarlet <- na.omit(data_scarlet)

#remove temp
name_col <- colnames(data_OD)
data_OD <- data_OD[,-2]
data_scarlet <- data_scarlet[,-2]
name_col <- name_col[-2]

#convertir les colonnes en numerique si character
data_OD[, 2:ncol(data_OD)] <- sapply(data_OD[, 2:ncol(data_OD)], as.numeric)
data_scarlet[, 2:ncol(data_scarlet)] <- sapply(data_scarlet[, 2:ncol(data_scarlet)], as.numeric)


#Calculate blank mean
meanblanck_OD <- as.data.frame(rowMeans(data_OD[grep("Blanc", names(data_OD), ignore.case=T)]))
meanblanck_scarlet <- as.data.frame(rowMeans(data_scarlet[grep("Blanc", names(data_scarlet), ignore.case=T)]))

#remove blank
for (i in 2:ncol(data_OD)){
  data_OD[,i] <- data_OD[,i] - meanblanck_OD[,1]
  data_scarlet[,i] <- data_scarlet[,i] - meanblanck_scarlet[,1]
}

colnames(data_OD) <- name_col
colnames(data_scarlet) <- name_col

#OD filter
export_data <- data.frame(matrix(NA,    # Create empty data frame
                                 nrow = 8,
                                 ncol = 0))
rownames(export_data) <- c("growth_rate","R²_GR","Scarlet_rate","R²_scarlet","mean_OD","mean_scarlet","ratio_scarletGR","ratio_meanScarletOD")
#####
for(x in 2:ncol(data_OD)){


 data_filter <- which.min(abs(data_OD[,x]-OD))
 if(data_filter< 4){
   data_filter = nrow(data_OD)/2
 }
 
 timelaps <- (data_OD[(data_filter-3):(data_filter+3),1])
 ODlaps <- na.omit(data_OD[(data_filter-3):(data_filter+3),x])
 scarletlaps <- na.omit(data_scarlet[(data_filter-3):(data_filter+3),x])

 blanc = FALSE
 
 for(r in 1:length(ODlaps)){
   if(ODlaps[r] <0){
     growth_rate = 0
     Scarlet_rate = 0
     mean_OD = 0
     mean_scarlet = 0
     Rsquare_scarlet = 0
     Rsquare_GR = 0
     ratio_meanScarletOD = 0
     Rsquare_GR = 0
     ratio_scarletGR = 0
     blanc = TRUE
     break
   }
 }
 if (blanc == FALSE) {
 sectime =c()
 sectimenorm = 0
 
 #time conversion
 for(y in 1:length(na.omit(timelaps))){
   sectime = c(sectime,abs(as.numeric(as.POSIXct(strptime(timelaps[y], "%Y-%m-%d %H:%M:%S")),tz = "GMT")))
   
 }
 df <- data.frame(sectime,ODlaps)
 diff_year = abs(df$sectime - lag(df$sectime))
 for(z in 2:length(sectime)){
   sectimenorm <- c(sectimenorm,sectimenorm[z-1]+diff_year[z])
 }
 
 #calcul
df <- data.frame(sectimenorm/3600,log(abs(ODlaps)),log(abs(scarletlaps)))
colnames(df) <- c("time","OD","scarlet")
df$OD[df$OD %in% c("-Inf")]<- 0
df$scarlet[df$scarlet %in% c("-Inf")]<- 0

df.lm <- lm(df$OD ~ df$time, data = df)
growth_rate <- summary(df.lm)$ coefficients[2]
if(growth_rate <0){
  growth_rate =0
}
Rsquare_GR <- summary(df.lm)$ r.squared
df.lm <- lm(df$scarlet ~ df$time, data = df)
Scarlet_rate <- summary(df.lm)$ coefficients[2]
if(Scarlet_rate <0){
  Scarlet_rate =0
}
Rsquare_scarlet <- summary(df.lm)$ r.squared
mean_OD <- mean(ODlaps)
mean_scarlet <- mean(scarletlaps)
ratio_scarletGR <- Scarlet_rate/growth_rate
if(growth_rate == 0){
  ratio_scarletGR =0
}
ratio_meanScarletOD <- mean_scarlet/mean_OD
 }
 cols <- which(grepl("Blanc|Dummy|blanc|dummy", names(export_data)))
 index <- c(2,4)
 export_data[index,cols] <- 1
#export data
export_data <- cbind(export_data,c(growth_rate,Rsquare_GR,Scarlet_rate,Rsquare_scarlet,mean_OD,mean_scarlet,ratio_scarletGR,ratio_meanScarletOD))

#assign colname
colnames(export_data)[ncol(export_data)]<- colnames(data_OD)[x]
}

plot <- as.data.frame(t(export_data))

name_row <- name_col[-1]
plot <- rownames_to_column(plot, var = "Souches")
plot$Souches <- name_row

# Filtrer les R²
plage_min <- 0.85

# Filtrer les données selon la plage minimale et maximale
plot <- plot[plot$`R²_GR` > plage_min,]
plot <- plot[plot$`R²_scarlet` > plage_min,]

filename_save<- paste(Sys.Date(),"_OD_",OD,"_export_",filename, sep = "")
write.xlsx(plot, file = filename_save ,sheetName="export", 
           col.names=TRUE, row.names=TRUE, append=FALSE)
if (exists("export_mixdata") == FALSE){
  export_mixdata <- c()
}
export_mixdata <- rbind(export_mixdata,plot)
}

filename_save<- paste(Sys.Date(),"_export_total_OD_",OD,".xlsx", sep = "")
write.xlsx(export_mixdata, file = filename_save ,sheetName="export", 
           col.names=TRUE, row.names=TRUE, append=FALSE)
}
setTxtProgressBar(pb, length(batch)*length(ODseq)+1)  # Sets the final progress bar 

close(pb) # Close the connection