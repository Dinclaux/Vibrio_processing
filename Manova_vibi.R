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

if (!require("MASS")){
  install.packages("MASS", dependencies = TRUE)
  library(MASS)
}

#########################################################
###              working directory                    ###
#########################################################


setwd("C:/Users/mdinclaux/Documents/Script/TWB_Processing")

rawdata <- read.xlsx("2022-09-30_Combinatoire VIBI_CONSTITUTIF.xlsx", header = TRUE, sheetName = "Bilan Constitutif")

#remove NA
rawdata <- rawdata[!is.na(rawdata[,3]),]

rawdata <- rawdata[,-1:-2]

rawdata$Promoteur <- factor(rawdata$Promoteur)
rawdata$RBS <- factor(rawdata$RBS)
rawdata$Terminateur <- factor(rawdata$Terminateur)

vibi.mod <- lm(cbind(growth_rate,Scarlet_rate,ratio_scarletGR)~Promoteur+Terminateur+RBS+Promoteur*Terminateur*RBS -1, data = rawdata)
options(max.print = 1000000)
sink("lm.txt")
print(summary(vibi.mod))
sink()


manoval_model <- manova(cbind(growth_rate,Scarlet_rate,ratio_scarletGR,ratio_meanScarletOD)~ Promoteur+RBS+Terminateur+Promoteur*RBS*Terminateur, data = rawdata)
summary(manoval_model)
summary.aov(manoval_model)
post_hoc <- lda( rawdata$RBS~cbind(rawdata$growth_rate,rawdata$Scarlet_rate,rawdata$ratio_scarletGR), CV=F)
summary(post_hoc)
dependent_vars <- cbind(rawdata$growth_rate,rawdata$Scarlet_rate)
independent_var <- cbind(rawdata$Promoteur,rawdata$RBS,rawdata$Terminateur)

plot_lda <- data.frame(rawdata[, "RBS"], lda = predict(post_hoc)$x)
ggplot(plot_lda) + geom_point(aes(x = lda.LD1, y = lda.LD2, colour = rawdata$RBS ), size = 4)


res.pca<-PCA(dependent_vars, scale.unit = TRUE, ncp = 6, graph = TRUE)
var <- get_pca_var(res.pca)
fviz_pca_biplot(res.pca,col.ind = rawdata$Promoteur, label = "var",addEllipses = TRUE)


MyResult.splsda <- splsda(rawdata$growth_rate, rawdata$Scarlet_rate)