###########Donnees microbiote###########
data_microbio <- read.table(file = 'hmp2_mgx_taxonomy.tsv', sep = '\t', header = TRUE)
colnames(data_microbio)

#On retourne pour que ce soit pareil que covars : 
#On veut les microbiotes en colonne et les indiv en lignes

#On retourne le data  
data_microbio_retourne <- t(data_microbio)
data_microbio_retourne<-as.data.frame(data_microbio_retourne)
colnames(data_microbio_retourne) <- data_microbio_retourne[1,]
data_microbio_retourne <- data_microbio_retourne[-1,]

#NA si cases vide 
data_microbio_retourne[data_microbio_retourne== "" ] <- NA 
sum(is.na(data_microbio_retourne)) #0

#Creation d'une colonne ExternalID pour la jointure avec covars 
data_micro_merge <- data_microbio_retourne 
data_micro_merge$External.ID <- row.names(data_microbio_retourne)

#Enlever le "_profile" pour pouvoir merger avec covars 
install.packages("tidyverse")
library(tidyverse)
data_micro_merge$External.ID <- str_replace(data_micro_merge$External.ID, "_profile", "")


##########1. Covariable + count des 25 microbiotes + abscence ou presence phenotype###########

###1.1 On merge data_micro_merge et covars sur external id

duplicate <- covar$External.ID[duplicated(covar$External.ID)] 
#2641 doublons de External ID dans covar

#On tri sur external id 
#covar <- covar[order(covar$External.ID), ]

duplicate2 <- data_micro_merge$External.ID[duplicated(data_micro_merge$External.ID)]
#Pas de doublons de External id dans les donnes de microbiote

#Merger les fichiers : avec les duplicates
merge_duplicate <- merge(covar, data_micro_merge, by = 'External.ID', sort = TRUE, all = FALSE)
#3826 indiv = 1638 de microbiote + 2188 doublons de covar et 940 variables
sum(is.na(merge_duplicate)) #1032
sum(is.na(merge_duplicate$Age.at.diagnosis)) #1032 : tous les NA sont dans cette variable 

#Merger les fichiers : sans les duplicates
duplicate3 <- which(duplicated(merge_duplicate$External.ID)) #trouver les lignes des duplicates
merge_sansdoublons <- merge_duplicate[-duplicate3,] #1638 lignes, 940 variables 

sum(is.na(merge_sansdoublons)) #440
sum(is.na(merge_sansdoublons$Age.at.diagnosis)) #tous les NA dans cette variable



#Stat sur jeu restreint 
#Voir si les distrib des covar sont les memes que dans l'echantillon entier 
table(merge_sansdoublons$sex)
table(merge_sansdoublons$sex) *100 / nrow(merge_sansdoublons)
table(merge_sansdoublons$Antibiotics)
table(merge_sansdoublons$Antibiotics) *100 / nrow(merge_sansdoublons)
table(merge_sansdoublons$Chemotherapy)
table(merge_sansdoublons$Chemotherapy) *100 / nrow(merge_sansdoublons)
table(merge_sansdoublons$diagnosis)
table(merge_sansdoublons$diagnosis) *100 / nrow(merge_sansdoublons)
summary(merge_sansdoublons$Age.at.diagnosis)

#En fonction du diagnosis 
#Categ
summary(glm(diagnosis~sex, data = merge_sansdoublons, family = "binomial")) #signif
summary(glm(diagnosis~Antibiotics, data = merge_sansdoublons, family = "binomial")) #signif
summary(glm(diagnosis~Chemotherapy, data = merge_sansdoublons, family = "binomial")) #signif 
#Continue 
by(merge_sansdoublons$Age.at.diagnosis, merge_sansdoublons$diagnosis, sd, na.rm = TRUE)
t.test(merge_sansdoublons$Age.at.diagnosis ~ merge_sansdoublons$diagnosis, var.equal = TRUE) #signif 
#pertinence avec autant de sujet ??

#et des variables entre elles 
summary(glm(sex~Antibiotics, data = merge_sansdoublons, family = "binomial")) #non signif
summary(glm(sex~Chemotherapy, data = merge_sansdoublons, family = "binomial")) #signif
t.test(merge_sansdoublons$Age.at.diagnosis ~ merge_sansdoublons$sex, var.equal = TRUE) #signif

summary(glm(Antibiotics~Chemotherapy, data = covar, family = "binomial")) #non signif 
t.test(covar$Age.at.diagnosis ~ covar$Antibiotics, var.equal = TRUE) #signif

t.test(covar$Age.at.diagnosis ~ covar$Chemotherapy, var.equal = TRUE) #non signif 


###1.2. Trouver les microbiotes les plus discriminants   

#Passer de char a num pour les count de microbio
for (i in c(9:940)) {
  merge_sansdoublons[,i] <- as.numeric(merge_sansdoublons[,i])
}

class(merge_sansdoublons)
head(str(merge_sansdoublons))

#Normalite des variables 
#log et +1 pour le pseudo count 
boxplot(log(merge_sansdoublons$`k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae`+1) ~ merge_sansdoublons$diagnosis) #OK

hist((merge_sansdoublons$k__Archaea))
hist(log(merge_sansdoublons$k__Archaea+1)) 
#pas de LN, on fera des wilcoxon 


#1: Choix des microbiotes independemment du phenotype
#Print moyennes/medianes des count des microbiotes et prendre les 25 premiers 

#Mediane de chaue colonne
#apply(merge_sansdoublons[,-c(1:8)],MARGIN = 2, median)
#microbiote_25median <- names(head(sort(apply(merge_sansdoublons[,-c(1:8)],MARGIN = 2, median), decreasing = TRUE), n.fam))

microbiote_25median <- data.frame()
k <- 1
for (i in c(9:940)) {
  median <- median(merge_sansdoublons[,i])
  microbiote_25median[k,1] <- colnames(merge_sansdoublons[i])
  microbiote_25median[k,2] <- median
  k <- k+1
}

n.fam = 25 
names(microbiote_25median) <- c("Species", "Median")
microbiote_25median <- microbiote_25median[order(microbiote_25median$Median, decreasing = TRUE), ] #on trie 
microbiote_25median <- microbiote_25median[1:n.fam,] #on prend les 25 premiers 

microbiote_5median <- microbiote_25median[1:5,] #5 premiers 
microbiote_2025median <- microbiote_25median[21:25,]
microbiote_15median <- microbiote_25median[1:15,]

#Moyenne de chaque colonnes
#apply(merge_sansdoublons[,-c(1:8)],MARGIN = 2, mean)
#microbiote_25mean <- names(head(sort(apply(merge_sansdoublons[,-c(1:8)],MARGIN = 2, mean), decreasing = TRUE), n.fam))

microbiote_25mean <- data.frame()
k <- 1
for (i in c(9:940)) {
  mean <- mean(merge_sansdoublons[,i])
  microbiote_25mean[k,1] <- colnames(merge_sansdoublons[i])
  microbiote_25mean[k,2] <- mean
  k <- k+1
}

names(microbiote_25mean) <- c("Species", "Mean")
microbiote_25mean <- microbiote_25mean[order(microbiote_25mean$Mean, decreasing = TRUE), ] #on trie
microbiote_25mean <- microbiote_25mean[1:n.fam,] #on prend les 25 


#2: Choix des microbiotes dependemment du phenotype
#Wilcoxon phenotype ~ species colonnes et prendre les 25 pvalue les plus petites 

microbiote_25wilcoxon <- data.frame()
k <- 1
for (i in c(9:940)) {
  test_wilcox <- wilcox.test(merge_sansdoublons[,i] ~ diagnosis, data = merge_sansdoublons)
  microbiote_25wilcoxon[k,1] <- colnames(merge_sansdoublons[i])
  microbiote_25wilcoxon[k,2] <- test_wilcox$p.value
  k <- k+1
  }

names(microbiote_25wilcoxon) <- c("Species", "pvalue")
microbiote_25wilcoxon <- microbiote_25wilcoxon[order(microbiote_25wilcoxon$pvalue), ]
microbiote_25wilcoxon <- microbiote_25wilcoxon[1:n.fam,]


#3. Recap deux methodes 

#Tres similaire 
microbiote_25medianbis
microbiote_25meanbis

microbiote_25wilcoxon #Assez different 
