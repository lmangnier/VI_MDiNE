library(mdine)

###########1. MICROBIOTE : Covariable + count des 25 microbiotes + abscence ou presence phenotype###########

###1. Mise en forme des donnees pour MDiNE

#1. Ajout de la colonne de ref 
#Prendre merge_sans_doublons + faire la somme en ligne de toutes les especes 
#hors celles qu'on etudie et faire une nouvelle colonne ref

###################25 premiers 

merge_sansdoublons$ref <- apply(merge_sansdoublons[, -c(1,2,3,4,5,6,7,8, 22, 128, 134, 135,136,137,276,383,384,
                                                        170,566,460,169,577,578,251,252,240,241,531,253,138,154,752,441)],1,sum,na.rm=TRUE)

###################5 premiers
merge_sansdoublons$ref <- apply(merge_sansdoublons[, -c(1,2,3,4,5,6,7,8, 22, 128, 134, 135,136)],1,sum,na.rm = TRUE)


###################5 derniers 20-25
merge_sansdoublons$ref <- apply(merge_sansdoublons[, -c(1,2,3,4,5,6,7,8, 253,138,154,752,441)],1,sum,na.rm = TRUE)

####################15 premiers 
merge_sansdoublons$ref <- apply(merge_sansdoublons[, -c(1,2,3,4,5,6,7,8, 22, 128, 134, 135,136,137,276,383,384,
                                                        170,566,460,169,577,578)],1,sum,na.rm=TRUE)



#2. Creation dtf pour MDiNE

#Especes selectionnes 
#covar_microbio <- merge_sansdoublons[, microbiote_25median[,"Species"]] #25 selectionnees 
#microbiote_25median[,"Species"]
#Covars
#covar_microbio <- merge_sansdoublons[, c(1:8)]

###############25 premiers 

#Lier l'ensemble : covars + especes selectionnees + ref 
covar_microbio <- merge_sansdoublons[, c(1,2,3,4,5,6,7,8, 22,128,134,135,136,137,276,383,384,
                                         170,566,460,169,577,578,251,252,240,241,531,253,138,154,752,441,941)]

which(colnames(merge_sansdoublons) == "k__Bacteria")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes|c__Clostridia")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_vulgatus")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Tannerellaceae")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Tannerellaceae|g__Parabacteroides")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae|g__Alistipes")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Roseburia")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Tannerellaceae|g__Parabacteroides|s__Parabacteroides_distasonis")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_caccae")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_ovatus")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Proteobacteria")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae")


#################5 premiers

covar_microbio5 <- merge_sansdoublons[, c(2,3,4,5,6,7,8,22,128,134,135,136,941)]

which(colnames(merge_sansdoublons) == "k__Bacteria")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales")
which(colnames(merge_sansdoublons) == "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae")


#################5 derniers 

covar_microbio5dernier <- merge_sansdoublons[, c(2,3,4,5,6,7,8, 253,138,154,752,441,941 )]

#################15 premiers 
covar_microbio15 <- merge_sansdoublons[, c(1,2,3,4,5,6,7,8, 22, 128, 134, 135,136,137,276,383,384,
                                            170,566,460,169,577,578,941)]


###2.MDiNE

#######################25 premiers
sum(is.na(covar_microbio)) #440 NA

#On garde que les lignes completes car model.matrix enleve les lignes avec NA
covar_microbio_complet <- covar_microbio[complete.cases(covar_microbio[,]),]
str(covar_microbio_complet) 

#X : #intercept + disease + covar (intercept + disease dans intro)
X <- model.matrix(~diagnosis+sex+Antibiotics+Chemotherapy+Age.at.diagnosis, data = covar_microbio_complet)

#Y : count OTU
Y <- covar_microbio_complet[,9:34]
Y <- as.matrix(Y)
sum(Y == 0)

#Condition du Y  
#any(floor(Y)!=Y) #pas de donnees entieres
#is.matrix(Y) #doit etre une matrice 
#which(any(Y<0)) #positive 

#Y doit etre entiers
for (i in c(9:34)) {
  covar_microbio_complet[,i] <- ceiling(covar_microbio_complet[,i])
}

for (i in c(1:26)) {
  Y[,i] <- ceiling(Y[,i])
}

sum(Y == 0)
sum(covar_microbio_complet == 0) #0 du diagnosis en plus 

md.fit <- mdine(Y=Y[index.random,], X = X[index.random,1:2], Z = X[index.random,2], mc.cores = 5, iter = 1000) 
index.random = sample(1:nrow(Y), 100, replace = FALSE) 

plot_networks(md.fit)

########################5 plus exprimes 

covar_microbio5_complet <- covar_microbio5[complete.cases(covar_microbio5[,]),]

X5 <- model.matrix(~diagnosis+sex+Antibiotics+Chemotherapy+Age.at.diagnosis, data = covar_microbio5_complet)

Y5 <- covar_microbio5_complet[,8:13]
Y5 <- as.matrix(Y5)

for (i in c(8:13)) {
  covar_microbio5_complet[,i] <- ceiling(covar_microbio5_complet[,i])
}

for (i in c(1:6)) {
  Y5[,i] <- ceiling(Y5[,i])
}

sum(Y5 == 0)
sum(X5 == 0)

sum(covar_microbio5_complet == 0) 

md.fit5 <- mdine(Y=Y5[index.random,], X = X5[index.random,1:3], Z = X5[index.random,2]) 
index.random = sample(1:nrow(Y5), 100, replace = FALSE) 

plot_networks(md.fit5) 


#####################5 moins exprimes 
covar_microbio5dernier_complet <- covar_microbio5dernier[complete.cases(covar_microbio5dernier[,]),]

X5dernier <- model.matrix(~diagnosis+sex+Antibiotics+Chemotherapy+Age.at.diagnosis, data = covar_microbio5dernier_complet)

Y5dernier <- covar_microbio5dernier_complet[,8:13]
Y5dernier <- as.matrix(Y5dernier)

for (i in c(1:6)) {
  Y5dernier[,i] <- ceiling(Y5dernier[,i])
}

md.fit5dernier <- mdine(Y=Y5dernier[index.random,], X = X5dernier[index.random,1:6], Z = X5dernier[index.random,2]) #sans covar
index.random = sample(1:nrow(Y5), 100, replace = FALSE) 

plot_networks(md.fit5dernier)

sum(Y5dernier == 0) 


##############15 premiers 
covar_microbio15_complet <- covar_microbio15[complete.cases(covar_microbio15[,]),]

X15 <- model.matrix(~diagnosis+sex+Antibiotics+Chemotherapy+Age.at.diagnosis, data = covar_microbio15_complet)

Y15 <- covar_microbio15_complet[,9:24]
Y15 <- as.matrix(Y15)

for (i in c(1:16)) {
  Y15[,i] <- ceiling(Y15[,i])
}

md.fit15 <- mdine(Y=Y15[index.random,], X = X15[index.random, c(1,2,3,4,6)], Z = X15[index.random,2]) 
index.random = sample(1:nrow(Y5), 500, replace = FALSE) 

plot_networks(md.fit15)

sum(Y15 == 0) 




###########2. METABOLITE : Covariable + count des 25 microbiotes + abscence ou presence des 25 metabolites###########


#1. Mise en forme des donnees pour MDiNE


