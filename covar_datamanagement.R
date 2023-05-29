###########Donnees covariables###########
data_covar <- read.csv('hmp2_metadata.csv')
#5533 obs, 490 var avec patient en lignes

class(data_covar)
str(data_covar)
colnames(data_covar)

#Selection covariable 
covar <- subset(data_covar, select = c(Project, External.ID, Participant.ID, sex , Antibiotics, Chemotherapy, Age.at.diagnosis, diagnosis))
str(covar)

#Diagnostic en binaire 
for (i in c("diagnosis")) {
  covar[,i] <- ifelse(covar[,i]=="CD",1,0)
}

#Tout en facteur 
for (i in c("Antibiotics", "sex", "diagnosis", "Chemotherapy" )) {
  covar[,i] <- as.factor(covar[,i])
}
str(covar)

#NA quand vides 
covar[covar== "" ] <- NA 
sum(is.na(covar)) #1384 NA dans age diagnosis 

###Stats descriptives : jeu covars complet 

##Univaries 
#Categorielle 
table(covar$sex) #equilibre
table(covar$sex) *100 / nrow(covar)
barplot(table(covar$sex))

table(covar$Antibiotics) #no
table(covar$Antibiotics) *100 / nrow(covar)
barplot(table(covar$Antibiotics))

table(covar$Chemotherapy) #no
table(covar$Chemotherapy) *100 / nrow(covar)
barplot(table(covar$Chemotherapy))

table(covar$diagnosis)
table(covar$diagnosis) *100 / nrow(covar)
barplot(table(covar$diagnosis))

#Continue 
summary(covar$Age.at.diagnosis)
sd(covar$Age.at.diagnosis, na.rm = TRUE)
hist(covar$Age.at.diagnosis) #pas de loi normale 


##Bivariees : 

#En fonction du diagnosis 
#Categorielle : reg log pour avoir des OR 
summary(glm(diagnosis~sex, data = covar, family = "binomial")) #signif
summary(glm(diagnosis~Antibiotics, data = covar, family = "binomial")) #signif
summary(glm(diagnosis~Chemotherapy, data = covar, family = "binomial")) #signif 

#Continue 
by(covar$Age.at.diagnosis, covar$diagnosis, sd, na.rm = TRUE)
t.test(covar$Age.at.diagnosis ~ covar$diagnosis, var.equal = TRUE) #signif 
#pertinence avec autant de sujet ?

#et des variables entre elles 
summary(glm(sex~Antibiotics, data = covar, family = "binomial")) #signif
summary(glm(sex~Chemotherapy, data = covar, family = "binomial")) #signif
t.test(covar$Age.at.diagnosis ~ covar$sex, var.equal = TRUE) #signif

summary(glm(Antibiotics~Chemotherapy, data = covar, family = "binomial")) #non signif 
t.test(covar$Age.at.diagnosis ~ covar$Antibiotics, var.equal = TRUE) #signif

t.test(covar$Age.at.diagnosis ~ covar$Chemotherapy, var.equal = TRUE) #non signif 


##################################Remarque
#correlation des covariables
#especes independantes les unes des autres 

##################################
