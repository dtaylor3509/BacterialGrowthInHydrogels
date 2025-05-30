#Relevant libraries
library(readxl)


RTypeFull <- read_excel("RTypeFull2.xlsx")
growthA = RTypeFull


# Correlations ------------------------------------------------------------

#Overall
cor.test(growthA$newstiffness,growthA$newval2,method = "spearman")

#For each species of bacteria
for(i in unique(growthA$species)){
  print(i)
  print(cor.test(subset(growthA,species==i)$newstiffness,subset(growthA,species==i)$newval2,method = "spearman"))
}

#For each species and hydrogel type
for(i in unique(growthA$species)){
  for(j in unique(growthA$temperature)){
    print(paste(i,j))
    print(cor.test(subset(growthA,species==i & temperature==j)$newstiffness,subset(growthA,species==i & temperature==j)$newval2,method = "spearman"))
  }
}

cor.test(growthA$newwaterloss,growthA$newval2,method = "spearman",data=growthA)

for(i in unique(growthA$species)){
  print(i)
  print(cor.test(subset(growthA,species==i)$newwaterloss,subset(growthA,species==i)$newval2,method = "spearman"))
}

for(i in unique(growthA$species)){
  for(j in unique(growthA$temperature)){
    print(paste(i,j))
    print(cor.test(subset(growthA,species==i & temperature==j)$newwaterloss,subset(growthA,species==i & temperature==j)$newval2,method = "spearman"))
  }
}



# Factor Analysis ---------------------------------------------------------

#Bacterial species
growthA$species = relevel(as.factor(growthA$species),ref = "E. Coli")

#Melting temperature/Hydrogel type
growthA$temperature = as.factor(growthA$temperature)

#Nutrient media
growthA$nm = relevel(as.factor(growthA$nm),ref = "TSB")

#Concentration
growthA$conc = as.factor(growthA$conc)

#Replicate count
growthA$rep = as.factor(growthA$rep)

#Indicator for gram-positive or gram-negative bacteria
growthA$gram = as.factor(growthA$gram)

growthA$EColi = as.factor(ifelse(growthA$species == 'E. Coli', 1, 0))
growthA$PFlu = as.factor(ifelse(growthA$species == 'P. Flu', 1, 0))
growthA$SAureus = as.factor(ifelse(growthA$species == 'S. Aureus', 1, 0))
growthA$BSubtilis = as.factor(ifelse(growthA$species == 'B. Subtilis', 1, 0))

summary(aov)

growth0s = subset(growthA,growthA$value == 0)
growth = subset(growthA,growthA$value != 0)  #using the strictly positive dataset


fact = lm(newval2 ~ species*temperature*nm*conc - species:temperature:nm:conc,data=growthA)
summary(aov(fact))

# Unused - extension work ---------------------------------------------------------

library("glmm")
library(glmnet)
library("gglasso")
library("glinternet")
library("grpreg")

#log and root analyses
growthLog = growth
growthLog$value = log(growth$value)
growthLog$logInv = log(growth$inv)

growthRoot = growthA
growthRoot$value = growth$value^0.5

#Creating the linear model for OLS estimates and model matrices



lmmodel = lm(value ~ PFlu*(temperature+nm+conc)^3 - PFlu:temperature:conc:nm - PFlu:temperature:conc
            + SAureus*(temperature+nm+conc)^3 - SAureus:temperature:conc:nm - SAureus:temperature:conc
            + BSubtilis*(temperature+nm+conc)^3 - BSubtilis:temperature:conc:nm - BSubtilis:temperature:conc
            , data=growthLog)

ggmatI = model.matrix(~ PFlu*(temperature+nm+conc)^3 - PFlu:temperature:conc:nm - PFlu:temperature:conc
                      + SAureus*(temperature+nm+conc)^3 - SAureus:temperature:conc:nm - SAureus:temperature:conc
                      + BSubtilis*(temperature+nm+conc)^3 - BSubtilis:temperature:conc:nm - BSubtilis:temperature:conc
                      , data=growthLog)
ggmat = ggmatI[,-1]

Icoefs = lmmodel$coefficients #The OLS estimates including the intercept
coefs = lmmodel$coefficients[-1] #The OLS estimates excluding the intercept
coefs[is.na(coefs)] = 0
p = length(coefs)

#Generating the group structure of the coefficients

groups = lmmodel$assign[-1]
groupv = factor(groups)
instancemat = model.matrix(~groupv+0)

#calculating adaptive weights
w = c()
bn = c()

for(i in 1:max(groups)){
  w[i] = norm(coefs*instancemat[,i],type = "2")^-1
  bn[i] = sum(instancemat[,i])
}

nLambda2 = 2000

#Using the 'grpreg' package to create the regularization path of the adaptive group lasso
#with the model matrix and group structure of the linear model and adaptive weights

Pgglass = grpreg(ggmat,growthLog$value,group=groups,group.multiplier =100*w,
                 family = "gaussian",intercept=TRUE,penalty = "grLasso",nlambda=nLambda2,log.lambda = TRUE,returnX = TRUE,lambda.min = 1e-5,eps = 1e-2)


nLambda = length(Pgglass$lambda)

#Plotting the regularization path
plot(Pgglass,log.l = TRUE,label = FALSE,norm=TRUE,xlab="log(100λ)")
abline(v = log(bestlambda)) #Used once the optimal lambda has been calculated

groupLabels = attr(terms(lmmodel),"term.labels")
leg = legend("topleft",legend = (groupLabels),text.width = 2,cex=0.44,col = rainbow(14),lty = rep(1,14),y.intersp=0.5)


#predictions & BIC
logPreds = matrix(0,nrow = 300,ncol=nLambda) #added
Preds = matrix(0,nrow = 300,ncol=nLambda) #added 
agResids = matrix(0,nrow = 300,ncol=nLambda)
RSS = c()
df = c()
BICpen = c()

BICpenAGL = matrix(0,nrow = nLambda,ncol=max(groups))
BICpenOLS = matrix(0,nrow = nLambda,ncol=max(groups))
df1 = c()
df2 = c()
agBIC = c()
agN = length(ggmat[,1])

#Algorithm to calculate the degrees of freedom for all models in the regularization path
#As defined by Yuan and Lin

for(i in 1:length(Pgglass$lambda)){
  logPreds[,i] = ggmatI%*%(coef(Pgglass))[,i] #[,i] is i'th model
  agResids[,i] = logPreds[,i] - growthLog$value
  RSS[i] = norm(agResids[,i],type = "2")^2
  
  for(j in 1:max(groups)){ #model i group j
    BICpenAGL[i,j] = norm(coef(Pgglass)[-1,i]*instancemat[,j],type = "2")
    BICpenOLS[i,j] = norm(coefs*instancemat[,j],type = "2")
  }
  df1[i] = length(BICpenAGL[i,][BICpenAGL[i,]>0])
  df2[i] = sum(BICpenAGL[i,]*(bn-1)/BICpenOLS[i,])
  agBIC[i] = log(RSS[i]/agN)+log(agN)*(df1[i] +df2[i])/agN
}

pres = c() #for 'presence', indicating the number of times a group's norm is non-zero
for(i in max(groups)){
  pres[i] = sum(BICpenAGL[,i] == 0)
}
groupLabels[order(pres)] #obtains order in which variables are added to the regularization path


bestindex = order(agBIC)[1]
bestindex2 = order(BIC(Pgglass))[1] #Using the BIC as defined by Breheny and Luang
bestlambda = Pgglass$lambda[agBIC == min(agBIC)]
bestlambda2 = Pgglass$lambda[BIC(Pgglass) == min(BIC(Pgglass))]
logofbestlambda = log(Pgglass$lambda[agBIC == min(agBIC)])
bestcoefs = coef(Pgglass)[,bestindex]

bestlambdaLog = bestlambda

#Plotting regularization parameter against BIC
plot(log(Pgglass$lambda*100),agBIC,type='l',ylab = "BIC",xlab = "log(λ)",col="blue")

#Calculates only one model for the optimal lambda
logbestModel =  grpreg(ggmat,growthLog$value,group=groups,group.multiplier =100*w, lambda = bestlambdaLog,
                                family = "gaussian",intercept=TRUE,penalty = "grLasso",nlambda=nLambda2,log.lambda = TRUE,returnX = TRUE,lambda.min = 1e-5,eps = 1e-2)

allLogpredictions = predict(logbestModel,X = ggmat)
logModelPredictions = exp(allLogpredictions)
logResids = logModelPredictions - exp(growthLog$value)
logRSS = sum(logResids^2)

#Variables indicating which groups are set to zero in the optimal model
removed = groupLabels[BICpenAGL[bestindex,] == 0]
kept = groupLabels[BICpenAGL[bestindex,] != 0]


