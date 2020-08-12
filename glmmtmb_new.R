
rm(list = ls())
  
  
  
#library(DHARMa)
library(usdm)
require(MuMIn)
require(glmmTMB)
#library(here)
library(ggplot2)
library(sjPlot)
require(parallel)
library(performance)
library(see)
library(randomForest)


dat<-read.csv("C:/Users/Maxwell Azali/Dropbox/Global Bleaching/Paper/Resistance/Data/Resistance data final.csv", stringsAsFactors = FALSE)

colnames(dat) [1] <- "unique.id"

# Use percent bleached and division as resistance metric
respVar<-dat[, c(1:4,10,11,19:21)]

PredictVar<-dat[,c(12:14,16:18)]

#Standardize continous predictor variables
data.std<-apply(X = PredictVar, MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))})

data.stds<-data.frame(cbind(respVar,data.std,dat[,c("Region","CoralProvince","Ecoregion","Province")]))

# Missing Hard coral values- 20 sites Australia(Ningaloo and Kimberlys); 1 site Reunion; M
data.stds[!complete.cases(data.stds),]

data.complete <- na.omit(data.stds)

# Tabulate spatial groupings
with(data.complete,table(CoralProvince))

with(data.complete,table(Province))

with(data.complete, table(Ecoregion))

with(data.complete, table(Region))

# Set reference levels for analysis
data.complete$Province <- as.factor(data.complete$Province)
data.complete$Province <- relevel(data.complete$Province, "Western Indian Ocean")

data.complete$CoralProvince <- as.factor(data.complete$CoralProvince)
data.complete$CoralProvince <- relevel(data.complete$CoralProvince, "Africa-India")

data.complete$Ecoregion <- as.factor(data.complete$Ecoregion)
data.complete$Ecoregion <- relevel(data.complete$Ecoregion, "East African Coral Coast")

data.complete$Region <- as.factor(data.complete$Region)
data.complete$Region <- relevel(data.complete$Region, "Non Coral triangle")

data.complete$pos <- numFactor(scale( data.complete$Longitude ), scale( data.complete$Latitude ))

data.complete$ID <- factor(rep(1, nrow(data.complete)))

summary(data.complete)

#PredictVar <- data.complete[,c("Hard.coral.perc.", "No..of.genera", "SST", "Kurtosis", "Skewness", "Region", "CoralProvince", "Ecoregion")]

#PredictVar <- data.complete[,c("Hard.coral.perc.", "No..of.genera", "SST", "Kurtosis", "Skewness", "Region")]

PredictVar <- data.complete[,c("Hard.coral.perc.", "No..of.genera", "SST", "Kurtosis", "Skewness", "CoralProvince")]

#### Functions 

##McClanahan et al. “Novel temperature patterns and geographic context critical to coral bleaching during the 2016 El Niño” 
##JosephMaina_2018_19
#this function evaluates the models and returns a model selection object
#for poisson models 
evalTextModel  <-  function(modelText,...) {
  #theMod <- model.sel(try(eval(parse(text = modelText))), extra =c(AIC, BIC, Cp,Nc = function(x) mean(ID)))
  #theMod <- model.sel(try(eval(parse(text = modelText))), extra =c(AIC, BIC, Cp,Nc = function(x) ModelRef(modelTextList, modelText)[[1]]))
  theMod <- try(eval(parse(text = modelText))) 
  tryErrorCheck  <-  all(class(theMod) %in% 'try-error')
  if(!tryErrorCheck)
    theMod
  #model.sel(theMod, rank.args = list(REML = FALSE), extra =c(AIC, BIC, Cp))
  #model.sel(theMod, rank=QAIC, rank.args=list(chat = deviance(theMod$mer) / df.residual(theMod$mer), REML = NULL), extra =c(AIC, BIC, Cp,OvDispPval = function(x) overdisp_fun(x$mer)[[4]]))
}

#arcsin squareroot transformation formula for percentages
trans.arcsine <- function(x){
  asin(sign(x) * sqrt(abs(x)))
}



checkNewVarAndAdd <-  function(originalData, targetData, newVar) {
  responseAndPredictorCheck  <-  sizesMatch(originalData, targetData)
  newVarNameDoesNotExist     <-  !(newVar %in% names(targetData))
  
  if(responseAndPredictorCheck & newVarNameDoesNotExist)
    originalData[[newVar]]
}

sizesMatch  <-  function(originalData, targetData) {
  nrow(originalData) == nrow(targetData)
}

runPredCombinations  <-  function(x, vec) {
  lapply(data.frame(combn(vec, x), stringsAsFactors=FALSE), identity)
}


##prepare the model text



prepareModelText  <-  function(Predictors, data) {
  ##dat           <-  data[, Predictors]
  dat <- subset(data, select = Predictors)
  dataClasses   <-  unique(sapply(dat, class))
  numOfClasses  <-  length(dataClasses)
  mixedClasses  <-  numOfClasses > 1
  if(mixedClasses) {
    fixStepNum   <-  paste0(Predictors[sapply(dat, class) == 'numeric'], collapse='+')
    fixStepCha   <-  paste0(Predictors[sapply(dat, class) == 'character'], collapse='+')
    fixStep      <-  paste0(fixStepNum, '+', fixStepCha)
  } else {
    if(dataClasses == 'numeric') {
      fixStep   <-  paste0(Predictors[sapply(dat, class) == 'numeric'], collapse='+')
      ranStep1  <-  paste0('(0+', Predictors ,  collapse='+')
    }
    if(dataClasses == 'character') {
      fixStep   <-  paste0(Predictors[sapply(dat, class) == 'character'], collapse='+')
    }
  }
  
  # Change predictors according to model
  #paste0('glmmTMB(formula = Res_CE_pbl_div~', fixStep,'+ exp(pos + 0 | ID)', ', family=Gamma(link=log)' ,', data=data.complete)')
 paste0('glmmTMB(formula = Res_DHM_pbl_div~', fixStep,'+ exp(pos + 0 | ID)', ', family=Gamma(link=log)' ,', data=data.complete)')
  
}

runModelText  <-  function(...) {
  modelText  <-  prepareModelText(...)
  evalTextModel(modelText)
}

###to test the list of predictors
#lapply(yourListOfPredictors, prepareModelText, data=yourDataset)to call this to the dataset

#getMaximumNOfCombs  <-  function(vecOfPredictors, combCeiling=3)##here adjust the number of variable combinations for fixed effects
getMaximumNOfCombs  <-  function(vecOfPredictors, combCeiling=10) {
  if(length(vecOfPredictors) >= combCeiling) {
    combCeiling
  } else {
    length(vecOfPredictors)
  }
}

##correlation matrix and table
corstarsl <- function(x){ 
  require(Hmisc) 
  require(xtable)
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
} 



### Analyses

vifPredCombinations  <-  list()
varnames<-colnames(PredictVar)#

maxCombs  <-  getMaximumNOfCombs(varnames)
for(j in 1:maxCombs) {
  vifPredCombinations  <-  append(runPredCombinations(j, varnames), vifPredCombinations)
}

##2. filter the combinations above with VIF<1.5
vifPredCombinations_new<- c()
for(con in vifPredCombinations){
  r <- subset(PredictVar, select = con)
  conClasses   <-  unique(sapply(r, class))
  numOfClasses  <-  length(conClasses)
  twoNCols <- ncol(r)==2
  numDf <- r[sapply(r,is.numeric)]
  zeroNumDf<-ncol(numDf)==0
  numeriCols<- ncol(numDf)>1
  
  if(length(con)<2){
    vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
  }else{
    if (zeroNumDf) { 
      vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
    }
    if (numOfClasses==2 && twoNCols) { 
      vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
    }
    if(numeriCols && max(vif(numDf)["VIF"])< 3){##vif cutoff
      vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
    }
    next
  }
}



##ns refers to models wirh no random slopes
modelText<-lapply(vifPredCombinations_new, prepareModelText, data.stds)

modelText.a<-unlist(modelText)
#modelText.b<-as.list(modelText.a)
##

#mods1<-mclapply(modelText.a[1:50], evalTextModel)
#mods2<-mclapply(modelText.a[1:47], evalTextModel)
#mods3<-mclapply(modelText.a[1:47], evalTextModel)
mods4<-mclapply(modelText.a[1:47], evalTextModel)

#modList.1<-mods1[!sapply(mods1, is.na(mods1))] 


modelSel1<-model.sel(mods1[1:47], REML=FALSE, rank = "AICc")
modelSel2<-model.sel(mods2, REML=FALSE, rank = "AICc")
modelSel3<-model.sel(mods3, REML=FALSE, rank = "AICc")
modelSel4<-model.sel(mods4, REML=FALSE, rank = "AICc")


#write.csv(modelSel1, file ="C:/Users/Maxwell Azali/Dropbox/Global Bleaching/Paper/Resistance/Data/Modelsel_Resistance_CE_Coraltriangle.csv" )

write.csv(modelSel2, file ="C:/Users/Maxwell Azali/Dropbox/Global Bleaching/Paper/Resistance/Data/Modelsel_Resistance_DHM_Coraltriangle.csv" )
write.csv(modelSel3, file ="C:/Users/Maxwell Azali/Dropbox/Global Bleaching/Paper/Resistance/Data/Modelsel_Resistance_CE_CoralProvince.csv" )
write.csv(modelSel4, file ="C:/Users/Maxwell Azali/Dropbox/Global Bleaching/Paper/Resistance/Data/Modelsel_Resistance_DHM_CoralProvince.csv" )

