rm(list = ls())
#setwd('/Users/josephmaina/Dropbox/Global Bleaching/Paper/Resistance/Data/')
setwd('~/Dropbox/Global Bleaching/Paper/Resistance/Data/')




#dat<-read.csv('/Users/josephmaina/Dropbox/Global Bleaching/Paper/Resistance/Data/Resistance_data_new.csv')

dat<-read.csv('~/Dropbox/Global Bleaching/Paper/Resistance/Data/Resistance_revised.csv')
#dat<-read.csv('~/Dropbox/Global Bleaching/Paper/Resistance/Data/Resistance_data_new.csv')

#source('functions_analyses_glmmTMB.R')


source("~/Dropbox/Global Bleaching/Paper/Resistance/Drafts General Info/CoralAdaptiveCapacityData/functions_analyses_azali_edit_glmmTMB.R")

#transformationa ccording to Daniel Zimprich
colnames(dat)
respVar<-dat[, c(6:8, 18:21)]

PredictVar<-dat[,c(9,10,11,12,13,15:17)]

#scaling <- function(y){ 
 # a<-min(y)
  #b<-max(y)
  #N<-length(y)
  #y_a<-y-a
  #b_a<-(b-a)
  #yscaled<-((y_a*(N-1))/(b_a*N))+(1/(2*N))
  #return(yscaled) 
#}

#resp.std<-scaling(respVar)

#resp.std<-apply(X = respVar, MARGIN = 2,
 #               FUN = function(x) scaling(x))

colnames(respVar)
  
##define predictors
colnames(dat)


data.std<-apply(X = PredictVar, MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))})

##




##combine variables for each predictor
#1. Create list with all possible combinations between predictors 
vifPredCombinations  <-  list()
varnames<-colnames(PredictVar)#

maxCombs  <-  getMaximumNOfCombs(varnames)
for(j in 1:maxCombs) {
  vifPredCombinations  <-  append(runPredCombinations(j, varnames), vifPredCombinations)
}

##2. filter the combinations above with VIF<1
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
    if(numeriCols && max(vif(numDf)["VIF"])<0){##vif cutoff
      vifPredCombinations_new <- c(vifPredCombinations_new, list(con))
    }
    next
  }
}

#data.stds<-cbind(resp.std,data.std)
data.stds<-data.frame(cbind(respVar,data.std))
#convert all to numeric
data.stds<-data.frame(apply(data.stds, 2, function(x) as.numeric(x)))
data.stds<-data.frame(cbind(data.stds,as.character(dat[["location"]]) ))
colnames(data.stds)[16]<-"location"
#colnames(data.stds)[1]<-"bleach_intensity"

library(glmmTMB)
modelText<-lapply(vifPredCombinations_new, prepareModelText, data.stds)
modList<-lapply(modelText, evalTextModel)
modList<-modList[!sapply(modList, is.null)] 
#save.image("modList_global.stress.RData")

#modelSel.clim.exp<-model.sel(modList, rank.args = list(REML = FALSE), extra =c(AIC, BIC, Ovd = function(x) overdisp_fun(x)[[4]]))

library(MuMIn)

#modelSel.GSMsub<-model.sel(modList, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
#write.csv(modelSel.GSMsub, 'modelSel.GSM Subtraction.csv')

#modelSel.cummulativeDHM<-model.sel(modList, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
#write.csv(modelSel.cummulativeDHM, 'modelSel.cummulativeDHM.csv')

#modelSel.GSMdiv<-model.sel(modList, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
#write.csv(modelSel.GSMdiv, 'modelSel.GSM Division.csv')

#modelSel.CEsub<-model.sel(modList, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
#write.csv(modelSel.CEsub, 'modelSel.CE Subtraction.csv')

#modelSel.CEdiv<-model.sel(modList, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
#write.csv(modelSel.CEdiv, 'modelSel.CE Division.csv')
#modelSel.DHMdiv<-model.sel(modList, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
#write.csv(modelSel.DHMdiv, 'modelSel.DHM Division.csv')

modelSel.DHMsub<-model.sel(modList, rank.args = list(REML = FALSE), extra =c(AIC, BIC))
write.csv(modelSel.DHMsub, 'modelSel.DHM Subtraction.csv')
#sjPlot
library(sjPlot)
#tab_model(modList, file = "GSM Subtraction.doc")
#tab_model(modList, file = "GSM Division.doc")
#tab_model(modList, file = "CE Subtraction.doc")
#tab_model(modList, file = "CE Division.doc")
#tab_model(modList, file = "DHM Division.doc")
tab_model(modList, file = "DHM Subtraction.doc")
