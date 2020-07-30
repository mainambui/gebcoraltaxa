# Adapted code for spatial glmm model
rm(list = ls())
library(DHARMa)

require(MuMIn)
require(glmmTMB)
#library(here)
library(ggplot2)
library(sjPlot)
require(parallel)
library(performance)
library(see)
library(randomForest)
#setwd('/Users/josephmaina/Dropbox/Global Bleaching/Paper/Resistance/Data/')
#setwd('~/Dropbox/Global Bleaching/Paper/Resistance/Data/')

#dat<-read.csv('/Users/josephmaina/Dropbox/Global Bleaching/Paper/Resistance/Data/Resistance_data_new.csv')

dat<-read.csv("C:/Users/Maxwell Azali/Dropbox/Global Bleaching/Paper/Resistance/Data/Resistance data final.csv")
colnames(dat) [1] <- "unique.id"
# Use percent bleached division with location as a random intercept
respVar<-dat[, c(1,10,11,19:21)]

PredictVar<-dat[,c(12:14,16:18)]

data.std<-apply(X = PredictVar, MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))})

data.stds<-data.frame(cbind(respVar,data.std,dat[,c("Region","CoralProvince","Ecoregion","Province")]))

with(data.stds,table(CoralProvince))

with(data.stds,table(Province))
unique(data.stds$Province)
data.stds$Province <- relevel(data.stds$Province, "Western Indian Ocean")

data.stds$pos <- numFactor(scale( data.stds$Longitude ), scale( data.stds$Latitude ))

data.stds$ID <- factor(rep(1, nrow(data.stds)))

data.complete <- na.omit(data.stds)
# Variogram correlation structure
library(gstat)

v_Res_CE <- variogram(Res_CE_pbl_div ~ 1, loc= ~Longitude+Latitude, data = data.complete)
plot(v_Res_CE)
best_vario <- fit.variogram(v_Res_CE, model = vgm(c("Exp", "Gau", "Mat", "Sph")))
best_vario
plot(v_Res_CE, model = best_vario)

v_Res_DHM <- variogram(Res_DHM_pbl_div ~ 1, loc= ~Longitude+Latitude, data = data.complete)
plot(v_Res_DHM)
best_vario1 <- fit.variogram(v_Res_DHM, model = vgm(c("Exp", "Gau", "Mat", "Sph")))
best_vario1
plot(v_Res_DHM, model = best_vario1)

# Parallel compute saves time  
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
clusterExport(clust, "data.complete")
invisible(clusterCall(clust, "library", "glmmTMB", character.only = TRUE))
# Specify global model 
# Note Gamma link log improves model fit by 20 AICc points

tmbMod1 <- glmmTMB(formula = Res_CE_pbl_div ~ Hard.coral.perc. + No..of.genera + SST + Skewness  + CoralProvince + Kurtosis + exp(pos + 0 | ID) , family=Gamma(link = "log"), data = data.complete )

tmbMod1$sdr$pdHess
summary(tmbMod1)
AICc(tmbMod1)

plot_model(tmbMod1)

#tmbMod1Vif <- check_collinearity(tmbMod1)
#plot(tmbMod1Vif)
#plot_model(tmbMod1, type = "eff", terms = c("CoralProvince"), sort.est = "re", show.intercept = TRUE)

#AICc(tmbMod1)




# Dredge but specify that models should not contain both SST and Kurtosis, due to high correlation
print(system.time(allmod<-pdredge(tmbMod1, cluster = clust, rank = AICc, subset = !(`cond(SST)` & `cond(Kurtosis)`))))


stopCluster(clust)
### Models delta AICc < 2
allmodsub <- subset(allmod, delta<2)

### Save dredged modes and best models <2 AICc
save(allmod, file="allmods.Rda")
save(allmodsub, file="allmodsub.Rda")


# Get top models delta AICc<2
m1 <-get.models(allmodsub, subset = 1)[[1]]
m2 <-get.models(allmodsub, subset = 2)[[1]]
m3 <-get.models(allmodsub, subset = 3)[[1]]
m4 <-get.models(allmodsub, subset = 4)[[1]]
m5 <-get.models(allmodsub, subset = 5)[[1]]

#check VIF of top models
check_collinearity(m2)
check_collinearity(m3)
check_collinearity(m4)
check_collinearity(m5)

# check model assumptions DHARMA

testResiduals( simulateResiduals(m1))
testResiduals( simulateResiduals(m2))
testResiduals( simulateResiduals(m3))
testResiduals( simulateResiduals(m4))
testResiduals( simulateResiduals(m5))

coefTable(m1,m2)

write.csv(allmod, file="Allmodels_ClimateExposure_Keith.csv")








#------------------------------------------------------
# Model with Spalding 2007 classification

#------------------------------------------------------
# Parallel dredge, saves time  
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))
clusterExport(clust, "data.complete")
invisible(clusterCall(clust, "library", "glmmTMB", character.only = TRUE))


# Specify global model 

# Note Gamma link log improves model fit by 20 AICc points

tmbMod2 <- glmmTMB(formula = Res_CE_pbl_div ~ Hard.coral.perc. + No..of.genera + SST + Skewness  + Province + Kurtosis + exp(pos + 0 | ID) , family=Gamma(link = "log"), data = data.complete)

tmbMod2$sdr$pdHess
summary(tmbMod2)
AICc(tmbMod2)

plot_model(tmbMod2)

#tmbMod2Vif <- check_collinearity(tmbMod2)
#plot(tmbMod1Vif)
#plot_model(tmbMod2, type = "eff", terms = c("CoralProvince"), sort.est = "re", show.intercept = TRUE)

#AICc(tmbMod2)






# Dredge but specify that models should not contain both SST and Kurtosis, due to high correlation
print(system.time(allmod2<-pdredge(tmbMod2, cluster = clust, rank = AICc, subset = !(`cond(SST)` & `cond(Kurtosis)`))))


stopCluster(clust)

### Models delta AICc < 2
allmodsub2 <- subset(allmod2, delta<2)

### Save dredged modes and best models <2 AICc
save(allmod2, file="allmods2.Rda")
save(allmodsub2, file="allmodsub2.Rda")


# Get top models delta AICc<2
m12 <-get.models(allmodsub2, subset = 1)[[1]]
update()

#check VIF of top models
check_collinearity(m12)
check_collinearity(m3)
check_collinearity(m4)
check_collinearity(m5)

# check model assumptions DHARMA

testResiduals( simulateResiduals(m12))


write.csv(allmod2, file="Allmodels_ClimateExposure_Spalding.csv")





#------------------------------------------------------
# DHM Models with Keith 2013 classification

#------------------------------------------------------
# Parallel dredge, saves time  
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))
clusterExport(clust, "data.complete")
invisible(clusterCall(clust, "library", "glmmTMB", character.only = TRUE))


# Specify global model 

# Note Gamma link log improves model fit by 20 AICc points

tmbMod3 <- glmmTMB(formula = Res_DHM_pbl_div ~ Hard.coral.perc. + No..of.genera + SST + Skewness  + CoralProvince + Kurtosis + exp(pos + 0 | ID) , family=Gamma(link = "log"), data = data.complete)

tmbMod3$sdr$pdHess
summary(tmbMod3)
AICc(tmbMod3)

plot_model(tmbMod3)

#tmbMod2Vif <- check_collinearity(tmbMod2)
#plot(tmbMod1Vif)
#plot_model(tmbMod2, type = "eff", terms = c("CoralProvince"), sort.est = "re", show.intercept = TRUE)

#AICc(tmbMod2)






# Dredge but specify that models should not contain both SST and Kurtosis, due to high correlation
print(system.time(allmod3<-pdredge(tmbMod3, cluster = clust, rank = AICc, subset = !(`cond(SST)` & `cond(Kurtosis)`))))


stopCluster(clust)

### Models delta AICc < 2
allmodsub3 <- subset(allmod3, delta<2)

### Save dredged modes and best models <2 AICc
save(allmod3, file="allmods3.Rda")
save(allmodsub3, file="allmodsub3.Rda")


# Get top models delta AICc<2
m13 <-get.models(allmodsub3, subset = 1)[[1]]
m23 <-get.models(allmodsub3, subset = 2)[[1]]
m33 <-get.models(allmodsub3, subset = 3)[[1]]

#check VIF of top models
check_collinearity(m13)
check_collinearity(m23)
check_collinearity(m33)


# check model assumptions DHARMA

testResiduals( simulateResiduals(m13))
testResiduals( simulateResiduals(m23))
testResiduals( simulateResiduals(m33))


write.csv(allmod3, file="Allmodels_DHM_Keith.csv")










#------------------------------------------------------
# DHM Models with Spalding 2007 classification

#------------------------------------------------------
# Parallel dredge, saves time  
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))
clusterExport(clust, "data.complete")
invisible(clusterCall(clust, "library", "glmmTMB", character.only = TRUE))


# Specify global model 

# Note Gamma link log improves model fit by 20 AICc points

tmbMod4 <- glmmTMB(formula = Res_DHM_pbl_div ~ Hard.coral.perc. + No..of.genera + SST + Skewness  + Province + Kurtosis + exp(pos + 0 | ID) , family=Gamma(link = "log"), data = data.complete)

tmbMod4$sdr$pdHess
summary(tmbMod4)
AICc(tmbMod4)

plot_model(tmbMod4)

#tmbMod2Vif <- check_collinearity(tmbMod2)
#plot(tmbMod1Vif)
#plot_model(tmbMod2, type = "eff", terms = c("CoralProvince"), sort.est = "re", show.intercept = TRUE)

#AICc(tmbMod2)






# Dredge but specify that models should not contain both SST and Kurtosis, due to high correlation
print(system.time(allmod4<-pdredge(tmbMod4, cluster = clust, rank = AICc, subset = !(`cond(SST)` & `cond(Kurtosis)`))))


stopCluster(clust)

### Models delta AICc < 2
allmodsub4 <- subset(allmod4, delta<2)

### Save dredged modes and best models <2 AICc
save(allmod4, file="allmods4.Rda")
save(allmodsub4, file="allmodsub4.Rda")


# Get top models delta AICc<2
m14 <-get.models(allmodsub4, subset = 1)[[1]]
m24 <-get.models(allmodsub4, subset = 2)[[1]]

#check VIF of top models
check_collinearity(m14)
check_collinearity(m24)


# check model assumptions DHARMA

testResiduals( simulateResiduals(m14))
testResiduals( simulateResiduals(m24))


write.csv(allmod4, file="Allmodels_DHM_Spalding.csv")



### Unused code
## Parallel computing tabmodel
num_cores <- (detectCores())		# All cores
num_cores

# Calculate the number of cores to use
num_cores <- (detectCores() - 1)	
num_cores

# Initiate cluster
cl <- makeCluster(5)	

# Load packages on the cluster named cl
clusterEvalQ(cl, library(glmmTMB)) 
clusterEvalQ(cl, library(MuMIn))
#clusterEvalQ(cl, library(sjPlot))
# Load the dataset named d on the cluster
clusterExport(cl, "data.stds")
clusterExport(cl, "allmodsub")

#tab_model(m1,m2,m3,m4,m5, show.aic = TRUE, show.ci = FALSE, show.se = TRUE,title = c("Resistance ClimateExposure KeithCoralProvince"), dv.labels = c("Delta AICc = 0", "Delta AICc = 0.20", "Delta AICc = 1.13","Delta AICc = 1.23","Delta AICc = 1.58"), file = "CE_Keith_topset.doc")

tab_model(m1, m2, m3, m4, m5, show.se = TRUE, show.aic = TRUE, show.aicc = TRUE, show.ci = FALSE,show.r2 = TRUE, collapse.ci = FALSE,  show.std = FALSE, show.stat = TRUE, show.icc = FALSE,string.se = "SE", string.p = "P", digits = 3, dv.labels = c("Model 1", "Model 2","Model 3","Model 4", "Model5"),col.order = c("est","se", "ci", "p"), file = "CE_Keith_topset.doc")

stargazer::stargazer(m1)
huxtable::huxreg(m1)
stopCluster(cl)
tab_model(m1)
#allmod<-lapply(pdredge(tmbMod1, cluster = clust, rank = AICc, subset = !(`cond(SST)` & `cond(Kurtosis)`)),eval)
