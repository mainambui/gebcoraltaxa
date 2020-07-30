rm(list = ls())
require(dplyr)
require(MuMIn)
require(glmmTMB)
#setwd('/Users/josephmaina/Dropbox/Global Bleaching/Paper/Resistance/Data/')
setwd('~/Dropbox/Global Bleaching/Paper/Resistance/Data/')

#dat<-read.csv('/Users/josephmaina/Dropbox/Global Bleaching/Paper/Resistance/Data/Resistance_data_new.csv')

dat<-read.csv('~/Dropbox/Global Bleaching/Paper/Resistance/Data/Resistance_revised.csv')

# Use percent bleached division with location as a random intercept
respVar<-dat[, c(2, 6:8, 18:21)]

PredictVar<-dat[,c(9,10,11,12,13,15:17)]

data.std<-apply(X = PredictVar, MARGIN = 2,FUN = function(x){(x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))})

data.stds<-data.frame(cbind(respVar,data.std))
data <-data.frame(cbind(respVar,PredictVar, dat[,c("Region")]))
colnames(data)[17] <- "Region"
#convert all to numeric
data.stds<-data.frame(apply(data.stds, 2, function(x) as.numeric(x)))
#data.stds<-data.frame(cbind(data.stds,as.character(dat[["location"]]) ))
#colnames(data.stds)[16]<-"location"
data.stds <- left_join(data.stds, dat[,c("unique.id", "location", "Region")]) 

# Normalize new resistance metric to 0-1
normalize <- function(x){(x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))}

data.stds$Res_DHM_norm <- normalize(data.stds$Res_DHM_pbl_div)
data.stds$Res_CE_norm <- normalize(data.stds$Res_CE_pbl_div)

data.stds$Res_DHM_norm[data.stds$Res_DHM_norm == 0] <- 0.0001
data.stds$Res_DHM_norm[data.stds$Res_DHM_norm == 1] <- 0.999

data.stds$Res_CE_norm[data.stds$Res_CE_norm == 0] <- 0.0001
data.stds$Res_CE_norm[data.stds$Res_CE_norm == 1] <- 0.999

data.stds <- na.omit(data.stds)
 
 ### Plots
 library(reshape2)
 library(dplyr)
  data %>%
    melt(id.vars=c("Res_CE_pbl_div", "CTN"), measure.vars = c(8,10:15)) %>%
    ggplot(aes(y = Res_CE_pbl_div, x= value, col=CTN)) + geom_point()+facet_wrap(.~variable, scales = "free_x") + theme_bw() + theme(legend.position = "top") + ylab("Resistance (Climate exposure)") + xlab("")
  ggsave("Resistance CE Bivariate.pdf", width=9, height = 9)  
  
  data %>%
    melt(id.vars=c("Res_DHM_pbl_div", "CTN"), measure.vars = c(8,10:15)) %>%
    ggplot(aes(y = Res_DHM_pbl_div, x= value, col=CTN)) + geom_point()+facet_wrap(.~variable, scales = "free_x") + theme_bw() + theme(legend.position = "top") + ylab("Resistance (Cumulative DHM)") + xlab("")
ggsave("Resistance DHM Bivariate.pdf", width=9, height = 9)  

## fitdistrplus
install.packages("fitdistrplus")
library(fitdistrplus)
dati <- data.stds$Res_DHM_norm
fits <- list(
  no = fitdist(dati,"norm", discrete = FALSE),
  lo = fitdist(dati,"logis", discrete = FALSE),
  ca = fitdist(dati,"cauchy", discrete = FALSE),
  we = fitdist(dati, "weibull", discrete = FALSE),
  un = fitdist(dati, "unif", discrete = FALSE),
  be = fitdist(dati, "beta", discrete = FALSE),
  ln = fitdist(dati, "lnorm", discrete = FALSE),
  ga = fitdist(dati, "gamma", discrete = FALSE)
  )
sapply(fits, function(i) i$loglik)
sapply(fits, function(i) i$aic)

descdist(data.stds$Res_DHM_norm)
descdist(data.stds$Res_DHM_pbl_div)
qqcomp(fitdist(data.stds$Res_DHM_norm, "beta"))

a <- data.stds %>%
  ggplot(aes(x=Res_DHM_pbl_div)) + geom_histogram(aes(y=..density..),bins = 40, color="black", fill="white") + geom_density() + theme_bw() + xlab("Resistance (Cumulative DHM)")
b <- data.stds %>%
  ggplot(aes(x=Res_CE_pbl_div)) + geom_histogram(aes(y=..density..),bins = 40, color="black", fill="white") + geom_density() + theme_bw() + xlab("Resistance (Climate exposure)")
library(patchwork)
pdf(file = "Resistance distribution.pdf", height = 6, width = 4)
a/b
dev.off()
### Cumulative DHM
global.model.identity <- glmmTMB(formula = log(Res_DHM_pbl_div) ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + (1|location) , family=gaussian(link = "identity"),  dispformula = ~1, data = data.stds )

AIC(global.model.identity)
res.identity <- simulateResiduals(global.model.identity)
plot(res.identity)
res$scaledResiduals
pdf(file = "Gaussian link identity_Random intercept.pdf", width = 12, height = 5)
plot(res.identity)
testResiduals(res.identity)
testSpatialAutocorrelation(res.identity)
dev.off()

global.model.log <- glmmTMB(formula = Res_DHM_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + (1|location) , family=gaussian(link = "log"),  dispformula = ~1, data = data.stds )
glm
AIC(global.model.log)
res.log <- simulateResiduals(global.model.log)

pdf(file = "Gaussian link log_Random intercept.pdf", width = 12, height = 5)
plot(res.log)
testResiduals(res.log)
testSpatialAutocorrelation(res.log)
dev.off()

anova(global.model.identity, global.model.log)

### Non Random
global.model.identity.nr <- glmmTMB(formula = Res_DHM_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness , family=gaussian(link = "identity"),  dispformula = ~1, data = data.stds )

AIC(global.model.identity.nr)
res.identity.nr <- simulateResiduals(global.model.identity.nr)

pdf(file = "Gaussian link identity_non Random intercept.pdf", width = 12, height = 5)
plot(res.identity.nr)
testResiduals(res.identity.nr)
dev.off()

global.model.log.nr <- glmmTMB(formula = Res_DHM_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + Region , family=gaussian(link = "log"),  dispformula = ~1, data = data.stds )
plot(global.model.log.nr)

AIC(global.model.log.nr)
res.log.nr <- simulateResiduals(global.model.log.nr)

pdf(file = "Gaussian link log_non Random intercept.pdf", width = 12, height = 5)
plot(res.log.nr)
testResiduals(res.log.nr)
testSpatialAutocorrelation(res.log.nr)
dev.off()
#family=beta_family(link = "logit")
#install.packages("DHARMa")
#library(DHARMa)

### Gamma link log
global.model.Glog.nr <- glmmTMB(formula = Res_DHM_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + (1|location), family=Gamma(link = "log"),  dispformula = ~1, data = data.stds )


res.Glog.nr <- simulateResiduals(global.model.Glog.nr)

pdf(file = "Gamma link log_Random intercept.pdf", width = 12, height = 5)
plot(res.Glog.nr)
testResiduals(res.Glog.nr)
dev.off()
testDispersion(res)

### Beta 
global.model.beta <- glmmTMB(formula = Res_DHM_norm ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + (1|location) , family= beta_family, data = data.stds )

res.beta <- simulateResiduals(global.model.beta)
pdf(file = "Beta link logit_Random intercept.pdf", width = 12, height = 5)
plot(res.beta)
testResiduals(res.beta)
dev.off()

global.model.beta.nr <- glmmTMB(formula = Res_DHM_norm ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness  , family= beta_family, data = data.stds )

res.beta.nr <- simulateResiduals(global.model.beta.nr)
pdf(file = "Beta link logit_non Random intercept.pdf", width = 12, height = 5)
plot(res.beta.nr)
testResiduals(res.beta.nr)
dev.off()


#global.model <- lme4::lmer(formula = Res_DHM_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + (1|location), data = data.stds, na.action = "na.fail" )

# Final models
global.model.log.nr <- glm(formula = Res_DHM_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + Region , family=gaussian(link = "log"), data = data.stds )

testResiduals(global.model.log.nr)
results <- dredge(global.model.log.nr, beta = "none", rank = AICc )
 mymodel <- model.avg(results, subset = delta <2 )
mod1 <- glmmTMB(formula = Res_DHM_pbl_div ~ Longitude + Skewness  +
           (1 | location), data = data.stds, ziformula = ~0, dispformula = ~1)
mod2 <- glmmTMB(formula = Res_DHM_pbl_div ~ Longitude + Skewness + SST +
                  (1 | location), data = data.stds, ziformula = ~0, dispformula = ~1)
 
### Climate exposure
global.model2 <- glmmTMB(formula = Res_CE_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness , family = gaussian(link = "log"), data = data.stds )

global.model22 <- glmmTMB(formula = Res_CE_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness  , family = gaussian(link = "log"), data = data.stds )

global.model23 <- glmmTMB(formula = Res_CE_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness +  (1|Region) , family = gaussian(link = "log"), data = data.stds )

anova(global.model22, global.model23)
res2 <- simulateResiduals(global.model2)
plot(res2)
testDispersion(res2)
testResiduals(res2)

global.model2.beta <- glmmTMB(formula = Res_CE_norm ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + Region , family= beta_family, data = data.stds )

res2.beta <- simulateResiduals(global.model2.beta)

plot(res2.beta)
testDispersion(res2.beta)
testResiduals(res2.beta)

results2 <- dredge(global.model2, beta = "none", rank = AICc )

mymodel2 <- model.avg(results2, subset = delta <2 )
subset(results2, cumsum(results2$weight) <= .95) 
subset(results2, 1/8 < weight/max(results2$weight)) 

mod3 <- glmmTMB(formula = Res_CE_pbl_div ~ Absolute.latitude + Longitude + Skewness  +
                  (1 | location), data = data.stds, ziformula = ~0, dispformula = ~1)

mod4 <- glmmTMB(formula = Res_CE_pbl_div ~ Absolute.latitude + Longitude +
                  (1 | location), data = data.stds, ziformula = ~0, dispformula = ~1)

mod5 <- glmmTMB(formula = Res_CE_pbl_div ~ Absolute.latitude + Kurtosis + Longitude +
                  (1 | location), data = data.stds, ziformula = ~0, dispformula = ~1)

mod6 <- glmmTMB(formula = Res_CE_pbl_div ~ Absolute.latitude + Kurtosis + Longitude + Skewness + (1 | location), data = data.stds, ziformula = ~0, dispformula = ~1)

tab_model( mod3, mod4, mod5, mod6, title = c("Climate exposure"), dv.labels = c("Delta AICc = 0", "Delta AICc = 0.8", "Delta AICc = 1.2", "Delta AICc = 1.8"), file = "Climate exposure top models.doc" )

library(sjPlot)
 tab_model( mod1, mod2, title = c("Cumulative DHM"), dv.labels = c("Delta AICc = 0", "Delta AICc = 1.8"), file = "Cumulative DHM top models.docx" )
 model.avg(results, subset = delta <2 )
 
 # Global stress model
 global.model3 <- glmmTMB(formula = Res_GSM_pbl_div ~ Longitude + Absolute.latitude + Hard.coral.perc. + No..of.genera + SST + Kurtosis + Skewness + (1|location), data = data.stds )

results3 <- dredge(global.model3, beta = "none", rank = AICc )
mymodel3 <- model.avg(results3, subset = delta <2 )

plot_model(global.model.beta)
 