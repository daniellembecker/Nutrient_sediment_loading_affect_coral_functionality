##Using purr to pull out effect sizes for data
##created by Danielle Becker 10/08/19
##edited by Danielle Becker 11/08/19

#clear dataframe#####
rm(list=ls())

##Install packages
# load packages
library(tidyverse)
library(png)
library(grid)
library(effsize)
library(parameters)
library(psycho)
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(lmerTest)
library(here)

#set wd
here()

#read in covariate data sheet
covariate.data <- read.csv("Effect_Sizes/Data/covariate.data.csv")

#standardize data to z scores for comparisons
covariate.data$Topt_scale <- scale(covariate.data$Topt, center = TRUE, scale = TRUE)
covariate.data$Pmax_scale<- scale(covariate.data$Pmax, center = TRUE, scale = TRUE)
covariate.data$lnc_scale <- scale(covariate.data$lnc, center = TRUE, scale = TRUE)
covariate.data$zoox_scale <- scale(covariate.data$zoox.per.cm2, center = TRUE, scale = TRUE)
covariate.data$chl_scale <- scale(covariate.data$chlA.ugcm2, center = TRUE, scale = TRUE)
covariate.data$afdw_scale <- scale(covariate.data$AFDW.mg.cm2., center = TRUE, scale = TRUE)
covariate.data$N.AT_scale <- scale(covariate.data$per.N.AT, center = TRUE, scale = TRUE)
covariate.data$N.ST_scale <- scale(covariate.data$N.ST, center = TRUE, scale = TRUE)
covariate.data$per.N.ST_scale <- scale(covariate.data$per.N.ST, center = TRUE, scale = TRUE)

#############################################################################################
#NUTRIENT GRADIENT
#make data sheet that pulls out the confidence intervals and effect sizes for extrinsic values
#run models for each parameter and then save the estimate and confidence interval summary table
Topt_scale.GP <- lmer(Topt_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
#a <- summary(Topt_scale.GP)$coefficients[2,1:2]
#confint(Topt_scale.GP, parm = "N.resid")
a <-model_parameters(Topt_scale.GP)
anova(Topt_scale.GP)

Topt_scale.R <- lmer(Topt_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="R")
#b <- summary(Topt_scale.R)$coefficients[2,1:2]
#confint(Topt_scale.R, parm = "N.resid")
b <- model_parameters(Topt_scale.R)
anova(Topt_scale.R)

Topt_scale.C <- lmer(Topt_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="C")
#c <- summary(Topt_scale.C)$coefficients[2,1:2]
#confint(Topt_scale.C, parm = "N.resid")
c <- model_parameters(Topt_scale.C)
anova(Topt_scale.C)

Pmax_scale.GP <- lmer(Pmax_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
#d <- summary(Pmax_scale.GP)$coefficients[2,1:2]
##confint(Pmax_scale.GP, parm = "N.resid", level = 0.90)
d <- model_parameters(Pmax_scale.GP)
anova(Pmax_scale.GP)

Pmax_scale.R <- lmer(Pmax_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="R")
#e <- summary(Pmax_scale.R)$coefficients[2,1:2]
#confint(Pmax_scale.R, parm = "N.resid", level = 0.90)
e <- model_parameters(Pmax_scale.R)
anova(Pmax_scale.R)

Pmax_scale.C <- lmer(Pmax_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="C")
#f <- summary(Pmax_scale.C)$coefficients[2,1:2]
#confint(Pmax_scale.C, parm = "N.resid", level = 0.90)
f <- model_parameters(Pmax_scale.C)
anova(Pmax_scale.C)

lnc_scale.GP <- lmer(lnc_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
#g <- summary(lnc_scale.GP)$coefficients[2,1:2]
#confint(lnc_scale.GP, parm = "N.resid", level = 0.90)
g <- model_parameters(lnc_scale.GP)
anova(lnc_scale.GP)

lnc_scale.R <- lmer(lnc_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="R")
#h <- summary(lnc_scale.R)$coefficients[2,1:2]
#confint(lnc_scale.R, parm = "N.resid", level = 0.90)
h <- model_parameters(lnc_scale.R)
anova(lnc_scale.R)

lnc_scale.C <- lmer(lnc_scale~N.resid + (1|site), data = covariate.data, subset = rate.type=="C")
#i <- summary(lnc_scale.C)$coefficients[2,1:2]
#confint(lnc_scale.C, parm = "N.resid", level = 0.90)
i <- model_parameters(lnc_scale.C)
anova(lnc_scale.C)

#rbind all sumary data for each param into one data sheet
extrin.dat <- rbind(a,b,c,d,e,f,g,h,i)

#make data table into a data frame
extrin.dat <- as.data.frame(extrin.dat)

#delete the rows with intercept information
extrin.dat<- extrin.dat[-c(1, 3, 5,7,9,11,13,15,17), ] 

#make rate type column with appropriate labels
extrin.dat$rate.type  <- c("GP", "R", "C", "GP", "R", "C","GP", "R", "C")

#make extrinsic param column with appropriate labels
extrin.dat$extrinsic.param  <- c("Topt_scale", "Topt_scale", "Topt_scale", "Pmax_scale", "Pmax_scale", "Pmax_scale", "b(Tc)", "b(Tc)", "b(Tc)")

#make environmental parma column with appropriate labels
extrin.dat$environ.param  <- c("N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid")

#make column with significant or not
extrin.dat$sig  <- c("non", "sig", "non", "sig", "non", "non", "sig", "non", "sig")

#make new data column combining labels for rate type and extrinsic param type
extrin.dat$combined <- paste(extrin.dat$extrinsic.param, extrin.dat$rate.type)

#rename coefficient column to estimat
colnames(extrin.dat)[2] <- "Estimate"

#change names of rate type params in data sheet
extrin.dat <- extrin.dat %>% 
  mutate(rate.type = recode(rate.type, `C` = "Net Calcification", `GP` = "Gross Photosynthesis", `R` = "Dark Respiration"))

#reorder data frame levels
extrin.dat$combined <- factor(extrin.dat$combined, c("b(Tc) C", "b(Tc) R", "b(Tc) GP", "Pmax_scale C", "Pmax_scale R", "Pmax_scale GP", "Topt_scale C", "Topt_scale R", "Topt_scale GP"))

 #reorder the legend labels
extrin.dat$rate.type <- factor(extrin.dat$rate.type, levels = c("Gross Photosynthesis", "Dark Respiration", "Net Calcification"))

#make plot of standardized effect sizes 
ggplot(extrin.dat, aes(x = Estimate, y = combined, col = rate.type, alpha = sig)) + #set aestetics for the graph
  theme_bw() + #chnage background theme
  geom_hline(yintercept = 9.35, lty = 1)  +
  geom_hline(yintercept = 6.5, lty = 1)  +
  geom_hline(yintercept = 3.5, lty = 1)  +
  theme(axis.title.x = element_text(size = 14, color = "black"), axis.text.x=element_text(size = 11, color = "black"), axis.text.y = element_text(size = 11, color = "black")) +
  geom_point(size = 3) + #increase the point size 
  geom_errorbarh(aes(x=Estimate, xmin=CI_low, xmax=CI_high), height = 0) + #graph estimates and st err for all values
  geom_vline(xintercept = 0, lty = 2)  + #add vertical line to graph and make it dashed 
  scale_y_discrete(breaks=c("Topt_scale R", "Topt_scale GP", "Topt_scale C", "Pmax_scale R", "Pmax_scale GP", "Pmax_scale C", "b(Tc) R", "b(Tc) GP", "b(Tc) C"), labels=c("R", "GP", "C", "R", "GP", "C", "R", "GP", "C")) +  #rename y axis tickmarks 
  scale_colour_manual(values = c('green4', 'blue', 'black')) +
  guides(alpha = FALSE) + #remove legend for alpha
  scale_alpha_discrete(range = c(0.3, 1.0)) + #set range for sig vs non sig alphas
  labs(col = "Rate Type") +
  xlab("Standardized Effect Sizes") + #rename x label
  ylab("")  #remove y label

ggsave(filename = "Effect_Sizes/Output/holobiont.effect.sizes.png", device = "png", width = 6, height = 5)


##############3#make data sheet that pulls out the confidence intervals and effect sizes######## 
#run models for each parameter and then save the estimate and CI summary table
sym <- lmer(zoox_scale~N.resid + (1|site), data = covariate.data)
#j <- summary(sym)$coefficients[2,1:2]
j <-model_parameters(sym)

chloro <- lmer(chl_scale~N.resid + (1|site), data = covariate.data)
#k <- summary(chloro)$coefficients[2,1:2]
k <-model_parameters(chloro)

AFDW <- lmer(afdw_scale~N.resid + (1|site), data = covariate.data)
#l <- summary(AFDW)$coefficients[2,1:2]
l <-model_parameters(AFDW)

N.AT <- lmer(N.AT_scale~N.resid + (1|site), data = covariate.data)
#m <- summary(N.AT)$coefficients[2,1:2]
m <-model_parameters(N.AT)

N.ST <- lmer(N.ST_scale~N.resid + (1|site), data = covariate.data)
#n <- summary(N.ST)$coefficients[2,1:2]
n <- model_parameters(N.ST)

per.N.ST <- lmer(per.N.ST_scale~N.resid + (1|site), data = covariate.data)
#n <- summary(N.ST)$coefficients[2,1:2]
o <- model_parameters(per.N.ST)

#rbind all sumary data for each param into one data sheet
intrin.dat <- rbind(j,k,l,m,n,o)

#make data table into a data frame
intrin.dat <- as.data.frame(intrin.dat)

#delete the rows with intercept information
intrin.dat<- intrin.dat[-c(1, 3, 5,7,9,11), ] 

#make rate type column with appropriate labels
intrin.dat$intrin.type  <- c("Endosymbiont Density", "Chlorophyll a Content", "Tissue Biomass", "% N Content Animal Tissue", "N Content Endosymbiont Tissue per cell", "% N content endosymbionts" )

#make environmental parma column with appropriate labels
intrin.dat$environ.param  <- c("N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid")

#make column with significant or not
intrin.dat$sig  <- c("sig", "sig", "non", "non", "sig", "sig")

#rename standard error column
colnames(intrin.dat)[2] <- "Estimate"

#reorder y axis labels
intrin.dat$Variable <- ordered(intrin.dat$intrin.type, levels = c("Endosymbiont Density", "Chlorophyll a Content", "N Content Endosymbiont Tissue per cell", "% N content endosymbionts", "% N Content Animal Tissue", "Tissue Biomass"))

#create a ggplot with all the effect sizes
ggplot(intrin.dat, aes(x = Estimate, y = Variable, alpha = sig)) + #set aestetics for the graph
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_text(size = 14), axis.text.x=element_text(size = 11), axis.text.y = element_text(size = 11)) +
  geom_point(size = 3) + #increase the point size 
  geom_errorbarh(aes(x=Estimate, xmin=CI_low, xmax=CI_high), height = 0) + #graph estimates and st err for all values
  geom_vline(xintercept = 0, lty = 2)  + #add vertical line to graph and make it dashed 
  scale_y_discrete(breaks=c("Tissue Biomass", "Endosymbiont Density", "Chlorophyll a Content",  "% N Content Animal Tissue", "N Content Endosymbiont Tissue per cell", "% N content endosymbionts"), labels=c("tissue biomass", "endosymbiont density", "chlorophyll a content", "coral % N content", " N content ST per cell", "% N symbionts")) +  #rename y axis tickmarks
  #scale_colour_manual(values = c('black', 'green4', 'blue', 'purple', 'darkred', 'red')) + #assign colors to legend
  guides(alpha = FALSE) + #remove legend for alpha
  scale_alpha_discrete(range = c(0.3, 1.0)) + #set range for sig vs non sig alphas
  xlab("Standardized Effect Sizes") + #rename x label
  ylab("")  #remove y label
  

ggsave(filename = "Effect_Sizes/Output/endo.coral.effect.sizes.png", device = "png", width = 5, height = 5)

######################################################################################
#running lmer models to compare influence of nutrient gradient on all params
#make data sheets with just GP, C, and R

Topt.GP <- lmer(Topt~N.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
anova(Topt.GP)
summary(Topt.GP)
tab_model(Topt.GP)

Topt.R <- lmer(Topt~N.resid + (1|site), data = covariate.data, subset = rate.type=="R")
anova(Topt.R)
summary(Topt.R)
tab_model(Topt.R)

Topt.C <- lmer(Topt~N.resid + (1|site), data = covariate.data, subset = rate.type=="C")
anova(Topt.C)
summary(Topt.C)
tab_model(Topt.C)

Pmax.GP <- lmer(Pmax~N.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
anova(Pmax.GP)
summary(Pmax.GP)
tab_model(Pmax.GP)

Pmax.R <- lmer(Pmax~N.resid + (1|site), data = covariate.data, subset = rate.type=="R")
anova(Pmax.R)
summary(Pmax.R)
tab_model(Pmax.R)

Pmax.C <- lmer(Pmax~N.resid + (1|site), data = covariate.data, subset = rate.type=="C")
anova(Pmax.C)
summary(Pmax.C)
tab_model(Pmax.C)

lnc.GP <- lmer(lnc~N.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
anova(lnc.GP)
summary(lnc.GP)
tab_model(lnc.GP)

lnc.R <- lmer(lnc~N.resid + (1|site), data = covariate.data, subset = rate.type=="R")
anova(lnc.R)
summary(lnc.R)
tab_model(lnc.R)

lnc.C <- lmer(lnc~N.resid + (1|site), data = covariate.data, subset = rate.type=="C")
anova(lnc.C)
summary(lnc.C)
tab_model(lnc.C)

sym <- lmer(zoox.per.cm2~N.resid + (1|site), data = covariate.data)
anova(sym)
tab_model(sym)

chloro <- lmer(chlA.ugcm2~N.resid + (1|site), data = covariate.data)
anova(chloro)
tab_model(chloro)

AFDW <- lmer(AFDW.mg.cm2.~N.resid + (1|site), data = covariate.data)
anova(AFDW)
tab_model(AFDW)

N.AT <- lmer(per.N.AT~N.resid + (1|site), data = covariate.data)
anova(N.AT)
tab_model(N.AT)

N.ST <- lmer(N.ST~N.resid + (1|site), data = covariate.data)
anova(N.ST)
tab_model(N.ST)

per.N.ST <- lmer(per.N.ST~N.resid + (1|site), data = covariate.data)
anova(per.N.ST)
tab_model(per.N.ST)

#organize the results for all models into distinct tables for nutrient comparisons
Topt_models <- tab_model(Topt.GP, Topt.R, Topt.C, pred.labels = c("Intercept", "% Nitrogen Content Residuals"),
                         dv.labels = c("Topt_GP", "Topt_R", "Topt_C"), string.ci = "Conf. Int (95%)",
                         string.p = "P-Value",file = "../../../Documents/CSUN/Thesis Defense/Tables/Topt_%N.doc")

Pmax_models <- tab_model(Pmax.GP, Pmax.R, Pmax.C, pred.labels = c("Intercept", "% Nitrogen Content Residuals"),
                         dv.labels = c("umax_GP", "umax_R", "umax_C"), string.ci = "Conf. Int (95%)",
                         string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/Pmax_%N.doc")

lnc_models <- tab_model(lnc.GP, lnc.R, lnc.C, pred.labels = c("Intercept",  "% Nitrogen Content Residuals"),
                        dv.labels = c("b(Tc)_GP", "b(Tc)_R", "b(Tc)_C"), string.ci = "Conf. Int (95%)",
                        string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/lnc_%N.doc")

intrinsic_models <- tab_model(sym, chloro, AFDW, N.AT, N.ST, per.N.ST,  collapse.ci = TRUE, pred.labels = c("Intercept",  "% Nitrogen Content Residuals"),
                              dv.labels = c("Endosymbiont Density", "Chlorophyll a Content", "Tissue Biomass", "Coral % N Content", "Endosymbiont N Content Per Cell", "Endosymbiont % N Content"), string.ci = "Conf. Int (95%)",
                              string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/intrinsic_%N.doc")


#############################################################################################
#SEDIMENTATION GRADIENT
#make data sheet that pulls out the confidence intervals and effect sizes for extrinsic values
#run models for each parameter and then save the estimate and confidence interval summary table
Topt_scale.GP.trap <- lmer(Topt_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
#a <- summary(Topt_scale.GP)$coefficients[2,1:2]
#confint(Topt_scale.GP, parm = "N.resid")
a <-model_parameters(Topt_scale.GP.trap)
anova(Topt_scale.GP.trap)

Topt_scale.R.trap <- lmer(Topt_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="R")
#b <- summary(Topt_scale.R)$coefficients[2,1:2]
#confint(Topt_scale.R, parm = "N.resid")
b <- model_parameters(Topt_scale.R.trap)
anova(Topt_scale.R.trap)

Topt_scale.C.trap <- lmer(Topt_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="C")
#c <- summary(Topt_scale.C)$coefficients[2,1:2]
#confint(Topt_scale.C, parm = "N.resid")
c <- model_parameters(Topt_scale.C.trap)
anova(Topt_scale.C.trap)

Pmax_scale.GP.trap <- lmer(Pmax_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
#d <- summary(Pmax_scale.GP)$coefficients[2,1:2]
##confint(Pmax_scale.GP, parm = "N.resid", level = 0.90)
d <- model_parameters(Pmax_scale.GP.trap)
anova(Pmax_scale.GP.trap)

Pmax_scale.R.trap <- lmer(Pmax_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="R")
#e <- summary(Pmax_scale.R)$coefficients[2,1:2]
#confint(Pmax_scale.R, parm = "N.resid", level = 0.90)
e <- model_parameters(Pmax_scale.R.trap)
anova(Pmax_scale.R.trap)

Pmax_scale.C.trap <- lmer(Pmax_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="C")
#f <- summary(Pmax_scale.C)$coefficients[2,1:2]
#confint(Pmax_scale.C, parm = "N.resid", level = 0.90)
f <- model_parameters(Pmax_scale.C.trap)
anova(Pmax_scale.C.trap)

lnc_scale.GP.trap <- lmer(lnc_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
#g <- summary(lnc_scale.GP)$coefficients[2,1:2]
#confint(lnc_scale.GP, parm = "N.resid", level = 0.90)
g <- model_parameters(lnc_scale.GP.trap)
anova(lnc_scale.GP.trap)

lnc_scale.R.trap <- lmer(lnc_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="R")
#h <- summary(lnc_scale.R)$coefficients[2,1:2]
#confint(lnc_scale.R, parm = "N.resid", level = 0.90)
h <- model_parameters(lnc_scale.R.trap)
anova(lnc_scale.R.trap)

lnc_scale.C.trap <- lmer(lnc_scale~trap.resid + (1|site), data = covariate.data, subset = rate.type=="C")
#i <- summary(lnc_scale.C)$coefficients[2,1:2]
#confint(lnc_scale.C, parm = "N.resid", level = 0.90)
i <- model_parameters(lnc_scale.C.trap)
anova(lnc_scale.C.trap)

#rbind all sumary data for each param into one data sheet
extrin.dat <- rbind(a,b,c,d,e,f,g,h,i)

#make data table into a data frame
extrin.dat <- as.data.frame(extrin.dat)

#delete the rows with intercept information
extrin.dat<- extrin.dat[-c(1, 3, 5,7,9,11,13,15,17), ] 

#make rate type column with appropriate labels
extrin.dat$rate.type  <- c("GP", "R", "C", "GP", "R", "C","GP", "R", "C")

#make extrinsic param column with appropriate labels
extrin.dat$extrinsic.param  <- c("Topt_scale", "Topt_scale", "Topt_scale", "Pmax_scale", "Pmax_scale", "Pmax_scale", "b(Tc)", "b(Tc)", "b(Tc)")

#make environmental parma column with appropriate labels
extrin.dat$environ.param  <- c("N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid")

#make column with significant or not
extrin.dat$sig  <- c("non", "non", "non", "sig", "non", "non", "non", "non", "non")

#make new data column combining labels for rate type and extrinsic param type
extrin.dat$combined <- paste(extrin.dat$extrinsic.param, extrin.dat$rate.type)

#rename coefficient column to estimat
colnames(extrin.dat)[2] <- "Estimate"

#change names of rate type params in data sheet
extrin.dat <- extrin.dat %>% 
  mutate(rate.type = recode(rate.type, `C` = "Net Calcification", `GP` = "Gross Photosynthesis", `R` = "Dark Respiration"))

#reorder data frame levels
extrin.dat$combined <- factor(extrin.dat$combined, c("b(Tc) C", "b(Tc) R", "b(Tc) GP", "Pmax_scale C", "Pmax_scale R", "Pmax_scale GP", "Topt_scale C", "Topt_scale R", "Topt_scale GP"))

#reorder the legend labels
extrin.dat$rate.type <- factor(extrin.dat$rate.type, levels = c("Gross Photosynthesis", "Dark Respiration", "Net Calcification"))

#make plot of standardized effect sizes 
ggplot(extrin.dat, aes(x = Estimate, y = combined, col = rate.type, alpha = sig)) + #set aestetics for the graph
  theme_bw() + #chnage background theme
  geom_hline(yintercept = 9.35, lty = 1)  +
  geom_hline(yintercept = 6.5, lty = 1)  +
  geom_hline(yintercept = 3.5, lty = 1)  +
  theme(axis.title.x = element_text(size = 14, color = "black"), axis.text.x=element_text(size = 11, color = "black"), axis.text.y = element_text(size = 11, color = "black")) +
  geom_point(size = 3) + #increase the point size 
  geom_errorbarh(aes(x=Estimate, xmin=CI_low, xmax=CI_high), height = 0) + #graph estimates and st err for all values
  geom_vline(xintercept = 0, lty = 2)  + #add vertical line to graph and make it dashed 
  scale_y_discrete(breaks=c("Topt_scale R", "Topt_scale GP", "Topt_scale C", "Pmax_scale R", "Pmax_scale GP", "Pmax_scale C", "b(Tc) R", "b(Tc) GP", "b(Tc) C"), labels=c("R", "GP", "C", "R", "GP", "C", "R", "GP", "C")) +  #rename y axis tickmarks 
  scale_colour_manual(values = c('green4', 'blue', 'black')) +
  guides(alpha = FALSE) + #remove legend for alpha
  scale_alpha_discrete(range = c(0.3, 1.0)) + #set range for sig vs non sig alphas
  labs(col = "Rate Type") +
  xlab("Standardized Effect Sizes") + #rename x label
  ylab("")  #remove y label


ggsave(filename = "Effect_Sizes/Output/holobiont.effect.sizes.trap.png", device = "png", width = 6, height = 5)


##############3#make data sheet that pulls out the confidence intervals and effect sizes######## 
#run models for each parameter and then save the estimate and CI summary table
sym.trap <- lmer(zoox_scale~trap.resid + (1|site), data = covariate.data)
#j <- summary(sym)$coefficients[2,1:2]
j <-model_parameters(sym.trap)

chloro.trap <- lmer(chl_scale~trap.resid + (1|site), data = covariate.data)
#k <- summary(chloro)$coefficients[2,1:2]
k <-model_parameters(chloro.trap)

AFDW.trap <- lmer(afdw_scale~trap.resid + (1|site), data = covariate.data)
#l <- summary(AFDW)$coefficients[2,1:2]
l <-model_parameters(AFDW.trap)

N.AT.trap <- lmer(N.AT_scale~trap.resid + (1|site), data = covariate.data)
#m <- summary(N.AT)$coefficients[2,1:2]
m <-model_parameters(N.AT.trap)

N.ST.trap <- lmer(N.ST_scale~trap.resid + (1|site), data = covariate.data)
#n <- summary(N.ST)$coefficients[2,1:2]
n <- model_parameters(N.ST.trap)

per.N.ST <- lmer(per.N.ST_scale~trap.resid + (1|site), data = covariate.data)
#n <- summary(N.ST)$coefficients[2,1:2]
o <- model_parameters(per.N.ST)

#rbind all sumary data for each param into one data sheet
intrin.dat <- rbind(j,k,l,m,n,o)

#make data table into a data frame
intrin.dat <- as.data.frame(intrin.dat)

#delete the rows with intercept information
intrin.dat<- intrin.dat[-c(1, 3, 5,7,9,11), ] 

#make rate type column with appropriate labels
intrin.dat$intrin.type  <- c("Endosymbiont Density", "Chlorophyll a Content", "Tissue Biomass", "% N Content Animal Tissue", "Tissue N Content Endosymbiont Tissue per cell", "%N endosymbionts")

#make environmental parma column with appropriate labels
intrin.dat$environ.param  <- c("N.resid", "N.resid", "N.resid", "N.resid", "N.resid", "N.resid")

#make column with significant or not
intrin.dat$sig  <- c("sig", "sig", "non", "sig", "sig", "sig")

#rename standard error column
colnames(intrin.dat)[2] <- "Estimate"

#reorder y axis labels
intrin.dat$Variable <- ordered(intrin.dat$intrin.type, levels = c("Endosymbiont Density", "Chlorophyll a Content", "Tissue N Content Endosymbiont Tissue per cell", "%N endosymbionts",  "% N Content Animal Tissue", "Tissue Biomass"))

#create a ggplot with all the effect sizes
ggplot(intrin.dat, aes(x = Estimate, y = Variable, alpha = sig)) + #set aestetics for the graph
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_text(size = 14), axis.text.x=element_text(size = 11), axis.text.y = element_text(size = 11)) +
  geom_point(size = 3) + #increase the point size 
  geom_errorbarh(aes(x=Estimate, xmin=CI_low, xmax=CI_high), height = 0) + #graph estimates and st err for all values
  geom_vline(xintercept = 0, lty = 2)  + #add vertical line to graph and make it dashed 
  scale_y_discrete(breaks=c("Tissue Biomass", "Endosymbiont Density", "Chlorophyll a Content",  "% N Content Animal Tissue", "Tissue N Content Endosymbiont Tissue per cell", "%N endosymbionts"), labels=c("tissue biomass", "endosymbiont density", "chlorophyll a content", "coral % N content", " N content ST per cell", "% N symbionts")) +  #rename y axis tickmarks
  #scale_colour_manual(values = c('black', 'green4', 'blue', 'purple', 'darkred', 'red')) + #assign colors to legend
  guides(alpha = FALSE) + #remove legend for alpha
  scale_alpha_discrete(range = c(0.3, 1.0)) + #set range for sig vs non sig alphas
  xlab("Standardized Effect Sizes") + #rename x label
  ylab("")  #remove y label


ggsave(filename = "Effect_Sizes/Output/endo.coral.effect.sizes.trap.png", device = "png", width = 5, height = 5)

######################################################################################
#running lmer models to compare influence of sedimentation gradient on all params
#make data sheets with just GP, C, and R

Topt.GP.trap <- lmer(Topt~trap.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
anova(Topt.GP.trap)
summary(Topt.GP.trap)
tab_model(Topt.GP.trap)

Topt.R.trap <- lmer(Topt~trap.resid + (1|site), data = covariate.data, subset = rate.type=="R")
anova(Topt.R.trap)
summary(Topt.R.trap)
tab_model(Topt.R.trap)

Topt.C.trap <- lmer(Topt~trap.resid + (1|site), data = covariate.data, subset = rate.type=="C")
anova(Topt.C.trap)
summary(Topt.C.trap)
tab_model(Topt.C.trap)

Pmax.GP.trap <- lmer(Pmax~trap.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
anova(Pmax.GP.trap)
summary(Pmax.GP.trap)
tab_model(Pmax.GP.trap)

Pmax.R.trap <- lmer(Pmax~trap.resid + (1|site), data = covariate.data, subset = rate.type=="R")
anova(Pmax.R.trap)
summary(Pmax.R.trap)
tab_model(Pmax.R.trap)

Pmax.C.trap <- lmer(Pmax~trap.resid + (1|site), data = covariate.data, subset = rate.type=="C")
anova(Pmax.C.trap)
summary(Pmax.C.trap)
tab_model(Pmax.C.trap)

lnc.GP.trap <- lmer(lnc~trap.resid + (1|site), data = covariate.data, subset = rate.type=="GP")
anova(lnc.GP.trap)
summary(lnc.GP.trap)
tab_model(lnc.GP.trap)

lnc.R.trap <- lmer(lnc~trap.resid + (1|site), data = covariate.data, subset = rate.type=="R")
anova(lnc.R.trap)
summary(lnc.R.trap)
tab_model(lnc.R.trap)

lnc.C.trap <- lmer(lnc~trap.resid + (1|site), data = covariate.data, subset = rate.type=="C")
anova(lnc.C.trap)
summary(lnc.C.trap)
tab_model(lnc.C.trap)

sym.trap <- lmer(zoox.per.cm2~trap.resid + (1|site), data = covariate.data)
anova(sym.trap)
tab_model(sym.trap)

chloro.trap <- lmer(chlA.ugcm2~trap.resid + (1|site), data = covariate.data)
anova(chloro.trap)
tab_model(chloro.trap)

AFDW.trap <- lmer(AFDW.mg.cm2.~trap.resid + (1|site), data = covariate.data)
anova(AFDW.trap)
tab_model(AFDW.trap)

N.AT.trap <- lmer(per.N.AT~trap.resid + (1|site), data = covariate.data)
anova(N.AT.trap)
tab_model(N.AT.trap)

N.ST.trap <- lmer(N.ST~trap.resid + (1|site), data = covariate.data)
anova(N.ST.trap)
tab_model(N.ST.trap)

per.N.ST.trap <- lmer(per.N.ST~trap.resid + (1|site), data = covariate.data)
anova(per.N.ST.trap)
tab_model(per.N.ST.trap)

#arrange models to make table of p estimates, CI values to add to paper, sedimentation models below
Topt_models.trap <- tab_model(Topt.GP.trap, Topt.R.trap, Topt.C.trap, pred.labels = c("Intercept", "Sedimentation Rate Residuals"),
          dv.labels = c("Topt_GP", "Topt_R", "Topt_C"), string.ci = "Conf. Int (95%)",
          string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/Topt_trap.doc")

Pmax_models.trap <- tab_model(Pmax.GP.trap, Pmax.R.trap, Pmax.C.trap, pred.labels = c("Intercept", "Sedimentation Rate Residuals"),
          dv.labels = c("umax_GP", "umax_R", "umax_C"), string.ci = "Conf. Int (95%)",
          string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/Pmax_trap.doc")
          
lnc_models.trap <- tab_model(lnc.GP.trap, lnc.R.trap, lnc.C.trap, pred.labels = c("Intercept", "Sedimentation Rate Residuals"),
           dv.labels = c("b(Tc)_GP", "b(Tc)_R", "b(Tc)_C"), string.ci = "Conf. Int (95%)",
          string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/lnc_trap.doc")

intrinsic_models.trap <- tab_model(sym.trap, chloro.trap, AFDW.trap, N.AT.trap, N.ST.trap, per.N.ST.trap, collapse.ci = TRUE, pred.labels = c("Intercept", "Sedimentation Rate Residuals"),
                        dv.labels = c("Endosymbiont Density", "Chlorophyll a Content", "Tissue Biomass", "Coral % N Content", "Endosymbiont N Content Per Cell", "Endosymbiont % N Content"), string.ci = "Conf. Int (95%)",
                        string.p = "P-Value", file = "../../../Documents/CSUN/Thesis Defense/Tables/intrinsic_trap.doc")

