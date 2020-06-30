##Practicing TPC

rm(list=ls())

##Install packages
# load packages
# load packages
library(nls.multstart)
library(sjPlot)
library(broom)
library(patchwork)
library(tidyverse)
library(nlstools)
library(dplyr)
library(car)
library(nls2)
library(car)
library(here)
library(lme4)
library(lmerTest)

#set wd
here()


#load data

photo.data <- read.csv("TPC_curves/Data/Photo.T.csv")
photo.data$X <- NULL
View(photo.data)
glimpse(photo.data)

#bring in calcification data sheet and rbind join to Photo.T 

calcification.data <- read.csv("TPC_curves/Data/calcification1.csv")
calcification.data$X <- NULL
mydata <- rbind(photo.data, calcification.data)
view(mydata)

write.csv(mydata, 'TPC_curves/Data/mydata.csv') 


mydata <- mydata %>% #filtering out NP, removes it from the list 
  filter(rate.type !="NP") 

#filter out certain outlier points for calcification and respo
mydata <- mydata %>%
  filter(!(fragment.ID == "PA2_C" & temp.Cat == 37)) %>%
  filter(!(fragment.ID == "PA9_C" & temp.Cat == 37)) %>%
  filter(!(fragment.ID == "PA2_D" & temp.Cat == 39)) %>%
  filter(!(fragment.ID == "PA52_C" & temp.Cat == 20)) %>%
  filter(!(fragment.ID == "PA54_C" & temp.Cat == 20)) %>%
  filter(!(fragment.ID == "PA2_D" & temp.Cat == 37))
  

mydata$log.rate <- log(mydata$umol.cm2.hr + 1)  #logging and adding 0.1 because a log of zero does not exist

# convert temp to K
mydata$K<-mydata$Temp.C + 273.15

# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp))) #units are eV/K, electrovolts/Kelvin
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}

#subset dataset
#  PA1 <- filter(mydata, individual.ID == 'PA1' & rate.type == "GP")
#  View(PA1)
# #
# # #run nls_multstart
# fit <- nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 26),
#                       data = PA1,
#                       iter = 500,
#                       start_lower = c(lnc = -2, E = 0.1, Eh = 0.2, Th = 285),
#                       start_upper = c(lnc = 3, E = 2, Eh = 5, Th = 330),
#                       supp_errors = 'Y',
#                       na.action = na.omit,
#                       lower = c(lnc = -2, E = 0, Eh = 0, Th = 0))
# 
#  fit


# fit over each set of groupings
fits <- mydata %>%
  group_by(., rate.type, fragment.ID) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 26.8),
                                                data = .x,
                                                iter = 1000,
                                                start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                                                start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                                supp_errors = 'Y',
                                                na.action = na.omit,
                                                lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))



#get r2, extract predit and pull out and join lists, shows o vs predcited and r2
#PredictedPA1_D <- predict(fits$fit[[1]])
#ObservedPA1_D <- mydata$umol.cm2.hr[mydata$fragment.ID == "PA1_D"]

#po <- lm(PredictedPA1_D ~ ObservedPA1_D)
#summary(po)
#plot(ObservedPA1_D,PredictedPA1_D)
#abline(po)
#legend("topleft", bty="n", legend=paste("r2 =", format(summary(po)$adj.r.squared, digits=4)))



# look at a single fit
summary(fits$fit[[1]])

# look at output object
select(fits, fragment.ID, data, fit)  

# get summary info
info <- fits %>%
  unnest_legacy(fit %>% map(glance))

# get params
params <- fits %>%
  unnest_legacy(fit %>% map(tidy))

#left join params with meta data file to have recovery block and site block in file
metadata <- read.csv(file="TPC_curves/Data/metadata.csv", header=T) #read in metadata file to add site block and recovery block
params <- left_join(params, metadata)

# get confidence intervals
CI <- fits %>% 
  unnest_legacy(fit %>% map(~ confint2(.x) %>%
                              data.frame() %>%
                              rename(., conf.A = X2.5.., conf.B = X97.5..))) %>%
  group_by(., fragment.ID) %>%
  mutate(., term = c('lnc', 'E', 'Eh', 'Th')) %>%
  ungroup()

#colnames(CI)[4:5]<-c("conf.low", "conf.high") # rename columns

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
  unnest_legacy(fit %>% map(augment))

#select(info, fragment.ID, logLik, AIC, BIC, deviance, df.residual)

# new data frame of predictions, do this to set a sequence to make a smooth curve with your prediction points
new_preds <- mydata %>%  do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 150), stringsAsFactors = FALSE)) #setting a specific sequence so you can have a smooth curve

# max and min for each curve
max_min <- mydata %>% group_by(fragment.ID) %>%
  summarise(., min_K = min(K), max_K = max(K)) %>%
  ungroup()

# create new predictions
preds2 <- fits %>%
  unnest_legacy(fit %>% map(augment, newdata = new_preds)) %>%
  merge(., max_min, by = "fragment.ID") %>%
  group_by(., fragment.ID) %>%
  filter(., K > unique(min_K) & K < unique(max_K)) %>%
  rename(., ln.rate = .fitted) %>%
  ungroup()

#left join preds2 with meta data file to have recovery block and site block in file
preds2 <- left_join(preds2, metadata)

#want to do ggplot where we look at individual curves for calcification rates of each fragment
mydataC <- mydata %>%
  filter(rate.type =="C")

preds2C<- preds2 %>%
  filter(rate.type =="C")

  ggplot() +
  geom_point(aes(K - 273.15, log.rate, col = rate.type), size = 2, mydataC) +
  geom_line(aes(K - 273.15, ln.rate, col = rate.type, group = fragment.ID), alpha = 0.5, preds2C) +
  facet_wrap(~ fragment.ID, labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = c('green4', 'blue', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression("Log " *mu*"mol (" *cm^-2 *hr^-1*")")) +
  xlab('Temperature (ºC)') +
  geom_hline(yintercept=0, color = "red") +
  theme(legend.position = c(0.91, 0.85))+
  labs(color = "Rate Type")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

  ggsave(filename = "TPC_curves/Output/calcindiv.curves.pdf", device = "pdf", width = 10, height = 10)
  
#want to do ggplot where we look at individual curves for respiration rates of each fragment
  mydataR <- mydata %>%
    filter(rate.type =="R")
  
  preds2R<- preds2 %>%
    filter(rate.type =="R")
  
  ggplot() +
    geom_point(aes(K - 273.15, log.rate, col = rate.type), size = 2, mydataR) +
    geom_line(aes(K - 273.15, ln.rate, col = rate.type, group = fragment.ID), alpha = 0.5, preds2R) +
    facet_wrap(~ fragment.ID, labeller = labeller(.multi_line = FALSE)) +
    scale_colour_manual(values = c('green4', 'blue', 'black')) +
    theme_bw(base_size = 12, base_family = 'Helvetica') +
    ylab(expression("Log " *mu*"mol (" *cm^-2 *hr^-1*")")) +
    xlab('Temperature (ºC)') +
    geom_hline(yintercept=0, color = "red") +
    theme(legend.position = c(0.91, 0.85))+
    labs(color = "Rate Type")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  ggsave(filename = "TPC_curves/Output/respindiv.curves.pdf", device = "pdf", width = 10, height = 10)
  
#want to do ggplot where we look at individual curves for photosynthesis rates of each fragment
  mydataGP <- mydata %>%
    filter(rate.type =="GP")
  
  preds2GP<- preds2 %>%
    filter(rate.type =="GP")
  
  ggplot() +
    geom_point(aes(K - 273.15, log.rate, col = rate.type), size = 2, mydataGP) +
    geom_line(aes(K - 273.15, ln.rate, col = rate.type, group = fragment.ID), alpha = 0.5, preds2GP) +
    facet_wrap(~ fragment.ID, labeller = labeller(.multi_line = FALSE)) +
    scale_colour_manual(values = c('green4', 'blue', 'black')) +
    theme_bw(base_size = 12, base_family = 'Helvetica') +
    ylab(expression("Log " *mu*"mol (" *cm^-2 *hr^-1*")")) +
    xlab('Temperature (ºC)') +
    geom_hline(yintercept=0, color = "red") +
    theme(legend.position = c(0.91, 0.85))+
    labs(color = "Rate Type")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  ggsave(filename = "TPC_curves/Output/photoindiv.curves.pdf", device = "pdf", width = 10, height = 10)

#make site.block and recoevry block a factor
  mydata$site.block <- as.factor(mydata$site.block)
  mydata$recovery.block <-as.factor(mydata$recovery.block)

  view(mydata)
    
#make set labels fro site block params to relabel in facet grid
labels <- c("1" = "Eastern", "2" = "Central", "3" = "Western")

# plot all P, R and C values in TPCs
ggplot() +
  geom_point(aes(K - 273.15, log.rate, col = rate.type), size = 2, mydata) +
  geom_line(aes(K - 273.15, ln.rate, col = rate.type, group = fragment.ID), alpha = 0.5, preds2) +
  scale_colour_manual(values = c('green4', 'blue', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression("Log " *mu*"mol (" *cm^-2 *hr^-1*")")) +
  xlab('Temperature (ºC)') +
  facet_grid(~ site.block, labeller = labeller(site.block=labels)) + #rename facet wrap headings
  theme(legend.position = "right") +
  labs(col = "Rate Type") +
  #labs(pch = "Region") + #rename legend title
  #guides(color = guide_legend(order = 1), pch = guide_legend(order = 2)) + #change order of legends 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  
ggsave(filename = "TPC_curves/Output/TPCcurves.pdf", device = "pdf", width = 40, height = 20)


# function for calculating Topt

get_topt <- function(E, Th, Eh){
  return((Eh*Th)/(Eh + (8.62e-05 *Th*log((Eh/E) - 1))))
}

# calc topts for all 
Topt_data <- params %>%
  dplyr::select(fragment.ID,  term, estimate, site, rate.type) %>%
  spread(term, estimate) %>%
  mutate(Topt = get_topt(E, Th, Eh)) %>% 
  left_join(.,metadata)
  #group_by(., rate.type, site, fragment.ID)

write.csv(Topt_data, 'TPC_curves/Data/Topt_data.csv') 


#get temerature back in celcius not K
Topt_data$Topt <- Topt_data$Topt - 273.15 

#drop NP
Topt_data$rate.type <- droplevels(Topt_data$rate.type)

#data summary for mean and se for topt per site 
data.summary <- Topt_data %>%
  group_by(site.letter, rate.type) %>% #tells to group by these two factors
  summarise(mean=mean(Topt), se=sd(Topt)/sqrt(n())) #calculates mean and s.e.
data.summary

data.summary$rate.type<- factor(data.summary$rate.type, levels = c("R", "GP", "C")) 

ggplot(data.summary, aes(x=site.letter, y=mean, col = rate.type, group=factor(site.letter))) +
  geom_point(position=position_dodge()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width= 0.4,position=position_dodge(0.9)) +
  labs(x = "Site", y = "Temperature (°C)") + 
  theme_bw()+
  theme(legend.text.align = 0) + #make legend text align left
  theme(legend.text=element_text(size=16), axis.text.x=element_text(color="black", size=16), plot.title = element_text(hjust =0.5, color = "black", size=22), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=18), axis.title.y = element_text(color="black", size=18),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

 
ggsave(filename = "TPC_curves/Output/Topt_graph.pdf", device = "pdf", width = 7, height = 5)


#data summary for mean and se for lnc per site
data.summary<-Topt_data %>%
  group_by(site.letter, rate.type) %>% #tells to group by these two factors
  summarise(mean=mean(lnc), se=sd(lnc)/sqrt(n())) #calculates mean and s.e.
data.summary

ggplot(data.summary, aes(x=site.letter, y=mean, col = rate.type, group=factor(site.letter))) +
  geom_point(position=position_dodge()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width= 0.4,position=position_dodge(0.9)) +
  labs(x = "Site", y = "Rate at a constant temperature (bTc)") + 
  theme_bw()+
  theme(legend.text.align = 0) + #make legend text align left
  theme(legend.text=element_text(size=16), axis.text.x=element_text(color="black", size=16), plot.title = element_text(hjust =0.5, color = "black", size=22), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=18), axis.title.y = element_text(color="black", size=18),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

ggsave(filename = "TPC_curves/Output/bTc_graph.pdf", device = "pdf", width = 8, height = 5)
    

#calculate Pmax values between sites

Pmax_data <- Topt_data %>%
  mutate(Pmax = schoolfield_high(lnc = lnc, E = E, Th = Th, Eh = Eh, temp = Topt + 273.15, Tc = 26)) %>% #add in factors that make up schoolfield function, reference topt to get pmax
  group_by(., rate.type, site, fragment.ID, site.letter)

write.csv(Pmax_data, 'TPC_curves/Data/Pmax_data.csv') # export all the uptake rates


#data summary for mean and se for pmax per site
data.summary<-Pmax_data %>%
  group_by(site.letter, rate.type) %>% #tells to group by these two factors
  summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.
data.summary

ggplot(data.summary, aes(x=site.letter, y=mean, col = rate.type, group=factor(site.letter))) +
  geom_point(position=position_dodge()) +
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width= 0.4,position=position_dodge(0.9)) +
  labs(x = "Site", y = "Maximum performance") + 
  theme_bw()+
  theme(legend.text.align = 0) + #make legend text align left
  theme(legend.text=element_text(size=16), axis.text.x=element_text(color="black", size=16), plot.title = element_text(hjust =0.5, color = "black", size=22), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=18), axis.title.y = element_text(color="black", size=18),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

ggsave(filename = "TPC_curves/Output/Pmax_graph.pdf", device = "pdf", width = 8, height = 5)



#join Pmax data to Topt to have all TPC params in one data sheet
TPC.data <- left_join(Pmax_data, Topt_data)

#rename fragment ID to individul ID
colnames(TPC.data)[colnames(TPC.data)=="fragment.ID"] <- "individual.ID"

#make fragment ID column with no _
TPC.data <- TPC.data %>%
  separate(individual.ID, c("fragment.ID")) 

#have to filter out PA39 and PA47 due to unreal symbiont data 
TPC.data <- TPC.data %>%
  filter(!(fragment.ID == "PA39"))  %>%
  filter(!(fragment.ID == "PA47"))
  

write.csv(TPC.data, 'Correlation_Matrix/Data/TPC.data.csv') 









