##Processing benthic survey data for Summer 2019
#Created by: Danielle Becker
#Created on: 20191122

#clear dataframe
rm(list=ls())

##Install packages
# load packages
library(tidyverse)
library(car)
library(dplyr)
library(here)

#set wd
here()

#load in meta data sheet
benthic.dat <- read_csv("Summer_2019/benthic.cover/Data/benthic.survey.csv")
view(benthic.dat)

#make all NA values in the data sheet read as 0's

benthic.dat <- benthic.dat %>%
  mutate_all(funs(ifelse(is.na(.), 0, .)))


#make P. acuta as a character
#as.character(benthic.dat$`Pocillopora acuta`)

#add the sum of p.acuta colonies per site 
benthic.dat.cal <- benthic.dat %>% 
  group_by(site.name) %>% 
  summarise(Pocillopora.acuta = sum(Pocillopora.acuta))

#divde p acuta colonies by 40 (quadrats) to get average per quadrat
benthic.dat.cal$per.quad <- benthic.dat.cal$Pocillopora.acuta/40

#multiply the quadrat value by percentage to get percent cover per site 
benthic.dat.cal$percent.cover <- benthic.dat.cal$per.quad*100

#add column with site.block
benthic.dat.cal$site.block  <- c(1, 3, 2, 2, 1, 3)

#add column with site
benthic.dat.cal$site  <- c("A", "B", "A", "B", "B", "A")

#load in the environmental data to compare p. acuta along a nutrient gradient
nut.dat <- read.csv("Summer_2019/parameter.connections/Data/covariate.data.csv")

#combine nut.dat and benthic cover
tot.dat <- left_join(benthic.dat.cal, nut.dat)

#NUTRIENT GRADIENT VS PERCENT COVER
#checking normality of variables
qqp(tot.dat$N.resid, "norm")

qqp(tot.dat$percent.cover, "norm")


#look at relationship of %N to p.acuta percent cover
ggplot(tot.dat, aes(x=N.resid, y=percent.cover)) + #set aesthetics for my x and y usuing means for x and y vars
  geom_point() +
  theme_bw()+ #overall theme for plot
  xlab(expression(bold("Percent Nitrogen (%)"))) + ylab(expression(bold(~bolditalic(Pocillopora~acuta)~ "Percent Cover"))) +
  theme(plot.title = element_text(face = "bold", color = "black", size = 14, hjust = 0.5)) + #adjust themes for plot title
  theme(axis.text.x=element_text(face="bold", color="black", size=12), axis.text.y=element_text(face="bold", color="black", size=12), axis.title.x = element_text(face="bold", color="black", size=14), axis.title.y = element_text(face="bold", color="black", size=14),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  geom_smooth(data=tot.dat, aes(x=N.resid, y=percent.cover), method = "lm") + #add lm relationship with raw data sheet, filter out for rate type represented 
  xlim(-0.095, 0.095)

ggsave(filename = "Summer_2019/benthic.cover/Output/percent.cover.p.acuta.png", device = "png", width = 7, height = 5)

#create linear model to test relationship
percent.cover.lm = lm(percent.cover ~ N.resid, data=tot.dat)


#calculate the r squared and p value
summary(percent.cover.lm) 


#SEDIMENTATION GRADIENT VS PERCENT COVER
##checking normality of variables
qqp(tot.dat$trap.resid, "norm")

qqp(tot.dat$percent.cover, "norm")


#look at relationship of %N to p.acuta percent cover
ggplot(tot.dat, aes(x=trap.resid, y=percent.cover)) + #set aesthetics for my x and y usuing means for x and y vars
  geom_point() +
  theme_bw()+ #overall theme for plot
  xlab(expression(bold("Sedimentation Rate (" *mg *~cm^-2~day^-1*")"))) +  ylab(expression(bold(~bolditalic(Pocillopora~acuta)~ "Percent Cover"))) +
  theme(plot.title = element_text(face = "bold", color = "black", size = 16, hjust = 0.5)) + #adjust themes for plot title
  theme(axis.text.x=element_text(face="bold", color="black", size=12), axis.text.y=element_text(face="bold", color="black", size=12), axis.title.x = element_text(face="bold", color="black", size=14), axis.title.y = element_text(face="bold", color="black", size=14),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + #adjust themes for chart x and y axis labels and axis tick mark labels
  geom_smooth(data=tot.dat, aes(x=trap.resid, y=percent.cover), method = "lm")  #add lm relationship with raw data sheet, filter out for rate type represented 


ggsave(filename = "Summer_2019/benthic.cover/Output/trap.resid.percent.cover.p.acuta.png", device = "png", width = 7, height = 5)

#create linear model to test relationship
percent.cover.lm = lm(percent.cover ~ trap.resid, data=tot.dat)

#calculate the r squared and p value
summary(percent.cover.lm) 

write.csv(tot.dat, file = "Summer_2019/benthic.cover/Data/percent.cover.csv")




