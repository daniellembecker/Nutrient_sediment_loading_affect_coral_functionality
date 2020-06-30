#Code for environmental characteristics by site
#Created By: Danielle Becker
#Created On: 06/29/2020
#clear list 

rm(list=ls())

#load libraries 

library(ggplot2)
library(patchwork)
library(PNWColors)
library(lme4)
library(lmerTest)
library(tidyverse)
library(here)

#wd
here()


#load data sheets
site.dat <- read_csv("Site_Characteristics/Data/site.characteristics.data.csv")
temp.dat <- read_csv("Site_Characteristics/Data/temp.dat.csv", col_types = cols(site.letter = col_character()))
light.dat <- read_csv("Site_Characteristics/Data/light.dat.csv")

#making three plots for site characteristics, one of percent cover data (biological), one for physical (temp/light), one for environemntal (nutrients) ()
###BIOLOGICAL###

#filter just biological parameters 
bio.dat <- site.dat %>%
  filter(type.param == "biological")

#rearrange data to specific order

x <- c("p.acuta.cover", "coral.cover", "algal.cover", "CCA.cover", "substrate.cover")

bio.dat <- bio.dat %>%
  mutate(parameter.measured =  factor(parameter.measured, levels = x)) %>%
  arrange(parameter.measured)


#use PNW color palette, build palette
pal=pnw_palette("Moth", type = "discrete")

#make stacked bar plot with % cover on the y axis
a <- ggplot(bio.dat, aes(fill=parameter.measured, y=mean, x=site.letter)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw()+ #overall theme for plot
  scale_fill_manual(values = pal, labels=c(expression(italic("Pocillopora acuta")), "Other Coral Species", "Algae", "CCA",  "Substrate")) + #change gradient color 
  xlab("Site") +  ylab("Benthic Cover (%)") +
  theme(legend.position="top") +
  labs(fill = "") +  
  theme(legend.direction = "horizontal") +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.text=element_text(size=16), axis.text.x=element_text(color="black", size=16), plot.title = element_text(hjust =0.5, color = "black", size=22), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=18), axis.title.y = element_text(color="black", size=18),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

ggsave(filename = "Site_Characteristics/Output/biological.plot.pdf", device = "pdf", width = 12, height = 12)


###PHYSICAL###

#temp data across sites summarised 
data.temp<-temp.dat %>%
  group_by(site.letter) %>% #tells to group by treatment
  summarise(mean=mean(Temp), se=sd(Temp)/sqrt(n())) #calculates mean and se
data.temp

#use PNW color palette, build palette
pal2=pnw_palette("Anemone", type = "discrete")

#make boxplot mean/var for temp and light seperate

b <- ggplot(temp.dat, aes(x = site.letter, y = Temp, fill = site.letter)) +
  geom_boxplot(fill = "cyan4") +
  theme_bw() + #overall theme for plot
  xlab("Site") +  ylab("Temperature (°C)") +
  labs(fill = "") + 
  theme(legend.position =  "none") +
  ylim(22,30) +
  theme(axis.text.x=element_text(color="black", size=16), plot.title = element_text(hjust =0.5, color = "black", size=22), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=18), axis.title.y = element_text(color="black", size=18),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

ggsave(filename = "Site_Characteristics/Output/temp.plot.pdf", device = "pdf", width = 12, height = 12)


#compare  light data across sites summarised 

data.light <-light.dat %>%
  group_by(site.letter) %>% #tells to group by treatment
  summarise(mean=mean(PFD), se=sd(PFD)/sqrt(n())) #calculates mean and se
data.light


c <- ggplot(light.dat, aes(x = site.letter, y = PFD, fill = site.letter)) +
  geom_boxplot(fill = "cyan3") +
  theme_bw() + #overall theme for plot
  xlab("Site") +  ylab(expression(atop("Photon Flux Density", (mu*mol~photons~m^{-2}~s^{-1})))) +
  theme(legend.position =  "none") +
  ylim(0,5000) +
  theme(axis.text.x=element_text(color="black", size=16), plot.title = element_text(hjust =0.5, color = "black", size=22), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=18), axis.title.y = element_text(color="black", size=16),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

ggsave(filename = "Site_Characteristics/Output/light.plot.pdf", device = "pdf", width = 12, height = 12)

###CHEMICAL###
#make plot with nutrient levels per site

#filter to only have watercoloumn nutrients
watercol.all <- site.dat%>%
  filter(type.param == "chemical")

#filter to just have DIN:DIP to place on other plot
DIN.DIP.watercol <- site.dat %>%
  filter(type.param == "ratio")

#use PNW color palette, build palette
pal.2=pnw_palette("Moth", n= 1, type = "discrete")

#make two bar plots for watercol and %N bar plots for water coloum nutrients

d <- ggplot(watercol.all, aes(x=site.letter, y = mean, fill = parameter.measured)) + 
  geom_bar(position=position_dodge(), stat="identity") + #determines the bar width
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se),width=.4,position=position_dodge(.9)) + #adds error bars
  theme_bw() + #overall theme for plot
  scale_fill_manual(values = pal, labels=c(expression(NO["3"]^~~{"-"}~+~NO["2"]^~~{"-"}), expression(NH["4"]^~~{"+"}), expression(PO["4"]^~~{"3-"}))) + #change gradient color 
  xlab("Site") +  ylab(expression(Concentration~(mu*mol~L^{-1}))) +
  labs(fill = "") +  
  scale_y_continuous(expand = c(0,0), limits = c(0,2.0)) +
  theme(legend.text.align = 0) + #make legend text align left
  theme(legend.text=element_text(size=16), axis.text.x=element_text(color="black", size=16), plot.title = element_text(hjust =0.5, color = "black", size=22), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=18), axis.title.y = element_text(color="black", size=18),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

ggsave(filename = "Site_Characteristics/Output/water.col.plot.pdf", device = "pdf", width = 12, height = 12)

e <- ggplot(DIN.DIP.watercol, aes(y=mean, x=site.letter)) + 
  geom_bar(position=position_dodge(), stat="identity", fill = "hotpink3") + #determines the bar width
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se),width=.4,position=position_dodge(.9)) + #adds error bars
  theme_bw() + #overall theme for plot
  xlab("Site") +  ylab("DIN:DIP") +
  labs(fill = "") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,6.8)) +
  theme(axis.text.x=element_text(color="black", size=16), plot.title = element_text(hjust =0.5, color = "black", size=22), axis.text.y=element_text(color="black", size=16), axis.title.x = element_text(color="black", size=18), axis.title.y = element_text(color="black", size=18),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #adjust themes for chart x and y axis labels and axis tick mark labels

ggsave(filename = "Site_Characteristics/Output/per.N.plot.pdf", device = "pdf", width = 12, height = 12)

###########################################################################################################
#make plots for all site condtions

figure <- a  / (b + c) / (d + e) +      #patchwork to combine plots
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 22, face = "bold"))   #edit the lettered text

figure


ggsave(filename = "Site_Characteristics/Output/site.char.pdf", device = "pdf", width = 15, height = 13)


###############################################################################################################
#source script for moorea map to combine with other plot later
source("Map/Scripts/MooreaMap.R")



combined.plots <- map / figure  +      #patchwork to combine plots
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 22, face = "bold"))   #edit the lettered text

combined.plots



ggsave("Site_Characteristics/Output/site.map.pdf", width = 16, height = 14)





