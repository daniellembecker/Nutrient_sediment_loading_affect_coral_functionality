###############################################################################################
##Making map of moorea with imagies and sites and cowplot of PI curves
##created by Danielle Becker 11/14/19
##edited by Danielle Becker 06/16/20


##Install packages
# load packages
library(viridis)
library(tidyverse)
library(lubridate)
library(data.table)
library(ggrepel)
library(ggpubr)
library(ggsn)
library(magick)
library(maps)
library(maptools)
library(ggplot2)
library(grid)
library(ggimage)
library(ggmap)
library(here)

#set wd
now()

#load data frame for nutrient and sedimentation data
mydata <- read_csv("Summer_2019/manuscript.figures/Data/map.data.csv")

#bounding box lowerleftlon, lowerleftlat, upperrightlon, upperrightlat
myLocation <- c(-149.9501,-17.60,-149.7399,-17.45)

#make map and set sources
myMap <- get_map(location=myLocation, source = "stamen", maptype = "terrain") 

#Get the bounding box for the map and reformat it into a data.frame for scalebar
bb <- attr(myMap, "bb")
bb2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))

#make name for sedimentation rate to use superscript in legend title
nameColor <- bquote(atop(Sedimentation~Rate~phantom(),
                         (mg~cm^-2~day^-1)))

#Add the scalebar to a ggmap, need ggsn package
map <- ggmap(myMap) + 
  labs(x = 'Longitude', y = 'Latitude') + 
  geom_point(aes(x = long, y = lat, size = sedimentation.rate, colour = N),  data = mydata) + 
  geom_text(data = mydata, aes(long, lat, label = map.letter), vjust = -1.5) + #add labels for each site 
  scale_size_continuous(limits = c(1.5, 3.5)) + #change range of sed rate scale in legend
  labs(colour = "% Nitrogen Content", size = nameColor) + #rename legend labels
  guides(colour = guide_colourbar(title.vjust = 3)) + #adjust N bar to move down from legend title to make more space
  scale_colour_viridis(direction = 1, option = "A", limits = c(0.65, 0.90)) + #change gradient color 
  scalebar(data = bb2, dist = 2, model  = "WGS84", transform = TRUE, dist_unit = "km", #data are bounding box of map, model gives dataum from google maps, transform recognizes as GPS points as decimal units, location sets location on map, anchor sets bounds of locatio non map
           location = "bottomleft", st.dist = 0.037, anchor = c( x = bb$ll.lon + 0.016, y = bb$ll.lat +0.014)) #add scale bar


north2(map, x=0.20, y=0.36, symbol=3) #add north symbol to map

dev.off




ggsave(filename = "Summer_2019/manuscript.figures/Output/Moorea.map.png", device = "png", width = 10, height = 8)











