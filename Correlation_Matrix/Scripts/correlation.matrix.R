##Using a correlation matrix to see general relationships with multiple plots between all parameters and data frames for environmental and physiological parameters
#correlation matrix is used to investigate the dependence between multiple variables at the same time. The result is a table containing the correlation coefficients between each variable and the others.

#clear dataframe
rm(list=ls())


##Install packages
# load packages
library(tidyverse)
library(PerformanceAnalytics)
library(reshape2)
library(viridis)
library(car)
library(GGally)
library(corrplot)
library(here)

#set wd
here()

#load in data files for parameters
physio.environ.data <- read.csv("Correlation_Matrix/Data/physio.environ.data.csv")
correlation.matrix.data <- read.csv("Correlation_Matrix/Data/TPC.data.csv")

#need to filter out PA39 and PA37 fragments as the symbiont counts were not representative of the sample, the sample had been clumpy
correlation.matrix.data <- correlation.matrix.data%>%
  filter(!(fragment.ID == "PA39"))  %>%
  filter(!(fragment.ID == "PA47"))

#combine parameters for sed and N with environ/physio data sheet to get overall data sheet of parameters
TPC.cor.data <- correlation.matrix.data %>%
  dplyr::select(fragment.ID, Topt, Pmax, rate.type) %>%
  #separate(fragment.ID, c("ID", "trial.type")) %>%
  #select(-c(trial.type)) %>%
  gather(variable,value, (Topt:Pmax))%>% #pivot longer can replace gather
  distinct() %>%
  unite(param,variable,rate.type) %>%
  spread(param,value) #pivot wider is new spread

#left join organized TPC data sheet with rest of the intrinsic physio params and environ params
all.params.data <- left_join(TPC.cor.data, physio.environ.data)

all.params.data$DIN <- all.params.data$NH4 + all.params.data$N.N 
all.params.data$N.P.ratio <- all.params.data$DIN / all.params.data$P

view(all.params.data)

#delete unnecessary columns from the new joined data sheet 
#could also use select to not hardcode, using a pipe to select specific environ variables
#or could make a vector with names and then call that to make a data frame
all.params.data <- all.params.data[, -c(1:11, 13:14, 16:17, 21:22, 25:29)] 

#rename various columns in data frame
names(all.params.data)[1] <- "Mean Light Intensity"
names(all.params.data)[2] <- "Mean Temperature"
names(all.params.data)[3] <- "Ammonium (NH4)"
names(all.params.data)[4] <- "Nitrate (NO3) and Nitrite (NO2)"
names(all.params.data)[5] <- "Phosphate (PO4)"
names(all.params.data)[6] <- "Sedimentation Rate"
names(all.params.data)[7] <- "Percent Nitrogen (N)"
names(all.params.data)[8] <- "DIN:DIP"

#reorder columns to clump nutrients together and environ char. together
all.params.data <- all.params.data[, c(1, 2, 6, 3, 4, 5, 7,8)]

view(all.params.data)


#correlation matrix code for all params 
#compute the correltaion matrix (correlation values organized into data frame)
cormat <- round(cor(all.params.data),4)
head(cormat)

#melt the correlation matrix means it reassembles data frame to be more effective to complete corr matrix
#to long format
melted_cormat <- melt(cormat)
head(melted_cormat)

#visulaize the correlation matrix in general
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

#save general correlation matrix
ggsave(filename = "Correlation_Matrix/Output/gen.correlation.matrix.pdf", device = "pdf", width = 20, height = 10)

# Get lower and upper triangle of the correlation matrix
#Note that, a correlation matrix has redundant information. We’ll use the functions below to set half of it to NA
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#apply upper tri calculation to graphc
upper_tri <- get_upper_tri(cormat)
upper_tri

#melt the correlation matrix
#melt the correlation data and drop the rows with NA values 
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#heatmap of correlation matrix
#negative correlations are in purple color and positive correlations in red
#scale_fill_gradient2 is used with the argument limit = c(-1,1) as correlation coefficients range from -1 to 1
#coord_fixed() : this function ensures that one unit on the x-axis is the same length as one unit on the y-axis
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "midnightblue", high = "firebrick4", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


#edited out because it rearranges the order, want the nutrient concentrations next to eachother in heatmap
#reorder the correlation matrix according to the correlation coefficient
#useful to identify the hidden pattern in the matrix
#hclust for hierarchical clustering order is used 
#helper function to reorder correlation matrix
#reorder_cormat <- function(cormat){
  
#use correlation between variables as distance
  #dd <- as.dist((1-cormat)/2)
  #hc <- hclust(dd)
 #cormat <-cormat[hc$order, hc$order]


#reorder the correlation data visualiztaion
#cormat <- reorder_cormat(cormat)
#upper_tri <- get_upper_tri(cormat)

#melt the correlation matrix
#melted_cormat <- melt(upper_tri, na.rm = TRUE)

#edit the number of sig figs listed in the value column
melted_cormat <- melted_cormat %>% mutate_at(vars(starts_with("value")), funs(round(., 2)))

# Create a ggheatmap with basic characteristics, etc. 
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "firebrick3", mid = "white", high = "dodgerblue3", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


# Print the heatmap
print(ggheatmap)

#add correlation coefficients to the heatmap
#geom_text() to add the correlation coefficients on the graph
#guides() to change the position of the legend title
#if else statement in melted data frame to quotes of black and white to adjust text color 
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 6) +
  theme(axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 18, face="bold", color="black"),
    axis.text.y = element_text(size = 18, face="bold", color="black"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0.8, 0),
    legend.title = element_text(size = 18, face="bold", color="black"),
    legend.text = element_text(size = 20, face="bold", color="black"),
    legend.position = c(0.48, 0.75),
    legend.direction = "horizontal") +
    guides(fill = guide_colorbar(barwidth = 12, barheight = 2, 
                               title.position = "top", title.hjust = 0.5, title.vjust = 1.0))

ggsave(filename = "Correlation_Matrix/Output/final.correlation.matrix.pdf", device = "pdf", width = 10, height = 10)







