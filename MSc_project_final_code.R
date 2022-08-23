rm(list=ls())#clear workspace

setwd("/Users/leahdenoun/")#set working directory

library(sf)
library(sp)

##DATA IMPORT AND MANAGEMENT##

#get dung removal data from SAFE repository
library(safedata)
set_safe_dir("leahdenoun:\\safedata")
show_record(3247492)#check data and worksheet options for 2011 using Zenodo number
dung_2011<- load_safe_data(3247492, "Dung_Removal")#import
show_record(3247494)##check data and worksheet options for 2015 using Zenodo number
dung_2015<-load_safe_data(3247494, "Dung_Removal")
dung_total<-rbind(dung_2011, dung_2015)#combine dung removal data both years 

#packages for data management 
library(dplyr)
library(tidyr)

#now that we have file with dung removal and points need to add above-ground carbon density (ACD)
#import ACD data 
library(readxl)
ab_CD<-  read_excel('/Users/leahdenoun/Downloads/Tom_swinfield_AGB_LIDAR.xlsx', sheet = 3)
#wrangle to get just sample point and ACD columns
ab_CD <- ab_CD%>%dplyr::select(Location, Numeric...7)%>%slice(5:1057)#extract columns and rows
#get first row as header
names(ab_CD) <- as.character(ab_CD[1,])
ab_CD <- ab_CD[-1,]#delete first row to obtain desired format
ab_CD$agb <- as.numeric(ab_CD$agb)#get ACD as numeric

#separate patch and station to then select based on station for each year
CD_2011<- ab_CD %>%
  filter(ID %in% dung_2011$Site)%>% rename(Site = ID)#select sites in our dataframe
dung_2011<- merge(dung_2011, CD_2011, all= TRUE)#get CD in big dataset 
CD_2015<- ab_CD %>%
  filter(ID %in% dung_2015$Site)%>% rename(Site = ID)#select sites in our dataframe
dung_2015<- merge(dung_2015, CD_2015, all= TRUE)#get CD in big dataset 

#need to format so Site column has numbers like in shp file 
dung_wrangle_2011 <- dung_2011%>%separate(Site, c("Patch", "Station"), "_")%>%separate(Date, c("Year",NA, NA), "-")%>%na.omit(Dung_removed)#at the end take out row with NA from OP1 patch 
dung_wrangle_2015 <- dung_2015%>%separate(Site, c("Patch", "Station"), "_")%>%separate(Date, c("Year",NA, NA), "-")%>%na.omit(Dung_removed) 

#import data shp file
locations<- st_read('/Users/leahdenoun/Downloads/SAFE_core_sampling_stations_UTM50N_WGS84/SAFE_core_sampling_stations_UTM50N_WGS84.shp')
locations<- subset(locations, Fract_Ord=="2")#subset to just have points in dung beetle dataset 
points_2011<- locations%>%dplyr::filter(Site!="LFE", Site!="VJR", Station!="742")#no VJR or LFE data from 2011 and one na so remove those

#need to select only rows with matching site 
#for 2011
removal_sf_2011<- merge(points_2011, dung_wrangle_2011, by.x= 'Station', by.y= 'Station' )
removal_sf_2011$agb <- as.numeric(removal_sf_2011$agb)#get CD from character to numeric format
removal_sp_2011<- as(removal_sf_2011, Class= "Spatial")#transform into a spatial object

#for 2015 different because there is VJR and LFE
removal_sf_2015<- merge(locations, dung_wrangle_2015, by.x= 'Station', by.y= 'Station' )
removal_sf_2015$agb <- as.numeric(removal_sf_2015$agb)
removal_sp_2015<- as(removal_sf_2015, Class= "Spatial")#transform into a spatial object

#reproject to borneo to get correct format to put into function
full_2011 <- st_as_sf(removal_sp_2011, crs = 29873)
full_2015 <- st_as_sf(removal_sp_2015, crs = 29873)

## SAMPLING POINT PAIRS TECHNIQUE ##

#make function to get summary dataframe from sf dataframe
summary_function <- function(x){
  # Calculate all pairs of polygons
  combns_x <- t(combn(length(x$geometry), 2))
  #functions to get variables 
  distance_x <- apply(combns_x, 1, function(y)
    st_distance(x$geometry[[y[1]]], x$geometry[[y[2]]]))#distance between points
  
  rem_x <- apply(combns_x, 1, function(y)
    mean(c(x$Dung_removed[[y[1]]], x$Dung_removed[[y[2]]]))/(sd(c(x$Dung_removed[[y[1]]], x$Dung_removed[[y[2]]]))))#inverse of coefficient of variation
  
  CD_x <- apply(combns_x, 1, function(y) 
    mean(c(x$agb[[y[1]]], x$agb[[y[2]]])))#mean ACD
  
  var_CD_x <- apply(combns_x, 1, function(y) 
    var(c(x$agb[[y[1]]], x$agb[[y[2]]]))) #variance ACD
  #make dataframe
  Summary_df <- cbind.data.frame(From=as.character(x$Station[combns_x[, 1]]), 
                                 To=as.character(x$Station[combns_x[, 2]]), 
                                 Distance=distance_x, 
                                 Stability=rem_x,
                                 Mean_CD=CD_x,
                                 Variance_CD=var_CD_x)
  return(Summary_df)
}

#run function for both years = all pairs of points compared regardless of site 
pairs_2011_summary <- summary_function(full_2011)
pairs_2011_summary$Year <- as.factor(2011)#add year to then include it in models 
pairs_2015_summary <- summary_function(full_2015)
pairs_2015_summary$Year <- as.factor(2015)

pairs_full <- rbind(pairs_2011_summary, pairs_2015_summary)
#see if mean and variance of biomass are correlated
cor.test(pairs_full$Mean_CD, pairs_full$Variance_CD)#they are correlated so just use one 

plot(log10(pairs_full$Distance), log10(pairs_full$Stability))#notice one big visible outlier
outliers_pairs <- subset(pairs_full, log10(Stability)>10)#one numbered outlier and 27 Inf values pruned_small <- subset(small, log10(Stability)<12)# version without outliers 
pruned_pairs <- subset(pairs_full, log10(Stability)<20)#exclude the infinite values 
pruned_pairs$Year <- as.factor(pruned_pairs$Year)#get year as a factor
hist(log10(pruned_pairs$Stability))#many very low values 
hist(log10(pruned_pairs$Distance))#few low values 
hist(log10(pruned_pairs$Variance_CD))#few low values 
plot(log10(pruned_pairs$Distance), log10(pruned_pairs$Stability))#much clearer can see potential trend of increasing stability with distance 

#save df used in analysis to avoid running everything when close rstudio
write.csv(pruned_pairs,"\\Users\\leahdenoun\\Downloads\\pruned_pairs.csv", row.names = FALSE)


# quantile regression
library(quantreg)

#model selection to see if eed to keep year as a covariate or not
pair_full <- rq(log10(Stability) ~ log10(Distance) * log10(Variance_CD) * Year, data = pruned_pairs, tau = 0.3)
pair_model_2 <- rq(log10(Stability) ~ log10(Distance) * log10(Variance_CD) + Year, data = pruned_pairs, tau =0.3)
anova.rq(pair_full, pair_model_2)#did individual models with quantiles 0.2-0.9
AIC.rq(pair_full)
AIC.rq(pair_model_2)
pruned_pairs$logD<- log10(pruned_pairs$Distance)
pruned_pairs$logV <- log10(pruned_pairs$Variance_CD)
pruned_pairs$logS <- log10(pruned_pairs$Stability) 
#model of this analysis 
pair_model <- rq(logS ~ logD * logV * Year, data = pruned_pairs, tau = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.7, 0.8, 0.9))
summary_pair_model<-summary(pair_model, se="boot") 
summary_pair_model
plot(summary_pair_model)


#getting pseudo r squared for each quantile regression model 
base_fit_pairs <- rq(log10(Stability)~1,tau=1:9/10, data=pruned_pairs)
rho <- function(u,tau=.5)u*(tau - (u < 0))
R1_pairs <- 1 - pair_model$rho/base_fit_pairs$rho #most explanation for 80th quantile 

#get AIC of each tau value model to see
AIC.rq(pair_model)#70th quantile model better 

#get logged values in data set for plot model and plotting (predict function specifically needs it)

#simple plot
pair_plot_model <- rq(logS ~ logD * logV * Year, data = pruned_pairs, tau = 0.8)#plotting only for one percentile 

#import packages to plot
library(sjPlot)
library(ggplot2)
#set theme and plot facet by year for mean of variance terciles: scatterplot and regression lines
set_theme(geom.outline.color = "antiquewhite4", 
          geom.outline.size = 1,
          axis.title.size = 1.5,
          axis.textsize = 1.2,
          legend.size = 1,
          legend.title.size = 1.2,
          legend.item.size = 1.3,
          geom.label.size = 3,
          base = theme_light())
plot_model(pair_plot_model, type="pred",terms=c("logD","logV[2.647,3.087,3.372]", "Year"), show.data =  TRUE, colors = c("black", "chocolate", "chartreuse4"),
           axis.title = c( "log10(Distance) (m)", "log10(Stability)"), legend.title = "log10(VCD)(t/ha)", title = "", dot.size = 0.7)+
  ylim(-0.5, 4)+
  theme(strip.text = element_text(size= 25, color = "black"))

#making a plot of regression estimates for each year 
#plot coefficients
pair_coefs <-  read_excel('/Users/leahdenoun/MSc_analysis_coefs.xlsx', sheet = 1)#manually inputted and summed estimates into excel sheet
pairs_coefs_f <- mutate(pair_coefs, quantile = quantile, intercept = as.numeric(intercept), log10D = as.numeric(log10D), log10VCD = as.numeric(log10VBIO), logDxlogVCD= as.numeric(logDxlogVBIO), Year=as.factor(Year))
#now in correct format 

#make four plots separately 
p_intercepts <-ggplot(data= pairs_coefs_f, aes(quantile, intercept, group=Year))+
  geom_point(aes(shape=Year), alpha=0.5, size = 2)+
  geom_line(aes(color = Year))+
  theme_light()+
  scale_x_continuous(breaks=seq(10, 90, 10))+
  labs(y="Intercept", x="Quantile (%)")+
  scale_color_manual(values = c("chocolate", "chartreuse4"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size= 16), legend.text = element_text(size=14))         


p_logD <-ggplot(data= pairs_coefs_f, aes(quantile, log10D, group=Year))+
  geom_point(aes(shape=Year), alpha=0.5, size = 2)+
  geom_line(aes(color = Year))+
  theme_light()+
  scale_x_continuous(breaks=seq(10, 90, 10))+
  labs(y="log10(Distance)", x="Quantile (%)")+
  scale_color_manual(values = c("chocolate", "chartreuse4"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size= 16), legend.text = element_text(size=14))         


p_logVCD <-ggplot(data= pairs_coefs_f, aes(quantile, log10VCD, group=Year))+
  geom_point(aes(shape=Year), alpha=0.5, size = 2)+
  geom_line(aes(color = Year))+
  theme_light()+
  scale_x_continuous(breaks=seq(10, 90, 10))+
  labs(y="log10(VCD)", x="Quantile (%)")+
  scale_color_manual(values = c("chocolate", "chartreuse4"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size= 16), legend.text = element_text(size=14))+
  ylim(-0.3, 0.2)


p_interaction <-ggplot(data= pairs_coefs_f, aes(quantile, logDxlogVCD, group=Year))+
  geom_point(aes(shape=Year), alpha=0.5, size = 2)+
  geom_line(aes(color = Year))+
  theme_light()+
  scale_x_continuous(breaks=seq(10, 90, 10))+
  labs(y="log10(Distance)*log10(VCD)", x="Quantile (%)")+
  scale_color_manual(values = c("chocolate", "chartreuse4"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size= 16), legend.text = element_text(size=14))+
  ylim(-0.04, 0.06)

#get plots in singular figure 
library(ggpubr)
figure_1 <- ggarrange(p_intercepts, p_logD, p_logVCD, p_interaction, 
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2, align= "v", 
                      common.legend = TRUE)
figure_1

## POINT GROUPINGS AT INCREASING SPATIAL SCALE TECHNIQUE ## 

#get coordinates in dataset to then aggregate based on distance 
clust_data_2011 <- full_2011 %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2])

clust_data_2015 <- full_2015 %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2])

#now function for all distances 
max(pairs_full$Distance)#find max area value for clusters: about 83775m 

clustering_function <- function(m){
  clusters_list <- list()
  for (j in seq(150, 83775, by = 150)){ #start at 150 since minimum distance between two points is 138 meters
    data_sub <-st_is_within_distance(m, dist = j)#cluster of points based on distance
    data_sub_list = list()
    for (i in seq(1:nrow(data_sub))){
      
      subset_x <- m %>% filter(row_number() %in% data_sub[[i]]) #subset cluster rows from big dataset
      rem_x <- mean(subset_x$Dung_removed/(sd(subset_x$Dung_removed))) #calculate stability of cluster
      CD_x <- mean(subset_x$agb)#mean biomass of cluster
      var_CD_x <- var(subset_x$agb)#variance in aboveground biomass of cluster 
      #make dataframe
      temp_clust_df <- cbind.data.frame(  Stability=rem_x,
                                          Mean_CD=CD_x,
                                          Variance_CD=var_CD_x,
                                          Points = nrow(subset_x))
      data_sub_list[[i]] <- temp_clust_df # add to sub list
    }
    clusters_data <- dplyr::bind_rows(data_sub_list)#join lists 
    clusters_data$Area <- j #add area column to know distance of aggregation in dataframe
    clusters_list[[j]] <- clusters_data# add to final list
  }
  return(clusters_list)
}

clust_fun_2011 <-clustering_function(clust_data_2011) #run function for 2011 data
clusters_temp_2011 <-dplyr::bind_rows(clust_fun_2011) #bind rows of all lists
clust_fun_2015 <- clustering_function(clust_data_2015)# run function for 2015 data
clusters_temp_2015 <-dplyr::bind_rows(clust_fun_2015) # bind rows 


clusters_final_2011 <-  clusters_temp_2011[!duplicated(clusters_temp_2011[c('Stability', 'Mean_CD', 'Variance_CD')]),] #takes out duplicates and keeps the one with lowest area value (first occurence)
clusters_final_cleaned_2011 <- na.omit(clusters_final_2011)#3 NA values so remove them 
clusters_final_cleaned_2011$Year <- 2011


clusters_final_2015 <-  clusters_temp_2015[!duplicated(clusters_temp_2015[c('Stability', 'Mean_CD', 'Variance_CD')]),] #takes out duplicates and keeps the one with lowest area value (first occurence)
clusters_final_cleaned_2015 <- na.omit(clusters_final_2015)#2 NA values taken out
clusters_final_cleaned_2015$Year <- 2015

clusters_final_cleaned_150m <- rbind(clusters_final_cleaned_2011, clusters_final_cleaned_2015)#join years
#save as file because took a while to compute the cluster_list 
write.csv(clusters_final_cleaned_150m,"\\Users\\leahdenoun\\Downloads\\clusters_final_cleaned_150m.csv", row.names = FALSE)

##QUANTILE REGRESSION FOR SECOND TECHNIQUE ##

#basic plots to get overall idea and see distributions 
hist(log10(clusters_final_cleaned_150m$Stability))
hist(log10(clusters_final_cleaned_150m$Area))
hist(log10(clusters_final_cleaned_150m$Variance_CD))
plot(log10(clusters_final_cleaned_150m$Area), log10(clusters_final_cleaned_150m$Stability))
clusters_final_cleaned_150m$Year <- as.factor(clusters_final_cleaned_150m$Year)#get year as factor

#see if mean CD and varCD correlated
cor.test(clusters_final_cleaned_150m$Mean_CD, clusters_final_cleaned_150m$Variance_CD)

#see if stability linked to number of points
cor.test(clusters_final_cleaned_150m$Stability, clusters_final_cleaned_150m$Points)

#model selection to see if year needs to be included as a covariate or not
clust_full <- rq(log10(Stability) ~ log10(Area) * log10(Variance_CD) * Year, data = clusters_final_cleaned_150m, tau = 0.1)
clust_model_2 <- rq(log10(Stability) ~ log10(Area) * log10(Variance_CD) + Year, data = clusters_final_cleaned_150m, tau =0.1)
anova.rq(clust_full, clust_model_2)#individually went through quantiles, better to keep year as a covariate 
AIC.rq(clust_full)
AIC.rq(clust_model_2)

clust_model <- rq(log10(Stability) ~ log10(Area) * log10(Variance_CD) * Year, data = clusters_final_cleaned_150m, tau = 1:9/10)
summary_clust_model<-summary(clust_model, se="boot")#get estimates using ootstrapping technique
summary_clust_model
plot(summary_clust_model)


# see pseudo r squared and AIC values for this analysis 
base_fit_clust <- rq(log10(Stability)~1,tau=1:9/10, data=clusters_final_cleaned_150m)
rho <- function(u,tau=.5)u*(tau - (u < 0))
R1_clust <- 1 - clust_model$rho/base_fit_clust$rho#most explained for 90th quantile

#get AIC of each tau value model to see
AIC.rq(clust_model)#70th quantile model better 

#plot
#get variables for the plot
clusters_final_cleaned_150m$logS <- log10(clusters_final_cleaned_150m$Stability)
clusters_final_cleaned_150m$logA <- log10(clusters_final_cleaned_150m$Area)
clusters_final_cleaned_150m$logVCD <- log10(clusters_final_cleaned_150m$Variance_CD)
clust_plot_model <- rq(logS ~ logA * logVCD * Year, data = clusters_final_cleaned_150m, tau = 0.9)

library(sjPlot)
library(ggplot2)
set_theme(geom.outline.color = "antiquewhite4", 
          geom.outline.size = 1,
          axis.title.size = 1.5,
          axis.textsize = 1.2,
          legend.size = 1,
          legend.title.size = 1.2,
          legend.item.size = 1.3,
          geom.label.size = 3,
          base = theme_light())
plot_model(clust_plot_model, type="pred",terms=c("logA","logVCD[2.647,3.087,3.372]", "Year"), show.data =  TRUE, colors = c("black", "chocolate", "chartreuse4"),#other palette option is met.brewer("Derain")
           axis.title = c( "log10(Distance) (m)", "log10(Stability)"), legend.title = "log10(VCD)(t/ha)", title = "", dot.size = 0.7)+
  theme(strip.text = element_text(size= 25, color = "black"))

#plot coefficient estimates
clust_coefs <-  read_excel('/Users/leahdenoun/MSc_analysis_coefs.xlsx', sheet = 2)#import excel sheet with manually entered and summed coefficients by year
clust_coefs_f <- mutate(clust_coefs, quantile = quantile, intercept = as.numeric(intercept), log10A = as.numeric(log10A), log10VCD = as.numeric(log10VBIO), logAxlogVCD= as.numeric(logAxlogVBIO), Year=as.factor(Year))
#get in right structure 

#make pots separately and arrange them
c_intercepts <-ggplot(data= clust_coefs_f, aes(quantile, intercept, group=Year))+
  geom_point(aes(shape=Year), alpha=0.5, size = 2)+
  geom_line(aes(color = Year))+
  theme_light()+
  scale_x_continuous(breaks=seq(10, 90, 10))+
  labs(y="Intercept", x="Quantile (%)")+
  scale_color_manual(values = c("chocolate", "chartreuse4"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size= 16), legend.text = element_text(size=14))+
  ylim(0, 4)


c_logA <-ggplot(data= clust_coefs_f, aes(quantile, log10A, group=Year))+
  geom_point(aes(shape=Year), alpha=0.5, size = 2)+
  geom_line(aes(color = Year))+
  theme_light()+
  scale_x_continuous(breaks=seq(10, 90, 10))+
  labs(y="log10(Distance)", x="Quantile (%)")+
  scale_color_manual(values = c("chocolate", "chartreuse4"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size= 16), legend.text = element_text(size=14))+
  ylim(-1.2, 0)


c_logVCD <-ggplot(data= clust_coefs_f, aes(quantile, log10VCD, group=Year))+
  geom_point(aes(shape=Year), alpha=0.5, size = 2)+
  geom_line(aes(color = Year))+
  theme_light()+
  scale_x_continuous(breaks=seq(10, 90, 10))+
  labs(y="log10(VCD)", x="Quantile (%)")+
  scale_color_manual(values = c("chocolate", "chartreuse4"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size= 16), legend.text = element_text(size=14))+ 
  ylim(-1, 0.1)


c_interaction <-ggplot(data= clust_coefs_f, aes(quantile, logAxlogVCD, group=Year))+
  geom_point(aes(shape=Year), alpha=0.5, size = 2)+
  geom_line(aes(color = Year))+
  theme_light()+
  scale_x_continuous(breaks=seq(10, 90, 10))+
  labs(y="log10(Distance)*log10(VCD)", x="Quantile (%)")+
  scale_color_manual(values = c("chocolate", "chartreuse4"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size= 16), legend.text = element_text(size=14))+
  ylim(-0.1, 0.3)


library(ggpubr)
figure_2 <- ggarrange(c_intercepts, c_logA, c_logVCD, c_interaction, 
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2, align= "v", 
                      common.legend = TRUE)

figure_2


