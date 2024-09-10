#### Revision for CEBP submission:
#### Luo et al. Where should the cancer control interventions target: 
#### A geospatial hotspot analysis for major cancer mortality 2018-2022 in the U.S.


# install.packages(c('sf', 'sfdep', 'spdep', 'sfdep', 'stringr', 'sfdep'))
require(data.table)
require(sf) 
require(sfdep)
require(spdep)
require(dplyr)
require(tidyr)
require(ggplot2) 
require(stringr)
# pkg for spatial reg
require(spaMM)
require(RSpectra)
library(spatialreg) 

# setwd('')
dd = readRDS(file='hotspot_Wonder_2022.RDS')
 

###  data dictionary  
# [1] "FIPS"         # county FIPS           
# [2] "pop2020M"     # 2020 census population male (R package tidycensus)
# [3] "pop2020F"     # 2020 census population female 
# [4] "pct_pov_2018"           
# [5] "med_age_2000"           
# [6] "med_age"                
# [7] "MFratio"                
## 2018 County Health Ranking cancer risk factors
## hcbips://www.countyhealthrankings.org/health-data/methodology-and-sources/data-documentation
# [8] "smoking"                
# [9] "obese"                  
# [10] "foodenv"                
# [11] "physicallyinactive"     
# [12] "exdrinking"             
# [13] "uninsured"              
# [14] "PCPrate"                
# [15] "mammography"            
# [16] "college"                
# [17] "unemployed"             
# [18] "avgdailyPM2.5"          
# [19] "County"       # county name and state             
# [20] "Cancer"       # 4 major cancer types            
# [21] "Deaths_2018_2022"       
# [22] "Population_2018_2022"   
# [23] "Deaths"                 
# [24] "pop"          # at-risk popu                
# [25] "death_rate"   # calculated annual mortality rates per 100k           
# [26] "suppressed"   # death count suppressed?          
# [27] "death_rate_age_adj_2000" 
# [28] "label_G_O"        # hotspot label, primary analysis          
# [29] "label_G_O_imp"    # hotspot label, sensitivity analysis: include imputed counties       
# [30] "label_G_O_age_adj_2000" # hotspot label, sensitivity analysis: age-adjust to 2000
# [31] "label_G_O_Rook"   # hotspot label, sensitivity analysis: use Rook contiguity neighbor      
# [32] "death_rate_Rook"        
# [33] "State" 
        

 
################ hotspot analysis: use Getis-Ord G statistics  ################
# calculate annual death_rate 2018_2022, per 100k
cancer.names = c('Colorectal','Lung','Breast','Prostate')
Cancer.names = c('colorectal','lung','breast','prostate')

## suppression % by cancer type
dd[,mean(suppressed),by=Cancer]

# only use US mainland for hotspot analysis 
cb <- st_read("../git/cb_2018_us_county_500k/cb_2018_us_county_500k.shp")
cb = cb %>% mutate(FIPS=as.numeric(GEOID))
cb = cb %>% filter(FIPS%in%dd$FIPS)
(nc = nrow(cb)) # N counties = 3108 (mainland US)

# color labels for hotspot clusters 
fc = c("red","orange","lightgrey","deepskyblue1","blue","white")

### Hotspot analysis
set.seed(123)
# some analysis code from 
# https://rpubs.com/heatherleeleary/hotspot_getisOrd_tut
res = list()  # all analysis results
for(ii in 1:4){ 
  ca = cancer.names[ii] 
  Ca = Cancer.names[ii] 
  cat(ca, '...\n')
  cb$suppressed = cb$label_G_O = cb$label_G_O_imp = cb$label_G_O_age_adj_2000 = cb$label_G_O_Rook = NULL
  ddi = dd[Cancer==ca, ]
  ddi$label_G_O = ddi$label_G_O_imp = ddi$label_G_O_age_adj_2000 = ddi$label_G_O_Rook = NULL
  cbi = ddi[data.table(FIPS=as.numeric(cb$GEOID)),,on='FIPS'] 
  cbi[,Mortality:=death_rate]
  tes_subset = cbind(cb, cbi[,-1])
  
  ## removed empty neighbor sets
  list_nb <- poly2nb(tes_subset, queen = T)
  empty_nb <- which(card(list_nb) == 0)
  tes_subset <- tes_subset[-empty_nb, ]
  empty_polygons <- tes_subset[empty_nb, ]
  # empty_polygons$NAME  # print neighborhood names # n=3
  # Identify neighbors with queen contiguity (edge/vertex touching)
  tes_nb <- poly2nb(tes_subset, queen = TRUE)
  # tes_nb <- poly2nb(tes_subset, queen = F) # Rook contiguity (only edge touching)
  
  # nb <- dnearneigh(tes_subset, d1 = 0, d2 = 0.4) #X# distance-based neighbor 
  # coo <- st_centroid(tes_subset)  # inv dist weight
  # dists <- nbdists(tes_nb, coo)
  # ids <- lapply(dists, function(x){1/x}) 
  # nbm <- nb2mat(tes_nb, glist = ids, style = "B")
  
  # Binary weighting assigns a weight of 1 to all neighboring features 
  # and a weight of 0 to all other features
  tes_w_binary <- nb2listw(tes_nb, style="B")
  # use county popu size as weights?
  # tes_w_size <- nb2listw(tes_nb, glist=tes_subset$pop, style="W")
  
  # adjacent matrix
  adj_matrix <- nb2mat(tes_nb, style="B")
  rownames(adj_matrix) <- colnames(adj_matrix) <- tes_subset$FIPS
  
  # spatial reg to impute those with suppressed mortality
  # random effects follow spatial correlation
  ddi0 = ddi[data.table(FIPS=tes_subset$FIPS),,on='FIPS']
  ddi0 = ddi0[suppressed!=T,]
  fitsp <- fitme(Deaths ~ pct_pov_2018+med_age+MFratio+
                   smoking+obese+foodenv+physicallyinactive+exdrinking+uninsured+
                   PCPrate+mammography+college+unemployed+avgdailyPM2.5+
                   adjacency(1|FIPS) + offset(log(pop)), 
                 adjMatrix = adj_matrix, data=ddi0, family='poisson')  
  res[[ca]]$fitsp = fitsp # summary(fitsp)$beta_table
  
  # calculate pearson residual-based pseudo R2
  eta =  cbind(1, as.matrix(ddi0[,c(4,6:18)])) %*% fitsp$fixef 
  fitsp0 <- fitme(Deaths ~ 1+ adjacency(1|FIPS) + offset(log(pop)), 
                  adjMatrix = adj_matrix, data=ddi0, family='poisson') 
  eta0 = rep(fitsp0$fixef, nrow(ddi0))
  res[[ca]]$pseudoR2_pearson = ddi0[,1-sum((Deaths-pop*exp(eta))^2/pop)/sum((Deaths-pop*exp(eta0))^2/pop)]
  
  # predict/impute Deaths counts of suppressed counties
  lphat = cbind(1, as.matrix(ddi[suppressed==T,c(4,6:18)])) %*% fitsp$fixef 
  Deaths.imp = ddi[suppressed==T,]$pop * exp(lphat)
  # how many imputed within suppression limit 10?
  res[[ca]]$imp.within.suppression = mean(Deaths.imp<10)  # 0.8600746
  Mortality.imp = pmax(pmin(Deaths.imp,9),1)/5/ddi[suppressed==T,]$pop*100000 
  ddi[suppressed==T,]$death_rate = Mortality.imp 
  dd[Cancer==ca&suppressed==T,death_rate:=Mortality.imp] 
  tes_subset$Mortality = ddi[,c('FIPS','death_rate')
  ][data.table(FIPS=tes_subset$FIPS),,on='FIPS']$death_rate
  
  # age-adjusted to 2000 popu! use med_age coef to calculate a ratio
  Deaths.adj.ratio = exp((ddi$med_age_2000-ddi$med_age) * fitsp$fixef['med_age']) 
  ddi$death_rate_age_adj_2000 = ddi$death_rate * Deaths.adj.ratio
  # ddi[,plot(death_rate, death_rate_age_adj_2000)]; abline(0,1)
  dd[Cancer==ca,death_rate_age_adj_2000:=death_rate * Deaths.adj.ratio]  
  tes_subset$Mortality_age_adj_2000 = ddi[,c('FIPS','death_rate_age_adj_2000')
  ][data.table(FIPS=tes_subset$FIPS),,on='FIPS']$death_rate_age_adj_2000
  
  # make a copy for convenience
  tes_subset$Mortality_2020 = tes_subset$Mortality 
  
  
  ### hotspot analysis, not age-adjusted, suppressed not shown
  tes_subset$Mortality = tes_subset$Mortality_2020
  # Calculate spatial lag 
  tes_lag <- lag.listw(tes_w_binary, tes_subset$Mortality)
  # Test for global G statistic of mortality
  res[[ca]]$Gtest = globalG.test(tes_subset$Mortality, tes_w_binary)
  # Identify neighbors, create weights, calculate spatial lag
  tes_nbs <- tes_subset |> 
    mutate(
      nb = st_contiguity(geometry),        # neighbors share border/vertex
      wt = st_weights(nb),                 # row-standardized weights
      tes_lag = st_lag(Mortality, nb, wt)  # calculate spatial lag of mortality
    ) 
  # Calculate the Gi using local_g_perm
  tes_hot_spots <- tes_nbs |> 
    mutate(
      Gi = local_g_perm(Mortality, nb, wt, nsim = 999)
      # nsim = number of Monte Carlo simulations (999 is default)
    ) |> 
    # The new 'Gi' column itself contains a dataframe, need to 'unnest' it
    unnest(Gi) 
  # Create a new data frame called 'tes_hot_spots"
  tes_hot_spots =
    tes_hot_spots |> #
    # with the columns 'gi' and 'p_folded_sim"
    # 'p_folded_sim' is the p-value of a folded permutation test
    # select(gi, p_folded_sim) |> 
    mutate(
      # Add a new column called "classification"
      classification = case_when(
        # Classify based on the following criteria:
        gi > 0 & p_folded_sim <= 0.01 ~ "Hot spot p<0.01",
        gi > 0 & p_folded_sim <= 0.05 ~ "Hot spot p<0.05",
        gi < 0 & p_folded_sim <= 0.01 ~ "Cold spot p<0.01",
        gi < 0 & p_folded_sim <= 0.05 ~ "Cold spot p<0.05",
        TRUE ~ "Not significant"
      )
    ) 
  ## attach back for plot with suppressed counties 
  cb$suppressed[cb$GEOID%in%tes_hot_spots$GEOID]= tes_hot_spots$suppressed
  cb$label_G_O[cb$GEOID%in%tes_hot_spots$GEOID] = tes_hot_spots$classification
  cb$label_G_O[is.na(cb$label_G_O) | cb$suppressed==T] = 'Suppressed'
  # attach label_G_O back
  # table(cb$label_G_O,useNA='ifany') 
  dd[Cancer==ca,label_G_O:=data.table(FIPS=cb$FIPS, label_G_O=cb$label_G_O)[ddi[,'FIPS'],,on='FIPS']$label_G_O] 
  # Convert hotspot clusters into a factor for easier plocbiing
  cb$label_G_O = factor(cb$label_G_O,
                        levels = c("Hot spot p<0.01","Hot spot p<0.05",
                                   "Not significant",
                                   "Cold spot p<0.05","Cold spot p<0.01",
                                   "Suppressed"))
  
  
  #S1# also plot suppressed...
  cb$label_G_O_imp[cb$GEOID%in%tes_hot_spots$GEOID] = tes_hot_spots$classification
  cb$label_G_O_imp[is.na(cb$label_G_O_imp)] = 'Not significant' # counties with empty neighbor
  # table(cb$label_G_O_imp,useNA='ifany') 
  dd[Cancer==ca,label_G_O_imp:=data.table(FIPS=cb$FIPS, label_G_O=cb$label_G_O_imp)[ddi[,'FIPS'],,on='FIPS']$label_G_O] 
  cb$label_G_O_imp = factor(cb$label_G_O_imp,
                            levels = c("Hot spot p<0.01","Hot spot p<0.05",
                                       "Not significant",
                                       "Cold spot p<0.05","Cold spot p<0.01")) 
  
  #S2# hotspot analysis using mortality rate, age-adjusted to 2000 popu
  tes_subset$Mortality = tes_subset$Mortality_age_adj_2000
  # Calculate spatial lag 
  tes_lag <- lag.listw(tes_w_binary, tes_subset$Mortality)
  # Test for global G statistic of mortality
  res[[ca]]$Gtest_age_adj_2000 = globalG.test(tes_subset$Mortality, tes_w_binary)
  # Identify neighbors, create weights, calculate spatial lag
  tes_nbs <- tes_subset |> 
    mutate(
      nb = st_contiguity(geometry),        # neighbors share border/vertex
      wt = st_weights(nb),                 # row-standardized weights
      tes_lag = st_lag(Mortality, nb, wt)  # calculate spatial lag of mortality
    ) 
  # Calculate the Gi using local_g_perm
  tes_hot_spots <- tes_nbs |> 
    mutate(
      Gi = local_g_perm(Mortality, nb, wt, nsim = 999)
      # nsim = number of Monte Carlo simulations (999 is default)
    ) |> 
    # The new 'Gi' column itself contains a dataframe, need to 'unnest' it
    unnest(Gi) 
  # Create a new data frame called 'tes_hot_spots"
  tes_hot_spots =
    tes_hot_spots |> #
    # with the columns 'gi' and 'p_folded_sim"
    # 'p_folded_sim' is the p-value of a folded permutation test
    # select(gi, p_folded_sim) |> 
    mutate(
      # Add a new column called "classification"
      classification = case_when(
        # Classify based on the following criteria:
        gi > 0 & p_folded_sim <= 0.01 ~ "Hot spot p<0.01",
        gi > 0 & p_folded_sim <= 0.05 ~ "Hot spot p<0.05",
        gi < 0 & p_folded_sim <= 0.01 ~ "Cold spot p<0.01",
        gi < 0 & p_folded_sim <= 0.05 ~ "Cold spot p<0.05",
        TRUE ~ "Not significant"
      )
    ) 
  ## attach back for plot with suppressed counties
  cb$label_G_O_age_adj_2000[cb$GEOID%in%tes_hot_spots$GEOID] = tes_hot_spots$classification
  cb$label_G_O_age_adj_2000[is.na(cb$label_G_O_age_adj_2000) | cb$suppressed==T] = 'Suppressed'
  # table(cb$label_G_O_age_adj_2000,useNA='ifany') 
  dd[Cancer==ca,label_G_O_age_adj_2000:=
       data.table(FIPS=cb$FIPS, label_G_O=cb$label_G_O_age_adj_2000)[ddi[,'FIPS'],,on='FIPS']$label_G_O] 
  cb$label_G_O_age_adj_2000 = factor(cb$label_G_O_age_adj_2000,
                                     levels = c("Hot spot p<0.01","Hot spot p<0.05",
                                                "Not significant",
                                                "Cold spot p<0.05","Cold spot p<0.01",
                                                "Suppressed"))
  
  
  #S3# sensitivity analysis using Rook contiguity in hotspot analysis
  tes_subset = cbind(cb, cbi[,-1])
  list_nb <- poly2nb(tes_subset, queen = F)
  empty_nb <- which(card(list_nb) == 0)
  tes_subset <- tes_subset[-empty_nb, ]
  empty_polygons <- tes_subset[empty_nb, ] 
  tes_nb <- poly2nb(tes_subset, queen = F) # Rook contiguity (only edge touching)
  tes_w_binary <- nb2listw(tes_nb, style="B")
  adj_matrix <- nb2mat(tes_nb, style="B")
  rownames(adj_matrix) <- colnames(adj_matrix) <- tes_subset$FIPS
  fitsp <- fitme(Deaths ~ pct_pov_2018+med_age+MFratio+
                   smoking+obese+foodenv+physicallyinactive+exdrinking+uninsured+
                   PCPrate+mammography+college+unemployed+avgdailyPM2.5+
                   adjacency(1|FIPS) + offset(log(pop)), 
                 adjMatrix = adj_matrix, data=tes_subset, family='poisson') 
  res[[ca]]$fitsp_Rook = fitsp # summary(fitsp)$beta_table
  lphat = cbind(1, as.matrix(ddi[suppressed==T,c(4,6:18)])) %*% fitsp$fixef 
  Deaths.imp = ddi[suppressed==T,]$pop * exp(lphat)
  Mortality.imp = pmax(pmin(Deaths.imp,9),1)/5/ddi[suppressed==T,]$pop*100000 
  ddi[suppressed==T,]$death_rate = Mortality.imp 
  dd[Cancer==ca&suppressed==T,death_rate_Rook:=Mortality.imp] 
  tes_subset$Mortality = ddi[,c('FIPS','death_rate')
  ][data.table(FIPS=tes_subset$FIPS),,on='FIPS']$death_rate
  tes_nbs <- tes_subset |> 
    mutate(
      nb = st_contiguity(geometry,queen=F), # neighbors share border 
      wt = st_weights(nb, style = "W"),    # row-standardized weights
      tes_lag = st_lag(Mortality, nb, wt)  # calculate spatial lag of mortality
    ) 
  tes_hot_spots <- tes_nbs |> 
    mutate(Gi = local_g_perm(Mortality, nb, wt, nsim = 999)) |> 
    unnest(Gi) 
  tes_hot_spots =
    tes_hot_spots |> 
    mutate(
      classification = case_when(
        gi > 0 & p_folded_sim <= 0.01 ~ "Hot spot p<0.01",
        gi > 0 & p_folded_sim <= 0.05 ~ "Hot spot p<0.05",
        gi < 0 & p_folded_sim <= 0.01 ~ "Cold spot p<0.01",
        gi < 0 & p_folded_sim <= 0.05 ~ "Cold spot p<0.05",
        TRUE ~ "Not significant"
      )
    )  
  cb$label_G_O_Rook[cb$GEOID%in%tes_hot_spots$GEOID] = tes_hot_spots$classification
  cb$label_G_O_Rook[is.na(cb$label_G_O_Rook) | cb$suppressed==T] = 'Suppressed'
  dd[Cancer==ca,label_G_O_Rook:=data.table(FIPS=cb$FIPS, label_G_O=cb$label_G_O_Rook)[ddi[,'FIPS'],,on='FIPS']$label_G_O] 
  cb$label_G_O_Rook = factor(cb$label_G_O_Rook,
                             levels = c("Hot spot p<0.01","Hot spot p<0.05",
                                        "Not significant",
                                        "Cold spot p<0.05","Cold spot p<0.01",
                                        "Suppressed"))
  # table(cb$label_G_O, cb$label_G_O_Rook)
  
  res[[ca]]$cb_to_plot = cb
}
#################### END: hotspot analysis ###################################


####################### Figures & Tables ######################################
library(usmap) 
# Map data for states and counties
states_data <- us_map("states") %>%filter(!abbr%in% c('AK','HI')) 
centroid_labels = usmapdata::centroid_labels("states")%>%filter(!abbr%in% c('AK','HI'))


### Figure 1: hotspot maps, US mainland
# pdf('hotspot_Wonder_2022_cancer4.pdf',width=14,height=10)
for(ii in 1:4){
  ca = cancer.names[ii]
  Ca = Cancer.names[ii]
  
  main.title=sprintf("Hotspots of %s cancer mortality (US mainland 2018-2022)",Ca)
  p0 = plot_usmap("counties", color = 'NA') +
    geom_sf(data = res[[ca]]$cb_to_plot,
            mapping = aes(geometry = geometry, fill = label_G_O),
            color = NA) +
    geom_sf(data = states_data,
            mapping = aes(geometry = geom),
            color = 'black',
            linewidth = 0.4,
            fill = NA) +
    geom_sf_text(data = centroid_labels,
                 mapping = aes(geometry = geom, label = abbr),
                 fontface = 'bold',
                 color = 'gray12',
                 size = 3) + 
    labs(fill = "Significant clusters", 
         title = main.title
    ) + scale_fill_manual(values=fc) +
    theme(legend.position.inside = c(0.04,0.2), 
          legend.title = element_text(size=16),
          legend.text = element_text(size=14),
          plot.title = element_text(size=18,hjust=0.4))+
    guides(fill=guide_legend(override.aes=list(colour="black",size=3)))
  print(p0)
}
dev.off()


### Figure S2: hotspot maps, US mainland, sensitivity analyses
# pdf('Figure S2 (sensitivity analyses).pdf',width=14,height=10)
for(ii in 1:4){
  ca = cancer.names[ii]
  Ca = Cancer.names[ii] 
  main.title=sprintf("Hotspots of %s cancer mortality (US mainland 2018-2022, age-adjust to 2000 population)",Ca)
  p0 = plot_usmap("counties", color = 'NA') +
    geom_sf(data = res[[ca]]$cb_to_plot,
            mapping = aes(geometry = geometry, fill = label_G_O_age_adj_2000),
            color = NA) +
    geom_sf(data = states_data,
            mapping = aes(geometry = geom),
            color = 'black',
            linewidth = 0.4,
            fill = NA) + 
    geom_sf_text(data = centroid_labels,
                 mapping = aes(geometry = geom, label = abbr),
                 fontface = 'bold',
                 color = 'gray12',
                 size = 3) + 
    labs(fill = "Significant clusters", 
         title = main.title
    ) + scale_fill_manual(values=fc) + 
    theme(legend.position.inside = c(0.04,0.2), 
          legend.title = element_text(size=16),
          legend.text = element_text(size=14),
          plot.title = element_text(size=18,hjust=0.4))+
    guides(fill=guide_legend(override.aes=list(colour="black",size=3)))
  print(p0)
  
  main.title=sprintf("Hotspots of %s cancer mortality (US mainland 2018-2022, include imputed data)",Ca)
  p0 = plot_usmap("counties", color = 'NA') +
    geom_sf(data = res[[ca]]$cb_to_plot,
            mapping = aes(geometry = geometry, fill = label_G_O_imp),
            color = NA) +
    geom_sf(data = states_data,
            mapping = aes(geometry = geom),
            color = 'black',
            linewidth = 0.4,
            fill = NA) + 
    geom_sf_text(data = centroid_labels,
                 mapping = aes(geometry = geom, label = abbr),
                 fontface = 'bold',
                 color = 'gray12',
                 size = 3) + 
    labs(fill = "Significant clusters", 
         title = main.title
    ) + scale_fill_manual(values=fc[1:5]) + 
    theme(legend.position.inside = c(0.04,0.2), 
          legend.title = element_text(size=16),
          legend.text = element_text(size=14),
          plot.title = element_text(size=18,hjust=0.4))+
    guides(fill=guide_legend(override.aes=list(colour="black",size=3)))
  print(p0)
  
  main.title=sprintf("Hotspots of %s cancer mortality (US mainland 2018-2022, use Rook contiguity neighbor)",Ca)
  p0 = plot_usmap("counties", color = 'NA') +
    geom_sf(data = res[[ca]]$cb_to_plot,
            mapping = aes(geometry = geometry, fill = label_G_O_Rook),
            color = NA) +
    geom_sf(data = states_data,
            mapping = aes(geometry = geom),
            color = 'black',
            linewidth = 0.4,
            fill = NA) + 
    geom_sf_text(data = centroid_labels,
                 mapping = aes(geometry = geom, label = abbr),
                 fontface = 'bold',
                 color = 'gray12',
                 size = 3) + 
    labs(fill = "Significant clusters", 
         title = main.title
    ) + scale_fill_manual(values=fc) + 
    theme(legend.position.inside = c(0.04,0.2), 
          legend.title = element_text(size=16),
          legend.text = element_text(size=14),
          plot.title = element_text(size=18,hjust=0.4))+
    guides(fill=guide_legend(override.aes=list(colour="black",size=3)))
  print(p0)
}
dev.off()


### Figure S3: direct mapping of mortality rates
# pdf('Figure S3 (direct mapping).pdf',width=14,height=10)
for(ii in 1:4){
  ca = cancer.names[ii]
  Ca = Cancer.names[ii]
  cb = res[[ca]]$cb_to_plot
  tt = dd[Cancer==ca,][data.table(FIPS=cb$FIPS),,on='FIPS']
  tt[suppressed==T,death_rate:=NA]
  cb$death_rate = tt$death_rate
  breaks = round(c(-0.01, quantile(cb$death_rate, seq(0.2,1,by=0.2), na.rm=T)+0.1),1)
  # breaks
  labels = c('Very low','Low','Medium','High','Very high')
  cb$label_5 = cut(cb$death_rate, breaks, include.lowest = T) 
  levels = paste(labels, levels(cb$label_5))
  cb$label_5 = levels[as.numeric(cb$label_5)]
  cb$label_5[is.na(cb$label_5)] = 'Suppressed'
  cb$label_5 = factor(cb$label_5, levels = c(rev(levels),'Suppressed'))
  table(cb$label_5, useNA='ifany')
  main.title=sprintf("Direct visualization of %s cancer mortality (US mainland 2018-2022)",Ca)
  p0 = plot_usmap("counties", color = 'NA') +
    geom_sf(data = cb,
            mapping = aes(geometry = geometry, fill = label_5),
            color = NA) +
    geom_sf(data = states_data,
            mapping = aes(geometry = geom),
            color = 'black',
            linewidth = 0.4,
            fill = NA) + 
    geom_sf_text(data = centroid_labels,
                 mapping = aes(geometry = geom, label = abbr),
                 fontface = 'bold',
                 color = 'gray12',
                 size = 3) + 
    labs(fill = "Mortality (per 100k)", 
         title = main.title
    ) + scale_fill_manual(values=fc) +
    # scale_fill_gradientn(n.breaks=5, colors=rev(fc[1:5]),na.value = fc[6] ) +
    # scale_fill_gradient2(low="blue",mid="lightgrey",high="red",
    #                      na.value="white", midpoint=median(cb$death_rate,na.rm=T)) +
    theme(legend.position.inside = c(0.04,0.2), 
          legend.title = element_text(size=16),
          legend.text = element_text(size=14),
          plot.title = element_text(size=18,hjust=0.4))+
    guides(fill=guide_legend(override.aes=list(colour="black",size=3)))
  print(p0)
}
dev.off()


### % counties have different hotspot clustering using Rook contiguity neighbor definition
dd[,.(n=sum(substr(label_G_O_Rook,1,8)!=substr(label_G_O,1,8)), 
      pd=mean(substr(label_G_O_Rook,1,8)!=substr(label_G_O,1,8))),by=Cancer]
# 1:     Breast    60 0.01930502
# 2: Colorectal    91 0.02927928
# 3:       Lung    69 0.02220077
# 4:   Prostate    56 0.01801802
# no counties are changed from hotspot to cold spot, or coldspot to hotspot...
# dd[Cancer==cancer.names[4],table(substr(label_G_O_Rook,1,8),substr(label_G_O,1,8))]


### Figure 2: avg mortality for each hotspot strata  
dd[,State:=substr(County, nchar(County)-2, nchar(County))]
dd[,hotspot_2020:=grepl('Hot spot',label_G_O)]
dd[,hotspot_2000:=grepl('Hot spot',label_G_O_age_adj_2000)]
dd$label_G_O = factor(dd$label_G_O,
                      levels = c("Hot spot p<0.01","Hot spot p<0.05",
                                 "Not significant",
                                 "Cold spot p<0.05","Cold spot p<0.01",
                                 "Suppressed"))
dd[,label_3:=label_G_O]
dd[grepl('Hot spot',label_G_O),label_3:='Hot spot']
dd[grepl('Cold spot',label_G_O),label_3:='Cold spot']
 
mort.strat = dd[,.(nc=.N,pop=sum(pop),
                   death_rate=sum(death_rate*pop)/sum(pop)),
                by=.(Cancer,label_3)][order(Cancer,label_3),]
names(mort.strat)[2] = 'label_G_O'
mort.strat = mort.strat[label_G_O!='Suppressed']
mort.strat[,label_G_O:=str_wrap(label_G_O, 4)]
mort.strat[,label_G_O:=factor(label_G_O, levels =c("Hot\nspot", "Not\nsignificant", "Cold\nspot"))]
mort.strat[,death_rate_sd:= sqrt(death_rate*(100000-death_rate)/pop)]
mort.strat[Cancer=='Colorectal',death_rate_overall:=16.1]
mort.strat[Cancer=='Lung', death_rate_overall:=41.2]
mort.strat[Cancer=='Breast', death_rate_overall:=25.9]
mort.strat[Cancer=='Prostate',death_rate_overall:=20.2]

mort.strat[,Cancer:=paste0(Cancer,' cancer')]
mort.strat[,Cancer:=factor(Cancer, levels = c("Colorectal cancer", "Lung cancer","Breast cancer","Prostate cancer"))]
# pdf('mortality_rate_by_hotspot_strata3.pdf',width=12,height=5)
ggplot(mort.strat, aes(x=label_G_O, y=death_rate, # shape=label_G_O,
                       color=label_G_O, group=Cancer)) +
  geom_errorbar(aes(ymin= death_rate-1.96* death_rate_sd, 
                    ymax= death_rate+1.96* death_rate_sd), 
                linewidth=1, width=0.2) +
  # geom_line() +
  geom_hline(aes(yintercept= death_rate_overall), linetype = "dashed") +
  geom_point(size=3) + 
  facet_wrap(. ~ Cancer, scales="free", ncol=4)+
  labs(title= '', x =  '', y = "Cancer mortality, per 100,000" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)
        , axis.ticks.x = element_blank()
        , strip.text.x = element_text(size = 18)
        , panel.grid.major = element_blank()
        , axis.title.y = element_text(size=20)
        , axis.text.x = element_text(size=16)  
        , axis.text.y = element_text(size=16) 
        , legend.position = "none" 
  )+
  scale_color_manual(breaks =levels(mort.strat$label_G_O), values=fc[c(1,3,5)]) # fc[c(1:5)] # 
dev.off()



### Table S1: spatial error regression output
for(ii in 1:4){
  ca = cancer.names[ii]
  tt = data.frame(summary(res[[ca]]$fitsp)$beta_table )
  tt$pv = round(pnorm(-abs(tt$t.value))*2, 3)
  tt[,1:2] = round(tt[,1:2]*100,2)
  tt = tt[,-3] 
  write.table(tt, file='spatial_error_reg.csv',sep=',',append = T)
}


### Table S2: all identified hotspot counties
write.csv(dd[grepl('Hot spot',label_G_O),c(1,19,20,21,24,25,28)][order(Cancer)], 
          file='hotspot_counties.csv')



## Figure S1 (high-res of Fig 1): hotspot maps, US mainland, with NCI cancer centers
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9074106/
# https://github.com/idblr/NCI_Cancer_Center_Catchment_Areas
# https://gis.cancer.gov/ncicatchment/
ncicc <- st_read("NCI_CancerCenter_Address_fall2023/NCI_CancerCenter_Address_fall2023.shp")
ncicc = ncicc[ncicc$state!='Hawaii',]
ncicc = ncicc[!grepl('St. Jude Children',ncicc$name),] # 63
# ncicc[ncicc$state=='Florida',]

# pdf('Figure S1 (hotspots with NCI cancer centers).pdf',width=14,height=10 )
for(ii in 1:4){
  ca = cancer.names[ii]
  Ca = Cancer.names[ii]
  
  main.title=sprintf("Hotspots of %s cancer mortality (US mainland 2018-2022)",Ca)
  p0 = plot_usmap("counties", color = 'NA') +
    geom_sf(data = res[[ca]]$cb_to_plot,
            mapping = aes(geometry = geometry, fill = label_G_O),
            color = NA) +
    geom_sf(data = states_data,
            mapping = aes(geometry = geom),
            color = 'black',
            linewidth = 0.4,
            fill = NA) +
    geom_sf(data = ncicc,
            mapping = aes(geometry = geometry ),
            shape = 16, col='black',
            size = 3) + 
    geom_sf_text(data = centroid_labels,
                 mapping = aes(geometry = geom, label = abbr),
                 fontface = 'bold',
                 color = 'gray12',
                 size = 3)+
    labs(fill = "Significant clusters", 
         title = main.title
    ) + scale_fill_manual(values=fc) + 
    theme(legend.position.inside = c(0.04,0.2), 
          legend.title = element_text(size=16),
          legend.text = element_text(size=14),
          plot.title = element_text(size=18,hjust=0.4))+
    guides(fill=guide_legend(override.aes=list(colour="black",size=3)))
  print(p0)
}
dev.off()

##################### END: Figures & Tables  ###################################
