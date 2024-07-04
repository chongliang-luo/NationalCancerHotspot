#### Luo et al. Where should the cancer control interventions target: 
#### A geospatial hotspot analysis for major cancer mortality 2016-2020 in the U.S.

### hotspot analysis: use Getis-Ord G statistics 
require(sf) 
require(sfdep)
require(spdep)
require(dplyr)
require(tidyr)
require(ggplot2)
require(data.table)
require(stringr)
  
## US county bound shape data
# https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
cb <- st_read("dat/cb_2018_us_county_500k/cb_2018_us_county_500k.shp")
dim(cb) # 3233   10  

## check state fips,
# https://www.fhwa.dot.gov/ohim/hpmsmanl/pdf/appa.pdf
# 2=Alaska 15=hawaii 60=American Samoa 66=Guam 69=Northern Marianas 72=puerto rico 78=Virgin Islands
table(cb$STATEFP) # 01 - 78 
fips49 = setdiff(names(table(cb$STATEFP))[1:51], c('02','15')) # only mainland US  


### Figure 1: hotspot analysis, and plots
cancer.names = c('crc','lung','breast','prostate')
# color labels for hotspot clusters 
fc = c("red","orange","lightgrey","deepskyblue1","blue","white")
set.seed(123)
pdf('NationalHotspot/git/hotspot_2020_cancer4_G_O.pdf',width=14,height=10)
# some analysis code from 
# https://rpubs.com/heatherleeleary/hotspot_getisOrd_tut
for(ii in 1:4){
  ## attach cancer mortality rate
  ca = cancer.names[ii]
  Ca = c('colorectal','lung','breast','prostate')[ii]
  mm = dd[cancer==ca, ]
  tt = mm[data.table(FIPS=as.numeric(cb$GEOID)),,on='FIPS']
  names(tt)[c(17,18)] = c('County', 'Mortality')
  tt = cbind(cb, tt[,-1])
  tes_subset = tt[!is.na(tt$Mortality),]
  # tes_subset = tt[!is.na(tt$pop2020M),] 
  
  ## removed empty neighbor sets
  list_nb <- poly2nb(tes_subset, queen = TRUE)
  empty_nb <- which(card(list_nb) == 0)
  tes_subset <- tes_subset[-empty_nb, ]
  empty_polygons <- tes_subset[empty_nb, ]
  # empty_polygons$NAME  # print neighborhood names 
  # Identify neighbors with queen contiguity (edge/vertex touching)
  tes_nb <- poly2nb(tes_subset, queen = TRUE)
  
  # Binary weighting assigns a weight of 1 to all neighboring features 
  # and a weight of 0 to all other features
  tes_w_binary <- nb2listw(tes_nb, style="B")
  # use county popu size as weights?
  # tes_w_size <- nb2listw(tes_nb, glist=tes_subset$pop2020M, style="W")
  
  # Calculate spatial lag 
  tes_lag <- lag.listw(tes_w_binary, tes_subset$Mortality)
  
  # Test for global G statistic of TreeEquity
  globalG.test(tes_subset$Mortality, tes_w_binary)
  
  # Identify neighbors, create weights, calculate spatial lag
  tes_nbs <- tes_subset |> 
    mutate(
      nb = st_contiguity(geometry),        # neighbors share border/vertex
      wt = st_weights(nb),                 # row-standardized weights
      tes_lag = st_lag(Mortality, nb, wt)  # calculate spatial lag of TreeEquity
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
  cb$label_G_O=cb$is_missing=cb$is_missing=NULL
  cb$label_G_O[cb$GEOID%in%tes_hot_spots$GEOID] = tes_hot_spots$classification
  cb$is_missing[cb$GEOID%in%tes_hot_spots$GEOID]= tes_hot_spots$is_missing
  cb$label_G_O[is.na(cb$label_G_O) | cb$is_missing==T] = 'Suppressed'
  # table(cb$label_G_O,useNA='ifany') 
  tt = data.table(FIPS=as.numeric(cb$GEOID), 
                  label_G_O=cb$label_G_O)[mm[,'FIPS'],,on='FIPS']
  tt[is.na(label_G_O),label_G_O:='Suppressed'] 
  dd[cancer==ca, label_G_O:=tt$label_G_O]
  
  # Convert hotspot clusters into a factor for easier plotting
  cb$label_G_O = factor(cb$label_G_O,
                        levels = c("Hot spot p<0.01","Hot spot p<0.05",
                                   "Not significant",
                                   "Cold spot p<0.05","Cold spot p<0.01",
                                   "Suppressed"))
  
  # plot US mainland
  p0 = cb |>  filter(STATEFP%in%fips49) |>
    # Visualize the results with ggplot2
    ggplot(aes(fill = label_G_O)) +
    geom_sf(color = "black", lwd = 0.1) +
    # scale_fill_brewer(type = "div", palette = 6) +
    theme_void() +
    labs(
      fill = "Significant clusters",
      title = sprintf("        Hotspots of %s cancer mortality (US 2016-2020)",Ca)
    )+
    scale_fill_manual(values=fc)
  
  print(p0) 
}
# dev.off()



### Fig 2: calculate avg mortality for each hotspot strata 
mort.strat = rbind(
  dd[cancer=='lung',.(cancer='Lung cancer',pop=sum(pop2020F+pop2020M),age_adjusted_death_rate=sum(age_adjusted_death_rate*(pop2020F+pop2020M),na.rm=T)/sum(pop2020F+pop2020M)),by=label_G_O],
  dd[cancer=='crc',.(cancer='Colorectal cancer',pop=sum(pop2020F+pop2020M),age_adjusted_death_rate=sum(age_adjusted_death_rate*(pop2020F+pop2020M),na.rm=T)/sum(pop2020F+pop2020M)),by=label_G_O],
  dd[cancer=='breast',.(cancer='Breast cancer (female)',pop=sum(pop2020F),age_adjusted_death_rate=sum(age_adjusted_death_rate*pop2020F,na.rm=T)/sum(pop2020F)),by=label_G_O],
  dd[cancer=='prostate',.(cancer='Prostate cancer (male)',pop=sum(pop2020M),age_adjusted_death_rate=sum(age_adjusted_death_rate*pop2020M,na.rm=T)/sum(pop2020M)),by=label_G_O]
)

mort.strat[,age_adjusted_death_rate_sd:= sqrt(age_adjusted_death_rate*(100000-age_adjusted_death_rate)/pop)]
mort.strat = mort.strat[label_G_O!='Suppressed']
mort.strat[,label_G_O:=str_wrap(label_G_O, 9)]
mort.strat[,label_G_O:=factor(label_G_O, levels = c("Hot spot\np<0.01", "Hot spot\np<0.05","Not\nsignificant", "Cold spot\np<0.05","Cold spot\np<0.01"))]
mort.strat[cancer=='Colorectal cancer',age_adjusted_death_rate_overall:=13.1]
mort.strat[cancer=='Lung cancer',age_adjusted_death_rate_overall:=35.0]
mort.strat[cancer=='Breast cancer (female)',age_adjusted_death_rate_overall:=19.6]
mort.strat[cancer=='Prostate cancer (male)',age_adjusted_death_rate_overall:=18.8]
mort.strat[,cancer:=factor(cancer, levels = c("Colorectal cancer", "Lung cancer","Breast cancer (female)","Prostate cancer (male)"))]

# pdf('NationalHotspot/mortality_rate_by_hotspot_strata.pdf',width=12,height=10)
ggplot(mort.strat, aes(x=label_G_O, y=age_adjusted_death_rate, # shape=label_G_O,
                       color=label_G_O, group=cancer)) +
  geom_errorbar(aes(ymin=age_adjusted_death_rate-1.96*age_adjusted_death_rate_sd, 
                    ymax=age_adjusted_death_rate+1.96*age_adjusted_death_rate_sd), 
                linewidth=1, width=0.2) +
  # geom_line() +
  geom_hline(aes(yintercept=age_adjusted_death_rate_overall), linetype = "dashed") +
  geom_point(size=3) + 
  facet_wrap(. ~ cancer, scales="free", ncol=2)+
  labs(title= '', x =  '', y = "Age adjusted death rate, per 100,000" ) +
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
  scale_color_manual(breaks =levels(mort.strat$label_G_O), values=fc[1:5])
# dev.off()




