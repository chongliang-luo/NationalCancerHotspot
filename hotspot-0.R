#### Luo et al. Where should the cancer control interventions target: 
#### A geospatial hotspot analysis for major cancer mortality 2016-2020 in the U.S.

## data clean, missing mortality imputation
require(data.table)
dd = readRDS(file='SCP2020_mortality.RDS')

###  data dictionary  
# [1] "FIPS"         # county FIPS           
# [2] "pop2020M"     # 2020 census population male (R package tidycensus)
# [3] "pop2020F"     # 2020 census population female      
## 2014 County Health Ranking cancer risk factors
## https://www.countyhealthrankings.org/health-data/methodology-and-sources/data-documentation
# [4] "pct_pov_2014"           
# [5] "smoking"                
# [6] "obese"                  
# [7] "foodenv"                
# [8] "physicallyinactive"     
# [9] "exdrinking"             
# [10] "uninsured"              
# [11] "PCPrate"                
# [12] "avgdailyPM2.5"          
# [13] "state_id"     # state ID          
# [14] "lat"          # county geo latitude          
# [15] "lng"          # county geo longitude         
# [16] "cancer"       # cancer type          
# [17] "State"        # county name   
## SCP 2016-2020 age-adjusted cancer mortality rate, 
## https://statecancerprofiles.cancer.gov/deathrates/index.php
# [18] "age_adjusted_death_rate" 
# [19] "avg_annual_count"       
# [20] "is_missing"             
# [21] "mort.mean"    # national average mortality
 

### Table 2: regress mortality rate on CHR risk factors to impute the missing mortality rates 
uni_multi_reg <- function(X, y, weights,family='gaussian', uni_sel_p=0.05){
  p = ncol(X)
  tt = matrix(NA,0,4)
  for(i in 1:p){
    fit = lm(y~X[,i],weights = weights)
    tt = rbind(tt, summary(fit)$coef[-1,])
  }
  idx = which(tt[,4]<uni_sel_p)
  fit = lm(y~as.matrix(X[,idx]),weights=weights)
  print(summary(fit)$r.squared)
  tt = rbind(tt, summary(fit)$coef[-1,])
  tt = round(tt,3)
  row.names(tt) = c(names(X), names(X)[idx])
  tt 
}

cancer.names = c('crc', 'lung', 'breast', 'prostate')
# ii = 'prostate'
fit4 = list()
for(ii in cancer.names){ 
  cat(ii, '...')
  ddi = dd[cancer==ii, ]
  ddi$pop = ddi$pop2020F + ddi$pop2020M
  if(ii=='breast') ddi$pop = ddi$pop2020F
  if(ii=='prostate') ddi$pop = ddi$pop2020M
  
  ddi[is_missing!=T,age_adjusted_death_count:=round(age_adjusted_death_rate*pop/100000)] 
  # counties with no missing: risk factors (pct_pov_2014 - avgdailyPM2.5), 
  fips.no.mis=ddi[apply(ddi[,c(4:12),with=F],1,function(a) !any(is.na(a)))==T]$FIPS  
  #  and outcome age_adjusted_death_rate 
  fips.no.mis2=ddi[apply(ddi[,c(4:12,23),with=F],1,function(a) !any(is.na(a)))==T]$FIPS  
  ddi = ddi[FIPS%in%fips.no.mis]    # 3141 -> 3134
  ddi2 = ddi[FIPS%in%fips.no.mis2]  # 2174 2750 1755 1576
  
  # popu weighted linear reg, R2 = 44.1, 70.8, 27.4, 15.8%
  fit4[[ii]] = uni_multi_reg(data.frame(ddi2[,c(4:12)]), ddi2$age_adjusted_death_rate, weights=ddi2$pop)
  fitl = ddi2[,lm(age_adjusted_death_rate~pct_pov_2014+smoking+obese+foodenv+physicallyinactive+exdrinking+uninsured+PCPrate+avgdailyPM2.5,weights = pop)]
  rh = predict(fitl,newdata=ddi[!FIPS%in%fips.no.mis2,], type = 'response') # predicted rate
  
  # ad-hoc cap age_adjusted_death_count at 5 (unadjusted avg_annual_count capped at 3)
  ddi[!FIPS%in%fips.no.mis2, age_adjusted_death_rate:=pmin(rh, 5/pop*100000)] 
  dd[cancer==ii&FIPS%in%fips.no.mis, age_adjusted_death_rate:=ddi$age_adjusted_death_rate]
}

# saveRDS(dd, file='NationalHotspot/git/SCP2020_mortality.RDS')




