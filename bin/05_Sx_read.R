#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 05- Processing of Sx Trial Data #


######## Process Sx trial data ########
######## Version 2.0- ht_pct added ####
######## Liam Gilson ##################
### Oct. 29, 2024 #####################


## SX data ingest and cleaning ##

library(dplyr)

############ Tree Data ###############

sx_trees<-read.csv("./Data/Sx_Trial/Sx_Growth_v21.csv",header = T)

# NAs are represented by ".", need to standardize

is.na(sx_trees)<-sx_trees == "."

# we have a bunch of numerics classified as chr because of these periods
sx_trees$HT3<-as.numeric(sx_trees$HT3)
sx_trees$HT6<-as.numeric(sx_trees$HT6)
sx_trees$HT10<-as.numeric(sx_trees$HT10)
sx_trees$HT16<-as.numeric(sx_trees$HT16)

# data also used "0" code for mortality- probably need to replace all zeros in the ht columns with NA

sx_trees$HT3[which(sx_trees$HT3==0)]<-NA
sx_trees$HT6[which(sx_trees$HT6==0)]<-NA
sx_trees$HT10[which(sx_trees$HT10==0)]<-NA
sx_trees$HT16[which(sx_trees$HT16==0)]<-NA
### Calculate increments ###


#Planted- spring 2005, measured Fall 2007,2010,2014,2020 (2022 partial measurements for some remote sensing project)
# 3-2007, 6,2010,10,2014,16,2020
# No increment from '05 to '07, we don't know seedling heights
sx_trees$Inc_2010<-sx_trees$HT6-sx_trees$HT3
sx_trees$Inc_2014<-sx_trees$HT10-sx_trees$HT6
sx_trees$Inc_2020<-sx_trees$HT16-sx_trees$HT10

# height percentiles
sx_trees <- sx_trees |> group_by(Site) |> mutate(pct_2007 = ecdf(HT3)(HT3),pct_2010 = ecdf(HT6)(HT6),pct_2014 = ecdf(HT10)(HT10),pct_2020 = ecdf(HT16)(HT16)) 
#sx_trees <- sx_trees |> group_by(Site) |> arrange(desc(DBH16)) |> mutate(BAL_ = cumsum(pi*(DIA/(2*12))^2))

# Need to collapse this into long format- one row per tree per increment

sx_clean<-sx_trees %>% tidyr::pivot_longer(cols=starts_with("Inc_"),names_to="year",values_to="dh",names_prefix="Inc_")

sx_clean$height.x <- case_when(
  sx_clean$year == 2010 ~ sx_clean$HT3,
  sx_clean$year == 2014 ~ sx_clean$HT6,
  sx_clean$year == 2020 ~ sx_clean$HT10
)

sx_clean$height.y <- case_when(
  sx_clean$year == 2010 ~ sx_clean$HT6,
  sx_clean$year == 2014 ~ sx_clean$HT10,
  sx_clean$year == 2020 ~ sx_clean$HT16
)

sx_clean$ht_pct.x <- case_when(
  sx_clean$year == 2010 ~ sx_clean$pct_2007,
  sx_clean$year == 2014 ~ sx_clean$pct_2010,
  sx_clean$year == 2020 ~ sx_clean$pct_2014
)

sx_clean$ht_pct.y <- case_when(
  sx_clean$year == 2010 ~ sx_clean$pct_2010,
  sx_clean$year == 2014 ~ sx_clean$pct_2014,
  sx_clean$year == 2020 ~ sx_clean$pct_2020
)


sx_trial <- sx_clean %>% select(Tree_ID,Site,Prov,year,height.x,height.y,ht_pct.x,ht_pct.y,dh)

sx_trial$dt<-case_when(sx_trial$year == 2010 ~ 3, sx_trial$year == 2014 ~ 4, sx_trial$year == 2020 ~ 6)

sx_trial$year<-as.numeric(sx_trial$year)

sx_trial$age <- sx_trial$year-2007


sx_trial$year.y <- sx_trial$year
sx_trial$year.x <- case_when(
  sx_trial$year == 2007 ~ 2005,
  sx_trial$year == 2010 ~ 2007,
  sx_trial$year == 2014 ~ 2010,
  sx_trial$year == 2020 ~ 2014
)

#saveRDS(sx_trial,file="./Data/Sx_Trial/sx_trial_v2.RDS")
# sx_provs<-read.csv("./Data/Sx_Trial/Sx_cna_prov_Normal_1961_1990SY.csv",header = T)
# sx_sites<-read.csv("./Data/Sx_Trial/Sx_cna_site_2007-2020SY.csv",header = T)
# 
# sx_ints <- intervals_cna()

