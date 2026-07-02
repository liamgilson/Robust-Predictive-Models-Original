#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 31- Process Prince George validation data #


######## Process PG trial data ########
######## Liam Gilson ##################


## pg data ingest and cleaning ##

library(dplyr)

############ Tree Data ###############

pg_trees<-read.csv("./Data/Pg_Trial/pg_trees.csv",header = T,
                   na.strings = ".")

# NAs are represented by ".", need to standardize
str(pg_trees)
# something slipped through in the ht5 and made it chr
# it's periods with trailing zeros...

is.na(pg_trees)<-pg_trees == ".    "

which.nonnum <- function(x) {
  badNum <- is.na(suppressWarnings(as.numeric(as.character(x))))
  which(badNum & !is.na(x))
}
# 
# badNum <- is.na(suppressWarnings(as.numeric(as.character(pg_trees$ht05))))
# which(badNum & !is.na(pg_trees$ht05))
# pg_trees$ht05[which(badNum & !is.na(pg_trees$ht05))]
is.na(pg_trees)<-pg_trees == ". "
pg_trees$ht05<-as.integer(pg_trees$ht05)
str(pg_trees)
#is.na(pg_trees)<-pg_trees == "."

pg_trees$Htfall03[which.nonnum(pg_trees$Htfall03)]
is.na(pg_trees$Htfall03)<-pg_trees$Htfall03 == ".  "
pg_trees$Htfall03<-as.integer(pg_trees$Htfall03)
str(pg_trees)
# all heights were cm

pg_trees$ht03<-pg_trees$Htfall03/100
pg_trees$ht05<-pg_trees$ht05/100
pg_trees$ht12<-pg_trees$ht12/100
pg_trees$ht17<-pg_trees$ht17/100
# pg_trees$ht05[which(badNum & !is.na(pg_trees$ht05))]


# so 999 in Htspring03 indicates "not planted", need to see what happens with mortality trees

# looking for dead code

test <- pg_trees |> filter(code05.1 == 7)
# remove unplanted
test <- pg_trees |> filter(ht03 == 9.99)

pg_trees <- pg_trees |> filter(ht03 != 9.99)



# okay it looks like dead trees had NA recorded for height, perfect



### Calculate increments ###


#Planted- spring 2005, measured Fall 2007,2010,2014,2020 (2022 partial measurements for some remote sensing project)
# 3-2007, 6,2010,10,2014,16,2020
# No increment from '05 to '07, we don't know seedling heights
pg_trees$Inc_2005<-pg_trees$ht05-pg_trees$ht03
pg_trees$Inc_2012<-pg_trees$ht12-pg_trees$ht05
pg_trees$Inc_2017<-pg_trees$ht17-pg_trees$ht12

# height percentiles
pg_trees <- pg_trees |> group_by(Site) |> mutate(pct_2003 = ecdf(ht03)(ht03),pct_2005 = ecdf(ht05)(ht05),pct_2012 = ecdf(ht12)(ht12),pct_2017 = ecdf(ht17)(ht17)) 
#pg_trees <- pg_trees |> group_by(Site) |> arrange(desc(DBH16)) |> mutate(BAL_ = cumsum(pi*(DIA/(2*12))^2))

# Need to collapse this into long format- one row per tree per increment

pg_clean<-pg_trees %>% tidyr::pivot_longer(cols=starts_with("Inc_"),names_to="year",values_to="dh",names_prefix="Inc_")

pg_clean$height.x <- case_when(
  pg_clean$year == 2005 ~ pg_clean$ht03,
  pg_clean$year == 2012 ~ pg_clean$ht05,
  pg_clean$year == 2017 ~ pg_clean$ht12
)

pg_clean$height.y <- case_when(
  pg_clean$year == 2005 ~ pg_clean$ht05,
  pg_clean$year == 2012 ~ pg_clean$ht12,
  pg_clean$year == 2017 ~ pg_clean$ht17
)

pg_clean$ht_pct.x <- case_when(
  pg_clean$year == 2005 ~ pg_clean$pct_2003,
  pg_clean$year == 2012 ~ pg_clean$pct_2005,
  pg_clean$year == 2017 ~ pg_clean$pct_2012
)

pg_clean$ht_pct.y <- case_when(
  pg_clean$year == 2005 ~ pg_clean$pct_2005,
  pg_clean$year == 2012 ~ pg_clean$pct_2012,
  pg_clean$year == 2017 ~ pg_clean$pct_2017
)

# I think at this point we should filter for only Class B stock

pg_clean <- pg_clean |> filter(Class == "B")

# record is tree ID, id is seedlot

pg_trial <- pg_clean %>% select(Record,Site,id,year,height.x,height.y,ht_pct.x,ht_pct.y,dh)
#pg_trial <- pg_clean %>% select(Tree_ID,Site,Prov,year,height.x,height.y,ht_pct.x,ht_pct.y,dh)


pg_trial$dt<-case_when(pg_trial$year == 2005 ~ 2, pg_trial$year == 2012 ~ 7, pg_trial$year == 2017 ~ 5)

pg_trial$year<-as.numeric(pg_trial$year)

pg_trial$age <- pg_trial$year-2003


pg_trial$year.y <- pg_trial$year
pg_trial$year.x <- case_when(
  pg_trial$year == 2005 ~ 2003,
  pg_trial$year == 2012 ~ 2005,
  pg_trial$year == 2017 ~ 2012
)

#saveRDS(pg_trial,file="./Data/pg_Trial/pg_trial.RDS")
# pg_provs<-read.csv("./Data/pg_Trial/pg_cna_prov_Normal_1961_1990SY.csv",header = T)
# pg_sites<-read.csv("./Data/pg_Trial/pg_cna_site_2007-2020SY.csv",header = T)
# 
# pg_ints <- intervals_cna()

#### Climate merge ####

#pg_trial <- readRDS(file="./Data/PG_Trial/pg_trial.RDS")

pg_provs<-read.csv("./Data/PG_Trial/prov_clim_Normal_1961_1990SY.csv",header = T)
pg_sites<-read.csv("./Data/PG_Trial/site_clim_2003-2017SY.csv",header = T)
pg_provs<-vpd_cna(pg_provs)
pg_sites<-vpd_cna(pg_sites)
# they used different abbreviations...
pg_sites$id1 <- case_when(
  pg_sites$id1 == "Church" ~ "Church",
  pg_sites$id1 == "Seebach" ~ "Seeb",
  pg_sites$id1 == "Hayden" ~ "Hayd"
)
pg_sites$id2 <- pg_sites$id1
pg_ints <- intervals_cna(pg_trial,pg_sites,"Site")

pg_clim <- inner_join(pg_trial,pg_ints,join_by(year.x,year.y,Site==site))

colnames(pg_provs)[1]<-"Prov"
pg_provs<-pg_provs[,-2]
colnames(pg_provs)[5:91]<-paste("P",colnames(pg_provs)[5:91],sep="_")
pg_clim <- pg_clim |> rename("Prov" = "id")
pg_clim <- inner_join(pg_clim,pg_provs,join_by(Prov))

#### format and transfer variables ####
pg_clim <- pg_clim |> rename("Tree_ID" = "Record")
pg_merge <- pg_clim |> select(Site,Tree_ID,Prov,Latitude,Longitude,Elevation,year.x,year.y,height.x,height.y,ht_pct.x,ht_pct.y,dt,dh,starts_with(c("mean_","min_","max_","P_")))
pg_merge$source <- "pg"
pg_merge$Type <- "Prov"
# check- heights may be cm
pg_merge$height.x[1:100]
max(pg_merge$height.x,na.rm=T)


mean_vars <- as.matrix(pg_merge |> ungroup() |> select(starts_with("mean_")))
p_vars <- as.matrix(pg_merge |> ungroup() |> select(starts_with("P_")))
difs <- mean_vars-p_vars
colnames(difs)<-gsub("mean_","T_",colnames(difs))
colnames(difs)

pg_merge <- cbind(pg_merge,difs)
rm(difs)

saveRDS(pg_merge,file="./Data/PG_Trial/pg_dat.RDS")

