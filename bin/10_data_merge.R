#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 10- Merging of all processed datasets #


################################################################################
# Final processing and merging of climate variables
################################################################################

# Dependencies:
# dis_functions.R
library(dplyr)
library(ggplot2)

#### FAIB ####

faib_cna_intervals<-read.csv(file="./Data/FAIB/cna_faib_1926-2019SY.csv")
faib_cna_intervals$id2<-faib_cna_intervals$id1
faib_cna_intervals<-vpd_cna(faib_cna_intervals)

faib_clim<- intervals_cna(faib_red,faib_cna_intervals,"plot_id")
faib<-inner_join(faib_red,faib_clim,join_by(year.x,year.y,plot_id==site))
faib_norms<-read.csv(file="./Data/FAIB/cna_faib_Normal_1961_1990SY.csv")
faib_norms<-faib_norms[,-2]
faib_norms<-vpd_cna(faib_norms)
colnames(faib_norms)[1]<-"plot_id"
colnames(faib_norms)[5:91]<-paste("P",colnames(faib_norms)[5:91],sep="_")
faib_full<-inner_join(faib,faib_norms,join_by(plot_id))

#### save and read processed file ####

#saveRDS(faib_full,file="./Data/FAIB/faib_climatev2.RDS")
faib_clim<-readRDS(file="./Data/FAIB/faib_climatev2.RDS")

#### AB Data ####

AB_data<-readRDS("./Data/AB_data/AB_datav2.RDS")

G103_clim<-read.csv(file="./Data/AB_data/cna_alberta_input_1980-2015SY.csv")
G103_clim$id2<-G103_clim$id1
G103_clim<-vpd_cna(G103_clim)

AB_clim_ints <- intervals_cna(AB_data,G103_clim,"TRIAL")

AB_clim <- inner_join(AB_data,AB_clim_ints,join_by(year.x,year.y,TRIAL==site))
# add provenances:

AB_prov <- read.csv(file="./Data/AB_data/ab_provs_cna_Normal_1961_1990SY.csv")
AB_prov<-AB_prov[,-2]
AB_prov<-vpd_cna(AB_prov)
colnames(AB_prov)[1]<-"ACC"
colnames(AB_prov)[5:91]<-paste("P",colnames(AB_prov)[5:91],sep="_")

AB_clim <- inner_join(AB_clim,AB_prov,join_by(ACC))
# After adding missing provs, only lose ACC=0 trees (not planted?) from 277A- acceptable

# save data:

#saveRDS(AB_clim,file="./Data/AB_data/ab_climv2.RDS")
AB_clim <- readRDS(file="./Data/AB_data/ab_climv2.RDS")
# checks/working code:
# we lose about 4000 obs with this inner join- must be obs with no provenance recorded?
# temp<-AB_clim %>% filter(ACC==38)
# sum(is.na(AB_data$ACC))
# data_list<-unique(AB_data$ACC)
# prov_list<-unique(AB_prov$ACC)
# sum(AB_data$ACC==0)
# which(!(data_list %in% prov_list))
# data_list[109]
# temp<-AB_data %>% filter(ACC==0)
# 3083 3081 3080 3079 0(277A)
# all these are real seedlots, except "0" from 277A, which is all NA, maybe not planted?
# added to cna file, with coords drawn from description pdfs

#### 410 Data ####

Data_410<-readRDS(file="./Data/410_csv/Data_410v2.csv")

clim_na_410<-read.csv(file="./Data/410_csv/Site_ClimNA_Input_1978-2016SY.csv")
clim_na_410 <- vpd_cna(clim_na_410)
# site in Data_410, and id2 for merge:
clim_410_ints<-intervals_cna(Data_410,clim_na_410,"site")
clim_410 <- inner_join(Data_410,clim_410_ints,join_by(year.x,year.y,site))

provs_410_clim <- read.csv(file="./Data/410_csv/prov_410_cna_Normal_1961_1990SY.csv")
provs_410_clim <- vpd_cna(provs_410_clim)
provs_410_clim<-provs_410_clim[,-2]
colnames(provs_410_clim)[1]<-"Prov"
colnames(provs_410_clim)[5:91]<-paste("P",colnames(provs_410_clim)[5:91],sep="_")


clim_410<-inner_join(clim_410,provs_410_clim,join_by(Prov))

#saveRDS(clim_410,file="./Data/410_csv/clim_410v2.RDS")
clim_410 <- readRDS(file="./Data/410_csv/clim_410v2.RDS")
# lost about ~1000 observations here- investigation follows:

# data_list<-unique(Data_410$Prov)
# prov_list<-unique(provs_410_clim$Prov)
# which(!(data_list %in% prov_list))
# data_list[113]
# temp<-Data_410 %>% filter(Prov==8609)
#### Sx Data ####
sx_trial <- readRDS(file="./Data/Sx_Trial/sx_trial_v2.RDS")

sx_provs<-read.csv("./Data/Sx_Trial/Sx_cna_prov_Normal_1961_1990SY.csv",header = T)
sx_sites<-read.csv("./Data/Sx_Trial/Sx_cna_site_2007-2020SY.csv",header = T)
sx_provs<-vpd_cna(sx_provs)
sx_sites<-vpd_cna(sx_sites)
sx_ints <- intervals_cna(sx_trial,sx_sites,"Site")

sx_clim <- inner_join(sx_trial,sx_ints,join_by(year.x,year.y,Site==site))

colnames(sx_provs)[2]<-"Prov"
sx_provs<-sx_provs[,-1]
colnames(sx_provs)[5:91]<-paste("P",colnames(sx_provs)[5:91],sep="_")
sx_clim <- inner_join(sx_clim,sx_provs,join_by(Prov))
# lose 10000 obs... not a problem, see below
#saveRDS(sx_clim,file="./Data/Sx_Trial/sx_clim_v2.RDS")
sx_clim <- readRDS(file="./Data/Sx_Trial/sx_clim_v2.RDS")
# data_list<-unique(sx_trial$Prov)
# prov_list<-unique(sx_provs$Prov)
# which(!(data_list %in% prov_list))
# data_list[36]
# temp<-sx_trial %>% filter(Prov==9999)
# all are 9999 trees- i.e., not planted. No worries

#### FIA Data ####

# fia_plot <- readRDS("./Data/FIA_DB/fia_plotv2.RDS")
# # convert elev to meters:
# fia_plot$ELEV<-fia_plot$ELEV*0.3048
# min(fia_plot$MEASYEAR,na.rm=T)
# max(fia_plot$MEASYEAR,na.rm=T)
# fia_early <- fia_plot %>% filter(MEASYEAR < 1998)
# fia_main <- fia_plot %>% filter(MEASYEAR >= 1998)
# fia_all_norms<-unique(fia_main[,-5])
# # fia early- 1983 to 1998
# # fia main- 1999 to 2023
# # will actually drop all early- earliest interval is 1998
# # colnames(fia_early)
# # fia_early<-fia_early[,-5]
# fia_main <- fia_main[,-5]
# # fia_early<-unique(fia_early)
# fia_main<-unique(fia_main)
# length(unique(fia_main$ID))
# # length(unique(fia_early$ID))
# length(unique(fia_all_norms$ID))

#write.csv(fia_early,file="./Data/FIA_DB/fia_early.csv")

#write.csv(fia_main,file="./Data/FIA_DB/fia_main_v2.csv")
#write.csv(fia_all_norms,file="./Data/FIA_DB/fia_all_norms_v2.csv")

#fia_clim_early<-read.csv(file="./Data/FIA_DB/fia_early_1983-1998SY.csv")
#fia_clim_main<-read.csv(file="./Data/FIA_DB/fia_main_v2_1998-2023SY.csv")
#fia_clim_annual <- rbind(fia_clim_early,fia_clim_main)
fia_cna_annual <- readRDS("./Data/FIA_DB/CNA_Andy/fia_corrected_annual.RDS")
fia_cna_annual<-vpd_cna(fia_cna_annual)
#rm(fia_clim_early)
#rm(fia_clim_main)

fia_dat <- readRDS(file="./Data/FIA_DB/fia_datv2.RDS")
# the column with the custom "ID" variable needs to be named "id2"
fia_cna_annual$ID.x<-fia_cna_annual$ID.y
colnames(fia_cna_annual)[2]<-"id2"
fia_cna_annual <- fia_cna_annual[,-92]
fia_cna_annual$id2<-as.numeric(fia_cna_annual$id2)

fia_ints<- intervals_cna_fia(fia_dat,fia_cna_annual,"ID")

fia_clim <- inner_join(fia_dat,fia_ints,join_by(year.x,year.y,ID==site))

#fia_norms<- read.csv(file="./Data/FIA_DB/fia_all_norms_v2_Normal_1961_1990SY.csv")
fia_norms <- readRDS("./Data/FIA_DB/CNA_Andy/fia_corrected_per.RDS")
# id1 was changed by csv conversion, could use "as.numeric" but want to preserve same form as ID in fia_dat
fia_norms$ID.x<-fia_norms$ID.y
fia_norms<-fia_norms[,-91]
fia_norms<-fia_norms[,-2]
colnames(fia_norms)[1]<-"ID"
# earlier version from fuzzed coordinates:
# fia_norms$X<-fia_norms$ID
# fia_norms<-fia_norms[,-2]
# colnames(fia_norms)[1]<-"ID"
colnames(fia_norms)
fia_norms<-vpd_cna(fia_norms)
# removing 1994 plot duplicated plot with different location
#fia_norms[c(2028,2292),]
#fia_norms<-fia_norms[-2028,]
colnames(fia_norms)[5:91]<-paste("P",colnames(fia_norms)[5:91],sep="_")


# unique(nchar(fia_norms$ID))
# test<-fia_norms
#fia_norms$ID <- ifelse(nchar(fia_norms$ID)==11,as.character(paste0(0,fia_norms$ID)),as.character(fia_norms$ID))
# unique(nchar(test$ID2))
# rm(test)
# unique(nchar(fia_clim$ID))

#simplest solution is to remove all dups
# dups<-fia_norms[duplicated(fia_norms$ID),]
# dup_list<-dups$ID
# # 
# fia_clim_red <- fia_clim |> filter(!ID %in% dup_list)
# fia_norms_red <- fia_norms |> filter(!ID %in% dup_list)

fia_clim2 <- inner_join(fia_clim,fia_norms,join_by(ID))
#saveRDS(fia_clim2,"./Data/FIA_DB/fia_clim2_v3.RDS")


#fia_clim2 <- inner_join(fia_clim,fia_norms,join_by(ID))

# fia_clim[16781,]$ID
# which(fia_norms$ID=="230102900361")
# fia_norms[c(2028,2292),]
# fia_plot %>% filter(ID == "230102900361")
# temp <- fia_clim_annual %>% filter(id2 == "230102900361")
# 
# fia_norms[12,]
# which(fia_clim$ID == "040200580406")
# temp <- head(fia_clim,n=5)
# 
# which(fia_norms$ID=="230301900717")
# fia_norms[c(1973,2457),]

# dups<-fia_norms[duplicated(fia_norms$ID),]
# dup_list<-dups$ID
# dup_table<- fia_norms %>% filter(ID %in% dup_list)
# okay so the 1994 version of this plot has a different location?
# all subsequent measurements are from later location
# remove fia_norms row 
# all NaN?
# key <- "ID"
# df1<-fia_dat
# df2<-fia_clim_annual
# issue is probable numeric coversion of IDs somewhere in the save csv and run through CNA process

# extract FIA plots for attempt to get real climate vars:


#### Cleanup/compatibility ####

# need to make sure units are all converted- standardizing to meters for height

## standard description ##

# Columns
 
#Site: site or stand (main plot for PSP) ID, needs to be unique for the plot, but the same across remeasurements
#Tree_ID: needs to be unique for the tree, across plots
#source: name of study or inventory, chr, even if numeric
#Prov: within study provenance ID, chr or numeric, "Local" for PSP data
#Type: type of study data, "PSP" or "Prov"
#Latitude, Longitude: Decimal, signed
#Elevation: Meters above SL
#year.x, year.y: first and last year of measurment interval
# height.x, height.y: height in meters at beginning and end of interval
# 

# start with FIA
# as long as we keep IDs/keys for all data, we can always get extra information
fia_merge<-fia_clim2 |> select(ID,TREE,Latitude,Longitude,Elevation,year.x,year.y,height.x,height.y,ht_pct.x,ht_pct.y,dt,dh,starts_with(c("mean_","min_","max_","P_")))
fia_merge$source<-"FIA"
fia_merge$Prov <- "Local"
fia_merge$Type <- "PSP"
# heights were in feet:
# fia_merge$height.x[1:100]
# fia_merge$dh[1:100]
# these appear to still be feet... 
# converted to meters with new code already
# fia_merge <- fia_merge |> mutate(height.x = height.x*0.3048,height.y=height.y*0.3048,dh=dh*0.3048)
fia_merge <- fia_merge |> rename(Site = ID, Tree_ID=TREE)

# Sx
sx_merge <- sx_clim |> select(Site,Tree_ID,Prov,Latitude,Longitude,Elevation,year.x,year.y,height.x,height.y,ht_pct.x,ht_pct.y,dt,dh,starts_with(c("mean_","min_","max_","P_")))
sx_merge$source <- "SX"
sx_merge$Type <- "Prov"
# check- heights may be cm
sx_merge$height.x[1:100]
max(sx_merge$height.x,na.rm=T)
# convert cm to meters:
sx_merge <- sx_merge |> mutate(height.x = height.x/100,height.y=height.y/100,dh=dh/100)

# 410
colnames(clim_410)
merge_410 <- clim_410 |> select(site,ID,Prov,Latitude,Longitude,Elevation,year.x,year.y,height.x,height.y,ht_pct.x,ht_pct.y,dt,dh,starts_with(c("mean_","min_","max_","P_")))
merge_410 <- merge_410 |> rename(Site = site, Tree_ID = ID)
merge_410$source <- "On_410"
merge_410$Type <- "Prov"
# heights could be in cm or dm
# merge_410$height.x[1:100] #meters, no adjustment needed
# max(merge_410$height.x,na.rm=T) # confirmed

# AB
colnames(AB_clim)
head(AB_clim[,1:10])
AB_merge <- AB_clim |> select(TRIAL,ID2,ACC,Latitude,Longitude,Elevation,year.x,year.y,height.x,height.y,ht_pct.x,ht_pct.y,dt,dh,starts_with(c("mean_","min_","max_","P_")))
AB_merge$source <- "AB"
AB_merge$Type <- "Prov"
AB_merge$height.x[1]

# this has been done earlier in processing now: 
# these steps should be moved to the earlier processing block, before RDS write, and applied to AB_clim
# AB_merge$height.x <- unlist(AB_merge$height.x)
# AB_merge$height.y <- unlist(AB_merge$height.y)
# AB Unit Correction

AB_merge <- AB_merge |> mutate(height.x = height.x/100,height.y=height.y/100,dh=dh/100)
#AB_merge$dh[1:100]
AB_merge <- AB_merge |> rename(Site=TRIAL,Tree_ID=ID2,Prov=ACC)
# FAIB
colnames(faib_clim)
# issues with last two VPD columns not having prefix
colnames(faib_clim)[214:215]#<-paste("P",colnames(faib_clim)[214:215],sep="_")
faib_merge <- faib_clim |> select(plot_id,tree_id,Latitude,Longitude,Elevation,year.x,year.y,height.x,height.y,ht_pct.x,ht_pct.y,dt,dh,starts_with(c("mean_","min_","max_","P_")))
faib_merge$source <- "faib"
faib_merge$Prov <- "Local"
faib_merge$Type <- "PSP"
#faib_merge$height.x[1:100] #meters- no adjustment needed
faib_merge <- faib_merge |> rename(Site=plot_id,Tree_ID = tree_id)
#### Merge/rbind ####

ht_dat <- do.call("rbind",list(fia_merge,AB_merge,sx_merge,faib_merge,merge_410))

#### Save Dataset ####


#saveRDS(ht_dat,"./Data/Modelling Datasets/ht_dat_v2.RDS")
ht_dat<-readRDS("./Data/Modelling Datasets/ht_dat_v3.RDS")

#### code to splice out and replace FIA data after new corrected points: ####

# # 492017 obs
# ht_dat <- ht_dat |> filter(source != "FIA")
# #492017-411387
# #80630- exactly the same as fia_clim_2
# ht_dat <- rbind(ht_dat,fia_merge)
# # 492037
# #saveRDS(ht_dat,"./Data/Modelling Datasets/ht_dat_v3.RDS")

#### Misc Code ####


#rm(test)
#test[sample(nrow(test),5),]
# lots of -9999 in some 410 data- look at 8072, 8326 D3, E1


# debugging of FIA intervals code
# df1<-fia_dat
# df2<-fia_clim_annual
# key<-"ID"
# ints<-unique(df1[,c(key,"year.x","year.y")])
# ints[1,]
# ints2<-ints[complete.cases(ints),]
# ints2[1,]
# temp_list <- list()
# is.na(df2) <- df2 == -9999
# #for (i in 1:nrow(ints2)){
# i<-1
#   int_years<-ints2[i,]$year.x:ints2[i,]$year.y
#   int_site <- ints2[i,key]
#   int_site
#   mean_clim<-as.data.frame.list(colMeans(df2 %>% filter(id2==as.numeric(int_site),Year %in% int_years) %>% select(!c(1:6))))
#   # no matches- NaN for everything.
#   df2 %>% filter(id2==as.numeric(int_site))
#   # 0 rows
#   df2$id2[1]
#   # stripped of leading zeros clearly
#   as.numeric(int_site)
#   int_site
#   df2 %>% filter(id2==20112246134)
#   fia_plot %>% filter(ID == "020112246134")
#   # it's in fia_plot
#   which(fia_clim_annual$id2 == 20112246134)
#   # (sprintf("%02d",test1$STATECD),sprintf("%02d",test1$STATECD),sprintf("%03d",test1$COUNTYCD),sprintf("%05d",test1$PLOT))
#   #   new/correct plot1$ID<-paste0(sprintf("%02d",plot1$STATECD),sprintf("%02d",plot1$UNITCD),sprintf("%03d",plot1$COUNTYCD),sprintf("%05d",plot1$PLOT))
#   # what if fia_clim still has erronious IDs?
#   # "020112246134" 2 state, 2 unit, 3 county, 5 plot, vs 2 state, 2 state, 3 county, 5 plot
#   which(df2$id2 == 20212246134)
#   # also doesn't work...
#   fia_dat %>% filter(ID == "020112246134")
#   sum(as.numeric(unique(fia_dat$ID)) %in%  fia_clim_annual$id2)
#   length(as.numeric(unique(fia_dat$ID)))
#   mean_clim <- mean_clim %>% rename_with( ~ paste("mean", .x, sep = "_"))
#   min_clim <- as.data.frame.list(sapply(df2 %>% filter(id2==as.numeric(int_site),Year %in% int_years) %>% select(ends_with("_wt")&!starts_with(c("DD","Eref","CMD","PAS"))),min)) %>%
#     rename_with( ~ paste("min", .x, sep = "_"))
#   list_max <- c("MAT","CMI","DD1040")
#   max_clim <- as.data.frame.list(sapply(df2 %>% filter(id2==as.numeric(int_site),Year %in% int_years) %>% select(ends_with("_sm")|all_of(list_max)),max)) %>%
#     rename_with( ~ paste("max", .x, sep = "_"))
#   clim<-cbind(mean_clim,min_clim,max_clim)
#   clim$site<-paste(int_site)
#   clim$year.x<-min(int_years)
#   clim$year.y<-max(int_years)
#   temp_list[[i]]<-clim
# }
# temp_df<-do.call(rbind,temp_list)
# return(temp_df)
