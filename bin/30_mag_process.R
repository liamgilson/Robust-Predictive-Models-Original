#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 40- Processing of MAGPlot data for external validation #


############## MAGPlot ########################################################
# for MAGPlot V. 1.2

# Dr. Liam Gilson
# May 2026

#### dependencies

library(dplyr)

#### Read-in data ####

magp_sites <- read.csv("./Data/MAGPlot/magp_sites.csv")
magp_trees <- read.csv("./Data/MAGPlot/magp_trees.csv")
magp_header <- read.csv("./Data/MAGPlot/magp_tree_header.csv")

colnames(magp_sites)
colnames(magp_trees)

unique(magp_sites$province)

qc_sites <- magp_sites %>% filter(province == "QC")

qc_ids <- qc_sites %>% pull(magp_site_id)

qc_trees <- magp_trees %>% filter(magp_site_id %in% qc_ids)
# temporary step to avoid loading magp_trees back into memory if something goes wrong
qc_trees_full <- qc_trees

rm(magp_trees)

qc_years <- magp_header %>% filter(magp_site_id %in% qc_ids) %>% select(magp_site_id,meas_num,meas_year,meas_month)
qc_years<- unique(qc_years)
qc_trees <- left_join(qc_trees,qc_years,by=join_by(magp_site_id,meas_num))
colnames(qc_trees)
# ecdf

# this is some kind of magplot error code, we remove all these trees
# -1 is not a real measurement, this would bias results, percentiles, etc.
# unclear why they don't use NULL or NA
qc_trees<-qc_trees[-which(qc_trees$height == -1),]

qc_trees <- qc_trees |> group_by(magp_site_id,meas_year)%>% mutate(ht_pct = ecdf(height)(height))
# dbh in cm
# BA not needed here
# qc_trees$BA <- 

qc_spruce <- qc_trees %>% filter(species_gs %in% c("PICE.GLA","PICE.ENG"), meas_est_height == "M")

#### Code for intervals: ####
# based on FAIB code:

tree_list<- unique(qc_spruce$tree_id)

mag_res<-list()
for (i in 1:length(tree_list)){
  tree<-tree_list[i]
  test <- qc_spruce %>% filter(tree_id == tree) %>% arrange(treemeas_num)
  #test <- faib_spruce %>% filter(SAMP_ID == test_tree)
  #len <- ifelse((nrow(test) %% 2)==0,nrow(test),nrow(test)-1)
  len<-nrow(test)
  full_len<-nrow(test)
  if(len >= 2){
    year.x<-test$meas_year[1:full_len-1]
    year.y<-test$meas_year[2:full_len]
    height.x<-test$height[1:full_len-1]
    height.y<-test$height[2:full_len]
    ht_pct.x<-test$ht_pct[1:full_len-1]
    ht_pct.y<-test$ht_pct[2:full_len]
    tree_id<-rep(tree,len-1)
    plot_id<-rep(test$magp_site_id[1],len-1)
    tree_num<-rep(test$tree_num[1],len-1)
    mag_res[[i]]<-cbind(plot_id,tree_id,tree_num,year.x,year.y,height.x,height.y,ht_pct.x,ht_pct.y)
  } else {next}
}

mag_proc<-do.call("rbind",mag_res)
mag_proc<-as.data.frame(mag_proc)
mag_proc <- mag_proc %>% mutate_at(c(1,3:9),as.numeric)

mag_proc$dh<-mag_proc$height.y-mag_proc$height.x
mag_proc$dt<-mag_proc$year.y-mag_proc$year.x

mag_proc %>% ggplot(aes(x=dh))+geom_histogram()

saveRDS(mag_proc,file="./Data/MAGplot/mag_red.RDS")


#### Climate NA section: ####

### Generate plot list for use in ClimateNA ####

spruce_sites_mag<-unique(mag_proc$plot_id)
mag_cna<-qc_sites %>% filter(magp_site_id %in% spruce_sites_mag)
cna_input_mag<-mag_cna %>% select(magp_site_id,latitude,longitude,elevation) %>% mutate(id1 = magp_site_id, id2=NA, lat=latitude,long=longitude, elev=elevation) %>%
  select(id1,id2,lat,long,elev)
# this still needs some edits before input to CNA, but faster to do them in excel
write.csv(cna_input_mag,file="./Data/MAGplot/cna_mag.csv")

# range
min(mag_proc$year.x)
max(mag_proc$year.y)
# 1970-2023
#### Read processed ClimateNA data, compute intervals, merge ####

mag_cna_intervals<-read.csv(file="./Data/MAGplot/cna_mag_1970-2023SY.csv")
mag_norms<-read.csv(file="./Data/MAGplot/cna_mag_Normal_1961_1990SY.csv")
mag_cna_intervals$id2<-mag_cna_intervals$id1
mag_cna_intervals<-vpd_cna(mag_cna_intervals)

# intervals
mag_clim<- intervals_cna(mag_proc,mag_cna_intervals,"plot_id")
mag_clim$site <- as.numeric(mag_clim$site)
mag<-inner_join(mag_proc,mag_clim,join_by(year.x,year.y,plot_id==site))

# normals

mag_norms<-mag_norms[,-2]
mag_norms<-vpd_cna(mag_norms)
colnames(mag_norms)[1]<-"plot_id"
colnames(mag_norms)[5:91]<-paste("P",colnames(mag_norms)[5:91],sep="_")

mag_full<-inner_join(mag,mag_norms,join_by(plot_id))

# transfers

mag_full$source <- "qc"
mag_full$Type <- "PSP"

mean_vars <- as.matrix(mag_full |> ungroup() |> select(starts_with("mean_")))
p_vars <- as.matrix(mag_full |> ungroup() |> select(starts_with("P_")))
difs <- mean_vars-p_vars
colnames(difs)<-gsub("mean_","T_",colnames(difs))
colnames(difs)

mag_full <- cbind(mag_full,difs)
rm(difs)


saveRDS(mag_full,file="./Data/MAGplot/mag_full.RDS")

#### Predictions ####

# models
# Best subsets, hypothetical pool
vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_CMD","max_VPD_max_sm","mean_MAT","mean_Eref")
# 50 variable blocks through all transfer variables
vars3<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_CMD" ,"max_VPD_max_sm", "T_DD_0_sp" , "Elevation" , "mean_Tave_sp")
# elastic net model
vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "Elevation", "max_DD_0_sm","max_CMI_sm","P_Tmin_sp","P_NFFD_sm")
# elastic net fed back into best subsets
vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_CMD","max_VPD_max_sm","mean_DD_0_at","P_DD_0_sp")
# "optimal" 1SE model from elastic net with blocked folds:
vars6<-c("height.x" , "ht_pct.x" , "dt" , "Type", "Elevation", "max_CMI_sm","P_NFFD_sm")
# stepwise variables
vars7<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","P_Tmin_sp","max_VPD_max_sm","mean_CMD","T_Tmax_wt")



#ridge_error(ht_tr,mag_full,vars1,0.01)
ridge_error(ht_tr,mag_full,vars2,0.01)
ridge_error(ht_tr,mag_full,vars3,0.01)
ridge_error(ht_tr,mag_full,vars4,0.01)
ridge_error(ht_tr,mag_full,vars5,0.01)
ridge_error(ht_tr,mag_full,vars6,0.01)
ridge_error(ht_tr,mag_full,vars7,0.01)

ridge_error(ht_tr,mag_full,vars2,0.01,type="bias")
ridge_error(ht_tr,mag_full,vars3,0.01,type="bias")
ridge_error(ht_tr,mag_full,vars4,0.01,type="bias")
ridge_error(ht_tr,mag_full,vars5,0.01,type="bias")
ridge_error(ht_tr,mag_full,vars6,0.01,type="bias")
ridge_error(ht_tr,mag_full,vars7,0.01,type="bias")

# variance estimate

var(ihs(mag_full$dh))


#### Without Elevation ####

#Best subsets, hypothetical pool
#formula5<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + Elevation + mean_CMD + max_VPD_max_sm"
# Elevation + max_VPD_max_sm + + mean_CMD + mean_MAT + mean_Eref
#vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_CMD","max_VPD_max_sm","mean_MAT","mean_Eref")
vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type" ,  "T_VPD_max_sm" , "T_MAT" , "max_VPD_max_sm" , "mean_CMD" ,  "mean_MAT")

# 50 variable blocks through all transfer variables
vars3<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_CMD" ,"max_VPD_max_sm", "T_Tmax_wt" , "P_NFFD_sm" , "P_Tmin_sp")
# elastic net model
#vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_NFFD_sm" , "max_CMI_sm" , "P_NFFD_sm" , "P_TD" , "P_bFFP")
vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "max_CMI_sm" , "P_Tmin_sm" , "P_NFFD_sm" , "P_TD" , "P_bFFP")
# elastic net fed back into best subsets
#vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_CMD","max_VPD_max_sm","mean_DD_0_at","P_DD_0_sp")
vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","mean_Tmax_at","mean_Eref_at","max_VPD_max_sm","mean_CMI_sm","T_Eref_sp")
# "optimal" 1SE model from elastic net with blocked folds:
# vars6<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_NFFD_sm", "mean_PAS_sp", "mean_EXT",
#            "max_DD_0_sm", "max_CMI_sm", "P_NFFD_sm", "P_PAS_sp", "P_CMD_at", "P_RH_at", "P_bFFP", "P_EXT", "T_RH_at")
vars6<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_NFFD_sm", "mean_PAS_sp","max_CMI_sm", 
         "P_NFFD_sm", "P_PAS_sp", "P_CMD_at", "P_RH_at", "P_bFFP", "P_EXT")
# stepwise variables
vars7<-c("height.x" , "ht_pct.x" , "dt" , "Type","mean_Tmax_at","mean_NFFD","mean_MAT","T_DD18_sm","max_CMI")


#### lambda table:

ridge_error(ht_tr,mag_full,vars3,0.0)
ridge_error(ht_tr,mag_full,vars3,0.0,type="bias")
ridge_error(ht_tr,mag_full,vars3,0.01)
ridge_error(ht_tr,mag_full,vars3,0.01,type="bias")
ridge_error(ht_tr,mag_full,vars3,0.22)
ridge_error(ht_tr,mag_full,vars3,0.22,type="bias")
ridge_error(ht_tr,mag_full,vars4,0.0)
ridge_error(ht_tr,mag_full,vars4,0.0,type="bias")
ridge_error(ht_tr,mag_full,vars4,0.01)
ridge_error(ht_tr,mag_full,vars4,0.01,type="bias")
ridge_error(ht_tr,mag_full,vars4,0.22)
ridge_error(ht_tr,mag_full,vars4,0.22,type="bias")
ridge_error(ht_tr,mag_full,vars7,0.0)
ridge_error(ht_tr,mag_full,vars7,0.0,type="bias")
ridge_error(ht_tr,mag_full,vars7,0.01)
ridge_error(ht_tr,mag_full,vars7,0.01,type="bias")
ridge_error(ht_tr,mag_full,vars7,0.22)
ridge_error(ht_tr,mag_full,vars7,0.22,type="bias")
