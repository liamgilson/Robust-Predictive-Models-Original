#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 03 Processing of BC Faib Data #


########################## FAIB Processing Code- V2 ############################

################################################################################

######### Dependencies #########

library(dplyr)
library(ggplot2)

########## Read-in ############

# read in data, size is ~800 mb, so this takes time

# Note- R reads this into RAM, so be careful!

#### Tree level information ####
faib_psp_tree <- read.csv("./Data/faib/faib_psp_tree1.csv", header=TRUE)
str(faib_psp_tree)
# all variables appear to be recognized correctly
na_ecdf <- function(x){if(all(is.na(x))){return(NA)}else{return(ecdf(x)(x))}}
#test<-c(1,2,3,NA)
faib_psp_tree <- faib_psp_tree |> filter (ld %in% c("L","I","V"))|> group_by(SAMP_ID,meas_yr) |> mutate(ht_pct = na_ecdf(height)) |> arrange(desc(baha)) |> mutate(BAL = cumsum(baha))
# need to check if this worked
#faib_psp_tree$SAMP_ID[1]
#test<-faib_psp_tree %>% filter(SAMP_ID == "4005149-PSP")

#### Plot level information ####
faib_plots<- read.csv("./Data/faib/sample.csv")
str(faib_plots)
# all variables appear to be recognized correctly

# we can filter to three tree species: SW, SX, SE

faib_spruce <- faib_psp_tree %>% filter(species %in% c("SW","SE","SX"))
# get rid of huge original data file
rm(faib_psp_tree)

# SAMP_ID is site, but there's also tree: tree_no
length(unique(faib_spruce$tree_no))
# tree_no isn't unique, need to make an actual, tree-level ID
faib_spruce$tree_id<-paste0(faib_spruce$SAMP_ID,faib_spruce$tree_no)

#### testing code
# faib_res<-list()
# test_tree<-faib_spruce$SAMP_ID[1]
# test <- faib_spruce %>% filter(SAMP_ID == test_tree)
# len <- ifelse((nrow(test) %% 2)==0,nrow(test),nrow(test)-1)
# full_len<-nrow(test)
# if(len >= 2){
# year.x<-test$meas_yr[(full_len-len+1):full_len-1]
# year.y<-test$meas_yr[(full_len-len+2):full_len]
# height.x<-test$ht_meas[(full_len-len+1):full_len-1]
# height.y<-test$ht_meas[(full_len-len+2):full_len]
# samp_id<-rep(test_tree,len-1)
# plot_num<-rep(test$plot_no[1],len-1)
# faib_res[[i]]<-cbind(samp_id,plot_num,year.x,year.y,height.x,height.y)
# } else {next}

#### Interval code ####

tree_list<- unique(faib_spruce$tree_id)

faib_res<-list()
for (i in 1:length(tree_list)){
  tree<-tree_list[i]
  test <- faib_spruce %>% filter(tree_id == tree) %>% arrange(meas_no)
  #test <- faib_spruce %>% filter(SAMP_ID == test_tree)
  #len <- ifelse((nrow(test) %% 2)==0,nrow(test),nrow(test)-1)
  len<-nrow(test)
  full_len<-nrow(test)
  if(len >= 2){
    year.x<-test$meas_yr[1:full_len-1]
    year.y<-test$meas_yr[2:full_len]
    height.x<-test$ht_meas[1:full_len-1]
    height.y<-test$ht_meas[2:full_len]
    ht_pct.x<-test$ht_pct[1:full_len-1]
    ht_pct.y<-test$ht_pct[2:full_len]
    BAL.x<-test$BAL[1:full_len-1]
    BAL.y<-test$BAL[2:full_len]
    tree_id<-rep(tree,len-1)
    plot_id<-rep(test$SAMP_ID[1],len-1)
    tree_num<-rep(test$tree_no[1],len-1)
    plot_num<-rep(test$plot_no[1],len-1)
    faib_res[[i]]<-cbind(plot_id,plot_num,tree_id,tree_num,year.x,year.y,height.x,height.y,ht_pct.x,ht_pct.y,BAL.x,BAL.y)
  } else {next}
}
faib_proc<-do.call("rbind",faib_res)

#### Clean-up ####

faib_df<-as.data.frame(faib_proc)
faib_red<-faib_df[complete.cases(faib_df[,c("height.x","height.y")]),]
faib_red$height.x<-as.numeric(faib_red$height.x)
faib_red$height.y<-as.numeric(faib_red$height.y)
faib_red$year.x<-as.numeric(faib_red$year.x)
faib_red$year.y<-as.numeric(faib_red$year.y)
faib_red$dh<-faib_red$height.y-faib_red$height.x
faib_red$dt<-faib_red$year.y-faib_red$year.x
  
faib_red %>% ggplot(aes(x=dh))+geom_histogram()

faib_plots[faib_plots$SAMP_ID==test_tree,]

saveRDS(faib_red,file="./Data/FAIB/faib_dat.RDS")
#### Generate plot list for use in ClimateNA ####

faib_cna_plots<-unique(faib_red$plot_id)
faib_cna<-faib_plots %>% filter(SAMP_ID %in% faib_cna_plots)
faib_cna<-faib_cna %>% select(SAMP_ID,latitude,longitude,elev,slope,aspect)
faib_cna<-unique(faib_cna)
cna_input_faib<-faib_cna %>% select(SAMP_ID,latitude,longitude,elev) %>% mutate(id1 = SAMP_ID, lat=latitude,long=longitude) %>%
  select(id1,lat,long,elev)
# this still needs some edits before input to CNA, but faster to do them in excel
write.csv(cna_input_faib,file="./Data/faib/cna_faib.csv")

#### Read processed ClimateNA data, compute intervals, merge ####

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

#saveRDS(faib_full,file="./Data/FAIB/faib_climateV2.RDS")
#faib_clim<-readRDS(file="./Data/FAIB/faib_climateV2.RDS")

# clean up

# rm(faib_cna_intervals)
# rm(faib)
# rm(faib_proc)
# rm(faib_res)
# rm(faib_plots)

# tests for open code (March 2026)
test<-faib_plots %>% filter(clstr_id == "4005149-PSP")
