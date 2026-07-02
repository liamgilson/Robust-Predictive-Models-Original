#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 32- Predict with FAIB data #

################################################################################
## Preds- PSP ##################################################################
################################################################################

faib_psp_tree <- read.csv("./Data/FAIB/faib_psp_tree1.csv", header=TRUE)

counts <- faib_psp_tree %>% count(SAMP_ID, species)
counts <- counts %>% tidyr::pivot_wider(id_cols=SAMP_ID,names_from=species,values_from=n)


test <- faib_psp_tree %>% filter(SAMP_ID=="59144 G000035", meas_yr==2001)
age_trees<-which(!is.na(test$age_tot))
mean(test$height[age_trees])
mean(test$age_tot,na.rm=T)
rm(faib_psp_tree)
rm(test)


head(ht_dat %>% filter(source == "faib"))
faib_plot <- ht_dat %>% filter(Site == "59144 G000035",year.y==2001)

faib585 <- read.csv("./Data/FAIB/pred_faib_13GCMs_ensemble_ssp585_2015-2100SY.csv",header = T)
faib370 <- read.csv("./Data/FAIB/pred_faib_13GCMs_ensemble_ssp370SY.csv",header = T)

faib585<-vpd_cna(faib585)
faib370<-vpd_cna(faib370)
# need to create df with year.x/year.y for 2020-2100, for each tree
# then run int function on this and projected climate data
# then merge with prov data
years.seq<-data.frame(year.x=seq(2015,2100,by=5)[1:17],year.y=seq(2015,2100,by=5)[2:18])
faib_test <-  ht_dat %>% filter(Site == "59144 G000035",year.y==2001,!is.na(dh))
faib_test<- faib_test[,c(1:13,129:218)]
n<-nrow(faib_test)
faib_test <- faib_test[rep(seq_len(n), each = 17), ]
faib_test$year.x <- rep(years.seq$year.x,n)
faib_test$year.y <- rep(years.seq$year.y,n)
faib_test$ht_pct.x<-faib_test$ht_pct.y
faib_test$ht.pct.y<-NA
faib_test$height.x<-faib_test$height.y
faib_test$height.y<-NA
faib_test$ht_pct.x<-ifelse(faib_test$year.x==2015,faib_test$ht_pct.x,NA)
faib_test$height.x<-ifelse(faib_test$year.x==2015,faib_test$height.x,NA)

# merge climate data
faib585$id2<-faib585$id1
faib_585_dat <- intervals_cna(faib_test,faib585,"Site")
faib370$id2<-faib370$id1
faib_370_dat <- intervals_cna(faib_test,faib370,"Site")

faib_pred<-inner_join(faib_test,faib_585_dat,join_by(year.x,year.y,Site==site))
#faib_pred<-inner_join(faib_test,faib_370_dat,join_by(year.x,year.y,Site==site))

# transfer variables
mean_vars <- as.matrix(faib_pred |> ungroup() |> select(starts_with("mean_")))
p_vars <- as.matrix(faib_pred |> ungroup() |> select(starts_with("P_")))
difs <- mean_vars-p_vars
colnames(difs)<-gsub("mean_","T_",colnames(difs))
colnames(difs)

faib_pred <- cbind(faib_pred,difs)
rm(difs)

# predictions

vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "T_DD_0_sp","T_RH","T_Tmax_wt","T_NFFD","T_Tave_sm","T_NFFD_wt")
mod4_unreg<-lm(ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type +  T_DD_0_sp+T_RH+T_Tmax_wt+T_NFFD+T_Tave_sm+T_NFFD_wt,data=ht_tr)
summary(mod4_unreg)

mod4_unreg_pred<-predict_pg(faib_pred,mod4_unreg,"test",n=17)
mod4_unreg_means <- mod4_unreg_pred %>% group_by(year.y) %>% summarize(mean_height=mean(height.y))
mod4_unreg_means  %>% ggplot(aes(x=year.y,y=mean_height))+geom_line()
mod4_unreg_means$model<-"Unregularized"

vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type","mean_PAS_sp","mean_EXT","mean_CMD","max_DD_0_sm",
         "max_CMI_sm" , "max_VPD_max_sm" , "P_Tave_sp" , "P_NFFD_sm" , "P_PAS_sp" , "P_CMD_at" , "P_RH_at", "P_EXT")
mod2 <-predict_pg_reg(faib_pred,"model2",vars2,lambda1=0.27,n=17)
mod2_means<- mod2 %>% group_by(year.y) %>% summarize(mean_height=mean(height.y))
mod2_means %>% ggplot(aes(x=year.y,y=mean_height))+geom_line()
mod2_means$model <- "Model 2"

vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "T_DD_0_sp","T_RH","T_Tmax_wt","T_NFFD","T_Tave_sm","T_NFFD_wt")
# test3 <- ridge_params(ht_tr,vars4,0.01)
# sigma3<-deviance(test3)/(nrow(ht_tr)-11)
mod4 <-predict_pg_reg(faib_pred,"test 2",vars4,lambda1=0.27,n=17)
mod4_means<- mod4 %>% group_by(year.y) %>% summarize(mean_height=mean(height.y))
mod4_means %>% ggplot(aes(x=year.y,y=mean_height))+geom_line()

mod4 %>% filter(Tree_ID==2485) %>% ggplot(aes(x=year.y,y=height.y))+geom_line()


mod4_means$model <- "Model 4"


vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","max_VPD_max_sm",
         "mean_CMD","mean_Eref_sm","T_RH_wt","mean_Tmax_at","mean_DD_0_at")
mod5 <-predict_pg_reg(faib_pred,"model 5",vars5,lambda1=0.2,n=17)
mod5_means <- mod5 %>% group_by(year.y) %>% summarize(mean_height=mean(height.y))
mod5_means$model<- "Model 5"

plottest<-rbind(mod4_unreg_means,mod2_means,mod4_means,mod5_means)

plottest %>% ggplot(aes(x=year.y,y=mean_height,color=model))+geom_line()


mod4_params<-ridge_params(ht_tr,vars4,0.27)
coef(mod4_params)

# New models

vars1<-c("height.x" , "ht_pct.x" , "dt" , "Type" ,  "T_VPD_max_sm" , "T_MAT" , "max_VPD_max_sm" , "mean_CMD" ,  "mean_MAT")
#formula5<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + Elevation + mean_CMD + max_VPD_max_sm"
vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_CMD","max_VPD_max_sm")
#formula7<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + T_DD_0_sp + T_Tmax_wt + T_TD + T_MAT"
vars3<-c("height.x" , "ht_pct.x" , "dt" , "Type", "T_DD_0_sp" , "T_Tmax_wt" , "T_TD" , "T_MAT")
#formula8<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_NFFD_sm + max_CMI_sm + P_NFFD_sm + P_TD + P_bFFP"
vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_NFFD_sm" , "max_CMI_sm" , "P_NFFD_sm" , "P_TD" , "P_bFFP")
#formula9<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + Elevation + mean_Eref_sm + max_VPD_max_sm + mean_CMD+P_Tave_sp+mean_Tmax_at +T_RH_wt"
vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_Eref_sm","max_VPD_max_sm","mean_CMD","P_Tave_sp","T_RH_wt","mean_Tmax_at")
# "optimal" 1SE model from elastic net with blocked folds:
# formula8a<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_NFFD_sm +  mean_PAS_sp + mean_EXT +
#max_DD_0_sm + max_CMI_sm + P_NFFD_sm + P_PAS_sp + P_CMD_at + P_RH_at + P_bFFP + P_EXT + T_RH_at"
vars6<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_NFFD_sm", "mean_PAS_sp", "mean_EXT",
         "max_DD_0_sm", "max_CMI_sm", "P_NFFD_sm", "P_PAS_sp", "P_CMD_at", "P_RH_at", "P_bFFP", "P_EXT", "T_RH_at")

mod1 <-predict_pg_reg(faib_pred,"model1",vars1,lambda1=0.001,n=17)
mod2 <-predict_pg_reg(faib_pred,"model2",vars2,lambda1=0.001,n=17)
mod3 <-predict_pg_reg(faib_pred,"model3",vars3,lambda1=0.001,n=17)
mod4 <-predict_pg_reg(faib_pred,"model4",vars4,lambda1=0.001,n=17)
mod5 <-predict_pg_reg(faib_pred,"model5",vars5,lambda1=0.001,n=17)
mod6 <-predict_pg_reg(faib_pred,"model6",vars6,lambda1=0.001,n=17)
mod7 <-predict_pg_reg(faib_pred,"model7",vars7,lambda1=0.001,n=17)
mods<-rbind(mod1,mod2,mod3,mod4,mod5,mod6,mod7)
mod_means<- mods%>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))


mod_means %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()

mods %>% filter(Tree_ID %in% head(unique(mods$Tree_ID),5)) %>% ggplot(aes(x=year.y,y=height.y,color=Model,group=interaction(Model,Tree_ID)))+geom_line()

############### What about historical time series? #############################

faibhist <- read.csv("./Data/FAIB/pred_faib_1901-2000SY.csv",header = T)
faibhist<-vpd_cna(faibhist)


# new dataframe for historic years

years.seq<-data.frame(year.x=seq(1905,1995,by=5)[1:19],year.y=seq(1910,2000,by=5)[1:19])
faib_test <-  ht_dat %>% filter(Site == "59144 G000035",year.y==2001,!is.na(dh))
faib_test<- faib_test[,c(1:13,129:218)]
n<-nrow(faib_test)
faib_test <- faib_test[rep(seq_len(n), each = 19), ]
faib_test$year.x <- rep(years.seq$year.x,n)
faib_test$year.y <- rep(years.seq$year.y,n)
faib_test$ht_pct.x<-faib_test$ht_pct.y
faib_test$ht.pct.y<-NA
faib_test$height.x<-faib_test$height.y
faib_test$height.y<-NA
faib_test$ht_pct.x<-ifelse(faib_test$year.x==1905,faib_test$ht_pct.x,NA)
faib_test$height.x<-ifelse(faib_test$year.x==1905,faib_test$height.x,NA)


# merge lines

faibhist$id2<-faibhist$id1
faib_hist_dat <- intervals_cna(faib_test,faibhist,"Site")

faib_pred<-inner_join(faib_test,faib_hist_dat,join_by(year.x,year.y,Site==site))

# transfer variables
mean_vars <- as.matrix(faib_pred |> ungroup() |> select(starts_with("mean_")))
p_vars <- as.matrix(faib_pred |> ungroup() |> select(starts_with("P_")))
difs <- mean_vars-p_vars
colnames(difs)<-gsub("mean_","T_",colnames(difs))
colnames(difs)

faib_pred <- cbind(faib_pred,difs)
rm(difs)


mod1 <-predict_pg_reg(faib_pred,"model1",vars1,lambda1=0.001,n=19)
mod2 <-predict_pg_reg(faib_pred,"model2",vars2,lambda1=0.001,n=19)
mod3 <-predict_pg_reg(faib_pred,"model3",vars3,lambda1=0.001,n=19)
mod4 <-predict_pg_reg(faib_pred,"model4",vars4,lambda1=0.001,n=19)
mod5 <-predict_pg_reg(faib_pred,"model5",vars5,lambda1=0.001,n=19)
mods<-rbind(mod1,mod2,mod3,mod4,mod5)
mod_means<- mods%>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))


mod_means %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()
