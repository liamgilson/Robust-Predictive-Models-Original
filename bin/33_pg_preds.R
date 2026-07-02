#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 33- Predictions with Prince George Trial #


################################################################################
# Predictions- PG test
################################################################################
# Liam Gilson

# packages
library(dplyr)
library(ggplot2)
library(glmnet)

# core functions

ihs <- function(x){log(x + sqrt(x^2 + 1))}
ihs_back <- function(x){exp((summary(x)$sigma^2)/2)}

intervals_cna<-function(df1,df2,key){
  ints<-unique(df1[,c(key,"year.x","year.y")])
  ints2<-ints[complete.cases(ints),]
  temp_list <- list()
  is.na(df2) <- df2 == -9999
  for (i in 1:nrow(ints2)){
    int_years<-ints2[i,]$year.x:ints2[i,]$year.y
    int_site <- ints2[i,key]
    mean_clim<-as.data.frame.list(colMeans(df2 %>% filter(id2==paste(int_site),Year %in% int_years) %>% select(!c(1:6))))
    mean_clim <- mean_clim %>% rename_with( ~ paste("mean", .x, sep = "_"))
    min_clim <- as.data.frame.list(sapply(df2 %>% filter(id2==paste(int_site),Year %in% int_years) %>% select(ends_with("_wt")&!starts_with(c("DD","Eref","CMD","PAS"))),min)) %>%
      rename_with( ~ paste("min", .x, sep = "_"))
    list_max <- c("MAT","CMI","DD1040")
    max_clim <- as.data.frame.list(sapply(df2 %>% filter(id2==paste(int_site),Year %in% int_years) %>% select(ends_with("_sm")|all_of(list_max)),max)) %>%
      rename_with( ~ paste("max", .x, sep = "_"))
    clim<-cbind(mean_clim,min_clim,max_clim)
    clim$site<-paste(int_site)
    clim$year.x<-min(int_years)
    clim$year.y<-max(int_years)
    temp_list[[i]]<-clim
  }
  temp_df<-do.call(rbind,temp_list)
  return(temp_df)
}

# pg data

pg_dat<-readRDS(file="./Data/PG_Trial/pg_dat.RDS")
pg_provs<-read.csv("./Data/PG_Trial/prov_clim_Normal_1961_1990SY.csv",header = T)
pg_585<-read.csv("./Data/PG_Trial/site_clim_13GCMs_ensemble_ssp585_2017-2100SY.csv",header = T)
pg_370<-read.csv("./Data/PG_Trial/site_clim_13GCMs_ensemble_ssp370_2017-2100SY.csv",header = T)
pg_585<-vpd_cna(pg_585)
pg_370<-vpd_cna(pg_370)
pg_provs<-vpd_cna(pg_provs)
# need to create df with year.x/year.y for 2020-2100, for each tree
# then run int function on this and projected climate data
# then merge with prov data
years.seq<-data.frame(year.x=seq(2020,2100,by=5)[1:16],year.y=seq(2020,2100,by=5)[2:17])
pg_test<-pg_dat %>% filter(Site=="Church",year.y==2017,!is.na(dh))
pg_test<-pg_test[,c(1:14,130:218)]
n<-nrow(pg_test)
pg_test <- pg_test[rep(seq_len(n), each = 16), ]
pg_test$year.x <- rep(years.seq$year.x,n)
pg_test$year.y <- rep(years.seq$year.y,n)
pg_test$ht_pct.x<-pg_test$ht_pct.y
pg_test$ht.pct.y<-NA
pg_test$height.x<-pg_test$height.y
pg_test$height.y<-NA
pg_test$ht_pct.x<-ifelse(pg_test$year.x==2020,pg_test$ht_pct.x,NA)
pg_test$height.x<-ifelse(pg_test$year.x==2020,pg_test$height.x,NA)

# merge climate data
pg_585$id2<-pg_585$id1
pg_585_dat <- intervals_cna(pg_test,pg_585,"Site")

pg_pred<-inner_join(pg_test,pg_585_dat,join_by(year.x,year.y,Site==site))

# transfer variables
mean_vars <- as.matrix(pg_pred |> ungroup() |> select(starts_with("mean_")))
p_vars <- as.matrix(pg_pred |> ungroup() |> select(starts_with("P_")))
difs <- mean_vars-p_vars
colnames(difs)<-gsub("mean_","T_",colnames(difs))
colnames(difs)

pg_pred <- cbind(pg_pred,difs)
rm(difs)


# prediction functions, one for lm, one for glmnet

# matrix setup


# results<-data.frame(TreeKey=rep(treelist$TreeKey,each=21),Year=rep(seq(2000,2100,by=5),78),DBH=NA,BAL=NA)
# temp <- treelist %>% mutate(DBH=DBH.x,Year=2000) %>% select(TreeKey,Year,DBH,BAL)
# results<-left_join(results,temp,join_by(TreeKey,Year))
# results<-results %>% rename(DBH=DBH.y,BAL=BAL.y) %>% select(TreeKey,Year,DBH,BAL)


# so this has to predict for each year, ca
predict_pg<-function(resmat,model,name,n=16){
  for (i in 1:n){
    year<-years.seq$year.y[i]
    tempdat<-resmat %>% filter(year.y==year)
    tempdat$dt<-5
    tempdat$dh<-ihs_back(model)*sinh(predict(model,newdata=tempdat))
    tempdat$height.y<-tempdat$height.x+tempdat$dh
    tempdat$ht_pct.y<-ecdf(tempdat$height.y)(tempdat$height.y)
    resmat[which(resmat$year.y==year),]$height.y<-tempdat$height.y
    resmat[which(resmat$year.y==year),]$ht_pct.y<-tempdat$ht_pct.y
    if (i < n){
      resmat[which(resmat$year.x==year),]$height.x<-tempdat$height.y
      resmat[which(resmat$year.x==year),]$ht_pct.x<-tempdat$ht_pct.y
    }
  }
  resmat$Model<-name
  return(resmat)
}

# Okay so a few notes here: functionally we need to fit the glmnet model object, so we need to create a clean model matrix with only the required variables
# lines 1-5 perform this, using hardcoded ht_tr to fit the model
# note especially converting any factor variables to numeric binary, or we get a "coerced NAs" error from glmnet
# we get sigma out of this model, so no need to provide it
# also note that we need to remove any factors from the validation (prediction) matrix or we get an error from predict.glmnet
# it is not known why the output of the ridge_params function, which should also be a glmnet model object, cannot be used in predict.glmnet
# also predict is producing a weird format, only a vector of values but a list of two names? needs to be put through "as.vector()" or it produces errors on assigment
predict_pg_reg<-function(resmat,name,vars,lambda1=0.01,n=16){
  for (i in 1:n){
    x<-ht_tr[,vars]
    x$Type <-ifelse(x$Type=="PSP",1,0)
    x$logheight<-log(x$height.x)
    y<- ihs(ht_tr$dh)
    mod_temp <- glmnet(x,y,alpha=0,lambda=lambda1)
    sigma<-deviance(mod_temp)/(nrow(ht_tr)-(length(vars)+1))
    year<-years.seq$year.y[i]
    tempdat<-resmat %>% filter(year.y==year)
    tempdat$dt<-5
    moddat <- tempdat[,vars]
    moddat$Type <- ifelse(moddat$Type=="PSP",1,0)
    moddat$logheight<-log(moddat$height.x)
    tempdat$dh<-exp(sigma/2)*sinh(predict(mod_temp,newx=as.matrix(moddat)))
    tempdat$height.y<-tempdat$height.x+as.vector(tempdat$dh)
    tempdat$ht_pct.y<-ecdf(tempdat$height.y)(tempdat$height.y)
    resmat[which(resmat$year.y==year),]$height.y<-tempdat$height.y
    resmat[which(resmat$year.y==year),]$ht_pct.y<-tempdat$ht_pct.y
    if (i < n){
      resmat[which(resmat$year.x==year),]$height.x<-tempdat$height.y
      resmat[which(resmat$year.x==year),]$ht_pct.x<-tempdat$ht_pct.y
    }
  }
  resmat$Model<-name
  return(resmat)
}

vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "T_DD_0_sp","T_RH","T_Tmax_wt","T_NFFD","T_Tave_sm","T_NFFD_wt")
mod4_unreg<-lm(ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type +  T_DD_0_sp+T_RH+T_Tmax_wt+T_NFFD+T_Tave_sm+T_NFFD_wt,data=ht_tr)
summary(mod4_unreg)

mod4_unreg_pred<-predict_pg(pg_pred,model,"test")
mod4_unreg_means <- test %>% group_by(year.y) %>% summarize(mean_height=mean(height.y))
mod4_unreg_means  %>% ggplot(aes(x=year.y,y=mean_height))+geom_line()
mod4_unreg_means$model<-"Unregularized"


vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type","mean_PAS_sp","mean_EXT","mean_CMD","max_DD_0_sm",
         "max_CMI_sm" , "max_VPD_max_sm" , "P_Tave_sp" , "P_NFFD_sm" , "P_PAS_sp" , "P_CMD_at" , "P_RH_at", "P_EXT")
mod2 <-predict_pg_reg(pg_pred,"model2",vars2,lambda1=0.27)
mod2_means<- mod2 %>% group_by(year.y) %>% summarize(mean_height=mean(height.y))
mod2_means %>% ggplot(aes(x=year.y,y=mean_height))+geom_line()
mod2_means$model <- "Model 2"

vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "T_DD_0_sp","T_RH","T_Tmax_wt","T_NFFD","T_Tave_sm","T_NFFD_wt")
# test3 <- ridge_params(ht_tr,vars4,0.01)
# sigma3<-deviance(test3)/(nrow(ht_tr)-11)
mod4 <-predict_pg_reg(pg_pred,"test 2",vars4,lambda1=0.27)
mod4_means<- mod4 %>% group_by(year.y) %>% summarize(mean_height=mean(height.y))
mod4_means %>% ggplot(aes(x=year.y,y=mean_height))+geom_line()

mod4 %>% filter(Tree_ID==2485) %>% ggplot(aes(x=year.y,y=height.y))+geom_line()


mod4_means$model <- "Model 4"


vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","max_VPD_max_sm",
         "mean_CMD","mean_Eref_sm","T_RH_wt","mean_Tmax_at","mean_DD_0_at")
mod5 <-predict_pg_reg(pg_pred,"model 5",vars5,lambda1=0.2)
mod5_means <- mod5 %>% group_by(year.y) %>% summarize(mean_height=mean(height.y))
mod5_means$model<- "Model 5"

plottest<-rbind(mod4_unreg_means,mod2_means,mod4_means,mod5_means)

plottest %>% ggplot(aes(x=year.y,y=mean_height,color=model))+geom_line()

# new models, new code

vars1<-c("height.x" , "ht_pct.x" , "dt" , "Type" ,  "T_VPD_max_sm" , "T_MAT" , "max_VPD_max_sm" , "mean_CMD" ,  "mean_MAT")
#formula5<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + Elevation + mean_CMD + max_VPD_max_sm"
vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_CMD","max_VPD_max_sm")
#formula7<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + T_DD_0_sp + T_Tmax_wt + T_TD + T_MAT"
vars3<-c("height.x" , "ht_pct.x" , "dt" , "Type", "T_DD_0_sp" + "T_Tmax_wt" + "T_TD" + "T_MAT")
#formula8<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_NFFD_sm + max_CMI_sm + P_NFFD_sm + P_TD + P_bFFP"
vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_NFFD_sm" , "max_CMI_sm" , "P_NFFD_sm" , "P_TD" , "P_bFFP")
#formula9<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + Elevation + mean_Eref_sm + max_VPD_max_sm + mean_CMD+P_Tave_sp+mean_Tmax_at +T_RH_wt"
vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_Eref_sm","max_VPD_max_sm","mean_CMD","P_Tave_sp","T_RH_wt","mean_Tmax_at")

mod1 <-predict_pg_reg(pg_pred,"model1",vars1,lambda1=0.001,n=16)
mod2 <-predict_pg_reg(pg_pred,"model2",vars2,lambda1=0.001,n=16)
mod3 <-predict_pg_reg(pg_pred,"model3",vars3,lambda1=0.001,n=16)
mod4 <-predict_pg_reg(pg_pred,"model4",vars4,lambda1=0.001,n=16)
mod5 <-predict_pg_reg(pg_pred,"model5",vars5,lambda1=0.001,n=16)
mod6 <-predict_pg_reg(pg_pred,"model6",vars6,lambda1=0.001,n=16)
mod7 <-predict_pg_reg(pg_pred,"model7",vars7,lambda1=0.001,n=16)
mods<-rbind(mod1,mod2,mod3,mod4,mod5,mod6,mod7)
mod_means<- mods %>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))


mod_means %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()

