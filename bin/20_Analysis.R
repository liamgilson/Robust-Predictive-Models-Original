#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 30- Analysis- Training/validation split, variable selection, internal validation #


#### Analysis of height data with validation split ####

# dependencies

library(dplyr)
library(ggplot2)
library(glmnet)

# for stepwise section:
library(stringr)

source("dis_functions.R")

# data read

ht_dat<-readRDS("./Data/Modelling Datasets/ht_dat_v4.RDS")
# new processed version with transfer variables and alberta cm measurements fixed
# should be able to just load and use
ht_dat_na <- ht_dat |> filter(!is.na(dh))

# probably 10% quantile, but could look at 30%, don't want to eliminate a whole data source
## This fails as it's almost entirely Ontario obs, with no rep. of Alberta trials
#quantile(ht_dat_na$mean_Tave_sm, 0.9)
#test <- ht_dat_na %>% filter(mean_Tave_sm > 17.9)
#test %>% count(source)

ht_dat_na %>% count(source)

# could do seperate quantiles for each study- as-is it's 25% of Ontario and no Alberta

ab_break <- ht_dat_na |> filter(source=="AB") |> pull(mean_Tave_sm) |> quantile(0.9)
fia_break <- ht_dat_na |> filter(source=="FIA") |> pull(mean_Tave_sm) |> quantile(0.9)
sx_break <- ht_dat_na |> filter(source=="SX") |> pull(mean_Tave_sm) |> quantile(0.9)
O410_break <- ht_dat_na |> filter(source=="On_410") |> pull(mean_Tave_sm) |> quantile(0.9)
faib_break <- ht_dat_na |> filter(source=="faib") |> pull(mean_Tave_sm) |> quantile(0.9)


temp_ab <- ht_dat_na |> filter(source=="AB") |> filter(mean_Tave_sm > ab_break)
temp_fia <- ht_dat_na |> filter(source=="FIA") |> filter(mean_Tave_sm > fia_break)
temp_sx <- ht_dat_na |> filter(source=="SX") |> filter(mean_Tave_sm > sx_break)
temp_410 <- ht_dat_na |> filter(source=="On_410") |> filter(mean_Tave_sm > O410_break)
temp_faib <- ht_dat_na |> filter(source=="faib") |> filter(mean_Tave_sm > faib_break)

ht_val <- do.call("rbind",list(temp_ab,temp_fia,temp_sx,temp_410,temp_faib))

temp_ab <- ht_dat_na |> filter(source=="AB") |> filter(mean_Tave_sm <= ab_break)
temp_fia <- ht_dat_na |> filter(source=="FIA") |> filter(mean_Tave_sm <= fia_break)
temp_sx <- ht_dat_na |> filter(source=="SX") |> filter(mean_Tave_sm <= sx_break)
temp_410 <- ht_dat_na |> filter(source=="On_410") |> filter(mean_Tave_sm <= O410_break)
temp_faib <- ht_dat_na |> filter(source=="faib") |> filter(mean_Tave_sm <= faib_break)

ht_tr <- do.call("rbind",list(temp_ab,temp_fia,temp_sx,temp_410,temp_faib))

rm(temp_ab)
rm(temp_fia)
rm(temp_sx)
rm(temp_410)
rm(temp_faib)

ht_tr %>% count(source)
ht_val %>% count(source)
#### Actual analysis ####

# folds generated here now, since we're going to base everything on internal (training set) CV
cv_folds_tr <- cv_climate(ht_tr,"mean_Tave_sm")

# exploratory best subsets- will this be the same as full training data?

test_vars<- c("mean_Tmax_sm","mean_Tmin_wt","mean_CMD","mean_DD1040","mean_VPD_max_sm","max_Tmax_sm","max_VPD_max_sm","mean_CMI",
              "T_CMI","T_DD1040","T_VPD_max_sm","T_Tmax_wt","T_Tmax_sm","T_CMD","T_MAT","T_MAP","mean_MAT","mean_MAP")
eq_vars<-paste0(test_vars,collapse="+")

mods_ex <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) + ht_pct.x + dt + Type +",eq_vars)),data=ht_tr,nvmax=20)
plot(mods_ex)

formula1<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type +  max_Tmax_sm"
formula2<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type +  max_Tmax_sm + mean_MAT"
formula3<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type +  T_VPD_max_sm + T_MAT + max_VPD_max_sm + mean_CMD + mean_MAT"
formula4<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type +  T_VPD_max_sm + max_VPD_max_sm + mean_CMD + T_DD1040"
cv_function(ht_tr,cv_folds_tr,formula1)
cv_function(ht_tr,cv_folds_tr,formula2)
cv_function(ht_tr,cv_folds_tr,formula3)
cv_function(ht_tr,cv_folds_tr,formula4)
# prefer formula 3 here

test_vars2<- c("mean_Tmax_sm","mean_Tmin_wt","mean_CMD","mean_DD1040","mean_VPD_max_sm","max_Tmax_sm","max_VPD_max_sm","mean_CMI","mean_MSP",
               "mean_Eref","T_Eref_sm",
               "T_CMI","T_CMD","T_DD1040","T_VPD_max_sm","T_MAT","T_MAP","T_DD_0","T_Tave_wt","T_Tmax_sm","mean_MAT","mean_MAP")
eq_vars2<-paste0(test_vars2,collapse="+")
mods_ex <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) + ht_pct.x + dt + Type +",eq_vars2)),data=ht_tr,nvmax=20)
plot(mods_ex)

formula5<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_CMD + max_VPD_max_sm"
formula6<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_CMD + max_VPD_max_sm + mean_Eref + T_DD_0"
formula6a<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_CMD + max_VPD_max_sm + mean_Eref + T_DD_0 + mean_MAT"
formula7<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + T_VPD_max_sm + T_MAT + max_VPD_max_sm + mean_CMD + mean_MAT"
#formula7<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + max_VPD_max_sm  + mean_CMD + mean_MAT + mean_Eref"
cv_function(ht_tr,cv_folds_tr,formula5)
cv_function(ht_tr,cv_folds_tr,formula6)
cv_function(ht_tr,cv_folds_tr,formula6a)
cv_function(ht_tr,cv_folds_tr,formula7)
# prefer formula 7 here

# iterative search approach
test_vars_transfer <- colnames(ht_tr)[219:249]
test_vars_transfer <- colnames(ht_tr)[250:270]
test_vars_transfer <- colnames(ht_tr)[271:291]
test_vars_transfer <- colnames(ht_tr)[292:305]
tr_vars<-paste0(test_vars_transfer,collapse="+")

mods_tr <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) +  ht_pct.x + dt + Type + mean_CMD + mean_VPD_max_sm +",tr_vars)),
                             data=ht_tr,nvmax=20)
plot(mods_tr)

best_var_tr <- c("T_DD_0_sp","T_Tmax_wt","T_Tave_sp","T_Tave_sm","T_Tmax_at","T_Rad_sp","T_Rad_sm",
                 "T_PAS_wt","T_Eref_wt","T_DD18_sm","T_TD","T_NFFD","T_CMI","T_RH","T_EXT","T_EMT")

#best_var_tr <- c("T_DD_0_sp","T_Tmax_wt","T_Tave_sp","T_Tave_sm","T_Tmax_at","T_Rad_sp","T_Rad_sm",
#                 "T_NFFD_sm","T_RH_wt","T_CMI_wt","T_RH","T_CMI","T_EMT","T_CMD","T_MAP","T_NFFD","T_bFFP")
#best_var_tr <- c("T_RH","T_NFFD","T_DD18","T_DD5","T_Eref","T_NFFD_wt","T_CMI_wt",
#                 "T_NFFD_sm","T_DD_0_sp","T_Rad_at","T_Rad_wt","T_Tave_sm","T_Tmax_wt")


tr_vars_best<-paste0(best_var_tr,collapse="+")
mods_tr <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) +  ht_pct.x + dt + Elevation + mean_CMD + mean_VPD_max_sm + mean_Eref + mean_MAP+ mean_MAT + mean_CMI + T_MAT+",tr_vars_best)),
                             data=ht_tr,nvmax=20)
plot(mods_tr)

formula7<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + T_DD_0_sp + T_Tmax_wt + T_TD + T_MAT"
cv_function(ht_tr,cv_folds_tr,formula7)
# worse than previous results

# "full" iterative search ##

# need to re-do this
# 14-100- means
# 101-128- min/max
# 129-215 P variables
# 219 to 305 T variables
# 5 variables leaves us groups of 45
# 86 vars per group
# two groups of 43

test_vars_transfer <- colnames(ht_tr)[14:57]
# mean_DD_0_sp + mean_Rad_wt + mean_PPT_sm + mean_Tave_wt + mean_NFFD_sm + mean_DD_18_sm + 
# + mean_Tmax_at + mean_Tmax_sp + mean_Tave_sp
test_vars_transfer <- colnames(ht_tr)[58:100]
# mean_Eref_sp + mean_CMD_wt + mean_RH_sm + mean_RH_at + mean_MWMT + mean_DD5+ mean_DD_0 + mean_FFP +
#  mean_VPD_max_sm + mean_CMD + mean_Eref + mean_VPD_av_sm
test_vars_transfer <- colnames(ht_tr)[101:128]
# min_Rad_wt + max_NFFD_sm + max_MAT + max_VPD_max_sm + max_CMI_sm + max_Tmax_sm
test_vars_transfer <- colnames(ht_tr)[129:172]
# P_Tmin_sp + P_Tmax_wt + P_Tave_wp + P_DD_0_wt + P_DD_18_wt + P_NFFD_sm + P_NFFD_sp + P_Tmin_at
test_vars_transfer <- colnames(ht_tr)[173:215]
# P_Eref_at + P_bFFP + P_DD5 + P_MAT + P_RH_at + P_VPD_av_sm + P_CMD + P_DD_18 + P_Eref_sm + P_CMI_sm
# transfer vars:
# best_var_tr <- c("T_DD_0_sp","T_Tmax_wt","T_Tave_sp","T_Tave_sm","T_Tmax_at","T_Rad_sp","T_Rad_sm",
#"T_PAS_wt","T_Eref_wt","T_DD18_sm","T_TD","T_NFFD","T_CMI","T_RH","T_EXT","T_EMT")
# 15

mods_tr <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) +  ht_pct.x + dt + Type + ",tr_vars)),
                             data=ht_tr,nvmax=20)
plot(mods_tr)

test_vars<-c("mean_DD_0_sp","mean_Rad_wt","mean_PPT_sm","mean_Tave_wt","mean_NFFD_sm","mean_DD_18_sm","mean_Tmax_at","mean_Tmax_sp","mean_Tave_sp",
             "mean_Eref_sp","mean_CMD_wt","mean_RH_sm","mean_RH_at","mean_MWMT","mean_DD5","mean_DD_0","mean_FFP","mean_VPD_max_sm","mean_CMD","mean_Eref","mean_VPD_av_sm",
             "min_Rad_wt","max_NFFD_sm","max_MAT","max_VPD_max_sm","max_CMI_sm","max_Tmax_sm",
             #"P_Tmin_sp","P_Tmax_wt","P_Tave_wt","P_DD_0_wt","P_DD_18_wt","P_NFFD_sm","P_NFFD_sp","P_Tmin_at",
             #"P_Eref_at","P_bFFP","P_DD5","P_MAT","P_RH_at","P_VPD_av_sm","P_CMD","P_DD_18","P_Eref_sm","P_CMI_sm",
             "T_DD_0_sp","T_Tmax_wt","T_Tave_sp","T_Tave_sm","T_Tmax_at","T_Rad_sp","T_Rad_sm","T_PAS_wt","T_Eref_wt","T_DD18_sm","T_TD","T_NFFD","T_CMI","T_RH","T_EXT","T_EMT")
# 64, need to cut 14
test_vars<-c("P_Tmin_sp","P_Tmax_wt","P_Tave_wt","P_DD_0_wt","P_DD_18_wt","P_NFFD_sm","P_NFFD_sp","P_Tmin_at",
             "P_Eref_at","P_bFFP","P_DD5","P_MAT","P_RH_at","P_VPD_av_sm","P_CMD","P_DD_18","P_Eref_sm","P_CMI_sm",
             "mean_Tave_wt","mean_NFFD_sm","mean_DD_18_sm","mean_Tmax_at","mean_Tmax_sp","mean_Tave_sp",
             "mean_Eref_sp","mean_CMD_wt","mean_RH_sm","mean_DD5","mean_DD_0","mean_FFP","mean_CMD","max_MAT","max_VPD_max_sm",
             "max_Tmax_sm","T_DD_0_sp","T_Tmax_wt","T_DD18_sm")

tr_vars<-paste0(test_vars,collapse="+")
mods_tr_iter <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) +  ht_pct.x + dt + Type + ",tr_vars)),
                                  data=ht_tr,nvmax=20,really.big=T)
mods_tr_iter2 <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) +  ht_pct.x + dt + Type + ",tr_vars)),
                                   data=ht_tr,nvmax=30)
mods_tr_iter3 <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) +  ht_pct.x + dt + Type + ",tr_vars)),
                                   data=ht_tr,nvmax=30)
plot(mods_tr_iter3)

formula8a<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + T_DD_0_sp + T_Tmax_wt + mean_Tmax_at + mean_CMD + max_MAT + max_VPD_max_sm"
formula8b<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_CMD + max_MAT + max_VPD_max_sm + mean_Tmax_at"
formula8c<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_CMD + max_VPD_max_sm + T_Tmax_wt + P_NFFD_sm + P_Tmin_sp"
formula8<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_CMD + max_VPD_max_sm + T_DD_0_sp + mean_Tave_sp"

cv_function(ht_tr,cv_folds_tr,formula8)
# elastic net
#### GLMNET-based elastic net ####
library(doParallel)
registerDoParallel(10)
x<-ht_tr[,c(8,10,12,14:215,218:305)]
which(!complete.cases(x))
x<-x[ , colSums(is.na(x))==0]
x$logheight<-log(x$height.x)
x$Type <-ifelse(x$Type=="PSP",1,0)
x<-as.matrix(x)
#x$int<-x$DBH.x*x$BAL
y<- ihs(ht_tr$dh)
#cv.glmnet(x,y,alpha=0,lambda = c(0.01,0.05,0.1),nfolds=4,parallel=TRUE)
elastic1<-cv.glmnet(x,y,nfolds=10,parallel=TRUE)
elastic2<-cv.glmnet(x,y,foldid= index_to_folds(cv_folds_tr),parallel=TRUE)
registerDoSEQ()

plot(elastic1)
plot(elastic2)
elastic1$lambda.1se
elastic2$lambda.1se
# might be worth noting that this approach suggests ~200 vars if taken at face value, and that's only based on a flawed 1SE interpretation
# error is still shrinking at all variables, this pretty much lacks robustness
# elbow is at ~16-22 variables, though
# looks like we're around 50 at ~exp(-5)
rownames(coef(elastic1, s=0.02,))[coef(elastic1, s = 0.02)[,1]!= 0]
rownames(coef(elastic2, s=0.038,))[coef(elastic2, s = 0.038)[,1]!= 0]
rownames(coef(elastic1, s=0.008,))[coef(elastic1, s = 0.008)[,1]!= 0]
rownames(coef(elastic1, s=0.05,))[coef(elastic1, s = 0.0045)[,1]!= 0]
# potential "just elastic net" model:
rownames(coef(elastic1, s=0.042,))[coef(elastic1, s = 0.042)[,1]!= 0]
rownames(coef(elastic2, s=0.04634947,))[coef(elastic2, s = 0.04634947)[,1]!= 0]
#formula8<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_NFFD_sm + max_CMI_sm + P_NFFD_sm + P_TD + P_bFFP"
formula8<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + max_DD_0_sm + max_CMI_sm + P_Tmin_sp + P_NFFD_sm"
cv_function(ht_tr,cv_folds_tr,formula8)
formula8a<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + max_DD_0_sm + max_CMI_sm + P_Tmin_sp + P_NFFD_sm + T_RH"
cv_function(ht_tr,cv_folds_tr,formula8a)
# mean_PAS_sp- 3.009
# T_DD18_sp - 3.0096
# T_RH_sp- 3.012
# T_RH - 3.0096

# could feed these back into best subsets and see if we can get a parsimonious model
elastic_vars <- rownames(coef(elastic_noel, s=0.0045,))[coef(elastic_noel, s = 0.0045)[,1]!= 0][c(5:40,42:46)]
#elastic_vars <- rownames(coef(elastic1, s=0.0045,))[coef(elastic1, s = 0.0045)[,1]!= 0][5:46]
tr_vars_en<-paste0(elastic_vars,collapse="+")
mods_tr_en <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + Type + log(height.x) +  ht_pct.x + dt +", tr_vars_en)),
                                data=ht_tr,nvmax=20)
plot(mods_tr_en)


formula10<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_NFFD_sm + T_DD_0_sm + P_RH_at + P_Eref_at + P_Tmin_at"
cv_function(ht_tr,cv_folds_tr,formula10)
formula10<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_Tmax_at +mean_Eref_at + mean_CMI_sm + max_VPD_max_sm + T_Eref_sp"

# formula 9 is the best overall CV performer, without regularization.

# I want to try a naive stepwise approach too
#14-215,219-305

# big issues here is radiation variables violate an assumption, we probably need to ensure these don't get selected by taking them out of the pool
# what are indicies: 
rad_ind<-which(str_detect(colnames(ht_tr),"[Rr]ad"))
ind<-c(14:215,219:305)
ind<-ind[! ind %in% rad_ind]

test_vars <- colnames(ht_tr)[ind]
tr_vars<-paste0(test_vars,collapse="+")
mods_step <- leaps::regsubsets(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) + ht_pct.x + dt + Type +",tr_vars)),data=ht_tr,nvmax=50,method="forward")
full_model <- lm(as.formula(paste("ihs(dh)  ~ height.x + log(height.x) + ht_pct.x + dt + Type +",tr_vars)),data=ht_tr)
start_model <- lm(ihs(dh)  ~ height.x + log(height.x) + ht_pct.x + dt + Type,data=ht_tr)
mods_step_v2 <- MASS::stepAIC(start_model,scope=list(lower=start_model,upper=full_model),direction="both")
#ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + T_DD18_sm + 
#mean_Tmax_at + max_MAT + mean_NFFD + mean_MAT + max_NFFD_sm
plot(mods_step)
coef(mods_step_v2,5)
# formula11<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_Tmax_at + mean_MAT + mean_NFFD + max_CMI + T_DD18_sm"
# cv_function(ht_tr,cv_folds_tr,formula11)
# formula12<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_Tmax_at + mean_MAT + mean_NFFD + mean_MAT + T_DD18_sm"
# cv_function(ht_tr,cv_folds_tr,formula12)
# Elevation + P_Tmin_sp + max_VPD_max_sm + mean_CMD + T_Tmax_wt
# c("height.x" , "ht_pct.x" , "dt" , "Type","mean_Tmax_at","mean_NFFD","mean_MAT","T_DD18_sm","max_CMI")
formula11<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_Tmax_at + mean_NFFD + mean_MAT + T_DD18_sm + max_CMI"
cv_function(ht_tr,cv_folds_tr,formula11)


# new, verifiable models for testing:
# model 3, original "guided" subsets
#ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type +  T_VPD_max_sm + T_MAT + max_VPD_max_sm + mean_CMD + mean_MAT"
# vars1<-c("height.x" , "ht_pct.x" , "dt" , "Type" ,  "T_VPD_max_sm" , "T_MAT" , "max_VPD_max_sm" , "mean_CMD" ,  "mean_MAT")


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

cv_function_reg(ht_tr,cv_folds_tr,vars1,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars2,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars3,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars4,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars5,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars6,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars7,0.01)/10

folds_tr <- index_to_folds(cv_folds_tr)
test_seq<-exp(seq(from=-7,to=1,by=0.1))
search_tr<-cv_reg_search2(ht_tr,vars1,folds_tr,n=10,test_seq)
search_tr<-cv_reg_search2(ht_tr,vars7,folds_tr,n=10,test_seq)

plot(search_tr)
search_tr$lambda.1se
search_tr$lambda.min

# 1 0.005/0.22
# 2 0.006/0.22
# 3 0.006/0.22
# 4 0.01/0.27
# 5 0.005/0.22
# 6 0.01/0.27
# 7 0.006/0.22

# at this point I think I really need an independent validation set
# Probably prince george geno trial
pg_merge<-readRDS(file="./Data/PG_Trial/pg_dat.RDS")

ridge_error(ht_tr,pg_merge,vars1,0.22)
ridge_error(ht_tr,pg_merge,vars2,0.24)
ridge_error(ht_tr,pg_merge,vars3,0.20)
ridge_error(ht_tr,pg_merge,vars4,0.27)
ridge_error(ht_tr,pg_merge,vars5,0.2)

# P_Ave_sp - doesn't really make sense
# max_VPD_max_sm mean_CMD mean_Eref_sm T_RH_wt mean_Tmax_at mean_DD_0_at
# folds for cv:
cv_folds_tr <- cv_climate(ht_tr,"mean_Tave_sm")

# # these are from selecting likely candidate vars, <50, then running through best subsets, then using CV
# vars1<-c("height.x" , "ht_pct.x" , "dt" , "Type" ,  "mean_MSP" , "T_DD_0_sp" , "T_Tmax_wt" , "T_RH" ,  "mean_CMD" , "mean_VPD_max_sm" , "mean_Eref" ,"T_DD1040" ,"mean_MAT")
# # Think this is full dataset elastic net, then into best subsets
# vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type","mean_PAS_sp","mean_EXT","mean_CMD","max_DD_0_sm",
#          "max_CMI_sm" , "max_VPD_max_sm" , "P_Tave_sp" , "P_NFFD_sm" , "P_PAS_sp" , "P_CMD_at" , "P_RH_at", "P_EXT")
# vars3<-c("height.x" , "ht_pct.x" , "dt" , "Type", "Elevation","mean_VPD_max_sm","mean_CMD","mean_DD_0_at","mean_Eref_sm","T_RH_wt")
# vars4<-c("height.x" , "ht_pct.x" , "dt" , "Type", "T_DD_0_sp","T_RH","T_Tmax_wt","T_NFFD","T_Tave_sm","T_NFFD_wt")
# vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","max_VPD_max_sm",
#          "mean_CMD","mean_Eref_sm","T_RH_wt","mean_Tmax_at","mean_DD_0_at")
# # max_VPD_max_sm mean_CMD mean_Eref_sm T_RH_wt mean_Tmax_at mean_DD_0_at


#cv_function_reg(ht_tr,cv_folds_tr,vars1,0.001)
cv_function_reg(ht_tr,cv_folds_tr,vars2,0.01)
cv_function_reg(ht_tr,cv_folds_tr,vars3,0.01)
cv_function_reg(ht_tr,cv_folds_tr,vars4,0.01)
cv_function_reg(ht_tr,cv_folds_tr,vars5,0.01)
cv_function_reg(ht_tr,cv_folds_tr,vars6,0.01)
cv_function_reg(ht_tr,cv_folds_tr,vars7,0.01)

folds_tr <- index_to_folds(cv_folds_tr)

search_tr<-cv_reg_search(ht_tr,vars4,folds_tr,n=10)
plot(search_tr)
search_tr$lambda.1se
search_tr$lambda.min


# fix lambda sequence
test_seq<-exp(seq(from=-7,to=1,by=0.1))
folds_tr<-index_to_folds(cv_folds_tr)
search_tr<-cv_reg_search2(ht_tr,vars1,folds_tr,n=10,test_seq)
search_tr<-cv_reg_search2(ht_tr,vars5,folds_tr,n=10,test_seq)

plot(search_tr)
search_tr$lambda.1se
search_tr$lambda.min
# min/1SE
# 1- 0.002/0.18
# 2- 0.001/0.272
# 3 - 0.001/0.20
# 4 - 0.007/0.27
# 5 - 0.0016/0.20

#0.2725318
#0.2018965
# ridge_error <- function(x_tr,x_val,vars,lambda)
ridge_error(ht_tr,ht_val,vars1,0.01)
ridge_error(ht_tr,ht_val,vars2,0.01)
ridge_error(ht_tr,ht_val,vars3,0.01)
ridge_error(ht_tr,ht_val,vars4,0.01)
ridge_error(ht_tr,ht_val,vars5,0.01)
ridge_error(ht_tr,ht_val,vars6,0.01)
ridge_error(ht_tr,ht_val,vars7,0.01)

ridge_error(ht_tr,ht_val,vars1,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars2,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars3,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars4,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars5,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars6,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars7,0.01,type="bias")

# at this point I think I really need an independent validation set
# Probably prince george geno trial
pg_merge<-readRDS(file="./Data/PG_Trial/pg_dat.RDS")

ridge_error(ht_tr,pg_merge,vars1,0.18)
ridge_error(ht_tr,pg_merge,vars2,0.27)
ridge_error(ht_tr,pg_merge,vars3,0.20)
ridge_error(ht_tr,pg_merge,vars4,0.27)
ridge_error(ht_tr,pg_merge,vars5,0.20)

#ridge_error(ht_tr,pg_merge,vars1,0.01)
ridge_error(ht_tr,pg_merge,vars2,0.01)
ridge_error(ht_tr,pg_merge,vars3,0.01)
ridge_error(ht_tr,pg_merge,vars4,0.01)
ridge_error(ht_tr,pg_merge,vars5,0.01)
ridge_error(ht_tr,pg_merge,vars6,0.01)
ridge_error(ht_tr,pg_merge,vars7,0.01)

ridge_error(ht_tr,pg_merge,vars1,0.01,type="bias")
ridge_error(ht_tr,pg_merge,vars2,0.01,type="bias")
ridge_error(ht_tr,pg_merge,vars3,0.01,type="bias")
ridge_error(ht_tr,pg_merge,vars4,0.01,type="bias")
ridge_error(ht_tr,pg_merge,vars5,0.01,type="bias")
ridge_error(ht_tr,pg_merge,vars6,0.01,type="bias")
ridge_error(ht_tr,pg_merge,vars7,0.01,type="bias")

# observations- including elevation kills P_ variables- so they're probably doing a lot of
# work in PSP data

# variance estimates
var(ihs(ht_tr$dh))
var(ihs(ht_val$dh))
var(ihs(pg_merge[!is.na(pg_merge$dh),]$dh))

# parameters
model2<-ridge_params(ht_tr,vars2,0.01)
model3<-ridge_params(ht_tr,vars3,0.01)
model4<-ridge_params(ht_tr,vars4,0.01)
model5<-ridge_params(ht_tr,vars5,0.01)
model6<-ridge_params(ht_tr,vars6,0.01)
model7<-ridge_params(ht_tr,vars7,0.01)
model2$beta
model3$beta
model4$beta
model5$beta
model6$beta
model7$beta

#### what happens if Type gets removed? ####
#vars6b<-c("height.x" , "ht_pct.x" , "dt" , "Elevation", "max_CMI_sm","P_NFFD_sm")
ridge_params_b <- function(x_tr,vars,lambda){
  x<-x_tr[,vars]
  x <- subset(x,select= -Type)
  x$logheight<-log(x$height.x)
  y<- ihs(x_tr$dh)
  mod_temp <- glmnet(x,y,alpha=0,lambda = lambda)
  return(mod_temp)
}
# this is a 1-1 replacement for the x_clean function: it will change results of all validation functions to exclude Type
# use carefully, and reset to dis_functions.R after use!
# x_clean<-function(x_tr,vars){
#   x<-x_tr[,vars]
#   #x$Type <-ifelse(x$Type=="PSP",1,0)
#   x <- subset(x,select= -Type)
#   x$logheight<-log(x$height.x)
#   y<- ihs(x_tr$dh)
#   na.list<-which(is.na(y))
#   na.list<-which(is.na(y))
#   if(any(na.list)){
#     x<-x[-na.list,]
#     y<-y[-na.list]}
#   return(list(x,y))
# }

# model6b<-ridge_params_b(ht_tr,vars6,0.01)
# model6b$beta
# model6$beta
# 
# vars5b<-c("height.x" , "ht_pct.x" , "dt" ,"Elevation","mean_CMD","max_VPD_max_sm","mean_DD_0_at","P_DD_0_sp")
# model5b<-ridge_params_b(ht_tr,vars5b,0.01)
# model5b$beta
# model5$beta

model7<-ridge_params_b(ht_tr,vars7,0.01)
model2<-ridge_params_b(ht_tr,vars2,0.01)
model3<-ridge_params_b(ht_tr,vars3,0.01)
model4<-ridge_params_b(ht_tr,vars5,0.01)
model6<-ridge_params_b(ht_tr,vars6,0.01)
model5<-ridge_params_b(ht_tr,vars4,0.01)
model7$beta
model2$beta
model3$beta
model4$beta
model6$beta
model5$beta

