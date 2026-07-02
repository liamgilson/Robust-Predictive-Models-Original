#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 41- Code for results tables for manscript #


#### Results Tables Compilation ####
# May 11, 2026
# Liam Gilson


#### Models (no elev) ####

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


#### Test-PSP (qc MAG data) ####

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

#### Validation error ####

#cv_function_reg(ht_tr,cv_folds_tr,vars1,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars2,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars3,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars4,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars5,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars6,0.01)/10
cv_function_reg(ht_tr,cv_folds_tr,vars7,0.01)/10

#ridge_error(ht_tr,ht_val,vars1,0.01)
ridge_error(ht_tr,ht_val,vars2,0.01)
ridge_error(ht_tr,ht_val,vars3,0.01)
ridge_error(ht_tr,ht_val,vars4,0.01)
ridge_error(ht_tr,ht_val,vars5,0.01)
ridge_error(ht_tr,ht_val,vars6,0.01)
ridge_error(ht_tr,ht_val,vars7,0.01)

#ridge_error(ht_tr,ht_val,vars1,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars2,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars3,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars4,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars5,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars6,0.01,type="bias")
ridge_error(ht_tr,ht_val,vars7,0.01,type="bias")

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


#### Graphics ####
custom_pallete <- c(wesanderson::wes_palette("Moonrise2")[2], wesanderson::wes_palette("Moonrise3")[c(5,4,3,2,1)])

c2_res2 <- read.csv("./Data/ch2_results_2.csv")

b<-c2_res2 %>% filter(Type=="bias") %>% mutate(Bias=abs(Result),Method=factor(Model,levels=
                                                                               c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
                                              Setting = factor(Setting,levels=c("Validation","Test-Provenance","Test-PSP"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Setting,y=Bias,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="Bias by Validation Setting",x="Validation Setting",y="Bias (m)")+scale_fill_manual(values=custom_pallete)

a<-c2_res2 %>% filter(Type=="MSE",Setting!="CV") %>% mutate(Method=factor(Model,levels=
                                                                                                     c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
                                             Setting = factor(Setting,levels=c("Validation","Test-Provenance","Test-PSP"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Setting,y=Result,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="MSE by Validation Setting",x="Validation Setting",y=expression(Mean~Squared~Error~(m^2)))+scale_fill_manual(values=custom_pallete)
library(wesanderson)
custom_pallete <- c(wes_palette("Moonrise2")[2], wes_palette("Moonrise3")[c(5,4,3,2,1)])
ggpubr::ggarrange(a,b,nrow=1,common.legend = TRUE,legend="bottom")
grid.arrange(a,b,nrow=1)
pdf("height_res1_v2.pdf",width=8,height=4)
dev.off()
