#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 40- Code for Manuscript Graphics #


## General Plotting Packages ##

library(ggplot2)
library(gridExtra)
library(ggpubr)

#### Data read in (if used independently)
#saveRDS(ht_tr,file="ht_tr_graphics.RDS")
#saveRDS(faib_pred,file="faib_pred.RDS")
#saveRDS(pg_pred,file="pg_pred.RDS")


####


samp1<-MASS::mvrnorm(100,c(1,1),Sigma=matrix(c(1,0.3,0.3,1),nrow=2,ncol=2,byrow=T))
samp1<-data.frame(x=samp1[,1],y=samp1[,2],set="Training")

samp2<-MASS::mvrnorm(100,c(7,7),Sigma=matrix(c(1,0.1,0.1,1),nrow=2,ncol=2,byrow=T))
samp2<-data.frame(x=samp2[,1],y=samp2[,2],set="Test")

samp<-rbind(samp1,samp2)

temp1 <- data.frame(int=summary(lm(y ~ x, data=samp1))$coef[1,1],slope=summary(lm(y ~ x, data=samp1))$coef[2,1],model="Training Only")
temp2 <- data.frame(int=summary(lm(y ~ x, data=samp))$coef[1,1],slope=summary(lm(y ~ x, data=samp))$coef[2,1],model="All Data")
temp3 <-  data.frame(int=(summary(lm(y ~ x, data=samp1))$coef[1,1]+summary(lm(y ~ x, data=samp))$coef[1,1])/2,
                     slope=(summary(lm(y ~ x, data=samp1))$coef[2,1]+summary(lm(y ~ x, data=samp))$coef[2,1])/2,model="Robust Model")
mods<-rbind(temp1,temp2,temp3)

samp %>% ggplot(aes(x=x,y=y,color=set))+geom_point()+geom_abline(data=mods,aes(intercept=int,slope=slope,linetype=model),size=1)


pdf("example_1.pdf",width=5,height=3.5)
samp %>% ggplot(aes(x=x,y=y,color=set))+geom_point()+geom_abline(data=mods,aes(intercept=int,slope=slope,linetype=model),size=1)
dev.off()

tiff("example_1.tiff", width = 5, height = 3.5, units = 'in', res = 500)
samp %>% ggplot(aes(x=x,y=y,color=set))+geom_point()+geom_abline(data=mods,aes(intercept=int,slope=slope,linetype=model),size=1)
dev.off()


#### Climate Envelope Graphic ####

pca_dat <- ht_dat %>% select(source,Type, starts_with("mean_"))
pca_dat <- unique(pca_dat)

pca_dat_clean <- within(pca_dat,rm("mean_MAR","mean_Rad_wt","mean_Rad_sp","mean_Rad_sm","mean_Rad_at"))
pca_dat_clean[which(rowSums(is.na(pca_dat_clean[,-1]))>0),]
pca_dat_clean <- pca_dat_clean[complete.cases(pca_dat_clean),]
pca_res <- prcomp(pca_dat_clean[,-c(1:2)],scale=TRUE)
autoplot(pca_res,data=pca_dat,color="Type")
pca_res_plot<-cbind(pca_dat_clean[,1:2],pca_res$x)

p1<-pca_res_plot %>% arrange(Type) %>% ggplot(aes(x=PC1,y=PC2,color=Type))+geom_point()+
  labs(title="Site Climate")+scale_color_manual(values=wesanderson::wes_palette("FantasticFox1")[c(3,1)],labels=c("Provenance","PSP"))


pca_dat_t <- ht_dat %>% select(source,Type, starts_with("T_"))
pca_dat_t <- unique(pca_dat_t)

pca_dat_t_clean <- within(pca_dat_t,rm("T_MAR","T_Rad_wt","T_Rad_sp","T_Rad_sm","T_Rad_at"))
pca_dat_t_clean[which(rowSums(is.na(pca_dat_t_clean[,-1]))>0),]
pca_dat_t_clean <- pca_dat_t_clean[complete.cases(pca_dat_t_clean),]
pca_res_t <- prcomp(pca_dat_t_clean[,-c(1:2)],scale=TRUE)
pca_res_t_plot<-cbind(pca_dat_t_clean[,1:2],pca_res_t$x)

p2<-pca_res_t_plot %>% arrange(desc(Type)) %>% ggplot(aes(x=PC1,y=PC2,color=Type))+geom_point()+
  labs(title="Transfer Climate")+scale_color_manual(values=wesanderson::wes_palette("FantasticFox1")[c(3,1)],labels=c("Provenance","PSP"))

ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend="bottom")

#tiff("coverage.tiff", width = 7, height = 3.5, units = 'in', res = 500)
pdf("coverage_v2.pdf", width = 7, height = 3.5)
ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend="bottom")
dev.off()

#### code for various numbers and statistics in data description ####

# counts by type:
ht_dat_na %>% group_by(Type) %>% count()
ht_dat %>% filter(Type=="PSP") %>% summarize(min_year = min(year.x, na.rm=T),max_year= max(year.y))

ht_dat %>% filter(Type=="Prov") %>% group_by(source) %>% summarize(provenances = length(unique(Prov)))
ht_dat %>% filter(Type=="Prov") %>% group_by(source) %>% summarize(provenances = length(unique(Site)))

ht_dat %>% filter(Type=="Prov") %>% summarize(min_year = min(year.x, na.rm=T),max_year= max(year.y))

# counts
sum(ht_dat_na[ht_dat_na$Type=="Prov",]$dh <= 0)
4227/287175
#### Prediction Plots ####

formula_psp<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_NFFD_sp + P_Eref_at + P_RH_wt + T_CMI_wt"

vars1<-c("height.x" , "ht_pct.x" , "dt" , "Type" ,  "mean_NFFD_sp" , "P_Eref_at" , "P_RH_wt" , "T_CMI_wt")

vars2<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_Eref_sm","max_VPD_max_sm","mean_CMD","P_Tave_sp","T_RH_wt","mean_Tmax_at")

vars3<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_NFFD_sm", "mean_PAS_sp", "mean_EXT",
         "max_DD_0_sm", "max_CMI_sm", "P_NFFD_sm", "P_PAS_sp", "P_CMD_at", "P_RH_at", "P_bFFP", "P_EXT", "T_RH_at")

mod1 <-predict_pg_reg(pg_pred,"PSP Only",vars1,lambda1=0.2,n=16)
mod2 <-predict_pg_reg(pg_pred,"2-Step",vars2,lambda1=0.2,n=16)
mod3 <-predict_pg_reg(pg_pred,"Elastic Net",vars3,lambda1=0.2,n=16)

mods<-rbind(mod1,mod2,mod3)
mod_means<- mods %>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))


a<- mod_means %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()

mod1 <-predict_pg_reg(pg_pred,"PSP Only",vars1,lambda1=0,n=16)
mod2 <-predict_pg_reg(pg_pred,"2-Step",vars2,lambda1=0,n=16)
mod3 <-predict_pg_reg(pg_pred,"Elastic Net",vars3,lambda1=0,n=16)

mods<-rbind(mod1,mod2,mod3)
mod_means<- mods %>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))

b<- mod_means %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()

ggarrange(a,b,nrow=1)

#### Predictions: All models ####

# new, verifiable models for testing:

# vars3<-c("height.x" , "ht_pct.x" , "dt" , "Type", "T_DD_0_sp" , "T_Tmax_wt" , "T_TD" , "T_MAT")
# # elastic net fed back into best subsets
# #formula9<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + Elevation + mean_Eref_sm + max_VPD_max_sm + mean_CMD+P_Tave_sp+mean_Tmax_at +T_RH_wt"
# vars5<-c("height.x" , "ht_pct.x" , "dt" , "Type","Elevation","mean_Eref_sm","max_VPD_max_sm","mean_CMD","P_Tave_sp","T_RH_wt","mean_Tmax_at")
# # "optimal" 1SE model from elastic net with blocked folds:
# # formula8a<-"ihs(dh) ~ height.x + log(height.x) + ht_pct.x + dt + Type + mean_NFFD_sm +  mean_PAS_sp + mean_EXT +
# #max_DD_0_sm + max_CMI_sm + P_NFFD_sm + P_PAS_sp + P_CMD_at + P_RH_at + P_bFFP + P_EXT + T_RH_at"
# vars6<-c("height.x" , "ht_pct.x" , "dt" , "Type", "mean_NFFD_sm", "mean_PAS_sp", "mean_EXT",
#          "max_DD_0_sm", "max_CMI_sm", "P_NFFD_sm", "P_PAS_sp", "P_CMD_at", "P_RH_at", "P_bFFP", "P_EXT", "T_RH_at")
# # stepwise variables
# vars7<-c("height.x" , "ht_pct.x" , "dt" , "Type","mean_Tmax_at","mean_MAT","mean_NFFD","mean_MAT","T_DD18_sm")


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

mod2 <-predict_pg_reg(pg_pred,"CV-Subsets",vars2,lambda1=0.01,n=16)
mod3 <-predict_pg_reg(pg_pred,"CV-Subsets (iter.)",vars3,lambda1=0.01,n=16)
mod4 <-predict_pg_reg(pg_pred,"Elastic Net",vars4,lambda1=0.01,n=16)
mod5 <-predict_pg_reg(pg_pred,"CV-Subsets (EN)",vars5,lambda1=0.01,n=16)
mod6 <-predict_pg_reg(pg_pred,"Elastic Net-1SE",vars6,lambda1=0.01,n=16)
mod7 <-predict_pg_reg(pg_pred,"Step AIC",vars7,lambda1=0.01,n=16)
#mods<-rbind(mod4,mod5,mod6,mod7)
mods<-rbind(mod2,mod3,mod4,mod5,mod6,mod7)
mod_means<- mods %>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))


# setup code for plots:
custom_pallete <- c(wesanderson::wes_palette("Moonrise2")[2], wesanderson::wes_palette("Moonrise3")[c(5,4,3,2,1)])
#mutate(Method=factor(Model,levels=c("Step-AIC","CV-Subsets","CV-Subsets-iter","EN-Subsets","Elastic-Net","EN-1SE")))

a<-mod_means %>% 
  mutate(Method=factor(Model,levels=c("Step AIC","CV-Subsets","CV-Subsets (iter.)","CV-Subsets (EN)","Elastic Net","Elastic Net-1SE"))) %>%
  ggplot(aes(x=year.y,y=mean_height,color=Method))+geom_line(linewidth=1)+labs(title="Church Road (PG Trial)",x="Year",y="Mean Height (m)")+
    scale_color_manual(values=custom_pallete)
 # vanilla colours version
a<-mod_means %>% 
  mutate(Method=factor(Model,levels=c("Step AIC","CV-Subsets","CV-Subsets (iter.)","CV-Subsets (EN)","Elastic Net","EN-1SE"))) %>%
  ggplot(aes(x=year.y,y=mean_height,color=Method))+geom_line(linewidth=1)+labs(title="Church Road (PG Trial)",x="Year",y="Mean Height (m)")

faib_pred<-readRDS(file="faib_pred.RDS")
mod2f <-predict_pg_reg(faib_pred,"CV-Subsets",vars2,lambda1=0.01,n=17)
mod3f <-predict_pg_reg(faib_pred,"CV-Subsets (iter.)",vars3,lambda1=0.01,n=17)
mod4f <-predict_pg_reg(faib_pred,"Elastic Net",vars4,lambda1=0.01,n=17)
mod5f <-predict_pg_reg(faib_pred,"CV-Subsets (EN)",vars5,lambda1=0.01,n=17)
mod6f <-predict_pg_reg(faib_pred,"Elastic Net-1SE",vars6,lambda1=0.01,n=17)
#mod7f <-predict_pg_rni(faib_pred,"EN-NI",vars6,lambda1=0.23,n=17)
mod7f <-predict_pg_reg(faib_pred,"Step AIC",vars7,lambda1=0.01,n=17)
#modsf<-rbind(mod4f,mod5f,mod6f,mod7f)
modsf<-rbind(mod2f,mod3f,mod4f,mod5f,mod6f,mod7f)
mod_meansf<- modsf%>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))
b<-mod_meansf %>% mutate(Method=factor(Model,levels=c("Step AIC","CV-Subsets","CV-Subsets (iter.)","CV-Subsets (EN)","Elastic Net","Elastic Net-1SE"))) %>%
  ggplot(aes(x=year.y,y=mean_height,color=Method))+geom_line(linewidth=1)+labs(title="Plot 59144 (FAIB)",x="Year",y="Mean Height (m)")+
  scale_color_manual(values=custom_pallete)
# vanilla colours version
b<-mod_meansf %>% mutate(Method=factor(Model,levels=c("Step AIC","CV-Subsets","CV-Subsets (iter.)","CV-Subsets (EN)","Elastic Net","EN-1SE"))) %>%
  ggplot(aes(x=year.y,y=mean_height,color=Method))+geom_line(linewidth=1)+labs(title="Plot 59144 (FAIB)",x="Year",y="Mean Height (m)")

ggarrange(a,b,nrow=1,common.legend=TRUE,legend="bottom")


pdf("fig1_v5.pdf", width = 7, height = 3.5)
ggarrange(a,b,nrow=1,common.legend=TRUE,legend="bottom")
dev.off()


#### Different Levels of lambda

modl1 <-predict_pg_reg(faib_pred,"0",vars4,lambda1=0.0,n=17)
modl2 <-predict_pg_reg(faib_pred,"0.01",vars4,lambda1=0.01,n=17)
modl3 <-predict_pg_reg(faib_pred,"0.22",vars4,lambda1=0.22,n=17)
#mod7f <-predict_pg_rni(faib_pred,"EN-NI",vars6,lambda1=0.23,n=17)
modl4 <-predict_pg_reg(faib_pred,"1",vars4,lambda1=1,n=17)
modsl<-rbind(modl1,modl2,modl3,modl4)
mod_meansl<- modsl%>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))
mod_meansl %>% ggplot(aes(x=year.y,y=mean_height,linetype=Model,color=Model))+geom_line()+
  labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)",linetype="Lambda")+
  theme(legend.position = "bottom")
mod_means_l2<-rbind(mod_meansl,gar_pred(68.9,0.01960819))

mod_means_l2 %>% mutate(empirical=as.logical(ifelse(Model=="Hu & Garcia",1,0))) %>%ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()+
  labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)",linetype="Lambda")+
  theme(legend.position = "bottom")

pdf("fig_lambda_v3.pdf", width = 5, height = 3.5)
mod_means_l2  %>%ggplot(aes(x=year.y,y=mean_height,linetype=Model,color=Model))+geom_line()+
  labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)",linetype="Lambda",color="Lambda")+
  scale_color_manual(values=c("black","black","black","black","red"))+
  theme(legend.position = "bottom")
dev.off()

mod_means_l2$Regularization <- case_when(
  mod_means_l2$Model=="0" ~"No regularization",
  mod_means_l2$Model=="0.01" ~"Minimum error",
  mod_means_l2$Model=="0.22" ~"1-SE error",
  mod_means_l2$Model=="1" ~"Over regularized",
  mod_means_l2$Model=="Hu & Garcia" ~"Hu & Garcia"
)

pdf("fig_lambda_v5.pdf", width = 5, height = 3.5)
mod_means_l2  %>%ggplot(aes(x=year.y,y=mean_height,color=Regularization))+geom_line()+
  labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)",linetype="Lambda",color="Lambda")+
  scale_color_manual(values=c("#440154FF","red","#3B528BFF","#21908CFF","#5DC863FF"))+
  theme(legend.position = "bottom")
dev.off()
png("fig_lambda_v3.png", width = 5, height = 3.5,res=300, units="in")
# color stuff:
#install.packages("viridis")
viridis::viridis(5)
#### Show error rates in FIA data ####

c<-ht_dat_na %>% filter(source=="FIA") %>% ggplot(aes(x=dh))+geom_histogram()+ xlim(-7,10)+
  geom_vline(xintercept = 0,colour="red")+labs(x="Height Increment",y="Count",title="Distribution of Increments in FIA Data")

pdf("fig3.pdf", width = 4, height = 4)
c
dev.off()



#### Height Growth Empirical EQ from Hu and Garcia ####

ht<-function(x,t){(283.9*x^0.5137)*(1-(1-(1.37/(283.9*x^0.5137))^0.5829)*exp(-x*(t-0.5)))^(1.71556)}
(283.9*x^0.5137)

# 62, 23.2

ht_solve<-function(x){(283.9*x^0.5137)*(1-(1-(1.37/(283.9*x^0.5137))^0.5829)*exp(-x*(62-0.5)))^(1.71556)}-23.2
install.packages("nleqslv")
library(nleqslv)
nleqslv(10,ht_solve)
0.01998504
# means from stand
ht_solve<-function(x){(283.9*x^0.5137)*(1-(1-(1.37/(283.9*x^0.5137))^0.5829)*exp(-x*(68.9-0.5)))^(1.71556)}-24.4
0.01960819

gar_pred <- function(t1,x){
  years<-seq(from=2020,to=2100,by=5)
  heights<-rep(NA,17)
  for (i in 1:17){
    t<-(i-1)*5+t1
    heights[i]<-ht(x,t)
  }
  res<-data.frame(Model="Hu & Garcia",year.y=years,mean_height=heights)
  return(res)
}
gar_pred(62,0.01998504)

# for lambda test:
# change labels for pres:
modl1$Model <- "No Reg."
modl2$Model <- "Min Error"
modl3$Model <- "1-SE"
modl4$Model <- "High Reg."
modsl<-rbind(modl1,modl2,modl3,modl4)
mod_meansl<- modsl%>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))
mod_means_l2<-rbind(mod_meansl,gar_pred(62,0.01998504))
mod_means_l2<-rbind(mod_meansl,gar_pred(68.9,0.01960819))
mod_means_l2 %>% ggplot(aes(x=year.y,y=mean_height,linetype=Model))+geom_line()+
  labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)",linetype="Lambda")+
  theme(legend.position = "bottom")

mod_means_l2 %>% ggplot(aes(x=year.y,y=mean_height,colour=Model))+geom_line()+
  labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)",colour="Lambda")+
  theme(legend.position = "bottom")+scale_color_hue(breaks=c("No Reg.","Min Error","1-SE","High Reg.","Hu & Garcia"))

pdf("fig_vic1.pdf", width = 7, height = 3.5)

dev.off()


mod4f <-predict_pg_reg(faib_pred,"Elastic Net",vars4,lambda1=0.23,n=17)
mod5f <-predict_pg_reg(faib_pred,"Subsets",vars5,lambda1=0.23,n=17)
mod6f <-predict_pg_reg(faib_pred,"EN-1SE",vars6,lambda1=0.23,n=17)
#mod7f <-predict_pg_rni(faib_pred,"EN-NI",vars6,lambda1=0.23,n=17)
mod7f <-predict_pg_reg(faib_pred,"Step AIC",vars7,lambda1=0.23,n=17)
modsf<-rbind(mod4f,mod5f,mod7f)
mod_meansf<- modsf%>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))
mod_meansf2<-rbind(mod_meansf,gar_pred(68.9,0.01960819))
mod_meansf2 %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()+labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)")

# victoria custom plot:

mod4 <-predict_pg_reg(pg_pred,"Elastic Net",vars4,lambda1=0.01,n=16)
mod5 <-predict_pg_reg(pg_pred,"Subsets",vars5,lambda1=0.01,n=16)
#mod6 <-predict_pg_reg(pg_pred,"EN-1SE",vars6,lambda1=0.01,n=16)
mod7 <-predict_pg_reg(pg_pred,"Step AIC",vars7,lambda1=0.01,n=16)
mods<-rbind(mod4,mod5,mod7)
mod_means<- mods %>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))


a<-mod_means %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()+labs(title="Prince George, SSP585",x="Year",y="Mean Height (m)")


mod4f <-predict_pg_reg(faib_pred,"Elastic Net",vars4,lambda1=0.01,n=17)
mod5f <-predict_pg_reg(faib_pred,"Subsets",vars5,lambda1=0.01,n=17)
#mod6f <-predict_pg_reg(faib_pred,"EN-1SE",vars6,lambda1=0.01,n=17)
#mod7f <-predict_pg_rni(faib_pred,"EN-NI",vars6,lambda1=0.23,n=17)
mod7f <-predict_pg_reg(faib_pred,"Step AIC",vars7,lambda1=0.01,n=17)
modsf<-rbind(mod4f,mod5f,mod7f)
mod_meansf<- modsf%>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))
b<-mod_meansf %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()+labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)")

ggarrange(a,b,nrow=1,common.legend=TRUE,legend="bottom")


pdf("fig_vic2.pdf", width = 7, height = 3.5)
ggarrange(a,b,nrow=1,common.legend=TRUE,legend="bottom")
dev.off()


# figure 1, manuscript version

# solve for pg empirical eq:
# temp <- pg_dat %>% filter(Site=="Church",year.y==2017)
# mean(temp$height.y,na.rm=T)
# rm(temp)
# ht_solve<-function(x){(283.9*x^0.5137)*(1-(1-(1.37/(283.9*x^0.5137))^0.5829)*exp(-x*(14-0.5)))^(1.71556)}-5
# nleqslv(10,ht_solve)
# 0.01073698


mod4 <-predict_pg_reg(pg_pred,"Elastic Net",vars4,lambda1=0.01,n=16)
mod5 <-predict_pg_reg(pg_pred,"EN-Subsets",vars5,lambda1=0.01,n=16)
mod6 <-predict_pg_reg(pg_pred,"EN-1SE",vars6,lambda1=0.01,n=16)
mod7 <-predict_pg_reg(pg_pred,"Stepwise-AIC",vars7,lambda1=0.01,n=16)
mods<-rbind(mod4,mod5,mod6,mod7)
mod_means<- mods %>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))
#mod_means<-rbind(mod_means,gar_pred(14,0.0167122))

a<-mod_means %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()+labs(title="Prince George, SSP585",x="Year",y="Mean Height (m)")


mod4f <-predict_pg_reg(faib_pred,"Elastic Net",vars4,lambda1=0.01,n=17)
mod5f <-predict_pg_reg(faib_pred,"EN-Subsets",vars5,lambda1=0.01,n=17)
mod6f <-predict_pg_reg(faib_pred,"EN-1SE",vars6,lambda1=0.01,n=17)
#mod7f <-predict_pg_rni(faib_pred,"EN-NI",vars6,lambda1=0.23,n=17)
mod7f <-predict_pg_reg(faib_pred,"Step AIC",vars7,lambda1=0.01,n=17)
modsf<-rbind(mod4f,mod5f,mod6f,mod7f)
mod_meansf<- modsf%>% group_by(Model,year.y) %>% summarize(mean_height=mean(height.y))
b<-mod_meansf %>% ggplot(aes(x=year.y,y=mean_height,color=Model))+geom_line()+labs(title="FAIB Plot 59144, SSP585",x="Year",y="Mean Height (m)")

ggarrange(a,b,nrow=1,common.legend=TRUE,legend="bottom")


pdf("fig1_v2.pdf", width = 7, height = 3.5)
ggarrange(a,b,nrow=1,common.legend=TRUE,legend="bottom")
dev.off()


#### Climate figures for defence ####

faib_base<-read.csv("./Data/FAIB/cna_faib_1926-2019SY.csv")
faib_base <- faib_base %>% select(Year,MAT,MWMT,MCMT) %>% mutate(year=Year-1925) %>% group_by(year) %>% summarise(across(everything(), mean)) %>% mutate(Climate="Historic")
#faib_base$year<-faib_base$Year-1925
faib_585 <- read.csv("./Data/FAIB/cna_faib_13GCMs_ensemble_ssp585_2011-2100Y.csv")
faib_370 <- read.csv("./Data/FAIB/cna_faib_13GCMs_ensemble_ssp370_2011-2100Y.csv")
faib_585 <- faib_585 %>% select(Year,MAT,MWMT,MCMT) %>% mutate(year=Year-2000) %>% group_by(year) %>% summarise(across(everything(), mean)) %>% mutate(Climate="")
faib_370 <- faib_370 %>% select(Year,MAT,MWMT,MCMT) %>% mutate(year=Year-2000) %>% group_by(year) %>% summarise(across(everything(), mean)) %>% mutate(Climate="SSP370")
faib_clim_dat <- rbind(faib_base,faib_370,faib_585)

str(faib_clim_dat)
faib_clim_dat %>% ggplot(aes(x=year,y=MAT,color=Climate))+geom_smooth()
faib_base %>% ggplot(aes(x=year,y=MAT))+geom_smooth()

faib_clim_dat %>% ggplot(aes(x=Year,y=MAT,color=Climate,fill=Climate))+geom_line(linewidth=1)+geom_ribbon(aes(ymin=MCMT,ymax=MWMT),alpha=.3, linetype=0)+
  labs(y="Temperature (C)",title="Historic and Projected (Shared Socioeconomic Pathways) Climate in British Columbia")+
  geom_segment(aes(x=1926,xend=2100,y=2.78,yend=2.78),linetype = "dotted",linewidth=0,color="black")
#geom_segment(aes(x = 2, xend = 4, y = 20, yend = 20))
pdf("clim_fig_v2.pdf", width = 9, height = 4)
dev.off()

faib_clim_dat %>% ggplot(aes(x=MAT,fill=Climate,color=Climate))+geom_histogram(bins=35,alpha=0.6,position = "identity")+theme(legend.position = "bottom")+
  scale_color_manual(values=wesanderson::wes_palette("Zissou1")[c(2,4,5)])+scale_fill_manual(values=wesanderson::wes_palette("Zissou1")[c(2,4,5)])+
  labs(title="Empirical Distribution of Mean Annual Temperature",x="Mean Annual Temperature (C)",y="Count")


faib_clim_dat %>% ggplot(aes(x=MAT,fill=Climate,color=Climate))+geom_histogram(bins=35,alpha=0.6,position = "identity",linewidth=1)+theme(legend.position = "bottom")+
  scale_color_manual(values=wesanderson::wes_palette("Moonrise3")[c(1,2,5)])+scale_fill_manual(values=wesanderson::wes_palette("Moonrise3")[c(1,2,5)])+
  labs(title="Empirical Distribution of Mean Annual Temperature",x="Mean Annual Temperature (C)",y="Count")

faib_clim_dat %>% ggplot(aes(x=MAT,fill=Climate,color=Climate))+geom_histogram(bins=35,alpha=0.6,position = "identity",linewidth=1)+theme(legend.position = "bottom")+
  scale_color_manual(values=wesanderson::wes_palette("FantasticFox1")[c(3,2,1)])+scale_fill_manual(values=wesanderson::wes_palette("FantasticFox1")[c(3,2,1)])+
  labs(title="Empirical Distribution of Mean Annual Temperature",x="Mean Annual Temperature (C)",y="Count")

pdf("clim_fig2.pdf",width=7,height=7)
dev.off()

#### Defence Tables ####
c2_res <- read.csv("./Data/ch2_results_1.csv")

b<-c2_res %>% filter(Type=="bias") %>% mutate(Bias=abs(Result),Method=factor(Model,levels=
  c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
  Setting = factor(Setting,levels=c("validation","external"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Setting,y=Bias,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="Bias by Validation Setting",x="Validation Setting",y="Bias (m)")+scale_fill_manual(values=custom_pallete)

a<-c2_res %>% filter(Type=="MSE",Setting %in% c("validation","external")) %>% mutate(Method=factor(Model,levels=
                                                                            c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
                                           Setting = factor(Setting,levels=c("CV","validation","external"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Setting,y=Result,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="MSE by Validation Setting",x="Validation Setting",y=expression(Mean~Squared~Error~(m^2)))+scale_fill_manual(values=custom_pallete)
library(wesanderson)
custom_pallete <- c(wes_palette("Moonrise2")[2], wes_palette("Moonrise3")[c(5,4,3,2,1)])
ggpubr::ggarrange(a,b,nrow=1,common.legend = TRUE,legend="bottom")
grid.arrange(a,b,nrow=1)
pdf("height_res1.pdf",width=8,height=4)
dev.off()

# elevation figure:
c2_res2 <- read.csv("./Data/ch2_results_2.csv")

a<-c2_res2 %>% filter(Type=="MSE",Setting=="CV") %>% mutate(Method=factor(Model,levels=
   c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
   Setting = factor(Setting,levels=c("CV","validation","external"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Elevation,y=Result,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="MSE: Cross-Validation",x="Elevation in Model?",y=expression(Mean~Squared~Error~(m^2)))+scale_fill_manual(values=custom_pallete)

b<-c2_res2 %>% filter(Type=="MSE",Setting=="validation") %>% mutate(Method=factor(Model,levels=
  c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
  Setting = factor(Setting,levels=c("CV","validation","external"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Elevation,y=Result,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="MSE: Validation Set",x="Elevation in Model?",y=expression(Mean~Squared~Error~(m^2)))+scale_fill_manual(values=custom_pallete)

c<-c2_res2 %>% filter(Type=="MSE",Setting=="external") %>% mutate(Method=factor(Model,levels=
                                                                                    c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
                                                                    Setting = factor(Setting,levels=c("CV","validation","external"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Elevation,y=Result,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="MSE: External Data",x="Elevation in Model?",y=expression(Mean~Squared~Error~(m^2)))+scale_fill_manual(values=custom_pallete)

ggpubr::ggarrange(a,b,c,nrow=1,common.legend = TRUE,legend="bottom")
pdf("elev_res1.pdf",width=8,height=4)
dev.off()


b<-c2_res2 %>% filter(Type=="bias",Setting=="validation") %>% mutate(Result=abs(Result),
  Method=factor(Model,levels=
                                                                                    c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
                                                                    Setting = factor(Setting,levels=c("CV","validation","external"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Elevation,y=Result,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="Bias: Validation Set",x="Elevation in Model?",y=expression(Mean~Squared~Error~(m^2)))+scale_fill_manual(values=custom_pallete)

c<-c2_res2 %>% filter(Type=="bias",Setting=="external") %>% mutate(Result=abs(Result),Method=factor(Model,levels=
                                                                                  c("Step-AIC"     ,   "CV-Subsets"   ,   "CV-Subsets-iter", "EN-Subsets", "Elastic-Net"  ,   "EN-1SE")),
                                                                  Setting = factor(Setting,levels=c("CV","validation","external"))) %>%
  mutate(Method=forcats::fct_recode(Method,"CV-Subsets (EN)"="EN-Subsets","CV-Subsets (iter)"="CV-Subsets-iter","Elastic-Net (1-SE)"="EN-1SE"))%>%
  ggplot(aes(x=Elevation,y=Result,fill=Method)) + geom_bar(stat="identity",position="dodge")+
  labs(title="Bias: External Data",x="Elevation in Model?",y=expression(Mean~Squared~Error~(m^2)))+scale_fill_manual(values=custom_pallete)



ggpubr::ggarrange(b,c,nrow=1,common.legend = TRUE,legend="bottom")
pdf("elev_res2.pdf",width=8,height=4)
dev.off()

# C3 figure:
c3_res <- read.csv("./Data/c3_res.csv")
a1<-c3_res %>% filter(Class=="Height") %>% mutate(Validation=factor(Validation,levels=c("Combined","PSP-Only","External")))%>%
  ggplot(aes(x=Validation,y=MSE,fill=Dataset))+geom_bar(stat="identity",position="dodge")+
  labs(title="MSE by Validation Setting",x="Validation Setting")+
  scale_fill_manual(values=wes_palette("FantasticFox1")[c(3,2,4)])
b1<-c3_res %>% filter(Class=="Height") %>% mutate(Validation=factor(Validation,levels=c("Combined","PSP-Only","External")))%>%
  ggplot(aes(x=Validation,y=Bias,fill=Dataset))+geom_bar(stat="identity",position="dodge")+
  labs(title="Bias by Validation Setting",x="Validation Setting")+
  scale_fill_manual(values=wes_palette("FantasticFox1")[c(3,2,4)])
pdf("sense_res1.pdf",width=8,height=4)
ggpubr::ggarrange(a1,b1,nrow=1,common.legend = TRUE,legend="bottom")
dev.off()

a2<-c3_res %>% filter(Class=="Mort") %>% mutate(Validation=factor(Validation,levels=c("Combined","PSP-Only","External")))%>%
  ggplot(aes(x=Validation,y=MSE,fill=Dataset))+geom_bar(stat="identity",position="dodge")+
  labs(title="Area Under the Curve by Validation Setting",x="Validation Setting",y="AUC")+
  scale_fill_manual(values=wes_palette("FantasticFox1")[c(3,2,4)])
b2<-c3_res %>% filter(Class=="Mort") %>% mutate(Validation=factor(Validation,levels=c("Combined","PSP-Only","External")))%>%
  ggplot(aes(x=Validation,y=Bias,fill=Dataset))+geom_bar(stat="identity",position="dodge")+
  labs(title="Phi by Validation Setting",x="Validation Setting",y="phi")
pdf("sense_res2.pdf",width=8,height=4)
ggpubr::ggarrange(a2,b2,nrow=1,common.legend = TRUE,legend="bottom")
dev.off()

pdf("sense_res2_v2.pdf",width=6,height=4)
a2
dev.off()                     


elasticpsp <- readRDS(file="./Data/Modelling Datasets/ht_sens_psp_en.RDS")
plot(elasticpsp)
test_seq<-exp(seq(from=-10,to=-1,by=0.1))


x<-ht_tr_psp[,c(8,10,12,14:215,219:305)]
# this removes all radiation variables
which(!complete.cases(x))
x<-x[ , colSums(is.na(x))==0]
x$logheight<-log(x$height.x)
cv_folds_tr <- cv_climate(x,"mean_Tave_sm")
#x$Type <-ifelse(x$Type=="PSP",1,0)
x<-as.matrix(x)
#x$int<-x$DBH.x*x$BAL
y<- ihs(ht_tr_psp$dh)
#cv.glmnet(x,y,alpha=0,lambda = c(0.01,0.05,0.1),nfolds=4,parallel=TRUE)
elastic_psp<-cv.glmnet(x,y,nfolds=10)
#elastic_psp_2<-cv.glmnet(x,y,foldid= index_to_folds(cv_folds_tr),parallel=TRUE)
plot(elastic_psp)
