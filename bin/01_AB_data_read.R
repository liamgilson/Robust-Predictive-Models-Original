#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 01- Processing of Alberta Provenance Trial Data #


############# AB Data Clean & Process ##########################################

library(stringr)
library(dplyr)
library(tidyr)

# G103 trial: Data read

# Functions G103
#(Trial) ID site proc Block tree height.x height.y year.x year.y dt dh
# need to filter unplanted/9999 trees out too
# conceptually start by making data wide, then long
# x<-G103B
# Year<-1980
# Age<-4

which2<-function(x,y){
  result<-rep(NA,length(x))
  for (i in 1:length(x)){
    result[i]<-which(y==x[i])
  }
  return(result)
}

proc_G103 <- function(x,Year,Age){
  x<-x[-(if(identical(which(x$ACC==9999), integer(0))) (nrow(x)+1) else which(x$ACC==9999)),]
  years<-str_extract(colnames(x),"(?<=HT)([0-9]{2})")
  years1<-unique(years[!is.na(years)])
  years2<-as.numeric(years1)-Age+Year
  n<-length(years1)-1
  for (i in 1:n){
    varname<-paste0("dh_",years2[i],"_",years2[i+1])
    x[[varname]]<-x[,which(years==years1[i+1])]-x[,which(years==years1[i])]
  }
  x <-pivot_longer(x,cols=starts_with("dh_"),names_to=c("year.x","year.y"),names_prefix = "dh_",names_sep="_",values_to="dh")
  x$ID2<-paste0(x$TRIAL,sprintf("%02d",x$ACC),sprintf("%04d",x$ID))
  x$year.x<-as.numeric(x$year.x)
  x$year.y<-as.numeric(x$year.y)
  x$dt<-x$year.y-x$year.x
  temp<-which2(years1[which2(x$year.x,years2)],years)
  result<-rep(NA,length(temp))
  for(i in 1:length(temp)){
    result[i]<-x[i,temp[i]]}
  x$height.x<-result
  temp<-which2(years1[which2(x$year.y,years2)],years)
  result<-rep(NA,length(temp))
  for(i in 1:length(temp)){
    result[i]<-x[i,temp[i]]}
  x$height.y<-result
#  x$ID <-paste0(x$site,x$Prov,sprintf("%02d",x$Block),sprintf("%02d",as.integer(x$Rep)))
  x <- x %>% select(any_of(c("ID2","TRIAL","ACC","REP","ROW","TREE","height.x","height.y","year.x","year.y","dh","dt")))
  x$height.x<-unlist(x$height.x)
  x$height.y<-unlist(x$height.y)
  x <- x %>% group_by(year.x) %>% mutate(ht_pct.x=ecdf(height.x)(height.x),ht_pct.y=ecdf(height.y)(height.y))
  return(x)
}

# G103B

G103B <- read.csv("./Data/AB_data/G103B.csv", na.strings=".")
G103B_proc<-proc_G103(G103B,Year=1980,Age=4)

G103C <- read.csv("./Data/AB_data/G103C.csv", na.strings=".")
G103C_proc<-proc_G103(G103C,Year=1980,Age=4)

G103D <- read.csv("./Data/AB_data/G103D.csv", na.strings=".")
G103D_proc<-proc_G103(G103D,Year=1981,Age=4)

G103E <- read.csv("./Data/AB_data/G103E.csv", na.strings=".")
G103E_proc<-proc_G103(G103E,Year=1981,Age=4)

G103F <- read.csv("./Data/AB_data/G103F.csv", na.strings=".")
G103F_proc<-proc_G103(G103F,Year=1981,Age=4)

G103G <- read.csv("./Data/AB_data/G103G.csv", na.strings=".")
G103G_proc<-proc_G103(G103G,Year=1982,Age=4)

G103H <- read.csv("./Data/AB_data/G103H.csv", na.strings=".")
G103H_proc<-proc_G103(G103H,Year=1982,Age=4)

G103RW <- read.csv("./Data/AB_data/G103RW.csv", na.strings=".")
G103RW_proc<-proc_G103(G103RW,Year=1982,Age=4)
# ACC is provenance code

G103<-rbind(G103B_proc,G103C_proc,G103D_proc,G103E_proc,G103F_proc,G103G_proc,G103H_proc,G103RW_proc)
G103$exp<-"G103"
saveRDS(G103,"./Data/AB_data/G103.RDS")

############# G276 ############
G276A <- read.csv("./Data/AB_data/G276A.csv", na.strings=".")
G276A_proc<-proc_G103(G276A,Year=1993,Age=3)
G276B <- read.csv("./Data/AB_data/G276B.csv", na.strings=".")
G276B_proc<-proc_G103(G276B,Year=1993,Age=3)
G276C <- read.csv("./Data/AB_data/G276C.csv", na.strings=".")
G276C_proc<-proc_G103(G276C,Year=1993,Age=3)
G276D <- read.csv("./Data/AB_data/G276D.csv", na.strings=".")
G276D_proc<-proc_G103(G276D,Year=1993,Age=3)
G276<-rbind(G276A_proc,G276B_proc,G276C_proc,G276D_proc)
G276$exp<-"G276"
saveRDS(G276,"./Data/AB_data/G276.RDS")

############ G277 #############

############# G276 ############
G277A <- read.csv("./Data/AB_data/G277A.csv", na.strings=".")
G277A_proc<-proc_G103(G277A,Year=1994,Age=3)
G277B <- read.csv("./Data/AB_data/G277B.csv", na.strings=".")
# Note: no REP in G277B
G277B$REP<-1
G277B_proc<-proc_G103(G277B,Year=1993,Age=3)
G277<-rbind(G277A_proc,G277B_proc)
G277$exp<-"G277"
saveRDS(G277,"./Data/AB_data/G277.RDS")

############# G347 ############

G347A <- read.csv("./Data/AB_data/G347A.csv", na.strings=".")
G347A_proc<-proc_G103(G347A,Year=1994,Age=3)
G347B <- read.csv("./Data/AB_data/G347B.csv", na.strings=".")
G347B_proc<-proc_G103(G347B,Year=1993,Age=3)
G347<-rbind(G347A_proc,G347B_proc)
G347$exp<-"G347"
saveRDS(G347,"./Data/AB_data/G347.RDS")

AB_data<-rbind(G103,G276,G277,G347)
#saveRDS(AB_data,"./Data/AB_data/AB_datav2.RDS")
AB_data<-readRDS("./Data/AB_data/AB_datav2.RDS")


sum(!is.na(AB_data$dh))




# x<-G277B
# Year<-1993
# Age<-3
# Debugging code

# years<-str_extract(colnames(G103B),"(?<=HT)([0-9]{2})")
# years1<-unique(years[!is.na(years)])
# #years1<-as.numeric(years1)-4+1980
# years2<-as.numeric(years1)-4+1980
# i<-1
# temp_x<-G103B
# i<-4
# varname<-paste0("dh_",years2[i],"_",years2[i+1])
# temp_x[[varname]]<-x[,which(years==years1[i+1])]-x[,which(years==years1[i])]
# temp_x <-pivot_longer(temp_x,cols=starts_with("dh_"),names_to=c("year.x","year.y"),names_prefix = "dh_",names_sep="_",values_to="dh")
# temp_x <- temp_x %>% select(ID,Trial,REP,ACC,ROW,TREE,height.x,height.y,year.x,year.y,dh)
# 
# G103B[[varname]]<-NA
# 
# temp<-which2(years1[which2(x$year.x,years2)],years)
# for(i in 1:length(temp)){
#   result[i]<-x[i,temp[i]]}
# which2(years1[which2(x$year.x,years2)],years)
