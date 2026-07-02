#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 02- Processing of Ontario Provenance Trial Data #


# 410 data read and cleanup
library(dplyr)
library(tidyr)
library(stringr)

# functions to speed things up
proc410<-function(x,Site,Year){
  x <- x %>% rename(height=starts_with("X19"))
  x$site <- Site
  x$year <- Year
  is.na(x)<-x == "."
  x$height<-as.numeric(x$height)/100
  x$ht_pct <- ecdf(x$height)(x$height)
  return(x)
}

ht_inc<-function(df1,df2){
  dup_rows<- which(duplicated(df1[,c("Block","SeedSource","Tree")]))
  if (!identical(dup_rows,integer(0))){
    df1<-df1[-dup_rows,]
    df2<-df2[-dup_rows,]}
  test<-left_join(df1,df2,by=c("Block","SeedSource","Tree","site"))
  test$dh<-test$height.y-test$height.x
  test$dt<-test$year.y-test$year.x
  colnames(test)[colnames(test)=="SeedSource"]<-"Prov"
  test$ID<-paste0(test$site,test$Prov,sprintf("%02d",test$Block),sprintf("%02d",test$Tree))
  test <- test %>% select(ID,site,Prov,Block,Tree,height.x,height.y,year.x,year.y,ht_pct.x,ht_pct.y,dt,dh)
  return(test)
}

proc_410_2000 <- function(x,Site){
  years<-str_extract(colnames(x),"(?<=X)([0-9]{4})(?=_Ht)")
  years1<-unique(years[!is.na(years)])
  stopifnot(length(years1)==2)
  y1<-min(years1)
  y2<-max(years1)
  x$year.x<-as.numeric(y1)
  x$year.y<-as.numeric(y2)
  colnames(x)[which(years==y1)]<-"height.x"
  colnames(x)[which(years==y2)]<-"height.y"
  x$site <- Site
  is.na(x)<-x == "."
  is.na(x)<-x == ""
  x$height.x<-as.numeric(x$height.x)
  x$height.y<-as.numeric(x$height.y)
  x$ht_pct.x <- ecdf(x$height.x)(x$height.x)
  x$ht_pct.y <- ecdf(x$height.y)(x$height.y)
  x$dh<-x$height.y-x$height.x
  x$dt<-x$year.y-x$year.x
  x$ID <-paste0(x$site,x$Prov,sprintf("%02d",x$Block),sprintf("%02d",as.integer(x$Rep)))
  colnames(x)[colnames(x)=="Rep"]<-"Tree"
  x <- x %>% select(ID,site,Prov,Block,Tree,height.x,height.y,year.x,year.y,ht_pct.x,ht_pct.y,dt,dh)
  return(x)
}

ht_inc_link<-function(df1,df2){
  df1<-df1 %>% select(ID,site,Prov,Block,Tree,height.y,year.y,ht_pct.y) %>% rename(height.x=height.y,year.x=year.y,ht_pct.x=ht_pct.y)
  df2<-df2 %>% select(ID,site,Prov,Block,Tree,height.x,year.x,ht_pct.x) %>% rename(height.y=height.x,year.y=year.x,ht_pct.y=ht_pct.x)
  df3<-inner_join(df1,df2,by=c("ID","site","Prov","Block","Tree"))
  df3$dt<-df3$year.y-df3$year.x
  df3$dh<-df3$height.y-df3$height.x
  return(df3)
}

ht_inc_link_2<-function(df1,df2){
  df1$ID<-paste0(df1$site,df1$SeedSource,sprintf("%02d",df1$Block),sprintf("%02d",df1$Tree))
  colnames(df1)[colnames(df1)=="SeedSource"]<-"Prov"
  df1<-df1 %>% select(ID,site,Prov,Block,Tree,height,ht_pct,year) %>% rename(height.x=height,ht_pct.x=ht_pct,year.x=year)
  df2<-df2 %>% select(ID,site,Prov,Block,Tree,height.x,ht_pct.x,year.x) %>% rename(height.y=height.x,ht_pct.y=ht_pct.x,year.y=year.x)
  df3<-inner_join(df1,df2,by=c("ID","site","Prov","Block","Tree"))
  df3$dt<-df3$year.y-df3$year.x
  df3$dh<-df3$height.y-df3$height.x
  return(df3)
}

process_410<-function(link1,link2,link3,link4,year=c(1980,1984,1988),site){
  first<-read.csv(link1)
  second<-read.csv(link2)
  third<-read.csv(link3)
  fourth<-read.csv(link4)
  first<-proc410(first,Site=site,Year=year[1])
  second<-proc410(second,Site=site,Year=year[2])
  third<-proc410(third,Site=site,Year=year[3])
  fourth<-proc_410_2000(fourth,Site=site)
  temp1<-ht_inc(first,second)
  temp2<-ht_inc(second,third)
  temp3<-ht_inc_link(temp2,fourth)
  temp4<-rbind(temp1,temp2,temp3,fourth)
  return(temp4)
}

process_410_2<-function(link1,link2,link4,year=c(1980,1984),site){
  first<-read.csv(link1)
  second<-read.csv(link2)
  fourth<-read.csv(link4)
  first<-proc410(first,Site=site,Year=year[1])
  second<-proc410(second,Site=site,Year=year[2])
  fourth<-proc_410_2000(fourth,Site=site)
  temp1<-ht_inc(first,second)
  temp3<-ht_inc_link(temp1,fourth)
  temp4<-rbind(temp1,temp3,fourth)
  return(temp4)
}


################################################################################
############################ Read data #########################################
################################################################################

# A1

A1_410<-process_410("./Data/410_csv/410A1_1980.csv","./Data/410_csv/410A1_1983.csv",
                         "./Data/410_csv/410A1_1988.csv","./Data/410_csv/410A1_2015.csv",
                         year=c(1980,1983,1988),site="A1")

# A2

A2_410<-process_410("./Data/410_csv/410A2_1978.csv","./Data/410_csv/410A2_1982.csv",
                         "./Data/410_csv/410A2_1988.csv","./Data/410_csv/410A2_2015.csv",
                         year=c(1978,1982,1988),site="A2")

# A3

A3_410<-process_410("./Data/410_csv/410A3_1978.csv","./Data/410_csv/410A3_1982.csv",
                         "./Data/410_csv/410A3_1988.csv","./Data/410_csv/410A3_2014.csv",
                         year=c(1978,1982,1988),site="A3")

# A4: manual with only two files

A4_1983 <- read.csv("./Data/410_csv/410A4_1983.csv")
A4_1983 <- proc410(A4_1983,Site="A4",Year=1983)
A4_2012 <- read.csv("./Data/410_csv/410A4_2012.csv")
A4_2012 <- proc_410_2000(A4_2012,Site="A4")
# special function for this case:
A4_temp1 <- ht_inc_link_2(A4_1983,A4_2012)
A4_410 <- rbind(A4_temp1,A4_2012)

# B1

B1_410<-process_410("./Data/410_csv/410B1_1979.csv","./Data/410_csv/410B1_1983.csv",
                    "./Data/410_csv/410B1_1988.csv","./Data/410_csv/410B1_2014.csv",
                    year=c(1979,1983,1988),site="B1")
# NAs coerced by trailing or leading spaces on periods, this is fine

# B2

# code for the situation of just three files, 2 in 80s, 1 2000s
B2_410<-process_410_2("./Data/410_csv/410B2_1983.csv","./Data/410_csv/410B1_1988.csv",
                      "./Data/410_csv/410B2_2014.csv",
                      year=c(1983,1988),site="B2")



# B3

B3_410<-process_410_2("./Data/410_csv/410B3_1984.csv","./Data/410_csv/410B3_1988.csv",
                      "./Data/410_csv/410B3_2014.csv",
                      year=c(1984,1988),site="B3")
# problem with Rep/Tree being classed as 'num' not 'int' for regex:
# x<-read.csv("./Data/410_csv/410B3_2014.csv")
# str(x)
# as.integer(x$Rep)
# paste0(x$site,x$Prov,sprintf("%02d",x$Block),sprintf("%02d",x$Rep))

# D1

D1_410<-process_410_2("./Data/410_csv/410D1_1981.csv","./Data/410_csv/410D1_1985.csv",
                      "./Data/410_csv/410D1.csv",
                      year=c(1981,1985),site="D1")

# D2

D2_410<-process_410_2("./Data/410_csv/410D2_1982.csv","./Data/410_csv/410D2_1985.csv",
                      "./Data/410_csv/410D2.csv",
                      year=c(1982,1985),site="D2")


# D3

# This is the one that was measured in both 2014 and 2015
D3_temp1<-read.csv("./Data/410_csv/410D3_1982.csv")
D3_temp2<-read.csv("./Data/410_csv/410D3_1985.csv")
D3_temp3<-read.csv("./Data/410_csv/410D3_2014.csv")
D3_temp1<-proc410(D3_temp1,Site="D3",Year=1982)
D3_temp2<-proc410(D3_temp2,Site="D3",Year=1985)
D3_temp4<-ht_inc(D3_temp1,D3_temp2)
# 2013 block 1/5, 2013 2,3,4
D3_2013 <- D3_temp3 %>% filter(Block %in% c(1,5))
D3_2013 <- D3_2013[,-c(10:11)]
D3_2013 <- proc_410_2000(D3_2013,Site="D3")
D3_2014 <- D3_temp3 %>% filter(Block %in% c(2:4))
D3_2014 <- D3_2014[,-c(8:9)]
D3_2014 <- proc_410_2000(D3_2014,Site="D3")
# try to merge both to get intervals from 1985
D3_2013_inc <- ht_inc_link(D3_temp4,D3_2013)
D3_2014_inc <- ht_inc_link(D3_temp4,D3_2014)

D3_410<-rbind(D3_temp4,D3_2013_inc,D3_2014_inc,D3_2013,D3_2014)

# E1
# same as A4 situation
E1_1984 <- read.csv("./Data/410_csv/410E1_1984.csv")
E1_1984 <- proc410(E1_1984,Site="E1",Year=1984)
E1_2012 <- read.csv("./Data/410_csv/410E1_2012.csv")
E1_2012 <- proc_410_2000(E1_2012,Site="E1")
# special function for this case:
E1_temp1 <- ht_inc_link_2(E1_1984,E1_2012)
which(duplicated(E1_1984[,1:3]))
E1_1984[2972,]
# 8028 duplicated in 4, 8039 dulplicated in 4, 8026 duplicated in 5
# with no unique code, there is no choice but to removed affected trees
# we cannot link them over time, other option is to remove whole site
E1_temp1 <- E1_temp1 %>% filter(!str_detect(ID,"E1803904|E1802804|E1802605"))
E1_410 <- rbind(E1_temp1,E1_2012)

# E2
# same as E1
E2_1987 <- read.csv("./Data/410_csv/410E2_1987.csv")
E2_1987 <- proc410(E2_1987,Site="E2",Year=1987)
# error from ". " character, still should be NA so coercion is fine here
E2_2012 <- read.csv("./Data/410_csv/410E2.csv")
E2_2012 <- proc_410_2000(E2_2012,Site="E2")
# 2001 and 2013 here

unique(E2_1987$SeedSource)
# there's non-numeric provenances in the 1987 data, we're just going to remove these
# they aren't in the provenance key, who knows what they are/were
# actually maybe coerce to NA to avoid competition calculation issues
E2_1987$SeedSource <- as.integer(E2_1987$SeedSource)
# special function for this case:
E2_temp1 <- ht_inc_link_2(E2_1987,E2_2012)
E2_410 <- rbind(E2_temp1,E2_2012)

# F1

# Note, "5yr" F1 measurements were in April 1987, so will be coded 1986
F1_410<-process_410_2("./Data/410_csv/410F1_1983.csv","./Data/410_csv/410F1_1987.csv",
                      "./Data/410_csv/410F1_2014.csv",
                      year=c(1983,1986),site="F1")

# F2

# messy here- Cornwall was planted in 1983, measured (1 yr) in April 84. No date on "5 yr" heights- probably April 1987, same as F1
F2_410<-process_410_2("./Data/410_csv/410F2_1983.csv","./Data/410_csv/410F2_1984.csv",
                      "./Data/410_csv/410F2.csv",
                      year=c(1983,1986),site="F2")

Data_410<-rbind(A1_410,A2_410,A3_410,A4_410,B1_410,B2_410,B3_410,D1_410,D2_410,D3_410,E1_410,E2_410,F1_410,F2_410)

################################################################################
############################Save and Read ######################################

#saveRDS(Data_410,file="./Data/410_csv/Data_410v2.RDS")
Data_410<-readRDS(file="./Data/410_csv/Data_410v2.RDS")

temp<-unique(Data_410[,c("site","year.x","year.y")])
write.csv(temp,file="./Data/410_csv/site_years.csv")


################################################################################
############################ Climate ###########################################

clim_na_410<-read.csv(file="./Data/410_csv/Site_ClimNA_Input_1978-2016SY.csv")






################################################################################
############################## Elevations/Provenances ##########################

provs_410 <- read.csv(file="./Data/410_csv/prov.csv")
provs_410$elev<-NA
# all longitudes are missing signs
provs_410$long<- -provs_410$long

library(httr)
library(jsonlite)
for (i in 1:nrow(provs_410)){
  lati<-provs_410[i,3]
  loni<-provs_410[i,2]
  res <- GET("http://geogratis.gc.ca/services/elevation/cdem/altitude?",query=list(lat=lati,lon=loni))
  res<-fromJSON(rawToChar(res$content))
  provs_410$elev[i]<-ifelse(is.null(res$altitude),NA,res$altitude)
}
which(is.na(provs_410$elev))
# from EPQS (provenances must be from US sites)
provs_410[12,4]<-346
provs_410[45,4]<-350
provs_410[56,4]<-448
provs_410[66,4]<-535
provs_410[78,4]<-420

# climateNA needs latitude and longitude in that order, and it screws up if they're swapped:
provs_410 <- provs410[,c(1,3,2,4)]
# needs "id1/id2" rename in excel, but otherwise should be fine for climateNA
write.csv(provs_410,file="./Data/410_csv/prov_410_cna.csv")

# res <- GET("http://geogratis.gc.ca/services/elevation/cdem/altitude?",query=list(lat=45.17,lon=-96.6))
# res<-fromJSON(rawToChar(res$content))
# lati<-provs_410[1,3]
# loni<-provs_410[1,2]
# res <- GET("http://geogratis.gc.ca/services/elevation/cdem/altitude?",query=list(lat=lati,lon=loni))
################################################################################
############################## Debugging #######################################

# years<-str_extract(colnames(A1_2015),"(?<=X)([0-9]{4})(?=_Ht)")
# years1<-unique(years[!is.na(years)])
# stopifnot(length(years1)==2)
# y1<-min(years1)
# y2<-max(years1)

# testing code

############### 410A1 ##########################################################
# A1_1980 <- read.csv("./Data/410_csv/410A1_1980.csv")
# #A1_1980 %>% rename(height=starts_with("X19"))
# A1_1980<-proc410(A1_1980,Site="A1",Year=1980)
# 
# A1_1983 <- read.csv("./Data/410_csv/410A1_1983.csv")
# A1_1983<-proc410(A1_1983,Site="A1",Year=1983)
# 
# test1<-ht_inc(A1_1980,A1_1983)
# 
# A1_1988 <- read.csv("./Data/410_csv/410A1_1988.csv")
# A1_1988<-proc410(A1_1988,Site="A1",Year=1988)
# A1_1988_undupe<-A1_1988
# A1_1988_undupe$SN<-NULL
# test2<-ht_inc(A1_1983,A1_1988)
# test2$SN<-NULL
# test3<-rbind(test1,test2)
# #test3$ID<-paste0(test3$site,test3$SeedSource,sprintf("%02d",test3$Block),sprintf("%02d",test3$Tree))
# 
# A1_2015 <- read.csv("./Data/410_csv/410A1_2015.csv")
# unique(A1_2015$Condition_2001)
# test4<-proc_410_2000(A1_2015,Site="A1")
# 
# which(duplicated(test3$ID))
# test5<-ht_inc_link(test2,test4)
# 
# A1_410<-rbind(test3,test5,test4)


# test<-full_join(A1_1980,A1_1983,by=c("Block","SeedSource","Tree"))
# A1_1980 %>% filter(Block==2,SeedSource==8058)
# which(duplicated(A1_1980[,1:3]))
# A1_1980[which(duplicated(A1_1980[,1:3])),]
# 
# # Duplicated data needs to be removed
# dup_rows<- which(duplicated(A1_1980[,1:3]))
# A1_1980 <- A1_1980[-dup_rows,]
# A1_1983 <- A1_1983[-dup_rows,]
# 
# test<-full_join(A1_1980,A1_1983,by=c("Block","SeedSource","Tree","site"))
# test$dh<-test$height.y-test$height.x
# test$dt<-test$year.y-test$year.x
# 
# hist(test$dh)
# which(duplicated(A1_1980[,c("Block","SeedSource","Tree")]))

##############A2
# A2_1978 <- read.csv("./Data/410_csv/410A2_1978.csv")
# #A1_1980 %>% rename(height=starts_with("X19"))
# A2_1978<-proc410(A2_1978,Site="A2",Year=1978)
# 
# A2_1982 <- read.csv("./Data/410_csv/410A2_1982.csv")
# colnames(A2_1982)
# A2_1982<-proc410(A2_1982,Site="A2",Year=1982)
# 
# A2_1<-ht_inc(A2_1978,A2_1982)
# 
# A2_1988 <- read.csv("./Data/410_csv/410A2_1988.csv")
# A2_1988<-proc410(A2_1988,Site="A2",Year=1988)
# 
# A2_2 <-ht_inc(A2_1982,A2_1988)
# 
# A2_2015 <- read.csv("./Data/410_csv/410A2_2015.csv")
# colnames(A2_2015)
# A2_3<-proc_410_2000(A2_2015,Site="A2")
# A2_4<-ht_inc_link(A2_2,A2_3)
# A2_410<-rbind(A2_1,A2_2,A2_4,A2_3)
# A2_410_test<-process_410("./Data/410_csv/410A2_1978.csv","./Data/410_csv/410A2_1982.csv",
#                          "./Data/410_csv/410A2_1988.csv","./Data/410_csv/410A2_2015.csv",
#                          year=c(1978,1982,1988),site="A2")

#### Debug after adding percentiles ####

# issue was that link function wasn't selecting/renaming new ht_pct vars
# link1<-"./Data/410_csv/410A1_1980.csv"
# link2<-"./Data/410_csv/410A1_1983.csv"
# link3<-"./Data/410_csv/410A1_1988.csv"
# link4<-"./Data/410_csv/410A1_2015.csv"
# year<-c(1980,1983,1988)
# site<-"A1"

# This problem was due to "integer(0)" returned from the duplicate test
# if we subset [-interger(0)] it removes the whole DF, so there needs to be a test before the subset command
# ht_inc was changed to address this, and tested on other sites

# df1<-D3_temp1
# df2<-D3_temp2
# 
# dup_rows<- which(duplicated(df1[,c("Block","SeedSource","Tree")]))
# if (!identical(dup_rows,integer(0))){
# df1<-df1[-dup_rows,]
# df2<-df2[-dup_rows,]}
# test<-left_join(df1,df2,by=c("Block","SeedSource","Tree","site"))
# test$dh<-test$height.y-test$height.x
# test$dt<-test$year.y-test$year.x
# colnames(test)[colnames(test)=="SeedSource"]<-"Prov"
# test$ID<-paste0(test$site,test$Prov,sprintf("%02d",test$Block),sprintf("%02d",test$Tree))
# test <- test %>% select(ID,site,Prov,Block,Tree,height.x,height.y,year.x,year.y,ht_pct.x,ht_pct.y,dt,dh)
