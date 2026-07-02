#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 04- Process FIA (USDA FS) data #


######## FIA Read ##########
######### V 2.0 ############
##### Liam Gilson ##########

#### Dependencies ####

library(DBI)
library(RSQLite)
library(ggplot2)
library(dplyr)

# statelist_test <- c("AZ","NM")
# 
# test<-FIA_reprocess(statelist_test)
# test3<-test[[1]][[2]]
# 
# test2<-unlist(test,recursive = F)
# test2 <- rbind(test[[1]][[2]],test[[2]][[2]])
# test3 <- lapply(test,`[[`,2)
# test3<-do.call("rbind",test3)

#### Read FIA databases and compile ####

statelist_fia <- c("AK","AZ","CA","CO","ID","ME","MI","MN","MT","NH","NM","NV","NY","OR","SD","UT","VT","WA","WY")
fia_temp1<-FIA_reprocess(statelist_fia)

fia_tree_list <- lapply(fia_temp1,`[[`,1)
fia_plot_list <- lapply(fia_temp1,`[[`,2)
fia_dat <- do.call("rbind",fia_tree_list)
fia_plot <- do.call("rbind",fia_plot_list)

#saveRDS(fia_dat,"./Data/FIA_DB/fia_datv2.RDS")
# careful not to overwrite fia_plot, as it should be identical to earlier version
#saveRDS(fia_plot,"./Data/FIA_DB/fia_plotv2.RDS")
# elevations would need to be re-done with slow API call in FIA_read.R


#### Stop here if updating FIA tree data #######################################



#### Code for getting elevations, getting plot info for comparisons code ####

missing_list <- which(is.na(fia_plot$ELEV))
# 112 missing elevations
library(httr)
library(jsonlite)
for (i in 1:length(missing_list)){
  index<-missing_list[i]
  lati<-fia_plot$LAT[index]
  loni<-fia_plot$LON[index]
  res <- GET("https://epqs.nationalmap.gov/v1/json?",query=list(x=loni,y=lati,units="Feet",wkid=4326,includeDate="False"))
  res<-fromJSON(rawToChar(res$content))
  fia_plot$ELEV[index]<-ifelse(is.null(res$value),NA,round(as.numeric(res$value),0))
}
which(is.na(fia_plot$ELEV))
# for some reason, obs 3781 doesn't work- manual entry returns 323, though
#fia_plot$ELEV[3781]<-323
#res<-GET("https://epqs.nationalmap.gov/v1/json?",query=list(x=loni,y=lati,units="Meters",wkid=4326,includeDate="False"))

# save again to avoid running the slow API calls again
#saveRDS(fia_plot,"./Data/FIA_DB/fia_plot.RDS")



#### Extract plots for passing to FIA ####

plot_cn <- unique(c(fia_dat$PLT_CN))
state_list<-c("AK","AZ","CA","CO","ID","ME","MI","MN","MT","NH","NM","NV","NY","OR","SD","UT","VT","WA","WY")

#FIA_OR <- dbConnect(RSQLite::SQLite(), "./Data/FIA_OR/SQLite_FIADB_OR.db")

test<-fia_plot_exr(plot_cn,state_list)

setdiff(test$CN,fia_dat$PLT_CN)

write.csv(test,file="./Data/FIA_DB/fia_plots.csv")


# debugging code
fia_plot[20129,]
fia_dat |> filter(ID == "330200700340")
fia_plot |> filter(ID == 330200700340)
min(fia_dat$year.x,na.rm=T)
# 1998 is earliest actual interval measurment

fia_plot[which(is.na(fia_plot$ELEV)),]
