#### Code for review of Gilson, Eskelson, Sattler, O'Neill ####
#### Copyright Liam Gilson, 2026 ####
# 00- Functions #
## Functions used throughout this codebase ##

# IHS functions

ihs <- function(x){log(x + sqrt(x^2 + 1))}
ihs_back <- function(x){exp((summary(x)$sigma^2)/2)}

# CNA calculate VPD function

vpd_cna<-function(df1){
  df1$VPD_max_sm <- (610.78*exp(17.2694*df1$Tmax_sm/(df1$Tmax_sm+237.3)))*(1-(df1$RH_sm/100))
  df1$VPD_av_sm <- (610.78*exp(17.2694*df1$Tave_sm/(df1$Tave_sm+237.3)))*(1-(df1$RH_sm/100))
  return(df1)
}

# failsafe version of ecdf for "all NA" scenario
na_ecdf <- function(x){if(all(is.na(x))){return(NA)}else{return(ecdf(x)(x))}}

# CNA data link function

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
# the output of this function needs to be joined to the list of increments by key, year.x, year.y


intervals_cna_fia<-function(df1,df2,key){
  ints<-unique(df1[,c(key,"year.x","year.y")])
  ints2<-ints[complete.cases(ints),]
  temp_list <- list()
  is.na(df2) <- df2 == -9999
  for (i in 1:nrow(ints2)){
    int_years<-ints2[i,]$year.x:ints2[i,]$year.y
    int_site <- ints2[i,key]
    mean_clim<-as.data.frame.list(colMeans(df2 %>% filter(id2==as.numeric(int_site),Year %in% int_years) %>% select(!c(1:6))))
    mean_clim <- mean_clim %>% rename_with( ~ paste("mean", .x, sep = "_"))
    min_clim <- as.data.frame.list(sapply(df2 %>% filter(id2==as.numeric(int_site),Year %in% int_years) %>% select(ends_with("_wt")&!starts_with(c("DD","Eref","CMD","PAS"))),min)) %>%
      rename_with( ~ paste("min", .x, sep = "_"))
    list_max <- c("MAT","CMI","DD1040")
    max_clim <- as.data.frame.list(sapply(df2 %>% filter(id2==as.numeric(int_site),Year %in% int_years) %>% select(ends_with("_sm")|all_of(list_max)),max)) %>%
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

# previous version had error in ID line, state listed twice, no unit
fia_process<-function(x){
  plot1 <- dbGetQuery(x, "SELECT CN,PLOT,PREV_PLT_CN,LAT,LON,ELEV,MEASYEAR,MEASMON,STATECD,UNITCD,COUNTYCD FROM PLOT;")
  tree1 <- dbGetQuery(x, "SELECT PLOT,PLT_CN,INVYR,SUBP,TREE,DIA,HT,HTCD,ACTUALHT,SPCD,STATUSCD, PREV_TRE_CN,CN,DAMTYP1 FROM TREE")
  tree2 <- tree1 %>% group_by(PLOT) %>% filter(any(SPCD %in% c(93,94))) %>% ungroup
  dbDisconnect(x)
  plot1$ID<-paste0(sprintf("%02d",plot1$STATECD),sprintf("%02d",plot1$UNITCD),sprintf("%03d",plot1$COUNTYCD),sprintf("%05d",plot1$PLOT))
  temp_yr <- plot1 %>% select(ID,CN,MEASYEAR)
  tree3<-merge(tree2,temp_yr,by.x="PLT_CN",by.y="CN")
  tree4 <- tree3 %>% filter(SPCD %in% c(93,94),STATUSCD == 1)
  tree4 <- tree4 %>% filter(PREV_TRE_CN != "NA")
  tree4$height.y <- tree4$HT
  tree4$height.x <- NA
  tree4$year.x<-NA
  tree4$year.y<-tree4$MEASYEAR
  for (i in 1:nrow(tree4)){
    key<-tree4$PREV_TRE_CN[i]
    index<-ifelse(identical(integer(0),x<-which(tree3$CN==key)),NA,x)
    tree4$height.x[i]<-ifelse(!is.na(index),tree3$HT[index],NA)
    tree4$year.x[i]<-ifelse(!is.na(index),tree3$MEASYEAR[index],NA)
  }
  tree4$dt<-tree4$year.y-tree4$year.x           
  tree4$dh<-tree4$height.y-tree4$height.x
  plotlist<-unique(tree4$ID)
  cna_vars <- plot1 %>% select(ID,LAT,LON,ELEV,MEASYEAR) %>% filter(ID %in% plotlist) 
  return(list(tree4,cna_vars))
}

fia_plot_exr<-function(plotlist,statelist){
  n <- length(statelist)
  res <- list()
  for (i in 1:n){
    state<-statelist[i]
    x <- dbConnect(RSQLite::SQLite(), paste0("./ Data/FIA_DB/SQLite_FIADB_",state,".db"))
    plot1 <- dbGetQuery(x, "SELECT CN,PLOT,PREV_PLT_CN,LAT,LON,ELEV,MEASYEAR,MEASMON,STATECD,UNITCD,COUNTYCD,INVYR,DESIGNCD FROM PLOT;")
    dbDisconnect(x)
    res[[i]] <- plot1 %>% select(CN,PLOT,PREV_PLT_CN,LAT,LON,ELEV,MEASYEAR,MEASMON,STATECD,UNITCD,COUNTYCD,INVYR,DESIGNCD) %>% filter(CN %in% plotlist)
  }
  return(do.call("rbind",res))
}

fia_process_v2<-function(x){
  plot1 <- dbGetQuery(x, "SELECT CN,PLOT,PREV_PLT_CN,LAT,LON,ELEV,MEASYEAR,MEASMON,STATECD,UNITCD,COUNTYCD FROM PLOT;")
  tree1 <- dbGetQuery(x, "SELECT PLOT,PLT_CN,INVYR,SUBP,TREE,DIA,HT,HTCD,ACTUALHT,SPCD,STATUSCD, PREV_TRE_CN,CN,DAMTYP1 FROM TREE")
  tree2 <- tree1 %>% group_by(PLT_CN) %>% filter(any(SPCD %in% c(93,94)),!all(is.na(HT))) %>% ungroup
  dbDisconnect(x)
  tree2 <- tree2 |> group_by(PLT_CN) |> mutate(ht_pct = ecdf(HT)(HT)) |> arrange(desc(DIA)) |> mutate(BAL = cumsum(pi*(DIA/(2*12))^2))
  plot1$ID<-paste0(sprintf("%02d",plot1$STATECD),sprintf("%02d",plot1$UNITCD),sprintf("%03d",plot1$COUNTYCD),sprintf("%05d",plot1$PLOT))
  temp_yr <- plot1 %>% select(ID,CN,MEASYEAR)
  tree3<-merge(tree2,temp_yr,by.x="PLT_CN",by.y="CN")
  # metric conversion here, only to tree3- HT, DIA, BAL
  tree3$HT<-0.3048*tree3$HT
  tree3$DIA <- 2.54*tree3$DIA
  tree3$BAL <- 0.092903*tree3$BAL
  #tree4 <- tree3 %>% filter(SPCD %in% c(93,94),STATUSCD == 1)
  tree4 <- tree3 %>% filter(SPCD %in% c(93,94),PREV_TRE_CN != "NA")
  # make variables for increment transformation:
  tree4$height.y <- tree4$HT
  tree4$height.x <- NA
  tree4$year.x<-NA
  tree4$year.y<-tree4$MEASYEAR
  tree4$dbh.x<-NA
  tree4$dbh.y<-tree4$DIA
  tree4$status.x<-NA
  tree4$status.y <- tree4$STATUSCD
  tree4$ht_pct.x <-NA
  tree4$ht_pct.y<-tree4$ht_pct
  tree4$BAL.x<-NA
  tree4$BAL.y<-tree4$BAL
  for (i in 1:nrow(tree4)){
    key<-tree4$PREV_TRE_CN[i]
    index<-ifelse(identical(integer(0),x<-which(tree3$CN==key)),NA,x)
    tree4$height.x[i]<-ifelse(!is.na(index),tree3$HT[index],NA)
    tree4$year.x[i]<-ifelse(!is.na(index),tree3$MEASYEAR[index],NA)
    tree4$dbh.x[i]<-ifelse(!is.na(index),tree3$DIA[index],NA)
    tree4$ht_pct.x[i]<-ifelse(!is.na(index),tree3$ht_pct[index],NA)
    tree4$BAL.x[i]<-ifelse(!is.na(index),tree3$BAL[index],NA)
    tree4$status.x[i]<-ifelse(!is.na(index),tree3$STATUSCD[index],NA)
  }
  tree4$dt<-tree4$year.y-tree4$year.x           
  tree4$dh<-tree4$height.y-tree4$height.x
  plotlist<-unique(tree4$ID)
  cna_vars <- plot1 %>% select(ID,LAT,LON,ELEV,MEASYEAR) %>% filter(ID %in% plotlist) 
  return(list(tree4,cna_vars))
}

FIA_reprocess <- function(statelist){
  n <- length(statelist)
  res <- list()
  for (i in 1:n){
    state<-statelist[i]
    x <- dbConnect(RSQLite::SQLite(), paste0("./Data/FIA_DB/SQLite_FIADB_",state,".db"))
    temp<-fia_process_v2(x)
    res[[i]] <- temp
  }
  return(res)
}

#### CV Functions ####

#### Cross Validation Functions ####
cv_function<-function(data,index,formula){
  n <- length(index)
  results<-rep(NA,n)
  for (i in 1:n){
    val_temp <- data[index[[i]],]
    train_temp <- data[-index[[i]],]
    mod_temp <- lm(formula,data=train_temp)
    results[i]<-sum((ihs(val_temp$dh) -predict(mod_temp,val_temp))^2,na.rm=T)/nrow(val_temp)
  }
  sum(results)
}

hv_function<-function(data,index,formula){
  n <- length(index)
  results<-rep(NA,n)
  test<-function(x,n){ifelse(x > 1,ifelse(x < n,return(c(x-1,x,x+1)),return(c(x-1,x))),return(c(x,x+1)))}
  for (i in 1:n){
    val_temp <- data[index[[i]],]
    train_temp <- data[-unlist(index[test(i,n)]),]
    mod_temp <- lm(formula,data=train_temp)
    results[i]<-sum((ihs(val_temp$dh) -predict(mod_temp,val_temp))^2,na.rm=T)/nrow(val_temp)
  }
  sum(results)
}

# for regularization (glmnet)

# data<-ht_dat_na
# index<-index_cv1
# variables<-vars1
# lambda<-0.1
# i<-1
# # nfit = as.matrix(cbind2(1, newx) %*% nbeta)
# newx<-xval
# object <- mod_temp
# a0 = t(as.matrix(object$a0))
# rownames(a0) = "(Intercept)"
# nbeta = methods::rbind2(a0, object$beta)
# cbind2(1, newx) %*% nbeta
# x<-train_temp[,variables]



cv_function_reg<-function(data,index,variables,lambda){
  n<-length(index)
  results<-rep(NA,n)
  data$Type <-ifelse(data$Type=="PSP",1,0)
  #data$Type<-as.factor(data$Type)
  for (i in 1:n){
    train_temp <- data[-index[[i]],]
    x<-train_temp[,variables]
    x$logheight<-log(x$height.x)
    #x$int<-x$DBH.x*x$BAL
    y<- ihs(train_temp$dh)
    val_temp <- data[index[[i]],]
    xval<- val_temp[,variables]
    xval$logheight<-log(xval$height.x)
    #xval$int<-xval$DBH.x*xval$BAL
    xval<-as.matrix(xval)
    mod_temp <- glmnet(x,y,alpha=0,lambda = lambda)
    results[i]<-sum((ihs(val_temp$dh) -predict(mod_temp,xval))^2,na.rm=T)/nrow(val_temp)
  }
  sum(results)
}

cv_reg_search <- function(data,variables,folds,n=10){
  train_temp <- data
  x<-train_temp[,variables]
  x$logheight<-log(x$height.x)
  y<- ihs(train_temp$dh)
  x$Type <-ifelse(x$Type=="PSP",1,0)
  x<-as.matrix(x)
  #cv.glmnet(x,y,alpha=0,nfolds=4,parallel=TRUE)
  #model<-ifelse(is.na(folds),cv.glmnet(x,y,nfolds=n),cv.glmnet(x,y,nfolds=n,foldid=folds))[[1]]
  #model <- cv.glmnet(x,y,nfolds=n,alpha=0)
  model <- cv.glmnet(x,y,nfolds=n,foldid=folds,alpha=0)
  return(model)
}

cv_reg_search2 <- function(data,variables,folds,n=10,lam){
  train_temp <- data
  x<-train_temp[,variables]
  x$logheight<-log(x$height.x)
  y<- ihs(train_temp$dh)
  x$Type <-ifelse(x$Type=="PSP",1,0)
  x<-as.matrix(x)
  #cv.glmnet(x,y,alpha=0,nfolds=4,parallel=TRUE)
  #model<-ifelse(is.na(folds),cv.glmnet(x,y,nfolds=n),cv.glmnet(x,y,nfolds=n,foldid=folds))[[1]]
  #model <- cv.glmnet(x,y,nfolds=n,alpha=0)
  model <- cv.glmnet(x,y,nfolds=n,foldid=folds,alpha=0,lambda = lam)
  return(model)
}

index_to_folds<-function(x){
  n<-length(x)
  n2<-length(unlist(x))
  folds<-rep(NA,n2)
  for (i in 1:n){
    folds[x[[i]]]<-i
  }
  return(folds)
}

index_to_folds2 <- function(x,y){
  n<-nrow(y)
  folds<-rep(NA,n)
  for (i in 1:10){
  folds[i]<-grep(paste0("\\b",i,"\\b"),x)}
  return(folds)
}

ridge_params <- function(x_tr,vars,lambda){
  x<-x_tr[,vars]
  x$Type <-ifelse(x$Type=="PSP",1,0)
  x$logheight<-log(x$height.x)
  y<- ihs(x_tr$dh)
  mod_temp <- glmnet(x,y,alpha=0,lambda = lambda)
  return(mod_temp)
}

# for penalized (penalized package)

cv_fun_pen<-function(data,index,formula,lambda){
  # make n variable, looking towards smaller species
  n<-length(index)
  results<-rep(NA,n)
  data$Type<-as.factor(data$Type)
  for (i in 1:n){
    val_temp <- data[index[[i]],]
    train_temp <- data[-index[[i]],]
    #lambda <- (lambda/2)*nrow(train_temp)
    # function fails with NA response, need to remove
    #penalized(DBH.ihs~DBH.x * BAL + dt + SMImean,data=tempdat,lambda1=0,lambda2=0.03,standardize = T)
    mod_temp <- penalized(as.formula(formula),data=train_temp,lambda1=0,lambda2=lambda,standardize=T)
    results[i]<-sum((ihs(val_temp$dh) -predict(mod_temp,data=val_temp)[,1])^2,na.rm=T)/nrow(val_temp)
  }
  sum(results)
}


# function to automate folds by summer mean temperature:
cv_climate<-function(x,var,n=10){
  #tempsort<-x %>% group_by(Site) %>% summarize(mean(MeanTempJuneToAugust))
  tempbreaks<-quantile(x[[var]], probs = seq(0.1, 0.9, by = 1/n))
  cv_groups_temp<-list()
  # for (i in 1:10) {
  #   cv_groups_temp[[i]] <- tempsort$Site[
  #     which(tempsort$`mean(MeanTempJuneToAugust)`>ifelse(i-1>0,tempbreaks[i-1],0) & tempsort$`mean(MeanTempJuneToAugust)`< ifelse(i<10,tempbreaks[i],20))]
  # }
  cv_folds_temp<-list()
  for (i in 1:n){
    cv_folds_temp[[i]]<-which(x[[var]] >=ifelse(i-1>0,tempbreaks[i-1],-5) & x[[var]]< ifelse(i<10,tempbreaks[i],25))
  } 
  #return(assign(paste0("cv_folds_","y"),cv_folds_temp))
  return(cv_folds_temp)
  #return(cv_groups_temp)
}


ridge_params <- function(x_tr,vars,lambda){
  x<-x_tr[,vars]
  x$Type <-ifelse(x$Type=="PSP",1,0)
  x$logheight<-log(x$height.x)
  y<- ihs(x_tr$dh)
  mod_temp <- glmnet(x,y,alpha=0,lambda = lambda)
  return(mod_temp)
}

ridge_error <- function(x_tr,x_val,vars,lambda,type="MSE"){
  dat <- x_clean(x_tr,vars)
  x<-dat[[1]]
  y<-dat[[2]]
  val<- x_clean(x_val,vars)
  xval<-val[[1]]
  yval<-val[[2]]
  mod_temp <- glmnet(x,y,alpha=0,lambda = lambda)
  xval<-as.matrix(xval)
  if (type == "MSE"){
  result<-sum((yval -predict(mod_temp,xval))^2,na.rm=T)/nrow(xval)
  } else if (type == "bias"){
    result<-sum((yval -predict(mod_temp,xval)),na.rm=T)/nrow(xval)
  } else {
    print("Error: type is not MSE or bias")
  }
  return(result)
}

x_clean<-function(x_tr,vars){
  x<-x_tr[,vars]
  x$Type <-ifelse(x$Type=="PSP",1,0)
  x$logheight<-log(x$height.x)
  y<- ihs(x_tr$dh)
  na.list<-which(is.na(y))
  na.list<-which(is.na(y))
  if(any(na.list)){
    x<-x[-na.list,]
    y<-y[-na.list]}
  return(list(x,y))
}
#head(x[var])

# cv_climate_vec<-function(x,var,n=10){
#   #tempsort<-x %>% group_by(Site) %>% summarize(mean(MeanTempJuneToAugust))
#   tempbreaks<-quantile(x[[var]], probs = seq(0.1, 0.9, by = 1/n))
#   cv_folds_temp<-rep(NA,nrow(x))
#   for (i in 1:n){
#     cv_folds_temp[which(x[[var]] >= ifelse(i-1>0,tempbreaks[i-1],-5) & x[[var]]< ifelse(i<10,tempbreaks[i],25))]<-i
#   } 
#   #return(assign(paste0("cv_folds_","y"),cv_folds_temp))
#   return(cv_folds_temp)
#   #return(cv_groups_temp)
# }

#### Prediction functions ####

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

#### No intercept functions: ####

cv_function_rni<-function(data,index,variables,lambda){
  n<-length(index)
  results<-rep(NA,n)
  data$Type <-ifelse(data$Type=="PSP",1,0)
  #data$Type<-as.factor(data$Type)
  for (i in 1:n){
    train_temp <- data[-index[[i]],]
    x<-train_temp[,variables]
    x$logheight<-log(x$height.x)
    #x$int<-x$DBH.x*x$BAL
    y<- ihs(train_temp$dh)
    val_temp <- data[index[[i]],]
    xval<- val_temp[,variables]
    xval$logheight<-log(xval$height.x)
    #xval$int<-xval$DBH.x*xval$BAL
    xval<-as.matrix(xval)
    mod_temp <- glmnet(x,y,alpha=0,lambda = lambda,intercept = FALSE)
    results[i]<-sum((ihs(val_temp$dh) -predict(mod_temp,xval))^2,na.rm=T)/nrow(val_temp)
  }
  sum(results)
}

predict_pg_rni<-function(resmat,name,vars,lambda1=0.01,n=16){
  for (i in 1:n){
    x<-ht_tr[,vars]
    x$Type <-ifelse(x$Type=="PSP",1,0)
    x$logheight<-log(x$height.x)
    y<- ihs(ht_tr$dh)
    mod_temp <- glmnet(x,y,alpha=0,lambda=lambda1,intercept=FALSE)
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
