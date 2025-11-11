# Import data and libraries ----

rm(list=ls())

library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(colorspace)
library(MASS)
library(foreign)
library(nnet)
library(margins)
library(lmtest)
library(glmmTMB)
library(writexl)
library(Hmisc)
library(socialmixr)
pacman::p_load(
  rio,           # import/export
  tidyverse,     # data mgmt and viz
  naniar,        # assess and visualize missingness
  mice           # missing data imputation
)

set.seed(4)

load("data_in/datamod.Rdata")
data_mod_children <- data_mod_children[data_mod_children$total_contacts_prol <= 100 & data_mod_children$total_contacts_prol!=0,] # censoring at 100
data_mod <- data_mod[data_mod$total_contacts_prol <= 100,] # censoring at 100
data_mod_tot <- rbind(data_mod, data_mod_children)

data_mod_tot$age_group_10 <- as.character(data_mod_tot$age_group_10)
data_mod_tot$age_group_10[data_mod_tot$respondent_age <= 9] <- '0-9 y'
data_mod_tot$age_group_10[data_mod_tot$respondent_age >= 10 & data_mod_tot$respondent_age < 20] <- '10-19 y'
data_mod_tot$age_group_10[data_mod_tot$age_group_10 == "18-29 y"] <- '20-29 y'

data_mod_tot <- data_mod_tot %>%
  group_by(EPID) %>%
  mutate(respondent_type = ifelse(dplyr::n_distinct(wave) > 1, "longitudinal", "panel")) %>%
  ungroup() %>%
  as.data.frame()

my_boot = function(x, times=10000) {
  # Bootstrap 95% CI
  cis = quantile(replicate(times, mean(sample(x, replace=TRUE))), probs=c(0.025,0.975))
  # Return results as a data frame
  data.frame(mean=mean(x), lower.ci=cis[1], upper.ci=cis[2])
}




# -----------------
# Table 1 main ----
varlist <- c("age_group_10", "respondent_gender", "income_threecat","hh_size","region_grouped_IT","vacc_covid_bin")

tot_table1 <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  data_temp <- data_mod_tot
  temp1 <- as.data.frame(table(data_temp[,var], data_temp[,"wave"], useNA="ifany"))
  temp2 <- as.data.frame(prop.table(table(data_temp[,var], data_temp[,"wave"], useNA="ifany"),2))
  temp1 <- spread(temp1, key = Var2, value = Freq)
  temp2 <- spread(temp2, key = Var2, value = Freq)
  temp <- cbind(temp1, temp2[2:3])
  colnames(temp) <- c("Var", "freq1", "freq2", "prop1", "prop2")
  temp$info1 <- paste0(temp$freq1, " (", round(temp$prop1*100, 1), "%)")
  temp$info2 <- paste0(temp$freq2, " (", round(temp$prop2*100, 1), "%)")
  temp <- temp %>% dplyr::select(Var, info1, info2)
  rm(temp1)
  rm(temp2)
  var2 <- sym(var)

    tdata <- data_temp %>%
    group_by(wave, !!var2) %>%
    do(my_boot(.$total_contacts_prol))
    
    tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")
    
  tdata <- tdata %>% dplyr::select(!!var2, wave, info)  
  tdata <- spread(tdata, key = wave, value = info)
  
  temp <- cbind(temp, tdata[,2:3])
  temp <- rbind(c(rep(NA, length(temp))), temp)

  temp <-  temp[,c(1,2,4,3,5)]
  tot_table1 <- rbind(tot_table1, temp)
}


## by wave
temp <- as.data.frame(table(data_mod_tot[,"wave"]))
temp <- t(as.matrix(temp$Freq))

tdata <- data_mod_tot %>%
  group_by(wave) %>%
  do(my_boot(.$total_contacts_prol))

tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")

tdata <- tdata %>% dplyr::select(wave, info)  
tdata <- spread(tdata, key = wave, value = info)

temp <- cbind("wave",temp, tdata)
temp <- rbind(c(rep(NA, length(temp))), temp)

wave_table1 <-  temp[,c(1,2,4,3,5)]
colnames(wave_table1) <- colnames(tot_table1)

# final table by wave
final_tablewave <- rbind(wave_table1, tot_table1)

# table total
table_totals <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  data_temp <- data_mod_tot
  temp1 <- as.data.frame(table(data_temp[,var], useNA="ifany"))
  temp2 <- as.data.frame(prop.table(table(data_temp[,var], useNA="ifany")))
  temp <- cbind(temp1, temp2[2])
  colnames(temp) <- c("Var", "freqtot", "proptot")
  temp$infotot <- paste0(temp$freqtot, " (", round(temp$proptot*100, 1), "%)")
  temp <- temp %>% dplyr::select(Var, infotot)
  rm(temp1)
  rm(temp2)
  var2 <- sym(var)
  
  tdata <- data_temp %>%
    group_by( !!var2) %>%
    do(my_boot(.$total_contacts_prol))
  
  tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")
  
  tdata <- tdata %>% dplyr::select(!!var2, info)
  
  temp <- cbind(temp, tdata[,2])
  temp <- rbind(c(rep(NA, length(temp))), temp)
  
  table_totals <- rbind(table_totals, temp)
}


## totals
temp <- as.data.frame(nrow(data_mod_tot))

tdata <- data_mod_tot %>%
  do(my_boot(.$total_contacts_prol))

tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")

tdata <- tdata %>% dplyr::select(info)  

temp <- cbind("wave",temp, tdata)
temp <- rbind(c(rep(NA, length(temp))), temp)

colnames(temp) <- colnames(table_totals)

## rbind
final_tabletot <- rbind(temp, table_totals)

## Final table 
final_table1 <- cbind(final_tablewave, final_tabletot)
final_table1 <- final_table1[,-c(6)]
write.table(final_table1, "data_out/tables/main_table1.txt",row.names = FALSE)


# -----------------
# Table 2 main ----
varlist <- c("education", "occupation_agg", "worked_in_presence", "schooled_in_presence")
data_mod$nowave <- 1
data_mod$schooled_in_presence <- ifelse(data_mod$presence_school=='In presence','Yes','No')
data_mod_children$schooled_in_presence <- ifelse(data_mod_children$presence_school=='In presence','Yes','No')

tot_table2 <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  
  if(var == "presence_school" || var == "schooled_in_presence"){
    data_temp <- data_mod %>% filter(occupation_agg == "Student" & weekend == "No" & respondent_age > 17)
  } else if(var == "worked_in_presence"){
    data_temp <- data_mod %>% filter(occupation_agg == "Employed" & weekend == "No" & respondent_age > 17)
  } else {
    data_temp <- data_mod %>% filter(respondent_age >= 17)
  }
  
  temp1 <- as.data.frame(table(data_temp[,var], data_temp[,"wave"], useNA="ifany"))
  temp2 <- as.data.frame(prop.table(table(data_temp[,var], data_temp[,"wave"], useNA="ifany"),2))
  temp1 <- spread(temp1, key = Var2, value = Freq)
  temp2 <- spread(temp2, key = Var2, value = Freq)
  temp <- cbind(temp1, temp2[2:3])
  colnames(temp) <- c("Var", "freq1", "freq2", "prop1", "prop2")
  temp$info1 <- paste0(temp$freq1, " (", round(temp$prop1*100, 1), "%)")
  temp$info2 <- paste0(temp$freq2, " (", round(temp$prop2*100, 1), "%)")
  temp <- temp %>% dplyr::select(Var, info1, info2)
  rm(temp1)
  rm(temp2)
  var2 <- sym(var)
  
  tdata <- data_temp %>%
    group_by(wave, !!var2) %>%
    do(my_boot(.$total_contacts_prol))
  
  tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")
  
  tdata <- tdata %>% dplyr::select(!!var2, wave, info)  
  tdata <- spread(tdata, key = wave, value = info)
  
  temp <- cbind(temp, tdata[,2:3])
  temp <- rbind(c(rep(NA, length(temp))), temp)
  
  temp <-  temp[,c(1,2,4,3,5)]
  tot_table2 <- rbind(tot_table2, temp)
}

## by wave
temp <- as.data.frame(table(data_mod[,"wave"]))
temp <- t(as.matrix(temp$Freq))

tdata <- data_mod %>%
  group_by(wave) %>%
  do(my_boot(.$total_contacts_prol))

tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")

tdata <- tdata %>% dplyr::select(wave, info)  
tdata <- spread(tdata, key = wave, value = info)

temp <- cbind("wave",temp, tdata)
temp <- rbind(c(rep(NA, length(temp))), temp)

wave_table2 <-  temp[,c(1,2,4,3,5)]
colnames(wave_table2) <- colnames(tot_table2)

## final table by wave
final_tablewave <- rbind(wave_table2, tot_table2)

# table total
table_totals <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  
  if(var == "presence_school" || var == "schooled_in_presence"){
    data_temp <- data_mod %>% filter(occupation_agg == "Student" & weekend == "No" & respondent_age > 17)
  } else if(var == "worked_in_presence"){
    data_temp <- data_mod %>% filter(occupation_agg == "Employed" & weekend == "No" & respondent_age > 17)
  } else {
    data_temp <- data_mod %>% filter(respondent_age >= 17)
  }
  
  temp1 <- as.data.frame(table(data_temp[,var], useNA="ifany"))
  temp2 <- as.data.frame(prop.table(table(data_temp[,var], useNA="ifany")))
  temp <- cbind(temp1, temp2[2])
  colnames(temp) <- c("Var", "freqtot", "proptot")
  temp$infotot <- paste0(temp$freqtot, " (", round(temp$proptot*100, 1), "%)")
  temp <- temp %>% dplyr::select(Var, infotot)
  rm(temp1)
  rm(temp2)
  var2 <- sym(var)
  
  tdata <- data_temp %>%
    group_by( !!var2) %>%
    do(my_boot(.$total_contacts_prol))
  
  tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")
  
  tdata <- tdata %>% dplyr::select(!!var2, info)
  
  temp <- cbind(temp, tdata[,2])
  temp <- rbind(c(rep(NA, length(temp))), temp)
  
  table_totals <- rbind(table_totals, temp)
}


## totals
temp <- as.data.frame(nrow(data_mod))

tdata <- data_mod %>%
  do(my_boot(.$total_contacts_prol))

tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")

tdata <- tdata %>% dplyr::select(info)  

temp <- cbind("wave",temp, tdata)
temp <- rbind(c(rep(NA, length(temp))), temp)

colnames(temp) <- colnames(table_totals)

## rbind
final_tabletot <- rbind(temp, table_totals)

## Final table 
final_table2 <- cbind(final_tablewave, final_tabletot)
final_table2 <- final_table2[,-c(6)]
write.table(final_table2, "data_out/tables/main_table2.txt",row.names = FALSE)


### --- Table TOTAL main numbers ----
varlist <- c("age_group_10", "respondent_gender", "income_threecat","hh_size")
tot_table1_no_wave <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  data_temp <- data_mod_tot
  temp1 <- as.data.frame(table(data_temp[,var], useNA="ifany"))
  temp2 <- as.data.frame(prop.table(table(data_temp[,var], useNA="ifany")))
  temp <- cbind(temp1, temp2[2])
  colnames(temp) <- c("Var", "freq", "prop")
  temp$info <- paste0(temp$freq, " (", round(temp$prop*100, 1), "%)")
  temp <- temp %>% dplyr::select(Var, info)
  var2 <- sym(var)
  tdata <- data_temp %>%
    group_by(!!var2) %>%
    do(my_boot(.$total_contacts_prol))
  tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")
  tdata <- tdata %>% dplyr::select(!!var2, info)  
  temp <- cbind(temp, tdata[,2])
  temp <- rbind(c(rep(NA, length(temp))), temp)
  temp <-  temp[,c(1,2,3)]
  tot_table1_no_wave <- rbind(tot_table1_no_wave, temp)
}
## -- Table2 --
varlist <- c("vacc_covid_bin", "d_time_since_covid","education", "occupation_agg", "worked_in_presence", "schooled_in_presence")
data_mod$nowave <- 1
data_mod$schooled_in_presence <- ifelse(data_mod$presence_school=='In presence','Yes','No')
tot_table2_no_wave <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  if(var == "presence_school" || var == "schooled_in_presence"){
    data_temp <- data_mod %>% filter(occupation_agg == "Student" & weekend == "No" & respondent_age > 17)
  } else if(var == "worked_in_presence"){
    data_temp <- data_mod %>% filter(occupation_agg == "Employed" & weekend == "No" & respondent_age > 17)
  } else {
    data_temp <- data_mod %>% filter(respondent_age >= 17)
    # data_temp['occupation_agg'] <- as.factor(ifelse(data_temp$respondent_age>15,data_temp$occupation_agg, ifelse(data_temp$respondent_age>5,"Student", ifelse(data_temp$kindergarten_05=="Yes","Student","Inactive"))))
                                        
  }
  temp1 <- as.data.frame(table(data_temp[,var], useNA="ifany"))
  temp2 <- as.data.frame(prop.table(table(data_temp[,var], useNA="ifany")))
  temp <- cbind(temp1, temp2[2])
  colnames(temp) <- c("Var", "freq", "prop")
  temp$info <- paste0(temp$freq, " (", round(temp$prop*100, 1), "%)")
  temp <- temp %>% dplyr::select(Var, info)
  var2 <- sym(var)
  tdata <- data_temp %>%
    group_by(!!var2) %>%
    do(my_boot(.$total_contacts_prol))
  tdata$info <- paste0(round(tdata$mean,1), " (", round(tdata$lower.ci,1), " - ", round(tdata$upper.ci,1), ") ")
  tdata <- tdata %>% dplyr::select(!!var2, info)  
  temp <- cbind(temp, tdata[,2])
  temp <- rbind(c(rep(NA, length(temp))), temp)
  temp <-  temp[,c(1,2,3)]
  tot_table2_no_wave <- rbind(tot_table2_no_wave, temp)
}

#### ----



# Table SI2 1 ----

varlist <- c("age_group_10", "respondent_gender", "income_threecat","hh_size","region_grouped_IT")
dir_indir_table1 <- NULL
for(i in 1:length(varlist)){
  
  var <- varlist[i]
  data_temp <- data_mod_tot
  var2 <- sym(var)
  
  tdata1 <- data_temp %>%
    group_by(wave, !!var2) %>%
    do(my_boot(.$total_contacts))
  
  tdata1$info <- paste0(round(tdata1$mean,1), " (", round(tdata1$lower.ci,1), " - ", round(tdata1$upper.ci,1), ") ")
  
  tdata1 <- tdata1 %>% dplyr::select(!!var2, wave, info)  
  tdata1 <- spread(tdata1, key = wave, value = info)
  
  tdata2 <- data_temp %>%
    group_by(wave, !!var2) %>%
    do(my_boot(.$c_sharedindoor))
  
  tdata2$info <- paste0(round(tdata2$mean,1), " (", round(tdata2$lower.ci,1), " - ", round(tdata2$upper.ci,1), ") ")
  
  tdata2 <- tdata2 %>% dplyr::select(!!var2, wave, info)  
  tdata2 <- spread(tdata2, key = wave, value = info)
  
  temp <- cbind(tdata1, tdata2[,2:3])
  colnames(temp) <- c("level", "w1dir", "w2dir", "w1ind", "w2ind")
  header <- as.data.frame(t(rep(NA, length(temp))))
  colnames(header) <- c("level", "w1dir", "w2dir", "w1ind", "w2ind")
  temp <- rbind(header, temp)
  
  temp <-  temp[,c(1,2,4,3,5)]
  dir_indir_table1 <- rbind(dir_indir_table1, temp)
  
}

## by wave 

tdata1 <- data_mod_tot %>%
  group_by(wave) %>%
  do(my_boot(.$total_contacts))

tdata1$info <- paste0(round(tdata1$mean,1), " (", round(tdata1$lower.ci,1), " - ", round(tdata1$upper.ci,1), ") ")

tdata1 <- tdata1 %>% dplyr::select(wave, info)  
tdata1 <- spread(tdata1, key = wave, value = info)

tdata2 <- data_mod_tot %>%
  group_by(wave) %>%
  do(my_boot(.$c_sharedindoor))

tdata2$info <- paste0(round(tdata2$mean,1), " (", round(tdata2$lower.ci,1), " - ", round(tdata2$upper.ci,1), ") ")

tdata2 <- tdata2 %>% dplyr::select(wave, info)  
tdata2 <- spread(tdata2, key = wave, value = info)

temp <- cbind("wave",tdata1, tdata2)
colnames(temp) <- c("level", "w1dir", "w2dir", "w1ind", "w2ind")
header <- as.data.frame(t(rep(NA, length(temp))))
colnames(header) <- c("level", "w1dir", "w2dir", "w1ind", "w2ind")
temp <- rbind(header, temp)

temp <-  temp[,c(1,2,4,3,5)]

## Final table
dir_indir_table1 <- rbind(temp, dir_indir_table1)
write.table(dir_indir_table1, "data_out/tables/SI2_table1.txt",row.names = FALSE)

# Table SI2 2 ----
varlist <- c("vacc_covid_bin", "covid",
             "education", "occupation_agg", "worked_in_presence", "presence_school")

dir_indir_table2 <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  
  if(var == "presence_school"){
    data_temp <- data_mod %>% filter(occupation_agg == "Student" & weekend == "No")
  } else if(var == "worked_in_presence"){
    data_temp <- data_mod %>% filter(occupation_agg == "Employed" & weekend == "No")
  } else {
    data_temp <- data_mod
  }
  
  var2 <- sym(var)
  
  tdata1 <- data_temp %>%
    group_by(wave, !!var2) %>%
    do(my_boot(.$total_contacts))
  
  tdata1$info <- paste0(round(tdata1$mean,1), " (", round(tdata1$lower.ci,1), " - ", round(tdata1$upper.ci,1), ") ")
  
  tdata1 <- tdata1 %>% dplyr::select(!!var2, wave, info)  
  tdata1 <- spread(tdata1, key = wave, value = info)
  
  tdata2 <- data_temp %>%
    group_by(wave, !!var2) %>%
    do(my_boot(.$c_sharedindoor))
  
  tdata2$info <- paste0(round(tdata2$mean,1), " (", round(tdata2$lower.ci,1), " - ", round(tdata2$upper.ci,1), ") ")
  
  tdata2 <- tdata2 %>% dplyr::select(!!var2, wave, info)  
  tdata2 <- spread(tdata2, key = wave, value = info)
  
  temp <- cbind(tdata1, tdata2[,2:3])
  colnames(temp) <- c("level", "w1dir", "w2dir", "w1ind", "w2ind")
  header <- as.data.frame(t(rep(NA, length(temp))))
  colnames(header) <- c("level", "w1dir", "w2dir", "w1ind", "w2ind")
  temp <- rbind(header, temp)
  
  temp <-  temp[,c(1,2,4,3,5)]
  dir_indir_table2 <- rbind(dir_indir_table2, temp)
  
}

## by wave

tdata1 <- data_mod %>%
  group_by(wave) %>%
  do(my_boot(.$total_contacts))

tdata1$info <- paste0(round(tdata1$mean,1), " (", round(tdata1$lower.ci,1), " - ", round(tdata1$upper.ci,1), ") ")

tdata1 <- tdata1 %>% dplyr::select(wave, info)  
tdata1 <- spread(tdata1, key = wave, value = info)

tdata2 <- data_mod %>%
  group_by(wave) %>%
  do(my_boot(.$c_sharedindoor))

tdata2$info <- paste0(round(tdata2$mean,1), " (", round(tdata2$lower.ci,1), " - ", round(tdata2$upper.ci,1), ") ")

tdata2 <- tdata2 %>% dplyr::select(wave, info)  
tdata2 <- spread(tdata2, key = wave, value = info)

temp <- cbind("wave",tdata1, tdata2)
colnames(temp) <- c("level", "w1dir", "w2dir", "w1ind", "w2ind")
header <- as.data.frame(t(rep(NA, length(temp))))
colnames(header) <- c("level", "w1dir", "w2dir", "w1ind", "w2ind")
temp <- rbind(header, temp)

temp <-  temp[,c(1,2,4,3,5)]

## Final table
dir_indir_table2 <- rbind(temp, dir_indir_table2)
write.table(dir_indir_table2, "data_out/tables/SI2_table2.txt",row.names = FALSE)

# Table SI3 ----
contacts <- read.csv("data_in/clean_contacts_proc_prolonged_proportional_setting.csv")
clean_respondent_info_proc_v2 <- read.csv("data_in/clean_respondent_info_proc.csv")
clean_respondent_info_proc_v2 <- clean_respondent_info_proc_v2 %>% dplyr::select(caseid, EPID)
contacts <- left_join(contacts, clean_respondent_info_proc_v2, by = "caseid")
rm(clean_respondent_info_proc_v2)

contacts$id <- paste0(contacts$EPID, "_", contacts$wave)
data_mod_tot$id <- paste0(data_mod_tot$EPID, "_", data_mod_tot$wave)
contacts <- contacts %>% filter(id %in% data_mod_tot$id)

#Home=”home”+ “conviventi”+”Homeguest”, Work= “work”, School=”school”, Leisure=”Leisure”+”Restaurant”, 
#Transport=”Transport”, Other=”Shopping”+”Otherindoor”+”Otheroutdoor”

contacts$loc_env <- contacts$location
contacts$loc_env <- ifelse(contacts$loc_env == "home" | contacts$loc_env == "conviventi" | contacts$loc_env == "homeguest" , "home", 
                           ifelse(contacts$loc_env == "leisure" | contacts$loc_env == "restaurant", "leisure",
                                  ifelse(contacts$loc_env == "shopping" | contacts$loc_env == "otherindoor" | contacts$loc_env == "otheroutdoor", "other", contacts$loc_env)))

varlist <- c("home", "work", "school", "leisure", "transport", "other")

location_table <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  
  var2 <- sym(var)
  
  ## Direct
  tdata1 <- contacts %>% 
    filter(is.na(is_soft) & loc_env == !!var) %>%
    group_by(wave, id) %>% 
    count() 
  
  temp <- data_mod_tot %>% filter(!(id %in% tdata1$id)) %>% dplyr::select(wave, id)
  temp$wave <- as.integer(temp$wave)
  tdata1 <- rbind(tdata1, temp)
  tdata1$n <- replace_na(tdata1$n, 0)
  print(nrow(tdata1))
  
  tdata1 <- tdata1  %>%
    group_by(wave) %>%
    do(my_boot(.$n))
  
  tdata1$info <- paste0(round(tdata1$mean,1), " (", round(tdata1$lower.ci,1), " - ", round(tdata1$upper.ci,1), ") ")
  
  tdata1 <- tdata1 %>% dplyr::select(wave, info)  
  tdata1 <- spread(tdata1, key = wave, value = info)
  
  ## Indirect
  tdata2 <- contacts %>% 
    filter(!is.na(is_soft) & loc_env == !!var) %>%
    group_by(wave, id) %>% 
    count()
  
    temp <- data_mod_tot %>% filter(!(id %in% tdata2$id)) %>% dplyr::select(wave, id)
    temp$wave <- as.integer(temp$wave)
    tdata2 <- rbind(tdata2, temp)
    tdata2$n <- replace_na(tdata2$n, 0)
    print(nrow(tdata2))
    
  tdata2 <- tdata2 %>%
    group_by(wave) %>%
    do(my_boot(.$n))
  
  tdata2$info <- paste0(round(tdata2$mean,1), " (", round(tdata2$lower.ci,1), " - ", round(tdata2$upper.ci,1), ") ")
  
  tdata2 <- tdata2 %>% dplyr::select(wave, info)  
  tdata2 <- spread(tdata2, key = wave, value = info)
  
  ## Total
  tdata3 <- contacts %>% 
    filter((is_soft == FALSE | is.na(is_soft)) & loc_env == !!var) %>%
    group_by(wave, id) %>% 
    count()
  
  temp <- data_mod_tot %>% filter(!(id %in% tdata3$id)) %>% dplyr::select(wave, id)
  temp$wave <- as.integer(temp$wave)
  tdata3 <- rbind(tdata3, temp)
  tdata3$n <- replace_na(tdata3$n, 0)
  print(nrow(tdata3))
  
  tdata3 <- tdata3 %>%
    group_by(wave) %>%
    do(my_boot(.$n))
  
  tdata3$info <- paste0(round(tdata3$mean,1), " (", round(tdata3$lower.ci,1), " - ", round(tdata3$upper.ci,1), ") ")
  
  tdata3 <- tdata3 %>% dplyr::select(wave, info)  
  tdata3 <- spread(tdata3, key = wave, value = info)
  
  temp <- cbind(var, tdata1, tdata2, tdata3)
  colnames(temp) <- c("location", "dirw1", "dirw2", "indirw1", "indirw2", "totw1", "totw2")
  
  location_table <- rbind(location_table, temp)
  
}

## by wave

## Direct
  tdata4 <- contacts %>% 
    filter(is.na(is_soft)) %>%
    group_by(wave, id) %>% 
    count()

    temp <- data_mod_tot %>% filter(!(id %in% tdata4$id)) %>% dplyr::select(wave, id)
    temp$wave <- as.integer(temp$wave)
    temp$id <- temp$id
    tdata4 <- rbind(tdata4, temp)
    tdata4$n <- replace_na(tdata4$n, 0)
    print(nrow(tdata4))

    tdata4 <- tdata4 %>%
      group_by(wave) %>%
      do(my_boot(.$n))
    
  tdata4$info <- paste0(round(tdata4$mean,1), " (", round(tdata4$lower.ci,1), " - ", round(tdata4$upper.ci,1), ") ")
  
  tdata4 <- tdata4 %>% dplyr::select(wave, info)  
  tdata4 <- spread(tdata4, key = wave, value = info)

## Indirect
  tdata5 <- contacts %>% 
    filter(!is.na(is_soft)) %>%
    group_by(wave, id) %>% 
    count()
  
  temp <- data_mod_tot %>% filter(!(id %in% tdata5$id)) %>% dplyr::select(wave, id)
  temp$wave <- as.integer(temp$wave)
  # temp$id <- as.double(temp$id)
  tdata5 <- rbind(tdata5, temp)
  tdata5$n <- replace_na(tdata5$n, 0)
  print(nrow(tdata5))
  
  tdata5 <- tdata5 %>%
    group_by(wave) %>%
    do(my_boot(.$n))
  
  tdata5$info <- paste0(round(tdata5$mean,1), " (", round(tdata5$lower.ci,1), " - ", round(tdata5$upper.ci,1), ") ")
  
  tdata5 <- tdata5 %>% dplyr::select(wave, info)  
  tdata5 <- spread(tdata5, key = wave, value = info)

## Total
  tdata6 <- contacts %>% 
    filter(is_soft == FALSE | is.na(is_soft)) %>%
    group_by(wave, id) %>% 
    count()
  
  temp <- data_mod_tot %>% filter(!(id %in% tdata6$id)) %>% dplyr::select(wave, id)
  temp$wave <- as.integer(temp$wave)
  # temp$id <- as.double(temp$id)
  tdata6 <- rbind(tdata6, temp)
  tdata6$n <- replace_na(tdata6$n, 0)
  print(nrow(tdata5))
  
  tdata6 <- tdata6 %>%
    group_by(wave) %>%
    do(my_boot(.$n))
  
  tdata6$info <- paste0(round(tdata6$mean,1), " (", round(tdata6$lower.ci,1), " - ", round(tdata6$upper.ci,1), ") ")
  
  tdata6 <- tdata6 %>% dplyr::select(wave, info)  
  tdata6 <- spread(tdata6, key = wave, value = info)

wave_table <- cbind("total", tdata4, tdata5, tdata6)
colnames(wave_table) <- c("location", "dirw1", "dirw2", "indirw1", "indirw2", "totw1", "totw2")

location_table <- rbind(location_table, wave_table)
location_table <-  location_table[,c(1,2,4,6,3,5,7)]
write.table(location_table, "data_out/tables/SI3_table.txt",row.names = FALSE)

# Paragraph contacts ----

# employed adults over weekdays, inpresence
data_mod %>% filter(weekend == "No" & occupation_agg == "Employed") %>% group_by(worked_in_presence) %>% do(my_boot(.$total_contacts_prol))
# student
data_mod <- subset(data_mod, select = -nowave)
for_main <- rbind(data_mod_children, data_mod[data_mod$occupation_agg == "Student",])
for_main %>% filter(sunday == "No") %>% group_by(presence_school) %>% do(round(my_boot(.$total_contacts_prol),1))
for_main %>% filter(weekend == "No") %>% group_by(presence_school) %>% do(round(my_boot(.$total_contacts_prol),1))
data_mod %>% filter(weekend == "No" & occupation_agg == "Student") %>% group_by(presence_school) %>% do(round(my_boot(.$total_contacts_prol),1))
# vaccinated
data_mod %>% group_by(d_vacc2) %>% do(round(my_boot(.$total_contacts_prol),1))

# Comparison comix-polymod ----

### POLYMOD ----

tot <- polymod
part <- tot$participants
cont <- tot$contacts

IT_part <- part %>% filter(participant_nationality == "IT")
IT_cont <- cont %>% filter(part_id %in% IT_part$part_id)

IT_part$part_id[!(IT_part$part_id %in% IT_cont$part_id)] #4 non riportano contatti, parrebbe

data <- IT_cont %>% group_by(part_id) %>% count() 
to_append <- data.frame(matrix(c(50000, 0, 50001, 0, 50002, 0, 50003, 0), byrow = TRUE, nrow=4))
colnames(to_append) <- c("part_id", "n")
data <- rbind(data, to_append)

to_merge <- IT_part %>% dplyr::select(part_id, part_age, part_gender, hh_size, dayofweek)
data <- merge(data, to_merge, by = "part_id")

# creating variables
data$part_age_group <- data$part_age
data$part_age_group <- ifelse(data$part_age >= 0 & data$part_age <= 4, "0-4", 
                              ifelse(data$part_age >= 5 & data$part_age <= 17, "5-17", 
                                     ifelse(data$part_age >= 18 & data$part_age <= 29, "18-29", 
                                            ifelse(data$part_age >= 30  & data$part_age <= 39, "30-39", 
                                                   ifelse(data$part_age >= 40 & data$part_age <= 49, "40-49", 
                                                          ifelse(data$part_age >= 50 & data$part_age <= 59, "50-59", 
                                                                 ifelse(data$part_age >= 60, "60+", data$part_age_group)))))))
data$part_age_group <- factor(data$part_age_group, levels = c("0-4", "5-17", "18-29", 
                                                              "30-39", "40-49", "50-59", "60+"))

data$part_gender <- factor(data$part_gender, levels = c(NA, "F", "M"))
data$hh_size <- ifelse(data$hh_size >= 6, "6+", data$hh_size)
data$weekday <- data$dayofweek
data$weekday <- ifelse(data$dayofweek %in% c(0,6), "weekend", "weekday")

# by relevant variables

varlist <- c("part_age_group", "part_gender", "hh_size", "weekday")

polymod_tab <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  var2 <- sym(var)
  temp1 <- data %>% group_by(!!var2) %>% count()
  temp2 <- data %>% group_by(!!var2) %>% do(round(my_boot(.$n),1))
  temp <- merge(temp1, temp2)
  colnames(temp) <- c("Category",  "n", "Mean", "Lower_CI", "Upper_CI")
  polymod_tab <- rbind(polymod_tab, temp)
}

all1 <- nrow(data)
all2 <- data %>% do(round(my_boot(.$n),1))
all <- as.data.frame(c("All", all1, all2))
colnames(all) <- c("Category",  "n", "Mean", "Lower_CI", "Upper_CI")
polymod_tab <- rbind(all, polymod_tab)
polymod_tab$info <- paste0(polymod_tab$Mean, " (", polymod_tab$Lower_CI, " - ", polymod_tab$Upper_CI, ") ")
polymod_tab[, c("Mean", "Lower_CI", "Upper_CI")] <- NULL

# by setting
IT_cont$location <- NA
for(i in 1:length(IT_cont$location)){
  if(IT_cont$cnt_home[i] == 1){
    IT_cont$location[i] <- "home"
  } else if(IT_cont$cnt_home[i] == 0 & IT_cont$cnt_work[i] == 1){
    IT_cont$location[i] <- "work"
  } else if(IT_cont$cnt_home[i] == 0 & IT_cont$cnt_work[i] == 0 & 
            IT_cont$cnt_school[i] == 1){
    IT_cont$location[i] <- "school"
  } else if(IT_cont$cnt_home[i] == 0 & IT_cont$cnt_work[i] == 0 &
            IT_cont$cnt_school[i] == 0 & IT_cont$cnt_leisure[i] ==1){
    IT_cont$location[i] <- "leisure"
  } else if(IT_cont$cnt_home[i] == 0 & IT_cont$cnt_work[i] == 0 &
            IT_cont$cnt_school[i] == 0 & IT_cont$cnt_leisure[i] ==0 &
            IT_cont$cnt_transport[i] == 1){
    IT_cont$location[i] <- "transport"
  }
}
IT_cont$location <- ifelse(is.na(IT_cont$location), "other_NA", IT_cont$location)

table_setting_polymod <- NULL
data <- IT_cont %>% group_by(part_id, location) %>% count()
locations <- c("home", "work", "school", "leisure", "transport", "other_NA")
for(i in locations){
  aux <- data %>% filter(location == i)
  obs <- nrow(aux)
  m_ID <- IT_part$part_id[!(IT_part$part_id %in% aux$part_id)]
  to_bind <- data.frame(matrix(c(m_ID, rep(i, length(m_ID)), rep(0, length(m_ID))), byrow =FALSE, ncol=3))
  colnames(to_bind) <- c("part_id", "location", "n")
  to_bind$part_id <- as.numeric(to_bind$part_id)
  to_bind$n <- as.numeric(to_bind$n)
  to_bind$location <- as.character(to_bind$location)
  aux <- rbind(aux, to_bind)
  
  temp <- aux %>% group_by(!!i) %>% do(round(my_boot(.$n),1))
  temp <- as.data.frame(c(temp, obs))
  colnames(temp) <- c("Setting", "Mean", "Lower_CI", "Upper_CI", "n")
  table_setting_polymod <- rbind(table_setting_polymod, temp)
}

table_setting_polymod$info <- paste0(table_setting_polymod$Mean, " (", table_setting_polymod$Lower_CI, " - ", table_setting_polymod$Upper_CI, ") ")
table_setting_polymod[, c("Mean", "Lower_CI", "Upper_CI")] <- NULL

rm(aux, temp, obs, to_bind, part, cont, data, IT_cont, IT_part, tot, to_append, to_merge, all)

### CoMix ----
# CODICE INUTILMENTE LUNGO, EVENTUALMENTE DA RENDERE PIU' EFFICIENTE

tot <- get_survey("https://doi.org/10.5281/zenodo.5041112")
# part <- tot$participants
### patching chr part_age data from comix
part_common <- read.csv("data_in/comix_it/CoMix_it_participant_common.csv")
part_extra <- read.csv("data_in/comix_it/CoMix_it_participant_extra.csv")
part_hh <- read.csv("data_in/comix_it/CoMix_it_hh_common.csv")
part_sd <- read.csv("data_in/comix_it/CoMix_it_sday.csv")
part <- merge(part_common,part_extra, by='part_id')
part <- merge(part,part_sd, by='part_id')
part <- merge(part,part_hh, by='hh_id')

part$part_age_original <- part$part_age
x <- part$part_age
x <- gsub("Under 1", "0-1", x)
ranges <- regmatches(x, regexec("^(\\d+)-(\\d+)$", x))
lower  <- as.numeric(sapply(ranges, `[`, 2))
upper  <- as.numeric(sapply(ranges, `[`, 3))
part$part_age <- ifelse(!is.na(upper), (lower + upper)/2, lower)
###

cont <- tot$contacts

IT_part <- part %>% filter(survey_round %in% c(4,7,8,9)) #taking only round A6 and A7
IT_cont <- cont %>% filter(part_id %in% IT_part$part_id)

zero <- IT_part$part_id[!(IT_part$part_id %in% IT_cont$part_id)] 

# aggiungo waves bambini
IT_part$survey_round <- ifelse(IT_part$survey_round == 4, 7, IT_part$survey_round)
IT_part$survey_round <- ifelse(IT_part$survey_round == 8, 9, IT_part$survey_round)

data <- IT_cont %>% group_by(part_id) %>% count() 
to_append <- data.frame(matrix(c(zero, rep(0, length(zero))), byrow = FALSE, nrow=length(zero)))
colnames(to_append) <- c("part_id", "n")
data <- rbind(data, to_append)
to_merge <- IT_part %>% dplyr::select(part_id, part_age, part_gender, hh_size, dayofweek)
data <- merge(data, to_merge, by = "part_id")

all1 <- nrow(data)
all2 <- data %>% do(round(my_boot(.$n),1))
table_stat <- as.data.frame(c("All", all1, all2))

#table_stat <- data.frame(t(as.matrix(c("All", n_tot_1, mean_tot_1,  n_tot_2, mean_tot_2))))
colnames(table_stat) <- c("Category", "n", "Mean", "Lower_CI", "Upper_CI" )

# harmonize variables
data$part_age <- ifelse(data$part_age > 60, "60+", 
                        ifelse(data$part_age == 0.5 | data$part_age == 2.5, "2", 
                               ifelse(data$part_age == 8 | data$part_age == 13.5 | 
                          data$part_age == 16.5, "11", 
                            ifelse(data$part_age == 18.5 | data$part_age == 22 | 
                            data$part_age == 39.5| data$part_age == 49.5|
                            data$part_age == 59.5, NA, data$part_age))))
data$part_age <- as.character(data$part_age)
data$part_age <- factor(data$part_age, levels = c("2", "11", "23.5", "34.5", "44.5", "54.5", "60+"), ordered = TRUE)

data$weekday <- data$dayofweek
data$weekday <- ifelse(data$dayofweek %in% c(0,6), "weekend", "weekday")

data$hh_size <- ifelse(data$hh_size >= 6, "6+", data$hh_size)

# create table
varlist <- c("part_age", "part_gender", "hh_size", "weekday")

comix_tab <- NULL
for(i in 1:length(varlist)){
  var <- varlist[i]
  var2 <- sym(var)
  temp1 <- data %>% group_by(!!var2) %>% count()
  temp2 <- data %>% group_by(!!var2) %>% do(round(my_boot(.$n),1))
  temp <- merge(temp1, temp2)
  colnames(temp) <- c("Category",  "n", "Mean", "Lower_CI", "Upper_CI")
  comix_tab <- rbind(comix_tab, temp)
}
comix_tab <- rbind(table_stat, comix_tab)
comix_tab$info <- paste0(comix_tab$Mean, " (", comix_tab$Lower_CI, " - ", comix_tab$Upper_CI, ") ")
comix_tab[, c("Mean", "Lower_CI", "Upper_CI")] <- NULL

# by setting
IT_cont$location <- NA
for(i in 1:length(IT_cont$location)){
  if(IT_cont$cnt_home[i] == TRUE){
    IT_cont$location[i] <- "home"
  } else if(IT_cont$cnt_home[i] == FALSE & IT_cont$cnt_work[i] == TRUE){
    IT_cont$location[i] <- "work"
  } else if(IT_cont$cnt_home[i] == FALSE & IT_cont$cnt_work[i] == FALSE & 
            IT_cont$cnt_school[i] == TRUE){
    IT_cont$location[i] <- "school"
  }  else if(IT_cont$cnt_home[i] == FALSE & IT_cont$cnt_work[i] == FALSE &
             IT_cont$cnt_school[i] == FALSE & IT_cont$cnt_leisure[i] == TRUE){
    IT_cont$location[i] <- "leisure"
  } else if(IT_cont$cnt_home[i] == FALSE & IT_cont$cnt_work[i] == FALSE &
            IT_cont$cnt_school[i] == FALSE & IT_cont$cnt_leisure[i] == FALSE &
            IT_cont$cnt_transport[i] == TRUE){
    IT_cont$location[i] <- "transport"
  }
}

IT_cont$location <- ifelse(is.na(IT_cont$location), "other_NA", IT_cont$location)
locations <- c("home", "work", "school", "leisure", "transport", "other_NA")
table_setting_comix <- NULL

data <- IT_cont %>% group_by(part_id, location) %>% count()

for(i in locations){
  aux <- data %>% filter(location == i)
  obs <- nrow(aux)
  m_ID <- IT_part$part_id[!(IT_part$part_id %in% aux$part_id)]
  to_bind <- data.frame(matrix(c(m_ID, rep(i, length(m_ID)), rep(0, length(m_ID))), byrow =FALSE, ncol=3))
  colnames(to_bind) <- c("part_id", "location", "n")
  to_bind$part_id <- as.numeric(to_bind$part_id)
  to_bind$n <- as.numeric(to_bind$n)
  to_bind$location <- as.character(to_bind$location)
  aux <- rbind(aux, to_bind)
  
  temp <- aux %>% group_by(!!i) %>% do(round(my_boot(.$n),1))
  temp <- as.data.frame(c(temp, obs))
  colnames(temp) <- c("Setting", "Mean", "Lower_CI", "Upper_CI", "n")
  table_setting_comix <- rbind(table_setting_comix, temp)
}

table_setting_comix$info <- paste0(table_setting_comix$Mean, " (", table_setting_comix$Lower_CI, " - ", table_setting_comix$Upper_CI, ") ")
table_setting_comix[, c("Mean", "Lower_CI", "Upper_CI")] <- NULL


### Our study ----
data_mod_tot$part_age <- ifelse(data_mod_tot$respondent_age >=0 & data_mod_tot$respondent_age < 5, "0-4", 
                         ifelse(data_mod_tot$respondent_age >= 5 & data_mod_tot$respondent_age < 18, "5-17",
                                ifelse(data_mod_tot$respondent_age >= 18 & data_mod_tot$respondent_age < 30, "18-29", 
                                       ifelse(data_mod_tot$respondent_age >= 30 & data_mod_tot$respondent_age < 40, "30-39",
                                              ifelse(data_mod_tot$respondent_age >= 40 & data_mod_tot$respondent_age < 50, "40-49",
                                                     ifelse(data_mod_tot$respondent_age >= 50 & data_mod_tot$respondent_age < 60, "50-59",
                                                     "60+"))))))
data_mod_tot$part_age <- factor(data_mod_tot$part_age, levels = c("0-4", "5-17", "18-29", "30-39", "40-49", "50-59", "60+"))


data_mod_tot$part_hh_size <- ifelse(data_mod_tot$hh_size_det >= 6, "6+", data_mod_tot$hh_size_det)
data_mod_tot$respondent_gender <- factor(data_mod_tot$respondent_gender, levels = c("Female", "Male"))

varlist <- c("part_age", "respondent_gender","part_hh_size", "weekend")
for(j in c(1,2)){
  temp_data <- data_mod_tot %>% filter(wave == j)
  
  all1 <- nrow(temp_data)
  all2 <- temp_data %>% do(round(my_boot(.$total_contacts_prol),1))
  all2$info <- paste0(all2$mean, " (", all2$lower.ci, " - ", all2$upper.ci, ") ")
  all2[, c("upper.ci", "lower.ci", "mean")] <- NULL
  all <- as.data.frame(c("All", all1, all2))
  colnames(all)<- c("Category", "n", "info")
  
  table_temp <- NULL
  name <- paste0("our_table_", j)
    for(i in 1:length(varlist)){
      var <- varlist[i]
      var2 <- sym(var)
      # tdata <- data_mod_tot %>%
      #   group_by(!!var2) %>%
      #   dplyr::summarise(count = n(), mean_c = round(mean(total_contacts_prol, na.rm=T), 1))
      # 
      # tdata <- cbind(tdata[1:(nrow(tdata)/2),], tdata[(nrow(tdata)/2+1):nrow(tdata), 3:4])
      # tdata <- tdata[,-1]
      # colnames(tdata) <- c("Var", "freq1", "mean1", "freq2", "mean2")
      
      temp1 <- temp_data %>% group_by(!!var2) %>% count()
      temp2 <- temp_data %>% group_by(!!var2) %>% do(round(my_boot(.$total_contacts_prol),1))
      temp <- merge(temp1, temp2)
      colnames(temp) <- c("Category", "n", "Mean", "Lower_CI", "Upper_CI")
      temp$info <- paste0(temp$Mean, " (", temp$Lower_CI, " - ", temp$Upper_CI, ") ")
      temp[, c("Mean", "Lower_CI", "Upper_CI")] <- NULL
    
      table_temp <- rbind(table_temp, temp)
    }
  table_temp <- rbind(all, table_temp)
  assign(name, table_temp)
}



our_table <- cbind(our_table_1, our_table_2)
our_table <- our_table[,-4]

# by setting
locations <- c("home", "work", "school", "leisure", "transport", "other")
table_setting_ourstudy <- NULL

data <- contacts %>% filter(is_soft == FALSE | is.na(is_soft)) %>% group_by(id, loc_env, wave) %>% count()

for(j in c(1,2)){
for(i in locations){
  aux <- data %>% filter(loc_env == i & wave == j)
  obs <- nrow(aux)
  wave <- j
  m_ID <- data_mod_tot$id[!(data_mod_tot$id %in% aux$id) & data_mod_tot$wave == j ]
  to_bind <- data.frame(matrix(c(m_ID, rep(i, length(m_ID)), rep(0, length(m_ID))), byrow =FALSE, ncol=3))
  colnames(to_bind) <- c("id", "loc_env", "n")
  to_bind$n <- as.numeric(to_bind$n)
  aux <- rbind(aux, to_bind)

  temp <- aux %>% group_by(!!i) %>% do(round(my_boot(.$n),1))
  temp <- as.data.frame(c(temp, obs, j))
  colnames(temp) <- c("Setting", "Mean", "Lower_CI", "Upper_CI", "n", "wave")
  table_setting_ourstudy <- rbind(table_setting_ourstudy, temp)
}
}
table_setting_ourstudy$info <- paste0(table_setting_ourstudy$Mean, " (", table_setting_ourstudy$Lower_CI, " - ", table_setting_ourstudy$Upper_CI, ") ")
table_setting_ourstudy[, c("Mean", "Lower_CI", "Upper_CI")] <- NULL
table_setting_ourstudy <- pivot_wider(table_setting_ourstudy, names_from = wave, values_from = c(n, info))
table_setting_ourstudy <- table_setting_ourstudy %>% dplyr::select(1,2,4,3,5)

## Tabelle totali ---
polymod_tab_fin <- rbind(polymod_tab[1:18,], rep(NA, 3), polymod_tab[19:20,])
our_table_fin <- rbind(our_table[1:8,], rep(NA,5), our_table[9:10,], rep(NA,5), our_table[11:16,], rep(NA,5), our_table[17:18,])
final_table <- cbind(polymod_tab_fin, comix_tab, our_table_fin)
final_table <- final_table[, -c(4,7)]

table_setting_fin <- cbind(table_setting_polymod,table_setting_comix, table_setting_ourstudy) 
table_setting_fin <- table_setting_fin[,-c(4,7)]

write.table(final_table, "data_out/tables/comixpolym_table1.txt",row.names = FALSE)
write.table(table_setting_fin, "data_out/tables/comixpolym_table2.txt",row.names = FALSE)



# ------------------------
# Table panel/long SI ----
varlist <- c("age_group_10", "respondent_gender", "income_threecat",
             "hh_size", "region_grouped_IT", "vacc_covid_bin")

make_block <- function(data, var) {
  # ---- COUNTS + PERC (denominator = wave x respondent_type) ----
  t3 <- table(data[[var]], data[["wave"]], data[["respondent_type"]], useNA = "ifany")
  cnt <- as.data.frame(t3)
  names(cnt) <- c("Var","wave","respondent_type","Freq")
  
  prp <- as.data.frame(prop.table(t3, margin = c(2,3)))
  names(prp) <- c("Var","wave","respondent_type","Prop")
  
  cntprp <- cnt %>%
    left_join(prp, by = c("Var","wave","respondent_type")) %>%
    mutate(info = paste0(Freq, " (", round(Prop*100, 1), "%)"),
           wave = as.character(wave)) %>%
    dplyr::select(Var, wave, respondent_type, info) %>%
    pivot_wider(names_from = c(wave, respondent_type), values_from = info, names_sep = "_")
  
  # ---- MEAN CONTACTS + 95% CI ----
  mn <- data %>%
    dplyr::group_by(wave, respondent_type, !!rlang::sym(var)) %>%
    do(my_boot(.$total_contacts_prol)) %>%
    ungroup() %>%
    mutate(info = paste0(round(mean,1), " (", round(lower.ci,1), " - ", round(upper.ci,1), ") "),
           wave = as.character(wave)) %>%
    transmute(Var = as.character(!!rlang::sym(var)), wave, respondent_type, info) %>%
    pivot_wider(names_from = c(wave, respondent_type), values_from = info, names_sep = "_")
  
  # ---- COMBINE & tidy names ----
  out <- cntprp %>% full_join(mn, by = "Var", suffix = c("_nperc", "_meanCI"))
  
  # Ensure all 8 expected columns exist (handles empty strata gracefully)
  need <- c("1_panel_nperc","1_longitudinal_nperc","2_panel_nperc","2_longitudinal_nperc",
            "1_panel_meanCI","1_longitudinal_meanCI","2_panel_meanCI","2_longitudinal_meanCI")
  for (nm in need) if (!nm %in% names(out)) out[[nm]] <- NA_character_
  
  # If pivot left unsuffixed names, normalize them
  norm <- function(df, old, new) if (old %in% names(df) && !(new %in% names(df))) dplyr::rename(df, !!new := dplyr::all_of(old)) else df
  out <- out %>%
    norm("1_panel", "1_panel_nperc") %>%
    norm("1_longitudinal", "1_longitudinal_nperc") %>%
    norm("2_panel", "2_panel_nperc") %>%
    norm("2_longitudinal", "2_longitudinal_nperc") %>%
    norm("1_panel.x", "1_panel_nperc") %>%
    norm("1_longitudinal.x", "1_longitudinal_nperc") %>%
    norm("2_panel.x", "2_panel_nperc") %>%
    norm("2_longitudinal.x", "2_longitudinal_nperc") %>%
    norm("1_panel.y", "1_panel_meanCI") %>%
    norm("1_longitudinal.y", "1_longitudinal_meanCI") %>%
    norm("2_panel.y", "2_panel_meanCI") %>%
    norm("2_longitudinal.y", "2_longitudinal_meanCI")
  
  # Final 8 columns (counts% then mean(CI))
  out %>%
    dplyr::select(Var,
                  `W1_panel_n%`              = `1_panel_nperc`,
                  `W1_panel_mean(CI)`        = `1_panel_meanCI`,
                  `W2_panel_n%`              = `2_panel_nperc`,
                  `W2_panel_mean(CI)`        = `2_panel_meanCI`,
                  `W1_longitudinal_n%`       = `1_longitudinal_nperc`,
                  `W1_longitudinal_mean(CI)` = `1_longitudinal_meanCI`,
                  `W2_longitudinal_n%`       = `2_longitudinal_nperc`,
                  `W2_longitudinal_mean(CI)` = `2_longitudinal_meanCI`)
}

final_table_8col <- dplyr::bind_rows(lapply(varlist, function(v) make_block(data_mod_tot, v)))
write.table(final_table_8col, "data_out/tables/panel_longitudinal_table1.txt",row.names = FALSE)

# ------------------------
