##################################
##################################
# PROCESSING DATA FOR STAT. MOD. #
##################################
##################################

# Import data and libraries #

rm(list=ls())

library(dplyr)
library(readr)
library(tidyr)
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
pacman::p_load(
  rio,           # import/export
  tidyverse,     # data mgmt and viz
  naniar,        # assess and visualize missingness
  mice           # missing data imputation
)
ego <- read.csv("data_in/clean_respondent_info_proc.csv")
contacts <- read.csv("data_in/clean_contacts_proc_prolonged_proportional_setting.csv")


# Prepare data ----
ego <- ego[!is.na(ego$isolation) & ego$isolation==4,] #then, keeping only people without restrictions

#ego$safe <- ifelse(!is.na(ego$perc_safe) & ego$perc_safe==1,1,
#                        ifelse(is.na(ego$perc_safe),NA,0))
#ego$zero <- ifelse(ego$total_contacts==0,1,0)

# #----------------------- Added 7th November 2025 Ele --------------------------#
# Dropped as already fixed in notebook 000.1 (Lore)
# ego$respondent_age <- ifelse(ego$respondent_sample == "recontact" & ego$wave == 2, ego$respondent_age + 1, ego$respondent_age)

#----------------------- Added 24th November 2023 Ele --------------------------#

# to_merge <- contacts %>% group_by(caseid) %>% count()
# colnames(to_merge) <- c("caseid", "total_contacts_prol")
# ego <- merge(ego, to_merge, by="caseid", all.x = T)
# ego$total_contacts_prol[is.na(ego$total_contacts_prol)] <- 0
# ego$total_contacts_nohh_prol <- ego$total_contacts_prol-ego$total_contacts_cohabitants

## Definition we had before December 10th ##
# ego$total_contacts_prol <- ego$total_contacts + ego$c_sharedindoor - ego$total_contacts_indoor
# ego$total_contacts_nohh_prol <- ego$total_contacts_prol-ego$outdoor_cohab

#-------------------------------------------------------------------------------#

# ------------------------Added December 10th 2023 Ele--------------------------#

#total_contacts_prol
ego$total_contacts_prol <- ego$total_contacts - ego$total_contacts_indoor + 
  pmax(ego$total_contacts_indoor, ego$c_sharedindoor)

# total_contacts_prol_soft
ego$total_contacts_prol_soft <- ego$total_contacts + ego$c_sharedindoor

#total_contacts_nohh_prol
#objective: tot_contacts_nohh_prol = contacts_nohh_outdoor + max(sharedindoor-(hh_size-1), total_contacts_nohh_indoor)

## contacts_nohh_outdoor= contatti outdoor - contatti hh outdoor
## = contatti totali - contatti indoor - contatti hh outdoor
ego$contacts_nohh_outdoor <- ego$total_contacts - ego$total_contacts_indoor - ego$outdoor_cohab

## c_sharedindoor_nohh 
ego$c_sharedindoor_nohh <- ego$c_sharedindoor - (ego$hh_size_det-1)

##total_contacts_nohh_indoor = contatti indoor - contatti hh indoor
## = contatti indoor - (contatti hh - contatti hh outdoor)
ego$contacts_nohh_indoor = ego$total_contacts_indoor - (ego$total_contacts_cohabitants - ego$outdoor_cohab)

## totale
ego$total_contacts_nohh_prol = ego$contacts_nohh_outdoor + pmax(ego$c_sharedindoor_nohh, ego$contacts_nohh_indoor)

ego$contacts_nohh_indoor = NULL
ego$c_sharedindoor_nohh = NULL
ego$contacts_nohh_outdoor = NULL

#-------------------------------------------------------------------------------#


ego$total_contacts_nohh <- ego$total_contacts-ego$total_contacts_cohabitants
ego$perc_risk <- ego$risky/ego$total_contacts_nohh
ego$perc_safe <- ego$masked_out/(ego$outdoor)
ego$perc_risk <- ifelse(is.nan(ego$perc_risk), NA, ego$perc_risk)
ego$risk <- ifelse(!is.na(ego$perc_risk) & ego$perc_risk==0,0,
                   ifelse(is.na(ego$perc_risk),NA,1))

ego$d_hh_size <- ifelse(ego$hh_size==1, "0", "1")

ego$d_children <- ifelse(ego$children_number >=1, "1", "0")
ego$d_senior65 <- ifelse(ego$senior65_cohabitant >=1, "1", "0")
ego$d_senior70 <- ifelse(ego$senior70_cohabitant >=1, "1", "0")
ego$d_senior75 <- ifelse(ego$senior75_cohabitant >=1, "1", "0")
ego$d_senior80 <- ifelse(ego$senior80_cohabitant >=1, "1", "0")

# replacing some values to zero for consistency: if household size=1, then no children, no cohabitants etc
ego$d_children<- ifelse(ego$hh_size_det ==1, 0, ego$d_children)
ego$d_senior65<- ifelse(ego$hh_size_det ==1, 0, ego$d_senior65)
ego$d_senior70<- ifelse(ego$hh_size_det ==1, 0, ego$d_senior70)
ego$d_senior75<- ifelse(ego$hh_size_det ==1, 0, ego$d_senior75)
ego$d_senior80<- ifelse(ego$hh_size_det ==1, 0, ego$d_senior80)
ego$chronic_comorb_cohab2<- ifelse(ego$hh_size_det ==1, "(-1,0]", ego$chronic_comorb_cohab2)


ego$income <- factor(ego$income)
levels(ego$income) = c("<1,000 eur","1,000-1,499 eur","1,500-1,999 eur", "2,000-2,499 eur","2,500-2,999 eur",
                       "3,000-3,499 eur", "3,500-3,999 eur", "4,000-4,999 eur", "5,000 eur","Doesn't answer")

ego$income_1000[is.na(ego$income_1000)] <- "Doesn't answer" 
ego$income_1000 = factor(ego$income_1000, levels =c("<1000","1000-1999","2000-2999","3000-3999","4000-4999",">5000","Doesn't answer"))
ego$income_1000 <- relevel(factor(ego$income_1000),ref= 1)


#----------------- Added 15th December ------------------------------#
ego$income_threecat = ifelse(ego$income=="<1,000 eur" | ego$income=="1,000-1,499 eur",1,
                                        ifelse(ego$income=="1,500-1,999 eur" | ego$income=="2,000-2,499 eur" | ego$income=="2,500-2,999 eur",2,
                                               ifelse(ego$income=="3,000-3,499 eur" | ego$income=="3,500-3,999 eur" | 
                                                        ego$income=="4,000-4,999 eur" | ego$income=="5,000 eur",3,
                                                      4)))
ego$income_threecat =as.factor(ego$income_threecat)
levels(ego$income_threecat)= c("<1,500 eur","1,500-2,999 eur",">3,000 eur","Doesn't answer")
#--------------------------------------------------------------------#

ego$perceived_income <- ifelse(ego$perceived_income ==933, NA, ego$perceived_income)
ego$perceived_income <- factor(ego$perceived_income)
ego$perceived_income <- relevel(factor(ego$perceived_income),ref= 3)


ego$occupation_agg <- ifelse(ego$occupation_agg==96,NA,
                             ego$occupation_agg)
ego$d_occupation <- ifelse((ego$occupation_agg == 1 |
                              ego$occupation_agg == 3) 
                           & !is.na(ego$occupation_agg), 1, 
                           ifelse(is.na(ego$occupation_agg), NA, 0))
ego$d_occupation <- relevel(factor(ego$d_occupation),ref= 2)
  
ego$occupation_agg2 <- ifelse(ego$occupation_agg == 1, 1,
                              ifelse(ego$occupation_agg == 3, 2,
                                     ifelse(ego$occupation_agg %in% c(2,4,5), 3,
                                            NA)))

ego$vacc_covid_bin <- as.factor(ifelse(ego$respondent_age<=5,4, ego$vacc_covid_bin)) #setting children to exempt
levels(ego$vacc_covid_bin) <- c("No doses", "One dose", "Two doses", "Three doses or more", "Exempt")

ego$vacc_covid_tot <- ifelse(is.na(ego$vacc_covid_1wave), 
                             ego$vacc_covid_2wave, 
                             ego$vacc_covid_1wave )
ego$vacc_covid_tot =factor(ego$vacc_covid_tot)

ego$vacc_covid_tot <- ifelse(ego$respondent_age<=5,1,ego$vacc_covid_tot) #setting children <5 to zero otherwise they will be dropped!
levels(ego$vacc_covid_tot) <- c("Non vaccinato per scelta",
                                "Esenzione per salute",
                                "Appuntamento prima dose",
                                "Prima dose",
                                "Seconda dose",
                                "Terza dose",
                                "Quarta dose",
                                "Quinta dose",
                                "Altro")

# ggplot(ego, aes(x= as.factor(vacc_covid_tot), y= total_contacts))+
#   geom_violin(draw_quantiles=c(0.25,0.5,0.75)) + ylim(0,20)

ego$vacc_covid_bin_2 <- ifelse(ego$vacc_covid_tot %in% c(1,2), 0,
                               ifelse(ego$vacc_covid_tot %in% c(3,4,5), 1,
                                      ifelse(ego$vacc_covid_tot %in% c(6,7,8), 2, NA )))
ego$vacc_covid_bin_2 <- relevel(factor(ego$vacc_covid_bin_2),ref= 2)
levels(ego$vacc_covid_bin_2) <- c("One/two doses", "No doses", "Three doses or more")


ego$d_vacc2 <- ifelse(ego$vacc_covid_tot <= 3 & ego$vacc_covid_tot != 999, 0, 
                      ifelse(ego$vacc_covid_tot == 999, NA, 1))
ego$d_vacc2 <- factor(ego$d_vacc2)
levels(ego$d_vacc2) <- c("No","Yes")
#ego$d_vacc2 <- relevel(factor(ego$d_vacc2),ref= 2)

ego$covid <- relevel(factor(ego$covid),ref= 2)

ego$time_since_covid_raw <- as.Date(ego$start_date_module2A)-as.Date(ego$covid_date)

ego$time_since_covid <- ifelse(ego$covid == 2, "N",
                                             ifelse(ego$covid == 1 & 
                                                      ego$time_since_covid_raw<120 & 
                                                      !is.na(ego$time_since_covid_raw), "Y_lessthan120days",
                                                    ifelse(ego$covid == 1 & 
                                                             ego$time_since_covid_raw>=120 &
                                                             !is.na(ego$time_since_covid_raw), "Y_120plusdaysago", 
                                                           ifelse(ego$covid == 1 & is.na(ego$time_since_covid_raw), "Y_nodate", NA
                                                           ))))
ego$d_time_since_covid <- ifelse(ego$time_since_covid=="Y_lessthan120days","Yes","No")

ego$age_group_10 <- ifelse(ego$respondent_age == 18 | ego$respondent_age == 19,
                           "20-29 y", ego$age_group_10)
ego$age_group_10 <- ifelse(ego$age_group_10 == "20-29 y", "18-29 y", ego$age_group_10)
ego$age_group_10 <- ifelse(ego$respondent_age<18, "0-17 y",ego$age_group_10)

ego$age_group_10 <- factor(ego$age_group_10, c("0-17 y","18-29 y", "30-39 y", "40-49 y",
                                               "50-59 y", "60-69 y", "70+ y"))

ego$age_group_25 <- cut(ego$respondent_age, breaks=c(-1,17,44,70,100), 
                        labels = c("0-17 y","18-44 y", "45-69 y", "70+ y"))

ego$refusal <- factor(ego$refusal)
ego$refusal <- relevel(ego$refusal, ref=2)
levels(ego$refusal) <- c("No","Yes","Doesn't know")


ego$other_vacc <- factor(ego$other_vacc)
ego$other_vacc <- relevel(ego$other_vacc, ref=2)
levels(ego$other_vacc) <- c("No","Yes","Doesn't know")

ego$kindergarten_05 <- ifelse(ego$kindergarten_05 ==2, 0,
                              ifelse(ego$kindergarten_05 ==1,1, NA))

# [presenza_work] In che modalità hai lavorato IERI? 
# <124>	Ho lavorato esclusivamente da remoto
# <125>	Ho lavorato del tutto o in parte in presenza
# <126>	Ieri non ho lavorato 
ego$worked_in_presence = ifelse(is.na(ego$presence_work) | ego$presence_work==1 | ego$presence_work==3,0,1) #set to 1 if "Ho lavorato del tutto o in parte in presenza"

#--------------------- Modified December 15th--------------------#
# [presenza_school] Che tipo di didattica hai effettuato IERI? 
# <127>	Esclusivamente da remoto (didattica a distanza/DAD)
# <128>	Del tutto o in parte in presenza
# <129>	Nessuna (ad es. la scuola/università era chiusa, ero malato/a) #TOSCRIPTING: FIXED

ego$presence_school <- factor(ego$presence_school)
levels(ego$presence_school) <- c("Remotely", "In presence","Didn't attend")
ego$presence_school = relevel(factor(ego$presence_school),ref= 2)
ego$presence_school[is.na(ego$presence_school)]="Didn't attend"

#---------------------------------------------------------------#

#create new variable for school or work
ego$attend_in_presence= ifelse(ego$worked_in_presence==1 | ego$presence_school=="In presence",1,0)
ego$attend_in_presence= factor(ego$attend_in_presence)
#ego$attend_in_presence = relevel(ego$attend_in_presence,ref= 2)
levels(ego$attend_in_presence)=c("No","Yes")

# [inf_sarscov]: Conosci personalmente almeno una persona che ha contratto, è stata ospedalizzata o è morta a causa di COVID-19
ego$knows_severe = ifelse(is.na(ego$inf_sarscov) | ego$inf_sarscov==4 | ego$inf_sarscov==1, 0, 1) # set to 1 if Conosco qualcuno che è stato ospedalizzato OR Conosco qualcuno che è morto


# ego$covid_hospital <- ifelse(ego$covid_hospital ==2, 0,
#                              ifelse(ego$covid_hospital ==1,1,NA)) # only 36 people hospitalized
# ego$covid_hospital[ego$covid==2] <- "No infection"


to_be_factor <- c("wave","region_grouped_IT", "respondent_gender", "education", 
                  "occupation_agg",  "d_occupation", "occupation_agg2", 
                  "sunday","weekend","worked_in_presence","presence_school",
                  "d_children", "d_senior65", "d_senior70", "d_senior75", "d_senior80",
                  "hh_size","d_hh_size",
                  "chronic_comorb_self2","chronic_comorb_cohab2",
                  "age_group_10","age_group_25", "covid", "vacc_covid_bin_2", #already factors with set reference levels, so this doesn't affect them
                  "d_vacc2", "vacc_covid_bin",
                  "other_vacc","refusal",
                  "kindergarten_05",
                  "perceived_income", "income","income_1000",
                  "d_time_since_covid","time_since_covid","knows_severe" #EPID removed for now
                   )
ego[to_be_factor] <- lapply(ego[to_be_factor], factor)

levels(ego$region_grouped_IT) <- c("North-West", "North-East", "Center", "South", "Islands")
levels(ego$respondent_gender) <- c("Male", "Female")
levels(ego$education) <- c("Elem/Middle school", "High school", "University") 
levels(ego$occupation_agg) <- c("Employed", "Home/family", "Student", "Retired", "Inactive")
levels(ego$d_occupation) <- c("Employed/Student","Other")
levels(ego$occupation_agg2) <- c("Employed", "Student", "Other (not working)")
levels(ego$sunday) <- c("No", "Yes")
levels(ego$weekend) <- c("No", "Yes")
levels(ego$worked_in_presence) <- c("No", "Yes")

levels(ego$d_children) <- c("No", "Yes")
levels(ego$d_senior65) <- c("No", "Yes")
levels(ego$d_senior70) <- c("No", "Yes")
levels(ego$d_senior75) <- c("No", "Yes")
levels(ego$d_senior80) <- c("No", "Yes")
levels(ego$d_hh_size) <- c("No cohabitants", "Has cohabitants")
levels(ego$chronic_comorb_self2) <- c("No", "Yes")
levels(ego$chronic_comorb_cohab2) <- c("No", "Yes")

levels(ego$covid) <- c("No", "Yes")
levels(ego$kindergarten_05) <- c("No", "Yes")
levels(ego$perceived_income) <- c("Medium","Very low", "Low", "High","Very high")
levels(ego$knows_severe)<- c("No", "Yes")
levels(ego$d_time_since_covid)<- c("No", "Yes")

##########################
# Variable labels to OUT #
##########################
var.labels = c(region_grouped_IT = "Macro-region of residence",
               respondent_gender = "Gender",
               education = "Highest education level achieved",
               occupation_agg = "Occupation",
               d_occupation = "Occupation",
               occupation_agg2 = "Occupation",
               sunday = "Contacts took place on Sunday",
               weekend = "Contacts took place over the weekend",
               worked_in_presence = "Worked in presence yesterday",
               attend_in_presence = "Went to work or schoold in presence yesterday",
               d_children = "Has minors in the household",
               d_senior65 = "Has people older than 65 in the household",
               d_senior70 = "Has people older than 70 in the household",
               d_senior75 = "Has people older than 75 in the household",
               d_senior80 = "Has people older than 80 in the household",
               d_hh_size = "Has at least a cohabitant",
               hh_size = "Household size (respondent included)",
               hh_size_det = "Household size (respondent included)",
               chronic_comorb_self2 = "Has at least a comorbidity",
               chronic_comorb_cohab2 = "At least one cohabitant has a comorbidity",
               covid = "Has been officially diagnosed with COVID-19 at least once",
               time_since_covid = "Time since last COVID infection",
               d_time_since_covid = "Had COVID-19 in the last 4 months",
               vacc_covid_bin = "Number of COVID-19 vaccine doses received",
               d_vacc2 = "Has received at least one dose of COVID-19 vaccine",
               vacc_covid_bin_2 = "Number of COVID-19 vaccine doses received",
               other_vacc = "Vaccinated against HPV, Meningo, Pneumo, HepA or Herpes zoster",
               refusal = "Refused vaccination at least once upon offer",
               perceived_income = "Perceived family wealth",
               income_1000 = "Family income",
               income = "Family income",
               income_threecat = "Family income",
               knows_severe = "Knows someone who was hospitalized/died due to COVID-19",
               age_group_10 = "Age group",
               age_group_25 = "Age group",
               wave = "Survey wave",
               EPID = "ID",
               
               respondent_age = "Age",
               kindergarten_05 = "Child goes to nursery",
               presence_school = "Attended school in presence",
               total_contacts = "Total contacts",
               total_contacts_nohh = "Total contacts (no household contacts)",
               indoor_nocohab = "Total indoor contacts (no household contacts)",
               indoor_mask_nocohab = "Total masked indoor contacts (no household contacts)",
               indoor_nomask_nocohab = "Total unmasked indoor contacts (no household contacts)",
               indoor_nocohab_notransp = "Total contacts indoor excl transport (no household contacts)",
               indoor_mask_nocohab_notransp = "Total masked indoor contacts excl transport (no household contacts)",
               outdoor_nocohab =  "Total outdoor contacts (no household contacts)",
               outdoor_mask_nocohab = "Total masked outdoor contacts (no household contacts)",
               outdoor_nomask_nocohab = "Total unmasked outdoor contacts (no household contacts)",
               risky = "Risky contacts (indoor, unmasked, no household)",
               total_contacts_prol = "Total contacts - prolonged",
               total_contacts_prol_soft = "Total contacts - prolonged soft",
               total_contacts_nohh_prol = "Total contacts - prolonged (no household contacts)",
               c_sharedindoor = "Total contacts - indirect", #added 29/7
               age = "Age (adult)"
)

data_mod <- ego[,colnames(ego) %in% names(var.labels)]
gg_miss_var(data_mod[names(var.labels)], show_pct = TRUE)

Hmisc::label(data_mod) =  as.list(var.labels[match(names(data_mod), names(var.labels))])


##############################################
##############################################
# SAVING STAT. MOD. CONTACTS AND Respondents #
##############################################
##############################################
# keep complete cases based on REGRESSORS that must be defined on everyone - so excluding nursery and all the contacts excl hh

var_reg = names(var.labels[1:37])
data_mod1 <-  dplyr::select(data_mod, var_reg)
data_mod1=data_mod1[complete.cases(data_mod1),] #5028 

data_mod2 = left_join(data_mod1,data_mod)

# Final datasets ----

# children
data_mod_children <- data_mod2[data_mod2$respondent_age <= 17, ]
data_mod_children <- arrange(data_mod_children,EPID,wave)

# adults 
data_mod <- data_mod2[data_mod2$respondent_age > 17, ]
data_mod <- arrange(data_mod,EPID,wave)


save(list=c("data_mod","data_mod_children","var_reg","var.labels"),file="data_in/datamod.RData")
write.csv(data_mod, "data_in/datamod.csv", row.names=FALSE)
write.csv(data_mod_children, "data_in/datamod_children.csv", row.names=FALSE)
