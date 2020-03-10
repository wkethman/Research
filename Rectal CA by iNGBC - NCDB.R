library(psych)
library(Rcmdr)
library(survey)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(twang)
library(haven)
library(survival)
library(dplyr)
library(survminer)
library(coxphw)
library(foreign)
library(jskm)

  dir <- "./Desktop/WK Desktop/Rectal CA - Appropriate  - NCDB/"
  NCDB <- as.data.frame(read_sav(paste(dir,"NCDB_PUF_RECTUM_Nov-03-2019 - Data/NCDBPUF_Rectum.0.2016.0.sav",sep = "")))

  #Clean data
#Variable reviewed - AGE, SEX, FACILITY_TYPE_CD, MEDICAID_EXPN_CODE, CDCC_TOTAL_BEST, PRIMARY_SITE, TNM_EDITION_NUMBER, SURG_DISCHARGE_DAYS, REASON_FOR_NO_SURGERY, RAD_LOCATION_OF_RX
#Need to review - REGIONAL_NODES_POSITIVE, REGIONAL_NODES_EXAMINED
#Notes:
  #Decided not to filter LYMPH_VASCULAR_INVASION because it wouldn't effect 30-/90-day outcomes and because 21,235 are unknown and would be filtered out
  #Call and figure out why 1,347 might have gotten 900 cGY boost dose
  
  NCDB <- subset(NCDB, YEAR_OF_DIAGNOSIS >= 2010 & YEAR_OF_DIAGNOSIS <= 2015) #Include only cases 2010-2015
  NCDB <- subset(NCDB, PUF_MULT_SOURCE == 0) #Include only those who were reported by one facility to NCDB
  
  #Create new histology type variable to better differentiate types
  attach(NCDB)
  NCDB$HISTOLOGY_TYPE[HISTOLOGY == 8140 | HISTOLOGY == 8144 | HISTOLOGY == 8201 | HISTOLOGY == 8210 | HISTOLOGY == 8211 | HISTOLOGY == 8213 | HISTOLOGY == 8220 | HISTOLOGY == 8261 | HISTOLOGY == 8262 | HISTOLOGY == 8263 | HISTOLOGY == 8490] <- 1 #Adenocarcinoma
  NCDB$HISTOLOGY_TYPE[HISTOLOGY == 8480 | HISTOLOGY == 8481] <- 2 #Mucinous adenocarcinoma
  NCDB$HISTOLOGY_TYPE[HISTOLOGY == 8936] <- 3 #GIST
  NCDB$HISTOLOGY_TYPE[HISTOLOGY == 8240 | HISTOLOGY == 8246 | HISTOLOGY == 8249] <- 4 #Neuroendocrine carcinoma
  detach(NCDB)
  
  NCDB$HISTOLOGY_TYPE <- replace(NCDB$HISTOLOGY_TYPE,is.na(NCDB$HISTOLOGY_TYPE), 5) 
  
  NCDB <- subset(NCDB, HISTOLOGY_TYPE == 1) #Filter to only Adenocarcinoma
  NCDB <- subset(NCDB, BEHAVIOR == 3) #Include only invasive disease
  NCDB <- subset(NCDB, RX_SUMM_TREATMENT_STATUS == 1) #Include only those who received treatment
  NCDB <- subset(NCDB, RX_HOSP_SURG_APPR_2010 != 9 & RX_HOSP_SURG_APPR_2010 != 0) #Exclude those who did not have surgery or if unknown
  
  #Filter TNM clinical staging
  NCDB$TNM_CLIN_T <- gsub("[[:space:]]", "", NCDB$TNM_CLIN_T) #Remove white spaces
  NCDB$TNM_CLIN_N <- gsub("[[:space:]]", "", NCDB$TNM_CLIN_N) #Remove white spaces
  NCDB$TNM_CLIN_M <- gsub("[[:space:]]", "", NCDB$TNM_CLIN_M) #Remove white spaces
  NCDB$TNM_CLIN_STAGE_GROUP <- gsub("[[:space:]]", "", NCDB$TNM_CLIN_STAGE_GROUP) #Remove white spaces
  
  NCDB <- subset(NCDB, TNM_CLIN_T == "c1" | TNM_CLIN_T == "c2" | TNM_CLIN_T == "c3" | TNM_CLIN_T == "c4" | TNM_CLIN_T == "c1A" | TNM_CLIN_T == "c1B" | TNM_CLIN_T == "c4A" | TNM_CLIN_T == "c4B") #Clinical T stage, excluding blank, not appplicable, cX, c0, pIS
  NCDB <- subset(NCDB, TNM_CLIN_N == "c0" | TNM_CLIN_N == "c1" | TNM_CLIN_N == "c2" | TNM_CLIN_N == "c1A" | TNM_CLIN_N == "c1B" | TNM_CLIN_N == "c1C" | TNM_CLIN_N == "c2A" | TNM_CLIN_N == "c2B") #Clinical N stage, excluding blank, not applicable, cX
  NCDB <- subset(NCDB, TNM_CLIN_M == "c0" | TNM_CLIN_M == "c1" | TNM_CLIN_M == "c1A" | TNM_CLIN_M == "c1B") #Clinical M stage, excluding blank
  NCDB <- subset(NCDB, TNM_CLIN_STAGE_GROUP != "88" & TNM_CLIN_STAGE_GROUP != "99" & TNM_CLIN_STAGE_GROUP != "")
  
  #Re-categorize T-stage
  attach(NCDB)
  NCDB$TNM_CLIN_T_COMB[TNM_CLIN_T == "c1" | TNM_CLIN_T == "c1A" | TNM_CLIN_T == "c1B"] <- 1 #T1
  NCDB$TNM_CLIN_T_COMB[TNM_CLIN_T == "c2"] <- 2 #T2
  NCDB$TNM_CLIN_T_COMB[TNM_CLIN_T == "c3"] <- 3 #T3
  NCDB$TNM_CLIN_T_COMB[TNM_CLIN_T == "c4" | TNM_CLIN_T == "c4A" | TNM_CLIN_T == "c4B"] <- 4 #T4
  detach(NCDB)
  
  #Re-categorize N-stage    
  attach(NCDB)
  NCDB$TNM_CLIN_N_COMB[TNM_CLIN_N == "c0"] <- 0 #N0
  NCDB$TNM_CLIN_N_COMB[TNM_CLIN_N == "c1" | TNM_CLIN_N == "c1A" | TNM_CLIN_N == "c1B" | TNM_CLIN_N == "c1C"] <- 1 #N1
  NCDB$TNM_CLIN_N_COMB[TNM_CLIN_N == "c2" | TNM_CLIN_N == "c2A" | TNM_CLIN_N == "c2B"] <- 2 #N2
  detach(NCDB)
  
  #Re-categorize M-stage  
  attach(NCDB)
  NCDB$TNM_CLIN_M_COMB[TNM_CLIN_M == "c0"] <- 0 #M0
  NCDB$TNM_CLIN_M_COMB[TNM_CLIN_M == "c1" | TNM_CLIN_M == "c1A" | TNM_CLIN_M == "c1B"] <- 1 #M1
  detach(NCDB)
  
  #T1-4, N1-2 clinical pre-operative disease, exclude metastatic disease
  NCDB <- subset(NCDB, TNM_CLIN_M_COMB == 0)

  #Cleaning data
  NCDB <- subset(NCDB, FACILITY_TYPE_CD <= 4) #Exclude other program types
  NCDB <- subset(NCDB, INSURANCE_STATUS < 9) #Exclude unknown insurance status
  NCDB <- subset(NCDB, RACE < 99) #Exclude unknown
  NCDB <- subset(NCDB, SPANISH_HISPANIC_ORIGIN < 9) #Exclude unknown
  
  #Create new race/hispanic origin variable
  attach(NCDB)
  NCDB$ETHNICITY[RACE == 1 & SPANISH_HISPANIC_ORIGIN == 0] <- 1 #Non-hispanic white
  NCDB$ETHNICITY[RACE == 2 & SPANISH_HISPANIC_ORIGIN == 0] <- 2 #Non-hispanic black
  NCDB$ETHNICITY[SPANISH_HISPANIC_ORIGIN == 0 & RACE >= 4 & RACE <= 97] <- 4 #Asian/Pacific islander
  NCDB$ETHNICITY[SPANISH_HISPANIC_ORIGIN >= 1] <- 3 #Hispanic
  NCDB$ETHNICITY[RACE == 3 | RACE == 98 & SPANISH_HISPANIC_ORIGIN == 0] <- 5 #Other
  detach(NCDB)
  
  #Re-categorize unplanned readmissions
  attach(NCDB)
  NCDB$UNPLAN_READMIT[READM_HOSP_30_DAYS == 0 | READM_HOSP_30_DAYS == 2] <- 0 #No unplanned readmissions
  NCDB$UNPLAN_READMIT[READM_HOSP_30_DAYS == 1 | READM_HOSP_30_DAYS == 3] <- 1 #Unplanned readmission occurred
  detach(NCDB)

  NCDB <- subset(NCDB, BEHAVIOR == 3) #Include only invasive disease
  NCDB <- subset(NCDB, GRADE != 9) #Differentiation, filter those that are unknown
  NCDB <- subset(NCDB, RX_SUMM_SURGICAL_MARGINS <= 3) #Exclude unknown or indeterminate margin status
  NCDB <- subset(NCDB, READM_HOSP_30_DAYS != 9) #Readmission unknown
  NCDB <- subset(NCDB, RX_SUMM_RADIATION <= 1) #Exclude patient whose radiation course was unknown or whom had atypical radiation administration (other than external beam) - All facilities
  NCDB <- subset(NCDB, RAD_TREAT_VOL == 0 | RAD_TREAT_VOL == 29) #0 - No radiation, 29 - pelvis
  NCDB <- subset(NCDB, RAD_REGIONAL_DOSE_CGY != 99999) #Exclude unknown regional dose
  NCDB <- subset(NCDB, RAD_BOOST_DOSE_CGY != 99999) #Exclude unknown regional dose
  NCDB <- subset(NCDB, RX_SUMM_SURGRAD_SEQ < 5) #Filter those who received either no radiation or neoadjuvant radiation, or intra-operative radiation
  NCDB <- subset(NCDB, REASON_FOR_NO_RADIATION < 8) #Exclude patients for which we don't know if radiation was given
  NCDB <- subset(NCDB, RX_SUMM_CHEMO < 88) #Exclude patients for which we don't know if chemotherapy was given
  NCDB <- subset(NCDB, RX_SUMM_SYSTEMIC_SUR_SEQ < 5) #Exclude systemic therapy where sequence was unknown or where given intra-operatively, this includes pre/post-operative chemo
  NCDB <- subset(NCDB, REGIONAL_NODES_EXAMINED >= 1 & REGIONAL_NODES_EXAMINED <= 90)
  NCDB <- subset(NCDB, PUF_30_DAY_MORT_CD != 9) #Exclude unknown mortality
  NCDB <- subset(NCDB, PUF_90_DAY_MORT_CD != 9) #Exclude unknown mortality
  NCDB <- subset(NCDB, !is.na(NCDB$SURG_DISCHARGE_DAYS)) 
  NCDB <- subset(NCDB, !is.na(NCDB$NO_HSD_QUAR_16)) 
  NCDB <- subset(NCDB, !is.na(NCDB$UR_CD_13)) 
  NCDB <- subset(NCDB, !is.na(NCDB$MED_INC_QUAR_16))
  
  NCDB$COMBINED_RAD_DOSE <- NCDB$RAD_REGIONAL_DOSE_CGY + NCDB$RAD_BOOST_DOSE_CGY 
  
  #Ideal radiation dose administered (irrespective of timing)
  attach(NCDB)
   NCDB$IDEAL_RAD_DOSE[COMBINED_RAD_DOSE == 0] <- 0 #No radiation received
   NCDB$IDEAL_RAD_DOSE[((COMBINED_RAD_DOSE >= 4500 & COMBINED_RAD_DOSE <= 5540) | COMBINED_RAD_DOSE == 2500)] <- 1 #Ideal dose received
  detach(NCDB)
   
  NCDB$IDEAL_RAD_DOSE <- replace(NCDB$IDEAL_RAD_DOSE, is.na(NCDB$IDEAL_RAD_DOSE), 2) #Not ideal dose
   
  #Ideal neoadjuvant radiation therapy administered
  attach(NCDB)
    NCDB$RAD_RECEIVED[IDEAL_RAD_DOSE == 0 | RX_SUMM_SURGRAD_SEQ == 3] <- 0 #No neoadjuvant radiation received
    NCDB$RAD_RECEIVED[((RX_SUMM_SURGRAD_SEQ == 2 | RX_SUMM_SURGRAD_SEQ == 4) & IDEAL_RAD_DOSE == 1)] <- 1 #Ideal standard or short-course dose given in neoadjuvant setting
  detach(NCDB)
  
  NCDB$RAD_RECEIVED <- replace(NCDB$RAD_RECEIVED, is.na(NCDB$RAD_RECEIVED), 2) #Received neoadjuvant but NOT ideal dose or timing
    
  #Re-categorize whether neoadjuvant chemotherapy was given (this assigns patients who only received adjuvant chemotherapy to the no chemotherapy group)
  attach(NCDB)
    NCDB$CHEMO_RECEIVED[RX_SUMM_SYSTEMIC_SUR_SEQ == 0 | RX_SUMM_SYSTEMIC_SUR_SEQ == 3] <- 0 #Neodjuvant chemotherapy NOT given
    NCDB$CHEMO_RECEIVED[RX_SUMM_SYSTEMIC_SUR_SEQ == 2 | RX_SUMM_SYSTEMIC_SUR_SEQ == 4] <- 1 #Neoadjuvant chemotherapy given
  detach(NCDB)
    
  #Categorize appropriate chemoradiation together
  attach(NCDB)
    NCDB$CHEMO_RAD_IDEAL[NCDB$RAD_RECEIVED == 0 & NCDB$CHEMO_RECEIVED == 0] <- 0 #No neoadjuvant chemoradiation given
    NCDB$CHEMO_RAD_IDEAL[NCDB$RAD_RECEIVED == 1 & NCDB$CHEMO_RECEIVED == 1] <- 1 #Appropriate neoadjuvant chemoradiation given
  detach(NCDB)
    
  NCDB$CHEMO_RAD_IDEAL <- replace(NCDB$CHEMO_RAD_IDEAL, is.na(NCDB$CHEMO_RAD_IDEAL), 2) #Received neoadjuvant chemotherapy but NOT ideal dose
  
  #Ideal neoadjuvant - 0, less than ideal neoadjuvant - 1  
  NCDB$IDEAL_TX <- ifelse(((NCDB$TNM_CLIN_T_COMB < 3) & (NCDB$TNM_CLIN_N_COMB == 0)),
                                ifelse(((NCDB$RAD_RECEIVED == 0) & (NCDB$CHEMO_RECEIVED == 0)), 0, 1),
                                ifelse(((NCDB$RAD_RECEIVED == 1) & (NCDB$CHEMO_RECEIVED == 1)), 0, 1))

  #Categorize margin status
  attach(NCDB)
  NCDB$MARGIN_STATUS[NCDB$RX_SUMM_SURGICAL_MARGINS == 0] <- 1 #Negative margin
  NCDB$MARGIN_STATUS[NCDB$RX_SUMM_SURGICAL_MARGINS > 0] <- 0 #Positive margin
  detach(NCDB)
  
  #Categorize nodal harvest status
  attach(NCDB)
  NCDB$NODE_STATUS[NCDB$REGIONAL_NODES_EXAMINED < 12] <- 0 #Poor nodal harvist
  NCDB$NODE_STATUS[NCDB$REGIONAL_NODES_EXAMINED >= 12] <- 1 #Good nodal harvist
  detach(NCDB)
  
  #Categorized oncologic status
  attach(NCDB)
  NCDB$ONC_SUCCESS[MARGIN_STATUS == 1 & REGIONAL_NODES_EXAMINED >= 12] <- 1 #Negative margins and nodes >=12 examined, success
  NCDB$ONC_SUCCESS[MARGIN_STATUS == 0 | REGIONAL_NODES_EXAMINED < 12] <- 0 #Positive margins and nodes <12 examined, unsuccessful
  detach(NCDB)
  
  #Categorize age  <50, 50-65, >65
  attach(NCDB)
  NCDB$AGE_CAT_50_65[NCDB$AGE < 50] <- 1 #Age <50
  NCDB$AGE_CAT_50_65[NCDB$AGE >= 50 & NCDB$AGE < 65] <- 2 #Age 50-65
  NCDB$AGE_CAT_50_65[NCDB$AGE >= 65] <- 3 #Age >65
  detach(NCDB)
  
  #Categorize age >65
  attach(NCDB)
  NCDB$AGE_CAT_65[NCDB$AGE < 65] <- 0 #Age <65
  NCDB$AGE_CAT_65[NCDB$AGE >= 65] <- 1 #Age >65
  detach(NCDB)
  
  attach(NCDB)
  NCDB$RURAL_METRO[NCDB$UR_CD_13 == 1 | NCDB$UR_CD_13 == 2 | NCDB$UR_CD_13 == 3] <- 1 #Metro
  NCDB$RURAL_METRO[NCDB$UR_CD_13 == 4 | NCDB$UR_CD_13 == 5 | NCDB$UR_CD_13 == 6 | NCDB$UR_CD_13 == 7] <- 2 #Urban
  NCDB$RURAL_METRO[NCDB$UR_CD_13 == 8 | NCDB$UR_CD_13 == 9] <- 3 #Rural
  detach(NCDB)
  
  attach(NCDB)
  NCDB$EARLY_LATE_STAGE[(NCDB$TNM_CLIN_T_COMB == 1 | NCDB$TNM_CLIN_T_COMB == 2) & NCDB$TNM_CLIN_N_COMB == 0] <- 1 #Early
  detach(NCDB)
  
  NCDB$EARLY_LATE_STAGE <- replace(NCDB$EARLY_LATE_STAGE, is.na(NCDB$EARLY_LATE_STAGE), 2) #Late
  
  attach(NCDB)
  NCDB$STAGE_IDEAL[NCDB$EARLY_LATE_STAGE == 1 & NCDB$IDEAL_TX == 0] <- 1 #Early stage ideal
  NCDB$STAGE_IDEAL[NCDB$EARLY_LATE_STAGE == 1 & NCDB$IDEAL_TX == 1] <- 2 #Early stage less than ideal
  NCDB$STAGE_IDEAL[NCDB$EARLY_LATE_STAGE == 2 & NCDB$IDEAL_TX == 0] <- 3 #Late stage ideal
  NCDB$STAGE_IDEAL[NCDB$EARLY_LATE_STAGE == 2 & NCDB$IDEAL_TX == 1] <- 4 #Late stage less than ideal
  detach(NCDB)
  
#Filter pathogic metastatic disease
  NCDB <- subset(NCDB, NCDB$TNM_PATH_M == "")
  
#STOP HERE FOR EITHER ANALYSIS

#Evaluate variables to see if they make sense
table(NCDB$RX_SUMM_TREATMENT_STATUS, useNA = 'ifany')
table(NCDB$PRIMARY_SITE, useNA = 'ifany')
table(NCDB$FACILITY_TYPE_CD, useNA = 'ifany')
table(NCDB$FACILITY_LOCATION_CD, useNA = 'ifany')
table(NCDB$NO_HSD_QUAR_16, useNA = 'ifany')
table(NCDB$UR_CD_13, useNA = 'ifany')
table(NCDB$MED_INC_QUAR_16, useNA = 'ifany')
table(NCDB$HISTOLOGY_TYPE, useNA = 'ifany') #Adenocarcinoma only
table(NCDB$RX_SUMM_SURG_PRIM_SITE, useNA = 'ifany') #Any facility
table(NCDB$RX_HOSP_SURG_PRIM_SITE, useNA = 'ifany') #Primary facility
table(NCDB$RX_HOSP_SURG_APPR_2010, useNA = 'ifany') #Surgical approach
table(NCDB$RX_SUMM_SURGICAL_MARGINS, useNA = 'ifany') #Surgical margins
table(NCDB$MARGIN_STATUS, useNA = 'ifany') #Surgical margins
table(NCDB$EARLY_LATE_STAGE, NCDB$TNM_CLIN_T_COMB)
table(NCDB$EARLY_LATE_STAGE, NCDB$TNM_CLIN_N_COMB)

#Propensity co-variates
summary(NCDB$AGE)
table(NCDB$AGE_CAT_50_65, useNA = 'ifany') #Surgical margins
table(NCDB$AGE_CAT_65, useNA = 'ifany') #Surgical margins
table(NCDB$CDCC_TOTAL_BEST, useNA = 'ifany')
table(NCDB$TNM_CLIN_T_COMB, useNA = 'ifany')
table(NCDB$TNM_CLIN_N_COMB, useNA = 'ifany')
table(NCDB$TNM_CLIN_M_COMB, useNA = 'ifany')
table(NCDB$GRADE, useNA = 'ifany')
table(NCDB$ETHNICITY, useNA = 'ifany')

#Grouping
table(NCDB$IDEAL_TX, useNA = 'ifany') #Appropriate treatment

#Outcomes
table(NCDB$MARGIN_STATUS, useNA = 'ifany')
summary(NCDB$REGIONAL_NODES_EXAMINED)
table(NCDB$ONC_SUCCESS, useNA = 'ifany')

table(NCDB$IDEAL_TX, NCDB$MARGIN_STATUS)
table(NCDB$IDEAL_TX, NCDB$NODE_STATUS)
table(NCDB$IDEAL_TX, NCDB$ONC_SUCCESS)
table(NCDB$IDEAL_TX, NCDB$PUF_90_DAY_MORT_CD)
table(NCDB$IDEAL_TX, NCDB$EARLY_LATE_STAGE)
table(NCDB$STAGE_IDEAL)

summary(NCDB$SURG_DISCHARGE_DAYS)
table(NCDB$UNPLAN_READMIT, useNA = 'ifany')
table(NCDB$PUF_30_DAY_MORT_CD, useNA = 'ifany')
table(NCDB$PUF_90_DAY_MORT_CD, useNA = 'ifany')

#Propensity for ideal treatment
var_names <- c("STAGE_IDEAL", "PUF_30_DAY_MORT_CD", "PUF_90_DAY_MORT_CD", "ONC_SUCCESS", "TNM_CLIN_T_COMB", "TNM_CLIN_N_COMB", "MARGIN_STATUS", "UNPLAN_READMIT", "SEX", "GRADE", "ETHNICITY", "CDCC_TOTAL_BEST", "FACILITY_TYPE_CD", "FACILITY_LOCATION_CD", "NO_HSD_QUAR_16", "INSURANCE_STATUS", "UR_CD_13", "MED_INC_QUAR_16", "AGE_CAT_50_65", "AGE_CAT_65")
NCDB[,var_names] <- lapply(NCDB[,var_names], factor)

NCDB$DX_LASTCONTACT_DEATH_MONTHS <- as.numeric(NCDB$DX_LASTCONTACT_DEATH_MONTHS)
NCDB$AGE <- as.numeric(NCDB$AGE)

#Survival based on ideal therapy or not
set.seed(1)
ps.NCDB <- ps(IDEAL_TX~AGE+SEX+FACILITY_TYPE_CD+CDCC_TOTAL_BEST+TNM_CLIN_T_COMB+TNM_CLIN_N_COMB+GRADE,
                  data=as.data.frame(NCDB),
                  n.trees = 3000,
                  stop.method = c("ks.mean"),
                  estimand = "ATE",
                  verbose = FALSE)
#plot(ps.NCDB, plots = 1)
#plot(ps.NCDB, plots = 2)
#plot(ps.NCDB, plots = 3)
#bal.table(ps.NCDB, digits = 3)
today <- Sys.Date()
file_name <- paste("Manuscript/Balance Table (Ideal) - ", format(today, format="%m%d%Y"), ".csv", sep = "")
write.csv(bal.table(ps.NCDB, digits = 3), file = paste(dir, file_name,sep = ""))
summary(ps.NCDB)

NCDB$weights <- get.weights(ps.NCDB, stop.method = "ks.mean")
design.ps <- svydesign(ids=~1, weights=~weights, data=NCDB)

#Outcomes
#Margin status
margin <- svyglm(MARGIN_STATUS ~ IDEAL_TX, design = design.ps, family=quasibinomial)
summary(margin)
exp(cbind(OR = coef(margin), confint(margin)))
#Nodes status
node <- svyglm(NODE_STATUS ~ IDEAL_TX, design = design.ps, family=quasibinomial)
summary(node)
exp(cbind(OR = coef(node), confint(node)))
#Lymph nodes examined
nodes_exam <- svyglm(REGIONAL_NODES_EXAMINED ~ IDEAL_TX, design = design.ps, family=gaussian())
summary(nodes_exam)
exp(cbind(OR = coef(nodes_exam), confint(nodes_exam)))
#Oncologic success
onc_success <- svyglm(ONC_SUCCESS ~ IDEAL_TX, design = design.ps, family=quasibinomial)
summary(onc_success)
exp(cbind(OR = coef(onc_success), confint(onc_success)))
#30-day mortality
#mort_30 <- svyglm(PUF_30_DAY_MORT_CD ~ IDEAL_TX, design = design.ps, family=quasibinomial)
#summary(mort_30)
#exp(cbind(OR = coef(mort_30), confint(mort_30)))
#90-day mortality
mort_90 <- svyglm(PUF_90_DAY_MORT_CD ~ IDEAL_TX, design = design.ps, family=quasibinomial)
summary(mort_90)
exp(cbind(OR = coef(mort_90), confint(mort_90)))

#Survival by ideal neoadjuvant treatment (weighted)
km.by.ideal <- svykm(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS == 0) ~ IDEAL_TX, design = design.ps, se = TRUE)
text.title <- element_text(size = 12, margin = margin(t = 0, r = 0, b = 20, l = 0, unit = "pt"), hjust = 0.5)
text.axis.x <- element_text(size = 10, margin = margin(t = 10, r = 0, b = 5, l = 0, unit = "pt"))
text.axis.y <- element_text(size = 10, margin = margin(t = 10, r = 10, b = 10, l = 5, unit = "pt"))
p.km <- svyjskm(km.by.ideal, ci=TRUE, ylims=c(0.5,1.0), xlims=c(0, 60), xlabs = "Time-to-event (months)", ylabs = "Survival probability", linecols="black", legendposition = c(0.3, 0.3), main="Survival analysis by ideal NCCN guideline-based care (weighted)", ystrataname = "", ystratalabs = c("iNGBC","non-iNGBC"), dashed=TRUE)
p.km + theme_bw() +
  theme(panel.grid.minor = element_blank(), plot.title = text.title, axis.text.y = text.axis.y, axis.text.x = text.axis.x, plot.margin = margin(t = 30, r = 20, b = 20, l = 20, unit = "pt"))
ggsave(filename = paste(dir,"Manuscript/Tables - Figures/Figure - 2.pdf",sep = ""), dpi = 300, width = 8, height = 6)

#Plot odds ratios
text.title <- element_text(size = 12, margin = margin(t = 0, r = 0, b = 20, l = 0, unit = "pt"), hjust = 0.5)
text.axis.x <- element_text(size = 10, margin = margin(t = 10, r = 0, b = 5, l = 0, unit = "pt"))
text.axis.y <- element_text(size = 10, margin = margin(t = 10, r = 10, b = 10, l = 0, unit = "pt"))
boxLabels = c("Negative margin", "Adequate nodal harvest", "Oncologic success", "90-day mortality")
yAxis = length(boxLabels):1
df <- data.frame(yAxis = length(boxLabels):1,
                 boxOdds = c(as.numeric(exp(coef(margin))[2]), as.numeric(exp(coef(node))[2]), as.numeric(exp(coef(onc_success))[2]), as.numeric(exp(coef(mort_90))[2])),
                 boxCILow = c(as.numeric(exp(confint(margin))[2,1]), as.numeric(exp(confint(node))[2,1]), as.numeric(exp(confint(onc_success))[2,1]), as.numeric(exp(confint(mort_90))[2,1])),
                 boxCIHigh = c(as.numeric(exp(confint(margin))[2,2]), as.numeric(exp(confint(node))[2,2]), as.numeric(exp(confint(onc_success))[2,2]), as.numeric(exp(confint(mort_90))[2,2])))
p.ideal <- ggplot(df, aes(x = boxOdds, y = yAxis))
p.ideal + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "black") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), plot.title = text.title, axis.text.y = text.axis.y, axis.text.x = text.axis.x, plot.margin = margin(t = 30, r = 20, b = 20, l = 20, unit = "pt")) +
  scale_y_continuous(breaks = yAxis, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,2,0.5)) +
  coord_trans(x = "log10") +
  labs(title = "Ideal NCCN guideline-based care", y = "", x = "Odds ratio (log scale)") +
  annotate("text", x = 0.73, y = 0.4, label = "Favors iNGBC") +
  annotate("text", x = 1.55, y = 0.4, label = "Favors non-iNGBC")
ggsave(filename = paste(dir,"Manuscript/Tables - Figures/Figure - 1.pdf",sep = ""), dpi = 300, width = 8, height = length(boxLabels))

#Survival for multiple categories
set.seed(1)
mnps.NCDB <- mnps(STAGE_IDEAL~AGE+SEX+FACILITY_TYPE_CD+CDCC_TOTAL_BEST+GRADE,
                  data=as.data.frame(NCDB),
                  n.trees = 3000,
                  stop.method = c("ks.mean"),
                  estimand = "ATE",
                  verbose = FALSE)
#plot(mnps.NCDB, plots = 1)
#plot(mnps.NCDB, plots = 2)
#plot(mnps.NCDB, plots = 3)
#bal.table(mnps.NCDB, digits = 3)
today <- Sys.Date()
file_name <- paste("Manuscript/Balance Table (Stage-Ideal) - ", format(today, format="%m%d%Y"), ".csv", sep = "")
write.csv(bal.table(mnps.NCDB, digits = 3), file = paste(dir, file_name,sep = ""))
summary(mnps.NCDB)

NCDB$weights_1 <- get.weights(mnps.NCDB, stop.method = "ks.mean")
design.mnps <- svydesign(ids=~1, weights=~weights_1, data=NCDB)

km.by.stage_ideal <- svykm(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS == 0) ~ STAGE_IDEAL, design = design.mnps, se = TRUE)
text.title <- element_text(size = 12, margin = margin(t = 0, r = 0, b = 20, l = 0, unit = "pt"), hjust = 0.5)
text.axis.x <- element_text(size = 10, margin = margin(t = 10, r = 0, b = 5, l = 0, unit = "pt"))
text.axis.y <- element_text(size = 10, margin = margin(t = 10, r = 10, b = 10, l = 5, unit = "pt"))
p.km <- svyjskm(km.by.stage_ideal, ci=TRUE, ylims=c(0.5,1.0), xlims=c(0, 60), xlabs = "Time-to-event (months)", ylabs = "Survival probability", linecols="black", legendposition = c(0.3, 0.3), main="Survival analysis by stage and ideal NCCN guideline-based care (weighted)", ystrataname = "", ystratalabs = c("Early stage - iNGBC","Early stage - non-iNGBC", "Later stage - iNGBC", "Later stage - non-iNGBC"), dashed=TRUE)
p.km + theme_bw() +
  theme(panel.grid.minor = element_blank(), plot.title = text.title, axis.text.y = text.axis.y, axis.text.x = text.axis.x, plot.margin = margin(t = 30, r = 20, b = 20, l = 20, unit = "pt"))
ggsave(filename = paste(dir,"Manuscript/Tables - Figures/Figure - 3.pdf",sep = ""), dpi = 300, width = 8, height = 6)

#Relevel before logistic regression
NCDB$SEX <- factor(NCDB$SEX, levels = c(1,2), labels = c("Male","Female"))
NCDB$SEX <- relevel(NCDB$SEX, ref = "Female")

NCDB$AGE_CAT_50_65 <- factor(NCDB$AGE_CAT_50_65, levels = c(1,2,3), labels = c("<50","50-65",">65"))
NCDB$AGE_CAT_50_65 <- relevel(NCDB$AGE_CAT_50_65, ref = "50-65") 

NCDB$ETHNICITY <- factor(NCDB$ETHNICITY, levels = c(1,2,3,4,5), labels = c("Non-hispanic white", "Non-hispanic black", "Hispanic", "Asian/Pacific islander", "Other"))
NCDB$ETHNICITY <- relevel(NCDB$ETHNICITY, ref = "Non-hispanic white") 

NCDB$FACILITY_LOCATION_CD <- factor(NCDB$FACILITY_LOCATION_CD, levels = c(1,2,3,4,5,6,7,8,9), labels = c("New England","Middle Atlantic","South Atlantic","East North Central","East South Central","West North Centra","West South Central","Mountain","Pacific"))
NCDB$FACILITY_LOCATION_CD <- relevel(NCDB$FACILITY_LOCATION_CD, ref = "New England") #Reference - New England

NCDB$FACILITY_TYPE_CD <- factor(NCDB$FACILITY_TYPE_CD, levels = c(1,2,3,4), labels = c("Community","Comprehensive community","Academic","Integrated network"))
NCDB$FACILITY_TYPE_CD <- relevel(NCDB$FACILITY_TYPE_CD, ref = "Academic") #Reference - Academic

NCDB$NO_HSD_QUAR_16 <- factor(NCDB$NO_HSD_QUAR_16, levels = c(1,2,3,4), labels = c(">17.6%","10.9-17.5%","6.3-10.8%","<6.3%"))
NCDB$NO_HSD_QUAR_16 <- relevel(NCDB$NO_HSD_QUAR_16, ref = "<6.3%") #Reference - <6.3% of population did not graduate from high school

NCDB$INSURANCE_STATUS <- factor(NCDB$INSURANCE_STATUS, levels = c(0,1,2,3,4), labels = c("Not insured","Private insurance","Medicaid","Medicare","Other government"))
NCDB$INSURANCE_STATUS <- relevel(NCDB$INSURANCE_STATUS, ref = "Private insurance") #Reference - Private insurance

NCDB$RURAL_METRO <- factor(NCDB$RURAL_METRO, levels = c(1,2,3), labels = c("Metro","Urban","Rural"))
NCDB$RURAL_METRO <- relevel(NCDB$RURAL_METRO, ref = "Metro") #Reference - Metro

NCDB$MED_INC_QUAR_16 <- factor(NCDB$MED_INC_QUAR_16, levels = c(1,2,3,4), labels = c("<40k","40-50K","50-63K",">63K"))
NCDB$MED_INC_QUAR_16 <- relevel(NCDB$MED_INC_QUAR_16, ref = ">63K") #Reference - >$63K

NCDB$EARLY_LATE_STAGE <- factor(NCDB$EARLY_LATE_STAGE, levels = c(1,2), labels = c("Early","Late"))
NCDB$EARLY_LATE_STAGE <- relevel(NCDB$EARLY_LATE_STAGE, ref = "Early") #Reference - Early

ideal <- glm(NCDB$IDEAL_TX ~ NCDB$SEX + NCDB$AGE_CAT_50_65 + NCDB$ETHNICITY + NCDB$FACILITY_LOCATION_CD + NCDB$FACILITY_TYPE_CD + NCDB$NO_HSD_QUAR_16 + NCDB$INSURANCE_STATUS + NCDB$RURAL_METRO + NCDB$EARLY_LATE_STAGE + NCDB$MED_INC_QUAR_16 + NCDB$MED_INC_QUAR_16 + NCDB$CDCC_TOTAL_BEST, family = binomial)
summary(ideal)
OR <- exp(cbind(OR = coef(ideal), confint(ideal)))