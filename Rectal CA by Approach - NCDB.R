# Effects of surgical approach in short and long-term outcomes in early stage rectal cancer: A Multicenter, Propensity Score-Weighted Cohort Study

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
library(ggplot2)
library(nnet)

  dir <- "./Desktop/WK Desktop/Rectal CA - Approach - NCDB/"
  NCDB <- as.data.frame(read_sav(paste(dir,"NCDB_PUF_RECTUM_Nov-03-2019 - Data/NCDBPUF_Rectum.0.2016.0.sav",sep = "")))

  #Clean data
#Variable reviewed - AGE, SEX, FACILITY_TYPE_CD, MEDICAID_EXPN_CODE, CDCC_TOTAL_BEST, PRIMARY_SITE, TNM_EDITION_NUMBER, SURG_DISCHARGE_DAYS, REASON_FOR_NO_SURGERY, RAD_LOCATION_OF_RX
#Need to review - REGIONAL_NODES_POSITIVE, REGIONAL_NODES_EXAMINED
#Notes:
  #Decided not to filter LYMPH_VASCULAR_INVASION because it wouldn't effect 30-/90-day outcomes and because 21,235 are unknown and would be filtered out
  
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
  
  NCDB <- subset(NCDB, (RX_SUMM_SURG_PRIM_SITE >= 30 & RX_SUMM_SURG_PRIM_SITE <= 60) | RX_SUMM_SURG_PRIM_SITE == 80) #Include partial proctectomy (30), coloanal (40), APR (50), total proct, NOS (60), proct, NOS (80) - exclude local and exent
  NCDB <- subset(NCDB, RX_SUMM_TREATMENT_STATUS == 1) #Include only those who received treatment
  NCDB <- subset(NCDB, RX_SUMM_SYSTEMIC_SUR_SEQ == 0) #No system therapy given
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
  
  #T1-2, N0 clinical pre-operatie disease
  NCDB <- subset(NCDB, TNM_CLIN_M_COMB == 0)
  NCDB <- subset(NCDB, TNM_CLIN_N_COMB == 0)
  NCDB <- subset(NCDB, TNM_CLIN_T_COMB == 1 | TNM_CLIN_T_COMB == 2)
  
  #Filter TNM pathologic staging
  NCDB$TNM_PATH_T <- gsub("[[:space:]]", "", NCDB$TNM_PATH_T) #Remove white spaces
  NCDB$TNM_PATH_N <- gsub("[[:space:]]", "", NCDB$TNM_PATH_N) #Remove white spaces
  NCDB$TNM_PATH_M <- gsub("[[:space:]]", "", NCDB$TNM_PATH_M) #Remove white spaces
  
  #Re-categorize ideal regional radiation - RAD_REGIONAL_DOSE_CGY (0, 4500-5040)
  attach(NCDB)
  NCDB$IDEAL_REG_RAD[RAD_REGIONAL_DOSE_CGY >= 4500 & RAD_REGIONAL_DOSE_CGY <= 5040 ] <- 1 #Ideal dose
  NCDB$IDEAL_REG_RAD[RAD_REGIONAL_DOSE_CGY < 4500 | RAD_REGIONAL_DOSE_CGY > 5040 ] <- 2 #Inappropriate dose
  NCDB$IDEAL_REG_RAD[RAD_REGIONAL_DOSE_CGY == 0] <- 0 #No radiation
  detach(NCDB)
  
  #Re-categorize ideal boost radiation - RAD_BOOST_DOSE_CGY (0, 540)
  attach(NCDB)
  NCDB$IDEAL_BOOST_RAD[RAD_BOOST_DOSE_CGY == 540] <- 1 #Ideal dose
  NCDB$IDEAL_BOOST_RAD[RAD_BOOST_DOSE_CGY != 540] <- 2 #Inappropriate dose
  NCDB$IDEAL_BOOST_RAD[RAD_BOOST_DOSE_CGY == 0] <- 0 #No radiation
  detach(NCDB)
  
  attach(NCDB)
  NCDB$RAD_RECEIVED[IDEAL_REG_RAD == 0 & IDEAL_BOOST_RAD == 0] <- 0 #No radiation
  NCDB$RAD_RECEIVED[IDEAL_REG_RAD == 1 | IDEAL_BOOST_RAD == 1] <- 1 #Ideal dose
  NCDB$RAD_RECEIVED[IDEAL_REG_RAD == 2 | IDEAL_BOOST_RAD == 2] <- 2 #Inappropriate radiation
  detach(NCDB)
  
  #Re-categorize whether chemotherapy was given
  attach(NCDB)
  NCDB$CHEMO_RECEIVED[RX_SUMM_CHEMO == 0 | RX_SUMM_CHEMO == 87 | RX_SUMM_CHEMO == 86 | RX_SUMM_CHEMO == 85 | RX_SUMM_CHEMO == 82] <- 0 #Chemotherapy not given
  NCDB$CHEMO_RECEIVED[RX_SUMM_CHEMO == 1 | RX_SUMM_CHEMO == 2 | RX_SUMM_CHEMO == 3] <- 1 #Chemotherapy received
  detach(NCDB)
  
  #No neoadjuvant chemoradiation
  NCDB <- subset(NCDB, RAD_RECEIVED == 0)
  NCDB <- subset(NCDB, CHEMO_RECEIVED == 0)
  
  #At this point we are looking at early-stage invasive rectal adenocarcinoma (2010-2015) treated at one facility with LAR/APR and who did not receive chemoradiation
  #Create table of volume frequency during the period (not annual) - this will be used as a marker for facility experience
  facility_vol <- table(NCDB$PUF_FACILITY_ID)
  NCDB$FACILITY_VOL <- as.numeric(facility_vol[NCDB$PUF_FACILITY_ID])
  
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

  #Create new approach variable
  attach(NCDB)
  NCDB$APPROACH[RX_HOSP_SURG_APPR_2010 == 5] <- 0 #Open
  NCDB$APPROACH[RX_HOSP_SURG_APPR_2010 == 3] <- 1 #Laparoscopic
  NCDB$APPROACH[RX_HOSP_SURG_APPR_2010 == 1] <- 2 #Robotic
  NCDB$APPROACH[RX_HOSP_SURG_APPR_2010 == 2 | RX_HOSP_SURG_APPR_2010 == 4] <- 3 #Converted to open
  detach(NCDB)
  
  NCDB <- subset(NCDB, APPROACH < 3)
  NCDB <- subset(NCDB, !is.na(NCDB$SURG_DISCHARGE_DAYS)) 
  NCDB <- subset(NCDB, GRADE != 9) #Differentiation, filter those that are unknown
  NCDB <- subset(NCDB, RX_SUMM_SURGICAL_MARGINS <= 3) #Exclude unknown or indeterminate margin status
  NCDB <- subset(NCDB, READM_HOSP_30_DAYS != 9) #Readmission unknown
  NCDB <- subset(NCDB, REGIONAL_NODES_EXAMINED >= 1 & REGIONAL_NODES_EXAMINED <= 90)
  NCDB <- subset(NCDB, PUF_30_DAY_MORT_CD != 9) #Exclude unknown mortality
  NCDB <- subset(NCDB, PUF_90_DAY_MORT_CD != 9) #Exclude unknown mortality
  
  #Re-categorize unplanned readmissions
  attach(NCDB)
  NCDB$UNPLAN_READMIT[READM_HOSP_30_DAYS == 0 | READM_HOSP_30_DAYS == 2] <- 0 #No unplanned readmissions
  NCDB$UNPLAN_READMIT[READM_HOSP_30_DAYS == 1 | READM_HOSP_30_DAYS == 3] <- 1 #Unplanned readmission occurred
  detach(NCDB)
  
  #Categorize margin status
  attach(NCDB)
  NCDB$MARGIN_STATUS[NCDB$RX_SUMM_SURGICAL_MARGINS == 0] <- 0 #Negative margin
  NCDB$MARGIN_STATUS[NCDB$RX_SUMM_SURGICAL_MARGINS > 0] <- 1 #Positive margin
  detach(NCDB)
  
  #Categorize node status
  attach(NCDB)
  NCDB$NODE_STATUS[NCDB$REGIONAL_NODES_EXAMINED < 12] <- 0 #Poor nodal harvist
  NCDB$NODE_STATUS[NCDB$REGIONAL_NODES_EXAMINED >= 12] <- 1 #Good nodal harvist
  detach(NCDB)
  
  #Categorized oncologic status
  attach(NCDB)
  NCDB$ONC_SUCCESS[MARGIN_STATUS == 0 & REGIONAL_NODES_EXAMINED >= 12] <- 1 #Negative margins and nodes >=12 examined, success
  NCDB$ONC_SUCCESS[MARGIN_STATUS == 1 | REGIONAL_NODES_EXAMINED < 12] <- 0 #Positive margins and nodes <12 examined, unsuccessful
  detach(NCDB)
  
  #Categorize age  <50, 50-65, >65
  attach(NCDB)
  NCDB$AGE_CAT_50_65[NCDB$AGE < 50] <- 1 #Age <50
  NCDB$AGE_CAT_50_65[NCDB$AGE >= 50 & NCDB$AGE < 65] <- 2 #Age 50-65
  NCDB$AGE_CAT_50_65[NCDB$AGE >= 65] <- 3 #Age >65
  detach(NCDB)
  
  #Categorize age <55, 55-65, 65-75, >75, mased on 1st quartile 54, median 64, 3rd - 75, min 20, max 90
  attach(NCDB)
  NCDB$AGE_CAT_55_65_75[NCDB$AGE < 55] <- 1 #Age <55
  NCDB$AGE_CAT_55_65_75[NCDB$AGE >= 55 & NCDB$AGE < 65] <- 2 #Age 55-65
  NCDB$AGE_CAT_55_65_75[NCDB$AGE >= 65 & NCDB$AGE < 75] <- 3 #Age 65-75
  NCDB$AGE_CAT_55_65_75[NCDB$AGE >= 75] <- 4 #Age >75
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
  
  var_names <- c("NO_HSD_QUAR_16", "MED_INC_QUAR_16", "AGE_CAT_55_65_75", "AGE_CAT_50_65", "AGE_CAT_65", "RURAL_METRO", "PUF_30_DAY_MORT_CD", "PUF_90_DAY_MORT_CD", "INSURANCE_STATUS", "TNM_CLIN_T_COMB", "TNM_PATH_T", "TNM_PATH_N", "TNM_PATH_M", "ONC_SUCCESS","MARGIN_STATUS","APPROACH","UNPLAN_READMIT", "SEX", "GRADE", "ETHNICITY","CDCC_TOTAL_BEST", "FACILITY_TYPE_CD", "FACILITY_LOCATION_CD")
  NCDB[,var_names] <- lapply(NCDB[,var_names], factor)
  NCDB$DX_LASTCONTACT_DEATH_MONTHS <- as.numeric(NCDB$DX_LASTCONTACT_DEATH_MONTHS)

#Evaluate variables to see if they make sense
summary(NCDB$AGE)
table(NCDB$SEX, useNA = 'ifany')
table(NCDB$ETHNICITY, useNA = 'ifany')
table(NCDB$INSURANCE_STATUS, useNA = 'ifany')
summary(NCDB$CDCC_TOTAL_BEST)
table(NCDB$TNM_CLIN_STAGE_GROUP, useNA = 'ifany')

t_table <- table(NCDB$APPROACH, NCDB$TNM_CLIN_T_COMB)
prop.table(t_table, 2)

table(NCDB$GRADE, useNA = 'ifany')
table(NCDB$HISTOLOGY_TYPE, useNA = 'ifany')

#Grouping
table(NCDB$APPROACH, useNA = 'ifany')

#Outcomes
table(NCDB$MARGIN_STATUS, useNA = 'ifany')
margin_status <- table(NCDB$APPROACH, NCDB$MARGIN_STATUS)
margin_status
prop.table(margin_status, 1)
summary(NCDB$REGIONAL_NODES_EXAMINED)
table(NCDB$NODE_STATUS, useNA = 'ifany')
node_status <- table(NCDB$APPROACH, NCDB$NODE_STATUS)
node_status
prop.table(node_status, 1)
table(NCDB$ONC_SUCCESS, useNA = 'ifany')
onc_success_table <- table(NCDB$APPROACH, NCDB$ONC_SUCCESS)
onc_success_table
prop.table(onc_success_table, 1)

summary(NCDB$SURG_DISCHARGE_DAYS)
table(NCDB$UNPLAN_READMIT, useNA = 'ifany')
table(NCDB$PUF_30_DAY_MORT_CD, useNA = 'ifany')
table(NCDB$PUF_90_DAY_MORT_CD, useNA = 'ifany')

#Propensity 
NCDB$APPROACH <- factor(NCDB$APPROACH, levels = c(0,1,2), labels = c("Open","Laparoscopic","Robotic"))
NCDB$APPROACH <- relevel(NCDB$APPROACH, ref = "Open")

set.seed(1)
mnps.NCDB <- mnps(APPROACH~AGE+SEX+INSURANCE_STATUS+ETHNICITY+CDCC_TOTAL_BEST+TNM_CLIN_T_COMB+GRADE+FACILITY_VOL,
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
file_name <- paste("Manuscript/Balance Table - ", format(today, format="%m%d%Y"), ".csv", sep = "")
write.csv(bal.table(mnps.NCDB, digits = 3), file = paste(dir, file_name,sep = ""))
summary(mnps.NCDB)

NCDB$weights <- get.weights(mnps.NCDB, stop.method = "ks.mean")
design.mnps <- svydesign(ids=~1, weights=~weights, data=NCDB)

#Margin status
margin <- svyglm(MARGIN_STATUS ~ as.factor(APPROACH), design = design.mnps, family=quasibinomial)
summary(margin)
exp(cbind(OR = coef(margin), confint(margin)))
#Nodes status
node <- svyglm(NODE_STATUS ~ as.factor(APPROACH), design = design.mnps, family=quasibinomial)
summary(node)
exp(cbind(OR = coef(node), confint(node)))
#Lymph nodes examined
nodes_exam <- svyglm(REGIONAL_NODES_EXAMINED ~ as.factor(APPROACH), design = design.mnps, family=gaussian())
summary(nodes_exam)
exp(cbind(OR = coef(nodes_exam), confint(nodes_exam)))
#Oncologic success
onc_success <- svyglm(ONC_SUCCESS ~ as.factor(APPROACH), design = design.mnps, family=quasibinomial)
summary(onc_success)
exp(cbind(OR = coef(onc_success), confint(onc_success)))
#Length of stay
los <- svyglm(SURG_DISCHARGE_DAYS ~ as.factor(APPROACH), design = design.mnps, family=gaussian())
summary(los)
exp(cbind(OR = coef(los), confint(los)))
#Readmissions
readm <- svyglm(UNPLAN_READMIT ~ as.factor(APPROACH), design = design.mnps, family=quasibinomial)
summary(readm)
exp(cbind(OR = coef(readm), confint(readm)))
#30-day mortality
mort_30 <- svyglm(PUF_30_DAY_MORT_CD ~ as.factor(APPROACH), design = design.mnps, family=quasibinomial)
summary(mort_30)
exp(cbind(OR = coef(mort_30), confint(mort_30)))
#90-day mortality
mort_90 <- svyglm(PUF_90_DAY_MORT_CD ~ as.factor(APPROACH), design = design.mnps, family=quasibinomial)
summary(mort_90)
exp(cbind(OR = coef(mort_90), confint(mort_90)))

#Plot odds ratios - Open vs. Laparoscopic
text.title <- element_text(size = 12, margin = margin(t = 0, r = 0, b = 20, l = 0, unit = "pt"), hjust = 0.5)
text.axis.x <- element_text(size = 10, margin = margin(t = 10, r = 0, b = 5, l = 0, unit = "pt"))
text.axis.y <- element_text(size = 10, margin = margin(t = 10, r = 10, b = 10, l = 0, unit = "pt"))
boxLabels = c("Margin status", "Node status", "Oncologic success", "Unplanned readmission", "90-day Mortality")
yAxis = length(boxLabels):1
df <- data.frame(yAxis = length(boxLabels):1,
                 boxOdds = c(as.numeric(exp(coef(margin))[2]), as.numeric(exp(coef(node))[2]), as.numeric(exp(coef(onc_success))[2]), as.numeric(exp(coef(readm))[2]), as.numeric(exp(coef(mort_90))[2])),
                 boxCILow = c(as.numeric(exp(confint(margin))[2,1]), as.numeric(exp(confint(node))[2,1]), as.numeric(exp(confint(onc_success))[2,1]), as.numeric(exp(confint(readm))[2,1]), as.numeric(exp(confint(mort_90))[2,1])),
                 boxCIHigh = c(as.numeric(exp(confint(margin))[2,2]), as.numeric(exp(confint(node))[2,2]), as.numeric(exp(confint(onc_success))[2,2]), as.numeric(exp(confint(readm))[2,2]), as.numeric(exp(confint(mort_90))[2,2])))
p.open_lap <- ggplot(df, aes(x = boxOdds, y = yAxis))
p.open_lap + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "black") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), plot.title = text.title, axis.text.y = text.axis.y, axis.text.x = text.axis.x, plot.margin = margin(t = 30, r = 20, b = 20, l = 20, unit = "pt")) +
  scale_y_continuous(breaks = yAxis, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,2,0.5)) +
  coord_trans(x = "log10") +
  labs(title = "Open vs. Laparoscopic Approach", y = "", x = "Odds ratio (log scale)") +
  annotate("text", x = 0.5, y = 0.25, label = "Favors laparoscopic") +
  annotate("text", x = 1.5, y = 0.25, label = "Favors open")
ggsave(filename = paste(dir,"Manuscript/Tables - Figures/Figure - 1a.pdf",sep = ""), dpi = 300, width = 8, height = length(boxLabels))

#Plot odds ratios - Open vs. Robotic
boxLabels = c("Margin status", "Node status", "Oncologic success", "Unplanned readmission", "90-day Mortality")
yAxis = length(boxLabels):1
df <- data.frame(yAxis = length(boxLabels):1,
                 boxOdds = c(as.numeric(exp(coef(margin))[3]), as.numeric(exp(coef(node))[3]), as.numeric(exp(coef(onc_success))[3]), as.numeric(exp(coef(readm))[3]), as.numeric(exp(coef(mort_90))[3])),
                 boxCILow = c(as.numeric(exp(confint(margin))[3,1]), as.numeric(exp(confint(node))[3,1]), as.numeric(exp(confint(onc_success))[3,1]), as.numeric(exp(confint(readm))[3,1]), as.numeric(exp(confint(mort_90))[3,1])),
                 boxCIHigh = c(as.numeric(exp(confint(margin))[3,2]), as.numeric(exp(confint(node))[3,2]), as.numeric(exp(confint(onc_success))[3,2]), as.numeric(exp(confint(readm))[3,2]), as.numeric(exp(confint(mort_90))[3,2])))
p.open_rob <- ggplot(df, aes(x = boxOdds, y = yAxis))
p.open_rob + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "black") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), plot.title = text.title, axis.text.y = text.axis.y, axis.text.x = text.axis.x, plot.margin = margin(t = 30, r = 20, b = 20, l = 20, unit = "pt")) +
  scale_y_continuous(breaks = yAxis, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0,2,0.5)) +
  coord_trans(x = "log10") +
  labs(title = "Open vs. Robotic Approach", y = "", x = "Odds ratio (log scale)") +
  annotate("text", x = 0.45, y = 0.25, label = "Favors robotic") +
  annotate("text", x = 2.0, y = 0.25, label = "Favors open")
ggsave(filename = paste(dir,"Manuscript/Tables - Figures/Figure - 1b.pdf",sep = ""), dpi = 300, width = 8, height = length(boxLabels))

#Evaluate predictors of open status compared to laparoscopic and robotic
multi_nom <- multinom(APPROACH ~ SEX + AGE_CAT_55_65_75 + SEX + ETHNICITY + CDCC_TOTAL_BEST + TNM_CLIN_T_COMB + FACILITY_TYPE_CD + FACILITY_LOCATION_CD + NO_HSD_QUAR_16 + INSURANCE_STATUS + RURAL_METRO + MED_INC_QUAR_16, data = NCDB)
summary(multi_nom)
z_multi_nom <- summary(multi_nom)$coefficients/summary(multi_nom)$standard.errors
p_multi_nom <- (1 - pnorm(abs(z_multi_nom), 0, 1))*2
exp(coef(multi_nom))
exp(confint(multi_nom))
# 
# #Open vs. laparoscopic
# text.title <- element_text(size = 12, margin = margin(t = 0, r = 0, b = 20, l = 0, unit = "pt"), hjust = 0.5)
# text.axis.x <- element_text(size = 10, margin = margin(t = 10, r = 0, b = 5, l = 0, unit = "pt"))
# text.axis.y <- element_text(size = 10, margin = margin(t = 10, r = 10, b = 10, l = 0, unit = "pt"))
# boxLabels = c("Sex", "Age 50-65", "Age >65", "Non-hispanic black", "Asian/Pacific Islander", "Hispanic", "Other", "CDCC Score 1", "CDCC Score 2", "CDCC Score >3", "T2", "")
# yAxis = length(boxLabels):1
# df <- data.frame(yAxis = length(boxLabels):1,
#                  boxOdds = c(as.numeric(exp(coef(margin))[2]), as.numeric(exp(coef(node))[2]), as.numeric(exp(coef(onc_success))[2]), as.numeric(exp(coef(readm))[2]), as.numeric(exp(coef(mort_90))[2])),
#                  boxCILow = c(as.numeric(exp(confint(margin))[2,1]), as.numeric(exp(confint(node))[2,1]), as.numeric(exp(confint(onc_success))[2,1]), as.numeric(exp(confint(readm))[2,1]), as.numeric(exp(confint(mort_90))[2,1])),
#                  boxCIHigh = c(as.numeric(exp(confint(margin))[2,2]), as.numeric(exp(confint(node))[2,2]), as.numeric(exp(confint(onc_success))[2,2]), as.numeric(exp(confint(readm))[2,2]), as.numeric(exp(confint(mort_90))[2,2])))
# p.open_lap <- ggplot(df, aes(x = boxOdds, y = yAxis))
# p.open_lap + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
#   geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
#   geom_point(size = 3.5, color = "black") +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(), plot.title = text.title, axis.text.y = text.axis.y, axis.text.x = text.axis.x, plot.margin = margin(t = 30, r = 20, b = 20, l = 20, unit = "pt")) +
#   scale_y_continuous(breaks = yAxis, labels = boxLabels) +
#   scale_x_continuous(breaks = seq(0,2,0.5)) +
#   coord_trans(x = "log10") +
#   labs(title = "Open vs. Laparoscopic Approach", y = "", x = "Odds ratio (log scale)") +
#   annotate("text", x = 0.5, y = 0.25, label = "Favors laparoscopic") +
#   annotate("text", x = 1.5, y = 0.25, label = "Favors open")
# ggsave(filename = paste(dir,"Manuscript/Tables - Figures/Figure - 1a.pdf",sep = ""), dpi = 300, width = 8, height = length(boxLabels))

#Evaluate predictors of 90-day mortality in only open cohort
NCDB_open <- subset(NCDB, APPROACH == 0 | APPROACH == "Open")

chi_var_names <- c("SEX", "AGE_CAT_55_65_75", "ETHNICITY", "CDCC_TOTAL_BEST", "TNM_CLIN_T_COMB", "GRADE", "FACILITY_TYPE_CD", "FACILITY_LOCATION_CD", "NO_HSD_QUAR_16", "INSURANCE_STATUS", "RURAL_METRO", "MED_INC_QUAR_16")
chi_tests_open <- lapply(NCDB_open[,chi_var_names], function(x) chisq.test(table(x, NCDB_open$PUF_90_DAY_MORT_CD)))

for (var in chi_var_names) {
  print(var)
  print(chi_tests_open[[var]])
  if(chi_tests_open[[var]]$p.value < 0.05) {
    print(table(NCDB_open$PUF_90_DAY_MORT_CD, NCDB_open[[var]]))
  }
}

mort_90_pred_open <- glm(NCDB_open$PUF_90_DAY_MORT_CD ~ NCDB_open$AGE_CAT_55_65_75 + NCDB_open$CDCC_TOTAL_BEST + NCDB_open$INSURANCE_STATUS, family = binomial)
summary(mort_90_pred_open)
exp(cbind(OR = coef(mort_90_pred_open), confint(mort_90_pred_open)))

# #Evaluate predictors of 90-day mortality in only laparoscopic cohort
# NCDB_lap <- subset(NCDB, APPROACH == 1)
# 
# chi_var_names <- c("SEX", "AGE_CAT_50_65", "ETHNICITY", "CDCC_TOTAL_BEST", "TNM_CLIN_T_COMB", "GRADE", "FACILITY_TYPE_CD", "FACILITY_LOCATION_CD", "NO_HSD_QUAR_16", "INSURANCE_STATUS", "RURAL_METRO", "MED_INC_QUAR_16")
# chi_tests_lap <- lapply(NCDB_lap[,chi_var_names], function(x) chisq.test(table(x, NCDB_lap$PUF_90_DAY_MORT_CD)))
# 
# for (var in chi_var_names) {
#   print(var)
#   print(chi_tests_lap[[var]])
#   if(chi_tests_lap[[var]]$p.value < 0.05) {
#     print(table(NCDB_lap$PUF_90_DAY_MORT_CD, NCDB_lap[[var]]))
#   }
# }
# 
# mort_90_pred_lap <- glm(NCDB_lap$PUF_90_DAY_MORT_CD ~ NCDB_lap$SEX + NCDB_lap$AGE_CAT_50_65 + NCDB_lap$ETHNICITY + NCDB_lap$CDCC_TOTAL_BEST + NCDB_lap$TNM_CLIN_T_COMB + NCDB_lap$GRADE + NCDB_lap$FACILITY_TYPE_CD + NCDB_lap$FACILITY_LOCATION_CD + NCDB_lap$NO_HSD_QUAR_16 + NCDB_lap$INSURANCE_STATUS + NCDB_lap$RURAL_METRO + NCDB_lap$MED_INC_QUAR_16, family = binomial)
# summary(mort_90_pred_lap)
# exp(cbind(OR = coef(mort_90_pred_lap), confint(mort_90_pred_lap)))
# 
# #Evaluate predictors of 90-day mortality in only robotic cohort
# NCDB_rob <- subset(NCDB, APPROACH == 2)
# 
# chi_var_names <- c("SEX", "AGE_CAT_50_65", "ETHNICITY", "CDCC_TOTAL_BEST", "TNM_CLIN_T_COMB", "GRADE", "FACILITY_TYPE_CD", "FACILITY_LOCATION_CD", "NO_HSD_QUAR_16", "INSURANCE_STATUS", "RURAL_METRO", "MED_INC_QUAR_16")
# chi_tests_rob <- lapply(NCDB_rob[,chi_var_names], function(x) chisq.test(table(x, NCDB_rob$PUF_90_DAY_MORT_CD)))
# 
# for (var in chi_var_names) {
#   print(var)
#   print(chi_tests_rob[[var]])
#   if(chi_tests_rob[[var]]$p.value < 0.05) {
#     print(table(NCDB_rob$PUF_90_DAY_MORT_CD, NCDB_rob[[var]]))
#   }
# }
# 
# mort_90_pred_rob <- glm(NCDB_rob$PUF_90_DAY_MORT_CD ~ NCDB_rob$SEX + NCDB_rob$AGE_CAT_50_65 + NCDB_rob$ETHNICITY + NCDB_rob$CDCC_TOTAL_BEST + NCDB_rob$TNM_CLIN_T_COMB + NCDB_rob$GRADE + NCDB_rob$FACILITY_TYPE_CD + NCDB_rob$FACILITY_LOCATION_CD + NCDB_rob$NO_HSD_QUAR_16 + NCDB_rob$INSURANCE_STATUS + NCDB_rob$RURAL_METRO + NCDB_rob$MED_INC_QUAR_16, family = binomial)
# summary(mort_90_pred_rob)
# exp(cbind(OR = coef(mort_90_pred_rob), confint(mort_90_pred_rob)))

#Survival by approach, need to exclude those who were upstaged after surgery
NCDB <- subset(NCDB, TNM_PATH_T == "p1" | TNM_PATH_T == "p2")
NCDB <- subset(NCDB, TNM_PATH_N == "p0")
NCDB <- subset(NCDB, TNM_PATH_M == "")

set.seed(1)
mnps.NCDB_2 <- mnps(APPROACH~AGE+SEX+INSURANCE_STATUS+ETHNICITY+CDCC_TOTAL_BEST+TNM_PATH_T+GRADE+FACILITY_VOL,
                  data=as.data.frame(NCDB),
                  n.trees = 3000,
                  stop.method = c("ks.mean"),
                  estimand = "ATE",
                  verbose = FALSE)

NCDB$weights_2 <- get.weights(mnps.NCDB_2, stop.method = "ks.mean")
design.mnps_2 <- svydesign(ids=~1, weights=~weights_2, data=NCDB)

#Survival by approach (weighted)
km.by.approach <- svykm(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS == 0) ~ APPROACH, design = design.mnps_2, se = TRUE)
text.axis.y <- element_text(size = 10, margin = margin(t = 10, r = 10, b = 10, l = 5, unit = "pt"))
p.km <- svyjskm(km.by.approach, ci=TRUE, ylims=c(0.6,1.0), xlims=c(0, 60), xlabs = "Time-to-event (months)", ylabs = "Survival probability", linecols="black", legendposition = c(0.2, 0.2), main="Survival analysis by surgical approach (weighted)", ystrataname = "Approach", ystratalabs = c("Open","Laparoscopic","Robotic"), dashed=TRUE)
p.km + theme_bw() +
  theme(panel.grid.minor = element_blank(), plot.title = text.title, axis.text.y = text.axis.y, axis.text.x = text.axis.x, plot.margin = margin(t = 30, r = 20, b = 20, l = 20, unit = "pt"))
ggsave(filename = paste(dir,"Manuscript/Tables - Figures/Figure - 2.pdf",sep = ""), dpi = 300, width = 8, height = 6)
#coxphfit <- coxph(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS == 0) ~ APPROACH, data=NCDB, weight=weights)
#summary(coxphfit)
#svycoxfit <- svycoxph(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS == 0) ~ APPROACH, design = design.mnps_2)
#summary(svycoxfit)
#plot(km.by.approach, pars=list(lty=c(1,2,3)), main = "Survival analysis - Surgical approach (Weighted)", xlab = "Months")
#legend("bottomleft",lty = c(1,2,3), legend = c("Laparoscopic", "Open", "Robotic"))
#res.cox <- coxph(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS == 0) ~ strata(APPROACH)+AGE+SEX+INSURANCE_STATUS+ETHNICITY+CDCC_TOTAL_BEST+TNM_CLIN_T_COMB+GRADE+FACILITY_VOL+MARGIN_STATUS+NODE_STATUS+ONC_SUCCESS, data=NCDB)
#summary(res.cox)
#ggadjustedcurves(res.cox, data=NCDB, size = 1)
