#################################################################################
###                                                                           ###
###    Balancing Weights Methods Review                                       ### 
###    Haedi Thelen, Todd Miano, Luke Keele                                   ###
###                                                                           ### 
###    This is the code for the companion paper "A primer on                  ### 
###    using balancing weights to estimate inverse probability weights"       ### 
###    We use a dataset extracted from an EHR that is described in:           ###
###    DOI: 10.34067/KID.0001432020                                           ###
###    This dataset contains several covariates and we can visualize          ###
###    the difference between not using weights, model-based weights          ###
###    and balancing weights. Here we clean the dataset prior to              ###
###    posting the datset  at                                                 ###
###                                                                           ###
#################################################################################

# Prepare data for Balance Weights Project
library(haven)
data <- read_dta("/Users/haedi/Library/CloudStorage/Box-Box/Projects/Summer Rotation/Balance weights dataset/nsaid_ras_dataset.dta")

#  prepare data
data <- data %>% 
  mutate(across(where(is.numeric), as.numeric)) %>%
  mutate(across(where(~ all(. %in% c(0, 1))), as.integer)) %>%
  mutate(presentation.ed = if_else(presentation == 1, 1, 0),
         presentation.icu = if_else(presentation == 2, 1, 0),
         presentation.or = if_else(presentation == 3, 1, 0),
         presentation.floor = if_else(presentation == 4,1,0),
         presentation.other = if_else(presentation == 5,1,0)) %>%
  mutate(periOp.no = if_else(periOp == 0, 1, 0),
         periOp.0 = if_else(periOp == 1, 1, 0),
         periOp.1 = if_else(periOp == 2, 1, 0),
         periOp.2 = if_else(periOp == 3,1,0),
         periOp.3 = if_else(periOp == 4,1,0)) %>%
  mutate(cancer.no = if_else(cancer == 0, 1, 0),
         cancer.noncomp = if_else(cancer == 1, 1, 0),
         cancer.metastatic = if_else(cancer == 2, 1, 0)) %>%
  mutate(dm.no = if_else(dm == 0, 1, 0),
         dm.noncomp = if_else(dm == 1, 1, 0),
         dm.comp = if_else(dm == 2, 1, 0)) %>%
  mutate(race.0 = if_else(race2 == 0, 1, 0),
         race.1 = if_else(race2 == 1, 1, 0),
         race.2 = if_else(race2 == 2, 1, 0)) %>%
  filter(ddiGroup %in% c(1, 2)) %>% # filter to patients with NSAID and either RASi or CCB
  mutate(nsaid_bp = if_else(ddiGroup == 1, 1, 0)) %>% #nsaid_bp, =1 if ddiGroup = 3 and 0 if ddiGroup = 4
  select('pain',"bp", "ddiGroup","nsaid_bp", "preAkiStatus", "admType", "sex", "age", "chf", 
         "metsol", "cpd", "hiv", "liver", "nonMet", "pvd", "cva", "mif", "valve", 
         "htn", "arry", "pCirc", "obese", "wtLoss", "fluid", "ckd", "afib", "osa", 
         "priorLos", "loopBase", "vancoBase", "wbcBase", "labnaBase", "hgbBase", 
         "labkBase", "labclBase", "platBase", "icuCurrent", "periOp.no", "periOp.0","periOp.1", "periOp.2", "periOp.3" 
         , "baseVentCurrent", 
         "bp", "ntxOther", "bBlocker", "abBlocker", "htnOther", "abxNTX", 
         "gramNegBroad", "gramNegNarrow", "location2", "presentation.ed", "presentation.icu",
         "presentation.or", "presentation.floor", "presentation.other", 
         "admissionFloor", "bmi2", "race2", 
         "race.0", "race.1", "race.2",
         "indexGFR", "dm", "dm.no", "dm.noncomp", "dm.comp",
         "cancer.no", "cancer.noncomp", "cancer.metastatic", "diurBase")

write.csv(data, "/Users/haedi/Library/CloudStorage/Box-Box/Data/BalWts/NSAID_RASClean.csv", row.names = F)
data.clean <- read.csv("/Users/haedi/Library/CloudStorage/Box-Box/Data/BalWts/NSAID_RASClean.csv")

