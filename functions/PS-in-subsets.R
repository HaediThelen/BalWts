# Calculate PS within the subsets of a given variable

# Make subsets of variable 
subsetter <- function(data, var){
  split.data <- split(data, data[[var]])
  return(split.data)
}



# Calculate PS in each dataset
calc.ps <- function(subset.list, covs, response){
  for (i in seq_along(subset.list)){
    # fit PS model
    psmod <- glm(reformulate(covs, response = response), 
                 family = binomial(),
                 data = subset.list[[i]])
    subset.list[[i]]$PS <- predict(psmod, type = "response")
    print(i)
  }
  return(subset.list)
}

# Function to calculate ps within subsets of given variable
  # covs are the covariates, with a -1
  # response is the trt variable
  # var is the variable that indicates the subsetting
  ps.subset <- function(data, covs, response, var) {
    splits <- subsetter(data, var)
    splits.ps <- calc.ps(splits, covs, response)
    # Rejoin the subsets 
    combined.splits <- do.call(rbind, splits.ps)
    return(combined.splits)
  }


order <- c("ddiID", "kEver","kStage","rrt","preAkiStatus","indexCre","ibu",
          "ibuDose","pain","admType","sex","year","age","chf",
          "metsol","cpd","hiv","liver","nonMet","pvd","cva","mif","valve","htn",
          "arry","pCirc","obese","wtLoss","fluid","ckd","afib","osa","indexDate",
          "endDate","crossEndTime", "crossOver","ddiDurTime","crossDurTime",
          "priorLos","hydralazineBase","hctzBase","loopBase", "vancoBase",
          "bactrimBase","pressBase","hydromBase","morphineBase","apapBase",
          "fentBase","ivPainBase","pcaBase","infusionBase","tubeBase","ibu2Base",
          "rasBase","metopBase", "wbcBase","labnaBase","hgbBase","labkBase",
          "labclBase","platBase","icuCurrent", "periOp","baseVentCurrent",
          "baseVentEver","pTime","pTime100","pTime1000","sup", "ntxOther",
          "abBlocker","htnOther","abxNTX","gramNegBroad","gramNegNarrow","center",
          "presentation","bmi","bmiCat","race","indexGFR","dm","cancer", "bmi.q",
          "admType.cat","race.white","race.black","race.other","center.hup",
          "center.presb", "center.pa","presentation.ed","presentation.icu",
          "presentation.or","presentation.floor", "presentation.other",
          "periOp.no","periOp.0","periOp.1","periOp.2","periOp.3","dm.no",
          "dm.noncomp","dm.comp", "cancer.no","cancer.noncomp","cancer.metastatic",
          "sup.no","sup.h2ra","sup.ppi","bmi.adm",  "age1","age2","age3",
          "indexGFR1","indexGFR2","indexGFR3","indexGFR4", "BW_ATE","BW_ATT",
          "PS","IPW", "SMR","BW_ATE_admType","PS_admType","IPW_admType",
          "BW_ATE_bmi.q","PS_bmi.q","IPW_bmi.q","BW_ATE_bmi.qXadmType",
          "PS_bmi.qXadmType", "IPW_bmi.qXadmType")




