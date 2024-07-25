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
###    and balancing weights.                                                 ###
###                                                                           ###
#################################################################################

#Load libraries
if (!require("foreign")) install.packages("foreign"); library(foreign)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("sandwich")) install.packages("sandwich"); library(sandwich)

#Install BalanceR package
if (!require("devtools")) install.packages("devtools"); library(devtools)
if (!require("balancer")) devtools::install_github("ebenmichael/balancer")

# Load functions
setwd("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/") #Path to repo
source("./functions/balance-plot.R") 
source("./functions/ess-function.R")

#################################################################################
#1. Load and Prepare Data
# Download data from  ############## Zenodo Posting TBD ##############
data <- read.csv("/Users/haedi/Library/CloudStorage/Box-Box/Data/BalWts/NSAID_RASClean.csv")  #Path to data
# Treatment variable: nsaid_bp, 1= treated (RASi), 0= control (CCB)

# List of covariates
covs <- c("admType", "sex", "age", "chf", 
            "cpd", "pvd", "cva", "mif", "valve", 
           "htn", "arry", "ckd", "afib", 
           "priorLos", "loopBase", "vancoBase", "icuCurrent", 
           "ntxOther", "bBlocker", "abBlocker", "htnOther", "abxNTX", 
           "gramNegBroad", "gramNegNarrow", 
            "bmi2", 
           "race.0", "race.1", "race.2",
           "indexGFR", "dm.no", "dm.noncomp", "dm.comp"
           , "diurBase")
covs_int <- c(covs, "-1") # include -1 for intercept in modeling steps

#################################################################################
#2. Estimate Model weights

#2a. Fit Propensity Score (PS) model
psmod <- glm(reformulate(covs_int, response = "nsaid_bp"), family = binomial(), data = data)
data$PS <- predict(psmod, type = "response")

#2b. Estimate Model-based Inverse Probability of Treatment Weights (IPTW)
data <- data %>%
  mutate(MW = ifelse(nsaid_bp ==1, 1/PS, 1/(1-PS)))

#2c. Evaluate balance
MW.balance <- bal.plots(data, "MW", 'nsaid_bp', covs) #try contrasts other than pain.. 

#2d. Calculate effective sample size (ESS)
MW.ess <- ess(data, "nsaid_bp", "MW")
summary(data$MW)
sd(data$MW)
sum(data$PS < 0.1)/nrow(data) # proportion of PS < 0.05
sum(data$PS > (1-0.1))/nrow(data) # proportion of PS > 0.95
# calculate proporiton of nsaid_bp = 1 with PS < 0.1
sum(data$nsaid_bp == 1 & data$PS < 0.1)/sum(data$nsaid_bp == 1)
sum(data$nsaid_bp == 0 & data$PS > 0.95)/sum(data$nsaid_bp == 0)


#2e. PS overlap plot
# prepare data for ggplot
ps.overlap.data <-data %>%
  select(nsaid_bp, PS) %>%
  mutate(nsaid_bp = as.factor(nsaid_bp))

# label the levels of the pain variable
levels(ps.overlap.data$nsaid_bp) <- c("Control", "Treated")

ggplot(ps.overlap.data, aes(x=PS, fill = factor(nsaid_bp)))+
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("blue", "red"), name = "Treatment Group", ) +
  labs(x = "Propensity Score", y = "Density") +
  ggtitle("Distribution of Propensity Score by Treatment") + 
  theme_minimal()

#################################################################################
#3. Estimate Balancing Weights
#3a. Prepare data for BalanceR
basis <- reformulate(covs_int) # prepare a formula object                       
X <- scale(model.matrix(as.formula(basis), data)) # prepare a scaled matrix 
# scaling is needed to calculate the weights, since they target a mean of 0
trt <- data$nsaid_bp
n <- nrow(data)

# estimate data-driven, initial hyperparameter
data.ctrl <- data %>% filter(nsaid_bp==0)
lambda.reg <- lm(reformulate(covs_int, response = "kEver"), data=data.ctrl)
var(lambda.reg$resid) # initial hyperparameter = 0.05

#3b. Estimate Balancing weights for multiple lambdas
lambdas <- c(0,0.0001, 0.01, 0.05, 0.5, 1, 5, 10, 50)

# For loop that estimates balance for different values of lambda
for (i in 1:length(lambdas)){
  print(i)
  out.pain.q <- multilevel_ate_qp(X, trt, rep(1,n), lambda = lambdas[i], lowlim = 0, uplim = 1,  
                                  verbose= TRUE, exact_global = TRUE, scale_sample_size = FALSE)
  data[[paste0("BW", lambdas[i])]] <- out.pain.q$weights
}

# Evaluate base case weights
summary(data$BW0.05)
sd(data$BW0.05)

#3c. Evaluate Balance
# View SMD plots and calculate percent bias reduction (PBR) for each value of lambda save to list PBRs
PBRs <- list()
for (i in 1:length(lambdas)){
  print(i)
  BW.balance <- bal.plots(data, paste0("BW", lambdas[i]), 'nsaid_bp', covs)
  PBRs[[i]] <- BW.balance$pbr
}

#convert list to dataframe with lambdas as column names
PBRs <- as.data.frame(PBRs)
colnames(PBRs) <- paste0("BW", lambdas)
PBRs

# Save the base case lambda for the balance plot
BW.balance <- bal.plots(data, "BW0.05", 'nsaid_bp', covs)

#3d. Calculate effective sample size for each value of lambda 
BW.ess <- data.frame(matrix(NA, nrow = 3, ncol = length(lambdas)))
rownames(BW.ess) <- c("Treatment", "Control", "Total")
colnames(BW.ess) <-  paste0("BW", lambdas)
for (i in 1:length(lambdas)){
    df <- ess(data, "nsaid_bp", paste0("BW", lambdas[i]))
    BW.ess[1,i] <- df[1]
    BW.ess[2,i] <- df[2]
    BW.ess[3,i] <- df[3]
}
BW.ess

#################################################################################
#4. Create overall Balance Plot using lambda = 0.05
#4a. Prepare a data frame combining unweighted, model-based, and balancing weights

# Relabel MW weights, and extract
MW.balance.data <- MW.balance$data 
levels(MW.balance.data$contrast) <- c("Model-Based Weights", "Unweighted") 
MW.balance.data.wt <- MW.balance.data %>% filter(contrast =="Model-Based Weights")

# Relabel BW weights
BW.balance.data <- BW.balance$data 
levels(BW.balance.data$contrast) <- c("Balancing Weights", "Unweighted")

# Add Model-based weights to BW weights
BW.bal.overall <- rbind(BW.balance.data, MW.balance.data.wt)
BW.bal.overall

# Change covariate names
rename_vector <- c(
  "admType" = "Admitted to Surgery",
  "sex" = "Sex",
  "age" = "Age",
  "chf" = "Congestive Heart Failure",
  "cpd" = "Chronic Pulmonary Disease",
  "pvd" = "Peripheral Vascular Disease",
  "cva" = "Cerebrovascular Accident",
  "mif" = "Myocardial Infarction",
  "valve" = "Valvular Disease",
  "htn" = "Hypertension",
  "arry" = "Arrhythmia",
  "ckd" = "Chronic Kidney Disease",
  "afib" = "Atrial Fibrillation",
  "priorLos" = "Prior Length of Stay",
  "loopBase" = "Loop Diuretic",
  "vancoBase" = "Vancomycin",
  "icuCurrent" = "Current ICU Stay",
  "ntxOther" = "Other Nephrotoxin",
  "bBlocker" = "Beta Blocker",
  "abBlocker" = "Alpha Beta Blocker",
  "htnOther" = "Other Hypertension Med",
  "abxNTX" = "Nephrotoxic Antibiotic",
  "gramNegBroad" = "Broad-spectram Antibiotic",
  "gramNegNarrow" = "Narrow-spectrum Antibiotic",
  "bmi2" = "BMI",
  "race.0" = "Race-White",
  "race.1" = "Race-Black",
  "race.2" = "Race-Other",
  "indexGFR" = "Baseline eGFR",
  "dm.no" = "Diabetes Mellitus - None",
  "dm.noncomp" = "Diabetes Mellitus - Noncomplicated",
  "dm.comp" = "Diabetes Mellitus - Complicated",
  "diurBase" = "Other Diuretic"
)
levels(BW.bal.overall$covariate) <- rename_vector[levels(BW.bal.overall$covariate)]

# Reorder Variables
order <- c(
  # Demographics
  "Age", "Sex", "BMI", "Race-White", "Race-Black", "Race-Other",
  
  # Clinical Factors
  "Current ICU Stay", "Prior Length of Stay", "Admitted to Surgery",
  
  # Comorbid conditions
  "Atrial Fibrillation", "Cerebrovascular Accident", "Chronic Kidney Disease", 
  "Chronic Pulmonary Disease", "Congestive Heart Failure","Arrhythmia", "Diabetes Mellitus - None",
  "Diabetes Mellitus - Noncomplicated", "Diabetes Mellitus - Complicated",
  "Hypertension", "Myocardial Infarction", "Peripheral Vascular Disease", 
  "Valvular Disease",
  
  # Medications
  "Alpha Beta Blocker", "Beta Blocker", "Broad-spectram Antibiotic", 
  "Diuretic", "Loop Diuretic", "Narrow-spectrum Antibiotic", 
  "Nephrotoxic Antibiotic", "Other Hypertension Med", "Other Nephrotoxin", 
  "Vancomycin",
  
  # Labs
  "Index GFR"
)
BW.bal.overall$covariate <- factor(BW.bal.overall$covariate, levels = order)


#4b. Plot
plot <- ggplot(data = BW.bal.overall, aes(x = std.dif, y = covariate, 
                                          color = factor(contrast),
                                          shape = factor(contrast))) +
  geom_vline(xintercept = c(-0.1, 0, 0.1)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                     name = "Contrast", 
                     breaks = c('Unweighted', "Model-Based Weights", "Balancing Weights")) +
  scale_shape_manual(values = c(16, 17, 15),  # Circle, Triangle, Square
                     name = "Contrast",
                     breaks = c('Unweighted', "Model-Based Weights", "Balancing Weights")) +
  xlab("Standardized Difference") + ylab("Covariates") + 
  ggtitle("Balance in Weighted Population") + 
  scale_y_discrete(limits = rev(levels(BW.bal.overall$covariate))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(plot)
 
 #################################################################################
#5. Create effective sample size table
og.n.tx <- table(data$nsaid_bp)[[2]]
og.n.cx<- table(data$nsaid_bp)[[1]]
og.n.total <- og.n.cx+og.n.tx
original.n <-c(og.n.tx, og.n.cx, og.n.total)
original.n <- as.data.frame(original.n)
rownames(original.n) <- c("Treatment", "Control", "Total")
colnames(original.n) <- "Original"
original.n

# Transpose and label MW.ess
MW.ess <- as.data.frame(t(MW.ess))
rownames(MW.ess) <- c("Treatment", "Control", "Total")
colnames(MW.ess) <- "Model-Based Weights"
MW.ess

# Compile ess.tab 
ess.tab <- round(cbind(original.n, MW.ess, BW.ess))
ess.tab

# Add MW-PBR to BW-PBRs
PBRs2 <- PBRs %>%
  mutate("Model-Based Weights" = round(MW.balance$pbr), Original = 0) %>%
  select(Original, "Model-Based Weights", everything())
  
row.names(PBRs2) <- c("PBR")
results <- rbind(ess.tab, PBRs2)
results <- results %>% arrange(row.names(results) != "PBR")
results 

write.csv(results, "/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/results/pbr.ess.tab.csv", row.names = T)


