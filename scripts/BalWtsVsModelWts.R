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
data <- read.csv("/Users/haedi/Library/CloudStorage/Box-Box/Data/BalWts/IbuAkiClean.csv")  #Path to data

# Treatment variable: pain, 1= treated (ibuprofen), 0= control (oxycodone)
# Outcome variable: kEver, 1= acute kindey Injury (AKI), 0= no AKI

# List of covariates
covs <- c("age", "sex", "admType", "presentation.ed", "presentation.icu",
          "presentation.or", "presentation.floor", "presentation.other", 
          "priorLos", "icuCurrent", "chf",  "dm.no", "dm.noncomp", "dm.comp", 
          "ckd", "indexGFR", "loopBase", "vancoBase") 

covs_int <- c(covs, "-1") # include -1 for intercept in modeling steps

#################################################################################
#2. Estimate Model weights

#2a. Fit Propensity Score (PS) model
psmod <- glm(reformulate(covs_int, response = "pain"), family = binomial(), data = data)
data$PS <- predict(psmod, type = "response")

#2b. Estimate Model-based Inverse Probability of Treatment Weights (IPTW)
data <- data %>%
  mutate(MW = ifelse(pain ==1, 1/PS, 1/(1-PS)))

#2c. Evaluate balance
MW.balance <- bal.plots(data, "MW", 'pain', covs)

#2d. Calculate effective sample size (ESS)
MW.ess <- ess(data, "pain", "MW")
summary(data$MW)
sd(data$MW)
sum(data$PS < 0.05)/nrow(data) # proportion of PS < 0.05

#2e. PS overlap plot
# prepare data for ggplot
ps.overlap.data <-data %>%
  select(pain, PS) %>%
  mutate(pain = as.factor(pain))

# label the levels of the pain variable
levels(ps.overlap.data$pain) <- c("Control", "Treated")

ggplot(ps.overlap.data, aes(x=PS, fill = factor(pain)))+
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("blue", "red"), name = "Treatment Group" , ) +
  labs(x = "Propensity Score", y = "Density") +
  ggtitle("Distribution of Propensity Score by Treatment") + 
  theme_minimal()

#################################################################################
#3. Estimate Balancing Weights
#3a. Prepare data for BalanceR
basis <- reformulate(covs_int) # prepare a formula object                       
X <- scale(model.matrix(as.formula(basis), data)) # prepare a scaled matrix 
# scaling is needed to calculate the weights, since they target a mean of 0
trt <- data$pain
n <- nrow(data)

# estimate data-driven, initial hyperparameter
data.ctrl <- data %>% filter(pain==0)
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
  BW.balance <- bal.plots(data, paste0("BW", lambdas[i]), 'pain', covs)
  PBRs[[i]] <- BW.balance$pbr
}

#convert list to dataframe with lambdas as column names
PBRs <- as.data.frame(PBRs)
colnames(PBRs) <- paste0("BW", lambdas)
PBRs

# Save the base case lambda for the balance plot
BW.balance <- bal.plots(data, "BW0.05", 'pain', covs)

#3d. Calculate effective sample size for each value of lambda 
BW.ess <- data.frame(matrix(NA, nrow = 3, ncol = length(lambdas)))
rownames(BW.ess) <- c("Treatment", "Control", "Total")
colnames(BW.ess) <-  paste0("BW", lambdas)
for (i in 1:length(lambdas)){
    df <- ess(data, "pain", paste0("BW", lambdas[i]))
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
levels(BW.bal.overall$covariate) <- 
  c("Admission Type Surgery", "Age", "Heart Failure", "CKD", 
    "Diabetes Mellitus - Complicated", "Diabetes Mellitus - None", 
    "Diabetes Mellitus - Uncomplicated", "ICU - Current", "eGFR", 
    "Loop Diuretic", "Admitted to ED", "Admitted to Medical Floor", 
    "Admitted to ICU", "Admitted to OR", "Admitted to Other", 
    "Prior Length of Stay", "Sex", "Vancomycin")

#4b. Plot
plot <- ggplot(data = BW.bal.overall, aes(x = std.dif, y = covariate, 
                                           shape = factor(contrast) 
                                           #, color = factor(contrast) # for color instead of shapes, toggle as needed this,  geom_point, scale_color_manual
)) +
  #geom_point(size = 2, shape = 16) +
  geom_point(size = 2, color = "black")+
  scale_shape_manual(name = "Contrast", values = c(5, 1, 15), 
                     breaks = c('Unweighted', "Model-Based Weights","Balancing Weights" )) +
  #scale_color_manual(name = "Contrast", values = c("darkred", "darkblue", "darkgreen")) +
  xlab("Standardized Difference") + ylab("Covariates") +
  ggtitle("Balance in Weighted Population") +
  scale_y_discrete(limits = rev(levels(BW.bal.overall$covariate))) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.1) +
  geom_vline(xintercept = -0.1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
#guides(shape = FALSE)

plot

#################################################################################
#5. Create effective sample size table
og.n.tx <- table(data$pain)[[2]]
og.n.cx<- table(data$pain)[[1]]
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


