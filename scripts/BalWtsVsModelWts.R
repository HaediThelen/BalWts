# Balancing Weights Methods Review

#Libraries
library(foreign)
library(balancer)
library(dplyr)
library(ggplot2)
library(sandwich)


#1. Load and Prepare Data
#1a. label the needed variables 
# Load Data
data <- read.csv("/Users/haedi/Library/CloudStorage/Box-Box/Data/BalWts/IbuAkiClean.csv")

# List of covs to study
covs <- c("age", "sex", "admType", 
          "presentation.ed", "presentation.icu",
          "presentation.or", "presentation.floor", "presentation.other", "priorLos", "icuCurrent",
          "chf",  "dm.no", "dm.noncomp", "dm.comp", "ckd", 
          "indexGFR", "loopBase",
          "vancoBase", "-1")

#2. Estimate Model weights
#2a. Fit PS model
psmod <- glm(reformulate(covs, response = "pain"), family = binomial(), data = data)
data$PS <- predict(psmod, type = "response")

#2b. Estimate weights
# Calculate Model-based IPTW (ATE) 
data <- data %>%
  mutate(MW = ifelse(pain ==1, 1/PS, 1/(1-PS)))

#2c. Evaluate balance
covs <- c("age", "sex", "admType", 
          "presentation.ed", "presentation.icu",
          "presentation.or", "presentation.floor", "presentation.other", "priorLos", "icuCurrent",
          "chf",  "dm.no", "dm.noncomp", "dm.comp", "ckd", 
          "indexGFR", "loopBase",
          "vancoBase")

source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/balance-plot.R")
MW.balance <- bal.plots(data, "MW", 'pain', covs)

#2d. Calculate effective sample size
source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/ess-function.R")
MW.ess <- ess(data, "pain", "MW")
summary(data$MW)
sd(data$MW)

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

#3. Estimate Balancing Weights
#3a. Prepare data for BalanceR
covs <- c("age", "sex", "admType", 
          "presentation.ed", "presentation.icu",
          "presentation.or", "presentation.floor", "presentation.other", "priorLos", "icuCurrent",
          "chf",  "dm.no", "dm.noncomp", "dm.comp", "ckd", 
          "indexGFR", "loopBase",
          "vancoBase", "-1")

basis <- reformulate(covs) # prepare a formula object                       
X <- scale(model.matrix(as.formula(basis), data)) # prepare a scaled matrix 
# scaling is needed to calculate the weights, since they target a mean of 0
trt <- data$pain
n <- nrow(data)

# estimate data-driven, initial hyperparameter
data.ctrl <- data %>% filter(pain==0)
lambda.reg <- lm(reformulate(covs, response = "kEver"), data=data.ctrl)
var(lambda.reg$resid) # initial hyperparameter = 0.05

#3b. Estimate Balancing weights for multiple lambdas
lambdas <- c(0,  0.01, 0.05, 0.5, 1, 5, 10, 50)

# For loop that estimates balance for different values of lambda
for (i in 1:length(lambdas)){
  print(i)
  out.pain.q <- multilevel_ate_qp(X, trt, rep(1,n), lambda = lambdas[i], lowlim = 0, uplim = 1,  
                                  verbose= TRUE, exact_global = TRUE, scale_sample_size = FALSE)
  data[[paste0("BW", lambdas[i])]] <- out.pain.q$weights
}

# Process
summary(data$BW0.05)
sd(data$BW0.05)

#3c. Evaluate Balance
covs <- c("age", "sex", "admType", 
          "presentation.ed", "presentation.icu",
          "presentation.or", "presentation.floor", "presentation.other", "priorLos", "icuCurrent",
          "chf",  "dm.no", "dm.noncomp", "dm.comp", "ckd", 
          "indexGFR", "loopBase",
          "vancoBase")

# View SMD plots and calculate PBR for each value of lambda save to list PBRs
source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/balance-plot.R")
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

# Save the initial lamba for the final balance plot
BW.balance <- bal.plots(data, "BW0.05", 'pain', covs)

#3d. Calculate effective sample size for each value of lambda 
source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/ess-function.R")
# create an empty dataframe to store the ESS for each lambda
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

#4. Create overall Balance Plot using lambda = 0.05
#4a. Prepare a data frame with just the balance we need
# Relabel MW weights, and extract
MW.balance.data <- MW.balance$data 
levels(MW.balance.data$contrast) <- c("Model-Based Weights", "Unweighted")
MW.balance.data.wt <- MW.balance.data %>% filter(contrast =="Model-Based Weights")
# relabel BW weights
BW.balance.data <- BW.balance$data 
levels(BW.balance.data$contrast) <- c("Balancing Weights", "Unweighted")
# Add Model-based weights to BW weights
BW.bal.overall <- rbind(BW.balance.data, MW.balance.data.wt)
BW.bal.overall

# Change covariate names
levels(BW.bal.overall$covariate) <- c("Admission Type Surgery", "Age", "Heart Failure",
                                       "CKD", "Diabetes Mellitus - Complicated",
                                       "Diabetes Mellitus - None", "Diabetes Mellitus - Uncomplicated",
                                       "ICU - Current", "eGFR", "Loop Diuretic", "Admitted to ED", "Admitted to Medical Floor",
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

#5. Create Effective Sample size table
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

# Add MWPBR to PBRs
PBRs2 <- PBRs %>%
  mutate("Model-Based Weights" = round(MW.balance$pbr), Original = 0) %>%
  select(Original, "Model-Based Weights", everything())
  
row.names(PBRs2) <- c("PBR")
results <- rbind(ess.tab, PBRs2)
results <- results %>% arrange(row.names(results) != "PBR")
results 

write.csv(results, "/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/results/pbr.ess.tab.csv", row.names = T)
