# Balancing Weights Methods Review
# 
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

#2. Estimate IPTW weights
#2a. Fit PS model
psmod <- glm(reformulate(covs, response = "pain"), family = binomial(), data = data)
data$PS <- predict(psmod, type = "response")

#2b. Estimate weights
# Calculate SMR (ATT)
data <- data %>%
  mutate(SMR = ifelse(pain ==1, 1, PS/(1-PS))) 
# Calculate IPWs (ATE)
data <- data %>%
  mutate(IPW = ifelse(pain ==1, 1/PS, 1/(1-PS)))

#2c. Evaluate balance
covs <- c("age", "sex", "admType", 
          "presentation.ed", "presentation.icu",
          "presentation.or", "presentation.floor", "presentation.other", "priorLos", "icuCurrent",
          "chf",  "dm.no", "dm.noncomp", "dm.comp", "ckd", 
          "indexGFR", "loopBase",
          "vancoBase")

source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/balance-plot.R")
IPW.balance <- bal.plots(data, "IPW", 'pain', covs)
SMR.balance <- bal.plots(data, "SMR", 'pain', covs)

#2d. Calculate effective sample size
source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/ess-function.R")
SMR.ess <- ess(data, "pain", "SMR")
IPW.ess <- ess(data, "pain", "IPW")

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

data.ctrl <- data %>% filter(pain==0)
lambda.reg <- lm(reformulate(covs, response = "kEver"), data=data.ctrl)
var(lambda.reg$resid)

#3b. Estimate Balancing weights for lambda = 0, 0.5, 1, 2 (maybe no?)
## process for ATT and save each in a new column

# Balance Weights ATT
# Calculate
out.pain <- multilevel_qp(X, trt, rep(1,n), lambda = 0.05, verbose= TRUE, 
                          exact_global = TRUE, scale_sample_size = FALSE)
# Process 
data$ATTwts <- pmax(out.pain$weights, 0) 
summary(data$ATTwts)
data$ATTwts[data$pain == 1] <- 1
summary(data$ATTwts)
sd(data$ATTwts)

# Balance Weights ATE
# Calculate 
out.pain.q <- multilevel_ate_qp(X, trt, rep(1,n), lambda = 1, lowlim = 0, uplim = 1,  
                                verbose= TRUE, exact_global = TRUE, scale_sample_size = FALSE)
# Process
data$ATEwts <- out.pain.q$weights
summary(data$ATEwts)
sd(data$ATEwts)

#3c. Evaluate Balance
covs <- c("age", "sex", "admType", 
          "presentation.ed", "presentation.icu",
          "presentation.or", "presentation.floor", "presentation.other", "priorLos", "icuCurrent",
          "chf",  "dm.no", "dm.noncomp", "dm.comp", "ckd", 
          "indexGFR", "loopBase",
          "vancoBase")

source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/balance-plot.R")
ATE.balance <- bal.plots(data, "ATEwts", 'pain', covs)
ATT.balance <- bal.plots(data, "ATTwts",'pain', covs)

#3d. Calculate effective sample size for each value of lambda 
source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/ess-function.R")
ATE.ess <- ess(data, "pain", "ATEwts")
ATT.ess <- ess(data, "pain", "ATTwts")



#4. Create overall Balance Plot 
#4a. Prepare a data frame with just the balance we need
# Relabel IPW weights, and extract
IPW.balance.data <- IPW.balance$data 
levels(IPW.balance.data$contrast) <- c("IPTW", "Unweighted")
IPW.balance.data.wt <- IPW.balance.data %>% filter(contrast =="IPTW")
# relabel ATE weights
ATE.balance.data <- ATE.balance$data 
levels(ATE.balance.data$contrast) <- c("Balancing Weights", "Unweighted")
# Add IPW weights to ATE weights
ATE.bal.overall <- rbind(ATE.balance.data, IPW.balance.data.wt)
ATE.bal.overall

# Change covariate names
levels(ATE.bal.overall$covariate) <- c("Admission Type Surgery", "Age", "Heart Failure",
                                       "CKD", "Diabetes Mellitus - Complicated",
                                       "Diabetes Mellitus - None", "Diabetes Mellitus - Uncomplicated",
                                       "ICU - Current", "eGFR", "Loop Diuretic", "Admitted to ED", "Admitted to Medical Floor",
                                       "Admitted to ICU", "Admitted to OR", "Admitted to Other", 
                                       "Prior Length of Stay", "Sex", "Vancomycin")

#4b. Plot
plot <- ggplot(data = ATE.bal.overall, aes(x = std.dif, y = covariate, 
                                     shape = factor(contrast) 
                                     #, color = factor(contrast) # for color instead of shapes, toggle as needed this,  geom_point, scale_color_manual
                                     )) +
  #geom_point(size = 2, shape = 16) +
  geom_point(size = 2, color = "black")+
  scale_shape_manual(name = "Contrast", values = c(5, 1, 15), 
                     breaks = c('Unweighted', "IPTW","Balancing Weights" )) +
  #scale_color_manual(name = "Contrast", values = c("darkred", "darkblue", "darkgreen")) +
  xlab("Standardized Difference") + ylab("Covariates") +
  ggtitle("Balance in Weighted Population") +
  scale_y_discrete(limits = rev(levels(ATE.bal.overall$covariate))) +
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

ess.tab <- round(rbind(original.n, IPW.ess, ATE.ess))
rownames(ess.tab) <- c("Unweighted", "IPTW", "Balancing Weights")
colnames(ess.tab) <- c("Treatment", "Control", "Total")
ess.tab

write.csv(ess.tab, "/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/results/ess.table.csv", row.names = T)



#6. Estimate ATT
#6a.Poisson model with IPTWs 

#6b.Poisson model with Balancing weights for each value of lambda  

