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

# estimate hyperparameter
data.ctrl <- data %>% filter(pain==0)
lambda.reg <- lm(reformulate(covs, response = "kEver"), data=data.ctrl)
var(lambda.reg$resid)

#3b. Estimate Balancing weights for lambda = 0, 0.5, 1, 2 (maybe no?)
# Balance Weights ATE
out.pain.q <- multilevel_ate_qp(X, trt, rep(1,n), lambda = 0.05, lowlim = 0, uplim = 1,  
                                verbose= TRUE, exact_global = TRUE, scale_sample_size = FALSE)
# Process
data$BW <- out.pain.q$weights
summary(data$BW)
sd(data$BW)

#3c. Evaluate Balance
covs <- c("age", "sex", "admType", 
          "presentation.ed", "presentation.icu",
          "presentation.or", "presentation.floor", "presentation.other", "priorLos", "icuCurrent",
          "chf",  "dm.no", "dm.noncomp", "dm.comp", "ckd", 
          "indexGFR", "loopBase",
          "vancoBase")

source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/balance-plot.R")
BW.balance <- bal.plots(data, "BW", 'pain', covs)

#3d. Calculate effective sample size for each value of lambda 
source("/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/functions/ess-function.R")
BW.ess <- ess(data, "pain", "BW")


#4. Create overall Balance Plot 
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

# sort ATE.bal.overall2 So it is reverse alphabetical by contrast
BW.bal.overall2 <- BW.bal.overall
# sort BW.bal.overall2 So it is reverse alphabetical by contrast
BW.bal.overall2 <- BW.bal.overall %>% 
  mutate(contrast = factor(contrast, levels = c("Balancing Weights", "Model-Based Weights", "Unweighted"))) %>%
  arrange(desc(contrast))

plot <- ggplot(data = BW.bal.overall2, 
               aes(x = std.dif, y = covariate, 
                  shape = factor(contrast), 
                  color = factor(contrast))) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.1) +
  geom_vline(xintercept = -0.1) +
  geom_point(size = 3, shape = 16, color = "black")+
  geom_point(size = 2, shape = 16) +
  scale_color_manual(name = "Contrast", values = c("#339933", "#FFCC00","#990000")) +
  xlab("Standardized Difference") + ylab("Covariates") +
  ggtitle("Balance in Weighted Population") +
  scale_y_discrete(limits = rev(levels(BW.bal.overall$covariate))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

plot

#5. Create Effective Sample size table
og.n.tx <- table(data$pain)[[2]]
og.n.cx<- table(data$pain)[[1]]
og.n.total <- og.n.cx+og.n.tx
original.n <-c(og.n.tx, og.n.cx, og.n.total)

ess.tab <- round(rbind(original.n, MW.ess, BW.ess))
rownames(ess.tab) <- c("Unweighted", "Model-Based Weights", "Balancing Weights")
colnames(ess.tab) <- c("Treatment", "Control", "Total")
ess.tab

write.csv(ess.tab, "/Users/haedi/Library/CloudStorage/Box-Box/Repos/Balwts/results/tweet.ess.table.csv", row.names = T)


############################################ Misc. ######
# Make plot comparing the box plots of ATEwts versus IPW
boxplot.data <- data %>%
  select(ATEwts, IPW)

bxplt <- ggplot(data = boxplot.data, aes(x = factor(1), y = ATEwts)) +
  geom_boxplot(aes(fill = "ATEwts"), alpha = 0.5) +
  geom_boxplot(aes(x = factor(2), y = IPW, fill = "IPW"), alpha = 0.5) +
  scale_fill_manual(values = c("ATEwts" = "#339933", "IPW" = "#FFCC00")) +
  xlab("Weighting Method") + ylab("Weight") +
  ggtitle("Comparison of ATEwts and IPW") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

bxplt

# make density plots of the same data
dens.plt <- ggplot(data = boxplot.data, aes(x = ATEwts, fill = "ATEwts")) +
  geom_density(alpha = 0.5) +
  geom_density(data = boxplot.data, aes(x = IPW, fill = "IPW"), alpha = 0.5) +
  scale_fill_manual(values = c("ATEwts" = "#339933", "IPW" = "#FFCC00")) +
  xlab("Weight") + ylab("Density") +
  ggtitle("Density Plot of ATEwts and IPW") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dens.plt


# estimate power of test for binary outcome, sample size of treated and control, and effect size
library(pwr)
#power for base ATEwts based on ESS
pwr.2p2n.test(n1 = 3903, n2 = 72467, h=0.05, sig.level = 0.05) 
# power for IPW based on ESS
pwr.2p2n.test(n1 = 4163, n2 = 70065, h=0.05, sig.level = 0.05)

