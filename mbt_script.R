#Herman Buschke Memory Binding Test - Czech Validation
##
##Authors: Ondrej Bezdicek, Jiri Motyl, Tomas Nikolai, Hana Horakova, Adela Fendrych Mazancova, Martin Vyhnalek, Jakub Hort, Robert Jech, Evzen Ruzicka, Herman Buschke
##Analysis by: Jiri Motyl, jiri.motyl@vfn.cz
##Corresponding Author: Ondrej Bezdicek, ondrej.bezdicek@vfn.cz
##Date of analysis: 2023-01-13
##R version 4.2.2 (2022-10-31 ucrt)


###Description of key variables

####Grouping Variables:
#####Disease (Main 4 groups: HC, PD-NC, PD-MCI, aMCI-AD)
#####Disease_NO (Main 4 groups coded cumerically: HC (01), PD-NC(02), PD-MCI(03), aMCI-AD(04))
#####Disease_amyloid (AD divided into 2 subgroups: HC, PD-NC, PD-MCI, aMCI-AD Aβ unconfirmed, aMCI-AD Aβ confirmed)
#####Disease_NO_amyloid (Main 4 groups coded cumerically: HC (01), PD-NC(02), PD-MCI(03), aMCI-AD Aβ unconfirmed (04), aMCI-AD Aβ confirmed(05))

####Demographic Variables:
#####Edu = Education
#####Age = Age
#####Sex = Sex

####Test Variables:
#####MMSE = MMSE
#####MBT_IR_CR_L1L2 = MBT Instant Recall List 1 + List 2
#####MBT_IR_TIP = MBT Instant Recall ,Total number of items cued recalled in the Paired condition on the MBT; 
#####MBT_IR_PIP = MBT Instant RecallThe num, ber of pairs cued recalled in the paired condition of MBT
#####MBT_IR_FR_30.120 = MBT Free Recall, Total number of Items recalled in the 2 minutes Free recall condition on the MBT
#####RAVLT15 = RAVLT Learning, Sum of scores ofrom trials 1-5
#####RAVLT30 = RAVLT Total Score in delayed recall after 30 minutes
#####TMTA = TMT A time
#####TMTB = TMT B time
#####StroopDots = Dots condition in Prague Stroop Test
#####StroopWords = Words condition in Prague Stroop Test
#####StroopColors = Colors condition in Prague Stroop Test



#necessary packages
install.packages("psych")
install.packages("tidyverse")
install.packages("rcompanion")
install.packages("PMCMRplus")
install.packages("WRS2")
install.packages("pROC")
install.packages("gt")
install.packages("janitor")
install.packages("papaja")
install.packages("jmv")



#necessary libraries
library(psych)
library(tidyverse)
library(rcompanion)
library(PMCMRplus)
library(WRS2)
library(pROC)
library(gt)
library(janitor)
library(papaja)
library(jmv)



#source csv file
MBT <- read_csv("mbt_221030.csv", 
                  na = "N/A")
spec(MBT)
View(MBT)

#variables
##variable types
MBT$Disease <- as.factor(MBT$Disease)
MBT$Disease_NO <- as.factor(MBT$Disease_NO)
MBT$Disease_NO <- factor(MBT$Disease_NO,
                         levels = c(1,2,3,4),
                         labels = c("HC", "PD-NC", "PD-MCI","AD-aMCI"))
MBT$Disease_NO_amyloid <- as.factor(MBT$Disease_NO_amyloid)
MBT$Disease_NO_amyloid <- factor(MBT$Disease_NO_amyloid,
                         levels = c(1,2,3,4, 5),
                         labels = c("HC", "PD-NC", "PD-MCI","AD-aMCI","AD-aMCI Aß conf."))

##new variables for the analyses - group comparisons
##aMCI-AD and aMCI-AD beta-amyloid confirmed as a single group
###PD-NC versus PD-MCI
MBT$Disease_NO_unclass <- unclass(MBT$Disease_NO)
MBT$PDNC_vs_PDMCI <- MBT$Disease_NO_unclass
MBT$PDNC_vs_PDMCI[MBT$PDNC_vs_PDMCI == 1] <- NA
MBT$PDNC_vs_PDMCI[MBT$PDNC_vs_PDMCI == 4] <- NA

###HC versus aMCI-AD
MBT$HC_vs_aMCIAD <- MBT$Disease_NO_unclass
MBT$HC_vs_aMCIAD[MBT$HC_vs_aMCIAD == 2] <- NA
MBT$HC_vs_aMCIAD[MBT$HC_vs_aMCIAD == 3] <- NA

###HC versus PD-MCI
MBT$HC_vs_PDMCI <- MBT$Disease_NO_unclass
MBT$HC_vs_PDMCI[MBT$HC_vs_PDMCI == 2] <- NA
MBT$HC_vs_PDMCI[MBT$HC_vs_PDMCI == 4] <- NA

###HC versus PD-NC
MBT$HC_vs_PDNC <- MBT$Disease_NO_unclass
MBT$HC_vs_PDNC[MBT$HC_vs_PDNC == 3] <- NA
MBT$HC_vs_PDNC[MBT$HC_vs_PDNC == 4] <- NA

###PD-MCI versus aMCI-AD
MBT$PDMCI_vs_aMCIAD <- MBT$Disease_NO_unclass
MBT$PDMCI_vs_aMCIAD[MBT$PDMCI_vs_aMCIAD == 1] <- NA
MBT$PDMCI_vs_aMCIAD[MBT$PDMCI_vs_aMCIAD == 2] <- NA

##new variables for the analyses - group comparisons
##aMCI-AD and aMCI-AD beta-amyloid confirmed as two separate groups
###aMCI-AD versus aMCI-AD amyloid confirmed
MBT$Disease_NO_amyl_unclass <- unclass(MBT$Disease_NO_amyloid)
MBT$aMCIAD_vs_aMCIADamyloid <- MBT$Disease_NO_amyl_unclass
MBT$aMCIAD_vs_aMCIADamyloid[MBT$aMCIAD_vs_aMCIADamyloid == 1] <- NA
MBT$aMCIAD_vs_aMCIADamyloid[MBT$aMCIAD_vs_aMCIADamyloid == 2] <- NA
MBT$aMCIAD_vs_aMCIADamyloid[MBT$aMCIAD_vs_aMCIADamyloid == 3] <- NA

###HC versus aMCI-AD amyloid confirmed
MBT$HC_vs_aMCIADamyloid <- MBT$Disease_NO_amyl_unclass
MBT$HC_vs_aMCIADamyloid[MBT$HC_vs_aMCIADamyloid == 2] <- NA
MBT$HC_vs_aMCIADamyloid[MBT$HC_vs_aMCIADamyloid == 3] <- NA
MBT$HC_vs_aMCIADamyloid[MBT$HC_vs_aMCIADamyloid == 4] <- NA

###PD-MCI versus aMCI-AD amyloid confirmed
MBT$PDMCI_vs_aMCIADamyloid <- MBT$Disease_NO_amyl_unclass
MBT$PDMCI_vs_aMCIADamyloid[MBT$PDMCI_vs_aMCIADamyloid == 1] <- NA
MBT$PDMCI_vs_aMCIADamyloid[MBT$PDMCI_vs_aMCIADamyloid == 2] <- NA
MBT$PDMCI_vs_aMCIADamyloid[MBT$PDMCI_vs_aMCIADamyloid == 4] <- NA

MBT$PDNC_vs_PDMCI <- as.factor(MBT$PDNC_vs_PDMCI)
MBT$HC_vs_aMCIAD <- as.factor(MBT$HC_vs_aMCIAD)
MBT$HC_vs_PDMCI <- as.factor(MBT$HC_vs_PDMCI)
MBT$HC_vs_PDNC <- as.factor(MBT$HC_vs_PDNC)
MBT$PDMCI_vs_aMCIAD <- as.factor(MBT$PDMCI_vs_aMCIAD)
MBT$aMCIAD_vs_aMCIADamyloid <- as.factor(MBT$aMCIAD_vs_aMCIADamyloid)
MBT$HC_vs_aMCIADamyloid <- as.factor(MBT$HC_vs_aMCIADamyloid)
MBT$PDMCI_vs_aMCIADamyloid <- as.factor(MBT$PDMCI_vs_aMCIADamyloid)

#Description of data
MBT_description <- describe(MBT, skew = T, IQR = T)
MBT_description_grouped <- describeBy(MBT, group = "Disease", skew = T, IQR = T)
MBT_description_grouped_amyl <- describeBy(MBT, group = "Disease_amyloid", skew = T, IQR = T)


##histograms
hist(MBT$MMSE)
hist(MBT$MBT_IR_CR_L1L2)
hist(MBT$MBT_IR_TIP)
hist(MBT$MBT_IR_PIP)
hist(MBT$MBT_IR_FR_30.120)
hist(MBT$RAVLT15)
hist(MBT$RAVLT30)
hist(MBT$TMTA)
hist(MBT$TMTB)
hist(MBT$StroopDots)
hist(MBT$StroopWords)
hist(MBT$StroopColors)

##QQ Plots
qqnorm(MBT$MMSE)
qqline(MBT$MMSE)

qqnorm(MBT$MBT_IR_CR_L1L2)
qqline(MBT$MBT_IR_CR_L1L2)

qqnorm(MBT$MBT_IR_TIP)
qqline(MBT$MBT_IR_TIP)

qqnorm(MBT$MBT_IR_PIP)
qqline(MBT$MBT_IR_PIP)

qqnorm(MBT$MBT_IR_FR_30.120)
qqline(MBT$MBT_IR_FR_30.120)

qqnorm(MBT$RAVLT15)
qqline(MBT$RAVLT15)

qqnorm(MBT$RAVLT30)
qqline(MBT$RAVLT30)

qqnorm(MBT$TMTA)
qqline(MBT$TMTA)

qqnorm(MBT$TMTB)
qqline(MBT$TMTB)

qqnorm(MBT$StroopDots)
qqline(MBT$StroopDots)

qqnorm(MBT$StroopWords)
qqline(MBT$StroopWords)

qqnorm(MBT$StroopColors)
qqline(MBT$StroopColors)


##transform skewed variables (skew >/< +1/-1)
###positive skew (log10)
MBT$TMTA_log10 <- log10(MBT$TMTA)
MBT$TMTB_log10 <- log10(MBT$TMTB)
MBT$StroopDots_log10 <- log10(MBT$StroopDots)
MBT$StroopWords_log10 <- log10(MBT$StroopWords)
MBT$StroopColors_log10 <- log10(MBT$StroopColors)
###negative skew - reflection and then (log10)
MBT$MMSE_reflectedlog10 <- log10(31 - MBT$MMSE)
MBT$MBT_IR_CR_L1L2_reflectedlog10 <- log10(33 - MBT$MBT_IR_CR_L1L2)
MBT$MBT_IR_TIP_reflectedlog10 <- log10(33 - MBT$MBT_IR_TIP)

hist(MBT$MMSE_reflectedlog10)
hist(MBT$MBT_IR_CR_L1L2_reflectedlog10)
hist(MBT$MBT_IR_TIP_reflectedlog10)

#violinplots

##complete AD
###MBT TIP
MBT_TIP_violinplot <- ggplot(MBT, aes(y=`MBT_IR_TIP`,x=`Disease_NO`, fill = `Disease_NO`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_TIP_violinplot+
  labs(y = "MBT TIP Total Score", x = NULL)

###MBT PIP
MBT_PIP_violinplot <- ggplot(MBT, aes(y=`MBT_IR_PIP`,x=`Disease_NO`, fill = `Disease_NO`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_PIP_violinplot+
  labs(y = "MBT PIP Total Score", x = NULL)

###MBT Free Recall
MBT_FR_violinplot <- ggplot(MBT, aes(y=`MBT_IR_FR_30.120`,x=`Disease_NO`, fill = `Disease_NO`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_FR_violinplot+
  labs(y = "MBT FR Total Score", x = NULL)

##beta-amyloid confirmed AD
###MBT TIP
MBT_TIP_violinplot <- ggplot(MBT, aes(y=`MBT_IR_TIP`,x=`Disease_NO_amyloid`, fill = `Disease_NO_amyloid`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_TIP_violinplot+
  labs(y = "MBT TIP Total Score", x = NULL)

###MBT PIP
MBT_PIP_violinplot <- ggplot(MBT, aes(y=`MBT_IR_PIP`,x=`Disease_NO_amyloid`, fill = `Disease_NO_amyloid`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_PIP_violinplot+
  labs(y = "MBT PIP Total Score", x = NULL)

###MBT Free Recall
MBT_FR_violinplot <- ggplot(MBT, aes(y=`MBT_IR_FR_30.120`,x=`Disease_NO_amyloid`, fill = `Disease_NO_amyloid`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_FR_violinplot+
  labs(y = "MBT FR Total Score", x = NULL)


#Basic statistics
##Basic display functions
###max 2 decimals
fn_2decimals <- function(x = Number) {
  papaja::printnum(x, digits = 2, big.interval=10)
}

###max 3 decimals
fn_3decimals <- function(x = Number) {
  papaja::printnum(x, gt1 = TRUE, digits = 3, big.interval=10)
}

###max 3 decimals (drop zero at the beginning)
fn_3decimals_no0 <- function(x = Number) {
  papaja::printnum(x, gt1 = FALSE, digits = 3, big.interval=10)
}

##Demographics - compare means / group comparisons
###function: Kruskal Wallis 
function_kw <- function(Variable, Grouping, Dataset) {
  fun01 <- kruskal.test(Variable ~ Grouping, data = Dataset)
  fun02 <- epsilonSquared(x = Variable, g = Grouping)
  fun03 <- dscfAllPairsTest(Variable ~ Grouping, data = Dataset)
  print(fun01)
  print(fun02)
  print(fun03)
}

fn_kw_basic <- function(Variable, Grouping, Dataset) {
  fun01 <- kruskal.test(Variable ~ Grouping, data = Dataset)
  fun02 <- epsilonSquared(x = Variable, g = Grouping)
  fun03 <- dscfAllPairsTest(Variable ~ Grouping, data = Dataset)
  kw_results <- c(fun01$p.value, unname(fun02[]))
  kw_results_p <- papaja::printnum(fun01$p.value, gt1 = TRUE, digits = 3)
  kw_results_ES <- papaja::printnum(unname(fun02[]), gt1 = TRUE, digits = 3)
  kw_results_final <- papaja::printnum(kw_results, gt1 = TRUE, digits = 3)
  return(kw_results_final)
}

fn_kw_p <- function(Variable, Grouping, Dataset) {
  fun01 <- kruskal.test(Variable ~ Grouping, data = Dataset)
  kw_results_p <- papaja::printnum(fun01$p.value, gt1 = TRUE, digits = 3)
  kw_results_p_updated <- ifelse(fun01$p.value < .001, "<0.001", kw_results_p)
  return(kw_results_p_updated)
}

fn_kw_p_stars <- function(Variable, Grouping, Dataset) {
  fun01 <- kruskal.test(Variable ~ Grouping, data = Dataset)
  mystars <- ifelse(fun01$p.value < .001, "***"
                    , ifelse(fun01$p.value < .01, "**"
                             , ifelse(fun01$p.value < .05, "*"
                                      , " ")))
  return(mystars)
}

fn_kw_es <- function(Variable, Grouping, Dataset) {
  fun02 <- epsilonSquared(x = Variable, g = Grouping)
  kw_results_es <- papaja::printnum(unname(fun02[]), gt1 = TRUE, digits = 3)
  return(kw_results_es)
}

###function: Welch's t-test 
fn_welch_t <- function(Variable, Grouping, Dataset) {
  fun01 <- t.test(Variable ~ Grouping, data = Dataset)
  return(fun01)
}

fn_wt_t <- function(Variable, Grouping, Dataset) {
  fun01 <- t.test(Variable ~ Grouping, data = Dataset)
  wt_results_t <- papaja::printnum(fun01[["statistic"]][["t"]], gt1 = TRUE, digits = 3)
  return(wt_results_t)
}

fn_wt_p <- function(Variable, Grouping, Dataset) {
  fun01 <- t.test(Variable ~ Grouping, data = Dataset)
  wt_results_p <- papaja::printnum(fun01$p.value, gt1 = TRUE, digits = 3)
  return(wt_results_p)
}

fn_wt_p_stars <- function(Variable, Grouping, Dataset) {
  fun01 <- t.test(Variable ~ Grouping, data = Dataset)
  mystars <- ifelse(fun01$p.value < .001, "***"
                    , ifelse(fun01$p.value < .01, "**"
                             , ifelse(fun01$p.value < .05, "*"
                                      , " ")))
  return(mystars)
}

fn_wt_t_stars <- function(Variable, Grouping, Dataset) {
  fun01 <- t.test(Variable ~ Grouping, data = Dataset)
  wt_results_t <- papaja::printnum(fun01[["statistic"]][["t"]], gt1 = TRUE, digits = 3)
  mystars <- ifelse(fun01$p.value < .001, "***"
                    , ifelse(fun01$p.value < .01, "**"
                             , ifelse(fun01$p.value < .05, "*"
                                      , " ")))
  return(paste(wt_results_t, mystars, sep=""))
}


###function: Chi-Squared
function_chi <- function(Variable, Grouping) {
  fun04 <- chisq.test(Variable, y = Grouping)
  fun05 <- cramerV(Variable, y = Grouping)
  print(fun04)
  print(fun05)
}

fn_chi_basic <- function(Variable, Grouping) {
  fun04 <- chisq.test(Variable, y = Grouping)
  fun05 <- cramerV(Variable, y = Grouping)
  chi_results <- c(fun04$p.value, unname(fun05))
  chi_results_final <- papaja::printnum(chi_results, gt1 = TRUE, digits = 3)
  return(chi_results_final)
}

fn_chi_p <- function(Variable, Grouping) {
  fun04 <- chisq.test(Variable, y = Grouping)
  chi_results_p <- papaja::printnum(fun04$p.value, gt1 = TRUE, digits = 3)
  chi_results_p_updated <- ifelse(fun04$p.value < .001, "<0.001", chi_results_p)
  return(chi_results_p_updated)
}

fn_chi_es <- function(Variable, Grouping) {
  fun05 <- cramerV(Variable, y = Grouping)
  chi_results_es <- papaja::printnum(unname(fun05), gt1 = TRUE, digits = 3)
  return(chi_results_es)
}

###Age
function_kw(MBT$Age, MBT$Disease_NO, MBT)
function_kw(MBT$Age, MBT$Disease_NO_amyloid, MBT)
posthoc_age <- dscfAllPairsTest(MBT$Age ~ MBT$Disease_NO, data = MBT)
####Single test for differences between bot AD-aMCI subgroups
fn_welch_t(MBT$Age, MBT$aMCIAD_vs_aMCIADamyloid, MBT)

###Education
function_kw(MBT$Edu, MBT$Disease_NO, MBT)
function_kw(MBT$Edu, MBT$Disease_NO_amyloid, MBT)
posthoc_edu <- dscfAllPairsTest(MBT$Edu ~ MBT$Disease_NO, data = MBT)
####Single test for differences between bot AD-aMCI subgroups
fn_welch_t(MBT$Edu, MBT$aMCIAD_vs_aMCIADamyloid, MBT)

###Sex
function_chi(MBT$Sex, MBT$Disease_NO)
table(MBT$Sex, MBT$Disease_NO)
function_chi(MBT$Sex, MBT$Disease_NO_amyloid)
table(MBT$Sex, MBT$Disease_NO_amyloid)


####Computes number of Females/Males in each group
function_sex_count <- function(data, patientgroup) {
  MBT_Disease <- filter(data, Disease_NO == patientgroup)
  no_females <- length(MBT_Disease$Sex[MBT_Disease$Sex=="female"])
  no_males <- length(MBT_Disease$Sex[MBT_Disease$Sex=="male"])
  paste(no_females, no_males, sep="/")
}

####Computes number of Females/Males in each group (AD devided)
function_sex_count_amyl <- function(data, patientgroup) {
  MBT_Disease <- filter(data, Disease_NO_amyloid == patientgroup)
  no_females <- length(MBT_Disease$Sex[MBT_Disease$Sex=="female"])
  no_males <- length(MBT_Disease$Sex[MBT_Disease$Sex=="male"])
  paste(no_females, no_males, sep="/")
}

####Sex - participants in each group
MBT_description_grouped$HC["ID*","n"]
function_sex_count(MBT, "HC")
MBT_description_grouped$"PD-NC"["ID*","n"]
function_sex_count(MBT, "PD-NC")
MBT_description_grouped$"PD-MCI"["ID*","n"]
function_sex_count(MBT, "PD-MCI")
MBT_description_grouped$"aMCI AD"["ID*","n"]
function_sex_count(MBT, "AD-aMCI")
MBT_description_grouped_amyl$"aMCI AD"["ID*","n"]
function_sex_count_amyl(MBT, "AD-aMCI")
MBT_description_grouped_amyl$"aMCI AD amyloid"["ID*","n"]
function_sex_count_amyl(MBT, "AD-aMCI Aß conf.")


###Table 01 Demographics:  Number of participants, means and SDs for Age and Education, and relevant tests (Kruskal-Walis, Chi-Squared, etc.). Divided by patient groups
rownames(MBT_description_grouped$HC)

fn_demographics_means <- function(Variable = "Age*", Group = MBT_description_grouped$"PD-NC"){
  paste(fn_2decimals(Group[Variable,"mean"])," ","(\u00b1",fn_2decimals(Group[Variable,"sd"]),")", sep = "")
}

fn_demographics_n <- function(Group = MBT_description_grouped$HC){
  as.character(Group["ID*","n"])
}

mbt_tab01 <- tribble(
  ~variable, ~HC_M, ~PDNC_M,  ~PDMCI_M, ~ADaMCI_ALL_M, ~p_4maingroups, ~ES_4maingroups, ~ADaMCI_unconf_M, ~ADaMCI_confamyl_M, ~p_5groups_amyl, ~ES_5groups_amyl,
  "n (All)", fn_demographics_n(MBT_description_grouped$HC), fn_demographics_n(MBT_description_grouped$"PD-NC"), fn_demographics_n(MBT_description_grouped$"PD-MCI"), fn_demographics_n(MBT_description_grouped$"aMCI AD"), "", "", fn_demographics_n(MBT_description_grouped_amyl$"aMCI AD"), fn_demographics_n(MBT_description_grouped_amyl$"aMCI AD amyloid"), "", "",
  "n (Females/Males)",function_sex_count(MBT, "HC"), function_sex_count(MBT, "PD-NC"), function_sex_count(MBT, "PD-MCI"), function_sex_count(MBT, "AD-aMCI"), fn_chi_p(MBT$Sex, MBT$Disease_NO), fn_chi_es(MBT$Sex, MBT$Disease_NO), function_sex_count_amyl(MBT, "AD-aMCI"), function_sex_count_amyl(MBT, "AD-aMCI Aß conf."), fn_chi_p(MBT$Sex, MBT$Disease_NO_amyloid), fn_chi_es(MBT$Sex, MBT$Disease_NO_amyloid),
  "Age", fn_demographics_means("Age*", MBT_description_grouped$HC),fn_demographics_means("Age*", MBT_description_grouped$"PD-NC"),fn_demographics_means("Age*", MBT_description_grouped$"PD-MCI"),fn_demographics_means("Age*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$Age,MBT$Disease_NO,MBT), fn_kw_es(MBT$Age,MBT$Disease_NO,MBT), fn_demographics_means("Age*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("Age*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$Age, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$Age, MBT$Disease_NO_amyloid, MBT),
  "Edu", fn_demographics_means("Edu*", MBT_description_grouped$HC),fn_demographics_means("Edu*", MBT_description_grouped$"PD-NC"),fn_demographics_means("Edu*", MBT_description_grouped$"PD-MCI"),fn_demographics_means("Edu*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$Edu,MBT$Disease_NO,MBT), fn_kw_es(MBT$Edu,MBT$Disease_NO,MBT), fn_demographics_means("Edu*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("Edu*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$Edu, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$Edu, MBT$Disease_NO_amyloid, MBT),
)

write.csv(mbt_tab01, "mbt_tab01.csv", na = "NA")


##Test scores of healthy controls and patient groups - compare means
function_kw(MBT$MMSE, MBT$Disease_NO, MBT)
function_kw(MBT$TMTA, MBT$Disease_NO, MBT)
function_kw(MBT$TMTB, MBT$Disease_NO, MBT)
function_kw(MBT$StroopDots, MBT$Disease_NO, MBT)
function_kw(MBT$StroopWords, MBT$Disease_NO, MBT)
function_kw(MBT$StroopColors, MBT$Disease_NO, MBT)
function_kw(MBT$RAVLT15, MBT$Disease_NO, MBT)
function_kw(MBT$RAVLT30, MBT$Disease_NO, MBT)
function_kw(MBT$MBT_IR_CR_L1,MBT$Disease_NO, MBT)
function_kw(MBT$MBT_IR_CR_L2, MBT$Disease_NO, MBT)
function_kw(MBT$MBT_IR_CR_L1L2, MBT$Disease_NO, MBT)
function_kw(MBT$MBT_IR_TIP, MBT$Disease_NO, MBT)
function_kw(MBT$MBT_IR_PIP, MBT$Disease_NO, MBT)
function_kw(MBT$MBT_IR_FR_30.120, MBT$Disease_NO, MBT)

posthoc_MBT_IR_CR_L1L2 <- function_kw(MBT$MBT_IR_CR_L1L2, MBT$Disease_NO, MBT)
posthoc_MBT_IR_TIP <- function_kw(MBT$MBT_IR_TIP, MBT$Disease_NO, MBT)
posthoc_MBT_IR_PIP <- function_kw(MBT$MBT_IR_PIP, MBT$Disease_NO, MBT)

###Table 02 Test scores of healthy controls and patients’ groups
mbt_tab02 <- tribble(
  ~variable, ~HC_M, ~PDNC_M,  ~PDMCI_M, ~ADaMCI_ALL_M, ~p_4maingroups, ~ES_4maingroups, ~ADaMCI_unconf_M, ~ADaMCI_confamyl_M, ~p_5groups_amyl, ~ES_5groups_amyl,
  "MMSE", fn_demographics_means("MMSE*", MBT_description_grouped$HC), fn_demographics_means("MMSE*", MBT_description_grouped$"PD-NC"),fn_demographics_means("MMSE*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("MMSE*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$MMSE,MBT$Disease_NO,MBT), fn_kw_es(MBT$MMSE,MBT$Disease_NO,MBT), fn_demographics_means("MMSE*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("MMSE*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$MMSE, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$MMSE, MBT$Disease_NO_amyloid, MBT),
  "TMT-A", fn_demographics_means("TMTA*", MBT_description_grouped$HC), fn_demographics_means("TMTA*", MBT_description_grouped$"PD-NC"),fn_demographics_means("TMTA*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("TMTA*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$TMTA,MBT$Disease_NO,MBT), fn_kw_es(MBT$TMTA,MBT$Disease_NO,MBT), fn_demographics_means("TMTA*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("TMTA*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$TMTA, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$TMTA, MBT$Disease_NO_amyloid, MBT),
  "TMT-B", fn_demographics_means("TMTB*", MBT_description_grouped$HC), fn_demographics_means("TMTB*", MBT_description_grouped$"PD-NC"),fn_demographics_means("TMTB*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("TMTB*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$TMTB,MBT$Disease_NO,MBT), fn_kw_es(MBT$TMTB,MBT$Disease_NO,MBT), fn_demographics_means("TMTB*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("TMTB*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$TMTB, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$TMTB, MBT$Disease_NO_amyloid, MBT),
  "PST-D", fn_demographics_means("StroopDots*", MBT_description_grouped$HC), fn_demographics_means("StroopDots*", MBT_description_grouped$"PD-NC"),fn_demographics_means("StroopDots*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("StroopDots*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$StroopDots,MBT$Disease_NO,MBT), fn_kw_es(MBT$StroopDots,MBT$Disease_NO,MBT), fn_demographics_means("StroopDots*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("StroopDots*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$StroopDots, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$StroopDots, MBT$Disease_NO_amyloid, MBT),
  "PST-W", fn_demographics_means("StroopWords*", MBT_description_grouped$HC), fn_demographics_means("StroopWords*", MBT_description_grouped$"PD-NC"),fn_demographics_means("StroopWords*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("StroopWords*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$StroopWords,MBT$Disease_NO,MBT), fn_kw_es(MBT$StroopWords,MBT$Disease_NO,MBT), fn_demographics_means("StroopWords*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("StroopWords*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$StroopWords, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$StroopWords, MBT$Disease_NO_amyloid, MBT),
  "PST-C", fn_demographics_means("StroopColors*", MBT_description_grouped$HC), fn_demographics_means("StroopColors*", MBT_description_grouped$"PD-NC"),fn_demographics_means("StroopColors*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("StroopColors*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$StroopColors,MBT$Disease_NO,MBT), fn_kw_es(MBT$StroopColors,MBT$Disease_NO,MBT), fn_demographics_means("StroopColors*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("StroopColors*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$StroopColors, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$StroopColors, MBT$Disease_NO_amyloid, MBT),
  "RAVLT T1-5", fn_demographics_means("RAVLT15*", MBT_description_grouped$HC), fn_demographics_means("RAVLT15*", MBT_description_grouped$"PD-NC"),fn_demographics_means("RAVLT15*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("RAVLT15*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$RAVLT15,MBT$Disease_NO,MBT), fn_kw_es(MBT$RAVLT15,MBT$Disease_NO,MBT), fn_demographics_means("RAVLT15*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("RAVLT15*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$RAVLT15, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$RAVLT15, MBT$Disease_NO_amyloid, MBT),
  "RAVLT-30 min", fn_demographics_means("RAVLT30*", MBT_description_grouped$HC), fn_demographics_means("RAVLT30*", MBT_description_grouped$"PD-NC"),fn_demographics_means("RAVLT30*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("RAVLT30*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$RAVLT30,MBT$Disease_NO,MBT), fn_kw_es(MBT$RAVLT30,MBT$Disease_NO,MBT), fn_demographics_means("RAVLT30*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("RAVLT30*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$RAVLT30, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$RAVLT30, MBT$Disease_NO_amyloid, MBT),
  "MBT CR L1", fn_demographics_means("MBT_IR_CR_L1*", MBT_description_grouped$HC), fn_demographics_means("MBT_IR_CR_L1*", MBT_description_grouped$"PD-NC"),fn_demographics_means("MBT_IR_CR_L1*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("MBT_IR_CR_L1*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$MBT_IR_CR_L1,MBT$Disease_NO,MBT), fn_kw_es(MBT$MBT_IR_CR_L1,MBT$Disease_NO,MBT), fn_demographics_means("MBT_IR_CR_L1*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("MBT_IR_CR_L1*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$MBT_IR_CR_L1, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$MBT_IR_CR_L1, MBT$Disease_NO_amyloid, MBT),
  "MBT CR L2", fn_demographics_means("MBT_IR_CR_L2*", MBT_description_grouped$HC), fn_demographics_means("MBT_IR_CR_L2*", MBT_description_grouped$"PD-NC"),fn_demographics_means("MBT_IR_CR_L2*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("MBT_IR_CR_L2*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$MBT_IR_CR_L2,MBT$Disease_NO,MBT), fn_kw_es(MBT$MBT_IR_CR_L2,MBT$Disease_NO,MBT), fn_demographics_means("MBT_IR_CR_L2*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("MBT_IR_CR_L2*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$MBT_IR_CR_L2, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$MBT_IR_CR_L2, MBT$Disease_NO_amyloid, MBT),
  "MBT CR L1+L2", fn_demographics_means("MBT_IR_CR_L1L2*", MBT_description_grouped$HC), fn_demographics_means("MBT_IR_CR_L1L2*", MBT_description_grouped$"PD-NC"),fn_demographics_means("MBT_IR_CR_L1L2*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("MBT_IR_CR_L1L2*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$MBT_IR_CR_L1L2,MBT$Disease_NO,MBT), fn_kw_es(MBT$MBT_IR_CR_L1L2,MBT$Disease_NO,MBT), fn_demographics_means("MBT_IR_CR_L1L2*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("MBT_IR_CR_L1L2*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$MBT_IR_CR_L1L2, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$MBT_IR_CR_L1L2, MBT$Disease_NO_amyloid, MBT),
  "MBT-TIP", fn_demographics_means("MBT_IR_TIP*", MBT_description_grouped$HC), fn_demographics_means("MBT_IR_TIP*", MBT_description_grouped$"PD-NC"),fn_demographics_means("MBT_IR_TIP*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("MBT_IR_TIP*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$MBT_IR_TIP,MBT$Disease_NO,MBT), fn_kw_es(MBT$MBT_IR_TIP,MBT$Disease_NO,MBT), fn_demographics_means("MBT_IR_TIP*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("MBT_IR_TIP*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$MBT_IR_TIP, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$MBT_IR_TIP, MBT$Disease_NO_amyloid, MBT),
  "MBT-PIP", fn_demographics_means("MBT_IR_PIP*", MBT_description_grouped$HC), fn_demographics_means("MBT_IR_PIP*", MBT_description_grouped$"PD-NC"),fn_demographics_means("MBT_IR_PIP*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("MBT_IR_PIP*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$MBT_IR_PIP,MBT$Disease_NO,MBT), fn_kw_es(MBT$MBT_IR_PIP,MBT$Disease_NO,MBT), fn_demographics_means("MBT_IR_PIP*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("MBT_IR_PIP*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$MBT_IR_PIP, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$MBT_IR_PIP, MBT$Disease_NO_amyloid, MBT),
  "MBT-FR 2min", fn_demographics_means("MBT_IR_FR_30.120*", MBT_description_grouped$HC), fn_demographics_means("MBT_IR_FR_30.120*", MBT_description_grouped$"PD-NC"),fn_demographics_means("MBT_IR_FR_30.120*", MBT_description_grouped$"PD-MCI"), fn_demographics_means("MBT_IR_FR_30.120*", MBT_description_grouped$"aMCI AD"), fn_kw_p(MBT$MBT_IR_FR_30.120,MBT$Disease_NO,MBT), fn_kw_es(MBT$MBT_IR_FR_30.120,MBT$Disease_NO,MBT), fn_demographics_means("MBT_IR_FR_30.120*", MBT_description_grouped_amyl$"aMCI AD"),fn_demographics_means("MBT_IR_FR_30.120*", MBT_description_grouped_amyl$"aMCI AD amyloid"), fn_kw_p(MBT$MBT_IR_FR_30.120, MBT$Disease_NO_amyloid, MBT), fn_kw_es(MBT$MBT_IR_FR_30.120, MBT$Disease_NO_amyloid, MBT)
)


##function: correlation table (thanks Igluca from Bavaria!!  [https://github.com/crsh/papaja/issues/210])
correlation_tab <- function(x, y, z, export=FALSE) { #x= corr.matrix, y=use, z=method
  
  r <-corr.test(x, use = y, method = z)$r	#taking just the correlation matrix; no N, or p
  p <-corr.test(x, use = y, method = z)$p	#taking the p*s
  
  #define notions for significance levels
  mystars <- ifelse(p < .001, "***"
                    , ifelse(p < .01, "**"
                             , ifelse(p < .05, "*"
                                      , " ")))
  
  #round r, define new matrix Rnew with the correlations from rnd and paste mystars
  rnd  <- papaja::printnum(r, gt1 = TRUE, digits = 3)  #round, drop leading 0 - Thanks CRSH!								                     
  Rnew <- matrix(paste(rnd, mystars, sep=""), ncol=ncol(rnd)) 
  
  #remove 1.0 correlations from diagonal  and set the strings
  diag(Rnew) <- ''		
  Rnew[upper.tri(Rnew)] <- ''								                	
  
  rownames(Rnew) <- paste(1:ncol(rnd), colnames(rnd), sep=" ")         #define number and name
  colnames(Rnew) <- paste(1:ncol(rnd), "", sep="") 			       #define number
  
  #fun-part: we trim the top half 
  Rnew[upper.tri(Rnew)] <- ''			
  Rnew
  
  Rnew <- cbind(round(describe(x)[,3:4],2), Rnew)		     #describe x, M sD - put them in the matrix
  colnames(Rnew)[1:2] <- c("M","SD")					      		#Beschriftung der neuen Spalten
  Rnew <- Rnew[,1:(ncol(Rnew)-1)]							        	#delete the last column (ugly)
  
  #export to clipboard
  
  if (export==TRUE){
    result<-write.table(Rnew
                        , "clipboard"
                        , sep=";"
                        , row.names=FALSE)
  }
  else result <- Rnew
  return(result)
  
}

##Correlation between MBT measures and demographics (Age, Years of Education)
MBT_Tab3_corr <- select(MBT, 
                        Age,
                        Edu,
                        MBT_IR_CR_L1,
                        MBT_IR_CR_L2,
                        MBT_IR_CR_L1L2,
                        MBT_IR_TIP,
                        MBT_IR_PIP,
                        MBT_IR_FR_30.120)

correlation_tab(MBT_Tab3_corr, "pairwise", "kendall", export = TRUE)
MBT_tab3_corr_results_csv <- correlation_tab(MBT_Tab3_corr, "pairwise", "kendall", export = FALSE)
rownames(MBT_tab3_corr_results_csv)   <- c("Age", "Education", "MBT IR-CR List 1", "MBT IR-CR List 2", "MBT IR-CR List 1+2", "MBT IR TIP", "MBT IR PIP", "MBT IR FR 2min")
colnames(MBT_tab3_corr_results_csv)   <- c("M", "SD","Age", "Education", "MBT IR-CR List 1", "MBT IR-CR List 2", "MBT IR-CR List 1+2", "MBT IR TIP", "MBT IR PIP")

MBT_tab3_corr_results_csv02 <- select(MBT_tab3_corr_results_csv,
                                      Age,
                                      Education)
MBT_tab3_corr_results_csv02 <- filter(MBT_tab3_corr_results_csv02, 
                                      !row_number() %in% c(1, 2))
MBT_tab03 <- mutate(MBT_tab3_corr_results_csv02, 
                                      Sex=c(
                                        fn_wt_t_stars(MBT$MBT_IR_CR_L1, MBT$Sex, MBT),
                                        fn_wt_t_stars(MBT$MBT_IR_CR_L2, MBT$Sex, MBT),
                                        fn_wt_t_stars(MBT$MBT_IR_CR_L1L2, MBT$Sex, MBT),
                                        fn_wt_t_stars(MBT$MBT_IR_TIP, MBT$Sex, MBT),
                                        fn_wt_t_stars(MBT$MBT_IR_PIP, MBT$Sex, MBT),
                                        fn_wt_t_stars(MBT$MBT_IR_FR_30.120, MBT$Sex, MBT)
                                        )
                                      )

write.csv(MBT_tab03, "MBT_tab03.csv", na = "NA")

##Correlation between MBT measures and other neuropsychological methods
mbt_supp_tab2_corr <- select(MBT, 
                        MBT_IR_CR_L1L2,
                        MBT_IR_TIP,
                        MBT_IR_PIP,
                        MBT_IR_FR_30.120,
                        MMSE,
                        TMTA,
                        TMTB,
                        StroopDots,
                        StroopWords,
                        StroopColors,
                        RAVLT15,
                        RAVLT30
                        )

mbt_supp_tab2_corr_results_csv <- correlation_tab(mbt_supp_tab2_corr, "pairwise", "kendall", export = FALSE)
rownames(mbt_supp_tab2_corr_results_csv)   <- c("MBT CR L1+L2", "MBT-TIP", "MBT-PIP", "MBT-FR 2min", "MMSE", "TMT-A", "TMT-B", "PST-D", "PST-W", "PST-C", "RAVLT 1-5", "RAVLT DR 30 min")
colnames(mbt_supp_tab2_corr_results_csv)   <- c("M", "SD","MBT CR L1+L2", "MBT-TIP", "MBT-PIP", "MBT-FR 2min", "MMSE", "TMT-A", "TMT-B", "PST-D", "PST-W", "PST-C", "RAVLT 1-5")

mbt_supp_tab2_corr_results_csv02 <- select(mbt_supp_tab2_corr_results_csv,
                                      3,
                                      4,
                                      5,
                                      6)
mbt_supp_tab2_corr_results_csv02 <- filter(mbt_supp_tab2_corr_results_csv02, 
                                      !row_number() %in% c(1))

write.csv(mbt_supp_tab2_corr_results_csv02, "MBT_supp_tab2.csv", na = "NA")


#ancova
set.seed(1234)
##four main groups
mbt_ancova_disease <- jmv::ancova(data = MBT,
            formula = MBT_IR_TIP ~ Disease_NO + Age + Edu,
            effectSize = c('eta', 'partEta'),
            homo = T,
            norm = T,
            qq = T,
            postHoc = ~ Disease_NO,
            postHocCorr = c('tukey','bonf'),
            emMeans = ~ Disease_NO,
            emmPlots = T,
            emmPlotData = T,
            emmPlotError = 'ci'
            )

mbt_ancova_disease$assump$homo
mbt_ancova_disease$assump$norm
mbt_ancova_disease$main
mbt_ancova_disease$postHoc
mbt_ancova_disease_df <- as.data.frame(mbt_ancova_disease$postHoc[[1]])
mbt_tab04a <- papaja::printnum(mbt_ancova_disease_df, gt1 = TRUE, digits = 3)
write.csv(mbt_tab04a, "MBT_tab04a.csv", na = "NA")

##AD divided into 2 subgroups (beta-amyloid confirmed/unconfirmed)
mbt_ancova_disease_amyl <- 
  jmv::ancova(data = MBT,
            formula = MBT_IR_TIP ~ Disease_NO_amyloid + Age + Edu,
            effectSize = c('eta', 'partEta'),
            homo = T,
            norm = T,
            qq = T,
            postHoc = ~ Disease_NO_amyloid ,
            postHocCorr = c('tukey','bonf'),
            emMeans = ~ Disease_NO_amyloid ,
            emmPlots = T,
            emmPlotData = T,
            emmPlotError = 'ci'
)

mbt_ancova_disease_amyl$assump$homo
mbt_ancova_disease_amyl$assump$norm
mbt_ancova_disease_amyl$main
mbt_ancova_disease_amyl$postHoc
mbt_ancova_disease_amyl_df <- as.data.frame(mbt_ancova_disease_amyl$postHoc[[1]])
mbt_tab04b <- papaja::printnum(mbt_ancova_disease_amyl_df, gt1 = TRUE, digits = 3)
write.csv(mbt_tab04b, "MBT_tab04b.csv", na = "NA")


#robust ancova
set.seed(1234)
##four main groups
MBT_robust_ancova_HCAD_age <- WRS2::ancova(MBT_IR_TIP ~ HC_vs_aMCIAD + Age, data=MBT, tr = 0.2)
MBT_robust_ancboot_HCAD_age <- WRS2::ancboot(MBT_IR_TIP ~ HC_vs_aMCIAD + Age, data=MBT, tr = 0.2, nboot=4000)

MBT_robust_ancova_HCPDNC_age <- WRS2::ancova(MBT_IR_TIP ~ HC_vs_PDNC + Age, data=MBT, tr = 0.2)
MBT_robust_ancboot_HCPDNC_age <- WRS2::ancboot(MBT_IR_TIP ~ HC_vs_PDNC + Age, data=MBT, tr = 0.2, nboot=4000)

MBT_robust_ancova_HCAD_age
MBT_robust_ancboot_HCAD_age


MBT_robust_ancova_HCPDNC_age
MBT_robust_ancboot_HCPDNC_age

##AD divided into 2 subgroups (beta-amyloid confirmed/unconfirmed)
MBT_robust_ancova_HCADamyl_age <- WRS2::ancova(MBT_IR_TIP ~ HC_vs_aMCIADamyloid + Age, data=MBT, tr = 0.2)

MBT_robust_ancova_HCADamyl_age


#ROC
ROC_MBT_PDAD <- roc(formula=MBT$PDMCI_vs_aMCIAD ~MBT$MBT_IR_CR_L1L2+MBT$MBT_IR_TIP+MBT$MBT_IR_PIP+MBT$RAVLT30+MBT$RAVLT15, levels = c(3,4), na.rm = TRUE, auc = TRUE, ci = TRUE, plot= TRUE)
ROC_MBT_PDADamyl <- roc(formula=MBT$PDMCI_vs_aMCIADamyloid ~MBT$MBT_IR_CR_L1L2+MBT$MBT_IR_TIP+MBT$MBT_IR_PIP+MBT$RAVLT30+MBT$RAVLT15, levels = c(3,5), na.rm = TRUE, auc = TRUE, ci = TRUE, plot= TRUE)
ROC_MBT_HCAD <- roc(formula=MBT$HC_vs_aMCIAD ~MBT$MBT_IR_CR_L1L2+MBT$MBT_IR_TIP+MBT$MBT_IR_PIP+MBT$RAVLT30+MBT$RAVLT15, levels = c(1,4), na.rm = TRUE, auc = TRUE, ci = TRUE, plot= TRUE)
ROC_MBT_HCADamyl <- roc(formula=MBT$HC_vs_aMCIADamyloid ~MBT$MBT_IR_CR_L1L2+MBT$MBT_IR_TIP+MBT$MBT_IR_PIP+MBT$RAVLT30+MBT$RAVLT15, levels = c(1,5), na.rm = TRUE, auc = TRUE, ci = TRUE, plot= TRUE)
ROC_MBT_HCPD <- roc(formula=MBT$HC_vs_PDMCI ~MBT$MBT_IR_CR_L1L2+MBT$MBT_IR_TIP+MBT$MBT_IR_PIP+MBT$RAVLT30+MBT$RAVLT15, levels = c(1,3), na.rm = TRUE, auc = TRUE, ci = TRUE, plot= TRUE)
ROC_MBT_PDNC_vs_PDMCI <- roc(formula=MBT$PDNC_vs_PDMCI~MBT$MBT_IR_CR_L1L2+MBT$MBT_IR_TIP+MBT$MBT_IR_PIP+MBT$RAVLT30+MBT$RAVLT15, levels = c(2,3), na.rm = TRUE, auc = TRUE, ci = TRUE, plot= TRUE)

##AUC CI
###HC versus AD-aMCI all
AUCCI_HCAD_MBT_TIP <- ci.auc(ROC_MBT_HCAD$'MBT$MBT_IR_TIP'$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCAD_MBT_PIP <- ci.auc(ROC_MBT_HCAD$'MBT$MBT_IR_PIP'$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCAD_RAVLT15 <- ci.auc(ROC_MBT_HCAD$'MBT$RAVLT15'$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCAD_RAVLT30 <- ci.auc(ROC_MBT_HCAD$'MBT$RAVLT30'$auc, conf.level=0.95, method=c("delong"))

ROC_MBT_HCAD$'MBT$MBT_IR_TIP'$auc
AUCCI_HCAD_MBT_TIP
AUCCI_HCAD_RAVLT15
ROC_MBT_HCAD$'MBT$RAVLT30'$auc
AUCCI_HCAD_RAVLT30


###HC versus AD-aMCI amyloid confirmed
AUCCI_HCADamyl_MBT_TIP <- ci.auc(ROC_MBT_HCADamyl$'MBT$MBT_IR_TIP'$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCADamyl_MBT_PIP <- ci.auc(ROC_MBT_HCADamyl$'MBT$MBT_IR_PIP'$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCADamyl_RAVLT15 <- ci.auc(ROC_MBT_HCADamyl$'MBT$RAVLT15'$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCADamyl_RAVLT30 <- ci.auc(ROC_MBT_HCADamyl$'MBT$RAVLT30'$auc, conf.level=0.95, method=c("delong"))

ROC_MBT_HCADamyl$'MBT$MBT_IR_TIP'$auc
AUCCI_HCADamyl_MBT_TIP
ROC_MBT_HCADamyl$'MBT$MBT_IR_PIP'$auc
AUCCI_HCADamyl_MBT_PIP
ROC_MBT_HCADamyl$'MBT$RAVLT15'$auc
AUCCI_HCADamyl_RAVLT15
ROC_MBT_HCADamyl$'MBT$RAVLT30'$auc
AUCCI_HCADamyl_RAVLT30


###PD-MCI versus AD-aMCI
AUCCI_PDAD_MBT_TIP <- ci.auc(ROC_MBT_PDAD$`MBT$MBT_IR_TIP`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDAD_MBT_PIP <- ci.auc(ROC_MBT_PDAD$`MBT$MBT_IR_PIP`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDAD_RAVLT15 <- ci.auc(ROC_MBT_PDAD$`MBT$RAVLT15`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDAD_RAVLT30 <- ci.auc(ROC_MBT_PDAD$`MBT$RAVLT30`$auc, conf.level=0.95, method=c("delong"))

ROC_MBT_PDAD$`MBT$MBT_IR_TIP`$auc
AUCCI_PDAD_MBT_TIP
ROC_MBT_PDAD$`MBT$MBT_IR_PIP`$auc
AUCCI_PDAD_MBT_PIP
ROC_MBT_PDAD$`MBT$RAVLT15`$auc
AUCCI_PDAD_RAVLT15
ROC_MBT_PDAD$`MBT$RAVLT30`$auc
AUCCI_PDAD_RAVLT30


###PD-MCI versus AD-aMCI beta-amyloid confirmed
AUCCI_PDADamyl_MBT_TIP <- ci.auc(ROC_MBT_PDADamyl$`MBT$MBT_IR_TIP`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDADamyl_MBT_PIP <- ci.auc(ROC_MBT_PDADamyl$`MBT$MBT_IR_PIP`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDADamyl_RAVLT15 <- ci.auc(ROC_MBT_PDADamyl$`MBT$RAVLT15`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDADamyl_RAVLT30 <- ci.auc(ROC_MBT_PDADamyl$`MBT$RAVLT30`$auc, conf.level=0.95, method=c("delong"))

ROC_MBT_PDADamyl$`MBT$MBT_IR_TIP`$auc
AUCCI_PDADamyl_MBT_TIP
ROC_MBT_PDADamyl$`MBT$MBT_IR_PIP`$auc
AUCCI_PDADamyl_MBT_PIP
ROC_MBT_PDADamyl$`MBT$RAVLT15`$auc
AUCCI_PDADamyl_RAVLT15
ROC_MBT_PDADamyl$`MBT$RAVLT30`$auc
AUCCI_PDADamyl_RAVLT30

###HC vs PD-MCI
AUCCI_HCPD_MBT_TIP <- ci.auc(ROC_MBT_HCPD$`MBT$MBT_IR_TIP`$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCPD_MBT_PIP <- ci.auc(ROC_MBT_HCPD$`MBT$MBT_IR_PIP`$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCPD_RAVLT15 <- ci.auc(ROC_MBT_HCPD$`MBT$RAVLT15`$auc, conf.level=0.95, method=c("delong"))
AUCCI_HCPD_RAVLT30 <- ci.auc(ROC_MBT_HCPD$`MBT$RAVLT30`$auc, conf.level=0.95, method=c("delong"))

ROC_MBT_HCPD$`MBT$MBT_IR_TIP`$auc
AUCCI_HCPD_MBT_TIP
ROC_MBT_HCPD$`MBT$MBT_IR_PIP`$auc
AUCCI_HCPD_MBT_PIP
ROC_MBT_HCPD$`MBT$RAVLT15`$auc
AUCCI_HCPD_RAVLT15
ROC_MBT_HCPD$`MBT$RAVLT30`$auc
AUCCI_HCPD_RAVLT30

###PD-NC vs PD-MCI
AUCCI_PDncPDmci_MBT_TIP <- ci.auc(ROC_MBT_PDNC_vs_PDMCI$`MBT$MBT_IR_TIP`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDncPDmci_MBT_PIP <- ci.auc(ROC_MBT_PDNC_vs_PDMCI$`MBT$MBT_IR_PIP`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDncPDmci_RAVLT15 <- ci.auc(ROC_MBT_PDNC_vs_PDMCI$`MBT$RAVLT15`$auc, conf.level=0.95, method=c("delong"))
AUCCI_PDncPDmci_RAVLT30 <- ci.auc(ROC_MBT_PDNC_vs_PDMCI$`MBT$RAVLT30`$auc, conf.level=0.95, method=c("delong"))

ROC_MBT_PDNC_vs_PDMCI$`MBT$MBT_IR_TIP`$auc
AUCCI_PDncPDmci_MBT_TIP
ROC_MBT_PDNC_vs_PDMCI$`MBT$MBT_IR_PIP`$auc
AUCCI_PDncPDmci_MBT_PIP
ROC_MBT_PDNC_vs_PDMCI$`MBT$RAVLT15`$auc
AUCCI_PDncPDmci_RAVLT15
ROC_MBT_PDNC_vs_PDMCI$`MBT$RAVLT30`$auc
AUCCI_PDncPDmci_RAVLT30

###AUC CI Table
fn_auc <- function(aucci_value) {
  fn01 <- fn_3decimals_no0(aucci_value[2])
  fn02 <- fn_3decimals_no0(aucci_value[1])
  fn03 <- fn_3decimals_no0(aucci_value[3])
  paste(fn01," (",fn02," - ",fn03,")", sep = "")
}

fn_auc(AUCCI_PDncPDmci_MBT_TIP)

mbt_tab05 <- tribble(
  ~variable, ~HC_vs_aMCIAD, ~HC_vs_aMCIADamyl, ~PDMCI_vs_aMCIAD, ~PDMCI_vs_aMCIADamyl, ~HC_vs_PDMCI, ~PDNC_vs_PDMCI,
  "MBT-TIP", fn_auc(AUCCI_HCAD_MBT_TIP), fn_auc(AUCCI_HCADamyl_MBT_TIP), fn_auc(AUCCI_PDAD_MBT_TIP), fn_auc(AUCCI_PDADamyl_MBT_TIP), fn_auc(AUCCI_HCPD_MBT_TIP), fn_auc(AUCCI_PDncPDmci_MBT_TIP),
  "MBT-PIP", fn_auc(AUCCI_HCAD_MBT_PIP), fn_auc(AUCCI_HCADamyl_MBT_PIP), fn_auc(AUCCI_PDAD_MBT_PIP), fn_auc(AUCCI_PDADamyl_MBT_PIP), fn_auc(AUCCI_HCPD_MBT_PIP), fn_auc(AUCCI_PDncPDmci_MBT_PIP),
  "RAVLT 1-5", fn_auc(AUCCI_HCAD_RAVLT15), fn_auc(AUCCI_HCADamyl_RAVLT15), fn_auc(AUCCI_PDAD_RAVLT15), fn_auc(AUCCI_PDADamyl_RAVLT15), fn_auc(AUCCI_HCPD_RAVLT15), fn_auc(AUCCI_PDncPDmci_RAVLT15),
  "RAVLT-DR 30min", fn_auc(AUCCI_HCAD_RAVLT30), fn_auc(AUCCI_HCADamyl_RAVLT30), fn_auc(AUCCI_PDAD_RAVLT30), fn_auc(AUCCI_PDADamyl_RAVLT30), fn_auc(AUCCI_HCPD_RAVLT30), fn_auc(AUCCI_PDncPDmci_RAVLT30),
)
write.csv(mbt_tab05, "MBT_tab05.csv", na = "NA")


##ROC coordinates
###HC versus AD-aMCI
mbt_supp_tab3a <- coords(ROC_MBT_HCAD$`MBT$MBT_IR_TIP`, input="threshold", ret=c("threshold",
                                                                     "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

mbt_supp_tab3a <- fn_3decimals(mbt_supp_tab3a)
write.csv(mbt_supp_tab3a, "MBT_supp_tab3a.csv", na = "NA")

mbt_supp_tab4a <- coords(ROC_MBT_HCAD$`MBT$MBT_IR_PIP`, input="threshold", ret=c("threshold",
                                                                     "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

mbt_supp_tab4a <- fn_3decimals(mbt_supp_tab4a)
write.csv(mbt_supp_tab4a, "MBT_supp_tab4a.csv", na = "NA")

coords(ROC_MBT_HCAD$`MBT$RAVLT15`, input="threshold", ret=c("threshold",
                                                            "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

coords(ROC_MBT_HCAD$`MBT$RAVLT30`, input="threshold", ret=c("threshold",
                                                            "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

###HC versus AD-aMCI amyloid confirmed
mbt_supp_tab3b <- coords(ROC_MBT_HCADamyl$`MBT$MBT_IR_TIP`, input="threshold", ret=c("threshold",
                                                               "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

mbt_supp_tab3b <- fn_3decimals(mbt_supp_tab3b)
write.csv(mbt_supp_tab3b, "MBT_supp_tab3b.csv", na = "NA")

mbt_supp_tab4b <- coords(ROC_MBT_HCADamyl$`MBT$MBT_IR_PIP`, input="threshold", ret=c("threshold",
                                                               "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

mbt_supp_tab4b <- fn_3decimals(mbt_supp_tab4b)
write.csv(mbt_supp_tab4b, "MBT_supp_tab4b.csv", na = "NA")

coords(ROC_MBT_HCADamyl$`MBT$RAVLT15`, input="threshold", ret=c("threshold",
                                                            "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

coords(ROC_MBT_HCADamyl$`MBT$RAVLT30`, input="threshold", ret=c("threshold",
                                                            "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

###PD-MCI versus AD-aMCI
coords(ROC_MBT_PDAD$`MBT$MBT_IR_TIP`, input="threshold", ret=c("threshold",
                                                                     "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

coords(ROC_MBT_PDAD$`MBT$MBT_IR_PIP`, input="threshold", ret=c("threshold",
                                                                     "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

coords(ROC_MBT_PDAD$`MBT$RAVLT15`, input="threshold", ret=c("threshold",
                                                            "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)

coords(ROC_MBT_PDAD$`MBT$RAVLT30`, input="threshold", ret=c("threshold",
                                                            "specificity", "sensitivity","youden"),
       drop=TRUE, best.method=c("youden"),
       best.weights=c(1, 0.5), transpose = FALSE, as.matrix=TRUE)


##ROC plots
###function to add legend outside of a plot
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

###PD-MCI versus AD-aMCI
par(oma=c(2, 0, 0, 0), xaxs="i")
roc(MBT$PDMCI_vs_aMCIAD, MBT$MBT_IR_TIP, plot=TRUE, percent=TRUE, col="#1c61b6", lwd=2, ci = T, print.auc=TRUE)
plot.roc(MBT$PDMCI_vs_aMCIAD, MBT$MBT_IR_PIP, percent=TRUE, col="#008600", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=45)
plot.roc(MBT$PDMCI_vs_aMCIAD, MBT$RAVLT15, percent=TRUE, col="pink", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=40)
plot.roc(MBT$PDMCI_vs_aMCIAD, MBT$RAVLT30, percent=TRUE, col="red", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=35)
title(main="AUC: PD-MCI vs. aMCI AD", line = 2.5)
add_legend("bottom", legend=c("MBT TIP", "MBT PIP", "RAVLT 1-5", "RAVLT Delayed Recall"), lwd=2, 
           col=c("#1c61b6", "#008600", "pink", "red"),
           horiz=T, bty='n', cex=0.8)

###PD-MCI versus AD-aMCI beta-amyloid confirmed
par(oma=c(2, 0, 0, 0), xaxs="i")
roc(MBT$PDMCI_vs_aMCIADamyloid, MBT$MBT_IR_TIP, plot=TRUE, percent=TRUE, col="#1c61b6", lwd=2, ci = T, print.auc=TRUE)
plot.roc(MBT$PDMCI_vs_aMCIADamyloid, MBT$MBT_IR_PIP, percent=TRUE, col="#008600", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=45)
plot.roc(MBT$PDMCI_vs_aMCIADamyloid, MBT$RAVLT15, percent=TRUE, col="pink", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=40)
plot.roc(MBT$PDMCI_vs_aMCIADamyloid, MBT$RAVLT30, percent=TRUE, col="red", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=35)
title(main="AUC: PD-MCI vs. aMCI AD Aß conf.", line = 2.5)
add_legend("bottom", legend=c("MBT TIP", "MBT PIP", "RAVLT 1-5", "RAVLT Delayed Recall"), lwd=2, 
           col=c("#1c61b6", "#008600", "pink", "red"),
           horiz=T, bty='n', cex=0.8)

###HC versus AD-aMCI
par(oma=c(2, 0, 0, 0), xaxs="i")
roc(MBT$HC_vs_aMCIAD, MBT$MBT_IR_TIP, plot=TRUE, percent=TRUE, col="#1c61b6", lwd=2, ci = T, print.auc=TRUE)
plot.roc(MBT$HC_vs_aMCIAD, MBT$MBT_IR_PIP, percent=TRUE, col="#008600", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=45)
plot.roc(MBT$HC_vs_aMCIAD, MBT$RAVLT15, percent=TRUE, col="pink", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=40)
plot.roc(MBT$HC_vs_aMCIAD, MBT$RAVLT30, percent=TRUE, col="red", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=35)
title(main="AUC: HC vs. aMCI AD", line = 2.5)
add_legend("bottom", legend=c("MBT TIP", "MBT PIP", "RAVLT 1-5", "RAVLT Delayed Recall"), lwd=2, 
           col=c("#1c61b6", "#008600", "pink", "red"),
           horiz=T, bty='n', cex=0.8)

###HC versus AD-aMCI beta-amyloid confirmed
par(oma=c(2, 0, 0, 0), xaxs="i")
roc(MBT$HC_vs_aMCIADamyloid, MBT$MBT_IR_TIP, plot=TRUE, percent=TRUE, col="#1c61b6", lwd=2, ci = T, print.auc=TRUE)
plot.roc(MBT$HC_vs_aMCIADamyloid, MBT$MBT_IR_PIP, percent=TRUE, col="#008600", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=45)
plot.roc(MBT$HC_vs_aMCIADamyloid, MBT$RAVLT15, percent=TRUE, col="pink", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=40)
plot.roc(MBT$HC_vs_aMCIADamyloid, MBT$RAVLT30, percent=TRUE, col="red", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=35)
title(main="AUC: HC vs. aMCI AD Aß conf.", line = 2.5)
add_legend("bottom", legend=c("MBT TIP", "MBT PIP", "RAVLT 1-5", "RAVLT Delayed Recall"), lwd=2, 
           col=c("#1c61b6", "#008600", "pink", "red"),
           horiz=T, bty='n', cex=0.8)

###HC versus PD-MCI
par(oma=c(2, 0, 0, 0), xaxs="i")
roc(MBT$HC_vs_PDMCI, MBT$MBT_IR_TIP, plot=TRUE, percent=TRUE, col="#1c61b6", lwd=2, ci = T, print.auc=TRUE)
plot.roc(MBT$HC_vs_PDMCI, MBT$MBT_IR_PIP, percent=TRUE, col="#008600", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=45)
plot.roc(MBT$HC_vs_PDMCI, MBT$RAVLT15, percent=TRUE, col="pink", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=40)
plot.roc(MBT$HC_vs_PDMCI, MBT$RAVLT30, percent=TRUE, col="red", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=35)
title(main="AUC: HC vs. PD-MCI", line = 2.5)
add_legend("bottom", legend=c("MBT TIP", "MBT PIP", "RAVLT 1-5", "RAVLT Delayed Recall"), lwd=2, 
           col=c("#1c61b6", "#008600", "pink", "red"),
           horiz=T, bty='n', cex=0.8)

###PD-NC versus PD-MCI
par(oma=c(2, 0, 0, 0), xaxs="i")
roc(MBT$PDNC_vs_PDMCI, MBT$MBT_IR_TIP, plot=TRUE, percent=TRUE, col="#1c61b6", lwd=2, ci = T, print.auc=TRUE)
plot.roc(MBT$PDNC_vs_PDMCI, MBT$MBT_IR_PIP, percent=TRUE, col="#008600", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=45)
plot.roc(MBT$PDNC_vs_PDMCI, MBT$RAVLT15, percent=TRUE, col="pink", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=40)
plot.roc(MBT$PDNC_vs_PDMCI, MBT$RAVLT30, percent=TRUE, col="red", lwd=2, ci = T, print.auc=TRUE, add=TRUE, print.auc.y=35)
title(main="AUC: PD-NC vs. PD-MCI", line = 2.5)
add_legend("bottom", legend=c("MBT TIP", "MBT PIP", "RAVLT 1-5", "RAVLT Delayed Recall"), lwd=2, 
           col=c("#1c61b6", "#008600", "pink", "red"),
           horiz=T, bty='n', cex=0.8)


#violinplots

##complete AD
###MBT TIP
MBT_TIP_violinplot <- ggplot(MBT, aes(y=`MBT_IR_TIP`,x=`Disease_NO`, fill = `Disease_NO`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_TIP_violinplot+
  labs(y = "MBT TIP Total Score", x = NULL)


###MBT PIP
MBT_PIP_violinplot <- ggplot(MBT, aes(y=`MBT_IR_PIP`,x=`Disease_NO`, fill = `Disease_NO`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_PIP_violinplot+
  labs(y = "MBT PIP Total Score", x = NULL)


###MBT Free Recall
MBT_FR_violinplot <- ggplot(MBT, aes(y=`MBT_IR_FR_30.120`,x=`Disease_NO`, fill = `Disease_NO`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_FR_violinplot+
  labs(y = "MBT FR Total Score", x = NULL)


##beta-amyloid confirmed AD
###MBT TIP
MBT_TIP_violinplot <- ggplot(MBT, aes(y=`MBT_IR_TIP`,x=`Disease_NO_amyloid`, fill = `Disease_NO_amyloid`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_TIP_violinplot+
  labs(y = "MBT TIP Total Score", x = NULL)


###MBT PIP
MBT_PIP_violinplot <- ggplot(MBT, aes(y=`MBT_IR_PIP`,x=`Disease_NO_amyloid`, fill = `Disease_NO_amyloid`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_PIP_violinplot+
  labs(y = "MBT PIP Total Score", x = NULL)


###MBT Free Recall
MBT_FR_violinplot <- ggplot(MBT, aes(y=`MBT_IR_FR_30.120`,x=`Disease_NO_amyloid`, fill = `Disease_NO_amyloid`)) + 
  geom_violin(alpha=0.3)+
  geom_boxplot(width=0.1, alpha=0.6)+
  theme(legend.position="none") +
  scale_fill_grey()
MBT_FR_violinplot+
  labs(y = "MBT FR Total Score", x = NULL)



#MBT Normative data
##select HC only
MBT_HC <- filter(MBT, HC == 1)
view(MBT_HC)

HC_Cases_MBTTIP <- filter(MBT_HC, is.na(MBT_IR_TIP) == F)
HC_Cases_MBTPIP <- filter(MBT_HC, is.na(MBT_IR_PIP) == F)


#percentiles
quantile(MBT_HC$MBT_IR_TIP, c(.1, .2, .3, .4), na.rm = T, type = 8, names = T)

mbt_tip_percentiles <- ecdf(MBT_HC$MBT_IR_TIP)(c(1:32))
mbt_pip_percentiles <- ecdf(MBT_HC$MBT_IR_PIP)(c(1:16))

view(mbt_tip_percentiles)
view(mbt_pip_percentiles)

#normative data descriptive
describe(HC_Cases_MBTTIP$MBT_IR_TIP)
describe(HC_Cases_MBTTIP$Age)

describe(HC_Cases_MBTPIP$MBT_IR_PIP)
describe(HC_Cases_MBTPIP$Age)
