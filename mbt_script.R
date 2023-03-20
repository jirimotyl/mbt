#necessary packages
install.packages("psych")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("readxl")
install.packages("WRS2")
install.packages("pROC")


#necessary libraries
library(psych)
library(tidyverse)
library(ggplot2)
library(readxl)
library(WRS2)
library(pROC)


#source excel file
MBT <- read_excel("MBT_last.xlsx", 
                                                              na = "N/A")

#variable types

MBT$HCaMCIAD <- MBT$`HC vs aMCI AD`
MBT$PDMCIaMCIAD <- MBT$`PD-MCI vs. aMCI AD`
MBT$HCPDMCI <- MBT$`HC vs PD-MCI AD`
MBT$HCPDNC <- MBT$`HC vs PD-NC`

MBT$Disease <- as.factor(MBT$Disease)
MBT$`HC vs aMCI AD` <- as.factor(MBT$`HC vs aMCI AD`)
MBT$`PD-MCI vs. aMCI AD` <- as.factor(MBT$`PD-MCI vs. aMCI AD`)
MBT$`HC vs PD-MCI AD` <- as.factor(MBT$`HC vs PD-MCI AD`)
MBT$HCaMCIAD <- as.factor(MBT$HCaMCIAD)
MBT$PDMCIaMCIAD <- as.factor(MBT$PDMCIaMCIAD)
MBT$HCPDMCI <- as.factor(MBT$HCPDMCI)
MBT$HCPDNC <- as.factor(MBT$HCPDNC)
MBT$Disease_NO <- as.factor(MBT$Disease_NO)
MBT$Disease_NO <- factor(MBT$Disease_NO,
                    levels = c(1,2,3,4),
                    labels = c("HC", "PD-NC", "PD-MCI","AD-aMCI"))
MBT$Disease_NO_unclass <- unclass(MBT$Disease_NO)
MBT$PDncPDmci <- MBT$Disease_NO_unclass
MBT$PDncPDmci[MBT$PDncPDmci == 1] <- NA
MBT$PDncPDmci[MBT$PDncPDmci == 4] <- NA

#select HC only
MBT_HC <- filter(MBT, HC == 1)
view(MBT_HC)

HC_Cases_MBTTIP <- filter(MBT_HC, is.na(MBT_IR_PairWords) == F)
HC_Cases_MBTPIP <- filter(MBT_HC, is.na(MBT_IR_PairPairs) == F)


#percentiles
quantile(MBT_HC$MBT_IR_PairWords, c(.1, .2, .3, .4), na.rm = T, type = 8, names = T)

mbt_tip_percentiles <- ecdf(MBT_HC$MBT_IR_PairWords)(c(1:32))
mbt_pip_percentiles <- ecdf(MBT_HC$MBT_IR_PairPairs)(c(1:16))

view(mbt_tip_percentiles)
view(mbt_pip_percentiles)

#normative data descriptive
describe(HC_Cases_MBTTIP$MBT_IR_PairWords)
describe(HC_Cases_MBTTIP$Age)

describe(HC_Cases_MBTPIP$MBT_IR_PairPairs)
describe(HC_Cases_MBTPIP$Age)

#just trials...
