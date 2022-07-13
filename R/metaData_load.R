#metadata
##############################################################################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(methylKit)
library(readxl)
library(lme4)
library(lmerTest)
library(readr)
library(ggeffects)
library(sjPlot)
library(qqman)
library(emmeans)
library(ggggeffects)
library(ggplot2)
library(genomation)
library(dplyr)



#path to directory containing the coverage files

dataPath = "/data/SBCS-EizaguirreLab/Eri/data/coverage/coverage_files/"

# list of paths to all coverage files
temp = list.files(path = dataPath,
                  pattern = ".sam.bismark.cov.gz",
                  full.names = T)

#path to meta data
metData = readxl::read_xlsx("/data/SBCS-EizaguirreLab/Eri/data/coverage/tables/Kost$metData$sex_trt <- paste(metData$Sex,metData$trtG1G2)
metData$sex_trt <- as.factor(metData$sex_trt)
metData$trtG1G2_NUM <- as.numeric(as.factor(metData$trtG1G2))

#read in metaData
fullMetadata4Eri <- read_csv("/data/SBCS-EizaguirreLab/Eri/data/coverage/tables/fullMetadata4Eri_2.csv")
#read metadata with alignment % for all samples (not present in fullMetadata4Eri)
align_for_full<-readxl::read_xlsx("/data/SBCS-EizaguirreLab/Eri/data/coverage/tables/135_samples_align_second.xlsx")
align_for_full<-as.data.frame(align_for_full)
names(align_for_full)[1]<- "Sample_Name"

mapping<-readxl::read_xlsx("/data/SBCS-EizaguirreLab/Eri/data/coverage/reads_mapping_135.xlsx")
mapping<-as.data.frame(mapping)

#merge dataframes
fullMetadata4Eri <- merge(fullMetadata4Eri,mapping,by="Sample_Name")

fullMetadata4Eri <- merge(fullMetadata4Eri,align_for_full,by="Sample_Name")
fullMetadata4Eri$Sex<- as.factor(fullMetadata4Eri$Sex)
fullMetadata4Eri$`% Aligned`<-fullMetadata4Eri$`% Aligned` * 100
fullMetadata4Eri <- merge(fullMetadata4Eri,mapping)

#get mapping % for samples
fullMetadata4Eri$mappingPerc_all_chr <- (fullMetadata4Eri$no_of_mapped_reads/fullMetadata4Eri$no_of_reads)*100

#mapping % only across autosomes
fullMetadata4Eri$mappingPerc_autosomes <- (fullMetadata4Eri$no_of_mapped_reads_no_sex_Chr/fullMetadata4Eri$no_of_reads)*100
#select samples which are only in G2 (offsring)
fullMetadata4Eri<-fullMetadata4Eri[fullMetadata4Eri[ ,26] == "O", ]

#calculate BCI (Body condition Index) for down the line models
fullMetadata4Eri$BCI <- residuals(lmer(Wnettofin ~ Slfin * Sex + (1|brotherPairID), data=fullMetadata4Eri))

fullMetadata4Eri$residuals_BCI_sex<-residuals(lmer(data = fullMetadata4Eri,BCI ~ Sex + (1|brotherPairID)))

fullMetadata4Eri$trtG1G2_NUM <- as.factor(fullMetadata4Eri$trtG1G2)



#convert meta data to data frame (optional)
metData <- data.frame(metData)

#convert treamtent to a variable and give variables values
#(1 = C), (2 = TC), (3 = TT), (4 = T), (5 = CC), (6 = CT)
metData$treatment_NUM <- as.numeric(as.factor(metData$trtG1G2))
metData$Sex<-as.factor(metData$Sex)

save(fullMetadata4Eri,metData, file = "metadata.RData")
