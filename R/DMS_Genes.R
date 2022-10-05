if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(methylKit)
library(genomation)
library(rtracklayer)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(purrr)
library(forcats)
library(tidyr)
library(AnnotationDbi)
library(GSEABase)
library(Category)
library(GO.db)
library(GOstats)
library(knitr)
library(org.Hs.eg.db)
library(rentrez)
library(devtools)
library(goEnrichment)
library(ggVennDiagram)
library(grDevices)

#####################################################################################################################################################################################################################
#Load Data
######################################################################################################################################################################################################################

load("R_data/myMetaData.RDS ")

annotBed12= readTranscriptFeatures("R_data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.transdec.bed12",remove.unusual = FALSE, up.flank = 1500, down.flank = 500)

getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}
for (i in 1:length(annotBed12)){
  annotBed12[[i]]$name <- getName(annotBed12[[i]]$name)
}

annotBed12$exon$name

GYgynogff = read.table("R_data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")


names(GYgynogff)[1]<-"chrom"
names(GYgynogff)[2]<-"length"

annotGff3 = readGFF("R_data/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.gff")
DMS_Sex <- readRDS("~/R_data/DMS_Sex.RDS")
DMS_Trt <- readRDS("~/R_data/DMS_Trt.RDS")
DMS_Sex_Trt <- readRDS("~/R_data/DMS_Sex_Trt.RDS")
######################################################################################################################################################################################################################
#Where are the DMS' located in context to the genome (exons, introns, promoters,intergenic)
######################################################################################################################################################################################################################

myDiff25p_Sex=getMethylDiff(DMS_Sex,difference=25,qvalue=0.01)
myDiff25p_Trt=getMethylDiff(DMS_Trt,difference=25,qvalue=0.01)
myDiff25p_Sex_Trt=getMethylDiff(DMS_Sex_Trt,difference=25,qvalue=0.01)
myDiff25p_Sex_Trt_sd_2 = myDiff25p_Sex_Trt[myDiff25p_Sex_Trt$meth.diff>25 + 2*sd(abs(myDiff25p_Sex_Trt$meth.diff)),]


#comparing CpG context between males and females
SexANN = annotateWithGeneParts(as(myDiff25p_Sex,"GRanges"),annotBed12)
SexANN_F = annotateWithGeneParts(as(myDiff25p_Sex[myDiff25p_Sex@treatment %in% 1],"GRanges"),annotBed12)
SexANN_M = annotateWithGeneParts(as(myDiff25p_Sex[myDiff25p_Sex@treatment %in% 2],"GRanges"),annotBed12)

#are there any regions which contain more sex-associated DMS' than expected

z_F_s<- SexANN_F@num.precedence
z_M_s<- SexANN_M@num.precedence


for_chi_Sex<- as.data.frame(cbind(z_F_s,z_M_s))
chisq<- chisq.test(for_chi_Sex)
#p-value = 0.2038
chisq

par(mfrow = c(1,2))

SexANN_F_context<- plotTargetAnnotation(SexANN_F,precedence=TRUE, main="SexANN_F", cex.legend = 0.5, border="white")
SexANN_M_context<- plotTargetAnnotation(SexANN_M,precedence=TRUE, main="SexANN_M", cex.legend = 0.5, border="white")

#comparing CpG context between the 4 treatment groups
TrtAnn = annotateWithGeneParts(as(myDiff25p_Trt,"GRanges"),annotBed12)
TrtAnn_E_control = annotateWithGeneParts(as(myDiff25p_Trt[myDiff25p_Trt@treatment %in% 1],"GRanges"),annotBed12)
TrtAnn_E_exposed = annotateWithGeneParts(as(myDiff25p_Trt[myDiff25p_Trt@treatment %in% 2],"GRanges"),annotBed12)
TrtAnn_NE_control = annotateWithGeneParts(as(myDiff25p_Trt[myDiff25p_Trt@treatment %in% 3],"GRanges"),annotBed12)
TrtAnn_NE_exposed = annotateWithGeneParts(as(myDiff25p_Trt[myDiff25p_Trt@treatment %in% 4],"GRanges"),annotBed12)

#are there any regions which contain more treatment-associated DMS' than expected

z_CC_s<- TrtAnn_NE_control@num.precedence
z_CE_s<- TrtAnn_NE_exposed@num.precedence
z_EE_s<- TrtAnn_E_exposed@num.precedence
z_EC_s<- TrtAnn_E_control@num.precedence



for_chi_Trt<- as.data.frame(cbind(z_CC_s,z_CE_s,z_EE_s,z_EC_s))
chisq_trt<-chisq.test(for_chi_Trt)
#p-value = 0.1468
chisq_trt

par(mfrow = c(2,2))

TrtAnn_E_control_Context<-plotTargetAnnotation(TrtAnn_E_control,precedence=TRUE, main="TrtAnn_E_control", cex.legend = 0.5, border="white")
TrtAnn_E_exposed_Context<-plotTargetAnnotation(TrtAnn_E_exposed,precedence=TRUE, main="TrtAnn_E_exposed", cex.legend = 0.5, border="white")
TrtAnn_NE_control_Context<-plotTargetAnnotation(TrtAnn_NE_control,precedence=TRUE, main="TrtAnn_NE_control", cex.legend = 0.5, border="white")
TrtAnn_NE_exposed_Context<-plotTargetAnnotation(TrtAnn_NE_exposed,precedence=TRUE, main="TrtAnn_NE_exposed", cex.legend = 0.5, border="white")


#comparing CpG context between 4 treatment groups' interaction with sex
Sex_Trt_ANN = annotateWithGeneParts(as(myDiff25p_Sex_Trt,"GRanges"),annotBed12)
Sex_Trt_ANN_F_EC = annotateWithGeneParts(as(myDiff25p_Sex_Trt_sd_2[myDiff25p_Sex_Trt_sd_2@treatment %in% 1],"GRanges"),annotBed12)
Sex_Trt_ANN_M_EC = annotateWithGeneParts(as(myDiff25p_Sex_Trt_sd_2[myDiff25p_Sex_Trt_sd_2@treatment %in% 2],"GRanges"),annotBed12)
Sex_Trt_ANN_F_EE = annotateWithGeneParts(as(myDiff25p_Sex_Trt_sd_2[myDiff25p_Sex_Trt_sd_2@treatment %in% 3],"GRanges"),annotBed12)
Sex_Trt_ANN_M_EE = annotateWithGeneParts(as(myDiff25p_Sex_Trt_sd_2[myDiff25p_Sex_Trt_sd_2@treatment %in% 4],"GRanges"),annotBed12)
Sex_Trt_ANN_F_CC = annotateWithGeneParts(as(myDiff25p_Sex_Trt_sd_2[myDiff25p_Sex_Trt_sd_2@treatment %in% 5],"GRanges"),annotBed12)
Sex_Trt_ANN_M_CC = annotateWithGeneParts(as(myDiff25p_Sex_Trt_sd_2[myDiff25p_Sex_Trt_sd_2@treatment %in% 6],"GRanges"),annotBed12)
Sex_Trt_ANN_F_CE = annotateWithGeneParts(as(myDiff25p_Sex_Trt_sd_2[myDiff25p_Sex_Trt_sd_2@treatment %in% 7],"GRanges"),annotBed12)
Sex_Trt_ANN_M_CE = annotateWithGeneParts(as(myDiff25p_Sex_Trt_sd_2[myDiff25p_Sex_Trt_sd_2@treatment %in% 8],"GRanges"),annotBed12)

#are there any regions which contain more Sex*treatment-associated DMS' than expected

z_CC_M<- Sex_Trt_ANN_M_CC@num.precedence
z_CC_F<- Sex_Trt_ANN_F_CC@num.precedence

z_CE_M<- Sex_Trt_ANN_M_CE@num.precedence
z_CE_F<- Sex_Trt_ANN_F_CE@num.precedence

z_EE_M<- Sex_Trt_ANN_M_EE@num.precedence
z_EE_F<- Sex_Trt_ANN_F_EE@num.precedence

z_EC_M<- Sex_Trt_ANN_M_EC@num.precedence
z_EC_F<- Sex_Trt_ANN_F_EC@num.precedence


for_chi_Trt_Sex<- as.data.frame(cbind(z_CC_M,z_CC_F,z_CE_M,z_CE_F,
                                  z_EE_M,z_EE_F,z_EC_M,z_EC_F))

for_chi_Trt_Sex_EC<- as.data.frame(cbind(z_CC_M,z_CC_F))
chisq_Trt_Sex<-chisq.test(for_chi_Trt_Sex_EC)
#p-value = 0.4015
chisq_Trt_Sex

#Vizualisation of DMS' in contect to genomic regions

par(mfrow = c(2,4))

Sex_Trt_ANN_F_EC_CON<-plotTargetAnnotation(Sex_Trt_ANN_F_EC,precedence=TRUE, main="F_EC", cex.legend = 0.4, border="white")
Sex_Trt_ANN_M_EC_CON<-plotTargetAnnotation(Sex_Trt_ANN_M_EC,precedence=TRUE, main="M_EC", cex.legend = 0.5, border="white")
Sex_Trt_ANN_F_EE_CON<-plotTargetAnnotation(Sex_Trt_ANN_F_EE,precedence=TRUE, main="F_EE", cex.legend = 0.5, border="white")
Sex_Trt_ANN_M_EE_CON<-plotTargetAnnotation(Sex_Trt_ANN_M_EE,precedence=TRUE, main="M_EE", cex.legend = 0.5, border="white")
Sex_Trt_ANN_F_CC_CON<-plotTargetAnnotation(Sex_Trt_ANN_F_CC,precedence=TRUE, main="F_CC", cex.legend = 0.5, border="white")
Sex_Trt_ANN_m_CC_CON<-plotTargetAnnotation(Sex_Trt_ANN_M_CC,precedence=TRUE, main="M_CC", cex.legend = 0.5, border="white")
Sex_Trt_ANN_F_CE_CON<-plotTargetAnnotation(Sex_Trt_ANN_F_CE,precedence=TRUE, main="F_EC", cex.legend = 0.5, border="white")
Sex_Trt_ANN_M_CE_CON<-plotTargetAnnotation(Sex_Trt_ANN_M_CE,precedence=TRUE, main="M_EC", cex.legend = 0.5, border="white")

Sex_Trt_ANN_F_EC_CON


######################################################################################################################################################################################################################
#what is the dispersion of the DMS' in terms of position and chromosome and how does the genome-wide methylation pattern differ between treatments
######################################################################################################################################################################################################################


annot_Sex <- as.data.frame(SexANN@members)
Sex_Diff_man_Plot<-makeManhattanPlots(DMSfile = myDiff25p_Sex, annotFile = annot_Sex, GYgynogff = GYgynogff, 
                                      mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of Sex Diff")
Sex_Diff_man_Plot

outliers_Sex = which(abs(myDiff25p_Sex$meth.diff)>25 + 2*sd(abs(myDiff25p_Sex$meth.diff)))

outliers_annot_Sex <- as.data.frame(SexANN_M@members)[outliers_Sex,]

Sex_biggest_Diff_man_Plot<-makeManhattanPlots(DMSfile = myDiff25p_Sex[outliers_Sex,], annotFile = outliers_annot_Sex, GYgynogff = GYgynogff, 
                                              mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of biggest Sex Diff")

Sex_biggest_Diff_man_Plot

annot_Trt <- as.data.frame(TrtAnn@members)
Trt_Diff_man_Plot<-makeManhattanPlots(DMSfile = myDiff25p_Trt, annotFile = annot_Trt, GYgynogff = GYgynogff, 
                                      mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of Trt Diff")
Trt_Diff_man_Plot


annot_Sex_Trt <- as.data.frame(Sex_Trt_ANN@members)
Sex_Trt_Diff_man_Plot<-makeManhattanPlots(DMSfile = myDiff25p_Sex_Trt, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of Sex+Trt Diff")
Sex_Trt_Diff_man_Plot


pos_of_interest <- paste(myDiff25p_Sex_Trt[["chr"]],myDiff25p_Sex_Trt[["start"]])

unite_for_DMS_Sex$start
#there should be 12559 CpGs
length(myDiff25p_Sex_Trt$chr)
table(unite_for_DMS_Sex_Trt$chr == myDiff25p_Sex_Trt$chr & 
        unite_for_DMS_Sex_Trt$start == myDiff25p_Sex_Trt$start)

table(duplicated(paste(unite_for_DMS_Sex_Trt$chr,unite_for_DMS_Sex_Trt$start)))


table(paste(unite_for_DMS_Sex_Trt$chr, unite_for_DMS_Sex_Trt$start) %in% paste(myDiff25p_Sex_Trt$chr, myDiff25p_Sex_Trt$start))

unite_for_DMS_Sex_Trt$pos <- paste(unite_for_DMS_Sex_Trt$chr, unite_for_DMS_Sex_Trt$start)
myDiff25p_Sex_Trt$pos <- paste(myDiff25p_Sex_Trt$chr, myDiff25p_Sex_Trt$start)
rows_to_keep<-which(unite_for_DMS_Sex_Trt$pos %in% myDiff25p_Sex_Trt$pos)
head(rows_to_keep)

methyl_unite_rows_selected<-methylKit::select(x = unite_for_DMS_Sex_Trt,i = rows_to_keep)

percMeth_for_trt_sex<-function(methobj,trt){
  obj2<-methylKit::reorganize(methobj, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% trt],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% trt])
  
  perc_obj2<-percMethylation(obj2)
  row_mean_trt_1 <-rowMeans(perc_obj2,na.rm = T)
  
}

F_E_control<-percMeth_for_trt(methyl_unite_rows_selected,trt = 1)
F_E_control_df<-data.frame(F_E_control)
F_E_control_df$chr <- obj2$chr
F_E_control_df$pos <- obj2$start
F_E_control_df$pos <- as.numeric(obj2$start)
colnames(F_E_control_df)[1] <- "avg.meth"

M_E_control<-percMeth_for_trt(methyl_unite_rows_selected,trt = 2)
F_E_exposed<-percMeth_for_trt(methyl_unite_rows_selected,trt = 3)
M_E_exposed<-percMeth_for_trt(methyl_unite_rows_selected,trt = 4)

F_NE_control<-percMeth_for_trt(methyl_unite_rows_selected,trt = 5)
F_NE_control_df<-data.frame(F_E_control)
F_NE_control_df$chr <- obj2$chr
F_NE_control_df$pos <- obj2$start
F_NE_control_df$pos <- as.numeric(obj2$start)
colnames(F_NE_control_df)[1] <- "avg.meth"



M_NE_control<-percMeth_for_trt(methyl_unite_rows_selected,trt = 6)
F_NE_exposed<-percMeth_for_trt(methyl_unite_rows_selected,trt = 7)
M_NE_exposed<-percMeth_for_trt(methyl_unite_rows_selected,trt = 8)

Sex_Trt_Diff_man_Plot<-makeManhattanPlots(DMSfile = F_E_control_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of Sex+Trt Diff")

Sex_Trt_Diff_man_Plot



par(mfrow=c(2,2))
plot1<-Sex_Trt_1
plot2<-Sex_Trt_2
plot3<-Sex_Trt_3
plot4<-Sex_Trt_4
plot5<-Sex_Trt_5
plot6<-Sex_Trt_6
plot7<-Sex_Trt_7
plot8<-Sex_Trt_8
grid.arrange(plot5,plot6,plot7,plot8, plot1, plot2,plot3,plot4,ncol=1,nrow=8)
#looking at Sex diif between treatments
grid.arrange(plot1, plot2, ncol=1)
grid.arrange(plot3, plot4, ncol=1)
grid.arrange(plot5, plot6, ncol=1)
grid.arrange(plot7, plot8, ncol=1)

#looking at trt differences
grid.arrange(plot_3,plot_4,plot_1,plot_2,ncol = 1)

unite_for_DMS_Sex_Trt<-readRDS("R_data/unite_for_DMS_Sex_Trt.RDS")
unite_for_DMS_Sex_Trt<- unite_for_DMS_Sex_Trt[!unite_for_DMS_Sex_Trt$chr %in% c("Gy_chrXIX","Gy_chrUn"),]

#####################################################################################################
#Reorganize Sex_trt into groups the 4 treatment groups with 2 subgroups sex
#####################################################################################################

EC<-methylKit::reorganize(unite_for_DMS_Sex_Trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(1, 2) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(1 2,)])
EC@treatment

EE<-methylKit::reorganize(unite_for_DMS_Sex_Trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(3, 4) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(3, 4)])
EE@treatment

CC<-methylKit::reorganize(unite_for_DMS_Sex_Trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(5, 6) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(5, 6)])
CC@treatment


CE<-methylKit::reorganize(unite_for_DMS_Sex_Trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(7, 8) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(7, 8)])
CE@treatment


#####################################################################################################
#Find the DMS' between sexes within each group
#####################################################################################################

EC_samples<-EC@sample.ids
EC_covariates <- as.data.frame(myMetaData$brotherPairID[ myMetaData$SampleID %in% EC_samples] )

ec_dif<-calculateDiffMeth(EC,covariates = EC_covariates,mc.cores = 3)
ec_dif$qvalue<-p.adjust(ec_dif$qvalue,"BH")

ec_dif<-readRDS("R_data/ec_diff.RDS")

ec_get_dif<-getMethylDiff(ec_dif,difference = 15,qvalue = 0.01)

saveRDS(ec_dif,file = "ec_diff.RDS")


EE_samples<-EE@sample.ids
EE_covariates <- as.data.frame(myMetaData$brotherPairID[ myMetaData$SampleID %in% EE_samples] )

ee_dif<-calculateDiffMeth(EE,covariates = EE_covariates,mc.cores = 3)
ee_dif$qvalue<-p.adjust(ee_dif$qvalue,"BH")

ee_get_dif<-getMethylDiff(ee_dif,difference = 15,qvalue = 0.01)

saveRDS(ee_dif,file = "ee_diff.RDS")

ee_dif<-readRDS("R_data/ee_diff.RDS")


CC_samples<-CC@sample.ids
CC_covariates <- as.data.frame(myMetaData$brotherPairID[ myMetaData$SampleID %in% CC_samples] )


cc_dif<-calculateDiffMeth(CC,covariates = covariates,mc.cores = 3)
cc_diff$qvalue<-p.adjust(cc_diff$qvalue,"BH")

cc_get_dif<-getMethylDiff(cc_diff,difference = 15,qvalue = 0.01)

saveRDS(cc_dif,file = "cc_diff.RDS")
cc_diff<-readRDS("R_data/cc_diff.RDS")



CE_samples<-CE@sample.ids
CE_covariates <- as.data.frame(myMetaData$brotherPairID[ myMetaData$SampleID %in% CE_samples] )

ce_dif<-calculateDiffMeth(CE,covariates = CE_covariates,mc.cores = 3)
ce_dif$qvalue<-p.adjust(ce_dif$qvalue,"BH")

ce_get_dif<-getMethylDiff(ce_dif,difference = 15,qvalue = 0.01)

ce_dif<-saveRDS(ce_dif,file = "ce_diff.RDS")
ce_dif<-readRDS("R_data/ce_diff.RDS")




#####################################################################################################
#Separate the Male and Female DMS'
#####################################################################################################

#DMS' with a meth.diff over 0 are female as they are the basis of the comparison
F_DMS_SEX_TRT_CC<-cc_get_dif[cc_get_dif$meth.diff > 0,]
M_DMS_SEX_TRT_CC<-cc_get_dif[cc_get_dif$meth.diff < 0,]

F_DMS_SEX_TRT_EC<-ec_get_dif[ec_get_dif$meth.diff > 0,]
M_DMS_SEX_TRT_EC<-ec_get_dif[ec_get_dif$meth.diff < 0,]

F_DMS_SEX_TRT_EE<-ee_get_dif[ee_get_dif$meth.diff > 0,]
M_DMS_SEX_TRT_EE<-ee_get_dif[ee_get_dif$meth.diff < 0,]

F_DMS_SEX_TRT_CE<-ce_get_dif[ce_get_dif$meth.diff > 0,]
M_DMS_SEX_TRT_CE<-ce_get_dif[ce_get_dif$meth.diff < 0,]

#####################################################################################################
#Compare overlap of DMS' within treatment groups and sexes
#####################################################################################################

listsex = list(Exposed_Exposed_M = paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start),
               Control_Exposed_M = paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start))
ggVennDiagram(listsex) + scale_fill_gradient(low = "blue", high = "red")

list_F=list(Control_Exposed_F = paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start),
            Exposed_Exposed_F = paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start))
ggVennDiagram(list_F) + scale_fill_gradient(low = "blue", high = "red")

#reorganise does not work on methyldiff obj
#ce_M_get_diff<-methylKit::reorganize(ce_get_dif, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 8],
#sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 8])

ce_M<-getData(ce_get_dif)

ce_M_get_diff<-ce_get_dif[,ce_get_dif@treatment %in% 8]
ce_F_get_diff<-ce_get_dif[,ce_get_dif@treatment %in% 7]

listx = list(ce_M_get_diff = paste(ce_M_get_diff$chr,ce_M_get_diff$start),
             ce_F_get_diff = paste(ce_F_get_diff$chr,ce_F_get_diff$start))
ggVennDiagram(listx) + scale_fill_gradient(low = "blue", high = "red")


listx = list(Exposed_Exposed = paste(ee_get_dif$chr,ee_get_dif$start),
             Control_Exposed = paste(ce_get_dif$chr,ce_get_dif$start),
             Control_Control = paste(cc_get_dif$chr,cc_get_dif$start),
             Exposed_Control = paste(ec_get_dif$chr,ec_get_dif$start))

pdf(file = "Figures_for_diss/venn_diagram.pdf",width = 10 ,height = 10)

ggVennDiagram(listx) + 
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_continuous(expand = expansion(mult = .2))
dev.off()


#############################################################################################
#2012 DMS' between Sex in ee trt
ee_get_dif_exclusive<-ee_get_dif[!paste(ee_get_dif$chr,ee_get_dif$start) %in% 
                                   c(paste(cc_get_dif$chr,cc_get_dif$start),
                                     paste(ce_get_dif$chr,ce_get_dif$start),
                                     paste(ec_get_dif$chr,ec_get_dif$start)),]
#860 DMS' are which are exclusive to ee trt are females
ee_get_dif_exclusive_F<-ee_get_dif_exclusive[ee_get_dif_exclusive$meth.diff>0]
#1152 DMS' are which are exclusive to ee trt are females
ee_get_dif_exclusive_M<-ee_get_dif_exclusive[ee_get_dif_exclusive$meth.diff<0]


#############################################################################################

#884 DMS' between Sex in ec trt
ec_get_dif_exclusive<-ec_get_dif[!paste(ec_get_dif$chr,ec_get_dif$start) %in% 
                                   c(paste(cc_get_dif$chr,cc_get_dif$start),
                                     paste(ce_get_dif$chr,ce_get_dif$start),
                                     paste(ee_get_dif$chr,ee_get_dif$start)),]
#379 DMS' are which are exclusive to ec trt are females
ec_get_dif_exclusive_F<-ec_get_dif_exclusive[ec_get_dif_exclusive$meth.diff>0]
#505 DMS' are which are exclusive to ec trt are males
ec_get_dif_exclusive_M<-ec_get_dif_exclusive[ec_get_dif_exclusive$meth.diff<0]
#############################################################################################
#1153 DMS' between Sex in ec trt
cc_get_dif_exclusive<-cc_get_dif[!paste(cc_get_dif$chr,cc_get_dif$start) %in% 
                                   c(paste(ee_get_dif$chr,ee_get_dif$start),
                                     paste(ce_get_dif$chr,ce_get_dif$start),
                                     paste(ec_get_dif$chr,ec_get_dif$start)),]
#618 DMS' are which are exclusive to ec trt are females
cc_get_dif_exclusive_F<-cc_get_dif_exclusive[cc_get_dif_exclusive$meth.diff>0]
#535 DMS' are which are exclusive to ec trt are males
cc_get_dif_exclusive_M<-cc_get_dif_exclusive[cc_get_dif_exclusive$meth.diff<0]
#############################################################################################
#1497 DMS' between Sex in ec trt
ce_get_dif_exclusive<-ce_get_dif[!paste(ce_get_dif$chr,ce_get_dif$start) %in% 
                                   c(paste(cc_get_dif$chr,cc_get_dif$start),
                                     paste(ee_get_dif$chr,ee_get_dif$start),
                                     paste(ec_get_dif$chr,ec_get_dif$start)),]
#737 DMS' are which are exclusive to ec trt are females
ce_get_dif_exclusive_F<-ce_get_dif_exclusive[ce_get_dif_exclusive$meth.diff>0]

#760 DMS' are which are exclusive to ec trt are males
ce_get_dif_exclusive_M<-ce_get_dif_exclusive[ce_get_dif_exclusive$meth.diff<0]

#############################################################################################


DMS_for_annot_ee_get_dif<-paste(ee_get_dif_exclusive$chr, ee_get_dif_exclusive$start,sep = " ")
DMS_for_annot_ec_get_dif<-paste(ec_get_dif_exclusive$chr, ec_get_dif_exclusive$start,sep = " ")
DMS_for_annot_cc_get_dif<-paste(cc_get_dif_exclusive$chr, cc_get_dif_exclusive$start,sep = " ")
DMS_for_annot_ce_get_dif<-paste(ce_get_dif_exclusive$chr, ce_get_dif_exclusive$start,sep = " ")

###############################################################################################

getAnnotationFun <- function(METHOBJ){
  A = annotateWithGeneParts(target = as(METHOBJ,"GRanges"), feature = annotBed12)
  # Heckwolf 2020: To be associated to a gene, the pop-DMS had to be either inside the gene or,
  # if intergenic, not further than 10 kb away from the TSS.
  rows2rm = which((A@dist.to.TSS$dist.to.feature>10000 |
                     A@dist.to.TSS$dist.to.feature< -10000) &
                    rowSums(A@members) %in% 0)
  if (is_empty(rows2rm)){
    METHOBJ2 = METHOBJ
  } else {
    METHOBJ2 = METHOBJ[-rows2rm,]
  }
  ## Re annotate the subsetted object
  B = annotateWithGeneParts(as(METHOBJ2,"GRanges"),annotBed12)
  ## Get genes associated
  C = getAssociationWithTSS(B)
  ## How many CpG per gene?
  #nCpG = table(C$feature.name)
  ## Get annotations for these genes
  subAnnot <- data.frame(subset(annotGff3, Name %in% C$feature.name))
  #subAnnot$nCpG = nCpG
  nCpGdf = data.frame(table(C$feature.name))
  names(nCpGdf) = c("Name", "nCpG")
  # Merge both by Name
  subAnnot = merge(subAnnot, nCpGdf)
  
  # Add extra info (nbr CpG per gene length, gene length, chrom name)
  subAnnot = subAnnot  %>% 
    mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb, chrom = seqid)
  
  
  return(subAnnot)
}


######################################################################################
df_ee <- data.frame(chr=sapply(strsplit(DMS_for_annot_ee_get_dif, " "), `[`, 1), 
                    start=sapply(strsplit(DMS_for_annot_ee_get_dif, " "), `[`, 2),
                    end=sapply(strsplit(DMS_for_annot_ee_get_dif, " "), `[`, 2)) 
# get annotation
anot_ee_get_dif <- getAnnotationFun(makeGRangesFromDataFrame(df_ee))

#816 Genes
nrow(anot_ee_get_dif)

anot_ee_get_dif_ <- anot_ee_get_dif %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

top_ee_genes = anot_ee_get_dif_[anot_ee_get_dif_$nCpGperGenekb>=1,]
######################################################################################
#### Genes asscoiated with ec 
df_ec <- data.frame(chr=sapply(strsplit(DMS_for_annot_ec_get_dif, " "), `[`, 1), 
                    start=sapply(strsplit(DMS_for_annot_ec_get_dif, " "), `[`, 2),
                    end=sapply(strsplit(DMS_for_annot_ec_get_dif, " "), `[`, 2)) 
# get annotation
anot_ec_get_dif <- getAnnotationFun(makeGRangesFromDataFrame(df_ec))

#389 Genes
nrow(anot_ec_get_dif)

anot_ec_get_dif_ <- anot_ec_get_dif %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

top_ec_genes = anot_ec_get_dif_[anot_ec_get_dif_$nCpGperGenekb>=1,]

######################################################################################
df_ce <- data.frame(chr=sapply(strsplit(DMS_for_annot_ce_get_dif, " "), `[`, 1), 
                    start=sapply(strsplit(DMS_for_annot_ce_get_dif, " "), `[`, 2),
                    end=sapply(strsplit(DMS_for_annot_ce_get_dif, " "), `[`, 2)) 
# get annotation
anot_ce_get_dif <- getAnnotationFun(makeGRangesFromDataFrame(df_ce))

#674 Genes
nrow(anot_ce_get_dif)

anot_ce_get_dif_ <- anot_ce_get_dif %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

top_ce_genes = anot_ce_get_dif_[anot_ce_get_dif_$nCpGperGenekb>=1,]

######################################################################################

df_cc <- data.frame(chr=sapply(strsplit(DMS_for_annot_cc_get_dif, " "), `[`, 1), 
                    start=sapply(strsplit(DMS_for_annot_cc_get_dif, " "), `[`, 2),
                    end=sapply(strsplit(DMS_for_annot_cc_get_dif, " "), `[`, 2)) 
# get annotation
anot_cc_get_dif <- getAnnotationFun(makeGRangesFromDataFrame(df_cc))

#453 Genes
nrow(anot_cc_get_dif)

anot_cc_get_dif_ <- anot_cc_get_dif %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

top_cc_genes = anot_cc_get_dif_[anot_cc_get_dif_$nCpGperGenekb>=1,]

######################################################################################

getGeneSummary <- function(topgenes){
  
  df = data.frame(GeneSymbol = sapply(topgenes$Note, function(x) sub("Similar to ", "", x) %>% str_extract(".*:") %>% str_remove(":") %>% toupper),  # extract uniprot symbol from note, then uppercase
                  seqid = topgenes$seqid) #add info: on how many BP? which chr? how many Cpg ar diffmeth?nCpGperGenekb = topgenes$nCpGperGenekb, 
  
  
  # Convert the uniprot gene names to entrez ids
  topGenesDF = unlist(mapIds(org.Hs.eg.db, keys = df$GeneSymbol, column = "ENTREZID", keytype = "SYMBOL")) %>% data.frame() 
  names(topGenesDF) = "ENTREZID" ; topGenesDF$GeneSymbol = rownames(topGenesDF) ; rownames(topGenesDF) = NULL
  # Retrieve gene summary & description
  genes = entrez_summary(db="gene", id=topGenesDF$ENTREZID)
  SummaDF = lapply(genes, function(x) x[["summary"]]) %>% unlist() %>% data.frame()
  names(SummaDF) = "Summary" ; SummaDF$ENTREZID = rownames(SummaDF) ; rownames(SummaDF) = NULL
  DescDF = lapply(genes, function(x) x[["description"]]) %>% unlist() %>% data.frame()
  names(DescDF) = "description" ; DescDF$ENTREZID = rownames(DescDF) ; rownames(DescDF) = NULL
  # Output complete table
  topGenesDF = merge(merge(DescDF, topGenesDF, all = T), SummaDF, all = T)
  
  
  
  # Add Number of brother pairs in which the gene is found
  topGenesDF = merge(topGenesDF, df)
}

top_cc_genes_df<-getGeneSummary(top_cc_genes)
top_ce_genes_df<-getGeneSummary(top_ce_genes)
top_ec_genes_df<-getGeneSummary(top_ec_genes)
top_ee_genes_df<-getGeneSummary(top_ee_genes)
#######################################################################################################
#are some of the genes shared between treatments? 
top_cc_genes_df$GeneSymbol
top_ce_genes_df$GeneSymbol
top_ee_genes_df$GeneSymbol
top_ec_genes_df$GeneSymbol

list_of_genes = list(CC = top_cc_genes_df$GeneSymbol,
                     CE = top_ce_genes_df$GeneSymbol,
                     EE = top_ee_genes_df$GeneSymbol,
                     EC = top_ec_genes_df$GeneSymbol)

pdf(file = "Figures_for_diss/gene_venn.pdf",width = 10 ,height = 10)

ggVennDiagram(list_of_genes) + scale_fill_gradient(low = "blue", high = "red")


dev.off()
ggVennDiagram(list_of_genes) + scale_fill_gradient(low = "blue", high = "red")

#little overlap with CC and CE having 2 genes in common and CE having 1 in common in EC
matching_genes<-top_ce_genes_df[top_ce_genes_df$GeneSymbol %in% 
                                  c(top_ec_genes_df$GeneSymbol,
                                    top_ee_genes_df$GeneSymbol,
                                    top_cc_genes_df$GeneSymbol),]

#"OR11A1" Olfactory receptors (GPCR)
#"RPL28"  Among its related pathways are Peptide chain elongation and Influenza Infection
#"ZBTB26", Predicted to be involved in regulation of transcription by RNA polymerase II
matching_genes$GeneSymbol

#######################################################################################################
gene_universe <- data.frame(
  subsetByOverlaps(GRanges(annotGff3), GRanges(unite_for_DMS_Sex_Trt))) %>% # subselect covered CpGs
  filter(lengths(Ontology_term)!=0) %>% # rm non existing GO terms
  filter(type %in% "gene")  %>% # keep all the 7416 genes with GO terms
  dplyr::select(c("Name", "Ontology_term")) %>%
  mutate(go_linkage_type = "IEA") %>% #NB: IEA but not necessarily true, it's from Interproscan after Maker. Sticklebacks (biomart) have 82701 IEA and 63 ISS.
  relocate("Ontology_term","go_linkage_type","Name") %>%
  unnest(Ontology_term) %>% # one GO per line (was a list before in this column)
  data.frame()

goFrame <- GOFrame(gene_universe, organism="Gasterosteus aculeatus")
goAllFrame <- GOAllFrame(goFrame)
gsc_universe <- GeneSetCollection(goAllFrame, setType = GOCollection())

makedfGO <- function(annot_trt_get_diff, gene_universe){
  #plotManhattanGenes(i)$anot4BP
  
  ## Create subuniverse:
  sub_universe <- gene_universe %>%
    subset(gene_universe$Name %in% unlist(annot_trt_get_diff$Parent))
  
  ## Run conditional hypergeometric test:
  runTestHypGeom <- function(sub_universe, onto){
    ## Constructing a GOHyperGParams objects or KEGGHyperGParams objects from a GeneSetCollection
    ## Then run hypergeometric test:
    GO_NO_fdr <- hyperGTest(GSEAGOHyperGParams(name="GO_set",
                                               geneSetCollection = gsc_universe,
                                               geneIds = as.vector(unique(sub_universe[["Name"]])), # gene ids for the selected gene set
                                               universeGeneIds = unique(gene_universe$Name),
                                               ontology = onto, # A string for GO to use ("BP", "CC", or "MF")
                                               pvalueCutoff = 0.05,
                                               conditional = TRUE, # see note above
                                               testDirection = "over")) # over represented GO terms
    
    
    
    # Use GOEnrich as a wrapper around GOStat for extra FDR comparison
    ## Does not solve all issues, but better than nothing. See: https://support.bioconductor.org/p/5571/
    GO_fdr <- joinGOEnrichResults(goEnrichTest(gsc=gsc_universe,
                                               gene.ids = as.vector(unique(sub_universe[["Name"]])),# genes in selected gene set
                                               univ.gene.ids = unique(gene_universe$Name),
                                               ontologies = onto, # A string for GO to use ("BP", "CC", or "MF")
                                               pvalue.cutoff = 0.05,
                                               cond = TRUE, # see note above
                                               test.dir = "over"),# over represented GO terms
                                  p.adjust.method = "fdr")
    
    
    return(list(GO_NO_fdr=GO_NO_fdr, GO_fdr=GO_fdr))
  }
  
  
  
  GO_MF <- runTestHypGeom(sub_universe = sub_universe, onto = "MF")
  GO_CC <- runTestHypGeom(sub_universe = sub_universe, onto = "CC")
  GO_BP <- runTestHypGeom(sub_universe = sub_universe, onto = "BP")
  
  # Get percentage of genes over reppresented in universe
  dfMFperc = GO_MF$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  dfCCperc = GO_CC$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  dfBPperc = GO_BP$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  
  # Add this information to FDR corrected table
  GO_MF_all = merge(GO_MF$GO_fdr, dfMFperc)
  GO_CC_all = merge(GO_CC$GO_fdr, dfCCperc)
  GO_BP_all = merge(GO_BP$GO_fdr, dfBPperc)
  
  # Merge the df MP and BP
  dfGO = rbind(GO_MF_all, GO_CC_all, GO_BP_all)
  
  dfGO = dfGO %>% mutate(Term = factor(x = GO.term, levels = GO.term))
  
  # Relabel GO group names
  dfGO$GO.category[dfGO$GO.category %in% "CC"]="Cellular components"
  dfGO$GO.category[dfGO$GO.category %in% "BP"]="Biological processes"
  dfGO$GO.category[dfGO$GO.category %in% "MF"]="Molecular functions"
  
  return(dfGO)
}

df_GO_cc<-makedfGO(annot_trt_get_diff = anot_cc_get_dif, gene_universe = gene_universe)
df_GO_ce<-makedfGO(annot_trt_get_diff = anot_ce_get_dif, gene_universe = gene_universe)
df_GO_ec<-makedfGO(annot_trt_get_diff = anot_ec_get_dif, gene_universe = gene_universe)
df_GO_ee<-makedfGO(annot_trt_get_diff = anot_ee_get_dif, gene_universe = gene_universe)

df_GO_cc$trtG1G2<-"CONTROL-CONTROL"
df_GO_ce$trtG1G2<-"CONTROL-EXPOSED"
df_GO_ec$trtG1G2<-"EXPOSED-CONTROL"
df_GO_ee$trtG1G2<-"EXPOSED-EXPOSED"

dfGO <- rbind(df_GO_cc, df_GO_ce, df_GO_ec, df_GO_ee)

dfGO_trt<-dfGO %>% ggplot(aes(x=trtG1G2,y = factor(GO.name))) +
  geom_point(aes(color = p.value.adjusted, size = genePercent)) +
  scale_color_gradient(name="adjusted\np-value", low = "red", high = "blue") +
  scale_size_continuous(name = "% of genes")+
  theme_bw() + ylab("") + xlab("Sex Differences within each treatment") +
  theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"), 
        legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"), 
        legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), legend.position="left",
        plot.margin = margin(1,1,1.5,1.2, "cm")) + # grey box for legend
  facet_grid((GO.category)~trtG1G2, scales="free",space = "free") 


pdf(file = "Figures_for_diss/dfGO_trt.pdf",width = 10 ,height = 25)

dfGO_trt

dev.off()

F_DMS_SEX_TRT_CC<-cc_get_dif[cc_get_dif$meth.diff > 0,]
M_DMS_SEX_TRT_CC<-cc_get_dif[cc_get_dif$meth.diff < 0,]

F_DMS_SEX_TRT_EC<-ec_get_dif[ec_get_dif$meth.diff > 0,]
M_DMS_SEX_TRT_EC<-ec_get_dif[ec_get_dif$meth.diff < 0,]

F_DMS_SEX_TRT_EE<-ee_get_dif[ee_get_dif$meth.diff > 0,]
M_DMS_SEX_TRT_EE<-ee_get_dif[ee_get_dif$meth.diff < 0,]

F_DMS_SEX_TRT_CE<-ce_get_dif[ce_get_dif$meth.diff > 0,]
M_DMS_SEX_TRT_CE<-ce_get_dif[ce_get_dif$meth.diff < 0,]


F_DMS_SEX_TRT_CC_exclusive<-F_DMS_SEX_TRT_CC[!paste(F_DMS_SEX_TRT_CC$chr,F_DMS_SEX_TRT_CC$start) %in% 
                                               c(paste(M_DMS_SEX_TRT_CC$chr,M_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_EC$chr,F_DMS_SEX_TRT_EC$start),
                                                 paste(M_DMS_SEX_TRT_EC$chr,M_DMS_SEX_TRT_EC$start),
                                                 paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start),
                                                 paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start),
                                                 paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start),
                                                 paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start)),]

M_DMS_SEX_TRT_CC_exclusive<-M_DMS_SEX_TRT_CC[!paste(M_DMS_SEX_TRT_CC$chr,M_DMS_SEX_TRT_CC$start) %in% 
                                               c(paste(F_DMS_SEX_TRT_CC$chr,F_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_EC$chr,F_DMS_SEX_TRT_EC$start),
                                                 paste(M_DMS_SEX_TRT_EC$chr,M_DMS_SEX_TRT_EC$start),
                                                 paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start),
                                                 paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start),
                                                 paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start),
                                                 paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start)),]

F_DMS_SEX_TRT_EC_exclusive<-F_DMS_SEX_TRT_EC[!paste(F_DMS_SEX_TRT_EC$chr,F_DMS_SEX_TRT_EC$start) %in% 
                                               c(paste(M_DMS_SEX_TRT_CC$chr,M_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_CC$chr,F_DMS_SEX_TRT_CC$start),
                                                 paste(M_DMS_SEX_TRT_EC$chr,M_DMS_SEX_TRT_EC$start),
                                                 paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start),
                                                 paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start),
                                                 paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start),
                                                 paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start)),]

M_DMS_SEX_TRT_EC_exclusive<-M_DMS_SEX_TRT_EC[!paste(M_DMS_SEX_TRT_EC$chr,M_DMS_SEX_TRT_EC$start) %in% 
                                               c(paste(M_DMS_SEX_TRT_CC$chr,M_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_EC$chr,F_DMS_SEX_TRT_EC$start),
                                                 paste(F_DMS_SEX_TRT_CC$chr,F_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start),
                                                 paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start),
                                                 paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start),
                                                 paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start)),]

F_DMS_SEX_TRT_EE_exclusive<-F_DMS_SEX_TRT_EE[!paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start) %in% 
                                               c(paste(M_DMS_SEX_TRT_CC$chr,M_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_EC$chr,F_DMS_SEX_TRT_EC$start),
                                                 paste(M_DMS_SEX_TRT_EC$chr,M_DMS_SEX_TRT_EC$start),
                                                 paste(F_DMS_SEX_TRT_CC$chr,F_DMS_SEX_TRT_CC$start),
                                                 paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start),
                                                 paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start),
                                                 paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start)),]

M_DMS_SEX_TRT_EE_exclusive<-M_DMS_SEX_TRT_EE[!paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start) %in% 
                                               c(paste(M_DMS_SEX_TRT_CC$chr,M_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_EC$chr,F_DMS_SEX_TRT_EC$start),
                                                 paste(M_DMS_SEX_TRT_EC$chr,M_DMS_SEX_TRT_EC$start),
                                                 paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start),
                                                 paste(F_DMS_SEX_TRT_CC$chr,F_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start),
                                                 paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start)),]

F_DMS_SEX_TRT_CE_exclusive<-F_DMS_SEX_TRT_CE[!paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start) %in% 
                                               c(paste(M_DMS_SEX_TRT_CC$chr,M_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_EC$chr,F_DMS_SEX_TRT_EC$start),
                                                 paste(M_DMS_SEX_TRT_EC$chr,M_DMS_SEX_TRT_EC$start),
                                                 paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start),
                                                 paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start),
                                                 paste(F_DMS_SEX_TRT_CC$chr,F_DMS_SEX_TRT_CC$start),
                                                 paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start)),]

M_DMS_SEX_TRT_CE_exclusive<-M_DMS_SEX_TRT_CE[!paste(M_DMS_SEX_TRT_CE$chr,M_DMS_SEX_TRT_CE$start) %in% 
                                               c(paste(M_DMS_SEX_TRT_CC$chr,M_DMS_SEX_TRT_CC$start),
                                                 paste(F_DMS_SEX_TRT_EC$chr,F_DMS_SEX_TRT_EC$start),
                                                 paste(M_DMS_SEX_TRT_EC$chr,M_DMS_SEX_TRT_EC$start),
                                                 paste(F_DMS_SEX_TRT_EE$chr,F_DMS_SEX_TRT_EE$start),
                                                 paste(M_DMS_SEX_TRT_EE$chr,M_DMS_SEX_TRT_EE$start),
                                                 paste(F_DMS_SEX_TRT_CE$chr,F_DMS_SEX_TRT_CE$start),
                                                 paste(F_DMS_SEX_TRT_CC$chr,F_DMS_SEX_TRT_CC$start)),]


DMS_for_annot_cc_get_dif_F<-paste(F_DMS_SEX_TRT_CC_exclusive$chr, F_DMS_SEX_TRT_CC_exclusive$start,sep = " ")
DMS_for_annot_cc_get_dif_M<-paste(M_DMS_SEX_TRT_CC_exclusive$chr, M_DMS_SEX_TRT_CC_exclusive$start,sep = " ")
DMS_for_annot_ec_get_dif_F<-paste(F_DMS_SEX_TRT_EC_exclusive$chr, F_DMS_SEX_TRT_EC_exclusive$start,sep = " ")
DMS_for_annot_ec_get_dif_M<-paste(M_DMS_SEX_TRT_EC_exclusive$chr, M_DMS_SEX_TRT_EC_exclusive$start,sep = " ")
DMS_for_annot_ee_get_dif_F<-paste(F_DMS_SEX_TRT_EE_exclusive$chr, F_DMS_SEX_TRT_EE_exclusive$start,sep = " ")
DMS_for_annot_ee_get_dif_M<-paste(M_DMS_SEX_TRT_EE_exclusive$chr, M_DMS_SEX_TRT_EE_exclusive$start,sep = " ")
DMS_for_annot_ce_get_dif_F<-paste(F_DMS_SEX_TRT_CE_exclusive$chr, F_DMS_SEX_TRT_CE_exclusive$start,sep = " ")
DMS_for_annot_ce_get_dif_M<-paste(M_DMS_SEX_TRT_CE_exclusive$chr, M_DMS_SEX_TRT_CE_exclusive$start,sep = " ")

##################################################################################################################################
df_cc_f <- data.frame(chr=sapply(strsplit(DMS_for_annot_cc_get_dif_F, " "), `[`, 1), 
                      start=sapply(strsplit(DMS_for_annot_cc_get_dif_F, " "), `[`, 2),
                      end=sapply(strsplit(DMS_for_annot_cc_get_dif_F, " "), `[`, 2)) 
# get annotation
anot_cc_get_dif_f <- getAnnotationFun(makeGRangesFromDataFrame(df_cc_f))

#253 Genes
nrow(anot_cc_get_dif_f)

anot_cc_get_dif_f <- anot_cc_get_dif_f %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

topgenes_cc_f<-anot_cc_get_dif_f[anot_cc_get_dif_f$nCpGperGenekb>1,]
topgenes_cc_f_summar<-getGeneSummary(topgenes_cc_f)
##################################################################################################################################

df_cc_m <- data.frame(chr=sapply(strsplit(DMS_for_annot_cc_get_dif_M, " "), `[`, 1), 
                      start=sapply(strsplit(DMS_for_annot_cc_get_dif_M, " "), `[`, 2),
                      end=sapply(strsplit(DMS_for_annot_cc_get_dif_M, " "), `[`, 2)) 
# get annotation
anot_cc_get_dif_m <- getAnnotationFun(makeGRangesFromDataFrame(df_cc_m))

#211 Genes
nrow(anot_cc_get_dif_m)

anot_cc_get_dif_m <- anot_cc_get_dif_m %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

topgenes_cc_m<-anot_cc_get_dif_m[anot_cc_get_dif_m$nCpGperGenekb>1,]
topgenes_cc_m_summar<-getGeneSummary(topgenes_cc_m)
##################################################################################################################################

df_ec_f <- data.frame(chr=sapply(strsplit(DMS_for_annot_ec_get_dif_F, " "), `[`, 1), 
                      start=sapply(strsplit(DMS_for_annot_ec_get_dif_F, " "), `[`, 2),
                      end=sapply(strsplit(DMS_for_annot_ec_get_dif_F, " "), `[`, 2)) 
# get annotation
anot_ec_get_dif_f <- getAnnotationFun(makeGRangesFromDataFrame(df_ec_f))

#198 Genes
nrow(anot_ec_get_dif_f)

anot_ec_get_dif_f <- anot_ec_get_dif_f %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

topgenes_ec_f<-anot_ec_get_dif_f[anot_ec_get_dif_f$nCpGperGenekb>1,]
topgenes_ec_f_summar<-getGeneSummary(topgenes_ec_f)
##################################################################################################################################

df_ec_m <- data.frame(chr=sapply(strsplit(DMS_for_annot_ec_get_dif_M, " "), `[`, 1), 
                      start=sapply(strsplit(DMS_for_annot_ec_get_dif_M, " "), `[`, 2),
                      end=sapply(strsplit(DMS_for_annot_ec_get_dif_M, " "), `[`, 2)) 
# get annotation
anot_ec_get_dif_m <- getAnnotationFun(makeGRangesFromDataFrame(df_ec_m))

#385 Genes
nrow(anot_ec_get_dif_m)

anot_ec_get_dif_m <- anot_ec_get_dif_m %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

topgenes_ec_m<-anot_ec_get_dif_m[anot_ec_get_dif_m$nCpGperGenekb>1,]
topgenes_ec_m_summar<-getGeneSummary(topgenes_ec_m)
##################################################################################################################################

df_ee_f <- data.frame(chr=sapply(strsplit(DMS_for_annot_ee_get_dif_F, " "), `[`, 1), 
                      start=sapply(strsplit(DMS_for_annot_ee_get_dif_F, " "), `[`, 2),
                      end=sapply(strsplit(DMS_for_annot_ee_get_dif_F, " "), `[`, 2)) 
# get annotation
anot_ee_get_dif_f <- getAnnotationFun(makeGRangesFromDataFrame(df_ee_f))

#385 Genes
nrow(anot_ee_get_dif_f)

anot_ee_get_dif_f <- anot_ee_get_dif_f %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

topgenes_ee_f<-anot_ee_get_dif_f[anot_ee_get_dif_f$nCpGperGenekb>1,]
topgenes_ee_f_summar<-getGeneSummary(topgenes_ee_f)
##################################################################################################################################

df_ee_m <- data.frame(chr=sapply(strsplit(DMS_for_annot_ee_get_dif_M, " "), `[`, 1), 
                      start=sapply(strsplit(DMS_for_annot_ee_get_dif_M, " "), `[`, 2),
                      end=sapply(strsplit(DMS_for_annot_ee_get_dif_M, " "), `[`, 2)) 
# get annotation
anot_ee_get_dif_m <- getAnnotationFun(makeGRangesFromDataFrame(df_ee_m))

#459 Genes
nrow(anot_ee_get_dif_m)

anot_ee_get_dif_m <- anot_ee_get_dif_m %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

topgenes_ee_m<-anot_ee_get_dif_m[anot_ee_get_dif_m$nCpGperGenekb>1,]
topgenes_ee_m_summar<-getGeneSummary(topgenes_ee_m)
##################################################################################################################################

df_ce_f <- data.frame(chr=sapply(strsplit(DMS_for_annot_ce_get_dif_F, " "), `[`, 1), 
                      start=sapply(strsplit(DMS_for_annot_ce_get_dif_F, " "), `[`, 2),
                      end=sapply(strsplit(DMS_for_annot_ce_get_dif_F, " "), `[`, 2)) 
# get annotation
anot_ce_get_dif_f <- getAnnotationFun(makeGRangesFromDataFrame(df_ce_f))

#324 Genes
nrow(anot_ce_get_dif_f)

anot_ce_get_dif_f <- anot_ce_get_dif_f %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

topgenes_ce_f<-anot_ce_get_dif_f[anot_ce_get_dif_f$nCpGperGenekb>1,]
topgenes_ce_f_summar<-getGeneSummary(topgenes_ce_f)
##################################################################################################################################

df_ce_m <- data.frame(chr=sapply(strsplit(DMS_for_annot_ce_get_dif_M, " "), `[`, 1), 
                      start=sapply(strsplit(DMS_for_annot_ce_get_dif_M, " "), `[`, 2),
                      end=sapply(strsplit(DMS_for_annot_ce_get_dif_M, " "), `[`, 2)) 
# get annotation
anot_ce_get_dif_m <- getAnnotationFun(makeGRangesFromDataFrame(df_ce_m))

#363 Genes
nrow(anot_ce_get_dif_m)

anot_ce_get_dif_m <- anot_ce_get_dif_m %>% 
  mutate(chr = gsub("Gy_chr", "", seqid), chrom_nr = seqid %>% deroman(), chrom_order=factor(chrom_nr) %>% 
           as.numeric()) %>% arrange(chrom_order)  %>% 
  mutate(geneLengthkb = (end - start)/1000, nCpGperGenekb = nCpG/geneLengthkb)

topgenes_ce_m<-anot_ce_get_dif_m[anot_ce_get_dif_m$nCpGperGenekb>1,]
topgenes_ce_m_summar<-getGeneSummary(topgenes_ce_m)
##################################################################################################################################



df_GO_cc_f<-makedfGO(annot_trt_get_diff = anot_cc_get_dif_f, gene_universe = gene_universe)
df_GO_cc_m<-makedfGO(annot_trt_get_diff = anot_cc_get_dif_m, gene_universe = gene_universe)
df_GO_ec_f<-makedfGO(annot_trt_get_diff = anot_ec_get_dif_f, gene_universe = gene_universe)
df_GO_ec_m<-makedfGO(annot_trt_get_diff = anot_ec_get_dif_m, gene_universe = gene_universe)
df_GO_ee_f<-makedfGO(annot_trt_get_diff = anot_ee_get_dif_f, gene_universe = gene_universe)
df_GO_ee_m<-makedfGO(annot_trt_get_diff = anot_ee_get_dif_m, gene_universe = gene_universe)
df_GO_ce_f<-makedfGO(annot_trt_get_diff = anot_ce_get_dif_f, gene_universe = gene_universe)
df_GO_ce_m<-makedfGO(annot_trt_get_diff = anot_ce_get_dif_m, gene_universe = gene_universe)
df_GO_cc_f$trtG1G2<-"CONTROL-CONTROL_female"
df_GO_cc_m$trtG1G2<-"CONTROL-CONTROL_male"
df_GO_ec_f$trtG1G2<-"EXPOSED-CONTROL_female"
df_GO_ec_m$trtG1G2<-"EXPOSED-CONTROL_male"
df_GO_ee_f$trtG1G2<-"EXPOSED-EXPOSED_female"
df_GO_ee_m$trtG1G2<-"EXPOSED-EXPOSED_male"
df_GO_ce_f$trtG1G2<-"CONTROL-EXPOSED_female"
df_GO_ce_m$trtG1G2<-"CONTROL-EXPOSED_male"

dfGO_sex <- rbind(df_GO_cc_f,df_GO_cc_m, df_GO_ec_f, df_GO_ec_m, df_GO_ee_f,df_GO_ee_m,df_GO_ce_f,df_GO_ce_m)

dfGO_Sex_trt<-dfGO_sex %>% ggplot(aes(x=trtG1G2,y = factor(GO.name))) +
  geom_point(aes(color = p.value.adjusted, size = genePercent)) +
  scale_color_gradient(name="adjusted\np-value", low = "red", high = "blue") +
  scale_size_continuous(name = "% of genes")+
  theme_bw() + ylab("") + xlab("Sex Differences within each treatment") +
  theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"), 
        legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"), 
        legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), legend.position="left",
        plot.margin = margin(1,1,1.5,1.2, "cm")) + # grey box for legend
  facet_grid((GO.category)~trtG1G2, scales="free",space = "free") 

pdf(file = "Figures_for_diss/dfGO_Sex_trt.pdf",width = 12 ,height = 45)

dfGO_Sex_trt

dev.off()
