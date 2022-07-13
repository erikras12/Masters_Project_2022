DMS_Sex<-readRDS("R_data/DMS_Sex.RDS")
DMS_Trt<-readRDS("R_data/DMS_Trt.RDS")
DMS_Sex_Trt<-readRDS("R_data/DMS_Sex_Trt.RDS")
load("R_data/myMetaData.RDS")


#read annotation files
annotBed12= readTranscriptFeatures("R_shit/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.transdec.bed12",remove.unusual = FALSE, up.flank = 1500, down.flank = 500)
GYgynogff = read.table("R_shit/Gy_allnoM_rd3.maker_apocrita.noseq_corrected_chromoAndLength.txt")
names(GYgynogff)[1]<-"chrom"
names(GYgynogff)[2]<-"length"
annotGff3 = readGFF("R_shit/Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.gff")

#obtain the significant difference with a specified methyaltion difference
myDiff25p_Sex=getMethylDiff(DMS_Sex,difference=25,qvalue=0.01)
myDiff25p_Trt=getMethylDiff(DMS_Trt,difference=25,qvalue=0.01)
myDiff25p_Sex_Trt=getMethylDiff(DMS_Sex_Trt,difference=45,qvalue=0.01)

#find where the DMS' are in contect to the genome (Promoter, intron, exon, intergenic)
SexANN = annotateWithGeneParts(as(myDiff25p_Sex,"GRanges"),annotBed12)
TrtAnn = annotateWithGeneParts(as(myDiff25p_Trt,"GRanges"),annotBed12)
Sex_Trt_ANN = annotateWithGeneParts(as(myDiff25p_Sex_Trt,"GRanges"),annotBed12)

#pie chart visualising information above
Sex_specific_CpG_Context<-plotTargetAnnotation(SexANN,precedence=TRUE, main="Sex", cex.legend = 0.5, border="white")
Sex_specific_CpG_Context<-plotTargetAnnotation(TrtAnn,precedence=TRUE, main="Sex", cex.legend = 0.5, border="white")
Sex_specific_CpG_Context<-plotTargetAnnotation(Sex_Trt_ANN,precedence=TRUE, main="Sex", cex.legend = 0.5, border="white")


#manhattan plot for DMS between the sexes in G2 
annot_Sex <- as.data.frame(SexANN@members)
Sex_Diff_man_Plot<-makeManhattanPlots(DMSfile = myDiff25p_Sex, annotFile = annot_Sex, GYgynogff = GYgynogff, 
                                      mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of Sex Diff")
Sex_Diff_man_Plot


#look at the DMS with the largest differences and plot them 
outliers_Sex = which(abs(myDiff25p_Sex$meth.diff)>25 + 2*sd(abs(myDiff25p_Sex$meth.diff)))

outliers_annot_Sex <- as.data.frame(SexANN@members)[outliers_Sex,]

Sex_biggest_Diff_man_Plot<-makeManhattanPlots(DMSfile = myDiff25p_Sex[outliers_Sex,], annotFile = outliers_annot_Sex, GYgynogff = GYgynogff, 
                                              mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of biggest Sex Diff")

Sex_biggest_Diff_man_Plot

#manhattan Plot showing DMS patterns accross the genome
#the data for trt treatment had more than 2 factors which meant that it was unable to be visualised via a manhattan Plot same with the Sex_trt data 

#load unite object for Trt and Sex_trt
unite_for_DMS_Trt_w_sx<-readRDS("R_data/unite_for_DMS_Trt.RDS")
unite_for_DMS_Trt<- unite_for_DMS_Trt_w_sx[!unite_for_DMS_Trt_w_sx$chr %in% c("Gy_chrXIX","Gy_chrUn"),]

unite_for_DMS_Sex_Trt_w_sx <- readRDS("R_data/unite_for_DMS_Sex_Trt.RDS")
unite_for_DMS_Sex_Trt<- unite_for_DMS_Sex_Trt_w_sx[!unite_for_DMS_Sex_Trt_w_sx$chr %in% c("Gy_chrXIX","Gy_chrUn"),]

#subset the unite object for both treatments by the samples which show significant differences in their respective mydiff objects
methyl_unite_rows_selected_trt<-subset_unite(myDiff25p_Trt,unite_for_DMS_Trt)
methyl_unite_rows_selected_Sex_trt<-subset_unite(myDiff25p_Sex_Trt,unite_for_DMS_Sex_Trt)
#these object can be used to again subset 
##################################################################################################################################################
##### subsetting the unite for G1G2 by the DMS for G1G2 treatment into its 4 factors
##################################################################################################################################################

#run percmet_for_trt to get average meth for each CpG (in this case for trt 1 = e_control)
e_control<-percMeth_for_trt(methyl_unite_rows_selected_trt,trt = 1)
e_control_df<-data.frame(e_control)
e_control_obj<-methylKit::reorganize(methyl_unite_rows_selected_trt, treatment = myMetaData$trtG1G2_NUM[myMetaData$trtG1G2_NUM %in% 1],
                              sample.ids = myMetaData$SampleID[myMetaData$trtG1G2_NUM %in% 1])
#add chr column and pos column for manhatton plot
e_control_df$chr <- e_control_obj$chr
e_control_df$pos <- e_control_obj$start
#name the methylation column avg.meth
colnames(e_control_df)[1] <- "avg.meth"


e_exposed<-percMeth_for_trt(methyl_unite_rows_selected_trt,trt = 2)
e_exposed_df<-data.frame(e_exposed)
e_exposed_obj<-methylKit::reorganize(methyl_unite_rows_selected_trt, treatment = myMetaData$trtG1G2_NUM[myMetaData$trtG1G2_NUM %in% 2],
                              sample.ids = myMetaData$SampleID[myMetaData$trtG1G2_NUM %in% 2])
e_exposed_df$chr <- e_exposed_obj$chr
e_exposed_df$pos <- e_exposed_obj$start
colnames(e_exposed_df)[1] <- "avg.meth"



ne_control<-percMeth_for_trt(methyl_unite_rows_selected_trt,trt = 3)
ne_control_df<-data.frame(ne_control)
ne_control_obj<-methylKit::reorganize(methyl_unite_rows_selected_trt, treatment = myMetaData$trtG1G2_NUM[myMetaData$trtG1G2_NUM %in% 3],
                              sample.ids = myMetaData$SampleID[myMetaData$trtG1G2_NUM %in% 3])
ne_control_df$chr <- ne_control_obj$chr
ne_control_df$pos <- ne_control_obj$start
colnames(ne_control_df)[1] <- "avg.meth"



ne_exposed<-percMeth_for_trt(methyl_unite_rows_selected_trt,trt = 4)
ne_exposed_df<-data.frame(ne_exposed)
ne_exposed_obj<-methylKit::reorganize(methyl_unite_rows_selected_trt, treatment = myMetaData$trtG1G2_NUM[myMetaData$trtG1G2_NUM %in% 4],
                              sample.ids = myMetaData$SampleID[myMetaData$trtG1G2_NUM %in% 4])
ne_exposed_df$chr <- ne_exposed_obj$chr
ne_exposed_df$pos <- ne_exposed_obj$start
colnames(ne_exposed_df)[1] <- "avg.meth"


#we can observe these data points on a manhattan plot now showing the percentage methylation of each DMS for each treatment

annot_Trt <- as.data.frame(TrtAnn@members)

##PLOT FOR TRT 1 (EXPOSED G1 - CONTROL G2)
TRT_1_mnplot<-makeManhattanPlots_mult_trt(DMSfile = e_control_df, annotFile = annot_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of e_control individuals")

TRT_2_mnplot<-makeManhattanPlots_mult_trt(DMSfile = e_exposed_df, annotFile = annot_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of e_exposed individuals")

TRT_3_mnplot<-makeManhattanPlots_mult_trt(DMSfile = ne_control_df, annotFile = annot_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of ne_control individuals")

TRT_4_mnplot<-makeManhattanPlots_mult_trt(DMSfile = ne_exposed_df, annotFile = annot_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of ne_exposed individuals")

grid.arrange(TRT_1_mnplot, TRT_2_mnplot,TRT_3_mnplot,TRT_4_mnplot, ncol=1)



## Preparing Sex_trt groups for manhatton Plot
### Exposed fater group (M&F)

F_E_control<-percMeth_for_trt_sex(methyl_unite_rows_selected_Sex_trt,trt = 1)



F_E_control_df<-data.frame(F_E_control)
F_E_control_obj<-methylKit::reorganize(methyl_unite_rows_selected_Sex_trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 1],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 1])
F_E_control_df$chr <- F_E_control_obj$chr
F_E_control_df$pos <- F_E_control_obj$start
colnames(F_E_control_df)[1] <- "avg.meth"



M_E_control<-percMeth_for_trt_sex(methyl_unite_rows_selected_Sex_trt,trt = 2)
M_E_control_df<-data.frame(M_E_control)
M_E_control_obj<-methylKit::reorganize(methyl_unite_rows_selected_Sex_trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 2],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 2])
M_E_control_df$chr <- M_E_control_obj$chr
M_E_control_df$pos <- M_E_control_obj$start
colnames(M_E_control_df)[1] <- "avg.meth"


F_E_exposed<-percMeth_for_trt_sex(methyl_unite_rows_selected_Sex_trt,trt = 3)
F_E_exposed_df<-data.frame(F_E_exposed)
F_E_exposed_obj<-methylKit::reorganize(methyl_unite_rows_selected_Sex_trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 3],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 3])
F_E_exposed_df$chr <-F_E_exposed_obj$chr
F_E_exposed_df$pos <-F_E_exposed_obj$start
colnames(F_E_exposed_df)[1] <- "avg.meth"


M_E_exposed<-percMeth_for_trt_sex(methyl_unite_rows_selected_Sex_trt,trt = 4)
M_E_exposed_df<-data.frame(M_E_exposed)
M_E_exposed_obj<-methylKit::reorganize(methyl_unite_rows_selected_Sex_trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 4],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 4])
M_E_exposed_df$chr <-M_E_exposed_obj$chr
M_E_exposed_df$pos <-M_E_exposed_obj$start
colnames(M_E_exposed_df)[1] <- "avg.meth"

annot_Sex_Trt <- as.data.frame(Sex_Trt_ANN@members)

##PLOT FOR TRT 1 Female(EXPOSED G1 - CONTROL G2)
Sex_Trt_1<-makeManhattanPlots_mult_trt(DMSfile = F_E_control_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of e_control Females")


##PLOT FOR TRT 2 male(EXPOSED G1 - CONTROL G2)
Sex_Trt_2<-makeManhattanPlots_mult_trt(DMSfile = M_E_control_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of e_control Males")


##PLOT FOR TRT 3 female(EXPOSED G1 - exposed G2)
Sex_Trt_3<-makeManhattanPlots_mult_trt(DMSfile =F_E_exposed_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of e_exposed Females")


##PLOT FOR TRT 4 male(EXPOSED G1 - Exposed G2)
Sex_Trt_4<-makeManhattanPlots_mult_trt(DMSfile = M_E_exposed_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of e_exposed Males")



### Control father group (M&F)



F_NE_control<-percMeth_for_trt_sex(methyl_unite_rows_selected_Sex_trt,trt = 5)
F_NE_control_df<-data.frame(F_NE_control)
F_NE_control_obj<-methylKit::reorganize(methyl_unite_rows_selected_Sex_trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 5],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 5])
F_NE_control_df$chr <- F_NE_control_obj$chr
F_NE_control_df$pos <- F_NE_control_obj$start
colnames(F_NE_control_df)[1] <- "avg.meth"


M_NE_control<-percMeth_for_trt_sex(methyl_unite_rows_selected_Sex_trt,trt = 6)
M_NE_control_df<-data.frame(M_NE_control)
M_NE_control_obj<-methylKit::reorganize(methyl_unite_rows_selected_Sex_trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 6],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 6])
M_NE_control_df$chr <- M_NE_control_obj$chr
M_NE_control_df$pos <- M_NE_control_obj$start
colnames(M_NE_control_df)[1] <- "avg.meth"


F_NE_exposed<-percMeth_for_trt_sex(methyl_unite_rows_selected_Sex_trt,trt = 7)
F_NE_exposed_df<-data.frame(F_NE_exposed)
F_NE_exposed_obj<-methylKit::reorganize(methyl_unite_rows_selected_Sex_trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 7],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 7])
F_NE_exposed_df$chr <-F_NE_exposed_obj$chr
F_NE_exposed_df$pos <-F_NE_exposed_obj$start
colnames(F_NE_exposed_df)[1] <- "avg.meth"


M_NE_exposed<-percMeth_for_trt_sex(methyl_unite_rows_selected_Sex_trt,trt = 8)
M_NE_exposed_df<-data.frame(M_NE_exposed)
M_NE_exposed_obj<-methylKit::reorganize(methyl_unite_rows_selected_Sex_trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% 8],
                              sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% 8])
M_NE_exposed_df$chr <-M_NE_exposed_obj$chr
M_NE_exposed_df$pos <-M_NE_exposed_obj$start
colnames(M_NE_exposed_df)[1] <- "avg.meth"




##PLOT FOR TRT 5 female(NOT EXPOSED G1 - CONTROL G2)
Sex_Trt_5<-makeManhattanPlots_mult_trt(DMSfile = F_NE_control_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of NE_control Females")




##PLOT FOR TRT 6 male(NOT EXPOSED G1 - CONTROL G2)
Sex_Trt_6<-makeManhattanPlots_mult_trt(DMSfile = M_NE_control_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of NE_control Males")



##PLOT FOR TRT 7 female(NOT EXPOSED G1 - EXPOSED G2)
Sex_Trt_7<-makeManhattanPlots_mult_trt(DMSfile = F_NE_exposed_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of NE_exposed Females")





##PLOT FOR TRT 7 male(NOT EXPOSED G1 - EXPOSED G2)
Sex_Trt_8<-makeManhattanPlots_mult_trt(DMSfile = M_NE_exposed_df, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of methylation values of NE_exposed males")




grid.arrange(Sex_Trt_1,Sex_Trt_2,Sex_Trt_3,Sex_Trt_4,Sex_Trt_5,Sex_Trt_6,Sex_Trt_7,Sex_Trt_8, ncol=1)



###################################################################################################################################################################
##### Venn #######################
###################################################################################################################################################################

#trying to identify how many DMS' treatments have in common 


#subset unite for DMS SEX and Trt into the 4 treatment groups 
CC<-methylKit::reorganize(unite_for_DMS_Sex_Trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(5, 6) ],
                                                      sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(5, 6)])
CC@treatment

EE<-methylKit::reorganize(unite_for_DMS_Sex_Trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(3, 4) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(3, 4)])
EE@treatment

CE<-methylKit::reorganize(unite_for_DMS_Sex_Trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(7, 8) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(7, 8)])
CE@treatment

EC<-methylKit::reorganize(unite_for_DMS_Sex_Trt, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(1, 2) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(1, 2)])
EC@treatment



CC<-methylKit::reorganize(unite_g2_in_5_no_sxChr, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(5, 6) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(5, 6)])
CC@treatment

CE<-methylKit::reorganize(unite_g2_in_5_no_sxChr, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(7, 8) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(7, 8)])
CE@treatment

EE<-methylKit::reorganize(unite_g2_in_5_no_sxChr, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(3, 4) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(3, 4)])
EE@treatment


EC<-methylKit::reorganize(unite_g2_in_5_no_sxChr, treatment = myMetaData$Sex_trt[myMetaData$Sex_trt %in% c(1, 2) ],
                          sample.ids = myMetaData$SampleID[myMetaData$Sex_trt %in% c(1, 2)])
EC@treatment


#create covariates df to take into account genetic variation
CC_samples<-CC@sample.ids
CC_covariates <- as.data.frame(myMetaData$brotherPairID[ myMetaData$SampleID %in% CC_samples] )


cc_dif<-calculateDiffMeth(CC,covariates = covariates,mc.cores = 3)

cc_get_dif<-getMethylDiff(cc_diff,difference = 25,qvalue = 0.01)

saveRDS(cc_dif,file = "cc_diff.RDS")
cc_diff<-readRDS("R_data/cc_diff.RDS")



CE_samples<-CE@sample.ids
CE_covariates <- as.data.frame(myMetaData$brotherPairID[ myMetaData$SampleID %in% CE_samples] )

ce_dif<-calculateDiffMeth(CE,covariates = CE_covariates,mc.cores = 3)

ce_get_dif<-getMethylDiff(ce_dif,difference = 25,qvalue = 0.01)

ce_dif<-saveRDS(ce_dif,file = "ce_diff.RDS")
ce_dif<-readRDS("R_data/ce_diff.RDS")


EE_samples<-EE@sample.ids
EE_covariates <- as.data.frame(myMetaData$brotherPairID[ myMetaData$SampleID %in% EE_samples] )

ee_dif<-calculateDiffMeth(EE,covariates = EE_covariates,mc.cores = 3)

ee_get_dif<-getMethylDiff(ee_dif,difference = 25,qvalue = 0.01)

saveRDS(ee_dif,file = "ee_diff.RDS")

ee_dif<-readRDS("R_data/ee_diff.RDS")


EC_samples<-EC@sample.ids
EC_covariates <- as.data.frame(myMetaData$brotherPairID[ myMetaData$SampleID %in% EC_samples] )

ec_dif<-calculateDiffMeth(EC,covariates = EC_covariates,mc.cores = 3)

ec_get_dif<-getMethylDiff(ec_dif,difference = 25,qvalue = 0.01)


saveRDS(ec_dif,file = "ec_diff.RDS")
ec_dif<-readRDS("R_data/ec_diff.RDS")


library(ggVennDiagram)
listx = list(Exposed_Exposed = paste(ee_get_dif$chr,ee_get_dif$start),
             Control_Exposed = paste(ce_get_dif$chr,ce_get_dif$start),
             Exposed_Control = paste(ec_get_dif$chr,ec_get_dif$start),
             Control_Control = paste(cc_get_dif$chr,cc_get_dif$start))
ggVennDiagram(listx) + scale_fill_gradient(low = "blue", high = "red")


list_ce_ec = list(Control_Exposed = paste(ce_get_dif$chr,ce_get_dif$start),
                  Exposed_Control = paste(ec_get_dif$chr,ec_get_dif$start))
ggVennDiagram(list_ce_ec) + scale_fill_gradient(low = "blue", high = "red")


list_ee_ec = list(Exposed_Exposed = paste(ee_get_dif$chr,ee_get_dif$start),
                 Exposed_Control = paste(ec_get_dif$chr,ec_get_dif$start))
ggVennDiagram(list_ee_ec) + scale_fill_gradient(low = "blue", high = "red")
