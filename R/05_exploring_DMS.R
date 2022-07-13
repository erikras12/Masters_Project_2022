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




annot_Trt <- as.data.frame(TrtAnn@members)
Trt_Diff_man_Plot<-makeManhattanPlots(DMSfile = myDiff25p_Trt, annotFile = annot_Trt, GYgynogff = GYgynogff, 
                                      mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of Trt Diff")
Trt_Diff_man_Plot


annot_Sex_Trt <- as.data.frame(Sex_Trt_ANN@members)
Sex_Trt_Diff_man_Plot<-makeManhattanPlots(DMSfile = myDiff25p_Sex_Trt, annotFile = annot_Sex_Trt, GYgynogff = GYgynogff, 
                                          mycols = c("red", "grey", "black", "green"), mytitle = "Manhattan plot of Sex+Trt Diff")
Sex_Trt_Diff_man_Plot

unite_for_DMS_Sex_Trt<-readRDS("R_shit/unite_for_DMS_Sex_Trt.RDS")

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


###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
myDiff25p_Sex$
  
  library(ggVennDiagram)
x = list(DMS_Sex_COMP = paste(myDiff25p_Sex$chr,myDiff25p_Sex$start),
         DMS_Sex_Trt_COMP = paste(myDiff25p_Sex_Trt$chr,myDiff25p_Sex_Trt$start),
         DMS_Trt_COMP = paste(myDiff25p_Trt$chr,myDiff25p_Trt$start))
ggVennDiagram(x) + scale_fill_gradient(low = "blue", high = "red")

#correction for difference in statring data sets
a = paste(DMS_Sex$chr,DMS_Sex$start)
b = paste(DMS_Trt$chr,DMS_Trt$start)
c = paste(DMS_Sex_Trt$chr,DMS_Sex_Trt$start)

pos_int<-Reduce(intersect, list(a,b,c))


a1 = paste(myDiff25p_Sex$chr,myDiff25p_Sex$start)
b1 = paste(myDiff25p_Trt$chr,myDiff25p_Trt$start)
c1 = paste(myDiff25p_Sex_Trt$chr,myDiff25p_Sex_Trt$start)



x = list(DMS_Trt_COMP = b1[b1 %in% pos_int],
         DMS_Sex_Trt_COMP = c1[c1 %in% pos_int])
ggVennDiagram(x) + scale_fill_gradient(low = "blue", high = "red")





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
