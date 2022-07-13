###DMS/DMR analysis###############################################################
##################################################################################
load("/data/SBCS-EizaguirreLab/Eri/data/coverage/metadata.RData")
load("/data/SBCS-EizaguirreLab/Eri/data/coverage/for_DMS.RData")

#combine trt with sex into new column giving you 8 different levels of sex_trt
myMetaData$Sex_trt<- paste(myMetaData$trtG1G2,myMetaData$Sex)

#factor of 8
myMetaData$Sex_trt <-as.numeric(as.factor(myMetaData$Sex_trt))



unite_G2_in_14_nosxchr@treatment <- myMetaData$Sex_trt
G2_cov_no_sxChr@treatment <- myMetaData$trtG1G2_NUM

covariates <- as.data.frame(myMetaData$brotherPairID)


############################################################################################
dataPath_111 = "/data/SBCS-EizaguirreLab/Eri/data/coverage/111_samples/"

# list of paths to all coverage files
temp_111 = list.files(path = dataPath_111,
                  pattern = ".sam.bismark.cov.gz",
                  full.names = T)

myObj_sex_trt = methylKit::methRead(as.list(temp_111),
                            mincov = 10,
                            pipeline = "bismarkCoverage",
                            sample.id = as.list(myMetaData$ID),
                            assembly = "Gynogen_pchrom_assembly_all",
                            treatment = myMetaData$Sex_trt,
                            context = "CpG")

#unite_g2_in_5=methylKit::unite(myObj_sex_trt,min.per.group = 5L, destrand=FALSE,mc.cores = 3)


table(myObj_sex_trt@treatment)




#filter myObj to remove bias from PCA
unite_g2_in_5_filt = methylKit::filterByCoverage(myObj_sex_trt,lo.count=10,lo.perc=NULL,
                                                hi.count=NULL,hi.perc=99.9)
#normalise coverage between samples
normalised_filtered_unite_g2_in_5= methylKit::normalizeCoverage(unite_g2_in_5_filt)

all_samples_normalised_filtered_unite_g2_in_5 = methylKit::unite(normalised_filtered_unite_g2_in_5, min.per.group = 5L ,mc.cores = 3)

unite_g2_in_5_no_sxChr = all_samples_normalised_filtered_unite_g2_in_5[!all_samples_normalised_filtered_unite_g2_in_5$chr %in% c("Gy_chrXIX","Gy_chrUn"),]

saveRDS(all_samples_normalised_filtered_unite_g2_in_5, file = "unite_for_DMS_Sex_Trt.RDS")

unite_g2_in_5_no_sxChr = unite_for_DMS_Sex_Trt[!unite_for_DMS_Sex_Trt$chr %in% c("Gy_chrXIX","Gy_chrUn"),]

myDiff_2=calculateDiffMeth(unite_for_DMS_Sex_Trt, covariates = covariates, mc.cores = 3)

saveRDS(myDiff_2, file = "DMS_Sex_Trt.RDS")

myDiff25p_2=getMethylDiff(myDiff_2,difference=25,qvalue=0.01)

annotBed12= readTranscriptFeatures("Gy_allnoM_rd3.maker_apocrita.noseq_corrected.gff.streamlined_for_AGAT.CURATED.transdec.bed12",remove.unusual = FALSE, up.flank = 1500, down.flank = 500)

diffAnn_PAR = annotateWithGeneParts(as(myDiff25p_2,"GRanges"),annotBed12)
diffAnn_PAR

unite_for_DMS_Sex_Trt<-readRDS("unite_for_DMS_Sex_Trt.RDS")

unite_for_DMS_Sex_Trt_no_sxchr = unite_for_DMS_Sex_Trt[!unite_for_DMS_Sex_Trt$chr %in% c("Gy_chrXIX","Gy_chrUn"),]


tile_sex_trt<-tileMethylCounts(unite_for_DMS_Sex_Trt_no_sxchr,win.size=100,step.size=100,cov.bases = 10)
tile_sex_trt_diff<-calculateDiffMeth(tile_sex_trt,covariates = covariates,mc.cores = 3)
tile_sex_trt_15<-getMethylDiff(tile_sex_trt_diff,difference = 15,qvalue = 0.01)

table(tile_sex_trt_15$chr)
table(tile_sex_15$chr)
table(tile_trt_15$chr)

save(tile_sex_trt_15,tile_sex_15,tile_trt_15,file = "DMR_objects.RData")
#########################################################################################################################################################################

myObj_sex_ = methylKit::methRead(as.list(temp_111),
                                    mincov = 10,
                                    pipeline = "bismarkCoverage",
                                    sample.id = as.list(myMetaData$ID),
                                    assembly = "Gynogen_pchrom_assembly_all",
                                    treatment = myMetaData$Sex,
                                    context = "CpG")
table(myObj_sex_@treatment)

unite_g2_in_5_sex=methylKit::unite(myObj_sex_,min.per.group = 25L, destrand=FALSE,mc.cores = 3)






#filter myObj to remove bias from PCA
unite_g2_in_5_filt_sex = methylKit::filterByCoverage(myObj_sex_,lo.count=10,lo.perc=NULL,
                                                 hi.count=NULL,hi.perc=99.9)
#normalise coverage between samples
normalised_filtered_unite_g2_in_5_sex= methylKit::normalizeCoverage(unite_g2_in_5_filt_sex)

all_samples_normalised_filtered_unite_g2_in_25_sex = methylKit::unite(normalised_filtered_unite_g2_in_5_sex, min.per.group = 5L ,mc.cores = 3)

saveRDS(all_samples_normalised_filtered_unite_g2_in_25_sex, file = "unite_for_DMS_Sex.RDS")

unite_g2_in_5_no_sxChr_sex = all_samples_normalised_filtered_unite_g2_in_25_sex[!all_samples_normalised_filtered_unite_g2_in_25_sex$chr %in% c("Gy_chrXIX","Gy_chrUn"),]


myDiff_3=calculateDiffMeth(unite_g2_in_5_no_sxChr_sex, covariates = covariates, mc.cores = 3)

saveRDS(myDiff_3, file = "DMS_Sex.RDS")

myDiff25p_3=getMethylDiff(myDiff_3,difference=25,qvalue=0.01)

diffAnn_PAR_3 = annotateWithGeneParts(as(myDiff25p_3,"GRanges"),annotBed12)
diffAnn_PAR_3

unite_for_DMS_Sex<-readRDS("unite_for_DMS_Sex.RDS")

unite_for_DMS_Sex_nosx<-unite_for_DMS_Sex[!unite_for_DMS_Sex$chr %in% c("Gy_chrXIX","Gy_chrUn"),]


tile_sex<-tileMethylCounts(unite_for_DMS_Sex_nosx,win.size=100,step.size=100,cov.bases = 10)
tile_sex_diff<-calculateDiffMeth(tile_sex,covariates = covariates,mc.cores = 3)





tile_sex_15<-getMethylDiff(tile_sex_diff,difference = 15,qvalue = 0.01)

###########################################################################################################################################################################


myMetaData$trtG1G2_NUM <-as.numeric(as.factor(myMetaData$trtG1G2_NUM))

myObj_trt_ = methylKit::methRead(as.list(temp_111),
                                 mincov = 10,
                                 pipeline = "bismarkCoverage",
                                 sample.id = as.list(myMetaData$ID),
                                 assembly = "Gynogen_pchrom_assembly_all",
                                 treatment = myMetaData$trtG1G2_NUM,
                                 context = "CpG")
table(myObj_trt_@treatment)

unite_g2_in_5_trt=methylKit::unite(myObj_trt_,min.per.group = 14L, destrand=FALSE,mc.cores = 3)






#filter myObj to remove bias from PCA
unite_g2_in_14_filt_trt = methylKit::filterByCoverage(myObj_trt_,lo.count=10,lo.perc=NULL,
                                                     hi.count=NULL,hi.perc=99.9)
#normalise coverage between samples
normalised_filtered_unite_g2_in_14_trt= methylKit::normalizeCoverage(unite_g2_in_14_filt_trt)

all_samples_normalised_filtered_unite_g2_in_14_trt = methylKit::unite(normalised_filtered_unite_g2_in_14_trt, min.per.group = 14L ,mc.cores = 3)

saveRDS(all_samples_normalised_filtered_unite_g2_in_14_trt, file = "unite_for_DMS_Trt.RDS")

unite_g2_in_14_no_sxChr_trt = all_samples_normalised_filtered_unite_g2_in_14_trt[!all_samples_normalised_filtered_unite_g2_in_14_trt$chr %in% c("Gy_chrXIX","Gy_chrUn"),]

myDiff_4=calculateDiffMeth(unite_g2_in_14_no_sxChr_trt, covariates = covariates, mc.cores = 3)

saveRDS(myDiff_4, file = "DMS_Trt.RDS")

myDiff25p_4=getMethylDiff(myDiff_3,difference=25,qvalue=0.01)


save(myDiff_2,myDiff_3,myDiff_4, file = "DMS.RData")

unite_for_DMS_trt<-readRDS("unite_for_DMS_Trt.RDS")

unite_for_DMS_trt_nosxchr = unite_for_DMS_trt[!unite_for_DMS_trt$chr %in% c("Gy_chrXIX","Gy_chrUn"),]


tile_trt<-tileMethylCounts(unite_for_DMS_trt_nosxchr,win.size=100,step.size=100,cov.bases = 10)
tile_trt_diff<-calculateDiffMeth(tile_trt,covariates = covariates,mc.cores = 3)

tile_trt_15=getMethylDiff(tile_trt_diff,difference=15,qvalue=0.01)




