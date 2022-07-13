###correcting for coverage bias########################
load("filtered_and_normalised_obj.RData")

G2_cov = methylKit::reorganize(normalised_filtered.myObj_cl,
                               sample.ids=metData$ID[!metData$ID %in% IDs &
                                                       metData$Generat %in% "O"],
                               treatment=metData$treatment_NUM[!metData$ID %in% IDs &
                                                                 metData$Generat %in% "O"])
G2_cov_unite= methylKit::unite(G2_cov,mc.cores = 8)

G2_cov_no_sxChr = G2_cov_unite[!G2_cov_unite$chr %in% c("Gy_chrXIX","Gy_chrUn"),]

unite_G2_in_14 = methylKit::unite(G2_cov,
                                  min.per.group=14L, mc.cores=3)# try with 8 cores

unite_G2_in_14_nosxchr = unite_G2_in_14[!unite_G2_in_14$chr %in% c("Gy_chrXIX","Gy_chrUn"),]

perc_meth_G2<-methylKit::percMethylation(G2_cov_no_sxChr)
meth_Cs = colSums(perc_meth_G2>=70 & !is.na(perc_meth_G2))
sum(meth_Cs)
Covered_Cs = colSums(!is.na(perc_meth_G2))
sum(Covered_Cs)
res__2<-residuals(lm(meth_Cs~Covered_Cs))


#G2_cov_no_offs_sxChr = G2_cov_no_sxChr[!G2_cov_no_sxChr$ID %in% IDs & G2_cov_no_sxChr$Generat %in% "O",]
percMethMat = methylKit::percMethylation(G2_cov_no_sxChr)
percMethDF = data.frame(SampleID = colnames(percMethMat),
                        ## number of methylated sites
                        Number_of_methylated_CpGs = colSums(percMethMat>=70 & !is.na(percMethMat)),#512493
                        ## number of sites covered in this sample
                        Number_of_covered_CpGs = colSums(!is.na(percMethMat)), #1015735
                        ## number of sites NOT covered in this sample
                        Nbr_NOTcoveredCpG = colSums(is.na(percMethMat)))

percMethDF$RMS_coveredCpG = percMethDF$Number_of_methylated_CpGs / percMethDF$Number_of_covered_CpGs

myMetaData = merge(fullMetadata4Eri, percMethDF)

myMetaData$res_No_meth_CpG_over_Number_of_covered_CpGs= residuals(
  lm(myMetaData$Number_of_methylated_CpGs ~ myMetaData$Number_of_covered_CpGs))

save(myMetaData,unite_G2_in_14_nosxchr,file = "for_DMS.RData")

cor.test(myMetaData$Number_of_methylated_CpGs, myMetaData$Number_of_covered_CpGs, method = "spearman")

plot(myMetaData$Number_of_methylated_CpGs, myMetaData$Number_of_covered_CpGs)
plot(myMetaData$RMS_coveredCpG, myMetaData$Number_of_covered_CpGs)
plot(myMetaData$res_No_meth_CpG_over_Number_of_covered_CpGs, myMetaData$Number_of_covered_CpGs)
cor.test(myMetaData$Number_of_methylated_CpGs, myMetaData$Number_of_covered_CpGs, method = "spearman")



cor.test(myMetaData$res_No_meth_CpG_over_Number_of_covered_CpGs,myMetaData$Number_of_covered_CpGs, method = "spearman")


#################################################################################################################
###Is there a difference in mapping efficeny between Males and Females that could introduce bias in our results
#################################################################################################################
#make two dataframes containing male and female samples
female_perc_align<-fullMetadata4Eri[fullMetadata4Eri$Sex == "F",]
male_perc_align<-fullMetadata4Eri[fullMetadata4Eri$Sex == "M",]

take the mapping effeciceny from both male and females
male_map_eff<-male_perc_align$mappingPerc_all_chr
female_map_eff<-female_perc_align$mappingPerc_all_chr

#see if there is a sig difference between the % mapping
wilcox.test(female_map_eff, male_map_eff,alternative = "two.sided")
#there is (this could introduce bias so we have to elimate this bias)
boxplot(male_map_eff, female_map_eff ,names=c("Male_map_eff","Female_map_eff"))

#we are mapping to a female genome and this leads to increased mapping especially at chromosome 19 (Female Sex chr)
#so we remove these chromosomesand look at the mapping %
male_map_eff_auto<-male_perc_align$mappingPerc_autosomes
female_map_eff_auto<-female_perc_align$mappingPerc_autosomes

wilcox.test(female_map_eff_auto, male_map_eff_auto,alternative = "two.sided")
#there is no sig difference now
boxplot(male_map_eff_auto, female_map_eff_auto ,names=c("Male_map_eff","Female_map_eff"))
