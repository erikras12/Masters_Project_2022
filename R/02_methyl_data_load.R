###methylobjects#####
#################################################################################
#create methylRawList of all samples

load("metadata.RData")

fullMetadata4Eri$trtG1G2_NUM <- as.numeric(as.factor(fullMetadata4Eri$trtG1G2))
#read coverage files into a methylbase object
myObj = methylKit::methRead(as.list(temp),
                            mincov = 10,
                            pipeline = "bismarkCoverage",
                            sample.id = as.list(metData$ID),
                            assembly = "Gynogen_pchrom_assembly_all",
                            treatment = metData$trtG1G2_NUM,
                            context = "CpG")

#Get some methylation and coverage stats to see if cleaning up is required
methylKit::getMethylationStats(myObj[[2]],plot=TRUE,both.strands=FALSE)
methylKit::getCoverageStats(myObj[[2]],plot=TRUE,both.strands=FALSE)

#unite all samples into one data form
meth=methylKit::unite(myObj, destrand=FALSE)

#cluster analysis on united object to look for outliers

methylKit::clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

pc= methylKit::PCASamples(meth,obj.return = TRUE, adj.lim=c(1,1))


#removing family 12 due to no offspring
family12 = c("S74",  "S75",  "S126", "S127")

#sample 12 removed due to poor methylation results
IDs= c("S12",  "S22",  "S110", "S118", "S142", family12)

#samples we want to keep
samples_remain=metData$ID[!metData$ID %in% IDs]
#treatment of said samples
treatment_of_rem_samples = metData$treatment_NUM[!metData$ID %in% IDs]

#remove samples is IDs
myObj_cl=reorganize(myObj, sample.ids=samples_remain,
                    treatment=treatment_of_rem_samples)

myObj_cl=methylKit::reorganize(
  myObj,
  sample.ids=metData$ID[!metData$ID %in% IDs],
  treatment=metData$treatment_NUM[!metData$ID %in% IDs])

##############################################################################
#working with samples I want to use
##############################################################################

#filter myObj to remove bias from PCA
filtered.myObj_cl = methylKit::filterByCoverage(myObj_cl,lo.count=10,lo.perc=NULL,
                                                hi.count=NULL,hi.perc=99.9)
#normalise coverage between samples
normalised_filtered.myObj_cl= methylKit::normalizeCoverage(filtered.myObj_cl)
save(normalised_filtered.myObj_cl, file = "filtered_and_normalised_obj.RData")

###now unite the remaining samples
#all_samples = methylKit::unite(normalised_filtered.myObj_cl,mc.cores = 8)


#cluster analysis on the remaining samples to see the variation in coverage now
#that outliers are removed
pc1 = methylKit::PCASamples(all_samples,obj.return = TRUE, adj.lim=c(1,1))

methylKit::clusterSamples(all_samples, dist="correlation", method="ward", plot=TRUE)
