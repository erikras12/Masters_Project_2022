###models####################################################################
##########################################################################
load("/data/SBCS-EizaguirreLab/Eri/data/coverage/metadata.RData")
load("/data/SBCS-EizaguirreLab/Eri/data/coverage/for_DMS.RData")

fullMetadata4Eri$residuals_BCI_sex<-residuals(lmer(data = fullMetadata4Eri,BCI ~ Sex + (1|brotherPairID)))

BCIII<-lmer(BCI ~ Sex * trtG1G2 + (1|brotherPairID),data = myMetaData)
anova(BCIII)
step(BCIII)


plop<-lmer(data = fullMetadata4Eri,res__2~trtG1G2 * residuals_BCI_sex + (1|brotherPairID))
plop0<-lmer(data = fullMetadata4Eri,res_2~ 1 + (1|brotherPairID))

######################################################################################################
plop<-lmer(data = fullMetadata4Eri,res__2~ Sex * trtG1G2 + (1|brotherPairID))
summary(plop)
anova(plop)
step(plop)
plop2<-lmer(data = myMetaData,res__2 ~ trtG1G2 + (1 | brotherPairID))

emmeans(plop2, list(pairwise ~ trtG1G2), adjust = "tukey")

res_B_S_trt <- lmer(BCI ~ Sex * residuals_BCI_sex + (1|brotherPairID),data = myMetaData)
anova(res_B_S_trt)
step(res_B_S_trt)

myMetaData$res_trt_sex_BCI<-residuals(lm(residuals_BCI_sex~PAT * outcome ,data = fullMetadata4Eri))

model_for_chris<-lmer(residuals_BCI_sex~PAT * outcome + (1|brotherPairID),data = myMetaData)

#residuals are taken in this way when the independent factors are actually dependent
#to observe whether they are independent you could do a chi squared test to see if they correlate or a simple lm with anova analyis 
#to see if they significnalt alter eachother

myMetaData$res_trt_sex_BCI <- residuals(lm(residuals_BCI_sex ~ PAT + outcome + PAT:outcome,data = myMetaData))

#################################################################################################################
###Whose infection status has more of an impact on offspring methylation, parents or offspring?
#################################################################################################################

trt_sex_BCI_mod<-lmer(res_No_meth_CpG_over_Number_of_covered_CpGs ~ Sex * PAT* outcome +(1|brotherPairID),data = myMetaData)
anova(trt_sex_BCI_mod)
plopppp<-ggpredict(trt_sex_BCI_mod)
summary(trt_sex_BCI_mod)
step(trt_sex_BCI_mod)
#################################################################################################################
###### FINAL MODEL ##############################################################################################
#################################################################################################################

trt_sex_BCI_mod<-lmer(res_No_meth_CpG_over_Number_of_covered_CpGs ~ Sex + PAT + (1|brotherPairID),data = myMetaData)
anova(trt_sex_BCI_mod)
summary(trt_sex_BCI_mod)
final_model_pred<-ggpredict(trt_sex_BCI_mod,terms = c("Sex"))
p<-plot(final_model_pred)
p
p +
  # Scales and theme
  scale_color_manual(values = c("deeppink", "blue")) +
  theme_bw() +
  labs(x = "Sex", y = "Predicted methylation level") +
  theme(axis.text=element_text(size=12),
  axis.title=element_text(size=14,face="bold"))

