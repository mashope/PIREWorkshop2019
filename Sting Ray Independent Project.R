#useful commands for exploring data
library(tidyverse) #load tidyverse so ggplot will actually work
str(bkr_sjr_old)     #show columns in data
summary(bkr_sjr_old) #show summary stats
table(bkr_sjr_old$SITE) #number of observations by site
by(bkr_sjr_old$CYP, bkr_sjr_old$SITE, mean) #mean CYP by Site
by(bkr_sjr_old$GST, bkr_sjr_old$SITE, mean) #mean GST by Site
by(bkr_sjr_old$UGT, bkr_sjr_old$SITE, mean) #mean UGT by Site
by(bkr_sjr_old$LPO, bkr_sjr_old$SITE, mean, na.rm = TRUE) #mean LPO by Site (with LPO NA removed)
pairs(bkr_sjr_old) #matrix scatterplot to visualize correlation

cor(bkr_sjr_old$CYP, bkr_sjr_old$GST) #correlation between CYP and GST
cor(bkr_sjr_old$CYP, bkr_sjr_old$UGT) #correlation between CYP and UGT
cor(bkr_sjr_old$CYP, bkr_sjr_old$LPO, use = "complete.obs") #correlation between CYP and LPO (with LPO NA removed)
cor(bkr_sjr_old$GST, bkr_sjr_old$UGT) #correlation between GST and UGT
cor(bkr_sjr_old$GST, bkr_sjr_old$LPO, use = "complete.obs") #correlation between GST and LPO (with LPO NA removed)
cor(bkr_sjr_old$UGT, bkr_sjr_old$LPO, use = "complete.obs") #correlation between UGT and LPO (with LPO NA removed)


#UNUSED

#table(bkr_sjr_old$CYP, bkr_sjr_old$SITE) #number of observations
#range(bkr_sjr_old$GST) #min and max
#rank(iris$Sepal.Length) #return rank of each ovservation 


#Make scatterplot of CYP to GST    
CYP2GST <- c("CYP", "GST")
CYP2GST <- bkr_sjr_old[CYP2GST] 

GC <- ggplot(data = CYP2GST)  + aes(x = CYP, y = GST) + geom_point()
geom_smooth(method = "lm")
show(GC)

GC + geom_smooth()
GC_TITLE <- GC + geom_smooth(method = "lm")
GC_TITLE + ggtitle("Stingray Biomarkers")

#Make scatterplot of UGT to LPO    
UGT2LPO <- c("UGT", "LPO")
UGT2LPO <- bkr_sjr_old[UGT2LPO] 
  
UL <- ggplot(data = UGT2LPO)  + aes(x = UGT, y = LPO) + geom_point()
geom_smooth(method = "lm")
show(UL)

UL + geom_smooth()
UL_TITLE <- UL + geom_smooth(method = "lm")
UL_TITLE + ggtitle("Stingray Biomarkers")
#Set xlim to cut extra space
UL_TITLE + xlim(0, 60) +ggtitle("Stingray Biomarkers")
#STATS FOR DAYS

###################################################################################################################################################
#Linear Model UGT to LPO
linear_model_UGT2LPO <- lm(UGT ~ LPO, data = bkr_sjr_old)
#Get p-values
summary(linear_model_UGT2LPO, xlim(80))  

#Linear Model CYP to GST
linear_model_CYP2GST <- lm(CYP ~ GST, data = bkr_sjr_old)
#get Df, SumSq, MeanSq, FValue, and Pr(>F)
anova(linear_model_CYP2GST)
#Get Residuals, R-squared (adjusted and multiple), F-statistic, and p-values
summary(linear_model_CYP2GST)  

#Linear Model SITE to GST
linear_model_SITE2GST <- lm(GST ~ SITE, data = bkr_sjr_old)
#get Df, SumSq, MeanSq, FValue, and Pr(>F)
anova(linear_model_SITE2GST)


ray_stats <- bkr_sjr_old %>%
  group_by(SITE) %>%
  summarise(
    meanGST = mean(GST),
    meanUGT = mean(UGT),
    meanLPO = mean(LPO, na.rm = TRUE), 
    meanCYP = mean(CYP))
ray_stats
