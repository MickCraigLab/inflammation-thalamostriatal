library(tidyverse)
library(lme4)
library(lmerTest)
library(ggpubr) 
library(ggplot2)
library(dplyr)
library(emmeans)
library(Hmisc)
library(corrplot)
library(scales)
library(factoextra)
library(ggfortify)
library(car)
library(glmnet) # for elastic net model
library(lm.beta) # needed for looking at loadings in the linear models of depression & fatigue against individual cytokines

#Set working directory - use comments to select either Windows or Mac environments
#setwd("C:/Users/mickc/OneDrive - University of Glasgow/Glasgow_OneDrive/Projects/Centaur")
setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Glasgow_OneDrive/Projects/Centaur")

#Import clinical psychiatry phenotype
psych_data<-read.csv("Centaur_clinical_phenotype_cleaned.csv",header=TRUE)

#Import biomarker data
#Note biomarker data in wide format
biomarker_data<-read.csv("Centaur_rawECL_dataset_cleaned.csv", header=TRUE, na.strings="#DIV/0!")

#Filter out entries with no biomarker data
psych_data_cleaned<-filter(psych_data, CK_data=="Y")
biomarker_data_cleaned<-filter(biomarker_data, CK_data=="Y")


#Replace missing values with NA - ignore, changed code to filter on import
#biomarker_data_cleaned[biomarker_data_cleaned=="#DIV/0!"] <-NA

#Pull out depression & fatigure scores for correlation matrix
psych_subset<-subset(psych_data_cleaned, select = c(Study_ID, Depression_T.Score,PROMIS_fatigue_total))

#dropping missing values
biomarker_data_cleaned2<- biomarker_data_cleaned%>% select_if(~ !any(is.na(.)))


#Append depression & fatigue scores to new data frame for correlation matrix; create second copy for stats
data_for_corr<-left_join(biomarker_data_cleaned2,psych_subset, by=c("Study_ID"))
psych_biomarker<-left_join(biomarker_data_cleaned2,psych_subset, by=c("Study_ID"))

#Drop categorical variables for participant ID, cytokine & fMRI
data_for_corr<-subset(data_for_corr, select = -c(Study_ID,CK_data, fMRI_data))

#Pull out cytokes to scale
cytokines<-subset(data_for_corr, select = -c(Depression_T.Score, PROMIS_fatigue_total))
cytokines_scaled<-scale(cytokines)

#merge back togther
data_for_corr_scaled<-cbind(cytokines_scaled,data_for_corr[,c("Depression_T.Score", "PROMIS_fatigue_total")])
psych_biomarker_scaled<-cbind(psych_biomarker[,c("Study_ID")],cytokines_scaled,psych_biomarker[,c("Depression_T.Score", "PROMIS_fatigue_total")])

#Run correlations, dropping columns with missing values

correlation_matrix_data<-cor(data_for_corr_scaled, use="complete.obs", method="spearman")

#Test_for_normality
normality_test<-do.call(rbind, lapply(data_for_corr_scaled, function(x) shapiro.test(x)[c("statistic", "p.value")]))


#Even log-normalising, data ain't normal so have to use non-parametric

# Function to generate p value matrix from rcorr test - got this from the internet in 2025, need to find source
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
# This came from https://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

correlation_matrix_data.p<-cor.mtest(correlation_matrix_data)

#Adjust p values to account for multiple comparisons
correlation_matrix_data.padj<-matrix(p.adjust(correlation_matrix_data.p, method = "fdr"), nrow = nrow(correlation_matrix_data.p))
colnames(correlation_matrix_data.padj) <- colnames(correlation_matrix_data.p)
rownames(correlation_matrix_data.padj) <- rownames(correlation_matrix_data.p)

#Draw correlation plot
Correlation_plot<-corrplot(correlation_matrix_data, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.5,
         p.mat=correlation_matrix_data.padj, sig.level = 0.05, insig = "blank")



#Individual correlations
Dep_BAFF.plot<-psych_biomarker%>%
  ggplot(aes(x=Depression_T.Score, y=BAFF))+
  geom_point(colour="black", size = 1)+
  theme_classic2()


#Determine whether we can model depression & fatigue scores with markers that come up as significant in the exploratory analysis
#Depression positively correlates with BAFF, EPO, GROalpha, and leptin, and negatively correlates with Eotaxin and PYY.
#Fatigue positively correlates with EPO, IL1RA, IL2RA, Mip1alpha and leptin, and negatively correlates with Eotaxin and IL17D.

Depression_lm<-lm(Depression_T.Score ~ BAFF  + EPO + CXCL1 + Leptin + Eotaxin + PYY ,data=psych_biomarker_scaled)
summary(Depression_lm)

Fatigue_lm<-lm(PROMIS_fatigue_total ~ EPO + IL.1Ra + IL.2Ra + CCL3 + Leptin + Eotaxin + IL.17D, data=psych_biomarker_scaled)
summary(Fatigue_lm)

#Check variance inflation factor
vif(Depression_lm)
vif(Fatigue_lm)

#Check residuals
par(mfrow = c(2,2))
#plot(Depression_lm)

#Check residuals
par(mfrow = c(2,2))
#plot(Fatigue_lm)


lm.beta(Depression_lm)
lm.beta(Fatigue_lm)

#Elastic net penalised regression


# Predictor matrix
x <- as.matrix(cytokines_scaled)

# Outcomes for both factors
y1 <- psych_biomarker_scaled$Depression_T.Score
y2 <- psych_biomarker_scaled$PROMIS_fatigue_total


#Run for depression
set.seed(123)

cv_enet_depression <- cv.glmnet(x, y1,
                     alpha = 0.5,   # Elastic Net (0 = ridge, 1 = lasso)
                     standardize = TRUE)

# Plot CV curve
#plot(cv_enet_depression)
     

lambda_min <- cv_enet_depression$lambda.min
lambda_1se <- cv_enet_depression$lambda.1se

#Elastic net, non-zero coefficients are significant
coef(cv_enet_depression, s = "lambda.min")

#Lasso method
cv_lasso_dep <- cv.glmnet(x, y1, alpha = 1)
coef(cv_lasso_dep, s = "lambda.min")

#Ridge method
cv_ridge_dep <- cv.glmnet(x, y1, alpha = 0)
coef(cv_ridge_dep, s = "lambda.min")

#Run for fatigue
set.seed(123)

cv_enet_fatigue <- cv.glmnet(x, y2,
                                alpha = 0.5,   # Elastic Net (0 = ridge, 1 = lasso)
                                standardize = TRUE)

# Plot CV curve
#plot(cv_enet_fatigue)


lambda_min <- cv_enet_fatigue$lambda.min
lambda_1se <- cv_enet_fatigue$lambda.1se

#Elastic net, non-zero coefficients are significant
coef(cv_enet_fatigue, s = "lambda.min")

#Lasso method
cv_lasso_fatigue <- cv.glmnet(x, y2, alpha = 1)
coef(cv_lasso_fatigue, s = "lambda.min")

#Ridge method
cv_ridge_fatigue <- cv.glmnet(x, y2, alpha = 0)
coef(cv_ridge_fatigue, s = "lambda.min")

#Now try dimensional reduction
#PCA regression
pca <- prcomp(cytokines_scaled)
summary(pca)

PCA_lm_dep<-lm(psych_biomarker_scaled$Depression_T.Score ~ pca$x[,1] + pca$x[,2])
summary(PCA_lm_dep)

PCA_lm_fat<-lm(psych_biomarker_scaled$PROMIS_fatigue_total ~ pca$x[,1] + pca$x[,2])
summary(PCA_lm_fat)

#Using all of the cytokine data in an unbiased way gives us no significant influences on depression or fatigue - this tells us that no single cytokine
#is driving this phenotype - useful to know.

#PCA regression using just the cocktail of cytokines that come up from the correlation matrix for depression to extract the shared inflammatory signal
pca_depression <- prcomp(psych_biomarker_scaled[, c("BAFF","EPO","CXCL1","Leptin","Eotaxin","PYY")], scale = TRUE)

summary(pca_depression)

pca_lm_depression<-lm(psych_biomarker_scaled$Depression_T.Score ~ pca_depression$x[,1])
summary(pca_lm_depression)

#Check loading factors for PCA
pca_depression$rotation[,1]

#PCA regression for fatigue again to extract shared inflammatory signal
pca_fatigue <- prcomp(psych_biomarker_scaled[, c("EPO", "IL.1Ra", "IL.2Ra", "CCL3", "Leptin", "Eotaxin", "IL.17D")], scale = TRUE)
summary(pca_fatigue)

pca_lm_fatigue<-lm(psych_biomarker_scaled$PROMIS_fatigue_total ~ pca_fatigue$x[,1])
summary(pca_lm_fatigue)

pca_fatigue$rotation[,1]

#Next, let's combine both sets of symptoms and cytokines into a multivariate linear model
pca_combined<- prcomp(psych_biomarker_scaled[, c("BAFF", "EPO", "CXCL1", "IL.1Ra", "IL.2Ra", "CCL3", "Leptin", "Eotaxin", "IL.17D", "PYY")], scale = TRUE)

combined_lm <- lm(cbind(Depression_T.Score, PROMIS_fatigue_total) ~ pca_combined$x[,1], data = psych_biomarker_scaled)
summary(manova(combined_lm)) #Is overall model significant? Yes

summary.aov(combined_lm) #Combined model perhaps not as good a fit as individuals but let's explore anyway

#Try partial least squares with both sets of cytokines together - commented  out as pls breaks the correlation plot package 
#library(pls)

#combined_pls_model <- plsr(cbind(Depression_T.Score, PROMIS_fatigue_total) ~  BAFF + EPO + CXCL1 + Leptin + Eotaxin + PYY + IL.1Ra + IL.2Ra + 
                  # CCL3 + IL.17D, data = psych_biomarker_scaled, scale = TRUE, validation = "CV")

#summary(combined_pls_model) #I *think* this tells us that the model is not robust enough to have predictive power

depression_loading<-pca_depression$rotation[,1]
fatigue_loading<-pca_fatigue$rotation[,1]

#Data wrangling to visualise loading for figure
#depression loading needs multiplied by -1 to match directionality of fatigue loading, so that markers having positive effect on depression score are 
#shown as positve on the graph
depression_adjust<-depression_loading*-1

Depression_loading_frame<-tibble(Psych_feature="Depression", Cytokine=names(depression_adjust), Loading=depression_adjust)
Fatigue_loading_frame<-tibble(Psych_feature="Fatigue", Cytokine=names(fatigue_loading), Loading=fatigue_loading)

Merged_loading<-bind_rows(Depression_loading_frame, Fatigue_loading_frame)%>%
  mutate(bar_colour=ifelse(Loading>0, "pos","neg"))%>%
  arrange(Loading)%>%
  mutate(Cytokine=fct_inorder(Cytokine))

LoadingColour<- c(
  "pos"="#2B73EB",
  "neg"="#C6403D"
)

#Plot PCA loadings
Loading_plot<-ggplot(data=Merged_loading, aes(x=Loading, y=Cytokine, fill=bar_colour))+
  geom_col()+
  scale_x_continuous(limits=c(-0.75,0.75))+
  scale_fill_manual(values=LoadingColour)+
  facet_wrap(~Psych_feature)+
  theme_bw()+
  labs(x="PCA Loading", y="Blood marker")+
  rremove("legend")
 
ggsave(file="Loading_plot.pdf", plot=Loading_plot, height=4, width=8)

#Plot PCA loadings
#Pull out PCA scores for first component for each participant
Dep_PCA_scores<-as.data.frame(pca_depression$x)
Fat_PCA_scores<-as.data.frame(pca_fatigue$x)

#Merge scores into dataframes
All_depression_data<-cbind(psych_biomarker_scaled,Dep_PCA_scores)
All_fatigue_data<-cbind(psych_biomarker_scaled,Fat_PCA_scores)

#Plot correlation between 
Dep_PCA.plot<-All_depression_data%>%
  ggplot(aes(x=Depression_T.Score, y=(PC1*-1)))+ #multiple by negative 1 to match figure 1A - directionality of PCA is artibtrary to this matches biology better
  geom_point(colour="black", size = 1)+
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE)+
  theme_classic2()+
  labs(x="PROMIS Depression Score", y="Principal component 1")

Fat_PCA.plot<-All_fatigue_data%>%
  ggplot(aes(x=PROMIS_fatigue_total, y=PC1))+
  geom_point(colour="black", size = 1)+
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE)+
  theme_classic2()+
  labs(x="PROMIS Fatigue Score", y="Principal component 1")


