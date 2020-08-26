## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(HealthyOralMicrobiome)
library(vegan)

## ----load_data-----------------------------------------------------------
#Load in metadata
Metadata <- HealthyOralMicrobiome::Healthy_Metadata
dim(Metadata)

W_unifrac <- HealthyOralMicrobiome::get_beta_div("w_unifrac", "validation", Metadata)
Bray_curt <- HealthyOralMicrobiome::get_beta_div("bray_curtis", "validation", Metadata)
dim(W_unifrac[["Distance"]])

## ---- eval=F-------------------------------------------------------------
#  weighted_sig_feats <- c("PM_WAIST_AVG", "PM_WAIST_HIP_RATIO", "PM_STANDING_HEIGHT_AVG", "PM_BIOIMPED_WEIGHT",
#                          "A_SDC_AGE_CALC", "A_SDC_GENDER", "A_SLE_LIGHT_EXP", "NUT_VEG_DAY_QTY",
#                          "REFINED_GRAIN_SERVINGS_DAY_QTY", "NUTS_SEEDS_SERVINGS_PER_DAY", "SALT_SEASONING", "PM_BIOIMPED_FFM")
#  
#  W_Valid_res <- list()
#  sample_size <- vector()
#  
#  for(feat in weighted_sig_feats){
#      feat_df <- W_unifrac[["Metadata"]][complete.cases(W_unifrac[["Metadata"]][,feat]),]
#      feat_dist <- W_unifrac[["Distance"]][rownames(feat_df), rownames(feat_df)]
#      sample_size[[feat]] <- dim(feat_dist)[1]
#  
#      W_Valid_res[[feat]] <-  adonis2(feat_dist ~ feat_df[,"Extraction_Number"] + feat_df[,feat], permutations = 1000, parallel = 20, by="margin")
#  }
#  sample_size
#  
#  bray_sig_feats <- c("A_HS_DENTAL_VISIT_LAST", "PM_BIOIMPED_BMI", "PM_WAIST_AVG", "PM_WAIST_HIP_RATIO",
#                      "PM_STANDING_HEIGHT_AVG", "PM_BIOIMPED_WEIGHT", "A_SDC_AGE_CALC", "A_SDC_GENDER",
#                      "A_SLE_LIGHT_EXP", "NUT_VEG_DAY_QTY", "NUT_JUICE_DAY_QTY", "REFINED_GRAIN_SERVINGS_DAY_QTY",
#                      "PM_BIOIMPED_FFM")
#  
#  Bray_Valid_res <- list()
#  Bray_sample_size <- vector()
#  
#  for(feat in bray_sig_feats){
#  
#    feat_df <- Bray_curt[["Metadata"]][complete.cases(Bray_curt[["Metadata"]][,feat]),]
#    feat_dist <- Bray_curt[["Distance"]][rownames(feat_df), rownames(feat_df)]
#    Bray_sample_size[[feat]] <- dim(feat_dist)[1]
#  
#    Bray_Valid_res[[feat]] <- adonis2(feat_dist ~ feat_df[,"Extraction_Number"] + feat_df[,feat], permutations = 1000, parallel = 20, by="margin")
#  
#  }
#  
#  

## ------------------------------------------------------------------------
Valid_W_unifrac_res <- HealthyOralMicrobiome::Valid_W_unifrac_res

Valid_w_uni_p <- vector()

weighted_sig_feats <- c("PM_WAIST_AVG", "PM_WAIST_HIP_RATIO", "PM_STANDING_HEIGHT_AVG", "PM_BIOIMPED_WEIGHT",
                        "A_SDC_AGE_CALC", "A_SDC_GENDER", "A_SLE_LIGHT_EXP", "NUT_VEG_DAY_QTY",
                        "REFINED_GRAIN_SERVINGS_DAY_QTY", "NUTS_SEEDS_SERVINGS_PER_DAY", "SALT_SEASONING", "PM_BIOIMPED_FFM")


for(feat in weighted_sig_feats){
  
  Valid_w_uni_p[[feat]] <- Valid_W_unifrac_res[[feat]]$`Pr(>F)`[2]
  
}


Valid_w_uni_p
Valid_sig_w_uni_p <- Valid_w_uni_p[which(Valid_w_uni_p < 0.05)]
### a few were not recovered...

valid_sig_w_uni_feats <- names(Valid_w_uni_p[which(Valid_w_uni_p < 0.05)])
valid_sig_w_uni_feats

Valid_w_uni_r2 <- vector()

for(feat in weighted_sig_feats){
  
  Valid_w_uni_r2[[feat]] <- Valid_W_unifrac_res[[feat]]$R2[2]
  
}

Valid_sig_w_uni_r2 <- Valid_w_uni_r2[which(Valid_w_uni_p < 0.05)]
Valid_sig_w_uni_r2

Valid_table_res_df <- data.frame("feature"=names(Valid_w_uni_p),
                                 "p"=Valid_w_uni_p,
                                 "R2"=Valid_w_uni_r2,
                                 stringsAsFactors = F)

Valid_table_sig_df <- data.frame("feature"=valid_sig_w_uni_feats,
                                 "p"=Valid_sig_w_uni_p,
                                 "R2"=Valid_sig_w_uni_r2,
                                 stringsAsFactors = F)
### make barplot for valid sig feats
Valid_table_sig_df$feature <- c("Waist Hip Ratio", "Height", "Weight", "Age", "Sex", "Fat Free Mass")
Valid_table_sig_df$Type <- c("Anthropometric", "Anthropometric", "Anthropometric", "Age", "Sex", "Anthropometric")

library(ggplot2)
library(cowplot)
theme_set(theme_classic())

Valid_table_sig_df$color <- "black"


for(i in 1:length(Valid_table_sig_df$Type)){
  
  if(Valid_table_sig_df$Type[i]=="Age"){
    Valid_table_sig_df$color[i] <- "Red"
  }else if(Valid_table_sig_df$Type[i]=="Anthropometric"){
    Valid_table_sig_df$color[i] <- "Purple"
  }else if(Valid_table_sig_df$Type[i]=="Diet"){
    Valid_table_sig_df$color[i] <- "Blue"
  }else if(Valid_table_sig_df$Type[i]=="Life Style"){
    Valid_table_sig_df$color[i] <- "Orange"
  }else if(Valid_table_sig_df$Type[i]=="Sex"){
    Valid_table_sig_df$color[i] <- "Pink"
  }else if(Valid_table_sig_df$Type[i]=="Life Style"){
    Valid_table_sig_df$color[i] <- "Green"
  }
  
}


### sort features by R2 value
Valid_table_sig_df$feature <- factor(Valid_table_sig_df$feature, levels = Valid_table_sig_df$feature[order(Valid_table_sig_df$R2, decreasing = T)])


effect_size_plot <- ggplot(Valid_table_sig_df, aes(x=feature, y=R2, fill=color)) + geom_bar(stat="identity") + coord_flip() +
  facet_grid(row=as.factor(Valid_table_sig_df$Type), scales="free", space="free") + ggtitle("Validation Weighted UniFrac Effect Sizes") + 
  theme(plot.title = element_text(hjust=0.5)) + theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  guides(fill=guide_legend(title="Type")) + scale_fill_identity(guide="legend", labels=c("Sex", "Anthropometric", "Age")) + xlab("Feature")
 
effect_size_plot


### Biplot

## We will need to remove some samples from the distance matrix that were NA...
Valid_sig_metadata <- W_unifrac[["Metadata"]][,valid_sig_w_uni_feats]
Valid_sig_metadata_comp <- Valid_sig_metadata[complete.cases(Valid_sig_metadata),]
dim(Valid_sig_metadata_comp)

Valid_sig_w_dist <- W_unifrac[["Distance"]][rownames(Valid_sig_metadata_comp), rownames(Valid_sig_metadata_comp)]
dim(Valid_sig_w_dist)

Valid_weighted_uni_pcoa <- cmdscale(Valid_sig_w_dist, k=2, eig=T)

#convert gender to numberic for complete_df
df_env <- Valid_sig_metadata_comp
table(df_env$A_SDC_GENDER)
df_env$A_SDC_GENDER <- as.numeric(as.character(df_env$A_SDC_GENDER))
df_env$A_SDC_GENDER
env_fit <- envfit(Valid_weighted_uni_pcoa ~ ., data=df_env[,valid_sig_w_uni_feats])


comp1 <- Valid_weighted_uni_pcoa$eig[1]/sum(Valid_weighted_uni_pcoa$eig)
comp1
comp2  <- Valid_weighted_uni_pcoa$eig[2]/sum(Valid_weighted_uni_pcoa$eig)
comp2

barplot(Valid_weighted_uni_pcoa$eig[1:20])
#save for composition barplot later...


theme_set(theme_cowplot())
Valid_plot_data_weight <- data.frame(PC1=Valid_weighted_uni_pcoa$points[,1], PC2=Valid_weighted_uni_pcoa$points[,2])

#get vectors

vec.sp.df <- as.data.frame(env_fit$vectors$arrows*sqrt(Valid_sig_w_uni_r2)*4)
rownames(vec.sp.df)
vec.sp.df$species <- c("Waist Hip Ratio", "Height", "Weight", "Age", "Sex", "Fat Free Mass")

library(ggrepel)
Valid_weighted_uni_gg_bi <- ggplot(data = Valid_plot_data_weight, aes(PC1, PC2)) + geom_point(colour="black", alpha=0.5) + coord_fixed() +
  geom_segment(data=vec.sp.df, aes(x=0, xend=Dim1, y=0, yend=Dim2), arrow = arrow(length = unit(0.5, "cm")), colour="red") +
  geom_text_repel(data=vec.sp.df, aes(x=Dim1, y=Dim2+.015, label=species), size=4) + xlab("PC1 46.34%") + ylab("PC2 19.01%") + ggtitle("Valid Weighted UniFrac")

Valid_weighted_uni_gg_bi

## ------------------------------------------------------------------------
Valid_Bray_res <- HealthyOralMicrobiome::Valid_Bray_res

Valid_Bray_p <- vector()

bray_sig_feats <- c("A_HS_DENTAL_VISIT_LAST", "PM_BIOIMPED_BMI", "PM_WAIST_AVG", "PM_WAIST_HIP_RATIO",
                    "PM_STANDING_HEIGHT_AVG", "PM_BIOIMPED_WEIGHT", "A_SDC_AGE_CALC", "A_SDC_GENDER", 
                    "A_SLE_LIGHT_EXP", "NUT_VEG_DAY_QTY", "NUT_JUICE_DAY_QTY", "REFINED_GRAIN_SERVINGS_DAY_QTY",
                    "PM_BIOIMPED_FFM")

for(feat in bray_sig_feats){
  
  Valid_Bray_p[[feat]] <- Valid_Bray_res[[feat]]$`Pr(>F)`[2]
  
}

Valid_Bray_sig_feats <- names(Valid_Bray_p[which(Valid_Bray_p < 0.05)])
Valid_Bray_sig_feats

Valid_Bray_sig_p <- Valid_Bray_p[which(Valid_Bray_p < 0.05)]

Valid_Bray_r2 <- vector()

for(feat in bray_sig_feats){
  Valid_Bray_r2[[feat]] <- Valid_Bray_res[[feat]]$R2[2]

}

Valid_Bray_sig_r2 <- Valid_Bray_r2[which(Valid_Bray_p < 0.05)]
Valid_Bray_sig_r2


Valid_Bray_res_Tab <- data.frame("feature"=names(Valid_Bray_p),
                                 "p"=Valid_Bray_p,
                                 "R2"=Valid_Bray_r2)


Valid_Bray_sig_res_Tab <- data.frame("feature"=Valid_Bray_sig_feats,
                                    "p"=Valid_Bray_sig_p,
                                    "R2"=Valid_Bray_sig_r2)

Valid_Bray_sig_res_Tab$feature <- c("Waist Size", "Waist Hip Ratio", "Height", "Weight", "Age", "Sex", "Fat Free Mass")
Valid_Bray_sig_res_Tab$color <- "black"
Valid_Bray_sig_res_Tab$Type <- c("Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "Age", "Sex", "Anthropometric")



library(ggplot2)
library(cowplot)
theme_set(theme_classic())

for(i in 1:length(Valid_Bray_sig_res_Tab$Type)){
  
  if(Valid_Bray_sig_res_Tab$Type[i]=="Age"){
    Valid_Bray_sig_res_Tab$color[i] <- "Red"
  }else if(Valid_Bray_sig_res_Tab$Type[i]=="Anthropometric"){
    Valid_Bray_sig_res_Tab$color[i] <- "Purple"
  }else if(Valid_Bray_sig_res_Tab$Type[i]=="Diet"){
    Valid_Bray_sig_res_Tab$color[i] <- "Blue"
  }else if(Valid_Bray_sig_res_Tab$Type[i]=="Life Style"){
    Valid_Bray_sig_res_Tab$color[i] <- "Orange"
  }else if(Valid_Bray_sig_res_Tab$Type[i]=="Sex"){
    Valid_Bray_sig_res_Tab$color[i] <- "Pink"
  }else if(Valid_Bray_sig_res_Tab$Type[i]=="Life Style"){
    Valid_Bray_sig_res_Tab$color[i] <- "Green"
  }
  
}
## Sort them by effect size
Valid_Bray_sig_res_Tab$feature <- factor(Valid_Bray_sig_res_Tab$feature, levels = Valid_Bray_sig_res_Tab$feature[order(Valid_Bray_sig_res_Tab$R2, decreasing = T)])


Valid_Bray_effect_size_plot <- ggplot(Valid_Bray_sig_res_Tab, aes(x=feature, y=R2, fill=color)) + geom_bar(stat="identity") + coord_flip() +
  facet_grid(row=as.factor(Valid_Bray_sig_res_Tab$Type), scales="free", space="free") + ggtitle("Validation Bray Curtis Effect Sizes") + 
  theme(plot.title = element_text(hjust=0.5)) + theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  guides(fill=guide_legend(title="Type")) + scale_fill_identity(guide="legend", labels=c("Sex", "Anthropometric", "Age"))
 
Valid_Bray_effect_size_plot

### Create biplot
## We will need to remove some samples from the distance matrix that were NA...
Valid_Bray_Metadata_Sig <- Bray_curt[["Metadata"]][,Valid_Bray_sig_feats]
Valid_Bray_Metadata_sig_comp <- Valid_Bray_Metadata_Sig[complete.cases(Valid_Bray_Metadata_Sig),]
dim(Valid_Bray_Metadata_sig_comp)
Valid_Bray_Dist_comp <- Bray_curt[["Distance"]][rownames(Valid_Bray_Metadata_sig_comp), rownames(Valid_Bray_Metadata_sig_comp)]

Valid_Bray_pcoa <- cmdscale(Valid_Bray_Dist_comp, k=2, eig=T)

#convert gender to numberic for complete_df
Bray_df_env <- Valid_Bray_Metadata_sig_comp
table(Bray_df_env$A_SDC_GENDER)
Bray_df_env$A_SDC_GENDER <- as.numeric(as.character(Bray_df_env$A_SDC_GENDER))
Bray_df_env$A_SDC_GENDER
Bray_env_fit <- envfit(Valid_Bray_pcoa ~ ., data=Bray_df_env[,Valid_Bray_sig_feats])


Bray_comp1 <- Valid_Bray_pcoa$eig[1]/sum(Valid_Bray_pcoa$eig)
Bray_comp1
Bray_comp2  <- Valid_Bray_pcoa$eig[2]/sum(Valid_Bray_pcoa$eig)
Bray_comp2

barplot(Valid_Bray_pcoa$eig[1:20])

theme_set(theme_cowplot())
Valid_plot_data_bray <- data.frame(PC1=Valid_Bray_pcoa$points[,1], PC2=Valid_Bray_pcoa$points[,2])

#get vectors

Bray_vec.sp.df <- as.data.frame(Bray_env_fit$vectors$arrows*sqrt(Valid_Bray_sig_r2)*4)
rownames(Bray_vec.sp.df)
Bray_vec.sp.df$species <- c("Waist Size", "Waist Hip Ratio", "Height", "Weight", "Age", "Sex", "Fat Free Mass")

library(ggrepel)
Valid_Bray_gg_bi <- ggplot(data = Valid_plot_data_bray, aes(PC1, PC2)) + geom_point(colour="black", alpha=0.5) + coord_fixed() +
  geom_segment(data=Bray_vec.sp.df, aes(x=0, xend=Dim1, y=0, yend=Dim2), arrow = arrow(length = unit(0.5, "cm")), colour="red") +
  geom_text_repel(data=Bray_vec.sp.df, aes(x=Dim1, y=Dim2+.015, label=species), size=4) + xlab("PC1 24.52%") + ylab("PC2 8.18%") + ggtitle("Valid Bray Curtis")

Valid_Bray_gg_bi


