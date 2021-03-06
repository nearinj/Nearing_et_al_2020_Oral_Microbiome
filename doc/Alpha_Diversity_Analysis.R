## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(HealthyOralMicrobiome)

## ---- Load_data----------------------------------------------------------

Metadata <- HealthyOralMicrobiome::Healthy_Metadata
dim(Metadata)


Metadata_alpha_Comp <- HealthyOralMicrobiome::get_alpha_divs(Metadata, "complete")


## ---- faiths_analysis----------------------------------------------------

features_to_test <- colnames(Metadata)[-42]
features_to_test

### test each feature
faiths_p <- vector()

for(feat in features_to_test){
  
  faiths_p[[feat]] <- HealthyOralMicrobiome::run_alpha_compare(df= Metadata_alpha_Comp, feature = feat, type="faiths")  
  
}

faiths_p
#apply fdr
faiths_q <- p.adjust(faiths_p, method = "fdr")
min(faiths_q)

### none pass alpha value of 0.1



## ---- shannon_analysis---------------------------------------------------

shannon_p <- vector()

for(feat in features_to_test){
  
  shannon_p[[feat]] <- HealthyOralMicrobiome::run_alpha_compare(df=Metadata_alpha_Comp, feature = feat, type="shannon")
}


shannon_q <- p.adjust(shannon_p, method = "fdr")
min(shannon_q)

## none pass fdr of 0.1


## ---- evenness_analysis--------------------------------------------------

evenness_p <- vector()

for(feat in features_to_test){
  
  evenness_p[[feat]] <- HealthyOralMicrobiome::run_alpha_compare(df=Metadata_alpha_Comp, feature=feat, type="evenness")
  
}

evenness_q <- p.adjust(evenness_p, method = "fdr")
min(evenness_q)

## none pass fdr of 0.1



## ---- richness_analysis--------------------------------------------------

richness_p <- vector()

for(feat in features_to_test){
  
  richness_p[[feat]] <- HealthyOralMicrobiome::run_alpha_compare(df=Metadata_alpha_Comp, feature=feat, type="richness")
  
}

richness_q <- p.adjust(richness_p, method="fdr")
min(richness_q)

## none pass alpha of 0.1



