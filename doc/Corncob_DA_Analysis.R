## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(HealthyOralMicrobiome)
library(corncob)
library(phyloseq)
library(dplyr)
library(tibble)
library(gplots)

## ---- eval=F, Load_in_data-----------------------------------------------
#  #Load in metadata
#  Metadata <- HealthyOralMicrobiome::Healthy_Metadata
#  dim(Metadata)
#  
#  
#  #load in ASV table
#  ASV_tab <- HealthyOralMicrobiome::Healthy_ASV_table
#  dim(ASV_tab)
#  
#  #load in Genus table
#  Genus_tab <- HealthyOralMicrobiome::Healthy_Genus_Table
#  dim(Genus_tab)
#  
#  ### Filter the samples
#  Comp_ASV_data <- Filt_samples(depth=5000, setType="complete", Taxa_table = ASV_tab, Metadata = Metadata, parallel = 20, cutoff_prop = 0.1)
#  Comp_Genus_data <- Filt_samples(depth=5000, setType="complete", Taxa_table = Genus_tab, Metadata = Metadata, parallel = 20, cutoff_prop = 0.1)
#  
#  ### finally we will scale the data so that the values will be comparable to one another
#  
#  
#  scaled_metadata_g<- Comp_Genus_data[["Metadata"]] %>% rownames_to_column('sample_name') %>%
#    mutate_if(is.numeric, scale) %>%
#    column_to_rownames('sample_name')
#  
#  scaled_metadata_a <- Comp_ASV_data[["Metadata"]] %>% rownames_to_column('sample_name') %>%
#    mutate_if(is.numeric, scale) %>%
#    column_to_rownames('sample_name')
#  
#  features_to_test <- c("A_SDC_AGE_CALC",
#                        "PM_STANDING_HEIGHT_AVG",
#                        "PM_WAIST_AVG",
#                        "PM_BIOIMPED_WEIGHT",
#                        "PM_WAIST_HIP_RATIO",
#                        "PM_BIOIMPED_FFM",
#                        "SALT_SEASONING",
#                        "NUTS_SEEDS_SERVINGS_PER_DAY",
#                        "NUT_VEG_DAY_QTY",
#                        "REFINED_GRAIN_SERVINGS_DAY_QTY",
#                        "A_SLE_LIGHT_EXP",
#                        "A_SDC_GENDER",
#                        "PM_BIOIMPED_BMI",
#                        "NUT_JUICE_DAY_QTY",
#                        "A_HS_DENTAL_VISIT_LAST")
#  

## ----eval=FALSE, fig.fullwidth=T,Run_CornCob_Tests-----------------------
#  ### create phyloseq object to analysis
#  ASVs_comp <- phyloseq::otu_table(Comp_ASV_data[['Taxa_Tab']], taxa_are_rows = T)
#  Meta_comp <- phyloseq::sample_data(scaled_metadata_a)
#  phylo_comp <- phyloseq::phyloseq(ASVs_comp, Meta_comp)
#  
#  features_to_test <- c("A_SDC_AGE_CALC",
#                        "PM_STANDING_HEIGHT_AVG",
#                        "PM_WAIST_AVG",
#                        "PM_BIOIMPED_WEIGHT",
#                        "PM_WAIST_HIP_RATIO",
#                        "PM_BIOIMPED_FFM",
#                        "SALT_SEASONING",
#                        "NUTS_SEEDS_SERVINGS_PER_DAY",
#                        "NUT_VEG_DAY_QTY",
#                        "REFINED_GRAIN_SERVINGS_DAY_QTY",
#                        "A_SLE_LIGHT_EXP",
#                        "A_SDC_GENDER",
#                        "PM_BIOIMPED_BMI",
#                        "NUT_JUICE_DAY_QTY",
#                        "A_HS_DENTAL_VISIT_LAST")
#  #list of features that were significant in beta diversity analysis.
#  
#  
#  cl <- parallel::makeCluster(15)
#  doParallel::registerDoParallel(cl)
#  foreach::getDoParWorkers()
#  `%dopar%` <- foreach::`%dopar%`
#  Comp_ASV_CC_Res <- foreach::foreach(i=1:15) %dopar% HealthyOralMicrobiome::run_corn_cob_analysis(features_to_test[i], phylo = phylo_comp)
#  parallel::stopCluster(cl)
#  ### Genus Analysis
#  
#  Genus_comp <- phyloseq::otu_table(Comp_Genus_data[["Taxa_Tab"]], taxa_are_rows = T)
#  Meta_comp <- phyloseq::sample_data(scaled_metadata_g)
#  phylo_comp_g <- phyloseq::phyloseq(Genus_comp, Meta_comp)
#  
#  
#  
#  cl <- parallel::makeCluster(15)
#  doParallel::registerDoParallel(cl)
#  foreach::getDoParWorkers()
#  Genus_corn_cob_res <- foreach::foreach(i=1:15) %dopar% HealthyOralMicrobiome::run_corn_cob_analysis(features_to_test[i], phylo=phylo_comp_g)
#  parallel::stopCluster(cl)
#  
#  

## ---- eval=T, fig.fullwidth=TRUE, fig.width=10, fig.height=10, Analysis_of_results----
## Genus analysis first
Genus_CC_res <- HealthyOralMicrobiome::Comp_Genus_CC_Res

features_to_test <- c("A_SDC_AGE_CALC",
                      "PM_STANDING_HEIGHT_AVG",
                      "PM_WAIST_AVG",
                      "PM_BIOIMPED_WEIGHT",
                      "PM_WAIST_HIP_RATIO",
                      "PM_BIOIMPED_FFM",
                      "SALT_SEASONING",
                      "NUTS_SEEDS_SERVINGS_PER_DAY",
                      "NUT_VEG_DAY_QTY",
                      "REFINED_GRAIN_SERVINGS_DAY_QTY",
                      "A_SLE_LIGHT_EXP",
                      "A_SDC_GENDER",
                      "PM_BIOIMPED_BMI",
                      "NUT_JUICE_DAY_QTY",
                      "A_HS_DENTAL_VISIT_LAST")

## Extract Genus level p-values for each feature tested
G_CC_scaled_p <- list()
for(i in 1:15){
  G_CC_scaled_p[[i]] <- Genus_CC_res[[i]]$p_fdr
  
}

G_CC_scaled_p_df <- do.call(rbind, G_CC_scaled_p)
G_CC_scaled_p_df <- data.frame(t(G_CC_scaled_p_df))
colnames(G_CC_scaled_p_df) <- features_to_test

#remove any taxa that weren't significant in atleast on feature
filt_G_CC_scaled_p_df <- G_CC_scaled_p_df[which(apply(G_CC_scaled_p_df, 1, function(x) {any(x <= 0.1)})),]
dim(filt_G_CC_scaled_p_df)

### leaves 45 taxa

keep_genus <- rownames(filt_G_CC_scaled_p_df)

## get coefficents for each remaining model
G_Corn_cob_coef <- data.frame()
i <- 0
for(feat in features_to_test){
  i <- i +1
  if(feat=="A_SDC_GENDER"){
    feat <- "A_SDC_GENDER2"
  }
  message(feat)
  feat_name <- paste("mu.",feat,sep="")
  for(j in 1:length(Genus_CC_res[[i]]$all_models)){
      G_Corn_cob_coef[j,i] <- Genus_CC_res[[i]]$all_models[[j]]$coefficients[,1][which(rownames(Genus_CC_res[[i]]$all_models[[j]]$coefficients)==feat_name)]
    
  }
  
}

rownames(G_Corn_cob_coef) <- rownames(Genus_CC_res[[1]]$data@otu_table)
colnames(G_Corn_cob_coef) <- features_to_test


#rename taxa
filt_G_Corn_cob_coef <- G_Corn_cob_coef[rownames(filt_G_CC_scaled_p_df),]
renamed_filt_G_Corn_cob_coef <- filt_G_Corn_cob_coef
rownames(renamed_filt_G_Corn_cob_coef) <- gsub(".*D_4", "D_4", rownames(renamed_filt_G_Corn_cob_coef))
rownames(renamed_filt_G_Corn_cob_coef) <- gsub("D_4__", "", rownames(renamed_filt_G_Corn_cob_coef))
rownames(renamed_filt_G_Corn_cob_coef) <- gsub("D_5__", "", rownames(renamed_filt_G_Corn_cob_coef))
rownames(renamed_filt_G_Corn_cob_coef) <- gsub(";", " ", rownames(renamed_filt_G_Corn_cob_coef))

## manually fix some rows...
rownames(renamed_filt_G_Corn_cob_coef)
rownames(renamed_filt_G_Corn_cob_coef)[29] <- "Veillonellaceae unclassified"
rownames(renamed_filt_G_Corn_cob_coef)[34] <- "Neisseriaceae unclassified"
rownames(renamed_filt_G_Corn_cob_coef)[41] <- "Mollicutes uncultured"


#set effect non-significant effect sizes to 0
set_na_filt_G_Corn_cob_coef <- renamed_filt_G_Corn_cob_coef
set_na_filt_G_Corn_cob_coef[filt_G_CC_scaled_p_df > 0.1] <- 0

#set column names
colnames(set_na_filt_G_Corn_cob_coef)
colnames(set_na_filt_G_Corn_cob_coef) <- c("Age", "Height", "Waist Size", "Weight",
                                           "Waist Hip Ratio", "Fat Free Mass", "Salt Usage", "Nut/Seed Servings",
                                           "Vegetable Servings", "Refined Grain Servings", "Sleeping Light Exposure", "Sex",
                                           "Body Mass Index", "Juice Servings", "Last Dental Visit")

#fix Column label for Sex
table(Genus_CC_res[[1]]$data@sam_data$A_SDC_GENDER)
set_na_filt_G_Corn_cob_coef$Sex <- -set_na_filt_G_Corn_cob_coef$Sex
colnames(set_na_filt_G_Corn_cob_coef)[12] <- "Male*"

#get phylums of each taxa
Phylums <- gsub("D_0__.*;D_1__", "", rownames(filt_G_Corn_cob_coef))
Phylums <- gsub(";D_2__.*", "", Phylums)
Phylums 


#set colors for each phylum
taxa_colors <- vector()
for(taxa in 1:length(Phylums)){
  
  if(Phylums[[taxa]] == "Actinobacteria"){
      taxa_colors[[taxa]] <- "Red"
  }else if(Phylums[[taxa]]=="Bacteroidetes"){
    taxa_colors[[taxa]] <- "Blue"
  }else if(Phylums[[taxa]]=="Cyanobacteria"){
    taxa_colors[[taxa]] <- "Orange"
  }else if(Phylums[[taxa]]=="Firmicutes"){
    taxa_colors[[taxa]] <- "Brown"
  }else if(Phylums[[taxa]]=="Fusobacteria"){
    taxa_colors[[taxa]] <- "Purple"
  }else if(Phylums[[taxa]]=="Proteobacteria"){
    taxa_colors[[taxa]] <- "#664163"
  }else if(Phylums[[taxa]]=="Spirochaetes"){
    taxa_colors[[taxa]] <- "Black"
  }else if(Phylums[[taxa]]=="Synergistetes"){
    taxa_colors[[taxa]] <- "Cyan"
  }else if(Phylums[[taxa]]=="Tenericutes"){
    taxa_colors[[taxa]] <- "#4a421f"
  }else
    print("error")
}

#set up color scale for heatmap
breaks=seq(-2, 2, length.out = 180) #41 values
#now add outliers
breaks=append(breaks, 10000)
breaks=append(breaks, -10000, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="blue",mid="white",high="brown")

#set plot dims
layout_mat <- rbind(c(0,0,4), c(2,1,3), c(0,0,0))
lwid=c(.1,2.75,1.2)
lhei=c(.2,3,.3)

#plot

gen_heat <- heatmap.2(as.matrix(set_na_filt_G_Corn_cob_coef),
          na.color="black", density.info = "none", symkey = F, trace="none", srtCol=35, dendrogram = "none", key.xlab = "log odds", key.ylab = "p Value", key=T, col=mycol, lmat = layout_mat, lwid = lwid, lhei = lhei, cexRow = 1.3, cexCol = 1.3, distfun=function(x) dist(x, method="manhattan"), colRow = taxa_colors, Rowv = F,
          rowsep = c(2, 10,29,30, 38, 39, 40), sepcolor = "black", keysize=0.2, key.par = list(cex=0.01, cex.axis=60, cex.main=60, mar=c(50,0,50,10), mgp=c(.01,50,1), tcl=-10))





#### Now lets make the ASV heatmap...


ASV_CC_res <- HealthyOralMicrobiome::Comp_ASV_CC_Res

ACC_scaled_p <- list()
for(i in 1:15){
  ACC_scaled_p[[i]] <- ASV_CC_res[[i]]$p_fdr
  
}

ACC_scaled_p_df <- do.call(rbind, ACC_scaled_p)
ACC_scaled_p_df <- data.frame(t(ACC_scaled_p_df))
colnames(ACC_scaled_p_df) <- features_to_test

#remove ASVs that were not significant in atleast one feature
filt_ACC_scaled_p_df <- ACC_scaled_p_df[which(apply(ACC_scaled_p_df, 1, function(x) {any(x <= 0.1)})),]
dim(filt_ACC_scaled_p_df)


ASV_to_keep <- rownames(filt_ACC_scaled_p_df)

ACorn_cob_coef <- data.frame()
i <- 0
for(feat in features_to_test){
  i <- i +1
  if(feat=="A_SDC_GENDER"){
    feat <- "A_SDC_GENDER2"
  }
  feat_name <- paste("mu.",feat,sep="")
  for(j in 1:length(ASV_CC_res[[i]]$all_models)){
      ACorn_cob_coef[j,i] <- ASV_CC_res[[i]]$all_models[[j]]$coefficients[,1][which(rownames(ASV_CC_res[[i]]$all_models[[j]]$coefficients)==feat_name)]
    
  }
  
}

rownames(ACorn_cob_coef) <- rownames(ASV_CC_res[[1]]$data@otu_table)
colnames(ACorn_cob_coef) <- features_to_test

filt_ACorn_cob_coef <- ACorn_cob_coef[rownames(filt_ACC_scaled_p_df),]
renamed_filt_ACorn_cob_coef <- filt_ACorn_cob_coef



set_na_filt_ACorn_cob_coef <- renamed_filt_ACorn_cob_coef
set_na_filt_ACorn_cob_coef[filt_ACC_scaled_p_df > 0.1] <- 0

colnames(set_na_filt_ACorn_cob_coef) <- c("Age", "Height", "Waist Size", "Weight", "Waist Hip Ratio", "Fat Free Mass", "Salt Usage", "Nut/Seed Servings",
                                          "Vegetable Servings", "Refined Grain Servings", "Sleeping Light Exposure", "Sex", "Body Mass Index", 
                                          "Juice Servings", "Last Dental Visit")
set_na_filt_ACorn_cob_coef$Sex <- -set_na_filt_ACorn_cob_coef$Sex
colnames(set_na_filt_ACorn_cob_coef)[12] <- "Male*"

ASV_taxa_labels <- HealthyOralMicrobiome::Taxa_Labels

taxa_asv <- match(rownames(set_na_filt_ACorn_cob_coef), rownames(ASV_taxa_labels))
taxa_asv

### only keep first 5 ASV characters
rownames(set_na_filt_ACorn_cob_coef) <- substring(rownames(set_na_filt_ACorn_cob_coef), 1, 5)
rownames(set_na_filt_ACorn_cob_coef) <- paste0("ASV-", rownames(set_na_filt_ACorn_cob_coef))
## okay now need to paste lowest taxaonomic assignment

taxa_rename <- ASV_taxa_labels
taxa_rename$Taxon <- gsub(".*D_5", "D_5", taxa_rename$Taxon)
taxa_rename$Taxon <- gsub("D_5__", "", taxa_rename$Taxon)
taxa_rename$Taxon <- gsub("D_6__", "", taxa_rename$Taxon)
taxa_rename$Taxon <- gsub(";", " ", taxa_rename$Taxon)
rownames(set_na_filt_ACorn_cob_coef) <- paste(rownames(set_na_filt_ACorn_cob_coef), taxa_rename$Taxon[taxa_asv], sep=" ")
#manually curate the rest
rownames(set_na_filt_ACorn_cob_coef) <- gsub(" bacterium", "", rownames(set_na_filt_ACorn_cob_coef))
rownames(set_na_filt_ACorn_cob_coef)[1] <- "ASV-0c2af Treponema 2 socranskii subsp. socranskii"
rownames(set_na_filt_ACorn_cob_coef)[3] <- "ASV-10d24 Haemophilus haemolyticus"
rownames(set_na_filt_ACorn_cob_coef)[4] <- "ASV-18820 Alloprevotella uncultured"
rownames(set_na_filt_ACorn_cob_coef)[7] <- "ASV-1e8e5 Oribacterium sinus"
rownames(set_na_filt_ACorn_cob_coef)[10] <- "ASV-2f0df Megasphaera micronuciformis"
rownames(set_na_filt_ACorn_cob_coef)[12] <- "ASV-43ec6 Capnocytophaga gingivalis"
rownames(set_na_filt_ACorn_cob_coef)[14] <- "ASV-4ca02 Selenomonas uncultured"
rownames(set_na_filt_ACorn_cob_coef)[15] <- "ASV-4fcb9 Prevotella uncultured"
rownames(set_na_filt_ACorn_cob_coef)[16] <- "ASV-50e51 Prevotella sp. oral taxon 299 str. F0039"
rownames(set_na_filt_ACorn_cob_coef)[24] <- "ASV-77909 Prevotella 7 denticola" 
rownames(set_na_filt_ACorn_cob_coef)[30] <- "ASV-9a2e6 Mogibacterium unclassified"
rownames(set_na_filt_ACorn_cob_coef)[31] <- "ASV-a0a4a Lachnospiraceae uncultured"
rownames(set_na_filt_ACorn_cob_coef)[33] <- "ASV-b8b0d Capnocytophaga granulosa" 
rownames(set_na_filt_ACorn_cob_coef)[35] <- "ASV-c891c Bacteroidetes F0058 uncultured"
rownames(set_na_filt_ACorn_cob_coef)[36] <- "ASV-ca12e Streptococcus anginosus subsp. anginosus"
rownames(set_na_filt_ACorn_cob_coef)[42] <- "ASV-eec30 Prevotella unclassified"


ASV_phylums <- gsub("D_0__.*;D_1__", "", ASV_taxa_labels$Taxon[taxa_asv])
ASV_phylums <- gsub("D_2__.*","", ASV_phylums)
ASV_phylums <- gsub(";", "", ASV_phylums)
ASV_phylums


order_ASV_phylums <- ASV_phylums[order(ASV_phylums)]
row_reorder <- set_na_filt_ACorn_cob_coef[order(ASV_phylums),]

taxa_colors <- vector()
for(taxa in 1:length(order_ASV_phylums)){
  
  if(order_ASV_phylums[[taxa]] == "Actinobacteria"){
      taxa_colors[[taxa]] <- "Red"
  }else if(order_ASV_phylums[[taxa]]=="Bacteroidetes"){
    taxa_colors[[taxa]] <- "Blue"
  }else if(order_ASV_phylums[[taxa]]=="Cyanobacteria"){
    taxa_colors[[taxa]] <- "Orange"
  }else if(order_ASV_phylums[[taxa]]=="Firmicutes"){
    taxa_colors[[taxa]] <- "Brown"
  }else if(order_ASV_phylums[[taxa]]=="Fusobacteria"){
    taxa_colors[[taxa]] <- "Purple"
  }else if(order_ASV_phylums[[taxa]]=="Proteobacteria"){
    taxa_colors[[taxa]] <- "#664163"
  }else if(order_ASV_phylums[[taxa]]=="Spirochaetes"){
    taxa_colors[[taxa]] <- "Black"
  }else if(order_ASV_phylums[[taxa]]=="Synergistetes"){
    taxa_colors[[taxa]] <- "Cyan"
  }else if(order_ASV_phylums[[taxa]]=="Tenericutes"){
    taxa_colors[[taxa]] <- "Gray"
  }else if(order_ASV_phylums[[taxa]]=="Epsilonbacteraeota"){
    taxa_colors[[taxa]] <- "#42256f"
  }else
    print("error")
    
}

reorder_set_na_filt_ACorn_cob_coef <- row_reorder[,c(12,6,2,8,10,7,11,3,13,4,14,9,5,1,15)]


layout_mat <- rbind(c(4,0,0), c(2,1,0), c(0,3,0))
lwid=c(2,3.5,.00001)
lhei=c(.2,3,.3)


breaks=seq(-2, 2, length.out = 200) #41 values
#now add outliers
breaks=append(breaks, 10000)
breaks=append(breaks, -10000, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="blue",mid="white",high="brown")

ASV_heat <- heatmap.2(as.matrix(reorder_set_na_filt_ACorn_cob_coef),na.color="black", density.info = "none", symkey = F, trace="none", srtCol=35, dendrogram = "none", key.xlab = "log odds", key.ylab = "p Value", key=T, col=mycol, lmat = layout_mat, lwid = lwid, lhei = lhei, cexRow = 1.3, cexCol = 1.3, Colv=F, distfun=function(x) dist(x, method="manhattan"), colRow = taxa_colors, Rowv=F, offsetRow = -56.5, sepcolor = "black", rowsep = c(7, 20, 35, 37,  41), keysize=0.2, key.par = list(cex=0.01, cex.axis=60, cex.main=60, mar=c(50,50,50,50), mgp=c(.01,50,1), tcl=-10))


### done

