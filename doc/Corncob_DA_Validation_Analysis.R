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

## ---- eval=F, load_data_in-----------------------------------------------
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
#  Valid_ASV_data <- Filt_samples(depth=5000, setType="validation", Taxa_table = ASV_tab, Metadata = Metadata, parallel = 20, cutoff_prop = 0.1)
#  Valid_Genus_data <- Filt_samples(depth=5000, setType="validation", Taxa_table = Genus_tab, Metadata = Metadata, parallel = 20, cutoff_prop = 0.1)
#  dim(Valid_ASV_data[["Taxa_Tab"]])
#  dim(Valid_Genus_data[["Taxa_Tab"]])
#  
#  ### finally we will scale the data so that the values will be comparable to one another
#  
#  
#  scaled_metadata_g<- Valid_Genus_data[["Metadata"]] %>% rownames_to_column('sample_name') %>%
#    mutate_if(is.numeric, scale) %>%
#    column_to_rownames('sample_name')
#  
#  scaled_metadata_a <- Valid_ASV_data[["Metadata"]] %>% rownames_to_column('sample_name') %>%
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
#  

## ---- eval=F, Run_Corncob_valid_tests------------------------------------
#  
#  ASVs_valid <- phyloseq::otu_table(Valid_ASV_data[['Taxa_Tab']], taxa_are_rows = T)
#  Meta_valid <- phyloseq::sample_data(scaled_metadata_a)
#  phylo_valid <- phyloseq::phyloseq(ASVs_valid, Meta_valid)
#  
#  
#  ## ASV Analysis
#  cl <- parallel::makeCluster(15)
#  doParallel::registerDoParallel(cl)
#  foreach::getDoParWorkers()
#  `%dopar%` <- foreach::`%dopar%`
#  Valid_ASV_CC_Res <- foreach::foreach(i=1:15) %dopar% HealthyOralMicrobiome::valid_corn_cob_analysis(features_to_test[i], scaled_metadata_a, taxa_tab = Valid_ASV_data[["Taxa_Tab"]])
#  parallel::stopCluster(cl)
#  
#  
#  ## Genus Analysis
#  
#  cl <- parallel::makeCluster(15)
#  doParallel::registerDoParallel(cl)
#  foreach::getDoParWorkers()
#  `%dopar%` <- foreach::`%dopar%`
#  Valid_Genus_CC_Res <- foreach::foreach(i=1:15) %dopar% HealthyOralMicrobiome::valid_corn_cob_analysis(features_to_test[i], scaled_metadata_g, taxa_tab = Valid_Genus_data[["Taxa_Tab"]])
#  parallel::stopCluster(cl)
#  
#  

## ---- eval=T, fig.width=10, fig.height=10, validation_heatmap------------
Genus_valid_CC_res <- HealthyOralMicrobiome::Valid_Genus_CC_Res
Genus_Comp_CC_res <- HealthyOralMicrobiome::Comp_Genus_CC_Res

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

Comp_scaled_Genus_p <- list()

for(i in 1:15){
  
  
  Comp_scaled_Genus_p[[i]] <- Genus_Comp_CC_res[[i]]$p_fdr
}


Comp_CC_p_df <- do.call(rbind, Comp_scaled_Genus_p)
Comp_CC_p_df <- data.frame(t(Comp_CC_p_df))
colnames(Comp_CC_p_df) <- features_to_test

#remove any taxa that weren't significant in atleast on feature
filt_Comp_scaled_p_df <- Comp_CC_p_df[which(apply(Comp_CC_p_df, 1, function(x) {any(x <= 0.1)})),]
dim(filt_Comp_scaled_p_df)

### get validation 
Genus_valid_cc_p <- list()
for(i in 1:15){
   Genus_valid_cc_p[[i]] <- Genus_valid_CC_res[[i]]$p
  
}

Genus_valid_cc_p_df <- do.call(rbind, Genus_valid_cc_p)
Genus_valid_cc_p_df <- data.frame(t(Genus_valid_cc_p_df))

#first only keep taxa that were significant in original analysis

filt_Genus_valid_cc_p_df <- Genus_valid_cc_p_df[rownames(filt_Comp_scaled_p_df),]
dim(filt_Genus_valid_cc_p_df)
### great

## grab coef for these taxa

Corn_cob_coef_valid_genus <- data.frame()
i <- 0
for(feat in features_to_test){
  i <- i +1
  if(feat=="A_SDC_GENDER"){
    feat <- "A_SDC_GENDER2"
  }
  feat_name <- paste("mu.",feat,sep="")
  for(j in 1:length(Genus_valid_CC_res[[i]]$all_models)){

    #set value to NA if p-value could not be calculated due to the low abundance of that organism (discriminant taxa)
    if(is.na(Genus_valid_CC_res[[i]]$all_models[[j]][1])){
      Corn_cob_coef_valid_genus[j,i] <- NA
      
    }else{
      Corn_cob_coef_valid_genus[j,i] <- Genus_valid_CC_res[[i]]$all_models[[j]]$coefficients[,1][which(rownames(Genus_valid_CC_res[[i]]$all_models[[j]]$coefficients)==feat_name)]
    }
  }
}

rownames(Corn_cob_coef_valid_genus) <- rownames(Genus_valid_CC_res[[1]]$data@otu_table)
colnames(Corn_cob_coef_valid_genus) <- features_to_test

Corn_cob_coef_valid_genus$A_SDC_GENDER <- -Corn_cob_coef_valid_genus$A_SDC_GENDER

filt_genus_coef <- Corn_cob_coef_valid_genus[rownames(filt_Genus_valid_cc_p_df),]
dim(filt_genus_coef)

### set coef to 0 if not significant in first analysis
identical(rownames(filt_genus_coef), rownames(filt_Comp_scaled_p_df))
identical(colnames(filt_genus_coef), colnames(filt_Comp_scaled_p_df))

filt_genus_coef[filt_Comp_scaled_p_df >= 0.1] <- 0

### set coef to 0 if valid p value is NA

filt_genus_coef[is.na(filt_Genus_valid_cc_p_df)] <- 0

### set coef to 0 if not significant in validation analysis

filt_genus_coef[filt_Genus_valid_cc_p_df >= 0.05] <- 0

renamed_filt_genus_coef <- filt_genus_coef
rownames(renamed_filt_genus_coef) <- gsub(".*D_4", "D_4", rownames(renamed_filt_genus_coef))
rownames(renamed_filt_genus_coef) <- gsub("D_4__", "", rownames(renamed_filt_genus_coef))
rownames(renamed_filt_genus_coef) <- gsub("D_5__", "", rownames(renamed_filt_genus_coef))
rownames(renamed_filt_genus_coef) <- gsub(";", " ", rownames(renamed_filt_genus_coef))

rownames(renamed_filt_genus_coef)[29] <- "Veillonellaceae unclassified"
rownames(renamed_filt_genus_coef)[34] <- "Neisseriaceae unclassified"
rownames(renamed_filt_genus_coef)[41] <- "Mollicutes uncultured"


### remove columns with all 0
remove_cols <- renamed_filt_genus_coef

remove_cols <- remove_cols[,colSums(remove_cols) != 0]
remove_cols

#set column names
colnames(remove_cols) <- c("Age", "Height", "Waist Hip Ratio", "Fat Free Mass", "Male*")

Phylums <- gsub("D_0__.*;D_1__", "", rownames(filt_genus_coef))
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

breaks=seq(-3, 3, length.out = 180) #41 values
#now add outliers
breaks=append(breaks, 10000)
breaks=append(breaks, -10000, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="blue",mid="white",high="brown")

layout_mat <- rbind(c(0,0,4), c(2,1,3), c(0,0,0))
lwid=c(.1,2.75,1.25)
lhei=c(.2,3,.3)

gen_heat <- heatmap.2(as.matrix(remove_cols),
          na.color="black", density.info = "none", symkey = F, trace="none", srtCol=35, dendrogram = "none", key.xlab = "log odds", key.ylab = "p Value", key=T, col=mycol, lmat = layout_mat, lwid = lwid, lhei = lhei, cexRow = 1.3, cexCol = 1.3, distfun=function(x) dist(x, method="manhattan"), colRow = taxa_colors, Rowv = F,
          rowsep = c(2, 10,29,30, 38, 39, 40), sepcolor = "black", keysize=0.2, key.par = list(cex=0.01, cex.axis=60, cex.main=60, mar=c(50,0,50,10), mgp=c(.01,50,1), tcl=-10), Colv=F)


## ---- eval=T, fig.width=10, fig.height=10, validation_ASV_heatmap--------

ASV_valid_CC_res <- HealthyOralMicrobiome::Valid_ASV_CC_Res
ASV_Comp_CC_res <- HealthyOralMicrobiome::Comp_ASV_CC_Res


ASV_Comp_p <- list()
for(i in 1:15){
  
  
  ASV_Comp_p[[i]] <- ASV_Comp_CC_res[[i]]$p_fdr
}

ASV_Comp_p_df <- do.call(rbind, ASV_Comp_p)
ASV_Comp_p_df <- data.frame(t(ASV_Comp_p_df))
colnames(ASV_Comp_p_df) <- features_to_test

ASV_filt_Comp_p_df <- ASV_Comp_p_df[which(apply(ASV_Comp_p_df, 1, function(x) {any(x <= 0.1)})),]
dim(ASV_filt_Comp_p_df)


ASV_Valid_p <- list()
for(i in 1:15){
  
  ASV_Valid_p[[i]] <- ASV_valid_CC_res[[i]]$p
  
}

ASV_Valid_p_df <- do.call(rbind, ASV_Valid_p)
ASV_Valid_p_df <- data.frame(t(ASV_Valid_p_df))
colnames(ASV_Valid_p_df) <- features_to_test

filt_ASV_Valid_p_df <- ASV_Valid_p_df[rownames(ASV_filt_Comp_p_df),]

Corn_cob_coef_valid_ASV <- data.frame()
i <- 0
for(feat in features_to_test){
  i <- i +1
  if(feat=="A_SDC_GENDER"){
    feat <- "A_SDC_GENDER2"
  }
  feat_name <- paste("mu.",feat,sep="")
  for(j in 1:length(ASV_valid_CC_res[[i]]$all_models)){

    #set value to NA if p-value could not be calculated due to the low abundance of that organism (discriminant taxa)
    if(is.na(ASV_valid_CC_res[[i]]$all_models[[j]][1])){
      Corn_cob_coef_valid_ASV[j,i] <- NA
      
    }else{
      Corn_cob_coef_valid_ASV[j,i] <- ASV_valid_CC_res[[i]]$all_models[[j]]$coefficients[,1][which(rownames(ASV_valid_CC_res[[i]]$all_models[[j]]$coefficients)==feat_name)]
    }
  }
}

rownames(Corn_cob_coef_valid_ASV) <- rownames(ASV_valid_CC_res[[1]]$data@otu_table)
colnames(Corn_cob_coef_valid_ASV) <- features_to_test

Corn_cob_coef_valid_ASV$A_SDC_GENDER <- -Corn_cob_coef_valid_ASV$A_SDC_GENDER

filt_ASV_coef <- Corn_cob_coef_valid_ASV[rownames(filt_ASV_Valid_p_df),]
dim(filt_ASV_coef)


### set coef to 0 if not significant in first analysis
identical(rownames(filt_ASV_coef), rownames(ASV_filt_Comp_p_df))
identical(colnames(filt_ASV_coef), colnames(ASV_filt_Comp_p_df))

filt_ASV_coef[ASV_filt_Comp_p_df >= 0.1] <- 0

### set coef to 0 if valid p value is NA

filt_ASV_coef[is.na(filt_ASV_Valid_p_df)] <- 0

### set coef to 0 if not significant in validation analysis

filt_ASV_coef[filt_ASV_Valid_p_df >= 0.05] <- 0

#remove columns that are all 0
ASV_remove_Cols <- filt_ASV_coef

ASV_remove_Cols <- ASV_remove_Cols[,colSums(ASV_remove_Cols) != 0]

ASV_taxa_labels <- HealthyOralMicrobiome::Taxa_Labels

taxa_asv <- match(rownames(ASV_remove_Cols), rownames(ASV_taxa_labels))
taxa_asv
renamed_ASV_coefs <- ASV_remove_Cols

rownames(renamed_ASV_coefs) <- substring(rownames(renamed_ASV_coefs), 1, 5)
rownames(renamed_ASV_coefs) <- paste0("ASV-", rownames(renamed_ASV_coefs))

taxa_rename <- ASV_taxa_labels
taxa_rename$Taxon <- gsub(".*D_5", "D_5", taxa_rename$Taxon)
taxa_rename$Taxon <- gsub("D_5__", "", taxa_rename$Taxon)
taxa_rename$Taxon <- gsub("D_6__", "", taxa_rename$Taxon)
taxa_rename$Taxon <- gsub(";", " ", taxa_rename$Taxon)

rownames(renamed_ASV_coefs) <- paste(rownames(renamed_ASV_coefs), taxa_rename$Taxon[taxa_asv], sep=" ")
rownames(renamed_ASV_coefs)

rownames(renamed_ASV_coefs) <- gsub(" bacterium", "", rownames(renamed_ASV_coefs))

rownames(renamed_ASV_coefs)[1] <- "ASV-0c2af Treponema 2 socranskii subsp. socranskii"
rownames(renamed_ASV_coefs)[3] <- "ASV-10d24 Haemophilus haemolyticus"
rownames(renamed_ASV_coefs)[4] <- "ASV-18820 Alloprevotella uncultured"
rownames(renamed_ASV_coefs)[7] <- "ASV-1e8e5 Oribacterium sinus"
rownames(renamed_ASV_coefs)[10] <- "ASV-2f0df Megasphaera micronuciformis"
rownames(renamed_ASV_coefs)[12] <- "ASV-43ec6 Capnocytophaga gingivalis"
rownames(renamed_ASV_coefs)[14] <- "ASV-4ca02 Selenomonas uncultured"
rownames(renamed_ASV_coefs)[15] <- "ASV-4fcb9 Prevotella uncultured"
rownames(renamed_ASV_coefs)[16] <- "ASV-50e51 Prevotella sp. oral taxon 299 str. F0039"
rownames(renamed_ASV_coefs)[24] <- "ASV-77909 Prevotella 7 denticola" 
rownames(renamed_ASV_coefs)[30] <- "ASV-9a2e6 Mogibacterium unclassified"
rownames(renamed_ASV_coefs)[31] <- "ASV-a0a4a Lachnospiraceae uncultured"
rownames(renamed_ASV_coefs)[33] <- "ASV-b8b0d Capnocytophaga granulosa" 
rownames(renamed_ASV_coefs)[35] <- "ASV-c891c Bacteroidetes F0058 uncultured"
rownames(renamed_ASV_coefs)[36] <- "ASV-ca12e Streptococcus anginosus subsp. anginosus"
rownames(renamed_ASV_coefs)[42] <- "ASV-eec30 Prevotella unclassified"

ASV_phylums <- gsub("D_0__.*;D_1__", "", ASV_taxa_labels$Taxon[taxa_asv])
ASV_phylums <- gsub("D_2__.*","", ASV_phylums)
ASV_phylums <- gsub(";", "", ASV_phylums)
ASV_phylums

order_ASV_phylums <- ASV_phylums[order(ASV_phylums)]
row_reorder <- renamed_ASV_coefs[order(ASV_phylums),]

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

layout_mat <- rbind(c(4,0,0), c(2,1,0), c(0,3,0))
lwid=c(2,3.5,.00001)
lhei=c(.2,3,.3)


breaks=seq(-2, 2, length.out = 200) #41 values
#now add outliers
breaks=append(breaks, 10000)
breaks=append(breaks, -10000, 0)
#create colour panel with length(breaks)-1 colours
mycol <- colorpanel(n=length(breaks)-1,low="blue",mid="white",high="brown")

colnames(row_reorder) <- c("Height", "Fat Free Mass", "Sleeping Light Exposure", "Male*")

ASV_heat <- heatmap.2(as.matrix(row_reorder),na.color="black", density.info = "none", symkey = F, trace="none", srtCol=35, dendrogram = "none", key.xlab = "log odds", key.ylab = "p Value", key=T, col=mycol, lmat = layout_mat, lwid = lwid, lhei = lhei, cexRow = 1.3, cexCol = 1.3, Colv=F, distfun=function(x) dist(x, method="manhattan"), colRow = taxa_colors, Rowv=F, offsetRow = -56.5, sepcolor = "black", rowsep = c(7, 20, 35, 37,  41), keysize=0.2, key.par = list(cex=0.01, cex.axis=60, cex.main=60, mar=c(50,50,50,50), mgp=c(.01,50,1), tcl=-10))


