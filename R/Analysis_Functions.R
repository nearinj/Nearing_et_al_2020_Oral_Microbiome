##### Contains all functions needed to run analysis found in vignettes


#' Filt_samples function
#'
#' This function takes a sequencing depth and whether you want to return the whole dataset, the completed cases dataset, or
#' the dataset that was used for validation.
#'
#'@param depth The depth that you want to filter samples to. Defaults to 5000
#'@param setType The dataset that type you want to return. Currently there are
#' three options "complete" returns the completed cases dataset (main cohort),
#' "validation" returns the validation cohort, "whole" returns all samples in the dataset.
#' @param Taxa_table The ASV or genus table that is to be filtered.
#' @param Metadata The metadata used for filtering.
#' @param cutoff_prop The cutoff proportion for fitlering rare taxa.
#' @param parallel The number of threads you want to use to filter your data.
#' @return Returns both the filtered metadata and taxa_table in a list. Filtering is done as follows:
#' Samples are filtered by read depth, rare taxa removed, and then samples filtered based on setType.
#'

Filt_samples <- function(depth=5000, setType="complete",Taxa_table, Metadata, cutoff_prop=0.1, parallel=1){
  ret_list <- list()

  if(cutoff_prop <= 0){
    int_taxa_tab <- Taxa_table

  }else{
    int_taxa_tab <- Taxa_table[,colSums(Taxa_table)>=depth]
  }


  int_taxa_tab_rare <- remove_rare_features(table = int_taxa_tab, parallel = parallel)

  int_metadata <- Metadata[colnames(int_taxa_tab_rare),]

  if(setType=="complete"){
      int_metadata2 <- int_metadata[complete.cases(int_metadata),]

  }else if(setType=="whole"){
    int_metadata2 <- int_metadata
  }else if(setType=="validation"){
    int_metadata2 <- int_metadata[!complete.cases(int_metadata),]
  }else{
    stop("Invalid SetType")
  }

  fin_taxa_tab <- int_taxa_tab_rare[,rownames(int_metadata2)]


  ret_list[["Metadata"]] <- int_metadata2
  ret_list[["Taxa_Tab"]] <- fin_taxa_tab
  if(!identical(rownames(int_metadata2), colnames(fin_taxa_tab))){
    stop("Something went wrong during filtering returning NULL")
  }
  return(ret_list)
}



#' get_alpha_divs Function
#' Takes in metadata and returns the assocaited alpha diversity values for each sample.
#' @param Metadata Input Metadata that is associated with alpha diversity
#' @param setType The set of alpha diversity values you want returned. Can either be "complete"
#' "validation" or "whole"
#'
#' @return Returns a dataframe containing the sample metadata and its related alpha diversity metrics
#'
#'
#'
#'
get_alpha_divs <- function(Metadata, setType="complete"){

  Faiths <- HealthyOralMicrobiome::Healthy_Faiths
  Shannon <- HealthyOralMicrobiome::Healthy_Shannon
  Evenness <- HealthyOralMicrobiome::Healthy_Evenness
  Richness <- HealthyOralMicrobiome::Healthy_Richness

  int_metadata <- Metadata[rownames(Faiths),]


  if(setType=="complete"){
    int_metadata2 <- int_metadata[complete.cases(Metadata),]
  }else if(setType=="whole"){
    int_metadata2 <- int_metadata
  }else if(setType=="validation"){
    int_metadata2 <- int_metadata[!complete.cases(Metadata),]
  }else{
    stop("Invalid SetType")
  }


  int_metadata2$faiths <- Faiths[rownames(int_metadata2),]
  int_metadata2$shannon <- Shannon[rownames(int_metadata2),]
  int_metadata2$evenness <- Evenness[rownames(int_metadata2),]
  int_metadata2$richness <- Richness[rownames(int_metadata2),]

  return(int_metadata2)
}


#' remove_rare_features function
#'
#' This function takes in a feature table and removes features that are found in less than the specififed cutoff
#' @param table The feature table that needs to be filtered
#' @param cutoff_pro The cutoff proportion used to identify rare features. Defaults to 0.1
#' @param parrellel The number of threads used to perform filtering. Defaults to 1.
#'
#' @return Returns the filtered feature table
#'

remove_rare_features <- function( table , cutoff_pro=0.1, parallel=1 ) {
  if(cutoff_pro==0){
    message("No filtering will be done due to cutoff_pro set to 0")
    return(table)
  }
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )
  if(parallel <= 1){
    for ( i in 1:nrow(table) ) {
      row_nonzero <- length( which( table[ i , ]  > 0 ) )
      if ( row_nonzero > cutoff ) {
        row2keep <- c( row2keep , i)
      }
    }
    return( table [ row2keep , , drop=F ])
  }else{
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(cl)
    message("Running with: ", foreach::getDoParWorkers(), " cores")

    res <- foreach::foreach(i = 1:nrow(table), .combine = c) %dopar% {
      row_nonzero <- length( which ( table[ i , ] > 0))
      if(row_nonzero > cutoff){
        i
      }else
        NULL
    }
    parallel::stopCluster(cl)
    return( table[res, , drop=F])

  }

}


#' run_corn_cob_analysis Function
#' Runs CornCob differential abundance analysis on the input feature and phyloseq object
#' Automatically runs DA of the feature while controlling for DNA extraction round and differiental variablility
#' @param feature Feature to be tested
#' @param phylo The phyloseq object that contains the metadata and ASV table to be tested.
#' @return Returns a list containing the DA results from Corncob
#'
#'
run_corn_cob_analysis <- function(feature, phylo){


  formula <- as.formula(paste("~ Extraction_Number",feature, sep="+"))

  null_formula <- as.formula("~ Extraction_Number")


  phi_formula <- as.formula(paste("~ Extraction_Number",feature, sep="+"))

  message(formula)

  message(null_formula)

  message(phi_formula)
  res <- corncob::differentialTest(formula = formula,
                          formula_null = null_formula,
                          phi.formula = phi_formula,
                          phi.formula_null = phi_formula,
                          data=phylo, test="Wald", fdr_cutoff = 0.1, boot=F)

  return(res)

}

#' run_corn_cob_analysis Function
#' Runs CornCob differential abundance analysis on the input feature and phyloseq object while removing samples that have NA data for each feature.
#' Automatically runs DA of the feature while controlling for DNA extraction round and differiental variablility
#' @param feature Feature to be tested
#' @param phylo The phyloseq object that contains the metadata and ASV table to be tested.
#' @return Returns a list containing the DA results from Corncob
#'
#'
valid_corn_cob_analysis <- function(feature, df, taxa_tab){

  ### need to get complete tab

  feat_df <- df[complete.cases(df[,feature]),]
  feat_taxa <- taxa_tab[,rownames(feat_df)]

  otu_tab <- phyloseq::otu_table(feat_taxa, taxa_are_rows = T)
  sam_data <- phyloseq::sample_data(feat_df)
  phylo <- phyloseq::phyloseq(otu_tab, sam_data)

  formula <- as.formula(paste("~ Extraction_Number",feature, sep="+"))

  null_formula <- as.formula("~ Extraction_Number")


  phi_formula <- as.formula(paste("~ Extraction_Number",feature, sep="+"))

  message(formula)

  message(null_formula)

  message(phi_formula)
  res <- corncob::differentialTest(formula = formula,
                          formula_null = null_formula,
                          phi.formula = phi_formula,
                          phi.formula_null = phi_formula,
                          data=phylo, test="Wald", fdr_cutoff = 0.1, boot=F)

  return(res)


}

#' get_beta_div Function
#' This function generates a list containing the rarified distance metric and its associated metadata table.
#' @param metric Can either be "w_unifrac" to get a Weighted UniFrac distance matrix or "bray_curtis" to get
#' a Bray Curtis dissimilarity matrix.
#' @param setType Choose between "complete" to get the original dataset. Whole to get all samples... although some will be removed as
#' they lack a depth of 5000 reads. Or validation to get the validation dataset.
#' @param The metadata table associated with the dataset.
#' @return Returns a list containing the "Metadata" and "Distance" matrix
#'
#'
#'
#'
#'

get_beta_div <- function(metric, setType, metadata){


  if(metric=="w_unifrac"){

    dis=HealthyOralMicrobiome::Healthy_W_Unifrac

  }else if(metric=="bray_curtis"){

    dis=HealthyOralMicrobiome::Healthy_Bray_Curt
  }else{
    stop("Please choose between w_unifrac or bray_curtis as your distance metric")
  }

  int_metadata <- metadata[rownames(dis),]

  if(setType=="complete"){
    int_metadata2 <- int_metadata[complete.cases(int_metadata),]

  }else if(setType=="whole"){
    int_metadata2 <- int_metadata
  }else if(setType=="validation"){
    int_metadata2 <- int_metadata[!complete.cases(int_metadata),]
  }else{
    stop("Invalid SetType")
  }

  ret_dis <- dis[rownames(int_metadata2), rownames(int_metadata2)]

  ret_list <- list()

  ret_list[["Metadata"]] <- int_metadata2
  ret_list[["Distance"]] <- ret_dis

  return(ret_list)

}


#' get_core_features function
#'
#' This function takes in a feature table and a percent cutoff and returns the rows that correspond to features
#' that are found in more than the cutoff proportion of samples.
#'
#' @param Feat_tab A feature table  with rownames as features and colnames and samples
#' @param cutoff A precent cutoff for the % of samples
#'
#' @return Returns the rows corresponding to features that are present in above the cutoff % of total samples
#'
#'
#'
#'
get_core_features <- function(Feat_tab, cutoff){
  #convert Feat_tab to presnece absence
  pre_Feat_tab <- Feat_tab
  pre_Feat_tab[pre_Feat_tab>0] <- 1
  #calc cuttoff_num

  cutoff_num <- cutoff*length(colnames(Feat_tab))


  core_rows <- apply(pre_Feat_tab, 1, function(Feat_tab) sum(Feat_tab) > cutoff_num)

  return(core_rows)

}


#' run_alpha_compare Function
#' This function takes in a metadata feature to test, a type of alpha diversity and the dataframe containing the metadata and the alpha div metric.
#' It then compares a baseline linear model that only contains DNA extraction and y ~ DNA.Extraction and a second model containing both DNA.Extraction
#' and the metadata feature using an ANOVA. Finally the p-value for the metadata feature from this anova is returned.
#'
#' @param df A dataframe containing the columns "faiths" "shannon" "evenness" and "richness" along with the metadata feature that will be tested.
#' @param feature The metadata feature that will be tested.
#' @param type The type of alpha diversity that you would like to test. The current options are "faiths", "richness", "evenness" and "shannon".
#'
#' @return Returns a p-value from an anova comparing a model that contains only dna extraction and a model that contains both dna extraction and the metadata
#' feature of interest.
#'
#'
#'
#'
run_alpha_compare <- function(df, feature, type){

  message("Running alpha div analysis for ", type, " on ", feature)
  #find NA's for that feature
  remove_feats <- which(is.na(df[,feature]) | df[,feature]==-1)

  #remove NA's
  if(length(remove_feats) != 0){
    df <- df[-remove_feats,]
  }
  if(type=="faiths"){
    faiths_base <- lm(df[,"faiths"] ~ df[,"Extraction_Number"])
    faiths_feat <- lm(df[,"faiths"] ~ df[,"Extraction_Number"] + df[,feature])
    anova_faiths <- anova(faiths_base, faiths_feat)
    return(anova_faiths$`Pr(>F)`[2])
  }else if(type=="richness"){
    richness_base <- lm(df[,"richness"] ~ df[,"Extraction_Number"])
    richness_feat <- lm(df[,"richness"] ~ df[,"Extraction_Number"] + df[,feature])
    anova_rich <- anova(richness_base, richness_feat)
    return(anova_rich$`Pr(>F)`[2])
  }else if(type=="shannon"){
    shannon_base <- lm(df[,"shannon"] ~ df[,"Extraction_Number"])
    shannon_feat <- lm(df[,"shannon"] ~ df[,"Extraction_Number"] + df[,feature])
    anova_shannon <- anova(shannon_base, shannon_feat)
    return(anova_shannon$`Pr(>F)`[2])
  }else if(type=="evenness"){
    evenness_base <- lm(df[,"evenness"] ~ df[,"Extraction_Number"])
    evenness_feat <- lm(df[,"evenness"] ~ df[,"Extraction_Number"] + df[,feature])
    anova_evenness <- anova(evenness_base, evenness_feat)
    return(anova_evenness$`Pr(>F)`[2])
  }else
    message("No type selected?!")
  return(NA)
}
