#' @title Integrate miRNA (MicroRNA) data.
#'
#' @description Add miRNA data to MAF file.
#'
#' @details For each variant, this function checks its genomic position against known miRNA target sites.
#' When an overlap is found, it appends the matched miRNA name(s) and target site details to the variant.
#'
#' @param maf Required. MAF data frame (required columns: Chromosome, Start_Position, End_Position).
#' @param mirna_target Data frame contain miRNA data (required columns: Chromosome, Start_Position, End_Position, miRNA)
#' 
#' @param projection The genome build projection for the variants you are working with (default is grch37)
#'
#' @return data frame in MAF format with mirna (and sites := length of the seed region) columns.
#'
#' @import dplyr tidyr Rcpp GAMBLR.data
#'
#' @export
#' 
#' \dontrun{
#' @examples
#' sample_df = target_targetscan (maf, mirna_target)
#'}



target_targetscan <- function(
    maf,
    mirna_target
){
    if (missing(mirna_target)){
      mirna_target <- GAMBLR.data::mirna_targetscan
    }
    if (!"sites" %in% names(mirna_target)){
      mirna_target$sites <- NA_integer_
    }  
    if (all(grepl("^chr", maf$Chromosome))) {
      mirna_target$Chromosome <- gsub("chr", "", mirna_target$Chromosome) 
      mirna_target$Chromosome <- paste0("chr", mirna_target$Chromosome)
    } else {
      # If there is a mix of prefixed and non-prefixed options
      maf$Chromosome <- gsub("chr", "", maf$Chromosome) 
      maf$Chromosome <- paste0("chr", maf$Chromosome)
      mirna_target$Chromosome <- gsub("chr", "", mirna_target$Chromosome) 
      mirna_target$Chromosome <- paste0("chr", mirna_target$Chromosome)
    }
    
    sourceCpp("./targetscan_core_code.cpp")
    
    targetscan_res = targetscan_rcpp(maf, mirna_target)
    
    final_maf = left_join(maf, targetscan_res, relationship = "many-to-many") %>% distinct()
    
    return(final_maf)
}
  