#' @title Get all Coding SSMs
#'
#' @description Retrieve all coding SSMs from GAMBL in MAF-like format, regardless of seq_type.
#'
#' @details Effectively retrieve coding SSM calls from one or all DNA seq_type. For additional optional arguments, see [GAMBLR.results::get_coding_ssm]

#' @param these_samples_metadata Supply a metadata table containing the sample/seq_type combinations you want. 
#' @param include_silent If set to TRUE, silent/synonymous mutations in the coding regions will also be returned. 
#' @param projection The desired genome build
#' 'grch37' or 'hg38' are allowed. Default is grch37
#'
#' @return A data frame containing all the MAF data columns (one row per mutation).
#'
#' @import dplyr tidyr GAMBLR.helpers GAMBLR.utils
#' @export
#'
#' @examples
#' library(GAMBLR.open)
#' my_meta = get_gambl_metadata(seq_type_filter = c("genome","capture"))
#' my_meta = check_and_clean_metadata(my_meta,duplicate_action="keep_first")
#' maf_all_seqtype = get_all_coding_ssm(my_meta)
#' 
#' table(maf_all_seqtype$maf_seq_type)
#' 
#' # most common mutations by gene and Variant_Classification
#' dplyr::group_by(maf_all_seqtype,
#'                 Hugo_Symbol,
#'                 Variant_Classification) %>% 
#'   dplyr::count() %>% 
#'   dplyr::arrange(desc(n))
get_all_coding_ssm = function(these_samples_metadata = NULL,
                              include_silent=FALSE,
                              projection = "grch37"){
  if(missing(these_samples_metadata)){
    these_samples_metadata = suppressMessages(
        get_gambl_metadata(seq_type_filter = c("genome","capture"))) %>%
        GAMBLR.helpers::check_and_clean_metadata(.,duplicate_action = "keep_first")
  }
  capture_ids = dplyr::filter(these_samples_metadata,seq_type=="capture") %>%
    pull(sample_id)
  genome_ids = dplyr::filter(these_samples_metadata,seq_type=="genome") %>%
    pull(sample_id)

  capture_maf = GAMBLR.data::sample_data[[projection]]$maf %>%
    dplyr::filter(Tumor_Sample_Barcode %in% capture_ids) %>%
    GAMBLR.utils::create_maf_data(.,projection) %>%
    mutate(.,maf_seq_type = "capture")
  genome_maf = GAMBLR.data::sample_data[[projection]]$maf  %>%
    dplyr::filter(Tumor_Sample_Barcode %in% genome_ids) %>%
    GAMBLR.utils::create_maf_data(.,projection) %>%
    mutate(.,maf_seq_type = "genome")

  if(length(capture_ids)>1 && length(genome_ids) > 1){
    merged_ssm = GAMBLR.utils::bind_genomic_data(capture_maf,genome_maf)
    return(merged_ssm)
  }else if(length(capture_ids)>1){
    return(capture_maf)
  }else if(length(genome_ids) > 1){
    return(genome_maf)
  }
}
