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
#' suppressPackageStartupMessages(library(GAMBLR.open))
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
    warning("No metadata supplied. Returning SSMs for all available samples. Supply these_samples_metadata to limit results to samples matching desired clinical features.")
    these_samples_metadata = suppressMessages(
        get_gambl_metadata(seq_type_filter = c("genome","capture"))) %>%
        GAMBLR.helpers::check_and_clean_metadata(.,duplicate_action = "keep_first")
  }
  these_samples_metadata = dplyr::filter(these_samples_metadata, seq_type != "mrna")
  capture_ids = dplyr::filter(these_samples_metadata,seq_type=="capture") %>%
    pull(sample_id)
  genome_ids = dplyr::filter(these_samples_metadata,seq_type=="genome") %>%
    pull(sample_id)

  # Return coding SSMs from the slms-3 pipeline, honouring `include_silent`.
  # (Historically this function applied no Pipeline or coding-class filter despite
  # its name, returning redundant rows from every pipeline; that is now fixed.)
  make_maf <- function(ids, seqtype){
    if(length(ids) == 0) return(NULL)
    GAMBLR.data::get_ssm_from_db(projection = projection,
                                 sample_ids = ids,
                                 tool_name = "slms-3",
                                 coding_only = TRUE,
                                 include_silent = include_silent) %>%
      GAMBLR.utils::create_maf_data(projection) %>%
      mutate(maf_seq_type = seqtype)
  }
  capture_maf = make_maf(capture_ids, "capture")
  genome_maf  = make_maf(genome_ids, "genome")

  if(!is.null(capture_maf) && !is.null(genome_maf)){
    return(GAMBLR.utils::bind_genomic_data(capture_maf, genome_maf))
  }else if(!is.null(capture_maf)){
    return(capture_maf)
  }else if(!is.null(genome_maf)){
    return(genome_maf)
  }
}
