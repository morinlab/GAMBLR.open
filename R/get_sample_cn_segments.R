#' @title Get Sample CN Segments.
#'
#' @description Get all segments for a single (or multiple) sample_id(s).
#'
#' @details This function returns CN segments. This works for single sample or multiple samples.
#' Specify the sample IDs you are interested in with `these_sample_ids` (as a vector of characters),
#' Or call this function with `these_samples_metadata` if you already have a metadata table subset to the sample IDs of interest.
#' If none of the above parameters are specified, the function will return CN segments for available samples (from get_gambl_metadata).
#' Note, this. function internally calls [id_ease] for dealing with sample IDs and metadata tables.
#'
#' @param these_sample_ids Optional, a vector of multiple sample_id (or a single sample ID as a string) that you want results for.
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' If not provided (and if `these_sample_ids` is not provided), the function will return all samples from the specified seq_type in the metadata.
#' @param projection Selected genome projection for returned CN segments. Default is "grch37".
#' @param this_seq_type Seq type for returned CN segments. Default is genome.
#' @param with_chr_prefix Set to TRUE to add a chr prefix to chromosome names. Default is FALSE.
#' @param streamlined Return a minimal output rather than full details. Default is FALSE.
#' @param verbose Set to FALSE to minimize the output to console. Default is TRUE. This parameter also dictates the verbosity of any helper function internally called inside the main function.
#' @param ... Any additional parameters.
#'
#' @return A data frame of segments for a specific or multiple sample ID(s).
#'
#' @import dplyr
#' @export
#'
#' @examples
#' #load pacakges
#' suppressPackageStartupMessages(library(dplyr))
#'
#' #get CN segments for one sample
#' dohh2_segs = get_sample_cn_segments(these_sample_ids = "DOHH-2",
#'                                     projection = "hg38",
#'                                     streamlined = TRUE)
#'
#' #get CN segments for DLBCL cell line
#' cell_line_meta = get_gambl_metadata() %>%
#'   dplyr::filter(cohort == "DLBCL_cell_lines")
#'
#' dlbcl_segs = get_sample_cn_segments(these_samples_metadata = cell_line_meta,
#'                                     streamlined = TRUE)
#'
get_sample_cn_segments = function(these_sample_ids = NULL,
                                  these_samples_metadata = NULL,
                                  projection = "grch37",
                                  this_seq_type = "genome",
                                  with_chr_prefix = FALSE,
                                  streamlined = FALSE,
                                  verbose = FALSE,
                                  ...){

  #warn/notify the user what version of this function they are using
  if (!isTRUE(getOption("GAMBLR.open.shown_cn_msg"))) {
    message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")
    options(GAMBLR.open.shown_cn_msg = TRUE)
  }

  #check if any invalid parameters are provided
  check_excess_params(...)

  # resolve samples of interest (id_ease retired)
  if(!is.null(these_samples_metadata)){
    metadata = dplyr::filter(these_samples_metadata, seq_type %in% this_seq_type)
  }else{
    metadata = get_gambl_metadata(seq_type_filter = this_seq_type)
    if(!is.null(these_sample_ids)){
      metadata = dplyr::filter(metadata, sample_id %in% these_sample_ids)
    }
  }

  sample_ids = metadata$sample_id

  #valid projections (kept static so we never load the multi-GB sample_data)
  valid_projections = c("grch37", "hg38")

  #return CN segments based on the selected projection
  if(projection %in% valid_projections){
    all_segs = GAMBLR.data::get_cn_segments_from_db(projection = projection,
                                                    sample_ids = sample_ids)
  }else{
    stop(paste("please provide a valid projection. The following are available:",
               paste(valid_projections,collapse=", ")))
  }

  #deal with chr prefixes
  if(!with_chr_prefix){
    all_segs = all_segs %>%
      dplyr::mutate(chrom = gsub("chr", "", chrom))
  }else{
    if(!grepl("chr", all_segs$chrom[1])){
      all_segs$chrom = paste0("chr", all_segs$chrom)
    }
  }

  if(streamlined){all_segs = dplyr::select(all_segs, ID, CN)}

  return(all_segs)
}
