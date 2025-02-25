#' @title Get Manta SVs
#'
#' @description Retrieve Manta SVs for one or many samples
#' 
#' @details Retrieve Manta SVs with additional VCF information to allow for
#' filtering of high-confidence variants.
#' To get SV calls for multiple samples, supply a metadata table via
#' `these_samples_metadata` that has been subset to only those samples.
#' The results will be restricted to the sample_ids within that data frame.
#' This function can also restrict the returned breakpoints within a genomic
#' region specified via `region` (in chr:start-end format).
#' Useful filtering parameters are also available, use `min_vaf` to set the
#' minimum tumour VAF for a SV to be returned and `min_score`
#' to set the lowest Manta somatic score for a SV to be returned.
#' In addition, the user can chose to return all variants, even
#' the ones not passing the filter criteria. To do so,
#' set `pass_filters = FALSE` (defaults to TRUE).
#' 
#' @param these_samples_metadata A metadata data frame to limit the
#' result to sample_ids within it
#' @param projection The projection genome build. Default is grch37.
#' @param min_vaf The minimum tumour VAF for a SV to be returned.
#' Default is 0.1.
#' @param min_score The lowest Manta somatic score for a SV to be returned.
#' Default is 40.
#' @param pass_filters If TRUE (default) only return SVs that are annotated with
#' PASS in the FILTER column. Set to FALSE to keep all variants,
#' regardless if they PASS the filters.
#' @param verbose Set to FALSE to minimize the output to console.
#' Default is TRUE. This parameter also dictates the verbose-ness of
#' any helper function internally called inside the main function.
#' @param region Specify a single region to fetch SVs anchored within
#' using the format "chrom:start-end"
#' @param chromosome DEPRECATED. Use `region` instead.
#' @param qstart DEPRECATED. Use `region` instead.
#' @param qend DEPRECATED. Use `region` instead.
#' @param pairing_status DEPRECATED.
#' @param these_sample_ids DEPRECATED. 
#' Subset your metadata and supply `these_samples_metadata`` instead.
#' @param ... Any additional parameters.
#' 
#' @export
#' @import dplyr
#' @examples
#' # lazily get every SV in the table with default quality filters
#' all_sv <- get_manta_sv()
#' head(all_sv)
#' 
#' # get all SVs for just one cohort
#' cohort_meta = suppressMessages(get_gambl_metadata()) %>% 
#'               dplyr::filter(cohort == "DLBCL_cell_lines")
#'
#' some_sv <- get_manta_sv(these_samples_metadata = cohort_meta, verbose=FALSE)
#' head(some_sv)
#' nrow(some_sv)
#' 
#' # get the SVs in a region around MYC
#' # WARNING: This is not the best way to find MYC SVs.
#' # Use annotate_sv on the full SV set instead.
#' myc_region_hg38 = "chr8:127710883-127761821"
#' myc_region_grch37 = "8:128723128-128774067"
#' 
#' hg38_myc_locus_sv <- get_manta_sv(region = myc_region_hg38,
#'                                 projection = "hg38",
#'                                 verbose = FALSE)
#' head(hg38_myc_locus_sv)
#' nrow(hg38_myc_locus_sv)
#' 
#' incorrect_myc_locus_sv <- get_manta_sv(region = myc_region_grch37,
#'                                 projection = "hg38",
#'                                 verbose = FALSE)
#' head(incorrect_myc_locus_sv)
#' nrow(incorrect_myc_locus_sv)
#' # The effect of specifying the wrong coordinate is evident
#' 
#' # Despite potentially being incomplete, we can nonetheless
#' # annotate these directly for more details
#' annotated_myc_hg38 = suppressMessages(
#'          annotate_sv(hg38_myc_locus_sv, genome_build = "hg38")
#' )
#' head(annotated_myc_hg38)
#' table(annotated_myc_hg38$partner)
#' # The usual MYC partners are seen here
#'
get_manta_sv = function(these_samples_metadata = NULL,
                        projection = "grch37",
                        region,
                        min_vaf = 0.1,
                        min_score = 40,
                        pass_filters = TRUE,
                        verbose = FALSE,
                        chromosome,
                        qstart,
                        qend,
                        pairing_status,
                        these_sample_ids = NULL,
                        ...){
  if (!missing(these_sample_ids) | !missing(chromosome) | !missing(qstart) | !missing(qend)) {
    stop("Parameters `these_sample_ids`, `chromosome`, `qstart` and `qend` deprecated and will be ignored.
    Please use `these_samples_metadata` and/or `region` instead.")
  }
  if(missing(these_samples_metadata)){
    these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(seq_type=="genome")
  }
  #warn/notify the user what version of this function they are using
  message("Using the bundled Manta SV (.bedpe) calls in GAMBLR.data...")
  
  #check if any invalid parameters are provided
  check_excess_params(...)
  
  #get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data), value = TRUE, invert = TRUE)
  
  metadata = these_samples_metadata
  
  sample_ids = metadata$sample_id
  
  #return manta SV based on the selected projection
  if(projection %in% valid_projections){
    manta_sv = GAMBLR.data::sample_data[[projection]]$bedpe %>% 
      dplyr::filter(tumour_sample_id %in% sample_ids)
  }else{
    stop(paste("please provide a valid projection.
    The following are available:",
               paste(valid_projections,collapse=", ")))
  }
  
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = startend[1]
    qend = startend[2]
  }
  
  manta_sv = manta_sv %>%
    dplyr::filter(VAF_tumour >= min_vaf,
                  SCORE >= min_score)
  
  if(verbose){
    no_manta = setdiff(metadata$sample_id, manta_sv$tumour_sample_id)
    
    if(length(no_manta) > 0){
      message(paste0("No Manta results found for ", length(no_manta), " samples..."))
      print(no_manta)
    }
  }
  
  #deal with chr prefixes based on the selected projection (if return is to be subset to regions...)
  if(!missing(region) || !missing(chromosome)){
    if(projection == "grch37"){
      if(grepl("chr", chromosome)){
        chromosome = gsub("chr", "", chromosome)
      }
    }else if(projection == "hg38"){
      if(!grepl("chr", chromosome)){
        chromosome = paste0("chr", chromosome)
      }
    }
    
    manta_sv = manta_sv %>%
      dplyr::filter((CHROM_A == chromosome & START_A >= qstart & START_A <= qend) | (CHROM_B == chromosome & START_B >= qstart & START_B <= qend))
  }
  
  if(verbose){
    message("\nThe following VCF filters are applied;")
    message(paste0("  Minimum VAF: ", min_vaf))
    message(paste0("  Minimum Score: ", min_score))
    message(paste0("  Only keep variants passing the quality filter: ", pass))
  }
  
  #PASS filter
  if(pass_filters){
    manta_sv = manta_sv %>%
      dplyr::filter(FILTER == "PASS")
  }
  
  #attach genome_build 
  manta_sv = create_genomic_data(manta_sv, projection)
  
  if(verbose){
    n_variants = nrow(manta_sv)
    unique_samples = unique(manta_sv$tumour_sample_id)
    message(paste0("Returning ", n_variants,
                   " variants from ", length(unique_samples), " sample(s)"))
    message("\nDone!")
  }
  return(manta_sv)
}
