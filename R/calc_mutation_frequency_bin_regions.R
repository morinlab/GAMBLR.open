#' @title Mutation counts across sliding windows for multiple regions.
#'
#' @description Obtain a long tidy or wide matrix of mutation counts across
#' sliding windows for multiple regions.
#'
#' @details This function takes a metadata table with `these_samples_metadata` 
#' parameter and internally calls `calc_mutation_frequency_bin_region` 
#' (that internally calls `get_ssm_by_regions`).
#' to retrieve mutation counts for sliding windows across one or more regions. 
#' May optionally provide any combination of a maf data frame, existing metadata,
#' or a regions data frame or named vector.
#'
#' @param regions_list Named vector of regions in the format 
#' c(name1 = "chr:start-end", name2 = "chr:start-end"). If neither `regions` nor 
#' `regions_bed` is specified, the function will use GAMBLR aSHM region information.
#' @param regions_bed Data frame of regions with four columns (chrom, start, end, name).
#' @param these_samples_metadata Metadata with at least sample_id column. 
#' If not providing a maf data frame, seq_type is also required.
#' @param these_sample_ids Vector of sample IDs. Metadata will be subset to
#' sample IDs present in this vector.
#' @param this_seq_type Optional vector of seq_types to include in heatmap. 
#' Default "genome". Uses default seq_type priority for samples with >1 seq_type.
#' @param maf_data Optional maf data frame. Will be subset to rows where 
#' Tumor_Sample_Barcode matches provided sample IDs or metadata table. 
#' If not provided, maf data will be obtained with get_ssm_by_regions().
#' @param region_padding Amount to pad the start and end coordinates by. Default 1000.
#' @param projection Genome build the function will operate in. Ensure this 
#' matches your provided regions and maf data for correct chr prefix handling. Default "grch37".
#' @param drop_unmutated Whether to drop bins with 0 mutations. If returning a 
#' matrix format, this will only drop bins with no mutations in any samples.
#' @param skip_regions Optional character vector of genes to exclude from the default aSHM regions.
#' @param only_regions Optional character vector of genes to include from the default aSHM regions.
#' @param slide_by Slide size for sliding window. Default 100.
#' @param window_size Size of sliding window. Default 500.
#' @param return_format Return format of mutations. Accepted inputs are "long" and 
#' "wide". Long returns a data frame of one sample ID/window per row. Wide returns 
#' a matrix with one sample ID per row and one window per column. Using the "wide" 
#' format will retain all samples and windows regardless of the drop_unmutated or 
#' min_count_per_bin parameters. Default wide.
#' @param ... Any additional parameters.
#' 
#' @return A table of mutation counts for sliding windows across one or more regions. May be long or wide.
#'
#' @import dplyr tidyr tibble
#' @export
#'
#' @examples

#'  #load metadata.
#'  my_meta = get_gambl_metadata()
#'  dlbcl_bl_meta = dplyr::filter(my_meta, pathology %in% c("DLBCL", "BL"))
#'
#'
#'  #get ashm regions
#'  some_regions = create_bed_data(grch37_ashm_regions,
#'                                 fix_names = "concat",
#'                                 concat_cols = c("gene","region"),
#'                                 sep="-")
#'  print(some_regions)
#'  mut_count_matrix <- calc_mutation_frequency_bin_regions(
#'    these_samples_metadata = dlbcl_bl_meta,
#'    regions_bed = some_regions
#'  )
#' dim(mut_count_matrix)
#' tail(mut_count_matrix[,c(1:10)])
calc_mutation_frequency_bin_regions <- function(regions_list = NULL,
                                                regions_bed = NULL,
                                                these_samples_metadata = NULL,
                                                these_sample_ids = NULL,
                                                this_seq_type = "genome",
                                                maf_data = NULL,
                                                projection = "grch37",
                                                region_padding = 1000,
                                                drop_unmutated = FALSE,
                                                skip_regions = NULL,
                                                only_regions = NULL,
                                                slide_by = 100,
                                                window_size = 500,
                                                return_format = "wide",
                                                ...){

  #check if any invalid parameters are provided
  check_excess_params(...)
  
  regions <- process_regions(regions_list = regions_list,
                             regions_bed = regions_bed,
                             region_padding = region_padding,
                             skip_regions = skip_regions,
                             only_regions = only_regions)
  
  regions_bed <- regions$regions_bed
  regions <- regions$regions_list
  
  if (
    (grepl("chr", regions_bed$chrom[1]) & projection == "grch37") |
    (!grepl("chr", regions_bed$chrom[1]) & projection == "hg38")
  ) {
    stop("chr prefixing status of provided regions and specified projection don't match. ")
  }
  # Harmonize metadata and sample IDs (id_ease retired)
  if(!is.null(these_samples_metadata)){
    metadata <- dplyr::filter(these_samples_metadata, seq_type %in% this_seq_type)
  }else{
    metadata <- get_gambl_metadata(seq_type_filter = this_seq_type)
    if(!is.null(these_sample_ids)){
      metadata <- dplyr::filter(metadata, sample_id %in% these_sample_ids)
    }
  }
  
  these_sample_ids <- metadata$sample_id

  # Pre-fetch all mutations across every requested (padded) region in a
  # single query, instead of letting the per-region loop below call
  # calc_mutation_frequency_bin_region() -> get_ssm_by_region() once per
  # region -- each a separate database round trip. For the full aSHM region
  # set (100+ regions) this was by far the largest cost in the GAMBLR.open
  # example suite (~140s of a ~240s total run). get_ssm_by_regions()
  # (plural) already consolidates a whole regions_bed into one query;
  # passing its result down as maf_data makes every per-region call below
  # take calc_mutation_frequency_bin_region()'s existing in-memory
  # subsetting path (cool_overlaps(), which already matches on inclusive
  # start/end boundaries -- same semantics either way, so this doesn't
  # change what gets counted) instead of querying the database again.
  #
  # process_regions() deliberately leaves regions_bed itself un-padded --
  # only `regions` (the per-region loop input, built from regions_bed a few
  # lines above) has region_padding applied -- so it has to be re-applied
  # here too, or this pre-fetch would miss mutations that fall inside the
  # padding but outside the raw region.
  if (is.null(maf_data)) {
    padded_regions_bed <- regions_bed %>%
      dplyr::mutate(start = start - region_padding, end = end + region_padding)
    maf_data <- get_ssm_by_regions(
      regions_bed = padded_regions_bed,
      these_samples_metadata = metadata,
      this_seq_type = unique(metadata$seq_type),
      streamlined = FALSE,
      projection = projection
    )
  }

  # Obtain sliding window mutation frequencies for all regions. Sequential
  # (not mclapply()) now that maf_data is pre-fetched: each iteration below
  # is just fast in-memory subsetting, not a database call, so there's
  # little left to parallelize -- and mclapply()'s fork-per-worker model
  # both multiplies peak memory use (each fork copies the parent's memory)
  # and is a poor fit for sharing one open DB connection across processes.
  dfs <- lapply(names(regions), function(x) {
    df <- calc_mutation_frequency_bin_region(
      region = regions[x],
      these_samples_metadata = metadata,
      maf_data = maf_data,
      projection = projection,
      drop_unmutated = drop_unmutated,
      slide_by = slide_by,
      window_size = window_size,
      min_count_per_bin = 0,
      return_count = TRUE,
      ...
    ) %>%
      dplyr::mutate(name = x)
    return(df)
  })
  
  all <- dplyr::bind_rows(dfs) %>%
    dplyr::distinct(bin, sample_id, .keep_all = TRUE)
  
  # If none of the samples are mutated, return the mutation frequency df and exit.
  if (max(all$mutation_count) == 0) {
    message("No mutations found in specified regions for specified samples. Exiting. ")
    return(all)
  }
  
  if (return_format == "wide") {
    # Convert mutation frequency table to a matrix
    all_wide <- all %>%
      dplyr::select(sample_id, mutation_count, bin) %>%
      pivot_wider(
        names_from = bin,
        values_from = mutation_count,
        values_fill = 0
      )
    return(all_wide)
  } else {
    return(all)
  }
}
