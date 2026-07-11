#' @title Get SSM By Regions.
#'
#' @description Efficiently retrieve all mutations across a range of genomic regions.
#'
#' @details This function internally calls get_ssm_by_region to retrieve SSM calls for the specified regions.
#'
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' @param this_seq_type The this_seq_type you want back, default is genome.
#' @param tool_name Optionally specify which tool to report variant from. The default is slms-3, also supports "publication" to return the exact variants as reported in the original papers.
#' @param regions_list A vector of regions in the chr:start-end format to restrict the returned SSM calls to.
#' @param regions_bed A data frame in BED format with the coordinates you want to retrieve (recommended).
#' This parameter can also accept an additional column with region names that will be added to the return if `use_name_column = TRUE`
#' @param streamlined If set to TRUE (default) only 3 columns will be kept in the returned data frame (start, sample_id and region_name).
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38), default is grch37.
#' @param use_name_column If your bed-format data frame has a name column (must be named "name") these can be used to name your regions.
#' @param verbose Set to TRUE to maximize the output to console. Default is TRUE.
#' This parameter also dictates the verbosity of any helper function internally called inside the main function.
#' @param ... Any additional parameters.
#'
#' @return Returns a data frame of variants in MAF-like format.
#'
#' @import tibble dplyr tidyr
#'
#' @export
#'
#' @examples
#' #basic usage, adding custom names from bundled ashm data frame
#' regions_bed = create_bed_data( GAMBLR.data::grch37_ashm_regions,
#'                           fix_names = "concat",
#'                           concat_cols = c("gene","region"),
#'                           sep="-")
#' 
#' my_meta = get_gambl_metadata()
#' # get a full MAF-format data frame for all aSHM regions on grch37 coordinates
#' ashm_maf = get_ssm_by_regions(regions_bed = regions_bed,
#'                               these_samples_metadata = my_meta,
#'                               streamlined = FALSE)
#' 
#' 
#'
#' one_region_maf = get_ssm_by_regions(regions_list = "2:136875000-136875097",
#'                          streamlined = FALSE,
#'                          projection = "grch37",
#'                          these_samples_metadata = my_meta)
#' \dontrun{
#' # This example fails, as it should
#' #ashm_maf = get_ssm_by_regions(regions_bed = regions_bed,
#' #                              these_samples_metadata = my_meta,
#' #                               projection="hg38")
#' # Error in get_ssm_by_regions(regions_bed = regions_bed, these_samples_metadata = my_meta,  : 
#' # requested projection: hg38 and genome_build of regions_bed: grch37 don't match
#' }
get_ssm_by_regions <- function(these_samples_metadata,
                               regions_list,
                               regions_bed,
                               this_seq_type = "genome",
                               streamlined = TRUE,
                               projection = "grch37",
                               verbose = FALSE,
                               use_name_column = FALSE,
                               tool_name = "slms-3",
                               ...) {

  # check provided projection
  # first, get valid projections
  valid_projections = c("grch37", "hg38")
  if (!projection %in% valid_projections) {
    stop("Please provide a valid projection. The following are available: ",
         paste(valid_projections, collapse = ", "), ".")
  }
  
  # check if any invalid parameters are provided
  check_excess_params(...)

  bed2region = function(x) {
    paste0(x[1], ":", as.numeric(x[2]), "-", as.numeric(x[3]))
  }

  if (missing(regions_list)) {
    if (!missing(regions_bed)) {
      if("bed_data" %in% class(regions_bed)){
        #confirm the genome builds match
        if(is.null(get_genome_build(regions_bed))){
          stop("something is wrong with regions_bed. No genome_build found!")
        }
        if(!get_genome_build(regions_bed)==projection){
          stop(paste("requested projection:",projection,"and genome_build of regions_bed:", get_genome_build(regions_bed), "don't match"))
        }
      }
      regions = apply(regions_bed, 1, bed2region)
    } else {
      if(projection == "grch37"){
        regions_bed = create_bed_data(grch37_ashm_regions,
                      genome_build = projection,
                      fix_names="concat",
                      concat_cols = c("gene","region"),
                      sep="-")
      }else if(projection == "hg38"){
        regions_bed = create_bed_data(grch37_ashm_regions,
                                      genome_build = projection,
                                      fix_names="concat",
                                      concat_cols = c("gene","region"),
                                      sep="-")
      }
      message(paste("defaulting to aSHM regions for ", projection))
      regions = apply(regions_bed, 1, bed2region)
    }
  } else {
    regions = regions_list
  }

  # Warn/notify the user what version of this function they are using
  if (!isTRUE(getOption("GAMBLR.open.shown_ssm_msg"))) {
    message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")
    options(GAMBLR.open.shown_ssm_msg = TRUE)
  }
    if (verbose) {
      print("Using the non-default engine for efficiency...")
    }

    # Build the regions data frame (Chromosome / Start_Position / End_Position / region)
    if(!missing(regions_bed) && "bed_data" %in% class(regions_bed)){
      regions_df = dplyr::select(regions_bed,1:4) %>%
        dplyr::rename(c("Chromosome"="chrom",
                        "Start_Position"="start",
                        "End_Position"="end",
                        "region"="name"))

    }else{
      regions_df <- as.data.frame(regions) %>%
        `names<-`("regions") %>%
        separate(
          regions,
          c("Chromosome", "Start_Position", "End_Position"),
          ":|-"
        ) %>%
        mutate(
          Start_Position = as.numeric(Start_Position),
          End_Position = as.numeric(End_Position),
          region = row_number()
        )
    }

    # Resolve samples of interest (id_ease retired)
    if(!missing(these_samples_metadata) && !is.null(these_samples_metadata)){
      metadata = dplyr::filter(these_samples_metadata, seq_type %in% this_seq_type)
    }else{
      metadata = get_gambl_metadata(seq_type_filter = this_seq_type)
    }

    # Pull only the mutations inside the requested regions from indexed SQL,
    # rather than the genome-wide MAF, then attribute each mutation to its region.
    # Widen by 1 bp so the (strict) SQL range is inclusive; cool_overlaps refines.
    sample_maf <- GAMBLR.data::get_ssm_from_db(
      projection = projection,
      sample_ids = metadata$sample_id,
      tool_name = tool_name,
      include_ashm = TRUE,
      regions = dplyr::transmute(regions_df,
                                 chrom = as.character(Chromosome),
                                 start = as.numeric(Start_Position) - 1,
                                 end   = as.numeric(End_Position) + 1)
    ) %>%
      dplyr::distinct(Tumor_Sample_Barcode, Chromosome,
                      Start_Position, End_Position, .keep_all = TRUE) %>%
      create_maf_data(projection) %>%
      mutate(maf_seq_type = this_seq_type)

    region_mafs <- cool_overlaps(
      sample_maf,
      regions_df
    ) %>%
      dplyr::rename_with(~ gsub(".x", "", .x, fixed = TRUE)) %>%
      dplyr::select(all_of(c(names(sample_maf), "region"))) %>%
      dplyr::group_split(region)
    maf_df = do.call(bind_rows, region_mafs)

    if(!use_name_column){
      maf_df = mutate(maf_df,region=paste0(Chromosome,":",Start_Position))
    }
    maf_df = dplyr::rename(maf_df,c("region_name"="region"))
    if(streamlined){
      
      maf_df = dplyr::select(maf_df,Start_Position,Tumor_Sample_Barcode,region_name) %>%
        dplyr::rename(c("sample_id"="Tumor_Sample_Barcode","start"="Start_Position"))
      
    }
    return(maf_df)
    

}