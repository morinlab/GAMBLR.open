#' @title Get SSM By Regions.
#'
#' @description Efficiently retrieve all mutations across a range of genomic regions.
#'
#' @details This function internally calls get_ssm_by_region to retrieve SSM calls for the specified regions.
#'
#' @param these_samples_metadata Optional, a metadata table (with sample IDs in a column) to subset the return to.
#' @param regions_list A vector of regions in the chr:start-end format to restrict the returned SSM calls to.
#' @param regions_bed A data frame in BED format with the coordinates you want to retrieve (recommended).
#' This parameter can also accept an additional column with region names that will be added to the return if `use_name_column = TRUE`
#' @param maf_data Use an already loaded MAF object of class "genomic_data" or "maf_data".
#' @param this_seq_type The this_seq_type you want back, default is genome.
#' @param streamlined If set to TRUE (default) only 3 columns will be kept in the returned data frame (start, sample_id and region_name).
#' @param projection Obtain variants projected to this reference (one of grch37 or hg38), default is grch37.
#' @param use_name_column If your bed-format data frame has a name column (must be named "name") these can be used to name your regions.
#' @param tool_name Optionally specify which tool to report variant from. The default is slms-3, also supports "publication" to 
#'  return the exact variants as reported in the original papers.
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
                               maf_data,
                               this_seq_type = "genome",
                               streamlined = TRUE,
                               projection = "grch37",
                               use_name_column = FALSE,
                               tool_name = "slms-3",
                               verbose = FALSE,
                               ...) {

  # check provided projection
  # first, get valid projections
  valid_projections = grep("meta", names(GAMBLR.data::sample_data),
                           value = TRUE, invert = TRUE)
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
    if(projection == "grch37"){
      regions = gsub(",", "", gsub("chr", "", regions_list))
    }else{
      regions = gsub(",", "", regions_list)
    }
  }

  if(!missing(maf_data)){
    # Warn/notify the user what version of this function they are using
    message("Using the supplied maf_data")

    # Confirm the genome build and projection match
    if(is.null(get_genome_build(maf_data))){
      stop("No genome_build found for maf_data! Is it a genomic_data or maf_data Object?")
    }
    if(!get_genome_build(maf_data)==projection){
      stop(paste("requested projection:",projection,"and genome_build of maf_data:", 
            get_genome_build(maf_data), "don't match"))
    }
    sample_maf <- maf_data %>%
      dplyr::filter(Tumor_Sample_Barcode %in% these_samples_metadata$sample_id)
    # Remove any existing region_name or name columns from this maf being made by
    # a previous call to get_ssm_by_regions. Only annotations from the current
    # regions_bed / regions_list will remain
    if("region_name" %in% names(sample_maf)){
      sample_maf <- sample_maf %>%
        dplyr::select(-region_name)
    }
    if("name" %in% names(sample_maf)){
      sample_maf <- sample_maf %>%
        dplyr::select(-name)
    }
  } else{
    # Warn/notify the user what version of this function they are using
    message("Using the bundled SSM calls (.maf) calls in GAMBLR.data...")
    if (verbose) {
      print("Using the non-default engine for efficiency...")
    }
    sample_maf <- get_ssm_by_samples(
      these_samples_metadata = these_samples_metadata,
      this_seq_type = this_seq_type,
      projection = projection,
      tool_name = tool_name
    )
  }

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
    # To remove duplicate rows in the case of the same region provided more than once
    maf_df <- unique(maf_df)

    if(streamlined){
      
      maf_df = dplyr::select(maf_df,Start_Position,Tumor_Sample_Barcode,region_name) %>%
        dplyr::rename(c("sample_id"="Tumor_Sample_Barcode","start"="Start_Position"))
      
    }
    return(maf_df)
    

}