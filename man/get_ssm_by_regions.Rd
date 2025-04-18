% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_ssm_by_regions.R
\name{get_ssm_by_regions}
\alias{get_ssm_by_regions}
\title{Get SSM By Regions.}
\usage{
get_ssm_by_regions(
  these_samples_metadata,
  regions_list,
  regions_bed,
  this_seq_type = "genome",
  streamlined = TRUE,
  projection = "grch37",
  verbose = FALSE,
  use_name_column = FALSE,
  tool_name = "slms-3",
  ...
)
}
\arguments{
\item{these_samples_metadata}{Optional, a metadata table (with sample IDs in a column) to subset the return to.}

\item{regions_list}{A vector of regions in the chr:start-end format to restrict the returned SSM calls to.}

\item{regions_bed}{A data frame in BED format with the coordinates you want to retrieve (recommended).
This parameter can also accept an additional column with region names that will be added to the return if \code{use_name_column = TRUE}}

\item{this_seq_type}{The this_seq_type you want back, default is genome.}

\item{streamlined}{If set to TRUE (default) only 3 columns will be kept in the returned data frame (start, sample_id and region_name).}

\item{projection}{Obtain variants projected to this reference (one of grch37 or hg38), default is grch37.}

\item{verbose}{Set to TRUE to maximize the output to console. Default is TRUE.
This parameter also dictates the verbosity of any helper function internally called inside the main function.}

\item{use_name_column}{If your bed-format data frame has a name column (must be named "name") these can be used to name your regions.}

\item{tool_name}{Optionally specify which tool to report variant from. The default is slms-3, also supports "publication" to return the exact variants as reported in the original papers.}

\item{...}{Any additional parameters.}
}
\value{
Returns a data frame of variants in MAF-like format.
}
\description{
Efficiently retrieve all mutations across a range of genomic regions.
}
\details{
This function internally calls get_ssm_by_region to retrieve SSM calls for the specified regions.
}
\examples{
#basic usage, adding custom names from bundled ashm data frame
regions_bed = create_bed_data( GAMBLR.data::grch37_ashm_regions,
                          fix_names = "concat",
                          concat_cols = c("gene","region"),
                          sep="-")

my_meta = get_gambl_metadata()
# get a full MAF-format data frame for all aSHM regions on grch37 coordinates
ashm_maf = get_ssm_by_regions(regions_bed = regions_bed,
                              these_samples_metadata = my_meta,
                              streamlined = FALSE)



one_region_maf = get_ssm_by_regions(regions_list = "2:136875000-136875097",
                         streamlined = FALSE,
                         projection = "grch37",
                         these_samples_metadata = my_meta)
\dontrun{
# This example fails, as it should
#ashm_maf = get_ssm_by_regions(regions_bed = regions_bed,
#                              these_samples_metadata = my_meta,
#                               projection="hg38")
# Error in get_ssm_by_regions(regions_bed = regions_bed, these_samples_metadata = my_meta,  : 
# requested projection: hg38 and genome_build of regions_bed: grch37 don't match
}
}
