#' @title Get SSM By Sample.
#'
#' @description Get the SSMs (i.e. load MAF) for a single sample.
#'
#' @details Alias for [get_ssm_by_samples] accepting a single sample ID via
#' `this_sample_id`, provided for API compatibility with GAMBLR.results.
#'
#' @param this_sample_id A single sample ID to retrieve mutations for.
#' @param these_samples_metadata Optional, a metadata table (with sample_id
#' column) to auto-subset the data to samples in that table before returning.
#' @param this_seq_type Default is genome.
#' @param projection The projection genome build. Supports hg38 and grch37.
#' @param ... Any additional parameters passed to [get_ssm_by_samples].
#'
#' @return data frame in MAF format.
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#'
#' maf = get_ssm_by_sample(this_sample_id = "DOHH-2")
#'
get_ssm_by_sample <- function(this_sample_id,
                              these_samples_metadata = NULL,
                              this_seq_type = "genome",
                              projection = "grch37",
                              ...) {
  get_ssm_by_samples(
    these_sample_ids = this_sample_id,
    these_samples_metadata = these_samples_metadata,
    this_seq_type = this_seq_type,
    projection = projection,
    ...
  )
}
