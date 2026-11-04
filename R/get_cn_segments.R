#' @title Get CN Segments.
#'
#' @description Retrieve all copy number segments from the GAMBL outputs
#'
#' @details This function merely loads and returns all the seg_data
#' available for a projection (genome build)
#' @param these_samples_metadata User must provide a metadata table to
#' restrict the data to the samples in your table.
#' The metadata also ensures the proper handling of duplicate sample_id
#' across seq_types and ensures the seq_type in the metadata faithfully
#' represents the seq_type of the data
#' @param projection Desired genome coordinate system for returned CN segments.
#' Default is "grch37".
#' @param this_seq_type Deprecated.
#' @param overlap_action What to do if a sample has overlapping segments,
#' which should never happen and would cause double counting downstream:
#' "warn" (default, naming the samples), "error", "resolve" (clean up with
#' [GAMBLR.utils::resolve_segment_overlaps]) or "none".
#' @param pipeline Optional: segments from published analyses instead of
#' the bundled GAMBL segments, e.g. "dlbclass" (DLBClass/Chapuy et al.) or
#' "schmitz_array" (Schmitz et al. arrays); see
#' [GAMBLR.data::published_cn_sources]. These cover capture samples of the
#' dlbcl_chapuy and dlbcl_schmitz cohorts, under GAMBL sample_ids. "all"
#' returns the bundled GAMBL segments (as pipeline "gambl") and every
#' published source side by side, e.g. for
#' GAMBLR.viz::plot_sample_cn_comparison. Default NULL: the bundled GAMBL
#' segments only.
#' @param source_rule How each sample's GAMBL segmentation is chosen, as in
#' GAMBLR.results::get_cn_segments: "cnv_master" (default), "qc_score",
#' "preference", or a
#' saved selection (a table, or the path of a TSV of one, with sample_id and
#' pipeline). GAMBLR.data bundles one GAMBL segmentation per sample (the
#' combined merge's), so every rule returns it: "qc_score" and "preference"
#' say so, and with a
#' saved selection, samples it lacks are dropped and those it assigns to
#' another pipeline get the bundled segments, with a message. Not used with
#' `pipeline`.
#' @param ... Additional parameters to be passed to the function.
#'
#' @return A data frame with CN segments for the specified region, with a
#' `pipeline` column naming each segment's source ("gambl" for the bundled
#' segments; published sources add source_id, purity and ploidy columns).
#'
#' @import dplyr
#' @export
#'
#' @examples
#' # Example for the capture samples:
#'
#' genome_metadata = get_gambl_metadata(seq_type_filter="genome")
#'
#' genome_segments_hg38 = get_cn_segments(
#'                              these_samples_metadata = genome_metadata,
#'                              projection="hg38")
#'
#' # the DLBClass (Chapuy et al.) segments of the Chapuy samples
#' chapuy = get_gambl_metadata(seq_type_filter = "capture") %>%
#'   dplyr::filter(cohort == "dlbcl_chapuy")
#' dlbclass_segments = get_cn_segments(chapuy, pipeline = "dlbclass")
#'
get_cn_segments = function(these_samples_metadata,
                           projection = "grch37",
                           pipeline = NULL,
                           source_rule = "cnv_master",
                           this_seq_type,
                           overlap_action = "warn", ...) {
  #warn/notify the user what version of this function they are using
  if (!isTRUE(getOption("GAMBLR.open.shown_cn_msg"))) {
    message("Using the bundled CN segments (.seg) calls in GAMBLR.data...")
    options(GAMBLR.open.shown_cn_msg = TRUE)
  }

  #check if any invalid parameters are provided
  check_excess_params(...)

  #valid projections (kept static so we never load the multi-GB sample_data)
  valid_projections = c("grch37", "hg38")

  metadata = these_samples_metadata

  #source_rule: GAMBLR.data has one GAMBL segmentation per sample, so only a
  #saved selection changes anything (the samples kept)
  if (!identical(source_rule, "cnv_master")) {
    if (!is.null(pipeline)) stop("give pipeline or source_rule, not both")
    if (is.character(source_rule) && length(source_rule) == 1 && source_rule %in% c("qc_score", "preference")) {
      message("GAMBLR.data bundles one GAMBL segmentation per sample (the combined merge's), ",
              "so source_rule = \"", source_rule, "\" returns it; choosing among pipelines needs GAMBLR.results")
    } else {
      sel = if (is.character(source_rule)) {
        if (!file.exists(source_rule)) stop("source_rule must be \"cnv_master\", \"qc_score\", \"preference\", ",
                                            "a selection table or the path of a saved one; not found: ", source_rule)
        utils::read.delim(source_rule, comment.char = "#", stringsAsFactors = FALSE)
      } else as.data.frame(source_rule)
      if (!"sample_id" %in% colnames(sel) && "ID" %in% colnames(sel)) sel$sample_id = sel$ID
      if (!all(c("sample_id", "pipeline") %in% colnames(sel))) {
        stop("source_rule selection needs sample_id and pipeline columns")
      }
      if ("selected" %in% colnames(sel)) sel = sel[sel$selected %in% TRUE, , drop = FALSE]
      lacking = !metadata$sample_id %in% sel$sample_id
      if (any(lacking)) {
        message("dropping ", length(unique(metadata$sample_id[lacking])),
                " sample(s) with no pipeline chosen by source_rule")
        metadata = metadata[!lacking, , drop = FALSE]
      }
      n_other = sum(sel$sample_id %in% metadata$sample_id & !sel$pipeline %in% "gambl")
      if (n_other > 0) {
        message(n_other, " sample(s) get the bundled GAMBL segments in place of the pipeline chosen in ",
                "source_rule (GAMBLR.data bundles one GAMBL segmentation per sample)")
      }
    }
  }

  sample_ids = metadata$sample_id
  if (!projection %in% valid_projections) {
    stop(paste("please provide a valid projection.",
               paste(valid_projections, collapse = ", ")))
  }
  #segments from published analyses (pipeline), and/or the bundled GAMBL ones
  published = NULL
  with_gambl = is.null(pipeline) || identical(pipeline, "all")
  if (!is.null(pipeline)) {
    if (!"get_published_cn_segments" %in% getNamespaceExports("GAMBLR.data")) {
      stop("pipeline needs a GAMBLR.data with get_published_cn_segments(); please update GAMBLR.data")
    }
    sources = GAMBLR.data::published_cn_sources()$source
    wanted = if (identical(pipeline, "all")) sources else pipeline
    unknown = setdiff(wanted, sources)
    if (length(unknown) > 0) {
      stop("unknown pipeline(s): ", paste(unknown, collapse = ", "),
           "; available: ", paste(c(sources, "all"), collapse = ", "))
    }
    published = GAMBLR.data::get_published_cn_segments(wanted, projection = projection,
                                                       sample_ids = sample_ids)
    published = dplyr::rename(published, pipeline = source)
    missing_ids = setdiff(sample_ids, published$ID)
    if (length(missing_ids) > 0) {
      message(length(missing_ids), " sample(s) have no ", paste(wanted, collapse = "/"),
              " segments (published sources cover capture samples of dlbcl_chapuy and dlbcl_schmitz)")
    }
  }
  if (with_gambl) {
    all_segs = GAMBLR.data::get_cn_segments_from_db(projection = projection,
                                                    sample_ids = sample_ids)
    all_segs$pipeline = rep("gambl", nrow(all_segs))
    if (!is.null(published)) {
      all_segs = dplyr::bind_rows(all_segs, published)
    }
  } else {
    all_segs = published
  }

  #ensure chr prefixes are there when necessary
  if(projection == "grch37") {
    if(grepl("chr", all_segs$chrom[1])) {
      all_segs = all_segs %>%
        dplyr::mutate(chrom = gsub("chr", "", chrom))
    }
  }else {
    if (!grepl("chr",all_segs$chrom[1])) {
      all_segs = all_segs %>%
        dplyr::mutate(chrom = paste0("chr", chrom))
    }
  }

  #segments of a sample (profile) should never overlap (per pipeline, when several)
  overlap_by = if ("pipeline" %in% colnames(all_segs)) c("ID", "pipeline") else NULL
  all_segs = GAMBLR.utils::handle_segment_overlaps(all_segs, action = overlap_action,
                                                   by = overlap_by)
  #return S3 class with CN segments and genome_build
  all_segs = create_seg_data(all_segs, projection)
  return(all_segs)
}
