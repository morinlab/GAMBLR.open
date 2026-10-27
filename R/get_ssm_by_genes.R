#' @title Get all available SSMs for an arbitrary set of genes
#'
#' @description Retrieve all SSMs for one or more genes using a
#' coordinate-based lookup, including non-coding variant classes (Intron,
#' 5'UTR, 3'UTR, flanks, etc.), unlike [get_coding_ssm] and
#' [get_all_coding_ssm].
#'
#' @details Each gene's coordinates come from
#' [GAMBLR.utils::gene_to_region], padded by `expand_by`, and its mutations
#' from [get_ssm_by_region] for each seq_type, keeping the rows whose
#' `Hugo_Symbol` is the gene or one of its aliases
#' ([GAMBLR.utils::expand_gene_aliases]; e.g. old and new histone gene
#' names). To keep only the mutations near each gene's transcription start
#' site, pass the result to [GAMBLR.utils::subset_ssm_by_tss].
#'
#' @param genes Gene symbol for one or more genes.
#' @param these_samples_metadata Optional metadata table (with sample_id and
#' seq_type columns) to restrict the samples. Default: every sample in
#' [get_gambl_metadata] of the seq_types in `this_seq_type`.
#' @param projection Obtain variants projected to this reference (one of
#' grch37 or hg38). Default grch37.
#' @param this_seq_type The seq_types to return mutations for. Default
#' c("genome", "capture").
#' @param expand_by Padding (bp) around each gene's coordinates for the
#' lookup, before keeping the gene's own mutations. Default 10000, so that
#' mutations attributed to the gene just outside its annotated span (e.g.
#' upstream, within VEP's 5 kb annotation window) aren't missed.
#' @param tool_name Variant caller to report, passed to [get_ssm_by_region].
#' Default "slms-3".
#' @param verbose Print a message for each gene. Default FALSE.
#' @param ... Any additional parameters.
#'
#' @return A maf_data object with all MAF columns, one row per mutation,
#' plus `maf_seq_type`.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' meta = get_gambl_metadata() %>%
#'   dplyr::filter(pathology == "DLBCL", seq_type == "genome") %>%
#'   head(40)
#' genes_maf = get_ssm_by_genes(c("BCL2", "PIM1", "EZH2"),
#'                              these_samples_metadata = meta)
#' dplyr::count(genes_maf, Hugo_Symbol, Variant_Classification)
#'
#' # non-coding mutations in the first 4 kb of each gene (LymphGen's Synon
#' # window)
#' synon = GAMBLR.utils::subset_ssm_by_tss(genes_maf, genes = c("BCL2", "PIM1"))
#' dplyr::count(synon, tss_gene)
get_ssm_by_genes = function(genes,
                            these_samples_metadata = NULL,
                            projection = "grch37",
                            this_seq_type = c("genome", "capture"),
                            expand_by = 10000,
                            tool_name = "slms-3",
                            verbose = FALSE,
                            ...){
  check_excess_params(...)
  if(!projection %in% c("grch37", "hg38")){
    stop("projection must be grch37 or hg38")
  }
  if(is.null(these_samples_metadata)){
    these_samples_metadata = suppressMessages(get_gambl_metadata())
  }
  these_samples_metadata = dplyr::filter(these_samples_metadata, seq_type %in% this_seq_type)

  all_ssms = list()
  for(gene in unique(genes)){
    gene_region = suppressMessages(GAMBLR.utils::gene_to_region(gene, projection = projection,
                                                                pad_length = expand_by))
    if(is.null(gene_region) || length(gene_region) == 0){
      warning("No coordinates found for gene: ", gene, ", skipping")
      next
    }
    if(verbose){
      message("Processing gene: ", gene, " with region: ", gene_region[1])
    }
    for(s_type in unique(these_samples_metadata$seq_type)){
      # gene_to_region() resolves the gene even when its MAF Hugo_Symbol is
      # another name (e.g. old/new histone gene names), so match all aliases
      ssms = suppressMessages(get_ssm_by_region(region = gene_region[1],
                                                these_samples_metadata = these_samples_metadata,
                                                this_seq_type = s_type,
                                                projection = projection,
                                                tool_name = tool_name,
                                                streamlined = FALSE))
      ssms = dplyr::filter(as.data.frame(ssms),
                           Hugo_Symbol %in% GAMBLR.utils::expand_gene_aliases(gene))
      if(nrow(ssms) > 0){
        all_ssms[[paste(gene, s_type, sep = "-")]] = dplyr::mutate(ssms, maf_seq_type = s_type)
      }
    }
  }
  all_ssm_maf = dplyr::bind_rows(all_ssms)
  if(nrow(all_ssm_maf) == 0){
    message("no mutations found for these genes and samples")
  }
  GAMBLR.utils::create_maf_data(all_ssm_maf, genome_build = projection)
}
