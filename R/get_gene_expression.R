#' @title Get Gene Expression.
#'
#' @description Get the expression of one or more genes for the GAMBL RNA-seq
#' samples linked to the bundled DNA samples.
#'
#' @details Returns DESeq2 variance-stabilised expression values. Values are
#' fetched from the GAMBLR.data expression source one gene at a time the first
#' time a gene is requested, then served from a local cache (see
#' [GAMBLR.data::get_expression_values()] for where the data come from, and
#' [GAMBLR.data::purge_expression_cache()] to clear the cache). The first
#' request for a gene is slower; later ones are fast.
#'
#' The function has the same name and main arguments as
#' `GAMBLR.results::get_gene_expression()`, which reads GAMBL's full
#' expression results for GSC users.
#'
#' RNA and DNA sample IDs often differ for the same biopsy. Use
#' `lazy_join = TRUE` to get `genome_sample_id` and `capture_sample_id`, the
#' DNA samples of the same biopsy, and join to DNA results on those columns.
#' When `these_samples_metadata` is given, an RNA sample is kept if its
#' `mrna_sample_id`, `genome_sample_id` or `capture_sample_id` is one of the
#' metadata's `sample_id`s, so DNA metadata from [get_gambl_metadata()] can be
#' passed directly.
#'
#' @param these_samples_metadata Optional data frame of sample metadata (e.g.
#'   from [get_gambl_metadata()]). Only RNA samples linked to its sample_ids
#'   are returned.
#' @param hugo_symbols One or more gene symbols. Cannot be used together with
#'   ensembl_gene_ids.
#' @param ensembl_gene_ids One or more Ensembl gene IDs. Cannot be used
#'   together with hugo_symbols.
#' @param all_genes Not available here (it would download all expression
#'   data); request the genes you need.
#' @param verbose Set to TRUE for more messages.
#' @param engine Accepted for compatibility with
#'   `GAMBLR.results::get_gene_expression()`; ignored.
#' @param format Either "wide" (default; one column per gene) or "long" (one
#'   row per sample and gene, with ensembl_gene_id, Hugo_Symbol and
#'   expression columns).
#' @param lazy_join If TRUE, add the genome_sample_id and capture_sample_id
#'   of the biopsy each RNA sample came from.
#' @param arbitrarily_pick A symbol can belong to more than one Ensembl gene.
#'   By default that is an error naming the genes, so you can request the one
#'   you want with ensembl_gene_ids; set this to TRUE to use the first of them
#'   (in genome order) instead.
#' @param HGNC Accepted for compatibility with
#'   `GAMBLR.results::get_gene_expression()`; ignored.
#' @param ... Any additional parameters.
#'
#' @return A data frame with mrna_sample_id, patient_id, protocol,
#'   ffpe_or_frozen and study (plus genome_sample_id and capture_sample_id
#'   with `lazy_join = TRUE`), then one column of expression per gene (named
#'   by symbol, or by Ensembl ID when genes were given as ensembl_gene_ids).
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' # needs an expression source, e.g.
#' # options(GAMBLR.data.expression_source = "<folder or URL>")
#' meta <- get_gambl_metadata()
#' myd88 <- get_gene_expression(meta, hugo_symbols = c("MYD88", "BCL2"),
#'                              lazy_join = TRUE)
#' }
get_gene_expression = function(these_samples_metadata,
                               hugo_symbols,
                               ensembl_gene_ids,
                               all_genes = FALSE,
                               verbose = FALSE,
                               engine = "grep",
                               format = "wide",
                               lazy_join = FALSE,
                               arbitrarily_pick = FALSE,
                               HGNC = FALSE,
                               ...){
  check_excess_params(...)
  if(!format %in% c("wide", "long")) stop("format must be \"wide\" or \"long\"")
  if(all_genes){
    stop("all_genes isn't available in GAMBLR.open (it would download all expression data). ",
         "Request the genes you need with hugo_symbols or ensembl_gene_ids.")
  }
  if(!missing(hugo_symbols) && !missing(ensembl_gene_ids)){
    stop("give hugo_symbols or ensembl_gene_ids, not both")
  }
  if(missing(hugo_symbols) && missing(ensembl_gene_ids)){
    stop("You must provide a vector of ensembl_gene_ids or hugo_symbols")
  }

  genes <- GAMBLR.data::expression_genes()
  by_symbol <- !missing(hugo_symbols)
  if(by_symbol){
    requested <- unique(as.character(hugo_symbols))
    found <- genes[genes$hgnc_symbol %in% requested, c("ensembl_gene_id", "hgnc_symbol"), drop = FALSE]
    multi <- unique(found$hgnc_symbol[duplicated(found$hgnc_symbol)])
    if(length(multi) > 0){
      if(!arbitrarily_pick){
        detail <- vapply(multi, function(s) paste0(s, " (", paste(found$ensembl_gene_id[found$hgnc_symbol == s],
                                                                  collapse = ", "), ")"), character(1))
        stop("these symbols belong to more than one Ensembl gene: ", paste(detail, collapse = "; "),
             ". Request the gene you want with ensembl_gene_ids, or set arbitrarily_pick = TRUE to use the first.")
      }
      if(verbose) message("using the first Ensembl gene for: ", paste(multi, collapse = ", "))
      found <- found[!duplicated(found$hgnc_symbol), , drop = FALSE]
    }
    missing_genes <- setdiff(requested, found$hgnc_symbol)
  }else{
    requested <- unique(sub("\\..*$", "", as.character(ensembl_gene_ids)))
    found <- genes[genes$ensembl_gene_id %in% requested, c("ensembl_gene_id", "hgnc_symbol"), drop = FALSE]
    missing_genes <- setdiff(requested, found$ensembl_gene_id)
  }
  if(length(missing_genes) > 0){
    warning("not in the expression data: ", paste(missing_genes, collapse = ", "), call. = FALSE)
  }
  if(nrow(found) == 0) stop("none of the requested genes are in the expression data")

  samples <- GAMBLR.data::expression_samples()
  if(!missing(these_samples_metadata)){
    ids <- unique(these_samples_metadata$sample_id)
    keep <- samples$mrna_sample_id %in% ids | samples$genome_sample_id %in% ids | samples$capture_sample_id %in% ids
    samples <- samples[keep, , drop = FALSE]
    if(verbose) message(nrow(samples), " RNA samples are linked to the samples in these_samples_metadata")
  }
  sample_cols <- c("mrna_sample_id", "patient_id",
                   if(lazy_join) c("genome_sample_id", "capture_sample_id"),
                   "protocol", "ffpe_or_frozen", "study")
  samples <- samples[, sample_cols, drop = FALSE]

  values <- GAMBLR.data::get_expression_values(found$ensembl_gene_id)
  values <- values[!is.na(values$expression), , drop = FALSE]
  values$Hugo_Symbol <- values$hgnc_symbol
  expression_long <- dplyr::inner_join(samples,
                                       values[, c("mrna_sample_id", "ensembl_gene_id", "Hugo_Symbol", "expression")],
                                       by = "mrna_sample_id")
  if(format == "long") return(expression_long)

  name_col <- if(by_symbol) "Hugo_Symbol" else "ensembl_gene_id"
  expression_long %>%
    dplyr::select(-dplyr::all_of(setdiff(c("Hugo_Symbol", "ensembl_gene_id"), name_col))) %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(name_col), values_from = "expression") %>%
    as.data.frame()
}

#' @title Check Gene Expression.
#'
#' @description List the RNA-seq samples available to [get_gene_expression()].
#'
#' @details Same name as `GAMBLR.results::check_gene_expression()`. One row
#'   per RNA sample (one per biopsy) linked to a bundled DNA sample.
#'
#' @param show_linkages Set to TRUE (default) to include genome_sample_id and
#'   capture_sample_id, the DNA samples of the same biopsy.
#' @param ... Any additional parameters.
#'
#' @return A data frame with mrna_sample_id, patient_id, protocol,
#'   ffpe_or_frozen and study (plus genome_sample_id and capture_sample_id).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rna <- check_gene_expression()
#' dplyr::count(rna, study)
#' }
check_gene_expression = function(show_linkages = TRUE, ...){
  check_excess_params(...)
  samples <- GAMBLR.data::expression_samples()
  if(!show_linkages) samples <- samples[, setdiff(colnames(samples), c("genome_sample_id", "capture_sample_id"))]
  samples
}
