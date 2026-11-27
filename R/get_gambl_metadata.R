#' @title Get GAMBL Metadata.
#'
#' @description Convenience function for loading the sample metadata.
#'
#' @details This bare bones function was developed to retrieve metadata for
#' non-GSC-users. Specify the seq type (`seq_type_filter`) for the samples you
#' want returned as the only argument.
#' It relies on the bundled metadata in this package.
#' Specify `case_set` argument to retreive samples from particular study.
#' Currently supported case_sets are: FL_Dreval (FL samples from Dreval et al),
#' DLBCL_Dreval (DLBCL samples from Dreval et al), FL-DLBCL-study (all samples
#' from Dreval et al), DLBCL_Arthur (all samples from Arthur et al study),
#' DLBCL_Hilton (all samples from Hilton et al DLBCL Trios study),
#' DLBCL_cell_lines (5 DLBCL cell lines), DLBCL_Chapuy (all samples from Chapuy
#' et al study), DLBCL_Schmitz (all samples from Schmitz et al study),
#' DLBCL_Reddy (all samples from Reddy et al study), DLBCL_Thomas (HTMCP DLBCLs
#' from Thomas et al study), BL_Thomas (BL samples from Thomas et al study)
#'
#' @param seq_type_filter Specify the seq type you want to return metadata for.
#' Default is "genome".
#' @param case_set Optionally specify study details to return samples from a
#' particular case set. See function description for supported case sets.
#' @param with_outcomes Set to TRUE to add the published patient-level
#' outcome columns of [get_gambl_outcomes] (e.g. CODE_OS, OS_YEARS,
#' CODE_PFS, PFS_YEARS, IPI), joined by patient_id. Patients without
#' released outcome data have NA. Default: FALSE.
#' @param ... Any additional parameters.
#'
#' @return A data frame with metadata, tailored for user without GSC access.
#'
#' \describe{
#'   \item{sample_id}{A unique identifier for the sample analyzed.}
#'   \item{seq_type}{The assay type used to produce this data (one of "genome","capture, "mrna", "promethION")}
#'   \item{study}{Which study/cohort block this sample was assembled from (e.g. "DLBCL_Arthur", "FL_Dreval"). Used internally by the `case_set` filters below.}
#'   \item{patient_id}{The anonymized unique identifier for this patient. For BC samples, this will be Res ID.}
#'   \item{pathology}{The diagnosis or pathology for the sample}
#'   \item{biopsy_res}{Biopsy-level result/classification identifier.}
#'   \item{genome_build}{The name of the genome reference the data were aligned to.}
#'   \item{pairing_status}{Matching status of the sample}
#'   \item{Tumor_Sample_Barcode}{Duplicate of sample_id for simplifying joins to MAF data frames}
#'   \item{cohort}{Name for a group of samples that were added together (usually from a single study), often in the format {pathology}_{cohort_descriptor}.}
#'   \item{COO_consensus}{Consensus call of COO between different sources.}
#'   \item{DHITsig_consensus}{Consensus call of DHIT signature (DHITsig, also called DZsig) status between different sources: POS, NEG or UNCLASS. DHITsig is defined within GCB DLBCL, but this column also has calls for ABC and COO-unclassified tumours, so comparing DHITsig groups with it across all DLBCL mostly compares GCB with ABC. Restrict to GCB tumours (\code{COO_consensus == "GCB"}) or use consensus_coo_dhitsig.}
#'   \item{consensus_coo_dhitsig}{Cell of origin and DHITsig combined: "ABC" for ABC DLBCLs, otherwise the DHITsig call ("DHITsigPos" or "DHITsigNeg"); NA outside DLBCL. DHITsig is a GCB signature, so compare DHITsig groups with this column rather than DHITsig_consensus, which also labels ABC tumours. Present in gambl_mutations.db builds made after 2026-11-26.}
#'   \item{EBV_status_inf}{Inferred EBV status of the tumor}
#'   \item{ffpe_or_frozen}{Whether the nucleic acids were extracted from a frozen or FFPE sample}
#'   \item{fl_grade}{Grade of FL samples}
#'   \item{bcl2_ba}{Result from breakapart FISH for BCL2 locus}
#'   \item{bcl2_cn}{Result from copy number FISH for BCL2 locus}
#'   \item{bcl6_ba}{Result from breakapart FISH for BCL6 locus}
#'   \item{bcl6_cn}{Result from copy number FISH for BCL6 locus}
#'   \item{myc_ba}{Result from breakapart FISH for MYC locus}
#'   \item{myc_cn}{Result from copy number FISH for MYC locus}
#'   \item{lymphgen}{LymphGen label}
#'   \item{lymphgen_cnv_noA53}{LymphGen label using model with CNV but excluding A53 class}
#'   \item{lymphgen_no_cnv}{LymphGen label using model without CNV}
#'   \item{lymphgen_with_cnv}{LymphGen label using model with CNV}
#'   \item{lymphgen_wright}{The LymphGen call for this sample from Wright et all (if applicable)}
#'   \item{normal_sample_id}{Sample id for normal tissue used in the analysis}
#'   \item{sex}{The biological sex of the patient (M, F or NA): one value per patient, from GAMBL's patient-level table.}
#'   \item{time_point}{Timing of biopsy in increasing alphabetical order (A = diagnosis, B = first relapse etc)}
#'   \item{TotalDuplicatedreads}{QC metric: total duplicated reads for this sample.}
#'   \item{TotalReads}{QC metric: total reads for this sample.}
#'   \item{TotalUniquelyMapped}{QC metric: total uniquely-mapped reads for this sample.}
#'   \item{TotalUnmappedreads}{QC metric: total unmapped reads for this sample.}
#'   \item{transformation}{Transformation status of the biopsy (BC Cancer patients; NA elsewhere): "de novo", "pre-transformation" (an indolent biopsy before the patient's first aggressive biopsy), "transformed" (an aggressive biopsy of a patient who transformed) or "post-transformation".}
#'   \item{genetic_subgroup}{Genetic subgroup from the source paper: Thomas et al. BL subgroups (DGG-BL, IC-BL, Q53-BL, DLBCL-A/B/C) and Dreval et al. FL (cFL, dFL); NA for other studies.}
#'   \item{MeanCorrectedCoverage}{QC metric: mean coverage after duplicate removal.}
#'   \item{n_ssm}{Number of SLMS-3 somatic SNV and indel calls. Not comparable between genomes and exomes.}
#'   \item{burden_per_mb}{Mutation burden: calls per Mb of the territory (burden_territory). Includes germline variants for unmatched samples (pairing_status); NA for samples under 15x MeanCorrectedCoverage.}
#'   \item{burden_snv_per_mb}{As burden_per_mb, SNVs only.}
#'   \item{burden_indel_per_mb}{As burden_per_mb, indels only.}
#'   \item{burden_territory}{Region the burden was measured over: "genome" (all calls over the non-N genome), "capture_targets" (calls in the capture kit's targets) or "coding_fallback" (coding calls over the protein-coding CDS, for exomes of unknown kit).}
#'   \item{burden_mb}{Size of that territory in Mb.}
#'   \item{burden_kit}{Capture kit whose targets were the territory, for "capture_targets".}
#' }
#'
#' The exact column set is controlled by
#' `data-raw/public_sample_meta_columns.txt` in GAMBLR.data -- see that
#' file if a column is missing here after a rebuild.
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' #return metadata for genome samples (here, the parameter is redundant because 
#' # 'genome' is the default)
#' genome_meta = get_gambl_metadata(seq_type_filter = "genome")
#'
#' #return metadata for capture samples.
#' capture_meta = get_gambl_metadata(seq_type_filter = "capture")
#'
#' #If you want metadata for genome and capture samples you can provide a vector of seq types
#' all_meta = get_gambl_metadata(seq_type_filter = c("genome", "capture"))
#'
#' dplyr::group_by(all_meta,cohort,seq_type) %>% 
#'     dplyr::count()
#'
get_gambl_metadata = function(
    seq_type_filter = c("genome","capture"),
    case_set,
    with_outcomes = FALSE,
    ...
){

    #check if any invalid parameters are provided
    check_excess_params(...)

    if (!isTRUE(getOption("GAMBLR.open.shown_metadata_msg"))) {
      message("Using the bundled sample_meta table in GAMBLR.data...")
      options(GAMBLR.open.shown_metadata_msg = TRUE)
    }
    con <- GAMBLR.data::gambl_mutations_db()
    metadata <- dplyr::tbl(con, "sample_meta") %>%
            dplyr::filter(seq_type %in% seq_type_filter) %>%
            dplyr::collect()


    if(!missing(case_set)){

        # pre-defined case sets
        if(case_set == "FL_Dreval"){
            metadata <- metadata %>%
                dplyr::filter(study == "FL_Dreval", pathology == "FL")
        }else if(case_set == "DLBCL_Dreval"){
            metadata <- metadata %>%
                dplyr::filter(study == "FL_Dreval", pathology == "DLBCL")
        }else if(case_set == "FL-DLBCL-study"){
            metadata <- metadata %>%
                dplyr::filter(study == "FL_Dreval")
        }else if(case_set == "DLBCL_Arthur"){
            metadata <- metadata %>%
                dplyr::filter(study == "DLBCL_Arthur")
        }else if(case_set == "DLBCL_Hilton"){
            metadata <- metadata %>%
                dplyr::filter(study == "DLBCL_Hilton")
        }else if(case_set == "DLBCL_cell_lines"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "DLBCL_cell_lines")
        }else if(case_set == "DLBCL_Chapuy"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "dlbcl_chapuy")
        }else if(case_set == "DLBCL_Schmitz"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "dlbcl_schmitz")
        }else if(case_set == "DLBCL_Reddy"){
            metadata <- metadata %>%
                dplyr::filter(cohort == "dlbcl_reddy")
        }else if(case_set == "BL_Thomas"){
            metadata <- metadata %>%
                dplyr::filter(study == "BL_Thomas")
        }else if(case_set == "DLBCL_Thomas"){
            metadata <- metadata %>%
                dplyr::filter(study == "DLBCL_Thomas")
        }else{
            message(paste("case set", case_set, "not available"))
            return()
        }
    }

    if(with_outcomes){
        # by patient_id only, so no other shared column can become a join key
        outcomes <- get_gambl_outcomes()
        outcomes <- dplyr::select(outcomes, -dplyr::any_of(setdiff(colnames(metadata), "patient_id")))
        metadata <- dplyr::left_join(metadata, outcomes, by = "patient_id")
    }

    #ensure only unique rows are returned
    return(unique(metadata))
}
