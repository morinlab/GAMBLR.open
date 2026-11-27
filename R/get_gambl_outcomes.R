#' @title Get GAMBL Outcomes.
#'
#' @description Get the published patient-level outcome data bundled in
#' GAMBLR.data.
#'
#' @details Reads the `patient_outcomes` table of the GAMBLR.data database
#' ([GAMBLR.data::gambl_mutations_db]): one row per patient, only for studies
#' that released outcome data, and only the endpoints each released:
#' \itemize{
#'   \item BL_Thomas (Thomas et al.): OS and PFS
#'   \item dlbcl_chapuy (Chapuy et al. 2018): OS, PFS and IPI
#'   \item dlbcl_schmitz (Schmitz et al. 2018): OS and PFS
#'   \item DLBCL_Hilton (the relapsed DLBCL tumour evolution study): OS
#' }
#' A patient in more than one of these studies has the endpoints of all of
#' them. Other patients have no row. Event codes are 1 for an event and 0 for
#' censored; times are from diagnosis. A censored time of 0 (no follow-up) is
#' NA.
#'
#' The function has the same name and arguments as
#' `GAMBLR.results::get_gambl_outcomes()`, which returns the full GAMBL
#' outcome table for GSC users. To add outcomes to sample metadata, use
#' `get_gambl_metadata(with_outcomes = TRUE)`, or join by `patient_id` only.
#'
#' @param patient_ids Optional vector of patient IDs to return.
#' @param time_unit Return follow-up times in years ("year", the default),
#'   months ("month") or days ("day", 365.25 per year).
#' @param censor_cbioportal Set to TRUE to add the OS and PFS status columns
#'   cBioPortal expects (OS_STATUS, DFS_STATUS and, with
#'   `time_unit = "month"`, DFS_MONTHS). Default is FALSE.
#' @param complete_missing Set to TRUE to add a row, censored at time 0, for
#'   every patient in `patient_ids` that has no outcome data. Patients with
#'   data keep their missing values as NA. Default is FALSE.
#' @param from_flatfile Accepted for compatibility with
#'   `GAMBLR.results::get_gambl_outcomes()`; ignored.
#' @param ... Any additional parameters.
#'
#' @return A data frame with one row per patient: patient_id and, depending
#'   on the endpoints released, CODE_OS, OS_YEARS (or OS_MONTHS/OS_DAYS),
#'   alive, CODE_PFS, PFS_YEARS (or PFS_MONTHS/PFS_DAYS) and IPI.
#'
#' \describe{
#'   \item{patient_id}{The anonymized unique identifier for this patient.}
#'   \item{CODE_OS}{Overall survival event: 1 = death from any cause, 0 = censored.}
#'   \item{OS_YEARS}{Overall survival time from diagnosis, in years.}
#'   \item{alive}{"alive" or "deceased" at last follow-up (from CODE_OS).}
#'   \item{CODE_PFS}{Progression-free survival event: 1 = progression or death, 0 = censored.}
#'   \item{PFS_YEARS}{Progression-free survival time from diagnosis, in years.}
#'   \item{IPI}{International Prognostic Index score, 0-5.}
#' }
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' outcomes <- get_gambl_outcomes()
#' dplyr::count(outcomes, alive)
#'
#' # sample metadata with outcomes
#' meta <- get_gambl_metadata(with_outcomes = TRUE)
#'
get_gambl_outcomes = function(patient_ids,
                              time_unit = "year",
                              censor_cbioportal = FALSE,
                              complete_missing = FALSE,
                              from_flatfile = TRUE,
                              ...){

  check_excess_params(...)
  if(!time_unit %in% c("year", "month", "day")){
    stop("time_unit must be one of \"year\", \"month\" or \"day\"")
  }

  con <- GAMBLR.data::gambl_mutations_db()
  if(!"patient_outcomes" %in% DBI::dbListTables(con)){
    stop("Your copy of the GAMBLR.data database has no patient_outcomes table (it predates outcome data). ",
         "Update it with GAMBLR.data::download_gambl_mutations_db().", call. = FALSE)
  }
  outcomes <- dplyr::tbl(con, "patient_outcomes") %>%
    dplyr::collect() %>%
    as.data.frame()

  if(!missing(patient_ids)){
    outcomes <- dplyr::filter(outcomes, patient_id %in% patient_ids)
    if(complete_missing){
      year_cols <- grep("_YEARS$", colnames(outcomes), value = TRUE)
      code_cols <- grep("^CODE_", colnames(outcomes), value = TRUE)
      added <- data.frame(patient_id = setdiff(unique(patient_ids), outcomes$patient_id))
      for(col in c(year_cols, code_cols)) added[[col]] <- rep(0, nrow(added))
      outcomes <- dplyr::bind_rows(outcomes, added)
    }
  }

  # follow-up times in months or days instead of years
  if(time_unit != "year"){
    year_cols <- grep("_YEARS$", colnames(outcomes), value = TRUE)
    factor <- if(time_unit == "month") 12 else 365.25
    suffix <- if(time_unit == "month") "_MONTHS" else "_DAYS"
    for(col in year_cols){
      outcomes[[sub("_YEARS$", suffix, col)]] <- outcomes[[col]] * factor
    }
    outcomes <- dplyr::select(outcomes, -dplyr::all_of(year_cols))
  }

  # censoring in the format cBioPortal expects
  if(censor_cbioportal){
    if("CODE_OS" %in% colnames(outcomes)){
      outcomes$OS_STATUS <- dplyr::case_when(outcomes$CODE_OS == 0 ~ "0:LIVING",
                                             outcomes$CODE_OS == 1 ~ "1:DECEASED")
    }
    if("CODE_PFS" %in% colnames(outcomes)){
      outcomes$DFS_STATUS <- dplyr::case_when(outcomes$CODE_PFS == 0 ~ "0:DiseaseFree",
                                              outcomes$CODE_PFS == 1 ~ "1:Recurred/Progressed")
    }
    if("PFS_MONTHS" %in% colnames(outcomes)) outcomes$DFS_MONTHS <- outcomes$PFS_MONTHS
  }

  return(outcomes)
}
