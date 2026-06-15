# Get GAMBL Metadata.

Convenience function for loading the sample metadata.

## Usage

``` r
get_gambl_metadata(seq_type_filter = c("genome", "capture"), case_set, ...)
```

## Arguments

- seq_type_filter:

  Specify the seq type you want to return metadata for. Default is
  "genome".

- case_set:

  Optionally specify study details to return samples from a particular
  case set. See function description for supported case sets.

- ...:

  Any additional parameters.

## Value

A data frame with metadata, tailored for user without GSC access.

- compression:

  Format of the original data used as input for our analysis pipelines
  (cram, bam or fastq)

- bam_available:

  Whether or not this file was available when last checked.

- patient_id:

  The anonymized unique identifier for this patient. For BC samples,
  this will be Res ID.

- sample_id:

  A unique identifier for the sample analyzed.

- seq_type:

  The assay type used to produce this data (one of "genome","capture,
  "mrna", "promethION")

- genome_build:

  The name of the genome reference the data were aligned to.

- cohort:

  Name for a group of samples that were added together (usually from a
  single study), often in the format pathology_cohort_descriptor.

- pathology:

  The diagnosis or pathology for the sample

- time_point:

  Timing of biopsy in increasing alphabetical order (A = diagnosis, B =
  first relapse etc)

- ffpe_or_frozen:

  Whether the nucleic acids were extracted from a frozen or FFPE sample

- COO_consensus:

  Consensus call of COO between different sources.

- DHITsig_consensus:

  Consensus call of DHIT signature status between different sources.

- EBV_status_inf:

  Inferred EBV status of the tumor

- lymphgen_no_cnv:

  LymphGen label using model without CNV

- lymphgen_with_cnv:

  LymphGen label using model with CNV

- lymphgen_cnv_noA53:

  LymphGen label using model with CNV but excluding A53 class

- lymphgen_wright:

  The LymphGen call for this sample from Wright et all (if applicable)

- fl_grade:

  Grade of FL samples

- normal_sample_id:

  Sample id for normal tissue used in the analysis

- pairing_status:

  Matching status of the sample

- lymphgen:

  LymphGen label

- molecular_BL:

  label of the sample according to the molecular BL classifier

- Tumor_Sample_Barcode:

  Duplicate of sample_id for simplifying joins to MAF data frames

- pathology_rank:

  Numeric rank for consistent ordering of samples by pathology

- hiv_status:

  HIV status of the sample

- age_group:

  Adult_BL or Pediatric_BL or Other, specific to the BLGSP study

- sex:

  The biological sex of the patient, if available. Allowable options: M,
  F, NA

## Details

This bare bones function was developed to retrieve metadata for
non-GSC-users. Specify the seq type (`seq_type_filter`) for the samples
you want returned as the only argument. It relies on the bundled
metadata in this package. Specify `case_set` argument to retreive
samples from particular study. Currently supported case_sets are:
FL_Dreval (FL samples from Dreval et al), DLBCL_Dreval (DLBCL samples
from Dreval et al), FL-DLBCL-study (all samples from Dreval et al),
DLBCL_Arthur (all samples from Arthur et al study), DLBCL_Hilton (all
samples from Hilton et al DLBCL Trios study), DLBCL_cell_lines (5 DLBCL
cell lines), DLBCL_Chapuy (all samples from Chapuy et al study),
DLBCL_Schmitz (all samples from Schmitz et al study), DLBCL_Reddy (all
samples from Reddy et al study), DLBCL_Thomas (HTMCP DLBCLs from Thomas
et al study), BL_Thomas (BL samples from Thomas et al study)

## Examples

``` r
#return metadata for genome samples (here, the parameter is redundant because 
# 'genome' is the default)
genome_meta = get_gambl_metadata(seq_type_filter = "genome")

#return metadata for capture samples.
capture_meta = get_gambl_metadata(seq_type_filter = "capture")

#If you want metadata for genome and capture samples you can provide a vector of seq types
all_meta = get_gambl_metadata(seq_type_filter = c("genome", "capture"))

dplyr::group_by(all_meta,cohort,seq_type) %>% 
    dplyr::count()
#> # A tibble: 21 × 3
#> # Groups:   cohort, seq_type [21]
#>    cohort             seq_type     n
#>    <chr>              <chr>    <int>
#>  1 BL_Adult           genome      91
#>  2 BL_Pediatric       genome     121
#>  3 BL_cell_lines      genome      22
#>  4 DLBCL_Gascoyne     genome      21
#>  5 DLBCL_GenomeCanada genome      59
#>  6 DLBCL_HTMCP        genome      43
#>  7 DLBCL_ICGC         genome      84
#>  8 DLBCL_LSARP_Trios  capture     12
#>  9 DLBCL_LSARP_Trios  genome     142
#> 10 DLBCL_Marra        genome      38
#> # ℹ 11 more rows
```
