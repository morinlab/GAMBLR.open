# Get SSM By Patients.

Get MAF-format data frame for more than one patient.

## Usage

``` r
get_ssm_by_patients(
  these_patient_ids,
  these_samples_metadata,
  projection = "grch37",
  this_seq_type = "genome",
  tool_name = "slms-3",
  this_study,
  verbose = FALSE,
  ...
)
```

## Arguments

- these_patient_ids:

  A vector of patient IDs that you want results for. The user can also
  use a metadata table that has been subset to the patient IDs of
  interest (see `these_samples_metadata`).

- these_samples_metadata:

  A metadata subset to contain the rows corresponding to the patients of
  interest. If the vector of patient IDs is missing
  (`these_patient_ids`), this function will default to all patient IDs
  in the metadata table given to this parameter.

- projection:

  Obtain variants projected to this reference (one of grch37 or hg38).
  Default is grch37.

- this_seq_type:

  The seq type you want results for. Default is "genome".

- tool_name:

  Optionally specify which tool to report variant from. The default is
  slms-3, also supports "publication" to return the exact variants as
  reported in the original papers.

- this_study:

  Optionally specify first name of the author for the paper from which
  the variants should be returned for. This parameter can either be a
  vector of indexes (integer) or a vector of characters (matching
  columns in MAF).

- verbose:

  Set to FALSE to minimize the output to console. Default is TRUE. This
  parameter also dictates the verbosity of any helper function
  internally called inside the main function.

- ...:

  Any additional parameters.

## Value

A data frame with SSM calls for the selected patients in MAF format.

## Details

This function returns variants from a set of patients. This function
internally calls get_ssm_by_samples. Thus, the main contents of this
function is to wrangle the provided patient IDs, so that the
corresponding sample IDs can be provided to the internal call of
`get_ssm_by_samples`. This function expects either a vector of patient
IDs (`these_patients_ids`) or an already subset metadata table
(`these_samples_metadata`).

## Examples

``` r
# Lets find which patient_id occur more than once in the metadata first
my_ids = get_gambl_metadata(seq_type_filter = c("genome","capture")) %>%
             dplyr::group_by(patient_id) %>%
             dplyr::tally() %>%
             dplyr::filter(n>1) %>%
             dplyr::pull(patient_id)

#now let's get every SSM for all samples from these patients
patient_maf = get_ssm_by_patients(these_patient_ids = my_ids)
#> Patient IDs and metadata were provided, this function will resort to all available patient IDs in the provided metadata.
patient_maf %>% dplyr::group_by(Tumor_Sample_Barcode) %>% 
                dplyr::count() %>% head()
#> genomic_data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>   Tumor_Sample_Barcode   n
#> 1      00-14595_tumorA 476
#> 2      00-14595_tumorB 596
#> 3      00-14595_tumorC 679
#> 4      00-14595_tumorD 679
#> 5      00-15201_tumorA 208
#> 6      00-15201_tumorB 142
```
