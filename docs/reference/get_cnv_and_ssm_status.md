# Get CNV and coding SSM combined status

For each specified chromosome region (gene name), return status 1 if the
copy number (CN) state is non-neutral, *i.e.* different from 2, or if
the region contains any coding simple somatic mutation (SSM).

## Usage

``` r
get_cnv_and_ssm_status(
  genes_and_cn_threshs,
  these_samples_metadata,
  maf_df,
  seg_data,
  cn_matrix,
  only_cnv = "none",
  genome_build = "grch37",
  include_hotspots = TRUE,
  review_hotspots = TRUE,
  adjust_for_ploidy = TRUE,
  include_silent = FALSE,
  this_seq_type
)
```

## Arguments

- genes_and_cn_threshs:

  A data frame with columns "gene_id" and "cn_thresh". The "gene_id"
  column stores gene symbols (characters) which determine the regions to
  return CNV and/or coding SSM status. The "cn_thresh" column stores
  integers that mean the maximum or minimum CN states to return status 1
  (contains CNV) for its respective gene. If this integer is below 2
  (neutral CN state for diploids), it is taken as the maximum (gene
  consider as tumor suppressor); if above 2, it is the minimum
  (oncogene); if equal to 2, do not consider CNV to return status.

- these_samples_metadata:

  The metadata for samples of interest to be included in the returned
  matrix. Can be created with `get_gambl_metadata` function.

- maf_df:

  Optional data frame containing the coding variants for your samples
  (i.e. output from `get_all_coding_ssm`)

- seg_data:

  Optionally provide the function with a data frame of segments that
  will be used instead of the GAMBL flatfiles

- cn_matrix:

  Instead of seg_data, you can provide a matrix of CN values for the
  samples in the metadata. See
  [GAMBLR.utils::segmented_data_to_cn_matrix](https://rdrr.io/pkg/GAMBLR.utils/man/segmented_data_to_cn_matrix.html)
  for more information on how to create this matrix.

- only_cnv:

  A vector of gene names indicating the genes for which only CNV status
  should be considered, ignoring SSM status. Set this argument to "all"
  or "none" (default) to apply this behavior to all or none of the
  genes, respectively.

- genome_build:

  Reference genome build. Possible values are "grch37" (default) or
  "hg38".

- include_hotspots:

  Logical parameter indicating whether hotspots object should also be
  tabulated. Default is TRUE.

- review_hotspots:

  Logical parameter indicating whether hotspots object should be
  reviewed to include functionally relevant mutations or rare
  lymphoma-related genes. Default is TRUE.

- adjust_for_ploidy:

  Set to FALSE to disable scaling of CN values by the genome-wide
  average per sample

- include_silent:

  Set to TRUE if you want Synonymous mutations to also be considered

- this_seq_type:

  Deprecated

## Value

A data frame with CNV and SSM combined status.

## Details

The user can choose from which regions are intended to return only copy
number variation (CNV) status, only coding SSM status, or at least the
presence of one of them. This behavior is controlled by the arguments
`genes_and_cn_threshs` (column `cn_thresh`) and `only_cnv`.

This function internally calls the `get_cn_states`, `get_ssm_by_samples`
and `get_coding_ssm_status`functions. Therefore, many of its arguments
are assigned to these functions. If needed, see the documentation of
these functions for more information.

In the case of returning NA values, this is due to the `get_cn_segments`
function not being able to internally return any copy number segments
from the specified chromosome region.

## Examples

``` r
suppressPackageStartupMessages(library(GAMBLR.open))
# Get sample metadata including a mix of seq_type
all_types_meta = suppressMessages(get_gambl_metadata()) %>% 
            dplyr::filter(pathology == "BL")
dplyr::group_by(all_types_meta, seq_type) %>% 
     dplyr::summarize(n=dplyr::n())
#> # A tibble: 1 × 2
#>   seq_type     n
#>   <chr>    <int>
#> 1 genome     234

# For MYC and SYNCRIP, return CNV and SSM combined status; for MIR17HG, 
# return only CNV status; for CCND3 return only SSM status
genes_and_cn_threshs = data.frame(
  gene_id=c("MYC", "MIR17HG", "CCND3","ID3","DDX3X", "SYNCRIP"),
  cn_thresh=c(3, 3, 2, 2, 2, 1)
)

genome_cnv_ssm_status = suppressMessages(get_cnv_and_ssm_status(
                           genes_and_cn_threshs,
                           genome_build = "hg38",
                           dplyr::filter(all_types_meta,seq_type=="genome"),
                           only_cnv = "MIR17HG"))

print(dim(genome_cnv_ssm_status))    
#> [1] 234   1
head(genome_cnv_ssm_status)   
#>                           MYC
#> Akata                       1
#> BL2                         1
#> BL30                        1
#> BL41                        1
#> BL70                        1
#> BLGSP-71-06-00001-01A-11D   0
colSums(genome_cnv_ssm_status)
#> MYC 
#> 173 
       
all_path_meta = suppressMessages(get_gambl_metadata())
all_path_meta = check_and_clean_metadata(all_path_meta,duplicate_action="keep_first")
dplyr::group_by(all_types_meta, seq_type) %>% 
     dplyr::summarize(n=dplyr::n())             
#> # A tibble: 1 × 2
#>   seq_type     n
#>   <chr>    <int>
#> 1 genome     234
all_pathology_status = suppressMessages(get_cnv_and_ssm_status(
                           genes_and_cn_threshs,
                           all_path_meta,
                           genome_build = "grch37",
                           only_cnv = "MIR17HG"))

print(dim(all_pathology_status))   
#> [1] 443   1
head(all_pathology_status)
#>           MYC
#> 01-20260T   1
#> 02-13135T   1
#> 02-20170T   0
#> 02-22991T   1
#> 03-34157T   0
#> 04-24029T   0
colSums(all_pathology_status)
#> MYC 
#> 150 
```
