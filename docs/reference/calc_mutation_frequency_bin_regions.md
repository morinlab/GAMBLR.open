# Mutation counts across sliding windows for multiple regions.

Obtain a long tidy or wide matrix of mutation counts across sliding
windows for multiple regions.

## Usage

``` r
calc_mutation_frequency_bin_regions(
  regions_list = NULL,
  regions_bed = NULL,
  these_samples_metadata = NULL,
  these_sample_ids = NULL,
  this_seq_type = "genome",
  maf_data = NULL,
  projection = "grch37",
  region_padding = 1000,
  drop_unmutated = FALSE,
  skip_regions = NULL,
  only_regions = NULL,
  slide_by = 100,
  window_size = 500,
  return_format = "wide",
  ...
)
```

## Arguments

- regions_list:

  Named vector of regions in the format c(name1 = "chr:start-end", name2
  = "chr:start-end"). If neither `regions` nor `regions_bed` is
  specified, the function will use GAMBLR aSHM region information.

- regions_bed:

  Data frame of regions with four columns (chrom, start, end, name).

- these_samples_metadata:

  Metadata with at least sample_id column. If not providing a maf data
  frame, seq_type is also required.

- these_sample_ids:

  Vector of sample IDs. Metadata will be subset to sample IDs present in
  this vector.

- this_seq_type:

  Optional vector of seq_types to include in heatmap. Default "genome".
  Uses default seq_type priority for samples with \>1 seq_type.

- maf_data:

  Optional maf data frame. Will be subset to rows where
  Tumor_Sample_Barcode matches provided sample IDs or metadata table. If
  not provided, maf data will be obtained with get_ssm_by_regions().

- projection:

  Genome build the function will operate in. Ensure this matches your
  provided regions and maf data for correct chr prefix handling. Default
  "grch37".

- region_padding:

  Amount to pad the start and end coordinates by. Default 1000.

- drop_unmutated:

  Whether to drop bins with 0 mutations. If returning a matrix format,
  this will only drop bins with no mutations in any samples.

- skip_regions:

  Optional character vector of genes to exclude from the default aSHM
  regions.

- only_regions:

  Optional character vector of genes to include from the default aSHM
  regions.

- slide_by:

  Slide size for sliding window. Default 100.

- window_size:

  Size of sliding window. Default 500.

- return_format:

  Return format of mutations. Accepted inputs are "long" and "wide".
  Long returns a data frame of one sample ID/window per row. Wide
  returns a matrix with one sample ID per row and one window per column.
  Using the "wide" format will retain all samples and windows regardless
  of the drop_unmutated or min_count_per_bin parameters. Default wide.

- ...:

  Any additional parameters.

## Value

A table of mutation counts for sliding windows across one or more
regions. May be long or wide.

## Details

This function takes a metadata table with `these_samples_metadata`
parameter and internally calls `calc_mutation_frequency_bin_region`
(that internally calls `get_ssm_by_regions`). to retrieve mutation
counts for sliding windows across one or more regions. May optionally
provide any combination of a maf data frame, existing metadata, or a
regions data frame or named vector.

## Examples

``` r
 #load metadata.
 my_meta = get_gambl_metadata()
 dlbcl_bl_meta = dplyr::filter(my_meta, pathology %in% c("DLBCL", "BL"))


 #get ashm regions
 some_regions = create_bed_data(grch37_ashm_regions,
                                fix_names = "concat",
                                concat_cols = c("gene","region"),
                                sep="-")
 print(some_regions)
#> BED Data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>    chrom     start       end           name region regulatory_comment
#> 1      1   6661482   6662702     KLHL21-TSS    TSS               <NA>
#> 2      1  23885584  23885835        ID3-TSS    TSS               <NA>
#> 3      1  28832551  28836339       RCC1-TSS    TSS               <NA>
#> 4      1  31229012  31232011     LAPTM5-TSS    TSS               <NA>
#> 5      1 150550814 150552135    MCL1-intron intron               <NA>
#> 6      2  88904839  88909096 EIF2Ak3-intron intron               <NA>
#> 7      2  88925456  88927581    EIF2Ak3-TSS    TSS               <NA>
#> 8      2 232572640 232574297       PTMA-TSS    TSS               <NA>
#> 9      1 157669490 157671299      FCRL3-TSS    TSS               <NA>
#> 10     1 203274698 203275778    BTG2-intron intron    active_promoter
 mut_count_matrix <- calc_mutation_frequency_bin_regions(
   these_samples_metadata = dlbcl_bl_meta,
   regions_bed = some_regions
 )
dim(mut_count_matrix)
#> [1]   763 14408
tail(mut_count_matrix[,c(1:10)])
#> # A tibble: 6 × 10
#>   sample_id `1_6661382` `1_6661482` `1_6661582` `1_6661682` `1_6661782`
#>   <chr>           <int>       <int>       <int>       <int>       <int>
#> 1 SP59456             0           0           0           0           0
#> 2 SP59460             0           0           0           0           0
#> 3 SU-DHL-10           0           0           0           0           0
#> 4 SU-DHL-4            0           0           0           0           0
#> 5 Seraphina           0           0           0           0           0
#> 6 Thomas              0           0           0           0           0
#> # ℹ 4 more variables: `1_6661882` <int>, `1_6661982` <int>, `1_6662082` <int>,
#> #   `1_6662182` <int>
```
