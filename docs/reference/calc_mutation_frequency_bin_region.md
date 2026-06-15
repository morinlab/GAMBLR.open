# Calculate Mutation Frequency By Sliding Window.

Count the number of mutations in a sliding window across a region for
all samples.

## Usage

``` r
calc_mutation_frequency_bin_region(
  region,
  these_samples_metadata = NULL,
  these_sample_ids = NULL,
  this_seq_type = "genome",
  maf_data = NULL,
  projection = "grch37",
  slide_by = 100,
  window_size = 1000,
  return_format = "long",
  min_count_per_bin = 0,
  return_count = TRUE,
  drop_unmutated = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- region:

  A string describing a genomic region in the "chrom:start-end" format.
  The region must be specified in this format OR as separate chromosome,
  start_pos, end_pos arguments.

- these_samples_metadata:

  Optional data frame containing a sample_id column. If not providing a
  maf file, seq_type is also a required column.

- these_sample_ids:

  Optional vector of sample IDs. Output will be subset to IDs present in
  this vector.

- this_seq_type:

  Optional vector of seq_types to include in heatmap. Default is
  "genome". Uses default seq_type priority for samples with \>1
  seq_type.

- maf_data:

  Optional maf data frame. Will be subset to rows where
  Tumor_Sample_Barcode matches provided sample IDs or metadata table. If
  not provided, maf data will be obtained with get_ssm_by_regions().

- projection:

  Specify which genome build to use. Required. Default grch37.

- slide_by:

  Slide size for sliding window. Default 100.

- window_size:

  Size of sliding window. Default 1000.

- return_format:

  Return format of mutations. Accepted inputs are "long" and "wide".
  Long returns a data frame of one sample ID/window per row. Wide
  returns a matrix with one sample ID per row and one window per column.
  Using the "wide" format will retain all samples and windows regardless
  of the drop_unmutated or min_count_per_bin parameters.

- min_count_per_bin:

  Minimum counts per bin, default is 0. Setting this greater than 0 will
  drop unmutated windows only when return_format is long.

- return_count:

  Boolean statement to return mutation count per window (TRUE) or binary
  mutated/unmutated status (FALSE). Default is TRUE.

- drop_unmutated:

  Boolean for whether to drop windows with 0 mutations. Only effective
  with "long" return format.

- ...:

  Any additional parameters.

## Value

Either a matrix or a long tidy table of counts per window.

## Details

This function is called to return the mutation frequency for a given
region, either from a provided input maf data frame or from the GAMBL
maf data. Regions are specified with the `region` parameter.
Alternatively, the region of interest can also be specified by calling
the function with `chromosome`, `start_pos`, and `end_pos` parameters.
This function operates on a single region. To return a matrix of sliding
window counts over multiple regions, see
`calc_mutation_frequency_bin_regions`.

## Examples

``` r
myc_region = "8:128747680-128753674" 
myc_mut_freq = calc_mutation_frequency_bin_region(region = myc_region,
                                                  slide_by = 10,
                                                  window_size = 10000)
#> Using GAMBLR.data::get_ssm_by_region...
dplyr::arrange(myc_mut_freq,desc(mutation_count))
#> # A tibble: 612,000 × 3
#>    sample_id                 bin         mutation_count
#>    <chr>                     <chr>                <int>
#>  1 BLGSP-71-26-00399-01A-01E 8_128747680             41
#>  2 BLGSP-71-26-00399-01A-01E 8_128747690             41
#>  3 BLGSP-71-26-00399-01A-01E 8_128747700             41
#>  4 BLGSP-71-26-00399-01A-01E 8_128747710             41
#>  5 BLGSP-71-26-00399-01A-01E 8_128747720             41
#>  6 BLGSP-71-26-00399-01A-01E 8_128747730             41
#>  7 BLGSP-71-26-00399-01A-01E 8_128747740             41
#>  8 BLGSP-71-26-00399-01A-01E 8_128747750             41
#>  9 BLGSP-71-26-00399-01A-01E 8_128747760             41
#> 10 BLGSP-71-26-00399-01A-01E 8_128747770             41
#> # ℹ 611,990 more rows
```
