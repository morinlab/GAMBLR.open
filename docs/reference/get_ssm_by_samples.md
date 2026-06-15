# Get SSM By Samples.

Get the SSMs (i.e. load MAF) for a single sample or a collection of
samples.

## Usage

``` r
get_ssm_by_samples(
  these_sample_ids = NULL,
  these_samples_metadata = NULL,
  this_seq_type = "genome",
  projection = "grch37",
  tool_name = "slms-3",
  verbose = FALSE,
  ...
)
```

## Arguments

- these_sample_ids:

  A vector of one or more sample IDs that you want results for.

- these_samples_metadata:

  Optional, a metadata table (with sample_id column) to auto-subset the
  data to samples in that table before returning. If not provided and
  these_sample_ids is also not provided, the function will return SSM
  for all samples from the specified seq_type in the bundled metadata.

- this_seq_type:

  Default is genome.

- projection:

  The projection genome build. Supports hg38 and grch37.

- tool_name:

  Optionally specify which tool to report variant from. The default is
  slms-3, also supports "publication" to return the exact variants as
  reported in the original papers.

- verbose:

  Enable for debugging/noisier output.

- ...:

  Any additional parameters.

## Value

data frame in MAF format.

## Details

Retrieve a maf for a specific sample or a set of samples. Either specify
the sample IDs of interest with `these_sample_ids`. Or a metadata table
subset to the sample IDs of interest with `these_samples_metadata`.

## Examples

``` r
#Get genome-wide set of mutations from all DLBCL cell lines

# 1. get our metadata for the DLBCL cell lines
cell_line_meta = get_gambl_metadata() %>%
  dplyr::filter(cohort == "DLBCL_cell_lines")

# 2. get the SSMs for the DLBCL cell lines
dlbcl_maf = get_ssm_by_samples(these_samples_metadata = cell_line_meta)

# 3. have a look:
dlbcl_maf %>% dplyr::group_by(Tumor_Sample_Barcode) %>%
              dplyr::count()
#> genomic_data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>   Tumor_Sample_Barcode     n
#> 1               DOHH-2 22089
#> 2             OCI-Ly10 30051
#> 3              OCI-Ly3 31532
#> 4            SU-DHL-10 26855
#> 5             SU-DHL-4 32824
```
