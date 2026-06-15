# Get SSM By Region.

Retrieve all SSMs from the GAMBL database within a single genomic
coordinate range.

## Usage

``` r
get_ssm_by_region(
  these_sample_ids = NULL,
  these_samples_metadata = NULL,
  maf_data,
  chromosome,
  qstart,
  qend,
  region = "",
  streamlined = FALSE,
  projection = "grch37",
  this_seq_type = "genome",
  tool_name = "slms-3",
  this_study,
  verbose = FALSE,
  ...
)
```

## Arguments

- these_sample_ids:

  Optional, a vector of multiple sample_id (or a single sample ID as a
  string) that you want results for.

- these_samples_metadata:

  Optional, a metadata table (with sample IDs in a column) to subset the
  return to. If not provided (and if `these_sample_ids` is not
  provided), the function will return all samples from the specified
  seq_type in the metadata.

- maf_data:

  Optional data frame with mutations in MAF format. If user provides a
  maf, the function trusts that the user has already subset this to
  samples of interest, correct seq_type. i.e the following parameters
  are ignored; `these_samples_metadata`, `these_sample_ids`, and
  `this_seq_type`

- chromosome:

  The chromosome you are restricting to (with or without a chr prefix).

- qstart:

  Query start coordinate of the range you are restricting to.

- qend:

  Query end coordinate of the range you are restricting to.

- region:

  Region formatted like chrX:1234-5678 instead of specifying chromosome,
  start and end separately.

- streamlined:

  Return Start_Position and Tumor_Smaple_Barcode as the only two MAF
  columns. Default is FALSE.

- projection:

  Obtain variants projected to this reference (one of grch37 or hg38).

- this_seq_type:

  The seq_type you want back, default is genome.

- tool_name:

  Optionally specify which tool to report variant from. The default is
  slms-3, also supports "publication" to return the exact variants as
  reported in the original papers.

- this_study:

  Optionally specify first name of the author for the paper from which
  the variants should be returned for.

- verbose:

  Set to FALSE to prevent ANY message to be printed. In most cases, this
  parameter should be left to TRUE. The parameter was added to
  accommodate for noisy output when running this function in a loop for
  retrieving SSM for multiple regions get_ssm_by_regions.

- ...:

  Any additional parameters.

## Value

A data frame containing all mutations (MAF) in the specified region.

## Details

This function lets the user specify a region of interest for returning
SSM calls within that region. There are multiple ways a region can be
specified. For example, the user can provide the full region in a
"region" format (chr:start-end) to the `region` parameter. Or, the user
can provide chromosome, start and end coordinates individually with
`chr`, `start`, and `end` parameters.

## Examples

``` r
my_mutations = get_ssm_by_region(region = "chr8:128,723,128-128,774,067")

#specifying chromosome, start and end individually
my_mutations = get_ssm_by_region(chromosome = "8",
                                 qstart = 128723128,
                                 qend = 128774067)
```
