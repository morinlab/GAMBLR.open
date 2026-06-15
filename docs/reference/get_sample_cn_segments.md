# Get Sample CN Segments.

Get all segments for a single (or multiple) sample_id(s).

## Usage

``` r
get_sample_cn_segments(
  these_sample_ids = NULL,
  these_samples_metadata = NULL,
  projection = "grch37",
  this_seq_type = "genome",
  with_chr_prefix = FALSE,
  streamlined = FALSE,
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

- projection:

  Selected genome projection for returned CN segments. Default is
  "grch37".

- this_seq_type:

  Seq type for returned CN segments. Default is genome.

- with_chr_prefix:

  Set to TRUE to add a chr prefix to chromosome names. Default is FALSE.

- streamlined:

  Return a minimal output rather than full details. Default is FALSE.

- verbose:

  Set to FALSE to minimize the output to console. Default is TRUE. This
  parameter also dictates the verbosity of any helper function
  internally called inside the main function.

- ...:

  Any additional parameters.

## Value

A data frame of segments for a specific or multiple sample ID(s).

## Details

This function returns CN segments. This works for single sample or
multiple samples. Specify the sample IDs you are interested in with
`these_sample_ids` (as a vector of characters), Or call this function
with `these_samples_metadata` if you already have a metadata table
subset to the sample IDs of interest. If none of the above parameters
are specified, the function will return CN segments for available
samples (from get_gambl_metadata). Note, this. function internally calls
[id_ease](https://morinlab.github.io/GAMBLR.open/reference/id_ease.md)
for dealing with sample IDs and metadata tables.

## Examples

``` r
#load pacakges
suppressPackageStartupMessages(library(dplyr))

#get CN segments for one sample
dohh2_segs = get_sample_cn_segments(these_sample_ids = "DOHH-2",
                                    projection = "hg38",
                                    streamlined = TRUE)

#get CN segments for DLBCL cell line
cell_line_meta = GAMBLR.data::sample_data$meta %>%
  dplyr::filter(cohort == "DLBCL_cell_lines")

dlbcl_segs = get_sample_cn_segments(these_samples_metadata = cell_line_meta,
                                    streamlined = TRUE)
```
