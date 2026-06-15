# Get CN Segments.

Retrieve all copy number segments from the GAMBL outputs

## Usage

``` r
get_cn_segments(
  these_samples_metadata,
  projection = "grch37",
  this_seq_type,
  ...
)
```

## Arguments

- these_samples_metadata:

  User must provide a metadata table to restrict the data to the samples
  in your table. The metadata also ensures the proper handling of
  duplicate sample_id across seq_types and ensures the seq_type in the
  metadata faithfully represents the seq_type of the data

- projection:

  Desired genome coordinate system for returned CN segments. Default is
  "grch37".

- this_seq_type:

  Deprecated.

- ...:

  Additional parameters to be passed to the function.

## Value

A data frame with CN segments for the specified region.

## Details

This function merely loads and returns all the seg_data available for a
projection (genome build)

## Examples

``` r
# Example for the capture samples:

genome_metadata = get_gambl_metadata(seq_type_filter="genome")

genome_segments_hg38 = get_cn_segments(
                             these_samples_metadata = genome_metadata,
                             projection="hg38")

```
