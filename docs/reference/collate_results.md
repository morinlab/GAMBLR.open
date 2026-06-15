# Collate Results

Bring together collated results for a selection of gambl samples.

## Usage

``` r
collate_results(
  sample_table,
  these_samples_metadata,
  join_with_full_metadata = FALSE,
  seq_type_filter = c("genome", "capture"),
  ...
)
```

## Arguments

- sample_table:

  A vector of characters with sample IDs, or a data frame with sample
  IDs in a column (sample_id). If provided, this will overwrite any
  sample subsets provided these_samples_metadata.

- these_samples_metadata:

  A metadata table with sample IDs of interest. If not provided, the
  function will get metadata for all available samples. This parameter
  is intended to use in combination with `join_with_full_metadata`.

- join_with_full_metadata:

  Set to TRUE to horizontally expand metadata with QC results. Default
  is FALSE. If `these_samples_metadata` is provided, collated resutls
  will be added to this metadata table. If not provided, the function
  will join collated results with all available metadata in the
  specified seq_type (`seq_type_filter`).

- seq_type_filter:

  Filtering criteria for `get_gambl_metadata` if
  `these_samples_metadata` is not provided, default is genomes and
  captures.

- ...:

  Any additional parameters.

## Value

A data frame with collated results.

## Details

Currently, this function only gathers QC metrics (`mirage_metrics`) as
the only collated result. Potentially, in the future, additional
collated results can be added by this function as well.

## Examples

``` r
#load packages
suppressPackageStartupMessages(library(dplyr))

#return collated results for all available samples
all_collated = collate_results()
#> Using the bundled collated results in GAMBLR.data...

#return available collated results for a metadata subset
fl_collated = collate_results(
 these_samples_metadata = get_gambl_metadata(
   seq_type_filter = "genome") %>% 
   dplyr::filter(pathology == "FL"))

#horizontally expand a metadata subset with collated results
fl_meta_collated = collate_results(
 join_with_full_metadata = TRUE, 
 these_samples_metadata = get_gambl_metadata(
   seq_type_filter = "genome") %>% 
   dplyr::filter(pathology == "FL"))

#horizontally expand all available metadata with collated results
all_meta_collated = collate_results(join_with_full_metadata = TRUE)
```
