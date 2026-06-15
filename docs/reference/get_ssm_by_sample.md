# Get SSM By Sample.

Get the SSMs (i.e. load MAF) for a single sample.

## Usage

``` r
get_ssm_by_sample(
  this_sample_id,
  these_samples_metadata = NULL,
  this_seq_type = "genome",
  projection = "grch37",
  ...
)
```

## Arguments

- this_sample_id:

  A single sample ID to retrieve mutations for.

- these_samples_metadata:

  Optional, a metadata table (with sample_id column) to auto-subset the
  data to samples in that table before returning.

- this_seq_type:

  Default is genome.

- projection:

  The projection genome build. Supports hg38 and grch37.

- ...:

  Any additional parameters passed to
  [get_ssm_by_samples](https://morinlab.github.io/GAMBLR.open/reference/get_ssm_by_samples.md).

## Value

data frame in MAF format.

## Details

Alias for
[get_ssm_by_samples](https://morinlab.github.io/GAMBLR.open/reference/get_ssm_by_samples.md)
accepting a single sample ID via `this_sample_id`, provided for API
compatibility with GAMBLR.results.

## Examples

``` r
maf = get_ssm_by_sample(this_sample_id = "DOHH-2")
```
