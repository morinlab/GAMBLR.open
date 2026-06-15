# Get all Coding SSMs

Retrieve all coding SSMs from GAMBL in MAF-like format, regardless of
seq_type.

## Usage

``` r
get_all_coding_ssm(
  these_samples_metadata = NULL,
  include_silent = FALSE,
  projection = "grch37"
)
```

## Arguments

- these_samples_metadata:

  Supply a metadata table containing the sample/seq_type combinations
  you want.

- include_silent:

  If set to TRUE, silent/synonymous mutations in the coding regions will
  also be returned.

- projection:

  The desired genome build 'grch37' or 'hg38' are allowed. Default is
  grch37

## Value

A data frame containing all the MAF data columns (one row per mutation).

## Details

Effectively retrieve coding SSM calls from one or all DNA seq_type. For
additional optional arguments, see
[GAMBLR.results::get_coding_ssm](https://rdrr.io/pkg/GAMBLR.results/man/get_coding_ssm.html)

## Examples

``` r
suppressPackageStartupMessages(library(GAMBLR.open))
my_meta = get_gambl_metadata(seq_type_filter = c("genome","capture"))
my_meta = check_and_clean_metadata(my_meta,duplicate_action="keep_first")
#> Duplicate rows (keeping first occurrence) for 'sample_id' and 'seq_type' have been dropped.
maf_all_seqtype = get_all_coding_ssm(my_meta)

table(maf_all_seqtype$maf_seq_type)
#> 
#> capture  genome 
#>   57951 3185525 

# most common mutations by gene and Variant_Classification
dplyr::group_by(maf_all_seqtype,
                Hugo_Symbol,
                Variant_Classification) %>% 
  dplyr::count() %>% 
  dplyr::arrange(desc(n))
#> genomic_data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>    Hugo_Symbol Variant_Classification       n
#> 1      Unknown                    IGR 1306295
#> 2        PTPRD                 Intron   21916
#> 3        IGLL5                 Intron   12689
#> 4         PCLO                 Intron    6479
#> 5        UNC5D                 Intron    5997
#> 6         BCL6                 Intron    5579
#> 7        ROBO2                 Intron    4850
#> 8         BCL2                 Intron    4810
#> 9        PTPRK                 Intron    4798
#> 10      BRINP3                 Intron    4774
```
