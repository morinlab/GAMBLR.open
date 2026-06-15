# Annotate Hotspots.

Annotate MAF-like data frome with a hot_spot column indicating recurrent
mutations.

## Usage

``` r
annotate_hotspots(mutation_maf, ...)
```

## Arguments

- mutation_maf:

  A data frame in MAF format.

- ...:

  Any other parameter. These parameters will be ignored.

## Value

The same data frame with one additional column "hot_spot".

## Details

This function takes an already loaded MAF data frame with the
`mutation_maf` parameter.

## Examples

``` r
my_metadata = get_gambl_metadata()
#> Using the bundled metadata in GAMBLR.data...
all_coding_ssm = get_coding_ssm(these_samples_metadata = my_metadata,
                                projection = "grch37",
                                this_seq_type = "genome") %>%
                   dplyr::filter(Hugo_Symbol %in% c("EZH2",
                                 "MEF2B","MYD88","KMT2D")) %>%
                   dplyr::arrange(Hugo_Symbol)
#> Using the bundled SSM calls (.maf) calls in GAMBLR.data...
#> id_ease: WARNING! 1838 samples in the provided metadata were removed because their seq types are not the same as in the `seq_type` argument. Use `verbose = TRUE` to see their IDs.
#> after linking with metadata, we have mutations from 859 samples

hot_ssms = annotate_hotspots(all_coding_ssm)
hot_ssms %>% dplyr::filter(!is.na(hot_spot)) %>%
      dplyr::select(1:5,37,hot_spot)
#> genomic_data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome HGVSp_Short hot_spot
#> 1         EZH2              0      .     GRCh37          7     p.Y646H     TRUE
#> 2         EZH2              0      .     GRCh37          7     p.Y646S     TRUE
#> 3         EZH2              0      .     GRCh37          7     p.Y646C     TRUE
#> 4         EZH2              0      .     GRCh37          7     p.Y646N     TRUE
#> 5         EZH2              0      .     GRCh37          7     p.Y646N     TRUE
#> 6         EZH2              0      .     GRCh37          7     p.Y646N     TRUE
#> 7         EZH2              0      .     GRCh37          7     p.Y646S     TRUE
#> 8         EZH2              0      .     GRCh37          7     p.Y646S     TRUE
#> 9         EZH2              0      .     GRCh37          7     p.Y646N     TRUE
#> 10        EZH2              0      .     GRCh37          7     p.Y646N     TRUE
```
