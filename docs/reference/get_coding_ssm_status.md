# Get Coding SSM Status.

Tabulate mutation status (SSM) for a set of genes.

## Usage

``` r
get_coding_ssm_status(
  gene_symbols,
  these_samples_metadata,
  maf_data,
  include_hotspots = TRUE,
  keep_multihit_hotspot = FALSE,
  review_hotspots = TRUE,
  genes_of_interest = c("FOXO1", "MYD88", "CREBBP"),
  genome_build,
  include_silent = FALSE,
  include_silent_genes,
  ...
)
```

## Arguments

- gene_symbols:

  A vector of gene symbols for which the mutation status will be
  tabulated. If not provided, lymphoma genes will be returned by
  default.

- these_samples_metadata:

  The metadata for samples of interest to be included in the returned
  matrix. Only the column "sample_id" is required. If not provided, the
  example metadata is used as default.

- maf_data:

  data frame in maf format. Must be in the grch37 projection.

- include_hotspots:

  Logical parameter indicating whether hotspots object should also be
  tabulated. Default is TRUE.

- keep_multihit_hotspot:

  Logical parameter indicating whether to keep the gene annotation as
  mutated when the gene has both hot spot and non-hotspot mutation.
  Default is FALSE. If set to TRUE, will report the number of
  non-hotspot mutations instead of tabulating for just mutation
  presence.

- review_hotspots:

  Logical parameter indicating whether hotspots object should be
  reviewed to include functionally relevant mutations or rare
  lymphoma-related genes. Default is TRUE.

- genes_of_interest:

  A vector of genes for hotspot review. Currently only FOXO1, MYD88, and
  CREBBP are supported.

- genome_build:

  Reference genome build for the coordinates in the MAF file. The
  default is inferred from maf_data.

- include_silent:

  Logical parameter indicating whether to include silent mutations into
  coding mutations. Default is FALSE.

- include_silent_genes:

  Optionally, provide a list of genes for which the Silent variants to
  be considered. If provided, the Silent variants for these genes will
  be included regardless of the include_silent argument.

- ...:

  Any other parameter. These parameters will be ignored.

## Value

A data frame with tabulated mutation status.

## Details

This function takes a data frame (in MAF-like format) and converts it to
a binary one-hot encoded matrix of mutation status for either a set of
user-specified genes (via gene_symbols) or, if no genes are provided,
default to all lymphoma genes. The default behaviour is to assign each
gene/sample_id combination as mutated only if there is a protein coding
mutation for that sample in the MAF but this can be configured to use
synonymous variants in some (via include_silent_genes) or all (via
include_silent) genes. This function also has other filtering and
convenience parameters giving the user full control of the return. For
more information, refer to the parameter descriptions and examples.
Currently only the grch37 genome build is supported for hotspot
annotation and review for this version of the function.

## Examples

``` r
coding_tabulated_df = get_coding_ssm_status(
 maf_data = get_coding_ssm(),
 gene_symbols = c("EZH2","KMT2D","CREBBP","MYC")
)
#> after linking with metadata, we have mutations from 859 samples
#> Joining with `by = join_by(sample_id)`
#> annotating hotspots
#> Joining with `by = join_by(sample_id)`
#> CREBBPHOTSPOT
#> OK
#> FOXO1HOTSPOT
#> MYD88HOTSPOT

head(coding_tabulated_df)
#>                   sample_id KMT2D CREBBP EZH2 MYC CREBBPHOTSPOT FOXO1HOTSPOT
#> 1                     Akata     0      0    0   1             0            1
#> 2                       BL2     0      0    0   1             0            1
#> 3                      BL30     1      0    0   1             0            0
#> 4                      BL41     0      0    0   1             0            0
#> 5                      BL70     0      0    0   1             0            0
#> 6 BLGSP-71-06-00001-01A-11D     0      0    0   0             0            1
#>   MYD88HOTSPOT
#> 1            0
#> 2            1
#> 3            0
#> 4            0
#> 5            0
#> 6            0

#all lymphoma genes from bundled NHL gene list
coding_tabulated_df = get_coding_ssm_status(
                           maf_data = get_coding_ssm()
                      )
#> after linking with metadata, we have mutations from 859 samples
#> No gene_symbols provided, defaulting to all lymphoma genes.
#> Joining with `by = join_by(sample_id)`
#> annotating hotspots
#> Joining with `by = join_by(sample_id)`
#> CREBBPHOTSPOT
#> OK
#> FOXO1HOTSPOT
#> OK
#> MYD88HOTSPOT
#> OK
head(coding_tabulated_df[,c(1:10)])
#>                   sample_id HIST1H2BC CDKN2A KMT2D STAT6 SMEK1 CREBBP ZNF423
#> 1                     Akata         0      1     0     0     0      0      0
#> 2                       BL2         0      0     0     0     0      0      0
#> 3                      BL30         0      0     1     0     0      0      0
#> 4                      BL41         0      0     0     0     0      0      0
#> 5                      BL70         0      0     0     0     0      0      0
#> 6 BLGSP-71-06-00001-01A-11D         0      0     0     0     0      0      0
#>   IRF8 BCL2
#> 1    0    0
#> 2    0    0
#> 3    0    0
#> 4    0    0
#> 5    0    0
#> 6    0    0
if (FALSE) { # \dontrun{
#this example would fail because hg38 is not supported by this function (yet)
coding_tabulated_df = get_coding_ssm_status(maf_data=
                        get_coding_ssm(projection = "hg38"))
# Error in get_coding_ssm_status(maf_data = get_coding_ssm(projection = "hg38")) : 
# Currently only grch37 projection (hg19 genome build) is supported.
} # }
```
