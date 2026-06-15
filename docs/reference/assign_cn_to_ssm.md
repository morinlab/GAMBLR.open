# Assign CN to SSM.

Annotate mutations with their copy number information.

## Usage

``` r
assign_cn_to_ssm(
  these_samples_metadata,
  maf_data,
  seg_data,
  projection,
  coding_only = FALSE,
  assume_diploid = FALSE,
  include_silent = FALSE,
  ...
)
```

## Arguments

- these_samples_metadata:

  Metadata table with one or more rows to specify the samples to
  process.

- maf_data:

  A data frame of mutations in MAF format or maf_data object (e.g. from
  `get_coding_ssm` or `get_ssm_by_sample`).

- seg_data:

  A data frame of segmented copy number data or seg_data object

- projection:

  Specified genome projection that returned data is relative to. This is
  only required when it cannot be inferred from maf_df or seg_df (or
  they are not provided).

- coding_only:

  Optional. Set to TRUE to restrict to only variants in coding space
  Default is to work with genome-wide variants.

- assume_diploid:

  Optional, this parameter annotates every mutation as copy neutral.
  Default is FALSE.

- include_silent:

  Logical parameter indicating whether to include silent mutations in
  coding space. Default is FALSE. This parameter only makes sense if
  `coding_only` is set to TRUE.

- ...:

  Any additional parameters.

## Value

A list containing a data frame (MAF-like format) with three extra
columns: - log.ratio is the log ratio from the seg file (NA when no
overlap). - LOH - CN (the rounded absolute copy number estimate of the
region based on log.ratio, NA when no overlap was found).

## Details

This function takes a metadata table and returns all mutations for the
samples in that metadata. Each mutation is annotated with the local copy
number state of each mutated site. The user can specify if only coding
mutations are of interest. To do so, set `coding_only = TRUE`. When
necessary, this function relies on `get_ssm_by_samples` and
`get_cn_segments` to obtain the required data.

## Examples

``` r
# long-handed way
# 1. get some metadata for a collection of samples
some_meta = get_gambl_metadata() %>%
        dplyr::filter(study=="FL_Dreval",
        grepl("SP",sample_id))
# 2. Get the SSMs for these samples

ssm_genomes_grch37 = get_coding_ssm(projection = "grch37",
                                  these_samples_metadata = some_meta)
#> after linking with metadata, we have mutations from 182 samples
# peek at the results
ssm_genomes_grch37 %>% dplyr::select(1:8)
#> genomic_data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position
#> 1       FBXO11              0      .     GRCh37          2       48050367
#> 2       ZNF608              0      .     GRCh37          5      124079933
#> 3       ZNF608              0      .     GRCh37          5      124079935
#> 4     HIST1H1E              0      .     GRCh37          6       26156823
#> 5         PIM1              0      .     GRCh37          6       37138950
#> 6        CCND3              0      .     GRCh37          6       41903710
#> 7         SGK1              0      .     GRCh37          6      134491506
#> 8          MYC              0      .     GRCh37          8      128750681
#> 9          MYC              0      .     GRCh37          8      128750878
#> 10         MYC              0      .     GRCh37          8      128750962
#>    End_Position Strand
#> 1      48050368      +
#> 2     124079933      +
#> 3     124079935      +
#> 4      26156823      +
#> 5      37138950      +
#> 6      41903710      +
#> 7     134491506      +
#> 8     128750681      +
#> 9     128750878      +
#> 10    128750962      +

# 3. Lazily let this function obtain the corresponding seg_data
# for the right genome_build
cn_list = assign_cn_to_ssm(some_meta,ssm_genomes_grch37)
#> Using the bundled CN segments (.seg) calls in GAMBLR.data...
#> Running in default mode of any...

cn_list$maf %>% dplyr::select(1:8,log.ratio,CN)
#> genomic_data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position
#> 1       FBXO11              0      .     GRCh37          2       48050367
#> 2       ZNF608              0      .     GRCh37          5      124079933
#> 3       ZNF608              0      .     GRCh37          5      124079935
#> 4     HIST1H1E              0      .     GRCh37          6       26156823
#> 5         PIM1              0      .     GRCh37          6       37138950
#> 6        CCND3              0      .     GRCh37          6       41903710
#> 7         SGK1              0      .     GRCh37          6      134491506
#> 8          MYC              0      .     GRCh37          8      128750681
#> 9          MYC              0      .     GRCh37          8      128750878
#> 10         MYC              0      .     GRCh37          8      128750962
#>    End_Position Strand log.ratio CN
#> 1      48050368      + 0.0000000  2
#> 2     124079933      + 0.0000000  2
#> 3     124079935      + 0.0000000  2
#> 4      26156823      + 0.0000000  2
#> 5      37138950      + 0.0000000  2
#> 6      41903710      + 0.0000000  2
#> 7     134491506      + 0.0000000  2
#> 8     128750681      + 0.1008947  2
#> 9     128750878      + 0.1008947  2
#> 10    128750962      + 0.1008947  2
if (FALSE) { # \dontrun{
# This wouldn't work because the hg38 seg_data is not bundled
ssm_genomes_hg38 = get_coding_ssm(projection = "hg38",
                                  these_samples_metadata = some_meta)
cn_list = assign_cn_to_ssm(some_meta,ssm_genomes_hg38)

# Easiest/laziest way:
cn_list = assign_cn_to_ssm(projection = "grch37")


cn_list$maf %>% dplyr::group_by(Tumor_Sample_Barcode,CN) %>%
  dplyr::count()
} # }
```
