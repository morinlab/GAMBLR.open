# Get Manta SVs

Retrieve Manta SVs for one or many samples

## Usage

``` r
get_manta_sv(
  these_samples_metadata = NULL,
  projection = "grch37",
  region,
  min_vaf = 0.1,
  min_score = 40,
  pass_filters = TRUE,
  verbose = FALSE,
  chromosome,
  qstart,
  qend,
  pairing_status,
  these_sample_ids = NULL,
  ...
)
```

## Arguments

- these_samples_metadata:

  A metadata data frame to limit the result to sample_ids within it

- projection:

  The projection genome build. Default is grch37.

- region:

  Specify a single region to fetch SVs anchored within using the format
  "chrom:start-end"

- min_vaf:

  The minimum tumour VAF for a SV to be returned. Default is 0.1.

- min_score:

  The lowest Manta somatic score for a SV to be returned. Default is 40.

- pass_filters:

  If TRUE (default) only return SVs that are annotated with PASS in the
  FILTER column. Set to FALSE to keep all variants, regardless if they
  PASS the filters.

- verbose:

  Set to FALSE to minimize the output to console. Default is TRUE. This
  parameter also dictates the verbose-ness of any helper function
  internally called inside the main function.

- chromosome:

  DEPRECATED. Use `region` instead.

- qstart:

  DEPRECATED. Use `region` instead.

- qend:

  DEPRECATED. Use `region` instead.

- pairing_status:

  DEPRECATED.

- these_sample_ids:

  DEPRECATED. Subset your metadata and supply \`these_samples_metadata“
  instead.

- ...:

  Any additional parameters.

## Details

Retrieve Manta SVs with additional VCF information to allow for
filtering of high-confidence variants. To get SV calls for multiple
samples, supply a metadata table via `these_samples_metadata` that has
been subset to only those samples. The results will be restricted to the
sample_ids within that data frame. This function can also restrict the
returned breakpoints within a genomic region specified via `region` (in
chr:start-end format). Useful filtering parameters are also available,
use `min_vaf` to set the minimum tumour VAF for a SV to be returned and
`min_score` to set the lowest Manta somatic score for a SV to be
returned. In addition, the user can chose to return all variants, even
the ones not passing the filter criteria. To do so, set
`pass_filters = FALSE` (defaults to TRUE).

## Examples

``` r
# lazily get every SV in the table with default quality filters
all_sv <- get_manta_sv()
#> Using the bundled Manta SV (.bedpe) calls in GAMBLR.data...
head(all_sv)
#> genomic_data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>   CHROM_A   START_A     END_A CHROM_B   START_B     END_B
#> 1       1 161658631 161658631       3  16509907  16509907
#> 2       1 161663959 161663959       9  37363320  37363320
#> 3       1 161663959 161663959       9  37363320  37363320
#> 4      11  65267283  65267283      14 106110907 106110907
#> 5      11  65267422  65267422      14 106110905 106110905
#> 6      14 106029118 106029118      18  60994852  60994852
#>                         manta_name SCORE STRAND_A STRAND_B tumour_sample_id
#> 1         MantaBND:21171:0:1:0:0:0   133        +        +         FL2002T1
#> 2        MantaBND:206628:0:1:0:0:0   122        +        +  09-15842_tumorA
#> 3        MantaBND:195941:0:1:0:0:0   151        +        +  09-15842_tumorB
#> 4      MantaBND:152220:0:1:0:0:0:0    88        +        -        15-38154T
#> 5      MantaBND:152220:0:1:0:0:0:0   135        -        +        15-38154T
#> 6 MantaBND:7:134969:682549:0:1:0:0   112        +        -  19-16466_tumorA
#>   normal_sample_id VAF_tumour  DP pair_status FILTER
#> 1          FL2002N      0.331 127     matched   PASS
#> 2  09-15842_normal      0.281 196     matched   PASS
#> 3  09-15842_normal      0.364 187     matched   PASS
#> 4        15-38154N      0.150 167     matched   PASS
#> 5        15-38154N      0.290 169     matched   PASS
#> 6  19-16466_normal      0.247  77     matched   PASS

# get all SVs for just one cohort
cohort_meta = suppressMessages(get_gambl_metadata()) %>% 
              dplyr::filter(cohort == "DLBCL_cell_lines")

some_sv <- get_manta_sv(these_samples_metadata = cohort_meta, verbose=FALSE)
head(some_sv)
#> genomic_data Object
#> Genome Build: grch37 
#> Showing first 10 rows:
#>   CHROM_A   START_A     END_A CHROM_B   START_B     END_B
#> 1      14 106329465 106329465      18  60793497  60793497
#> 2       3  72551428  72551428       3  72551541  72551541
#> 3       3 189255439 189255439       3 189255440 189255440
#> 4       8 128748200 128748200      14 106114286 106114286
#>                    manta_name SCORE STRAND_A STRAND_B tumour_sample_id
#> 1   MantaBND:194451:1:2:0:0:0   103        +        -           DOHH-2
#> 2    MantaDEL:39707:0:0:0:0:0    45        +        -        SU-DHL-10
#> 3    MantaINS:48107:0:0:0:0:0    64        +        -         SU-DHL-4
#> 4 MantaBND:135279:0:1:0:0:0:0    84        +        -           DOHH-2
#>   normal_sample_id VAF_tumour DP pair_status FILTER
#> 1        14-11247N      0.290 69   unmatched   PASS
#> 2        14-11247N      1.000 16   unmatched   PASS
#> 3        14-11247N      0.686 70   unmatched   PASS
#> 4        14-11247N      0.700 20   unmatched   PASS
nrow(some_sv)
#> [1] 4

# get the SVs in a region around MYC
# WARNING: This is not the best way to find MYC SVs.
# Use annotate_sv on the full SV set instead.
myc_region_hg38 = "chr8:127710883-127761821"
myc_region_grch37 = "8:128723128-128774067"

hg38_myc_locus_sv <- get_manta_sv(region = myc_region_hg38,
                                projection = "hg38",
                                verbose = FALSE)
head(hg38_myc_locus_sv)
#> genomic_data Object
#> Genome Build: hg38 
#> Showing first 10 rows:
#>   CHROM_A   START_A     END_A CHROM_B   START_B     END_B
#> 1    chr8 127716025 127716934   chr14 105862581 105863164
#> 2    chr8 127716523 127716523   chr14 105862757 105862757
#> 3    chr8 127718148 127718148   chr14 105860256 105860256
#> 4    chr8 127718150 127718150   chr14 105860564 105860564
#> 5    chr8 127720983 127720983   chr14 105859144 105859144
#> 6    chr8 127721754 127721755   chr14 105863234 105863235
#>                   manta_name SCORE STRAND_A STRAND_B          tumour_sample_id
#> 1 MantaBND:85417:0:1:0:1:0:0   112        +        - BLGSP-71-06-00084-01A-01D
#> 2 MantaBND:85417:0:1:0:1:0:0   173        -        + BLGSP-71-06-00084-01A-01D
#> 3 MantaBND:90480:0:1:0:0:0:0   152        +        - BLGSP-71-08-00036-01A-01D
#> 4 MantaBND:90480:0:1:0:0:0:0   163        -        + BLGSP-71-08-00036-01A-01D
#> 5  MantaBND:640107:0:1:2:1:0   129        +        - BLGSP-71-23-00408-01A-01E
#> 6   MantaBND:85092:0:1:0:0:0   112        -        + BLGSP-71-08-00023-01A-01D
#>            normal_sample_id VAF_tumour  DP pair_status FILTER
#> 1 BLGSP-71-06-00084-99A-01D      0.373  83     matched   PASS
#> 2 BLGSP-71-06-00084-99A-01D      0.373 153     matched   PASS
#> 3 BLGSP-71-08-00036-10A-01D      0.390 172     matched   PASS
#> 4 BLGSP-71-08-00036-10A-01D      0.325 237     matched   PASS
#> 5 BLGSP-71-06-00286-99A-01D      0.320 100   unmatched   PASS
#> 6 BLGSP-71-08-00023-12A-01D      0.241 199     matched   PASS
nrow(hg38_myc_locus_sv)
#> [1] 200

incorrect_myc_locus_sv <- get_manta_sv(region = myc_region_grch37,
                                projection = "hg38",
                                verbose = FALSE)
head(incorrect_myc_locus_sv)
#> genomic_data Object
#> Genome Build: hg38 
#> Showing first 10 rows:
#>  [1] CHROM_A          START_A          END_A            CHROM_B         
#>  [5] START_B          END_B            manta_name       SCORE           
#>  [9] STRAND_A         STRAND_B         tumour_sample_id normal_sample_id
#> [13] VAF_tumour       DP               pair_status      FILTER          
#> <0 rows> (or 0-length row.names)
nrow(incorrect_myc_locus_sv)
#> [1] 0
# The effect of specifying the wrong coordinate is evident

# Despite potentially being incomplete, we can nonetheless
# annotate these directly for more details
annotated_myc_hg38 = suppressMessages(
         annotate_sv(hg38_myc_locus_sv, genome_build = "hg38")
)
head(annotated_myc_hg38)
#>    chrom1    start1      end1 chrom2    start2      end2   name score strand1
#>    <char>     <num>     <num> <char>     <num>     <num> <char> <num>  <char>
#> 1:      8 127716025 127716934     14 105862581 105863164      .   112       +
#> 2:      8 127716523 127716523     14 105862757 105862757      .   173       -
#> 3:      8 127718148 127718148     14 105860256 105860256      .   152       +
#> 4:      8 127718150 127718150     14 105860564 105860564      .   163       -
#> 5:      8 127720983 127720983     14 105859144 105859144      .   129       +
#> 6:      8 127721754 127721755     14 105863234 105863235      .   112       -
#>    strand2          tumour_sample_id   gene partner  fusion
#>     <char>                    <char> <char>  <char>  <char>
#> 1:       - BLGSP-71-06-00084-01A-01D    MYC     IGH IGH-MYC
#> 2:       + BLGSP-71-06-00084-01A-01D    MYC     IGH IGH-MYC
#> 3:       - BLGSP-71-08-00036-01A-01D    MYC     IGH IGH-MYC
#> 4:       + BLGSP-71-08-00036-01A-01D    MYC     IGH IGH-MYC
#> 5:       - BLGSP-71-23-00408-01A-01E    MYC     IGH IGH-MYC
#> 6:       + BLGSP-71-08-00023-01A-01D    MYC     IGH IGH-MYC
table(annotated_myc_hg38$partner)
#> 
#> IGH 
#> 198 
# The usual MYC partners are seen here
```
