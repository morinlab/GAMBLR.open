# Get ASHM Count Matrix.

Prepare a matrix with one row per sample and one column per region using
a set of hypermutated regions.

## Usage

``` r
get_ashm_count_matrix(
  regions_bed,
  these_samples_metadata,
  this_seq_type,
  projection = "grch37"
)
```

## Arguments

- regions_bed:

  A bed file with one row for each region.

- these_samples_metadata:

  This is used to complete your matrix. All GAMBL samples will be used
  by default. Provide a data frame with at least sample_id for all
  samples if you are using non-GAMBL data.

- this_seq_type:

  The seq type to return results for. Only used if no metadata is
  provided with these_samples_metadata.

- projection:

  Which genome build to use for the mutations (must match the coordinate
  system your regions to avoid a nonsense result)

## Value

matrix

## Details

Values are the number of mutations in that patient in the region.

## Examples

``` r
regions_bed = create_bed_data(GAMBLR.data::grch37_ashm_regions,
                              fix_names="concat",
                              concat_cols=c("gene","region"),
                              sep="-")
my_meta = get_gambl_metadata() %>% dplyr::filter(pathology=="DLBCL")
shm_matrix <- get_ashm_count_matrix(
     regions_bed = regions_bed,
     this_seq_type = "genome",
     these_samples_metadata = my_meta
)
head(shm_matrix[,c(1:12)])
#>           AICDA-TSS BACH2-TSS BCL11A-TSS BCL2-TSS BCL2-intron BCL6-Intergenic-1
#> 02-13135T         0         1          0       20           0                 0
#> 02-20170T         1         1          0        0           0                 0
#> 02-22991T         0         0          0        0           0                 0
#> 04-24937T         0         2         15        4           0                 0
#> 04-28140T         0         0          0        0           0                 0
#> 04-29264T         0         0          0        0           0                 0
#>           BCL6-Intergenic-2 BCL6-Intergenic-3 BCL6-Intergenic-4
#> 02-13135T                 0                 1                 0
#> 02-20170T                 0                 2                 0
#> 02-22991T                 0                 0                 1
#> 04-24937T                 0                 1                 0
#> 04-28140T                 0                 0                 0
#> 04-29264T                 0                 2                 2
#>           BCL6-Intergenic-5 BCL6-TSS BCL7A-TSS
#> 02-13135T                 0        4         3
#> 02-20170T                 0        2         0
#> 02-22991T                 0        1         0
#> 04-24937T                 0       12         2
#> 04-28140T                 0        0         0
#> 04-29264T                 1        3         3
if (FALSE) { # \dontrun{
#this example should fail because the regions_bed is not hg38
shm_matrix <- get_ashm_count_matrix(regions_bed=regions_bed,
                            this_seq_type = "genome",
                            these_samples_metadata = my_meta,
                            projection = "hg38")
# Error in get_ashm_count_matrix(
# Your projection argument does not match the genome_build of regions_bed
} # }
# format the name column to include the coordinates instead of the gene
regions_bed = create_bed_data(GAMBLR.data::hg38_ashm_regions,
                           fix_names="concat",
                           concat_cols=c("chr_name","hg38_start","hg38_end"),
                           sep="-")

 matrix_hg38 <- get_ashm_count_matrix(regions_bed=regions_bed,
                                      this_seq_type = "genome",
                                      these_samples_metadata = my_meta,
                                      projection = "hg38")
print(dim(matrix_hg38))
#> [1] 2312  126
print(head(matrix_hg38[,c(1:8)]))
#>           chr1-203305570-203306650 chr1-226732862-226740184
#> 02-13135T                        0                        0
#> 02-20170T                        0                        0
#> 02-22991T                        0                        0
#> 04-24937T                       21                        0
#> 04-28140T                        0                        0
#> 04-29264T                        1                        0
#>           chr1-226733387-226740281 chr1-28506039-28509827
#> 02-13135T                        0                      0
#> 02-20170T                        0                      0
#> 02-22991T                        0                      0
#> 04-24937T                        0                      0
#> 04-28140T                        0                      0
#> 04-29264T                        0                      0
#>           chr1-30756165-30759164 chr11-102317439-102319346
#> 02-13135T                      0                         0
#> 02-20170T                      0                         1
#> 02-22991T                      0                         0
#> 04-24937T                      0                         4
#> 04-28140T                      0                         0
#> 04-29264T                      0                         0
#>           chr11-111377353-111379499 chr11-118883749-118885805
#> 02-13135T                         0                         0
#> 02-20170T                         0                         0
#> 02-22991T                         0                         0
#> 04-24937T                         1                         0
#> 04-28140T                         0                         0
#> 04-29264T                         1                         0
```
