# Summarizing mutation status as a feature matrix

## Tutorial: assemble feature matrix

In 2023, the Morin Lab developed a classifier of Follicular Lymphoma
(FL) predictive of histological transformation to more aggressive DLBCL.
This study was published in [Blood
(2023).](https://doi.org/10.1182/blood.2022018719) How was the binary
feature matrix assembled for that machine learning model?

In this quick tour we will show how GAMBLR.data resources and
functionality of GAMBLR.open can help you to generate such matrix.

``` r
# Load packages
library(GAMBLR.open)
library(tidyverse)
```

This tutorial explores how to obtain metadata for the samples, simple
somatic mutations in maf format, and auto-magically transform it into a
binary matrix of features.

To simplify this tutorial, we will only include a small subset of
features (not the whole set as was used in the original paper), but this
example will be able to illustrate the functionality and highlight the
main steps of the process.

### Obtain metadata

In the [previous
tutorial](https://morinlab.github.io/GAMBLR.open/articles/getting_started.md),
we have already explored the function
[`get_gambl_metadata()`](https://morinlab.github.io/GAMBLR.open/reference/get_gambl_metadata.md).
We will use it to retreive the metadata:

``` r
metadata <- get_gambl_metadata() %>%
    filter(
        study == "FL_Dreval"
    )
```

We can now see that the samples from the FL study are in the cohort
`FL_Dreval`. Let’s explore these samples more:

``` r
table(metadata$pathology)
```

    COMFL DLBCL    FL
       21   209   213 

For illustration purposes, let’s take 10 samples from the FL study: 5 FL
and 5 DLBCL:

``` r
# Only filter for the samples from FL study
metadata <- metadata %>%
    filter(pathology %in% c("FL", "DLBCL")) %>%
    group_by(pathology) %>%
    arrange(sample_id) %>%
    slice_tail(n = 5) %>%
    ungroup

# How does our metadata looks like now?
str(metadata)
```

    tibble [10 × 30] (S3: tbl_df/tbl/data.frame)
     $ age_group           : chr [1:10] "Other" "Other" "Other" "Other" ...
     $ bam_available       : logi [1:10] TRUE TRUE TRUE TRUE TRUE TRUE ...
     $ cohort              : chr [1:10] "DLBCL_ICGC" "DLBCL_ICGC" "DLBCL_ICGC" "DLBCL_ICGC" ...
     $ compression         : chr [1:10] "cram" "cram" "cram" "cram" ...
     $ COO_consensus       : chr [1:10] "ABC" "GCB" NA "ABC" ...
     $ DHITsig_consensus   : chr [1:10] "DHITsigNeg" "DHITsigNeg" "NA" "DHITsigNeg" ...
     $ EBV_status_inf      : chr [1:10] NA NA NA NA ...
     $ ffpe_or_frozen      : chr [1:10] "frozen" "frozen" "frozen" "frozen" ...
     $ fl_grade            : chr [1:10] NA NA NA NA ...
     $ genetic_subgroup    : chr [1:10] "dFL" "dFL" "dFL" "dFL" ...
     $ genome_build        : chr [1:10] "hs37d5" "hs37d5" "hs37d5" "hs37d5" ...
     $ hiv_status          : chr [1:10] NA NA NA NA ...
     $ lymphgen            : chr [1:10] "Other" "Other" "BN2" "MCD-COMP" ...
     $ lymphgen_cnv_noA53  : chr [1:10] "Other" "Other" "BN2" "BN2/MCD" ...
     $ lymphgen_no_cnv     : chr [1:10] "MCD" "Other" "BN2" "BN2/MCD" ...
     $ lymphgen_with_cnv   : chr [1:10] "A53" "A53" "BN2" "BN2/MCD" ...
     $ lymphgen_wright     : chr [1:10] NA NA NA NA ...
     $ molecular_BL        : chr [1:10] "non-BL" "non-BL" NA "non-BL" ...
     $ normal_sample_id    : chr [1:10] "SP59410" "SP59446" "SP59450" "SP59454" ...
     $ pairing_status      : chr [1:10] "matched" "matched" "matched" "matched" ...
     $ pathology           : chr [1:10] "DLBCL" "DLBCL" "DLBCL" "DLBCL" ...
     $ pathology_rank      : num [1:10] 19 19 19 19 19 15 15 15 15 15
     $ patient_id          : chr [1:10] "DO27833" "DO27851" "DO27853" "DO27855" ...
     $ reference_PMID      : num [1:10] 37084389 37084389 37084389 37084389 37084389 ...
     $ sample_id           : chr [1:10] "SP59412" "SP59448" "SP59452" "SP59456" ...
     $ seq_type            : chr [1:10] "genome" "genome" "genome" "genome" ...
     $ sex                 : chr [1:10] "M" "F" "M" "F" ...
     $ study               : chr [1:10] "FL_Dreval" "FL_Dreval" "FL_Dreval" "FL_Dreval" ...
     $ time_point          : chr [1:10] "A" "A" "A" "A" ...
     $ Tumor_Sample_Barcode: chr [1:10] "SP59412" "SP59448" "SP59452" "SP59456" ...

Our subsetting worked and we can now proceed with matrix assembly.

### Generate feature matrix

We will now use the metadata from the previous step to assemble the
binary feature matrix.

#### Return the simple somatic mutations

First, lets return the data frame with simple somatic mutations in maf
format and store it in a variable. Technically this step is not strictly
necessary, as each function used below will be able to retreive it for
you when the maf data is not provided, but we will advocate for good
practice here and have the maf data stored in a designated variable. As
we will be utilizing both coding and non-coding mutations, we will take
advantage of returning the variants per sample (without necessarily
restricting to coding-only variants).

``` r
# Obtain simple somatic mutations
maf <- get_ssm_by_samples(
    these_samples_metadata = metadata
)
```

Wow, that was super easy and blazing fast! Did it even work?? Let’s
confirm:

``` r
# How many samples are present in maf?
length(unique(maf$Tumor_Sample_Barcode))
```

    [1] 10

``` r
# Are all these samples the ones we are interested in and requested with metadata?
sort(unique(maf$Tumor_Sample_Barcode)) == sort(metadata$Tumor_Sample_Barcode)
```

     [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

``` r
# What are the mutations in the maf? Are they just coding?
table(maf$Variant_Classification)
```

                   3'Flank                  3'UTR                5'Flank
                        29                     20                    127
                     5'UTR        Frame_Shift_Del        Frame_Shift_Ins
                        61                      7                      1
                       IGR           In_Frame_Del                 Intron
                        46                      2                   1502
         Missense_Mutation      Nonsense_Mutation                    RNA
                        78                      9                     15
                    Silent          Splice_Region            Splice_Site
                        29                      6                      8
    Translation_Start_Site
                         1 

We can see from the above outputs that we got somatic mutations for all
requested samples and the maf contains both coding and non-coding
variants. We can now proceed to the next steps.

#### Coding mutations as features

The feature matrix in the FL study contained coding mutation status at
selected genes denoted as 1 when the mutation was present and 0 when
there was no mutation. In addition to that, mutation hotspots at some
genes were also taken into account. To complicate things more,
annotation of hotspots was performed differently depending on the
specific gene and mutational effect. Specifically, the *CREBBP* missemse
mutations were considered hotspots when they occured in the KAT domain,
but not outside of it. The mutations in *FOXO1* were considered hotspots
when the AA change was at M1.

GAMBLR.data provides a one-stop shop to achieve this level of details
out-of-the-box in a simple and convenient way, and is directly used by
the `get_coding_ssm_status` function. This function has a logical
arguments `include_hotspots` to separate regular mutations from the ones
occurring at hotspots, and `review_hotspots`, which will handle the
specific cases we described above in an automated way. Both of these
arguments are `TRUE` by default, so you do not need to toggle them
separately, but in this example we will specify them explicitly just to
illustrate that this will happen during the function call. To keep the
example clean and concise, we will also annotate the mutation status
only for a few of selected genes.

``` r
# Specify genes to annotate
our_genes <- c(
    "CREBBP", "MYD88", "RRAGC",
    "PIM1", "BCL2", "BCL6"
)

# Generate binary matrix
coding_matrix <- get_coding_ssm_status(
    gene_symbols = our_genes,
    these_samples_metadata = metadata,
    maf_data = maf,
    include_hotspots = TRUE,
    review_hotspots = TRUE
)

coding_matrix
```

       sample_id PIM1 CREBBP RRAGC MYD88 BCL6 BCL2 CREBBPHOTSPOT MYD88HOTSPOT
    1    SP59412    1      1     0     0    0    0             0            0
    2    SP59448    0      0     0     0    0    0             0            0
    3    SP59452    1      0     0     1    1    0             0            0
    4    SP59456    0      0     0     0    0    0             0            0
    5    SP59460    0      0     0     0    1    1             0            1
    6    SP59420    0      1     0     0    0    0             0            0
    7    SP59424    0      0     0     0    0    0             1            0
    8    SP59432    0      0     1     0    0    0             0            0
    9    SP59436    0      1     0     0    0    0             0            0
    10   SP59464    0      0     1     0    0    0             0            0

We can see that in this example there is only one sample with hotspot
mutation in *CREBBP*, SP59424 (annotated as `1` in the column
`CREBBPHOTSPOT`). Let’s sanity check this annotation for illustration
purposes:

``` r
crebbp_hotspot_mutation <- maf %>%
    filter(
        Tumor_Sample_Barcode == "SP59424",
        Hugo_Symbol == "CREBBP"
    ) %>%
    select(Chromosome, Start_Position, Variant_Classification)

crebbp_hotspot_mutation
```

    genomic_data Object
    Genome Build: grch37
    Showing first 10 rows:
      Chromosome Start_Position Variant_Classification
    1         16        3788650      Missense_Mutation

This mutation falls within KAT domain and is a missense variant, so it
is indeed correct to be annotated as hotspot. How can we check it does
fall within KAT domain? This is easy with GAMBLR.data!

``` r
# GAMBLR.data has the regions to be considered as hotspots
hotspot_regions_grch37
```

                   chrom     start       end
    CREBBP            16   3785000   3791000
    EZH2               7 148508764 148506238
    NOTCH1             9 139391455 139391455
    NOTCH2             1 120459150 120459150
    CD79B_trunc       17  62007172  62007172
    CD79B_NONtrunc    17  62006800  62006800

``` r
# Now check that CREBBP mutations
between(
    crebbp_hotspot_mutation$Start_Position,
    hotspot_regions_grch37["CREBBP", "start"],
    hotspot_regions_grch37["CREBBP", "end"]
)
```

    [1] TRUE

We have now generated matrix of coding mutations and SSM hotspots in a
binary format and are ready to proceed to the next step.

#### aSHM mutations as features

The FL study also annotated the non-coding mutations at selected aSHM
targets as features of binary matrix. Again, all the necessarily means
to do it in a simple step are provided by GAMBLR.data. The regions
targeted by aSHM [are
available](https://morinlab.github.io/GAMBLR.data/resources/bundled_data.html#ashm-regions)
in the GAMBLR.data. Since they are always complemented with new regions
once they are identified, the latest and most comprehensive version is
always available by referring to `{projection}_ashm_regions`. However,
this use case is a perfect example to demonstarate how to operate on
versioned iterations of the aSHM target list, as the list has been
updated since the time the study was published and at the time of
publication version 0.2 was used. We can refer to it directly by the
version number:

``` r
regions_bed <- somatic_hypermutation_locations_GRCh37_v0.2

head(regions_bed)
```

      chr_name hg19_start  hg19_end    gene     region regulatory_comment
    1     chr2   96808901  96811913   DUSP2   intron-1           enhancer
    2    chr17   56407732  56410140 TSPOAP1 intergenic           enhancer
    3    chr11  128339774 128345731    ETS1    introns           enhancer
    4    chr11  128388492 128394163    ETS1      TSS-2    active_promoter
    5     chr6   31548325  31550717     LTB   intron-1           enhancer
    6     chr3   32020518  32024930 OSBPL10      TSS-1    active_promoter

We will perform some simple modifications to it to make our experience
better and only select few regions for illustrative purposes:

``` r
our_regions <- c(
    "BCL6-TSS",
    "BCL7A-TSS",
    "RHOH-TSS",
    "ZFP36L1-TSS"
)

regions_bed <- create_bed_data(
    regions_bed,
    fix_names = "concat",
    concat_cols = c("gene","region"),
    sep = "-",
    genome_build = "grch37"
)

regions_bed <- regions_bed %>%
    filter(name %in% our_regions)

regions_bed
```

    genomic_data Object
    Genome Build: grch37
    Showing first 10 rows:
      chrom     start       end        name region regulatory_comment
    1     3 187458526 187464632    BCL6-TSS    TSS               <NA>
    2    12 122456912 122464036   BCL7A-TSS    TSS    poised_promoter
    3    14  69257848  69259739 ZFP36L1-TSS    TSS    active_promoter
    4     4  40193105  40204231    RHOH-TSS    TSS    active_promoter

Now we can see whether or not there are any mutations within these
regions:

``` r
ashm_matrix <- cool_overlaps(
    maf,
    regions_bed,
    columns2 = c("chrom", "start", "end")
) %>%
    group_by(Tumor_Sample_Barcode, name) %>%
    summarize(n = n()) %>%
    pivot_wider(
        id_cols = Tumor_Sample_Barcode,
        names_from = name,
        values_from = n,
        values_fill = 0
    ) %>%
    column_to_rownames("Tumor_Sample_Barcode")
ashm_matrix
```

            BCL6-TSS ZFP36L1-TSS BCL7A-TSS RHOH-TSS
    SP59412        4           1         0        0
    SP59420        1           0         1        0
    SP59424        0           0         0        1
    SP59436        1           0         0        0
    SP59448        6           4        14        3
    SP59452       18           2        12       28
    SP59456        2           6         0        9
    SP59460        1           1         1        2
    SP59464        1           0         0        0

We now calculated the number of mutations in each region for each
sample. The original study used pathology-specific values per region to
convert these counts to binary, but here we will use an arbitrary cutoff
of 5 to binarize these counts. We will also convert rownames to column
so it will be easier for us to unite all matrices into single one at the
later steps:

``` r
ashm_matrix[ashm_matrix <= 5] = 0
ashm_matrix[ashm_matrix > 5] = 1
ashm_matrix <- ashm_matrix %>%
    rownames_to_column("sample_id")

ashm_matrix
```

      sample_id BCL6-TSS ZFP36L1-TSS BCL7A-TSS RHOH-TSS
    1   SP59412        0           0         0        0
    2   SP59420        0           0         0        0
    3   SP59424        0           0         0        0
    4   SP59436        0           0         0        0
    5   SP59448        1           0         1        0
    6   SP59452        1           0         1        1
    7   SP59456        0           1         0        1
    8   SP59460        0           0         0        0
    9   SP59464        0           0         0        0

We have now generated matrix of non-coding mutations in a binary format
and are ready to unite all matrices together.

#### Unite into full feature matrix

We can now combine both coding and non-coding mutations into single
matrix:

``` r
feature_matrix <- left_join(
    coding_matrix,
    ashm_matrix
) %>%
    mutate(
        across(
            where(is.numeric), ~ replace_na(.x, 0)
        )
    )

feature_matrix
```

       sample_id PIM1 CREBBP RRAGC MYD88 BCL6 BCL2 CREBBPHOTSPOT MYD88HOTSPOT
    1    SP59412    1      1     0     0    0    0             0            0
    2    SP59448    0      0     0     0    0    0             0            0
    3    SP59452    1      0     0     1    1    0             0            0
    4    SP59456    0      0     0     0    0    0             0            0
    5    SP59460    0      0     0     0    1    1             0            1
    6    SP59420    0      1     0     0    0    0             0            0
    7    SP59424    0      0     0     0    0    0             1            0
    8    SP59432    0      0     1     0    0    0             0            0
    9    SP59436    0      1     0     0    0    0             0            0
    10   SP59464    0      0     1     0    0    0             0            0
       BCL6-TSS ZFP36L1-TSS BCL7A-TSS RHOH-TSS
    1         0           0         0        0
    2         1           0         1        0
    3         1           0         1        1
    4         0           1         0        1
    5         0           0         0        0
    6         0           0         0        0
    7         0           0         0        0
    8         0           0         0        0
    9         0           0         0        0
    10        0           0         0        0

That’s it!

*Happy GAMBLing!*

      /$$$$$$     /$$$$$$    /$$      /$$   /$$$$$$$    /$$        .:::::::
     /$$__  $$   /$$__  $$  | $$$    /$$$  | $$__  $$  | $$        .::    .::
    | $$  \__/  | $$  \ $$  | $$$$  /$$$$  | $$  \ $$  | $$        .::    .::
    | $$ /$$$$  | $$$$$$$$  | $$ $$/$$ $$  | $$$$$$$   | $$   <-   .: .::
    | $$|_  $$  | $$__  $$  | $$  $$$| $$  | $$__  $$  | $$        .::  .::
    | $$  \ $$  | $$  | $$  | $$\  $ | $$  | $$  \ $$  | $$        .::    .::
    |  $$$$$$/  | $$  | $$  | $$ \/  | $$  | $$$$$$$/  | $$$$$$$$  .::      .::
     \______/   |__/  |__/  |__/     |__/  |_______/   |________/
     ~GENOMIC~~~~~~~~~~~~~OF~~~~~~~~~~~~~~~~~B-CELL~~~~~~~~~~~~~~~~~~IN~~~~~~
     ~~~~~~~~~~~~ANALYSIS~~~~~~MATURE~~~~~~~~~~~~~~~~~~~LYMPHOMAS~~~~~~~~~~R~ 
