# Genetic classification of tumours

In another tutorial, we saw at how GAMBLR.open can be used to generate
binary matrices representing the presence/absence of mutations and other
genetic features (e.g. common CNVs, SVs). In this tutorial, we will
review how the binary matrix generation can be utilized to classify
samples according to different classification systems and tumour labels
can be inferred to assign tumours into genetic subgroups.

> **Tip**
>
> We would not need to assemble the matrix separately - the `classify_`
> family of GAMBLR functions will directly support the maf/seg/bedpe
> data and handle matrix construction for you.

The GAMBLR contains [a collection of trained
models](https://github.com/morinlab/GAMBLR.predict) and functions to
pre-format inputs for these models. It contains classifiers of
[Burkitt](https://doi.org/10.1182/blood.2022016534) and
[Follicular](https://doi.org/10.1182/blood.2022018719) lymphomas
originally published, as well as reproduced of DLBCL classification by
the groupings of [Chapuy et
al](https://doi.org/10.1038/s41591-018-0016-8), [Lacy et
al](https://doi.org/10.1182/blood.2019003535), and [Runge et
al](https://doi.org/10.1111/bjh.17132).

## Load packages and data

First, we will load required packages. Only `GAMBLR.open` and
`tidyverse` are needed - they both will load for you all required
functionality.

``` r
# Load packages
library(GAMBLR.open)
library(tidyverse)
```

Next, we will obtain the data. We would need the metadata, SSM in maf
format, CNV in seg format, and SV in bedped format. For the
demonstration purposes, we will use the data bundled with GAMBLR to show
the classification functionality. Here, we will subset data to a small
number of samples and will illustrate the required formatting and
minimal required information.

### Metadata

The metadata is a data frame which contains the column `sample_id`
listing sample ids for tumours to be classified, and column `pathology`
which will dictate sliding threshold for aSHM site annotation when
necessary. The column `seq_type` is not required for classification
purposes, but will be kept in this tutorial for the purpose of
subsequent retreival of mutation data.

``` r
# Load metadata
metadata <- get_gambl_metadata() %>%
    filter(seq_type == "genome") %>%
    filter(pathology %in% c("FL", "DLBCL", "BL")) %>%
    group_by(sample_id) %>%
    slice_head() %>%
    ungroup %>%
    filter(!study %in% c("DLBCL_Arthur", "DLBCL_Thomas"))

metadata %>%
    count(pathology)
```

| pathology |   n |
|:----------|----:|
| BL        | 234 |
| DLBCL     | 334 |
| FL        | 219 |

``` r
# Demonstrate the required columns in metadata
metadata <- metadata %>%
    select(sample_id, pathology, seq_type) 

head(metadata) 
```

| sample_id       | pathology | seq_type |
|:----------------|:----------|:---------|
| 00-14595_tumorC | DLBCL     | genome   |
| 00-14595_tumorD | FL        | genome   |
| 00-15201_tumorA | DLBCL     | genome   |
| 00-15201_tumorB | DLBCL     | genome   |
| 00-26427_tumorA | DLBCL     | genome   |
| 00-26427_tumorC | DLBCL     | genome   |

### SSM (maf format)

The key data to be provided for the genetic subgroup classification is
SSM in maf format. Only a small subset of columns is required to be used
as input information, and here we would demonstrate the expected format
and minimal required information for SSM.

``` r
maf <- get_ssm_by_samples(
    these_samples_metadata = metadata
) %>% as.data.frame

# Only these columns are required
maf_columns <- c(
    "Hugo_Symbol", "NCBI_Build",
    "Chromosome", "Start_Position", "End_Position",
    "Variant_Classification", "HGVSp_Short",
    "Tumor_Sample_Barcode"
)

maf <- maf %>%
    select(all_of(maf_columns))

head(maf) 
```

| Hugo_Symbol | NCBI_Build | Chromosome | Start_Position | End_Position | Variant_Classification | HGVSp_Short | Tumor_Sample_Barcode |
|:---|:---|:---|---:|---:|:---|:---|:---|
| BTBD3 | GRCh37 | 20 | 11867899 | 11867899 | 5'Flank | NA | FL1011T1 |
| MTOR | GRCh37 | 1 | 11170344 | 11170344 | Intron | NA | FL1011T1 |
| SPEN | GRCh37 | 1 | 16228491 | 16228491 | Intron | NA | FL1011T1 |
| JAK1 | GRCh37 | 1 | 65416629 | 65416629 | Intron | NA | FL1011T1 |
| NOTCH2 | GRCh37 | 1 | 120460767 | 120460767 | Intron | NA | FL1011T1 |
| DCAF6 | GRCh37 | 1 | 167981016 | 167981016 | Intron | NA | FL1011T1 |

### CNV (seg format)

Several classifiers ([Chapuy et
al](https://doi.org/10.1038/s41591-018-0016-8), [Lacy et
al](https://doi.org/10.1182/blood.2019003535), and [Runge et
al](https://doi.org/10.1111/bjh.17132)) require a CNV information to be
provided, and here we would demonstrate the expected format and minimal
required information for CNV in standard seg format.

``` r
seg <- get_cn_segments(
    these_samples_metadata = metadata
) %>% as.data.frame

# Only these columns are required
seg_columns <- c(
    "ID", "chrom", "start", "end", "LOH_flag", "log.ratio"
)

seg <- seg %>%
    select(all_of(seg_columns))

head(seg) 
```

| ID        | chrom |     start |       end | LOH_flag | log.ratio |
|:----------|:------|----------:|----------:|---------:|----------:|
| 02-13135T | 1     |     10001 |    762600 |        0 | 0.0000000 |
| 02-13135T | 1     |    762601 | 121500000 |        0 | 0.0000000 |
| 02-13135T | 1     | 142600000 | 161506889 |        0 | 0.0000000 |
| 02-13135T | 1     | 161506890 | 161652716 |        0 | 0.0000000 |
| 02-13135T | 1     | 161652717 | 162110568 |        0 | 0.7279926 |
| 02-13135T | 1     | 162110569 | 162111399 |        0 | 0.0000000 |

### SV (bedpe format)

Several classifiers require SV data to be provided, and here we would
demonstrate the expected format and minimal required information for SV
in standard bed format.

``` r
bedpe <- get_manta_sv(
    these_samples_metadata = metadata
) %>% as.data.frame

# Only these columns are required
bedpe_columns <- colnames(bedpe)[1:11]

bedpe <- bedpe %>%
    select(all_of(bedpe_columns))

head(bedpe) 
```

| CHROM_A | START_A | END_A | CHROM_B | START_B | END_B | manta_name | SCORE | STRAND_A | STRAND_B | tumour_sample_id |
|:---|---:|---:|:---|---:|---:|:---|---:|:---|:---|:---|
| 1 | 161658631 | 161658631 | 3 | 16509907 | 16509907 | MantaBND:21171:0:1:0:0:0 | 133 | \+ | \+ | FL2002T1 |
| 1 | 161663959 | 161663959 | 9 | 37363320 | 37363320 | MantaBND:206628:0:1:0:0:0 | 122 | \+ | \+ | 09-15842_tumorA |
| 1 | 161663959 | 161663959 | 9 | 37363320 | 37363320 | MantaBND:195941:0:1:0:0:0 | 151 | \+ | \+ | 09-15842_tumorB |
| 11 | 65267283 | 65267283 | 14 | 106110907 | 106110907 | MantaBND:152220:0:1:0:0:0:0 | 88 | \+ | \- | 15-38154T |
| 11 | 65267422 | 65267422 | 14 | 106110905 | 106110905 | MantaBND:152220:0:1:0:0:0:0 | 135 | \- | \+ | 15-38154T |
| 14 | 106029118 | 106029118 | 18 | 60994852 | 60994852 | MantaBND:7:134969:682549:0:1:0:0 | 112 | \+ | \- | 19-16466_tumorA |

> **Note**
>
> The column `manta_name` can be any other column from the tool of
> choice that generated bedpe output and reports the unique id of the SV
> event, and is expected to be unique for the SV event.

## Classify DLBCL

Here, we will explore the reproduced DLBCL classification by the
groupings of [Chapuy et al](https://doi.org/10.1038/s41591-018-0016-8),
[Lacy et al](https://doi.org/10.1182/blood.2019003535), and [Runge et
al](https://doi.org/10.1111/bjh.17132). The classification algorithm can
be easily controlled by switching the argument `method` of the
`classify_dlbcl()` function.

> **Important**
>
> The DLBCL classifiers were not released with the publications, and the
> models provided here are the best attempt on reproducing the original
> models. While the model described by Chapuy et al is reproduced with
> \> 92% accuracy, the Lacy and HMRN classifier results need to be taken
> with caution as the accuracy of the bundled model is only 80%.

### Chapuy et al. (C0-C5 clusters)

The paper published originally in 2018 by [Chapuy et
al](https://doi.org/10.1038/s41591-018-0016-8) was the first attempt to
use systematic approach and classify DLBCL patients into genetic
subgroupings using genomic information.

> **Note**
>
> This model classifies tumours according to the original 2018
> publication, and is not directly related to the 2024 DLBCLass model.

Here is how these subgroupings can be recapitulated with GAMBLR
functionality:

``` r
predictions <- classify_dlbcl(
    these_samples_metadata = metadata %>%
        filter(pathology == "DLBCL"),
    maf_data = maf,
    seg_data = seg,
    sv_data = bedpe
)
```

The `classify_` family of functions by default return both the assembled
matrix and predictions. The model reports the confidence of each
tumour’s subgoup vote, as well as final label.

``` r
# assembled matrix
predictions$matrix[1:5,1:5] 
```

|                 | ACTB | B2M | BCL10 | BCL11A | BCL2 |
|:----------------|-----:|----:|------:|-------:|-----:|
| 00-14595_tumorC |    0 |   0 |     0 |      0 |    2 |
| 00-15201_tumorA |    0 |   0 |     0 |      0 |    0 |
| 00-15201_tumorB |    0 |   0 |     0 |      0 |    0 |
| 00-26427_tumorA |    0 |   0 |     0 |      2 |    0 |
| 00-26427_tumorC |    0 |   0 |     0 |      0 |    0 |

``` r
# subgroup assignment
predictions$predictions %>%
 slice_head(n=10)
```

| sample_id | C0 | C1 | C2 | C3 | C4 | C5 | Chapuy_cluster |
|:---|---:|---:|---:|---:|---:|---:|:---|
| 00-14595_tumorC | 0 | 0.1060348 | 0.0000000 | 0.5110284 | 0.2784030 | 0.1045338 | C3 |
| 00-15201_tumorA | 0 | 0.0931391 | 0.0386996 | 0.1889158 | 0.3859610 | 0.2932845 | C4 |
| 00-15201_tumorB | 0 | 0.2842088 | 0.0914615 | 0.1265154 | 0.0730396 | 0.4247747 | C5 |
| 00-26427_tumorA | 0 | 0.0797927 | 0.0670896 | 0.0457110 | 0.1034915 | 0.7039151 | C5 |
| 00-26427_tumorC | 0 | 0.1152963 | 0.0000000 | 0.0000000 | 0.3449964 | 0.5397073 | C5 |
| 01-14774_tumorA | 0 | 0.3303927 | 0.0349260 | 0.1965204 | 0.3047769 | 0.1333841 | C1 |
| 01-14774_tumorB | 0 | 0.3278894 | 0.0070903 | 0.1219750 | 0.4249661 | 0.1180793 | C4 |
| 01-16433_tumorA | 0 | 0.0000000 | 0.0000000 | 0.8286825 | 0.1713175 | 0.0000000 | C3 |
| 01-16433_tumorB | 0 | 0.2873270 | 0.2551277 | 0.2676392 | 0.1094943 | 0.0804118 | C1 |
| 01-23117_tumorA | 0 | 0.2100450 | 0.0000000 | 0.0463066 | 0.2689441 | 0.4747043 | C5 |

``` r
count(predictions$predictions, Chapuy_cluster) 
```

| Chapuy_cluster |   n |
|:---------------|----:|
| C0             |   2 |
| C1             |  20 |
| C2             |  74 |
| C3             | 112 |
| C4             |  47 |
| C5             |  79 |

### Lacy et al.

The paper published originally in 2020 by [Lacy et
al](https://doi.org/10.1182/blood.2019003535) did not release the
classification algorithm, but we recapitulated it by training random
forest model with 80% accuaracy. Here is how these subgroupings can be
recapitulated with GAMBLR functionality:

``` r
predictions <- classify_dlbcl(
    these_samples_metadata = metadata %>%
        filter(pathology == "DLBCL"),
    maf_data = maf,
    seg_data = seg,
    sv_data = bedpe,
    method = "lacy"
)
```

Similar to the Chapuy et al predictions, we can see both the constructed
matrix and the confidence of each tumour’s subgoup vote, as well as
final label.

``` r
# assembled matrix
predictions$matrix[1:5,1:5]
```

|                 | BCL11A_amp | MALT1_amp | REL_amp | XPO1_amp | CD58_OR_del |
|:----------------|-----------:|----------:|--------:|---------:|------------:|
| 00-14595_tumorC |          0 |         0 |       0 |        0 |           0 |
| 00-15201_tumorA |          0 |         0 |       0 |        0 |           0 |
| 00-15201_tumorB |          0 |         0 |       0 |        0 |           0 |
| 00-26427_tumorA |          0 |         0 |       0 |        0 |           0 |
| 00-26427_tumorC |          0 |         0 |       0 |        0 |           1 |

``` r
# subgroup assignment
head(predictions$predictions) 
```

| sample_id       |  BCL2 | MYD88 | NOTCH2 | Other | SOCS1/SGK1 | TET2/SGK1 | Lacy_cluster |
|:----------------|------:|------:|-------:|------:|-----------:|----------:|:-------------|
| 00-14595_tumorC | 0.278 | 0.070 |  0.058 | 0.010 |      0.500 |     0.084 | SOCS1/SGK1   |
| 00-15201_tumorA | 0.220 | 0.102 |  0.148 | 0.070 |      0.440 |     0.020 | SOCS1/SGK1   |
| 00-15201_tumorB | 0.052 | 0.186 |  0.500 | 0.210 |      0.026 |     0.026 | NOTCH2       |
| 00-26427_tumorA | 0.020 | 0.310 |  0.194 | 0.470 |      0.000 |     0.006 | Other        |
| 00-26427_tumorC | 0.020 | 0.754 |  0.032 | 0.162 |      0.016 |     0.016 | MYD88        |
| 01-14774_tumorA | 0.262 | 0.064 |  0.196 | 0.016 |      0.398 |     0.064 | SOCS1/SGK1   |

``` r
count(predictions$predictions, Lacy_cluster) 
```

| Lacy_cluster |   n |
|:-------------|----:|
| BCL2         | 110 |
| MYD88        |  57 |
| NOTCH2       |  34 |
| Other        |  42 |
| SOCS1/SGK1   |  81 |
| TET2/SGK1    |  10 |

### HMRN

The paper published by [Runge et al](https://doi.org/10.1111/bjh.17132)
modified the original Lacy et al. subgoupings to more closely adhere the
genetic subgroups to the
[LymphGen](https://doi.org/10.1016/j.ccell.2020.03.015) classification
system. The difference between the Runge and Lacy methods is that
tumours with truncating *NOTCH1* mutation will be assigned to a separate
subgrouping category regardless of other genetic alterations present in
the same tumour sample. Here is how these subgroupings can be
recapitulated with GAMBLR functionality:

``` r
predictions <- classify_dlbcl(
    these_samples_metadata = metadata %>%
        filter(pathology == "DLBCL"),
    maf_data = maf,
    seg_data = seg,
    sv_data = bedpe,
    method = "hmrn"
)
```

Take a look at the resulting output

``` r
# assembled matrix
predictions$matrix[1:5,1:5]
```

|                 | BCL11A_amp | MALT1_amp | REL_amp | XPO1_amp | CD58_OR_del |
|:----------------|-----------:|----------:|--------:|---------:|------------:|
| 00-14595_tumorC |          0 |         0 |       0 |        0 |           0 |
| 00-15201_tumorA |          0 |         0 |       0 |        0 |           0 |
| 00-15201_tumorB |          0 |         0 |       0 |        0 |           0 |
| 00-26427_tumorA |          0 |         0 |       0 |        0 |           0 |
| 00-26427_tumorC |          0 |         0 |       0 |        0 |           1 |

``` r
predictions$matrix[1:5,1:5]
```

|                 | BCL11A_amp | MALT1_amp | REL_amp | XPO1_amp | CD58_OR_del |
|:----------------|-----------:|----------:|--------:|---------:|------------:|
| 00-14595_tumorC |          0 |         0 |       0 |        0 |           0 |
| 00-15201_tumorA |          0 |         0 |       0 |        0 |           0 |
| 00-15201_tumorB |          0 |         0 |       0 |        0 |           0 |
| 00-26427_tumorA |          0 |         0 |       0 |        0 |           0 |
| 00-26427_tumorC |          0 |         0 |       0 |        0 |           1 |

``` r
# subgroup assignment
head(predictions$predictions) 
```

| sample_id | BCL2 | MYD88 | NOTCH2 | Other | SOCS1/SGK1 | TET2/SGK1 | NOTCH1 | hmrn_cluster |
|:---|---:|---:|---:|---:|---:|---:|---:|:---|
| 00-14595_tumorC | 0.278 | 0.070 | 0.058 | 0.010 | 0.500 | 0.084 | 0 | SOCS1/SGK1 |
| 00-15201_tumorA | 0.220 | 0.102 | 0.148 | 0.070 | 0.440 | 0.020 | 0 | SOCS1/SGK1 |
| 00-15201_tumorB | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1 | NOTCH1 |
| 00-26427_tumorA | 0.020 | 0.310 | 0.194 | 0.470 | 0.000 | 0.006 | 0 | Other |
| 00-26427_tumorC | 0.020 | 0.754 | 0.032 | 0.162 | 0.016 | 0.016 | 0 | MYD88 |
| 01-14774_tumorA | 0.262 | 0.064 | 0.196 | 0.016 | 0.398 | 0.064 | 0 | SOCS1/SGK1 |

``` r
count(predictions$predictions, hmrn_cluster)
```

| hmrn_cluster |   n |
|:-------------|----:|
| BCL2         | 110 |
| MYD88        |  56 |
| NOTCH1       |   3 |
| NOTCH2       |  32 |
| Other        |  41 |
| SOCS1/SGK1   |  81 |
| TET2/SGK1    |  11 |

## Classify FL

Here, we will explore the classification of FL tumours into genetic
subgroups with differential propensity for transformation as originally
described [here](https://doi.org/10.1182/blood.2022018719). The
developed model can be easily accessed with `classify_fl()` function.

``` r
predictions <- classify_fl(
    these_samples_metadata = metadata %>%
        filter(pathology  %in% c("FL", "DLBCL")),
    maf_data = maf,
    output = "both"
)
```

Similar to DLBCL classifier, we can take a look at the assembled matrix,
as well as at the predictions and the confidence of each tumour’s vote:

``` r
# assembled matrix
predictions$matrix[1:10, 1:5]
```

|                 | ACTB | ARID1A | ATP6AP1 | ATP6V1B2 | B2M |
|:----------------|-----:|-------:|--------:|---------:|----:|
| 00-14595_tumorC |    0 |      1 |       0 |        0 |   0 |
| 00-14595_tumorD |    0 |      0 |       0 |        0 |   0 |
| 00-15201_tumorA |    0 |      0 |       0 |        0 |   0 |
| 00-15201_tumorB |    0 |      1 |       0 |        0 |   0 |
| 00-26427_tumorA |    0 |      0 |       0 |        0 |   0 |
| 00-26427_tumorC |    0 |      0 |       0 |        0 |   0 |
| 01-14774_tumorA |    0 |      0 |       0 |        0 |   1 |
| 01-14774_tumorB |    0 |      0 |       0 |        0 |   0 |
| 01-16433_tumorA |    0 |      0 |       0 |        0 |   0 |
| 01-16433_tumorB |    0 |      1 |       0 |        0 |   1 |

``` r
# subgroup assignment
head(predictions$predictions) 
```

| sample_id       |   cFL |   dFL | is_cFL |
|:----------------|------:|------:|:-------|
| 00-14595_tumorC | 0.052 | 0.948 | dFL    |
| 00-14595_tumorD | 0.038 | 0.962 | dFL    |
| 00-15201_tumorA | 0.000 | 1.000 | dFL    |
| 00-15201_tumorB | 0.078 | 0.922 | dFL    |
| 00-26427_tumorA | 0.018 | 0.982 | dFL    |
| 00-26427_tumorC | 0.014 | 0.986 | dFL    |

``` r
count(predictions$predictions, is_cFL)
```

| is_cFL |   n |
|:-------|----:|
| cFL    |  41 |
| dFL    | 510 |

## Classify BL

Here, we will explore the classification of BL tumours into genetic
subgroups as originally described in [Thomas et
al](https://doi.org/10.1182/blood.2022016534). The developed model
reported in this study can be easily accessed with `classify_bl()`
function.

``` r
predictions <- classify_bl(
    these_samples_metadata = metadata %>%
        filter(pathology  %in% c("BL", "DLBCL")),
    maf_data = maf
)
```

Similar to other classifiers, we can take a look at the assembled
matrix, as well as at the predictions and the confidence of each
tumour’s vote:

``` r
# assembled matrix
predictions$matrix[1:10, 1:5]
```

|                 | ARID1A | B2M | BCL10 | BCL11A | BCL2 |
|:----------------|-------:|----:|------:|-------:|-----:|
| 00-14595_tumorC |      1 |   0 |     0 |      0 |    1 |
| 00-15201_tumorA |      0 |   0 |     0 |      0 |    0 |
| 00-15201_tumorB |      1 |   0 |     0 |      0 |    0 |
| 00-26427_tumorA |      0 |   0 |     0 |      1 |    0 |
| 00-26427_tumorC |      0 |   0 |     0 |      0 |    0 |
| 01-14774_tumorA |      0 |   1 |     1 |      0 |    0 |
| 01-14774_tumorB |      0 |   0 |     0 |      0 |    0 |
| 01-16433_tumorA |      0 |   0 |     0 |      0 |    1 |
| 01-16433_tumorB |      1 |   1 |     0 |      0 |    1 |
| 01-23117_tumorA |      0 |   1 |     0 |      0 |    0 |

``` r
# subgroup assignment
head(predictions$predictions) 
```

| sample_id       | DGG-BL | DLBCL | IC-BL | Q53-BL | BL_subgroup |
|:----------------|-------:|------:|------:|-------:|:------------|
| 00-14595_tumorC |  0.116 | 0.816 | 0.068 |  0.000 | DLBCL       |
| 00-15201_tumorA |  0.052 | 0.930 | 0.002 |  0.016 | DLBCL       |
| 00-15201_tumorB |  0.296 | 0.648 | 0.050 |  0.006 | DLBCL       |
| 00-26427_tumorA |  0.036 | 0.964 | 0.000 |  0.000 | DLBCL       |
| 00-26427_tumorC |  0.000 | 1.000 | 0.000 |  0.000 | DLBCL       |
| 01-14774_tumorA |  0.072 | 0.894 | 0.018 |  0.016 | DLBCL       |

``` r
count(predictions$predictions, BL_subgroup)
```

| BL_subgroup |   n |
|:------------|----:|
| DGG-BL      | 106 |
| DLBCL       | 345 |
| IC-BL       | 108 |
| Q53-BL      |   9 |

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
