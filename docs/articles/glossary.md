# GAMBLR terms and concepts

There are several concepts that often come up when using any of the
GAMBLR packages. In particular, there are a few arguments that are
required by most functions and that have a very specific meaning in the
context of GAMBLR. This guide will help you understand the importance of
these and how they are used across GAMBLR.

### metadata

In it’s minimal form, this is a data frame with a set of 7 required
columns: `patient_id`, `Tumor_Sample_Barcode`, `sample_id`, `seq_type`,
`sex`, `cohort`, and `pathology`. In reality many of these variables
such as patient details (e.g. `sex`, `pathology`), convenience variables
such as `cohort` or `study` can contain NA values but those columns must
be present in the metadata. The main purpose of this data frame is to
provide a structure for the metadata that is used by numerous GAMBLR
functions. In general, it provides linkage between unique sample
identifiers and basic metadata fields that are needed by various
functions. Notably, the columns `Tumor_Sample_Barcode` and `sample_id`
are expected to be identical. In most cases, `sample_id` is used but the
alias `Tumor_Sample_Barcode` is required to be present mostly for ease
of linking MAF-type data with the metadata.

Here is an example from GAMBLR.data:

| sample_id | Tumor_Sample_Barcode | patient_id | seq_type | pathology | cohort | sex |
|:---|:---|:---|:---|:---|:---|:---|
| Reddy_2677T | Reddy_2677T | Reddy_2677 | capture | DLBCL | dlbcl_reddy | F |
| Reddy_3540T | Reddy_3540T | Reddy_3540 | capture | DLBCL | dlbcl_reddy | M |
| DLBCL10995T | DLBCL10995T | DLBCL10995 | capture | DLBCL | dlbcl_schmitz | NA |
| Reddy_3742T | Reddy_3742T | Reddy_3742 | capture | DLBCL | dlbcl_reddy | M |
| Reddy_2899T | Reddy_2899T | Reddy_2899 | capture | DLBCL | dlbcl_reddy | F |
| Reddy_3781T | Reddy_3781T | Reddy_3781 | capture | DLBCL | dlbcl_reddy | M |
| 10-36955_tumorC | 10-36955_tumorC | 10-36955 | genome | FL | DLBCL_LSARP_Trios | M |
| Reddy_2781T | Reddy_2781T | Reddy_2781 | capture | DLBCL | dlbcl_reddy | F |
| Reddy_3418T | Reddy_3418T | Reddy_3418 | capture | DLBCL | dlbcl_reddy | F |
| DO221129 | DO221129 | DO221129 | genome | DLBCL | NA | F |
| DLBCL10987T | DLBCL10987T | DLBCL10987 | capture | DLBCL | dlbcl_schmitz | NA |
| Reddy_3948T | Reddy_3948T | Reddy_3948 | capture | DLBCL | dlbcl_reddy | F |
| Reddy_2540T | Reddy_2540T | Reddy_2540 | capture | DLBCL | dlbcl_reddy | M |
| DLBCL-MC_F502_RFK-Tumor | DLBCL-MC_F502_RFK-Tumor | DLBCL-MC_F502_RFK | capture | DLBCL | dlbcl_chapuy | NA |
| Reddy_2232T | Reddy_2232T | Reddy_2232 | capture | DLBCL | dlbcl_reddy | F |

If you are analyzing your own data in GAMBLR and don’t have detailed
metadata, the bare minimum you could get away with is a data frame with
dummy values for every required column. Here’s an example that assumes
you are working entirely with exome (`capture`) data from a collection
of Burkitt lymphoma (BL) patients.

``` r
actual_ids = c("my_sample1","my_sample2","my_sample3","my_sample_1000")

metadata = data.frame(sample_id=actual_ids,
                      Tumor_Sample_Barcode=actual_ids,
                      patient_id = actual_ids,
                      seq_type = "capture",
                      pathology = "BL",
                      cohort = NA,
                      sex=NA)
metadata %>% kableExtra::kable()
```

| sample_id | Tumor_Sample_Barcode | patient_id | seq_type | pathology | cohort | sex |
|:---|:---|:---|:---|:---|:---|:---|
| my_sample1 | my_sample1 | my_sample1 | capture | BL | NA | NA |
| my_sample2 | my_sample2 | my_sample2 | capture | BL | NA | NA |
| my_sample3 | my_sample3 | my_sample3 | capture | BL | NA | NA |
| my_sample_1000 | my_sample_1000 | my_sample_1000 | capture | BL | NA | NA |

### genome build

GAMBLR functions and the data in GAMBLR.data all support two main
genomic coordinate systems. As GAMBL itself and the data bundled and
available through this package represents a large collection of samples
sequenced both locally and externally, we opted to support only two
genome build flavours. For somewhat nuanced reasons, coordinates using
`grch37` or `hg38` are supported. `grch37` uses the same numeric
coordinate system a NCBI `hg19`, but for the former, chromosome names
lack the “chr” prefix. In contrast, the `hg38` genome build is always
chr-prefixed and is in the same coordinate system as the `GRCh38` genome
build. If you are analyzing your own data in GAMBLR, obviously you will
have to know which coordinate system and genome build is appropriate.
Generally, you can make any results that correspond to an unsupported
genome build (e.g. `hg19` or `GRCh38`) by removing or adding the “chr”
prefix to the chromosome name.

``` r
# example of grch37 coordinates
head(GAMBLR.data::chromosome_arms_grch37)
```

      chromosome     start       end arm
    1          1     10000 121500000   p
    2          1 142600000 249250621   q
    3          2     10000  90500000   p
    4          2  96800000 243199373   q
    5          3     10000  87900000   p
    6          3  98300000 198022430   q

``` r
# example of hg38 coordinates
head(GAMBLR.data::chromosome_arms_hg38)
```

      chromosome     start       end arm
    1       chr1     10000 121700000   p
    2       chr1 143200000 248956422   q
    3       chr2         0  91800000   p
    4       chr2  96000000 242193529   q
    5       chr3         0  87800000   p
    6       chr3  98600000 198295559   q

### projection

You can think of `projection` as an alias for genome build that is used
in a very specific (and common) context within GAMBLR packages. When you
use a function that retrieves some sort of data on a coordinate system,
you must tell that function which coordinate system you want. Within
GAMBLR.data, there is data available relative to `grch37` or `hg38`.
Behind the scenes, results have all been *lifted over* or *projected* to
both builds using tools such as `liftOver` or `crossMap`. Hence, we
refer to the results relative to any given genome build as a
`projection`.
