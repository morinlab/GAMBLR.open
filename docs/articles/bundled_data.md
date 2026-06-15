# GAMBLR.data

The GAMBLR.data package, which is installed and loaded as part of
GAMBLR.open, comes with many different bundled objects with various
uses. These can be divided in the following categories:

## Somatic variants

- `sample_data` A list of data frames containing the metadata, simple
  somatic, copy number, and structural variants collected together from
  the supplemental tables of large sequencing studies of B-cell
  lymphomas.

> **Important**
>
> The data stored in `sample_data` is not meant to be accessed directly.
> You should use one of the various [`get_`
> functions](https://morinlab.github.io/GAMBLR.open/articles/reference/index.html#get-functions)
> in GAMBLR.open to retrieve the data for your use case.

## Curated gene lists

- `gene_blacklist` A tibble with gene symbols (Hugo) that fall within
  blacklisted regions of the genome. The genes in this data object
  represent common sequencing artifacts and are discarded during the
  data analysis.

``` r
# Look at 15 randomly sampled genes in the blacklist
slice_sample(gene_blacklist,n=15) 
```

| Gene       |
|:-----------|
| RP1L1      |
| KIR2DS2    |
| KRTAP3-2   |
| DMRTB1     |
| RYR1       |
| RP11       |
| KIR2DL1    |
| RYR3       |
| ZFAT-AS1_2 |
| MED15      |
| HSP90B2P   |
| KRTAP5-14P |
| KRTAP20-4  |
| OBSCN-AS1  |
| XIRP2      |

- `lymphoma_genes` A data frame with a manually curated set of genes
  commonly mutated in lymphomas. For each pathology, there is a column
  that indicates whether (TRUE) or not (FALSE) each gene is in either
  the Tier 1 or Tier 2 list for that pathology. For more granularity,
  you can also use the corresponding PATHOLOGY_Tier column
  e.g. DLBCL_Tier.

``` r
# Look at the first 15 DLBCL genes in Tier 1
dplyr::filter(lymphoma_genes,DLBCL_Tier == 1) %>% 
    slice_head(n=15) 
```

| Gene | MCL_Tier | MCL | FL_Tier | FL | DLBCL_Tier | DLBCL | BL_Tier | BL | ensembl_gene_id | hgnc_symbol | LymphGen | Reddy | Chapuy |
|:---|---:|:---|---:|:---|---:|:---|---:|:---|:---|:---|:---|:---|:---|
| ACTB | NA | FALSE | 1 | TRUE | 1 | TRUE | NA | FALSE | ENSG00000075624 | ACTB | TRUE | TRUE | TRUE |
| ACTG1 | NA | FALSE | 2 | TRUE | 1 | TRUE | NA | FALSE | ENSG00000267807 | ACTG1 | TRUE | FALSE | FALSE |
| ARID1A | NA | FALSE | 1 | TRUE | 1 | TRUE | 1 | TRUE | ENSG00000117713 | ARID1A | TRUE | TRUE | FALSE |
| ATM | 1 | TRUE | NA | FALSE | 1 | TRUE | NA | FALSE | ENSG00000149311 | ATM | FALSE | TRUE | FALSE |
| B2M | 2 | TRUE | 1 | TRUE | 1 | TRUE | 2 | TRUE | ENSG00000166710 | B2M | FALSE | TRUE | TRUE |
| BCL10 | NA | FALSE | 1 | TRUE | 1 | TRUE | NA | FALSE | ENSG00000142867 | BCL10 | TRUE | TRUE | TRUE |
| BCL2 | NA | FALSE | 1 | TRUE | 1 | TRUE | 2 | TRUE | ENSG00000171791 | BCL2 | TRUE | TRUE | TRUE |
| BCL6 | NA | FALSE | 1 | TRUE | 1 | TRUE | 1 | TRUE | ENSG00000113916 | BCL6 | TRUE | TRUE | TRUE |
| BCL7A | NA | FALSE | 1 | TRUE | 1 | TRUE | 1 | TRUE | ENSG00000110987 | BCL7A | FALSE | TRUE | FALSE |
| BIRC6 | NA | FALSE | NA | FALSE | 1 | TRUE | NA | FALSE | ENSG00000115760 | BIRC6 | FALSE | TRUE | FALSE |
| BRAF | NA | FALSE | NA | FALSE | 1 | TRUE | NA | FALSE | ENSG00000157764 | BRAF | FALSE | TRUE | TRUE |
| BTG1 | NA | FALSE | 2 | TRUE | 1 | TRUE | 2 | TRUE | ENSG00000133639 | BTG1 | TRUE | TRUE | TRUE |
| BTG2 | NA | FALSE | 2 | TRUE | 1 | TRUE | NA | FALSE | ENSG00000159388 | BTG2 | TRUE | TRUE | TRUE |
| BTK | NA | FALSE | 1 | TRUE | 1 | TRUE | NA | FALSE | ENSG00000010671 | BTK | FALSE | TRUE | FALSE |
| CARD11 | 1 | TRUE | 1 | TRUE | 1 | TRUE | 2 | TRUE | ENSG00000198286 | CARD11 | FALSE | TRUE | TRUE |

- `lymphoma_genes_comprehensive` An older (and outdated) curated list of
  genes reported as significantly mutated in the large lymphoma studies.
  Both Ensembl ID and Hugo Symbol are available as gene identifiers.
  This data contains annotations for the studies by Chapuy, Reddy,
  Wright (LymphGen), Lacy, as well as annotations for whether the gene
  is curated, reported as SMG in other_studies, or a target of aSHM.
  This is mostly included for backwards compatability and you should use
  `lymphoma_genes` instead.

## Coordinate-based resources

- `chromosome_arms_grch37`: A data frame with the chromosome arm
  coordinates with respect to the grch37 projection.
- `chromosome_arms_hg38` A data frame with the chromosome arm
  coordinates with respect to the hg38 projection.
- `grch37_gene_coordinates` A data frame of all gene coordinates with
  respect to grch37. Contains both Ensembl ID and Hugo Symbol as
  identifiers.
- `grch37_lymphoma_genes_bed` A data frame in the bed format for genes
  commonly associated with B-cell lymphomas. Coordinates are with
  respect to grch37.
- `grch37_oncogene` A data frame with the coordinates of lymphoma
  oncogenes relative to grch37. Used in mapping of the breakpoint
  coordinates.
- `grch37_partners` A data frame of translocation partners for oncogenes
  with coordinates relative to grch37.
- `hg38_gene_coordinates` A data frame of all gene coordinates with
  respect to hg38. Contains both Ensembl ID and Hugo Symbol as
  identifiers.
- `hg38_lymphoma_genes_bed` A data frame in the bed format for genes
  commonly associated with B-cell lymphomas. Coordinates are with
  respect to hg38.
- `hg38_oncogene` A data frame with the coordinates of lymphoma
  oncogenes relative to the hg38. Used in mapping of the breakpoint
  coordinates.
- `hg38_partners` A data frame of translocation partners for oncogenes
  with relative coordinates to hg38.
- `grch37_all_gene_coordinates` A data frame of protein-coding gene
  coordinates relative to grch37. Contains both Ensembl ID and Hugo
  Symbol as identifiers. Mainly here for backwards compatibility with
  earlier GAMBLR versions.
- `hotspot_regions_grch37` A data frame of mutation hotspot regions
  relative to grch37.
- `hotspot_regions_hg38` A data frame of mutation hotspot regions
  relative to hg38.
- `target_regions_grch37` A data frame with coordinates of the regions
  of the genome targeted by the whole exome sequencing panel Agilent V5
  (no UTR) relative to grch37.
- `target_regions_hg38` A data frame with coordinates of the regions of
  the genome targeted by the whole exome sequencing panel Agilent V5 (no
  UTR) relative to hg38.

## aSHM regions

- `grch37_ashm_regions` Aberrant somatic hypermutation (aSHM) regions
  relative to grch37. This object always by default refers to the most
  recent version of the aSHM regions.

``` r
# Peek at the first 15 rows
grch37_ashm_regions %>%
    slice_head(n=15) 
```

| chr_name | hg19_start |  hg19_end | gene    | region   | regulatory_comment |
|:---------|-----------:|----------:|:--------|:---------|:-------------------|
| chr1     |    6661482 |   6662702 | KLHL21  | TSS      | NA                 |
| chr1     |   23885584 |  23885835 | ID3     | TSS      | NA                 |
| chr1     |   28832551 |  28836339 | RCC1    | TSS      | NA                 |
| chr1     |   31229012 |  31232011 | LAPTM5  | TSS      | NA                 |
| chr1     |  150550814 | 150552135 | MCL1    | intron   | NA                 |
| chr2     |   88904839 |  88909096 | EIF2Ak3 | intron   | NA                 |
| chr2     |   88925456 |  88927581 | EIF2Ak3 | TSS      | NA                 |
| chr2     |  232572640 | 232574297 | PTMA    | TSS      | NA                 |
| chr1     |  157669490 | 157671299 | FCRL3   | TSS      | NA                 |
| chr1     |  203274698 | 203275778 | BTG2    | intron   | active_promoter    |
| chr1     |  226864857 | 226873452 | ITPKB   | intron   | weak_enhancer      |
| chr1     |  226920563 | 226927885 | ITPKB   | TSS      | active_promoter    |
| chr1     |  226921088 | 226927982 | ITPKB   | intron-1 | enhancer           |
| chr1     |  245023502 | 245029083 | HNRNPU  | TSS      | NA                 |
| chr2     |   60773789 |  60783486 | BCL11A  | TSS      | active_promoter    |

- `hg38_ashm_regions` Aberrant somatic hypermutation (aSHM) regions
  relative to hg38. This object always by default refers to the most
  recent version of the aSHM regions.

``` r
# Peek at the first 15 rows
hg38_ashm_regions %>%
    slice_head(n=15) 
```

| chr_name | hg38_start |  hg38_end | gene    | region   | regulatory_comment |
|:---------|-----------:|----------:|:--------|:---------|:-------------------|
| chr1     |    6601422 |   6602642 | KLHL21  | TSS      | NA                 |
| chr1     |   23559093 |  23559344 | ID3     | TSS      | NA                 |
| chr1     |   28506039 |  28509827 | RCC1    | TSS      | NA                 |
| chr1     |   30756165 |  30759164 | LAPTM5  | TSS      | NA                 |
| chr1     |  150578338 | 150579659 | MCL1    | intron   | NA                 |
| chr2     |   88605321 |  88609578 | EIF2Ak3 | intron   | NA                 |
| chr2     |   88625938 |  88628063 | EIF2Ak3 | TSS      | NA                 |
| chr2     |  231707930 | 231709587 | PTMA    | TSS      | NA                 |
| chr1     |  157699700 | 157701509 | FCRL3   | TSS      | NA                 |
| chr1     |  203305570 | 203306650 | BTG2    | intron   | active_promoter    |
| chr1     |  226677156 | 226685751 | ITPKB   | intron   | weak_enhancer      |
| chr1     |  226732862 | 226740184 | ITPKB   | TSS      | active_promoter    |
| chr1     |  226733387 | 226740281 | ITPKB   | intron-1 | enhancer           |
| chr1     |  244860200 | 244865781 | HNRNPU  | TSS      | NA                 |
| chr2     |   60546654 |  60556351 | BCL11A  | TSS      | active_promoter    |

## Other resources

- `colour_codes` A data frame with colour codes (HEX) arranged into
  different categories, groups.
- `dhitsig_genes_with_weights` A data frame with double hit signature
  genes (both as ensembl IDs and Hugo symbols) and importance scores.
- `gambl_metadata` A data frame with metadata for a collection of GAMBL
  samples. This represents a collection of whole genome, exome,
  targeted, RNA, and PrometION sequencing samples available as a data
  set known as GAMBL. This object rather serves an FYI purpose as not
  all samples listed here are published and bundled with GAMBLR.data.
- `hgnc2pfam.df` A dataset containing the mapping table between Hugo
  symbol, UniProt ID, and Pfam ACC. This dataset comes from the g3viz
  package and was obtained via this URL:
  https://github.com/morinlab/g3viz/tree/master/data
- `hotspots_annotations` Hotspot coordinates used in the feature
  annotation during matrix assembly of data for cFL classifier.
- `mirage_metrics` A data frame providing the data reported in the
  Supplemental Table of the MIRAGE study by [Dreval et al,
  2022](https://doi.org/10.1182/blood.2022016095)
- `mutation.table.df` A data frame providing the linkage between Variant
  Classification, Mutation_Class, and Short_Name for the simple somatic
  mutations.
- reddy_genes A data frame of the genes reported as significantly
  mutated by the study of [Reddy et al,
  2017](https://doi.org/10.1016/j.cell.2017.09.027)
- `wright_genes_with_weights` Wright genes with weight values from the
  study by [Scott et al, 2014](https://10.1182/blood-2013-11-536433).
