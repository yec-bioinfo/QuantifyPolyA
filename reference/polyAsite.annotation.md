# Annotate the poly(A) sites based on the genomic structures

Annotate the poly(A) sites based on the genomic structures

## Usage

``` r
polyAsite.annotation(polyA, gff, seq.levels=NA)
```

## Arguments

- polyA:

  A data.frame containing the info of clustered poly(A) sites - PAC.

- gff:

  A genome annotation file in GFF or GTF format, GTF format is
  recommended.

- seq.levels:

  A character vector of the names of Chromosomes/Contigs. This parameter
  is useful when the Chromosomes/Contigs names in genome file (used for
  short reads mapping) and genome annotation file are not consistent.

## Value

The polyA data.frame will be return with three additional columns
'gene_id', 'distance, 'type'.
