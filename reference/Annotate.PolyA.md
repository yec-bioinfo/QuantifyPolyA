# Annotate the PACs

Annotate the PACs based on existed genome annotations, and generate the
PAC table of all experimental samples.

## Usage

``` r
Annotate.PolyA(QpolyA,gff,seq.levels=NA)
```

## Arguments

- QpolyA:

  A QuantifyPolyA object containing clean poly(A) sites and PACs.

- gff:

  A genome annotation file in GFF or GTF format, GTF format is
  recommended.

- seq.levels:

  A character vector of the names of Chromosomes/Contigs. This parameter
  is useful when the Chromosomes/Contigs names in genome file (used for
  short reads mapping) and genome annotation file are not consistent.
