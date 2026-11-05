# Construct genomic ranges

Construct genomic ranges

## Usage

``` r
buildGenomicRanges(seqname,position,score,strand)
```

## Arguments

- seqname:

  A vector containing the sequence (Chromosomes/Contigs) names of
  poly(A) sites.

- position:

  A numeric vector containing the genomic positions of poly(A) sites.

- score:

  A numeric vector containing the numbers of reads supporting each
  poly(A) site.

- strand:

  A character vector containing the strand information of poly(A) sites.

## Value

A Granges object.
