# Quantify the APA dynamics between canonical and non-canonical PACs in a gene.

Quantify the APA dynamics between canonical and non-canonical PACs in a
gene.

## Usage

``` r
Quantify.CNCAPA(QpolyA,colData,contrast)
```

## Arguments

- QpolyA:

  A QuantifyPolyA object generated after poly(A) site clustering and
  annotation.

- colData:

  A data.frame object containing the information of experimental design.

- contrast:

  A two element vector specifying the two samples to be compared, e.g.
  contrast=c("condition","Brain","UHR").

## Value

A tibble object containing APA dynamics metrics of all APA genes.
