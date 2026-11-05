# Quantify APA dynamics among canonical PACs

Quantify APA dynamics among canonical PACs of a specific gene, canonical
PACs represent PACs located at 3' UTR.

## Usage

``` r
Quantify.CanonicalAPA(QpolyA,colData,contrast)
```

## Arguments

- QpolyA:

  A QuantifyPolyA object generated after poly(A) site clustering and
  annotation.

- colData:

  A data.frame object containing the information of experimental design.

- contrast:

  A three element vector specifying the two samples to be compared, e.g.
  contrast=c("condition","Brain","UHR").

## Value

A tibble object containing APA dynamics metrics of all APA genes.
