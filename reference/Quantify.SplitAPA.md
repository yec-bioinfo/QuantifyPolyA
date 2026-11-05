# Quantify APA dynamics among split PACs

Quantify APA dynamics among split PACs identified by the weighted
density peak clustering.

## Usage

``` r
Quantify.SplitAPA(QpolyA,colData,contrast)
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

A tibble object containing APA dynamics metrics of all split PACs.
