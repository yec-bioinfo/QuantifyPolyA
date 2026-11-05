# Quantify the APA dynamics across a whole gene

Quantify the APA dynamics across a whole gene

## Usage

``` r
Quantify.GeneAPA(QpolyA,colData,contrast)
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
