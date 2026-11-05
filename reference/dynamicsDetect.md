# Detect APA dynamics between two samples

Detect APA dynamics between two samples

## Usage

``` r
dynamicsDetect(sample_x, sample_y, sample_strand, sample_center)
```

## Arguments

- sample_x:

  A data.frame object with multiple rows and columns, each row
  represents a PAC, each column represents a replicate.

- sample_y:

  A data.frame object with multiple rows and columns, each row
  represents a PAC, each column represents a replicate.

- sample_strand:

  A vector containing the strand information of PACs.

- sample_center:

  A vector containing the coordinates (centers) information of PACs.

## Value

A data.frame object containing three columns 'pd', 'r', 'p.value'.
