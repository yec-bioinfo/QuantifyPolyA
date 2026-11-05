# Cluster poly(A) sites into groups by a certain distance iterative.

Cluster poly(A) sites into groups by a certain distance iterative.

## Usage

``` r
simpleCluster(points.gr, max.gapwidth = 24)
```

## Arguments

- points.gr:

  A GRanges object constructed from poly(A) sites.

- max.gapwidth:

  A cutoff to limit the maximum distance between two adjacent sites in a
  PAC.

## Value

A reduced GRanges object.
