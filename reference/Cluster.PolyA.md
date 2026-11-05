# Generate PACs

Cluster poly(A) sites into clusters using a weighted density peak
clustering algorithm.

## Usage

``` r
Cluster.PolyA(QpolyA, max.gapwidth = 24, mc.cores = 4)
```

## Arguments

- QpolyA:

  A QuantifyPolyA object containing clean poly(A) sites.

- max.gapwidth:

  A cutoff limiting the max distance between two adjacent poly(A) sites
  in a poly(A) site cluster (PAC), default value is 24.

- mc.cores:

  An integer indicating the number of processors/cores to perform the
  clustering, default value is 4.
