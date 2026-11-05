# Cluster the one-dimensional poly(A) sites into groups based on local density

Cluster the one-dimensional poly(A) sites into groups based on local
density

## Usage

``` r
findPeaks(sub_pos,sub_wts,min_delta=24)
```

## Arguments

- sub_pos:

  The genomic coordinates of poly(A) sites.

- sub_wts:

  The numbers of reads support each poly(A) site in 'sub_pos'.

- min_delta:

  A cutoff to limit the minimum distance between two centers.

## Value

A tibble object with five columns 'group', 'start', 'end', 'sum.wts',
'center'.
