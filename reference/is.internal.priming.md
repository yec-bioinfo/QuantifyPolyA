# Check whether the poly(A) sites are resulted from internal priming artifacts or not.

Check whether the poly(A) sites are resulted from internal priming
artifacts or not.

## Usage

``` r
is.internal.priming(polyA, fasta, flank_len = 15, win_size = 10, min_A = 8)
```

## Arguments

- polyA:

  A data.frame containing the information of raw poly(A) sites.

- fasta:

  A string specifies the location and name of genome file in FASTA
  format.

- flank_len:

  An integer specifying the neighborhood to search for template polyA.

- win_size:

  An integer specifying the size of the sliding window.

- min_A:

  An integer specifying the minimum number of base A in a window.

## Value

A logic vector indicates whether the poly(A) sites are resulted from
internal priming artifacts or not.
