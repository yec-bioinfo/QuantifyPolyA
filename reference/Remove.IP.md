# Remove IP

Remove poly(A) sites which potentially result from internal priming
events.

## Usage

``` r
Remove.IP(QpolyA, fasta, flank_len = 15, win_size = 10, min_A = 8)
```

## Arguments

- QpolyA:

  A QuantifyPolyA object containing all raw poly(A) sites.

- fasta:

  A string specifying the location and name of the corresponding genome
  file in FASTA format.

- flank_len:

  An integer specifying the neighborhood to search for template polyA.

- win_size:

  An integer specifying the size of the sliding window.

- min_A:

  An integer specifying the minimum number of base A in a window.

## Value

A QuantifyPolyA object containing poly(A) sites with IP removed.
