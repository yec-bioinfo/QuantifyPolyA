# Filter out PACs

Filter out PACs with low count of supporting reads across multiple
samples.

## Usage

``` r
Filter.PolyA(QpolyA,min_count = 10, min_sample = 1)
```

## Arguments

- QpolyA:

  A QuantifyPolyA object containing clean poly(A) sites and PAC table.

- min_count:

  A fixed value specifying the minimum read count in a PAC.

- min_sample:

  A fixed value specifying the minimum number of samples with min_count
  reads in a PAC.
