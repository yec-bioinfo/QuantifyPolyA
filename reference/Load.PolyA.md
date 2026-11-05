# Load raw poly(A) site data

Load the raw poly(A) site information extracted from 3' end sequencing
data.

## Usage

``` r
Load.PolyA(files,dir)
```

## Arguments

- files:

  A character vector specifying the names of files containing poly(A)
  site information, each file contains 4 columns 'seqnames', 'strand',
  'coord', 'score'.

- dir:

  A string specifying the directory of files containing poly(A) site
  information.

## Value

A QuantifyPolyA object containing all raw poly(A) site information.
