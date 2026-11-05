# Visualize poly(A) profiles

Visualize poly(A) profiles of a specific gene across multiple samples

## Usage

``` r
Visualize.PolyA(QpolyA,anno,gene.id,sample.names,coord.lim)
```

## Arguments

- QpolyA:

  A QuantifyPolyA object containing clean poly(A) sites and PAC table.

- anno:

  A Granges object of gene model info. It could be generated using the
  function 'import' of R package 'rtracklayer', e.g. gtf =
  import('./Mus_musculus.GRCm38.84.gtf')

- gene.id:

  A string indicating the gene to be visualized, e.g. gene.id =
  'ENSMUSG00000073490'.

- sample.names:

  A character vector specifying the samples to be visualized.

- coord.lim:

  Two numeric values, specifying the left and right limit of xaxis.
