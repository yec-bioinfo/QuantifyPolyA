# Class QuantifyPolyA.

Class `QuantifyPolyA.` defines a poly(A) dataset.

## Slots

- `sample_names`:

  A character vector storing names of samples.

- `pre.polyA`:

  A list storing un-clustered poly(A) sites.

- `polyA`:

  A data.frame storing PACs (Poly(A) Site Clusters) generated from the
  weighted density peak clustering.

- `simple.clusters`:

  A data.frame storing PACs generated from the k-nt iteratively
  clustering.

- `split.clusters`:

  A data.frame storing only the split PACs.
