# Package index

## Workflow

Main poly(A) quantification workflow

- [`Load.PolyA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Load.PolyA.md)
  : Load raw poly(A) site data
- [`Remove.IP()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Remove.IP.md)
  : Remove IP
- [`Cluster.PolyA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Cluster.PolyA.md)
  : Generate PACs
- [`Annotate.PolyA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Annotate.PolyA.md)
  : Annotate the PACs
- [`Filter.PolyA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Filter.PolyA.md)
  : Filter out PACs

## APA dynamics

- [`Quantify.SplitAPA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Quantify.SplitAPA.md)
  : Quantify APA dynamics among split PACs
- [`Quantify.CanonicalAPA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Quantify.CanonicalAPA.md)
  : Quantify APA dynamics among canonical PACs
- [`Quantify.CNCAPA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Quantify.CNCAPA.md)
  : Quantify the APA dynamics between canonical and non-canonical PACs
  in a gene.
- [`Quantify.GeneAPA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Quantify.GeneAPA.md)
  : Quantify the APA dynamics across a whole gene

## Visualization and statistics

- [`Visualize.PolyA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Visualize.PolyA.md)
  : Visualize poly(A) profiles
- [`DESeq2.PolyA()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/DESeq2.PolyA.md)
  : Perform basic differential expression analysis using DESeq2 package

## Classes and low-level functions

- [`QuantifyPolyA-class`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/QuantifyPolyA-class.md)
  : Class QuantifyPolyA.
- [`show(`*`<QuantifyPolyA>`*`)`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/show-methods.md)
  : Method show.
- [`is.internal.priming()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/is.internal.priming.md)
  : Check whether the poly(A) sites are resulted from internal priming
  artifacts or not.
- [`polyAsite.annotation()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/polyAsite.annotation.md)
  : Annotate the poly(A) sites based on the genomic structures
- [`buildGenomicRanges()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/buildGenomicRanges.md)
  : Construct genomic ranges
- [`simpleCluster()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/simpleCluster.md)
  : Cluster poly(A) sites into groups by a certain distance iterative.
- [`findPeaks()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/findPeaks.md)
  : Cluster the one-dimensional poly(A) sites into groups based on local
  density
- [`dynamicsDetect()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/dynamicsDetect.md)
  : Detect APA dynamics between two samples
- [`FU()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/FU.md)
  : Calculate the Pearson correlation coefficient of two PAC vectors.
- [`Motif.Search()`](https://yec-bioinfo.github.io/QuantifyPolyA/reference/Motif.Search.md)
  : Motif search
