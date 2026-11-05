#' Class QuantifyPolyA.
#'
#' Class \code{QuantifyPolyA.} defines a poly(A) dataset.
#'
#' @name QuantifyPolyA-class
#' @rdname QuantifyPolyA-class
#' @slot sample_names A character vector storing names of samples.
#' @slot pre.polyA A list storing un-clustered poly(A) sites.
#' @slot polyA A data.frame storing PACs (Poly(A) Site Clusters) generated from the weighted density peak clustering.
#' @slot simple.clusters A data.frame storing PACs generated from the k-nt iteratively clustering.
#' @slot split.clusters A data.frame storing only the split PACs.
#' @importFrom methods setClass setMethod
#' @exportClass QuantifyPolyA
#'
setClass("QuantifyPolyA",slots=list(sample_names="character",
                                    pre.polyA="list",
                                    polyA='data.frame',
                                    simple.clusters='data.frame',
                                    split.clusters ='data.frame'))
#' Method show.
#'
#' @rdname show-methods
#' @aliases show,ANY-method
#' @param object A QuantifyPolyA object.
#' @importFrom utils head
setMethod("show",
          "QuantifyPolyA",
          function(object) {
            cat('Number of samples:',length(object@sample_names),'\n')
            cat('Number of PAS in each sample:\n')
            print(sapply(object@pre.polyA, nrow))
            
            if (nrow(object@polyA)>0){
              print(head(object@polyA[, 1:min(11, ncol(object@polyA)), drop=F],2))
              cat(sprintf('%s...[%d x %d]%s','@polyA',nrow(object@polyA), ncol(object@polyA),'\n'))
            }
            
            if (nrow(object@simple.clusters)>0){
              print(head(object@simple.clusters[, 1:8, drop=F],2))
              cat(sprintf('%s...[%d x %d]%s','@simple.clusters',nrow(object@simple.clusters), ncol(object@simple.clusters),'\n'))
            }
            
            if (nrow(object@split.clusters)>0){
              print(head(object@split.clusters[ , , drop=F],2))
              cat(sprintf('%s...[%d x %d]%s','@split.clusters',nrow(object@split.clusters), ncol(object@split.clusters),'\n'))
            }
          }
)


##############################################
#           Load raw poly(A) data            #
##############################################
#' @title Load raw poly(A) site data
#' @description Load the raw poly(A) site information extracted from 3' end sequencing data.
#' @name Load.PolyA
#' @usage Load.PolyA(files,dir)
#' @param files A character vector specifying the names of files containing poly(A) site information, each file contains 4 columns 'seqnames', 'strand', 'coord', 'score'.
#' @param dir A string specifying the directory of files containing poly(A) site information.
#' @return A QuantifyPolyA object containing all raw poly(A) site information.
#' @importFrom tools file_path_sans_ext
#' @importFrom utils read.table
#' @importFrom methods new
#' @export
#'
Load.PolyA <- function(files,dir){
  # Check parameters.
  if (missing(files)){
    if (missing(dir)){
      stop("Parameter 'files' or 'dir' should be provided!")
    } else{
      if (!dir.exists(dir)) stop(paste('Directory',dir,'does not exist!'))
      files = list.files(dir,'\\.bed$',full.names = T)
      if (length(files) == 0) stop(paste('No bed files was found in directory',dir,'!'))
    }
  } else{
    for (file in files) {
      if (!file.exists(file)) stop(paste('Input file',file,'does not exist!'))
    }
  }
  
  # Load poly(A) file.
  alt_names = file_path_sans_ext(basename(files))
  pre.polyA = list()
  
  for (i in 1:length(files)) {
    pre.polyA[[alt_names[i]]] = read.table(file = files[i],sep = '\t',stringsAsFactors = FALSE)
    colnames(pre.polyA[[alt_names[i]]]) = c('seqnames','strand','coord','score')
  }
  print('Load polyA files finished!')
  QpolyA = new("QuantifyPolyA",sample_names = alt_names, pre.polyA = pre.polyA)
}


##############################################
#       Remove internal priming events       #
##############################################
#' @title Remove IP
#' @description Remove poly(A) sites which potentially result from internal priming events.
#' @name Remove.IP
#' @usage Remove.IP(QpolyA, fasta, flank_len = 15, win_size = 10, min_A = 8)
#' @param QpolyA A QuantifyPolyA object containing all raw poly(A) sites.
#' @param fasta A string specifying the location and name of the corresponding genome file in FASTA format.
#' @param flank_len An integer specifying the neighborhood to search for template polyA.
#' @param win_size An integer specifying the size of the sliding window.
#' @param min_A An integer specifying the minimum number of base A in a window.
#' @return A QuantifyPolyA object containing poly(A) sites with IP removed.
#' @export
#'
Remove.IP <- function(QpolyA,fasta,flank_len=15,win_size=10,min_A=8){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!file.exists(fasta)) stop(paste('Fasta file',fasta,'does not exist!'))
  if (!is.numeric(flank_len)) stop("'flank_len' is not a number!")
  if (flank_len<=0) stop("'flank_len' should be larger than 0, the default value is 15!")
  if (!is.numeric(win_size)) stop("'win_size' is not a number!")
  if (win_size<=0) stop("'win_size' should be larger than 0, the default value is 10!")
  if (!is.numeric(min_A)) stop("'min_A' is not a number!")
  if (min_A<=0) stop("'min_A' should be larger than 0, the default value is 8!")
  
  # Remove internal priming events.
  for (alt_name in QpolyA@sample_names) {
    is.IP = is.internal.priming(QpolyA@pre.polyA[[alt_name]],fasta = fasta,flank_len,win_size,min_A)
    print(paste(sum(is.IP),'internal priming events were found in sample', alt_name,'!'))
    QpolyA@pre.polyA[[alt_name]] = QpolyA@pre.polyA[[alt_name]][!is.IP,]
  }
  return(QpolyA)
}


##############################################
#          Cluster poly(A) sites             #
##############################################
#' @title Generate PACs
#' @description Cluster poly(A) sites into clusters using a weighted density peak clustering algorithm.
#' @name Cluster.PolyA
#' @usage Cluster.PolyA(QpolyA, max.gapwidth = 24, mc.cores = 4)
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites.
#' @param max.gapwidth A cutoff limiting the max distance between two adjacent poly(A) sites in a poly(A) site cluster (PAC), default value is 24.
#' @param mc.cores An integer indicating the number of processors/cores to perform the clustering, default value is 4.
#' @importFrom pbmcapply pbmcmapply
#' @importFrom dplyr bind_rows
#' @export
#'
Cluster.PolyA <- function(QpolyA, max.gapwidth = 24, mc.cores = 4){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.numeric(max.gapwidth)) stop("'max.gapwidth' is not a number!")
  if (max.gapwidth<=0) stop("'max.gapwidth' should be larger than 0, the default value is 24!")
  if (!is.numeric(mc.cores)) stop("'mc.cores' is not a number!")
  if (mc.cores<=0) stop("'mc.cores' should be larger than 0, the default value is 4!")
  
  # Declare
  seqnames = strand = coord = score = sum.wts = NULL
  
  # Combine data of different samples.
  polyA = bind_rows(QpolyA@pre.polyA) %>%
    group_by(seqnames,strand,coord) %>%
    summarise(score = sum(score),.groups='drop')
  
  # Simple clustering by distance.
  print('Simple clustering by distance!')
  points.gr = buildGenomicRanges(seqname = polyA$seqnames, position = polyA$coord,
                                 score = polyA$score, strand = polyA$strand)
  
  simple.clusters = simpleCluster(points.gr, max.gapwidth = max.gapwidth)
  simple.clusters.df = as.data.frame(simple.clusters)
  simple.clusters.df$split_label = NA
  
  # Re-clustering of polyA sites with large width using weighted density peak calling algorithm.
  print('Re-clustering by weighted density peak clustering algorithm!')
  idx = which(simple.clusters.df$width > max.gapwidth)
  pos = extractList(points.gr@ranges@start,simple.clusters$revmap[idx])
  wts = extractList(points.gr$score,simple.clusters$revmap[idx])
  
  split.clusters = pbmcmapply(findPeaks,pos,wts,SIMPLIFY = F, mc.cores = mc.cores)
  
  lens = sapply(split.clusters,nrow)
  split.clusters.df = bind_rows(split.clusters[lens!=1], .id = "split_label")
  
  split.clusters.df$seqnames = rep(simple.clusters.df$seqnames[idx[lens!=1]],times = lens[lens!=1])
  split.clusters.df$strand = rep(simple.clusters.df$strand[idx[lens!=1]],times = lens[lens!=1])
  split.clusters.df$width = split.clusters.df$end - split.clusters.df$start + 1
  split.clusters.df = dplyr::rename(split.clusters.df,score=sum.wts)
  
  # Generate final clusters.
  polyA = simple.clusters.df[-idx[lens!=1],c("seqnames","start","end","width","strand","score","center","split_label")]
  polyA = rbind(polyA,split.clusters.df[,c("seqnames","start","end","width","strand","score","center","split_label")])
  
  QpolyA@simple.clusters = simple.clusters.df
  QpolyA@split.clusters = split.clusters.df
  QpolyA@polyA = polyA
  return(QpolyA)
}


##############################################
#    Annotate and quantify polyA sites       #
##############################################
#' @title Annotate the PACs
#' @description Annotate the PACs based on existed genome annotations, and generate the PAC table of all experimental samples.
#' @name Annotate.PolyA
#' @usage Annotate.PolyA(QpolyA,gff,seq.levels=NA)
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites and PACs.
#' @param gff A genome annotation file in GFF or GTF format, GTF format is recommended.
#' @param seq.levels A character vector of the names of Chromosomes/Contigs. This parameter is useful when the Chromosomes/Contigs names in genome file (used for short reads mapping) and genome annotation file are not consistent.
#' @importFrom S4Vectors countQueryHits aggregate
#' @export
#'
Annotate.PolyA <- function(QpolyA,gff,seq.levels=NA){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!file.exists(gff)) stop(paste('Annotation file',gff,'does not exist!'))
  
  # Annotate polyA sites
  print('Annotate polyA sites!')
  polyA = polyAsite.annotation(QpolyA@polyA,gff,seq.levels)
  
  # Generate full table of polyA sites
  print('Generate full table of polyA sites!')
  polyA.gr = GRanges(seqnames = polyA$seqnames,ranges = IRanges(start = polyA$start, end = polyA$end),strand = polyA$strand)
  
  # Declare
  score = NULL
  
  for (name in QpolyA@sample_names) {
    points.gr = GRanges(seqnames = QpolyA@pre.polyA[[name]]$seqnames,
                        ranges = IRanges(start = QpolyA@pre.polyA[[name]]$coord, end = QpolyA@pre.polyA[[name]]$coord),
                        strand = QpolyA@pre.polyA[[name]]$strand,
                        score = QpolyA@pre.polyA[[name]]$score)
    
    hits = findOverlaps(polyA.gr,points.gr,ignore.strand=FALSE)
    agg = S4Vectors::aggregate(points.gr, hits, score=sum(score))
    polyA[,name] = 0
    polyA[countQueryHits(hits) > 0L,name] = agg$score
  }
  QpolyA@polyA = polyA
  return(QpolyA)
}


##############################################
#           Visualize polyA sites            #
##############################################
#' @title Visualize poly(A) profiles
#' @description Visualize poly(A) profiles of a specific gene across multiple samples
#' @name Visualize.PolyA
#' @usage Visualize.PolyA(QpolyA,anno,gene.id,sample.names,coord.lim)
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites and PAC table.
#' @param anno A Granges object of gene model info. It could be generated using the function 'import' of R package 'rtracklayer', e.g. gtf = import('./Mus_musculus.GRCm38.84.gtf')
#' @param gene.id A string indicating the gene to be visualized, e.g. gene.id = 'ENSMUSG00000073490'.
#' @param sample.names A character vector specifying the samples to be visualized.
#' @param coord.lim Two numeric values, specifying the left and right limit of xaxis.
#' @importFrom ggplot2 ggplot autoplot aes geom_point labs theme xlim facet_grid
#' @importFrom ggalt geom_lollipop
#' @importFrom ggbio tracks theme_tracks_sunset
#' @importFrom tidyr gather
#' @importFrom S4Vectors values values<-
#' @export
#'
Visualize.PolyA <- function(QpolyA,anno,gene.id,sample.names,coord.lim){
  options(warn = -1)
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (class(anno)!='GRanges') stop('anno is not a GRanges object!')
  if (missing(gene.id)) stop('gene.id is not set!')
  if (!(gene.id %in% anno@elementMetadata@listData$gene_id)) stop(paste(c(gene.id,' dose not exist, please check the gene id!'),collapse = ''))
  if (missing(sample.names)) {
    sample.names = QpolyA@sample_names
  }else{
    sample.names = QpolyA@sample_names[str_detect(QpolyA@sample_names,paste0(sample.names,collapse = '|'))]
  }
  
  # Declare
  gene_id = type = coord = cluster = center = seqnames = strand = score = NULL
  
  # Extract sites
  sub.polyA = subset(QpolyA@polyA,gene_id %in% gene.id & type != 'intergenic' & score>=10)
  sub.polyA = sub.polyA[order(sub.polyA$center),]
  sub.polyA$cluster = 1:nrow(sub.polyA)
  sub.sites = c()
  for (name in sample.names) {
    tmp = subset(QpolyA@pre.polyA[[name]],seqnames==sub.polyA$seqnames[1] & strand==sub.polyA$strand[1]&
                   coord <= max(sub.polyA$end) & coord >= min(sub.polyA$start))
    tmp$sample = name
    tmp$cluster = 0
    for (i in 1:nrow(sub.polyA)) {
      idx = which(tmp$coord >= sub.polyA$start[i] & tmp$coord <= sub.polyA$end[i])
      tmp$cluster[idx] = i
    }
    tmp = subset(tmp,cluster!=0)
    sub.sites = rbind(sub.sites,tmp)
  }
  
  # Extract annotation info
  gene = anno[which(anno@elementMetadata@listData$gene_id==gene.id&
                      (anno@elementMetadata@listData$type %in% c('gene')))];
  gene.name = gene$gene_name
  exon = anno[which(anno@elementMetadata@listData$gene_id==gene.id&
                      (anno@elementMetadata@listData$type %in% c('exon')))];
  colnames(values(exon))[4] = "model"
  values(exon)[4] = values(exon)[2]
  exon@elementMetadata@listData$tx_name = exon@elementMetadata@listData$transcript_name
  exon@elementMetadata@listData$tx_id = exon@elementMetadata@listData$transcript_id
  grl = S4Vectors::split(exon, exon$tx_id)
  
  # Generate figure
  if (missing(coord.lim)) {
    coord.lim = gene
  }
  p.isoform = ggbio::autoplot(grl) +
    ggbio::xlim(coord.lim)
  p.points = ggplot() +
    geom_lollipop(data= gather(sub.polyA,sample.names,key = "sample", value = "score"),
                  aes(x=center,y = score,fill=factor(cluster)),point.size=3,shape = 23,alpha=0.5)+
    geom_point(data=sub.sites,aes(x=coord,y=score,color=factor(cluster)))+
    facet_grid(sample~.,scales = "free") +
    labs(color='Cluster') +
    ggbio::xlim(coord.lim)
  #scale_y_log10()
  
  # Plot
  tracks("isoforms"= p.isoform,
         "poly(A) site" =p.points,
         heights = c(0.5,2),
         title=paste('Gene:  ',gene.id,'-',gene.name))+
    theme_tracks_sunset()+
    theme(legend.position = 'none',legend.justification=c(1,1))
}


##############################################
#       Filter low count polyA sites         #
##############################################
#' @title Filter out PACs
#' @description Filter out PACs with low count of supporting reads across multiple samples.
#' @name Filter.PolyA
#' @usage Filter.PolyA(QpolyA,min_count = 10, min_sample = 1)
#' @param QpolyA A QuantifyPolyA object containing clean poly(A) sites and PAC table.
#' @param min_count A fixed value specifying the minimum read count in a PAC.
#' @param min_sample A fixed value specifying the minimum number of samples with min_count reads in a PAC.
#' @export
#'
Filter.PolyA <- function(QpolyA,min_count = 10,min_sample = 1){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.numeric(min_count)) stop("'min_count' is not a number!")
  if (min_count<=0) stop("'min_count' should be larger than 0, the default value is 10!")
  if (!is.numeric(min_sample)) stop("'min_sample' is not a number!")
  if (min_sample<=0) stop("'min_sample' should be larger than 0, the default value is 1!")
  
  QpolyA@polyA = subset(QpolyA@polyA, rowSums(QpolyA@polyA[,QpolyA@sample_names]>=min_count) >= min_sample)
  rownames(QpolyA@polyA) = paste('PA',1:nrow(QpolyA@polyA),sep = '')
  return(QpolyA)
}


##############################################
#          Perform DESeq2 analysis           #
##############################################
#' Perform basic differential expression analysis using DESeq2 package
#' @name DESeq2.PolyA
#' @usage DESeq2.PolyA(QpolyA,colData)
#' @param QpolyA A QuantifyPolyA object.
#' @param colData A data.frame object containing the information of experimental design.
#' @importFrom DESeq2 DESeq counts vst DESeqDataSetFromMatrix
#' @importFrom grDevices colorRampPalette
#' @importFrom factoextra fviz_pca_ind
#' @importFrom ggplot2 coord_fixed xlab ylab
#' @importFrom uwot umap
#' @importFrom FactoMineR PCA
#' @importFrom SummarizedExperiment assay
#' @importFrom stats dist
#' @importFrom S4Vectors isEmpty
#' @export
#'
DESeq2.PolyA <- function(QpolyA,colData){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")
  
  # Declare
  condition = X1 = X2 = NULL
  
  # Build DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = QpolyA@polyA[,rownames(colData)],
                                colData = colData,
                                design = ~ condition)
  # Run DESeq
  dds <- DESeq(dds)
  
  # Normalized table
  polyA.normailized = cbind(QpolyA@polyA[,c(1:11)],counts(dds,normalized=TRUE))
  
  ### Data transformations and visualization
  ## Count data transformations
  vsd <- vst(dds, blind=FALSE)
  
  ## Data quality assessment by sample clustering and visualization
  # Heatmap of the sample-to-sample distances
  # library(pheatmap)
  # sampleDists <- dist(t(assay(vsd)))
  # sampleDistMatrix <- as.matrix(sampleDists)
  # rownames(sampleDistMatrix) <- vsd$condition
  # colnames(sampleDistMatrix) <- NULL
  # colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
  # sample.heatmap = pheatmap(sampleDistMatrix,
  #                           clustering_distance_rows=sampleDists,
  #                           clustering_distance_cols=sampleDists,
  #                           col=colors)
  
  # Principal component plot of the samples
  # library(FactoMineR)
  # library(factoextra)
  pcaData <- PCA(t(assay(vsd)), graph = FALSE)
  sample.pca = fviz_pca_ind(pcaData,
                            geom.ind = "point", # show points only (nbut not "text")
                            col.ind = factor(colData$condition), # color by groups
                            addEllipses = FALSE, # Concentration ellipses
                            legend.title = "Groups"
  )
  
  # tSNE plot of the samples
  # library(Rtsne) # Load package
  # set.seed(100) # Sets seed for reproducibility
  # tsne_out <- Rtsne(as.matrix(t(assay(vsd))),perplexity = 3) # Run TSNE
  #
  # sample.tsne = ggplot(data.frame(tsne_out$Y,condition=colData$condition), aes(X1, X2, color=condition)) +
  #   geom_point(size=3) +
  #   xlab(paste0("tSNE1")) +
  #   ylab(paste0("tSNE2")) +
  #   coord_fixed()
  
  # umap plot of the samples
  # library(uwot)
  # run UMAP algorithm
  umap_out <- umap(as.matrix(t(assay(vsd))), n_neighbors = 2, init = "spca")
  
  sample.umap = ggplot(data.frame(umap_out,condition=colData$condition), aes(X1, X2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("Dim1")) +
    ylab(paste0("Dim2")) +
    coord_fixed()
  
  return(list(DESeq2.Result = dds,
              PCA.Plot = sample.pca,
              UMAP.Plot = sample.umap))
}


#########################################################
#    Quantify APA dynamics among split polyA sites      #
#########################################################
#' @title Quantify APA dynamics among split PACs
#' @description Quantify APA dynamics among split PACs identified by the weighted density peak clustering.
#' @name Quantify.SplitAPA
#' @usage Quantify.SplitAPA(QpolyA,colData,contrast)
#' @param QpolyA A QuantifyPolyA object generated after poly(A) site clustering and annotation.
#' @param colData A data.frame object containing the information of experimental design.
#' @param contrast A three element vector specifying the two samples to be compared, e.g. contrast=c("condition","Brain","UHR").
#' @return A tibble object containing APA dynamics metrics of all split PACs.
#' @importFrom dplyr filter n do
#' @importFrom S4Vectors isEmpty
#' @export
#'
Quantify.SplitAPA <- function(QpolyA,colData,contrast){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")
  
  control = rownames(colData)[colData[,contrast[1]]==contrast[2]]
  if (isEmpty(control)) stop(paste('Sample name',contrast[2],'is not valid!'))
  treat = rownames(colData)[colData[,contrast[1]]==contrast[3]]
  if (isEmpty(treat)) stop(paste('Sample name',contrast[3],'is not valid!'))
  
  # Declare
  split_label = . = NULL
  
  # Select out split PACs
  split.Sites = QpolyA@polyA %>% filter(!is.na(split_label)) %>%
    group_by(split_label) %>%
    filter(n()>=2)
  
  res = split.Sites %>% group_by(split_label) %>% do(dynamicsDetect(.[,control],.[,treat],.$strand,.$center))
}


############################################################
#    Quantify APA dynamics among canonical polyA sites     #
############################################################
#' @title Quantify APA dynamics among canonical PACs
#' @description Quantify APA dynamics among canonical PACs of a specific gene, canonical PACs represent PACs located at 3' UTR.
#' @name Quantify.CanonicalAPA
#' @usage Quantify.CanonicalAPA(QpolyA,colData,contrast)
#' @param QpolyA A QuantifyPolyA object generated after poly(A) site clustering and annotation.
#' @param colData A data.frame object containing the information of experimental design.
#' @param contrast A three element vector specifying the two samples to be compared, e.g. contrast=c("condition","Brain","UHR").
#' @return A tibble object containing APA dynamics metrics of all APA genes.
#' @importFrom S4Vectors isEmpty
#' @export
#'
Quantify.CanonicalAPA <- function(QpolyA,colData,contrast){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")
  
  control = rownames(colData)[colData[,contrast[1]]==contrast[2]]
  if (isEmpty(control)) stop(paste('Sample name',contrast[2],'is not valid!'))
  treat = rownames(colData)[colData[,contrast[1]]==contrast[3]]
  if (isEmpty(treat)) stop(paste('Sample name',contrast[3],'is not valid!'))
  
  # Declare
  gene_id = type = . = NULL
  
  canonical.Sites = QpolyA@polyA %>% filter(type %in% c('3UTR','ext_3UTR')) %>%
    group_by(gene_id) %>%
    filter(n()>=2)
  
  res = canonical.Sites %>% group_by(gene_id) %>% do(dynamicsDetect(.[,control],.[,treat],.$strand,.$center))
}


################################################################################
#    Quantify APA dynamics between canonical and non-canonical polyA sites     #
################################################################################
#' Quantify the APA dynamics between canonical and non-canonical PACs in a gene.
#' @name Quantify.CNCAPA
#' @usage Quantify.CNCAPA(QpolyA,colData,contrast)
#' @param QpolyA A QuantifyPolyA object generated after poly(A) site clustering and annotation.
#' @param colData A data.frame object containing the information of experimental design.
#' @param contrast A two element vector specifying the two samples to be compared, e.g. contrast=c("condition","Brain","UHR").
#' @return A tibble object containing APA dynamics metrics of all APA genes.
#' @importFrom dplyr vars full_join mutate summarise_at
#' @importFrom tidyselect all_of
#' @importFrom S4Vectors isEmpty
#' @export
#'
Quantify.CNCAPA <- function(QpolyA,colData,contrast){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")
  
  control = rownames(colData)[colData[,contrast[1]]==contrast[2]]
  if (isEmpty(control)) stop(paste('Sample name',contrast[2],'is not valid!'))
  treat = rownames(colData)[colData[,contrast[1]]==contrast[3]]
  if (isEmpty(treat)) stop(paste('Sample name',contrast[3],'is not valid!'))
  
  # Declare
  type = gene_id = canonical.label = strand = center = . = NULL
  
  # Group PACs
  gene.Sites = QpolyA@polyA %>% filter(type != 'intergenic') %>%
    mutate(canonical.label = type %in% c('3UTR','ext_3UTR')) %>%
    group_by(gene_id,canonical.label)
  
  gene.Sites.Count = gene.Sites %>% summarise_at(vars(all_of(c(control,treat))),sum) %>%
    group_by(gene_id) %>%
    filter(n()>=2)
  
  gene.Sites.Info = gene.Sites %>% summarise(strand = strand[1],center = mean(center),.groups = 'drop') %>%
    group_by(gene_id) %>%
    filter(n()>=2)
  
  gene.Sites = full_join(gene.Sites.Info,gene.Sites.Count,by=c('gene_id','canonical.label'))
  
  res = gene.Sites %>% group_by(gene_id) %>% do(dynamicsDetect(.[,control],.[,treat],.$strand,.$center))
}


############################################################
#        Quantify APA dynamics across a whole gene         #
############################################################
#' Quantify the APA dynamics across a whole gene
#' @name Quantify.GeneAPA
#' @usage Quantify.GeneAPA(QpolyA,colData,contrast)
#' @param QpolyA A QuantifyPolyA object generated after poly(A) site clustering and annotation.
#' @param colData A data.frame object containing the information of experimental design.
#' @param contrast A two element vector specifying the two samples to be compared, e.g. contrast=c("condition","Brain","UHR").
#' @return A tibble object containing APA dynamics metrics of all APA genes.
#' @importFrom S4Vectors isEmpty
#' @export
#'
Quantify.GeneAPA <- function(QpolyA,colData,contrast){
  # Check parameters.
  if (class(QpolyA) != "QuantifyPolyA") stop(paste('QpolyA should be a QuantifyPolyA object!'))
  if (!is.data.frame(colData)) stop("'colData' is not a data.frame!")
  if (isEmpty(colData$condition)) stop("'colData' does not contain a 'condition' column!")
  if (!isEmpty(setdiff(rownames(colData),QpolyA@sample_names))) stop("Inconsistency between rownames of 'colData' and sample names!")
  
  control = rownames(colData)[colData[,contrast[1]]==contrast[2]]
  if (isEmpty(control)) stop(paste('Sample name',contrast[2],'is not valid!'))
  treat = rownames(colData)[colData[,contrast[1]]==contrast[3]]
  if (isEmpty(treat)) stop(paste('Sample name',contrast[3],'is not valid!'))
  
  # Declare
  type = gene_id = . = NULL
  
  gene.Sites = QpolyA@polyA %>% filter(type != 'intergenic') %>%
    group_by(gene_id) %>%
    filter(n()>=2)
  
  res = gene.Sites %>% group_by(gene_id) %>% do(dynamicsDetect(.[,control],.[,treat],.$strand,.$center))
}


##############################################
#          Motif search function             #
##############################################
#' @title Motif search
#' @description Search canonical poly(A) signals at the upstream of poly(A) sites.
#' @name Motif.Search
#' @param Sites A data.frame containing the information of poly(A) sites.
#' @param fasta A string specifies the location and name of genome file in FASTA format.
#' @return A list containing the types and genomic locations of signals.
#' @importFrom stringr str_extract_all str_locate_all
#' @export
#'
Motif.Search <- function(Sites,fasta){
  # Check parameters.
  if (!file.exists(fasta)) stop(paste('Fasta file',fasta,'does not exist!'))
  
  # Construct a data.frame to build a bed file for extracting sequences
  seq.bed = data.frame(
    chr = as.character(Sites$seqnames),
    start = Sites$start,
    end = Sites$end,
    strand = as.character(Sites$strand),
    stringsAsFactors = FALSE
  )
  
  # Check index file of fasta file
  if (!file.exists(paste0(fasta,'.fai'))) indexFa(fasta)
  
  # Sort
  idx = order(seq.bed$chr,seq.bed$start,decreasing = FALSE)
  seq.bed.sort = seq.bed[idx,]
  
  # Extract sequences
  seq.df = get.fasta(seq.bed.sort,fasta = fasta,verbose = FALSE,check.chr = FALSE)
  
  Sites$seq[idx] = seq.df$sequence
  Sites$index[idx] = seq.df$index
  
  signal.str = c()
  signal.pos = c()
  signal.str[Sites$strand=='+'] = str_extract_all(Sites$seq[Sites$strand=='+'],'.ATAAA|A.TAAA|AA.AAA|AAT.AA|AATA.A|AATAA.', simplify = F)
  signal.pos[Sites$strand=='+'] = str_locate_all(Sites$seq[Sites$strand=='+'],'.ATAAA|A.TAAA|AA.AAA|AAT.AA|AATA.A|AATAA.')
  signal.str[Sites$strand=='-'] = str_extract_all(Sites$seq[Sites$strand=='-'],'.TTATT|T.TATT|TT.ATT|TTT.TT|TTTA.T|TTTAT.', simplify = F)
  signal.pos[Sites$strand=='-'] = str_locate_all(Sites$seq[Sites$strand=='-'],'.TTATT|T.TATT|TT.ATT|TTT.TT|TTTA.T|TTTAT.')
  
  # num = sum(sapply(plus.signal, length)>0) +  sum(sapply(minus.signal, length)>0)
  #
  # minus.signal = DNAStringSet(unlist(minus.signal))
  # minus.signal = as.character(reverseComplement(minus.signal))
  #
  # logo = ggplot() + geom_logo( c(unlist(plus.signal),minus.signal)) + theme_logo()
  # annotate('rect', xmin = 0, xmax = 7, ymin = -0.05, ymax = 1.3, alpha = .1, col='black', fill='yellow')
  
  return(list(strings = signal.str, positions = signal.pos))
}


################################################
#       Detect dynamics of poly(A) sites       #
################################################
#' Detect APA dynamics between two samples
#' @name dynamicsDetect
#' @param sample_x A data.frame object with multiple rows and columns, each row represents a PAC, each column represents a replicate.
#' @param sample_y A data.frame object with multiple rows and columns, each row represents a PAC, each column represents a replicate.
#' @param sample_strand A vector containing the strand information of PACs.
#' @param sample_center A vector containing the coordinates (centers) information of PACs.
#' @return A data.frame object containing three columns 'pd', 'r', 'p.value'.
#' @importFrom stats chisq.test
#' @importFrom dplyr pull
#' @export
#'
dynamicsDetect <- function(sample_x,sample_y,sample_strand,sample_center){
  options(warn = -1)
  
  status = 1
  if(sample_strand[1]=='-') status = -1
  
  n_x = ncol(sample_x)
  n_y = ncol(sample_y)
  name_x = colnames(sample_x)
  name_y = colnames(sample_y)
  diff_x_y = data.frame(x = rep(name_x,each = n_y),
                        y = rep(name_y,times = n_x),
                        pd = rep(NA,n_x*n_y),
                        r = rep(NA,n_x*n_y),
                        p.value = rep(NA,n_x*n_y), stringsAsFactors = FALSE)
  
  for (x in 1:n_x) {
    ind.sum_x = sum(sample_x[,x])
    if(ind.sum_x==0) next
    ind.pet_x = sample_x[,x]/ind.sum_x
    
    for (y in 1:n_y) {
      ind.sum_y = sum(sample_y[,y])
      if(ind.sum_y==0) next
      ind.pet_y = sample_y[,y]/ind.sum_y
      
      diff_x_y[(x-1)*n_y+y,'pd'] = sum(abs(ind.pet_x - ind.pet_y))/2
      diff_x_y[(x-1)*n_y+y,'r'] = FU(pull(sample_x,x),pull(sample_y,y),sample_center)*status
      
      test.data = cbind(sample_x[,x],sample_y[,y])
      test.data = test.data[which(rowSums(test.data)>0),]
      if(dim(test.data)[1]==1) {diff_x_y[(x-1)*n_y+y,'p.value'] = 1}
      else {diff_x_y[(x-1)*n_y+y,'p.value'] = chisq.test(test.data)$p.value}
    }
  }
  #diff_x_y
  data.frame(pd = mean(diff_x_y$pd), r = mean(diff_x_y$r), p.value = mean(diff_x_y$p.value))
}


################################################
#                FU function                   #
################################################
#' Calculate the Pearson correlation coefficient of two PAC vectors.
#' @name FU
#' @param s1 A numeric vector of PAC usages of sample X.
#' @param s2 A numeric vector of PAC usages of sample Y.
#' @param score A numeric vector containing the genomic coordinates of PACs.
#' @return A numeric value in a range of [-1,1].
#' @export
#'
FU <- function(s1,s2,score) {
  p1=s1/sum(s1+s2)
  p2=s2/sum(s1+s2)
  p=rbind(p1,p2)
  u=c(1,2)
  v=score
  ubar=sum(p1)*u[1]+sum(p2)*u[2]
  vbar=sum((p1+p2)*v)
  
  fz=0
  for (i in 1:2) {
    for (j in 1:length(v)) {
      fz=fz+(u[i]-ubar)*(v[j]-vbar)*p[i,j]
    }
  }
  
  fm1=0
  for (i in 1:2) {
    fm1=fm1+(u[i]-ubar)^2*sum(p[i,])
  }
  
  fm2=0
  for (j in 1:length(v)) {
    fm2=fm2+(v[j]-vbar)^2*sum(p[,j])
  }
  
  r=fz/sqrt(fm1*fm2)
  return(r)
}


##############################################
#       Check internal priming events        #
##############################################
#' Check whether the poly(A) sites are resulted from internal priming artifacts or not.
#' @name is.internal.priming
#' @param polyA A data.frame containing the information of raw poly(A) sites.
#' @param fasta A string specifies the location and name of genome file in FASTA format.
#' @param flank_len An integer specifying the neighborhood to search for template polyA.
#' @param win_size An integer specifying the size of the sliding window.
#' @param min_A An integer specifying the minimum number of base A in a window.
#' @return A logic vector indicates whether the poly(A) sites are resulted from internal priming artifacts or not.
#' @importFrom dplyr left_join
#' @importFrom bedr get.fasta
#' @importFrom Rsamtools indexFa
#' @importFrom stringr str_detect str_split
#' @importFrom matrixStats rowCumsums
#' @export
#'
is.internal.priming <- function(polyA,fasta,flank_len=15,win_size=10,min_A=8){
  # Construct a data.frame to build a bed file for extracting sequences
  seq.bed = data.frame(
    chr = polyA$seqnames,
    start = polyA$coord-flank_len,
    end = polyA$coord+flank_len,
    strand = polyA$strand,
    stringsAsFactors = FALSE
  )
  
  # Check index file of fasta file
  if (!file.exists(paste0(fasta,'.fai'))) indexFa(fasta)
  
  # Get info of fasta file
  fai = read.table(file = paste0(fasta,'.fai'))
  colnames(fai) = c('chr','len','offset','linebase','linewidth')
  
  # Check for valid regions i.e. start < 0 or end > seq.len
  seq.bed = left_join(seq.bed,fai[,c(1,2)], by='chr')
  seq.bed$start[seq.bed$start<=0] = 1
  idx = seq.bed$end>seq.bed$len
  seq.bed$end[idx] = seq.bed$len[idx]
  
  # Sort
  idx = order(seq.bed$chr,seq.bed$start,decreasing = FALSE)
  seq.bed = seq.bed[idx,]
  
  # Extract sequences
  seq.df = get.fasta(seq.bed,fasta = fasta,verbose = FALSE,check.chr = FALSE)
  
  # Check internal priming events
  is.IP = rep(TRUE,nrow(seq.bed))
  plus.idx = which(seq.bed$strand=='+')
  minus.idx = which(seq.bed$strand=='-')
  
  is.IP[plus.idx] = str_detect(seq.df$sequence[plus.idx],'A{6,}')
  is.IP[minus.idx] = str_detect(seq.df$sequence[minus.idx],'T{6,}')
  
  plus.chars = str_split(seq.df$sequence[plus.idx],'',simplify = T)=='A'
  plus.chars = rowCumsums(plus.chars,na.rm=TRUE)
  plus.chars = rowSums(cbind(plus.chars[,win_size],plus.chars[,(win_size+1):(2*flank_len)] - plus.chars[,1:(2*flank_len-win_size)])>=min_A,na.rm = TRUE)
  
  is.IP[plus.idx] = is.IP[plus.idx]|plus.chars
  
  minus.chars = str_split(seq.df$sequence[minus.idx],'',simplify = T)=='T'
  minus.chars = rowCumsums(minus.chars,na.rm=TRUE)
  minus.chars = rowSums(cbind(minus.chars[,win_size],minus.chars[,(win_size+1):(2*flank_len)] - minus.chars[,1:(2*flank_len-win_size)])>=min_A,na.rm = TRUE)
  
  is.IP[minus.idx] = is.IP[minus.idx]|minus.chars
  
  # return result
  is.IP[idx] = is.IP
  return(is.IP)
}


##############################################
#          construct genome ranges           #
##############################################
#' Construct genomic ranges
#' @name buildGenomicRanges
#' @usage buildGenomicRanges(seqname,position,score,strand)
#' @param seqname A vector containing the sequence (Chromosomes/Contigs) names of poly(A) sites.
#' @param position A numeric vector containing the genomic positions of poly(A) sites.
#' @param score A numeric vector containing the numbers of reads supporting each poly(A) site.
#' @param strand A character vector containing the strand information of poly(A) sites.
#' @return A Granges object.
#' @export
#'
buildGenomicRanges <- function(seqname,position,score,strand = '*'){
  points.gr = GRanges(seqnames = seqname,ranges = IRanges(start = position,width = 1),
                      strand = strand,score = score)
  points.gr = sort(points.gr)
}


##############################################
#       simple clustering by distance        #
##############################################
#' Cluster poly(A) sites into groups by a certain distance iterative.
#' @name simpleCluster
#' @param points.gr A GRanges object constructed from poly(A) sites.
#' @param max.gapwidth A cutoff to limit the maximum distance between two adjacent sites in a PAC.
#' @return A reduced GRanges object.
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges extractList
#' @importFrom BiocGenerics start
#' @export
#'
simpleCluster <- function(points.gr,max.gapwidth=24){
  
  # Cluster points by distance
  range.gr = reduce(points.gr,min.gapwidth=max.gapwidth,with.revmap=T,ignore.strand=FALSE)
  
  # Sum the score
  range.gr$score = sum(extractList(points.gr$score,range.gr$revmap))
  
  # Select the center
  idx = BiocGenerics::which.max(extractList(points.gr$score,range.gr$revmap))
  range.gr$center = start(points.gr)[idx+c(0,cumsum(lengths(range.gr$revmap))[1:(length(range.gr$revmap)-1)])]
  
  # Return result
  return(range.gr)
}


##############################################
#             findPeaks function             #
##############################################
#' Cluster the one-dimensional poly(A) sites into groups based on local density
#' @name findPeaks
#' @usage findPeaks(sub_pos,sub_wts,min_delta=24)
#' @param sub_pos The genomic coordinates of poly(A) sites.
#' @param sub_wts The numbers of reads support each poly(A) site in 'sub_pos'.
#' @param min_delta A cutoff to limit the minimum distance between two centers.
#' @return A tibble object with five columns 'group', 'start', 'end', 'sum.wts', 'center'.
#' @importFrom dplyr tibble %>% group_by summarise
#' @importFrom outliers scores
#' @export
#'
findPeaks <- function(sub_pos,sub_wts,min_delta=24){
  
  # Declare
  group = pos = wts = NULL
  
  ND <- length(sub_pos)
  
  if(ND<=2){
    res <- tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)])
    return(res)
  }
  
  dc <- 0.5
  #paste('Computing Rho with gaussian kernel of radius: ',dc,collapse = '')
  
  rho = rep(0,ND)
  #
  # Gaussian kernel
  #
  for(i in 1:(ND-1)){
    for(j in (i+1):ND){
      tmp_dist = abs(sub_pos[i]-sub_pos[j])
      tmp_wts = sub_wts[i]*sub_wts[j]
      rho[i]=rho[i]+exp(-(tmp_dist/dc)^2)*sub_wts[j];
      rho[j]=rho[j]+exp(-(tmp_dist/dc)^2)*sub_wts[i];
    }
  }
  for(i in 1:ND){
    rho[i]=rho[i]+sub_wts[i];
  }
  
  maxd <- abs(sub_pos[1] - sub_pos[ND])
  rho_sorted = sort(rho,decreasing = T,index.return=T)
  ordrho <- rho_sorted$ix
  rho_sorted <- rho_sorted$x
  
  delta <- rep(-1.,ND)
  nneigh <- rep(0,ND)
  
  for(ii in 2:ND){
    delta[ordrho[ii]] = maxd
    for(jj in 1:(ii-1)){
      tmp_dist = abs(sub_pos[ordrho[ii]]-sub_pos[ordrho[jj]])
      if(tmp_dist<=delta[ordrho[ii]]){
        delta[ordrho[ii]] = tmp_dist
        nneigh[ordrho[ii]] = ordrho[jj]
      }
    }
  }
  delta[ordrho[1]]=max(delta)
  
  decision.data <- data.frame(rho=rho,delta=delta,rho.delta=rho*delta)
  outlier <- rho>max(rho)/2 & delta>max(delta)/2 & delta>min_delta
  
  tmp <- decision.data$rho.delta
  tmp[which(outlier)] <- 0
  outlier <- (outlier | scores(tmp, type="chisq", prob=0.99))& delta>min_delta
  
  if(!sum(outlier)){
    res <- tibble(group=1,start=min(sub_pos),end=max(sub_pos),sum.wts=sum(sub_wts),center=sub_pos[which.max(sub_wts)])
    return(res)
  }
  
  NCLUST = 0
  cl = rep(-1,ND)
  icl = c()
  for(i in 1:ND){
    if(outlier[i]){
      NCLUST = NCLUST + 1
      cl[i] = NCLUST
      icl[NCLUST] = i
    }
  }
  #paste('NUMBER OF CLUSTERS: ',NCLUST,collapse = '')
  
  # Assignation
  for(i in 1:ND){
    if (cl[ordrho[i]]==-1){
      cl[ordrho[i]] = cl[nneigh[ordrho[i]]];
    }
  }
  
  center <- rep(1,ND)
  center[icl] <- 2
  res <- data.frame(pos= sub_pos,wts = sub_wts,group = cl,center=center)
  res <- res %>% group_by(group) %>% summarise(start=min(pos),end=max(pos),sum.wts = sum(wts),center = pos[center==2],.groups = 'drop')
}


##############################################
#            annotation function             #
##############################################
#' Annotate the poly(A) sites based on the genomic structures
#' @name polyAsite.annotation
#' @usage polyAsite.annotation(polyA, gff, seq.levels=NA)
#' @param polyA A data.frame containing the info of clustered poly(A) sites - PAC.
#' @param gff A genome annotation file in GFF or GTF format, GTF format is recommended.
#' @param seq.levels A character vector of the names of Chromosomes/Contigs. This parameter is useful when the Chromosomes/Contigs names in genome file (used for short reads mapping) and genome annotation file are not consistent.
#' @return The polyA data.frame will be return with three additional columns 'gene_id', 'distance, 'type'.
#' @importFrom GenomicFeatures cds exons genes threeUTRsByTranscript fiveUTRsByTranscript
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom GenomicRanges GRanges findOverlaps follow distance
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors countQueryHits aggregate
#' @importFrom BiocGenerics width
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @export
#'
polyAsite.annotation <- function(polyA, gff, seq.levels=NA){
  
  # build granges object for polyA
  polyA.gr = GRanges(seqnames = polyA$seqnames,ranges = IRanges(start = polyA$start, end = polyA$end),strand = polyA$strand)
  
  # import annotation info from GFF/GTF file
  txdb = txdbmaker::makeTxDbFromGFF(gff)
  if(!sum(is.na(seq.levels))){
    seqlevels(txdb) = seq.levels
  }
  
  # extract specific regions
  genes.gr = genes(txdb) # gene
  threeUTRs.gr = unlist(threeUTRsByTranscript(txdb)) # 3' UTR
  fiveUTRs.gr = unlist(fiveUTRsByTranscript(txdb)) # 5' UTR
  cds.gr = cds(txdb) # CDS
  exons.gr = exons(txdb) # EXON
  
  threeUTR.mean.len = mean(width(threeUTRs.gr),na.rm = TRUE) # mean length of 3' UTR
  
  # annotation by gene
  # Declare
  gene_id = NULL
  hits.gene = findOverlaps(polyA.gr,genes.gr,ignore.strand=FALSE)
  agg <- S4Vectors::aggregate(genes.gr, hits.gene, gene_id=BiocGenerics::paste(gene_id,collapse=';'))
  polyA$gene_id[countQueryHits(hits.gene) > 0L] <- agg$gene_id
  
  # calculate the distance between polyA site and its most upstream gene (distance of polyA site with gene annotation in previous step was set to 0)
  hits.upstream = follow(polyA.gr,genes.gr,ignore.strand=FALSE)
  idx = !is.na(hits.upstream)
  polyA$distance[idx] = distance(polyA.gr[idx],genes.gr[hits.upstream[idx]],ignore.strand=FALSE) + 1
  polyA$distance[!is.na(polyA$gene_id)] = 0
  
  # set gene_id for polyA site without gene annotation in previous step
  idx = idx & is.na(polyA$gene_id) # choose elements without gene annotation and having a precedent gene
  polyA$gene_id[idx] = genes.gr$gene_id[hits.upstream[idx]]
  
  # set type to 'ext_3UTR' if polyA site locating in a intergenic region but within 'threeUTR.mean.len' bp of the downstream of a specific gene
  idx = polyA$distance <= threeUTR.mean.len*2 & polyA$distance > 0
  polyA$type[idx] = 'ext_3UTR'
  
  
  # start annotation
  # annotate 3' UTR
  idx = which(polyA$distance == 0)
  if(length(idx)){
    hits.threeUTR = findOverlaps(polyA.gr[idx],threeUTRs.gr,ignore.strand=FALSE)
    polyA$type[idx[unique(hits.threeUTR@from)]] = '3UTR'
  }
  
  # annotate 5' UTR
  idx = which(polyA$distance == 0 & is.na(polyA$type))
  if(length(idx)){
    hits.fiveUTR = findOverlaps(polyA.gr[idx],fiveUTRs.gr,ignore.strand=FALSE)
    polyA$type[idx[unique(hits.fiveUTR@from)]] = '5UTR'
  }
  
  # annotate CDS
  idx = which(polyA$distance == 0 & is.na(polyA$type))
  if(length(idx)){
    hits.cds = findOverlaps(polyA.gr[idx],cds.gr,ignore.strand=FALSE)
    polyA$type[idx[unique(hits.cds@from)]] = 'CDS'
  }
  
  # annotate exon
  idx = which(polyA$distance == 0 & is.na(polyA$type))
  if(length(idx)){
    hits.exon = findOverlaps(polyA.gr[idx],exons.gr,ignore.strand=FALSE)
    polyA$type[idx[unique(hits.exon@from)]] = 'exon'
  }
  
  # annotate intron
  idx = which(polyA$distance == 0 & is.na(polyA$type))
  if(length(idx)){
    polyA$type[idx] = 'intron'
  }
  
  # annotate intergenic
  idx = which(is.na(polyA$type))
  if(length(idx)){
    polyA$type[idx] = 'intergenic'
  }
  
  # return result
  return(polyA)
}
