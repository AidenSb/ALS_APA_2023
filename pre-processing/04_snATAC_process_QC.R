#Pre-procesing of snATAC-Seq data
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(monocle3)
library(cicero)
library(SeuratWrappers)
set.seed(1234)

setwd("/ALS_multiome")


# Define sample names and their corresponding metadata
sample_names <- c(paste0("CTRL", 1:6), paste0("sALSnoFTLD", 1:8),paste0("C9ALSFTLD", 1:3),paste0("C9ALSFTLD", 5:7), paste0("C9ALSnoFTLD", 1:3))
snATAC_diagnosis <- c(rep("control", 4), rep("C9ALS", 9), rep("sALS", 8))
snATAC_sex<-c("female","male","female","male","male","female","male","female","male","female","male","male",
              "male","male","male","female","male","male","female","female","female"
              )
snATAC_age<-c(50,59,72,72,59,57,65,88,43,72,37,50,60,60,58,59,72,57,71,66,88)

# Define the paths to your peaks.bed files
peak_files <- c(
  "/ALS_snATAC/CTRL1_snATAC/outs/peaks.bed",
  "/ALS_snATAC/CTRL2_snATAC/outs/peaks.bed",
  "/ALS_snATAC/CTRL3_snATAC/outs/peaks.bed",
  "/ALS_snATAC/CTRL4_snATAC/outs/peaks.bed",
  "/ALS_snATAC/sALSnoFTLD1_snATAC/outs/peaks.bed",
  "/ALS_snATAC/sALSnoFTLD2_snATAC/outs/peaks.bed",
  "/ALS_snATAC/sALSnoFTLD3_snATAC/outs/peaks.bed",
  "/ALS_snATAC/sALSnoFTLD4_snATAC/outs/peaks.bed",
  "/ALS_snATAC/sALSnoFTLD5_snATAC/outs/peaks.bed",
  "/ALS_snATAC/sALSnoFTLD6_snATAC/outs/peaks.bed",
  "/ALS_snATAC/sALSnoFTLD7_snATAC/outs/peaks.bed",
  "/ALS_snATAC/sALSnoFTLD8_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSFTLD1_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSFTLD2_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSFTLD3_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSFTLD5_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSFTLD6_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSFTLD7_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSnoFTLD1_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSnoFTLD2_snATAC/outs/peaks.bed",
  "/ALS_snATAC/C9ALSnoFTLD3_snATAC/outs/peaks.bed"
)

# Initialize an empty list to store GRanges objects
granges_list <- list()

# Loop through each file, read the peaks, and convert to GRanges
for (file_path in peak_files) {
  peaks <- read.table(file = file_path, col.names = c("chr", "start", "end"))
  granges_list[[file_path]] <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
}

# Combine all GRanges objects into a unified set
combined.peaks <- reduce(granges_list)

# Filter out peaks based on length criteria
filtered.peaks <- combined.peaks[width(combined.peaks) < 10000 & width(combined.peaks) > 20]

# Display the combined and filtered peaks
filtered.peaks

#Now 
sample_paths <- paste0("/ALS_snATAC/", sample_names, "_snATAC/outs/")

# Initialize list to store Seurat objects
seurat_objects <- list()

# Loop over samples
for (i in seq_along(sample_names)) {
  # Construct file paths
  peaks_file <- paste0(sample_paths[i], "peaks.bed")
  metadata_file <- paste0(sample_paths[i], "singlecell.csv")
  fragments_file <- paste0(sample_paths[i], "fragments.tsv.gz")
  
  # Read peaks and metadata
  peaks <- read.table(file = peaks_file, col.names = c("chr", "start", "end"))
  metadata <- read.table(file = metadata_file, sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)[-1,]
  metadata <- metadata[metadata$passed_filters > 500, ]
  
  # Convert DF to Gragges and use filtered.peaks
  granges <- makeGRangesFromDataFrame(filtered.peaks)
    
  # Create Signac Fragment Object
  frags <- CreateFragmentObject(path = fragments_file, cells = rownames(metadata))
  
  # Create FeatureMatrix as input
  counts <- FeatureMatrix(fragments = frags, features = combined.peaks, cells = rownames(metadata))
  
  # Combine to create Seurat Object
  assay <- CreateChromatinAssay(counts, fragments = frags)
  seurat_obj <- CreateSeuratObject(assays = list(ATAC = assay), meta.data = metadata, project = sample_names[i])
  
  # Store objects in list
  seurat_objects[[sample_names[i]]] <- seurat_obj
}


#Merge all samples for further processing
unint_ATAC <- Reduce(function(x, y) merge(x, y, add.cell.ids = sample_names, project = "ALS_multiome"), seurat_objects)


#Save file for future steps
saveRDS(unint_ATAC, file="/mnt/WORKHORSE/C9ALSFTLD_multiome/objects/snATAC/unint_ATAC.RDS")

unint_ATAC[["ATAC"]]
DefaultAssay(unint_ATAC) <- "ATAC"

granges(unint_ATAC)
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(unint_ATAC) <- annotations

# compute nucleosome signal score per cell
unint_ATAC <- NucleosomeSignal(object = unint_ATAC)

# compute TSS enrichment score per cell
unint_ATAC <- TSSEnrichment(object = unint_ATAC, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
unint_ATAC$pct_reads_in_peaks <- unint_ATAC$peak_region_fragments / unint_ATAC$passed_filters * 100
unint_ATAC$blacklist_ratio <- unint_ATAC$blacklist_region_fragments / unint_ATAC$peak_region_fragments
unint_ATAC$high.tss <- ifelse(unint_ATAC$TSS.enrichment > 2, 'High', 'Low')
unint_ATAC$nucleosome_group <- ifelse(unint_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

saveRDS(unint_ATAC,file="objects/snATAC/snATAC_unint_preQC.RDS")

#apply appropriate cut-offs based on violin plots from above QC data
unint_QC <- subset(
  x = unint_ATAC,
  subset = peak_region_fragments > 2000 &
    peak_region_fragments < 30000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    TSS.enrichment < 20
)
unint_QC
saveRDS(unint_QC,file="objects/snATAC/snATAC_unint_postQC.RDS")
rm(unint_ATAC)
gc()

unint_QC[["ATAC"]]
unint_QC <- FindVariableFeatures(unint_QC)
unint_QC <- RunTFIDF(unint_QC)
unint_QC <- FindTopFeatures(unint_QC, min.cutoff = 50)
unint_QC <- RunSVD(unint_QC)
unint_QC <- RunUMAP(unint_QC, reduction = 'lsi', dims = 2:30)
gc()

#perform batch correction using Harmony (Seurat wrapper) - batch effect evident on UMAP
snATAC <- RunHarmony(
  object = unint_QC,
  group.by.vars = 'sample',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)
gc()
snATAC <- RunUMAP(snATAC, reduction = "harmony", dims = 2:30)
snATAC <- FindNeighbors(object = snATAC, reduction = 'harmony', dims = 2:30)
snATAC <- FindClusters(object = snATAC, verbose = FALSE, algorithm = 4, method ="igraph")

#calculate Gene Activity Matrix
gene.activities <- GeneActivity(snATAC)
snATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)

snATAC <- NormalizeData(
  object = snATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(snATAC$nCount_RNA)
)

#transfer labels from snRNA
snRNA <- readRDS(file="objects/snRNA/snRNA_final.RDS")
snATAC <- FindVariableFeatures(snATAC, nfeatures=3000)

transfer.anchors <- FindTransferAnchors(
  reference = snRNA,
  query = snATAC,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = snRNA$subclass_clusters,
  weight.reduction = snATAC[['harmony']],
  dims = 2:30
)
snATAC <- AddMetaData(object = snATAC, metadata = predicted.labels)
Idents(snATAC) <- "predicted.id"
snATAC$subclass_clusters <- Idents(snATAC)
table(snATAC$subclass_clusters)

#Cutoff label transfer
snATAC <- snATAC[,snATAC$prediction.score.max >= 0.5]

DefaultAssay(snATAC) <- "ATAC"
Idents(snATAC) <- "subclass_clusters"
peaks <- CallPeaks(snATAC,
                   macs2.path="ENV/bin/macs2",
                   group.by="subclass_clusters")
peaks2 <- keepStandardChromosomes(peaks,
                                  pruning.mode="coarse")
peaks2 <- subsetByOverlaps(x=peaks2,
                             ranges=blacklist_hg38_unified,
                             invert=T)
macs2_counts <- FeatureMatrix(fragments=Fragments(snATAC), features=peaks2, cells=colnames(snATAC))

#save in case crash
saveRDS(macs2_counts, file = "objects/snATAC/snATAC_macs2counts.RDS")
saveRDS(peaks2, file = "objects/snATAC/snATAC_MACS2peaks_only.RDS")

DefaultAssay(snATAC) <- "ATAC"
Idents(snATAC) <- "subclass_clusters"
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

snATAC[["ATAC"]] <- NULL
snATAC[["ATAC"]] <- CreateChromatinAssay(counts = macs2_counts,
                                         fragments = Fragments(snATAC),
                                         annotation = annotations
)

#Run ChromVAR to calculate motif activity
cisBP<-readRDS(file="cisBP_human_pfms_2021.rds")

snATAC <- AddMotifs(snATAC,
                    genome=BSgenome.Hsapiens.UCSC.hg38,
                    pfm=cisBP)
snATAC <- RunChromVAR(object = snATAC,
                      genome = BSgenome.Hsapiens.UCSC.hg38)
snATAC <- RegionStats(snATAC,
                      genome=BSgenome.Hsapiens.UCSC.hg38)


genome<-seqlengths(snATAC)
genome.df<-data.frame("chr"=names(genome),"length"=genome)

snATAC.cds <- as.cell_data_set(snATAC)
snATAC.cic <- make_cicero_cds(snATAC.cds, reduced_coordinates = reducedDims(snATAC.cds)$UMAP)
conns <- run_cicero(snATAC.cic, genomic_coords = genome.df, sample_num=100)
head(conns)

ccans<-generate_ccans(conns)
head(ccans)

save(conns, file = "objects/snATAC/cicero/snATAC_cicero_conns.rda")
save(ccans, file = "objects/snATAC/cicero/snATAC_cicero_ccans.rda")

links <- ConnectionsToLinks(conns=conns, ccans=ccans)

saveRDS(links, file = "objects/snATAC/cicero/snATAC_cicero_links_unadded.RDS")

Links(snATAC) <- links

saveRDS(snATAC, file = "objects/snATAC/cicero/snATAC_cicero_links.RDS")