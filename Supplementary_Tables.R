#Export supplementary tables

library(Seurat)
library(Signac)
library(dplyr)
library(readr)
library(openxlsx)
set.seed(1234)

setwd("/working_dir")

# Supplementary Table 2 - table summary from DESeq2 results
directory <- "tables/DEGs"
csv_files <- list.files(directory, pattern="*.csv", full.names = TRUE)
wb <- createWorkbook()
for (file in csv_files) {
  data <- read_csv(file)
  sheet_name <- tools::file_path_sans_ext(basename(file))
  sheet_name <- gsub("_vs_control_signif_genes", "", sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
}
saveWorkbook(wb, "Supplementary_Data_Table2.xlsx", overwrite = TRUE)

#supplementary Table 3 (DARs)

library(Seurat)
library(Signac)
set.seed(1234)

snATAC_C9ALSFTLD <- readRDS(file="objects/snATAC_C9ALSFTLDvCTRL.RDS")
snATAC_C9ALSFTLD

fragments<-CreateFragmentObject(path="objects/fragments1.tsv.gz",cells=colnames(snATAC_C9ALSFTLD))
Fragments(snATAC_C9ALSFTLD)<-NULL
Fragments(snATAC_C9ALSFTLD)<-fragments

table(snATAC_C9ALSFTLD$subclass_DE)

class_types<-c("Oligodendrocytes","OPC","Astrocytes","Micro-PVM","Excitatory","Inhibitory")
Idents(snATAC_C9ALSFTLD)<-'class_clusters'
DefaultAssay(snATAC_C9ALSFTLD) <- 'ATAC'
combined_df <- data.frame()
for(i in 1:length(class_types)){
  
  cur_celltype <- class_types[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(snATAC_C9together, class_clusters == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers1 <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c("C9ALSFTLD"),
                         ident.2=c("control"),
                         test.use = "LR",
                         latent.vars = c("sex","peak_region_fragments","on_target_fragments","nCount_ATAC"),
                         logfc.threshold = 0.1,
                         #mean.fxn=rowMeans,
                         min.pct = 0.1,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='ATAC')
    markers<-subset(markers1,p_val_adj<0.01)
    cur_markers<-rownames(markers)
	closest_genes<-ClosestFeature(snATAC_C9ALSFTLD,regions=cur_markers)
	markers$genes<-closest_genes$gene_name
	markers$gene_biotype<-closest_genes$gene_biotype
	markers$type<-closest_genes$type
	markers$closest_region<-closest_genes$closest_region
	markers$query_region<-closest_genes$query_region
  markers$distance<-closest_genes$distance
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$class_clusters <- cur_celltype
  combined_df <- rbind(combined_df, markers)
}

snATAC_sALSnoFTLD <- readRDS(file="objects/snATAC_sALSnoFTLDvCTRL.RDS")
snATAC_sALSnoFTLD

fragments<-CreateFragmentObject(path="objects/fragments1.tsv.gz",cells=colnames(snATAC_sALSnoFTLD))
Fragments(snATAC_sALSnoFTLD)<-NULL
Fragments(snATAC_sALSnoFTLD)<-fragments

table(snATAC_sALSnoFTLD$subclass_DE)


class_types<-c("Oligodendrocytes","OPC","Astrocytes","Micro-PVM","Excitatory","Inhibitory")
Idents(snATAC_sALSnoFTLD)<-'class_clusters'
DefaultAssay(snATAC_sALSnoFTLD) <- 'ATAC'
combined_df <- data.frame()
for(i in 1:length(class_types)){
  
  cur_celltype <- class_types[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(snATAC_sALSnoFTLD, class_clusters == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers1 <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c("sALSnoFTLD"),
                         ident.2=c("control"),
                         test.use = "LR",
                         latent.vars = c("sex","peak_region_fragments","on_target_fragments","nCount_ATAC"),
                         logfc.threshold = 0.1,
                         #mean.fxn=rowMeans,
                         min.pct = 0.1,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='ATAC')
    markers<-subset(markers1,p_val_adj<0.01)
    cur_markers<-rownames(markers)
	  closest_genes<-ClosestFeature(snATAC_sALSnoFTLD,regions=cur_markers)
	  markers$genes<-closest_genes$gene_name
	  markers$gene_biotype<-closest_genes$gene_biotype
	  markers$type<-closest_genes$type
	  markers$closest_region<-closest_genes$closest_region
	  markers$query_region<-closest_genes$query_region
    markers$distance<-closest_genes$distance
    markers$diff <- markers$pct.1 - markers$pct.2
    markers$class_clusters <- cur_celltype
    combined_df <- rbind(combined_df, markers)
}

class.sALSnoFTLD.LR <- combined_df

save(class.sALSnoFTLD.LR,file="tables/class_sALSnoFTLD_LR.rda")
write.csv(class.sALSnoFTLD.LR, file='tables/class_sALSnoFTLD_LR.csv', quote=F)

snATAC_class.C9ALSFTLD.LR <- combined_df

C9ALS_exc_upper_DAR <- subset(diagnoses_exc_nodeep, diagnoses_exc_nodeep$diagnosis == "C9ALS")
sALS_exc_upper_DAR <- subset(diagnoses_exc_nodeep, diagnoses_exc_nodeep$diagnosis == "sALS")

C9ALS_exc_upper_DAR$diagnosis <- NULL
C9ALS_exc_upper_DAR$class_clusters <- "Excitatory_upper"
C9ALS_DAR <- rbind(class.C9together.LR,C9ALS_exc_upper_DAR)
C9ALS_DAR$diagnosis <- "C9ALS"
C9ALS_DAR$cluster <- paste0(C9ALS_DAR$diagnosis, "_", C9ALS_DAR$class_clusters)

sALS_exc_upper_DAR$diagnosis <- NULL
sALS_exc_upper_DAR$class_clusters <- "Excitatory_upper"
sALS_DAR <- rbind(snATAC_class_sALSnoFTLD_LR,sALS_exc_upper_DAR)
sALS_DAR$diagnosis <- "sALS"
sALS_DAR$cluster <- paste0(sALS_DAR$diagnosis, "_", sALS_DAR$class_clusters)
DARs <- rbind(C9ALS_DAR, sALS_DAR)

idents_DAR <- unique(DARs$cluster)
idents_DAR

wb_DAR <- createWorkbook()

for (x in idents_DAR) {
  output <- which(DARs$cluster == x)
  addWorksheet(wb_DAR, output)
  writeData(wb_DAR, x, output)
}
saveWorkbook(wb_DAR, file = "Supplementary_Data_Table3.xlsx", overwrite=T)

#Supplementary Table 4 put together manually from rGREAT output

#Supplementary Table 5
#chromVAR FAM motifs
load("/mnt/WORKHORSE/C9ALSFTLD_multiome/tables/chromvar/FAM_class_motif_markers.rda")
rownames(class_clusters_motif_markers) <- NULL
head(class_clusters_motif_markers)
idents_motifs <- unique(class_clusters_motif_markers$cluster)

wb_motifs <- createWorkbook()

for (x in idents_motifs) {
  TF_cluster <- class_clusters_motif_markers[class_clusters_motif_markers$cluster == x,]
  addWorksheet(wb_motifs,x)
  writeData(wb_motifs, x, TF_cluster)
}
saveWorkbook(wb_motifs, file = "Supplementary_Data_Table_5.xlsx", overwrite = T)

#Supplementary Table 6
#chromVAR motifs
load("/mnt/WORKHORSE/C9ALSFTLD_multiome/tables/snATAC_class_C9together_cutoff.rda")
load("/mnt/WORKHORSE/C9ALSFTLD_multiome/tables/snATAC_class_sALSnoFTLD_cutoff.rda")
load("/mnt/WORKHORSE/C9ALSFTLD_multiome/tables/snATAC_exc_nodeep_cutoff.rda")

C9ALS_exc_upper <- subset(exc.nodeep.chromvar.cutoff, exc.nodeep.chromvar.cutoff$diagnosis == "C9ALS")
sALS_exc_upper <- subset(exc.nodeep.chromvar.cutoff, exc.nodeep.chromvar.cutoff$diagnosis == "sALS")

C9ALS_exc_upper$diagnosis <- NULL
C9ALS_exc_upper$class_clusters <- "Excitatory_upper"
C9ALS_DE_motifs <- rbind(C9ALS.cutoff.df,C9ALS_exc_upper)
C9ALS_DE_motifs$diagnosis <- "C9ALS"
C9ALS_DE_motifs$cluster <- paste0(C9ALS_DE_motifs$diagnosis, "_", C9ALS_DE_motifs$class_clusters)

sALS_exc_upper$diagnosis <- NULL
sALS_exc_upper$class_clusters <- "Excitatory_upper"
sALS_DE_motifs <- rbind(sALS.cutoff.df,sALS_exc_upper)
sALS_DE_motifs$diagnosis <- "sALS"
sALS_DE_motifs$cluster <- paste0(sALS_DE_motifs$diagnosis, "_", sALS_DE_motifs$class_clusters)
DE_motifs <- rbind(C9ALS_DE_motifs, sALS_DE_motifs)
rownames(DE_motifs) <- NULL
idents_DE_motifs <- unique(DE_motifs$cluster)

wb_DE_motifs <- createWorkbook()

for (x in idents_DE_motifs) {
  DE_cluster <- DE_motifs[DE_motifs$cluster == x,]
  addWorksheet(wb_DE_motifs,x)
  writeData(wb_DE_motifs, x, DE_cluster)
}

saveWorkbook(wb_DE_motifs, file = "Supplementary_Data_Table_6.xlsx", overwrite = T)