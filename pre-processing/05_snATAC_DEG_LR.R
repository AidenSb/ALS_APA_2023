#Calculate DARs for Figures 2-3
library(Seurat)
library(Signac)
set.seed(1234)

snATAC_C9together <- readRDS(file="objects/snATAC_C9togethervCTRL.RDS")
snATAC_C9together

fragments<-CreateFragmentObject(path="objects/fragments.tsv.gz",cells=colnames(snATAC_C9together))
Fragments(snATAC_C9together)<-NULL
Fragments(snATAC_C9together)<-fragments

table(snATAC_C9together$subclass_DE)

class_types<-c("Oligodendrocytes","OPC","Astrocytes","Micro-PVM","Excitatory","Inhibitory")
Idents(snATAC_C9together)<-'class_clusters'
DefaultAssay(snATAC_C9together) <- 'ATAC'
combined_df <- data.frame()
for(i in 1:length(class_types)){
  
  cur_celltype <- class_types[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(snATAC_C9together, class_clusters == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$C9together
  print(unique(Idents(cur_seurat_obj)))
  
  markers1 <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c("C9ALS"),
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
	closest_genes<-ClosestFeature(snATAC_C9together,regions=cur_markers)
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

class.C9together.LR <- combined_df

save(class.C9together.LR,file="tables/class_C9together_LR.rda")
write.csv(class.C9together.LR, file='tables/class_C9together_LR.csv', quote=F)

library(Seurat)
library(Signac)
set.seed(1234)

snATAC_sALSnoFTLD <- readRDS(file="objects/snATAC_sALSnoFTLDvCTRL.RDS")
snATAC_sALSnoFTLD

fragments<-CreateFragmentObject(path="objects/fragments.tsv.gz",cells=colnames(snATAC_sALSnoFTLD))
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