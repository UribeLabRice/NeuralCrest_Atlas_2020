# Set up the 69hpf data for processing. -----------------------------------
library(dplyr)
library(Seurat)
library(reticulate)
rm(list = ls())

nc2.data <- Read10X(data.dir = "C:/Users/agh8/Desktop/RU_10XGenSC42/Sox10Enriched-69h_analysis/outs/filtered_gene_bc_matrices/DanioGRCz10/")
#nc2.data <- Read10X(data.dir = "C:/Users/agh8/Desktop/RU_10XGenSC42/Sox10Enriched-69h_analysis/outs/Ensb_filtered_gene_bc_matrices/DanioGRCz10/")
ProjectTitle <- "69HPF"
nc2 <- CreateSeuratObject(counts = nc2.data, project = ProjectTitle, min.cells = 3, min.features = 200)
nc2

pcnum <- 20
res <- 1.8
imagesize<- 4000
filepath<- "C:/Users/agh8/Desktop/Rice University/Uribe lab/Transcriptome/R programming/Seurat Plots/"

# NC2 Quality Control ---------------------------------------------------------
# Here we will filter cells that above 2500 UMI and under 200 UMI
#Also will remove cells with >5% Mitochondrial genes (indicates dying or damaged cell)

#Add Mitocohondrial QC to Seurat object
nc2[["percent.mt"]] <- PercentageFeatureSet(nc2, pattern = "^MT-")

#Add Cell Cycle QC to Seurat Object
nc2[["CellCycle"]] <- PercentageFeatureSet(nc2)

#Visualize QC metrics prior to thresholding
VlnPlot(nc2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1<- FeatureScatter(nc2, feature1= "nCount_RNA", feature2 = "percent.mt")
plot2<- FeatureScatter(nc2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
nc2 <- subset(nc2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# NC2 Normilization, scaling, and feature selection of Data -------------------

#normalizaing the data
nc2 <- NormalizeData(nc2, normalization.method = "LogNormalize", scale.factor = 10000)


#find variable features (eg. id most significant dimensions)
nc2 <- FindVariableFeatures(nc2, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nc2),10)

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(nc2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

#scaling data
all.genes <- rownames(nc2)
nc2 <- ScaleData(nc2, features = all.genes)

#Remove unwanted sources of variation
# nc2 <- ScaleData(nc2, vars.to.regress = "percent.mt") or refer to the Seurat V3 vignette on sctransform

#linear dimensional reduction (i.e. find the PCs)
nc2 <- RunPCA(nc2, features = VariableFeatures(object =nc2))
#Print the  top PC gene lists
print(nc2[["pca"]],dims = 1:pcnum, nfeatures = 5)

nc2qc.pcaVis <- VizDimLoadings(nc2,dims = 1:pcnum, reduction = "pca")
fname = paste(filepath,"QC PLOTS/","69hpf_PCAViz_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(VizDimLoadings(nc2,dims = 1:pcnum, reduction = "pca"))
dev.off()


nc2qc.pcaHM %<a-%{ DimHeatmap(nc2, dims = 1:pcnum, cells = 100, balanced = TRUE, fast = TRUE)}
fname = paste(filepath,"/QC PLOTS/","69hpf_PCAHM_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize*4, height = imagesize*5,units = "px",res = 300,pointsize = 24)
nc2qc.pcaHM
p1.save <-recordPlot()
dev.off()

nc2 <- JackStraw(nc2, reduction = "pca", num.replicate = 100,dims = pcnum)
nc2 <- ScoreJackStraw(nc2, dims = 1:pcnum)
nc2qc.jk<-JackStrawPlot(nc2, dims = 1:pcnum)
fname = paste(filepath,"QC PLOTS/","69hpf_JackStraw_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc2qc.jk)
dev.off()

nc2qc.EP<-ElbowPlot(nc2, ndims = pcnum)
fname = paste(filepath,"QC PLOTS/","69hpf_ElbowPlot_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc2qc.EP)
dev.off()


# NC2 Clustering of Cells -----------------------------------------------------


nc2 <- FindNeighbors(nc2, dims = 1:pcnum)
nc2 <- FindClusters(nc2, resolution = res)

head(Idents(nc2),5)


# NC2 Run TSEN ----------------------------------------------------------------
nc2 <- RunTSNE(object = nc2, dims = 1:pcnum, do.fast = TRUE )
nc2_tsne <- TSNEPlot(nc2, label = FALSE)
fname = paste(filepath,"69hpf_tsne_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc2_tsne)
dev.off()

nc2_tsne <- TSNEPlot(nc2, label = TRUE)
fname = paste(filepath,"69hpf_tsne_PC",pcnum,"+res",res,"_Labled.tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc2_tsne)
dev.off()


#proposed cluster labels
new.cluster.ids <- (c("cl0:Chondrogenic differentiating mesenchyme 1",
                      "cl1:Fin bud",
                      "cl2:Mesenchyme Migratory",                   
                      "cl3:Neural progenitor" ,
                      "cl4:Melanophore",
                      "cl5:Sympatho-enteric progenitor" ,
                      "cl6:Chondrogenic differentiating mesenchyme 2" ,
                      "cl7:Fin bud",
                      "cl8:Chondrogenic proliferative mesenchyme",
                      "cl9:Pigmented muscle" ,
                      "cl10:Neural progenitor",
                      "cl11:Otic epithelium" ,
                      "cl12:Enteric neuron"   ,                       
                      "cl13:Pigment Progenitor",
                      "cl14:Schwann Cells" ,
                      "cl15:Xanthophore",                              
                      "cl16:Iridiophore" ,
                      "cl17:Differentiating mesenchyme 1",
                      "cl18:Proliferative melanophores",      
                      "cl19:Muscle"                     ,              
                      "cl20:Chondrogenic mesenchyme migratory",         
                      "cl21:Unidentified"                      ,       
                      "cl22:Differentiating mesenchyme 2"))
names(new.cluster.ids) <- levels(nc2)
nc2 <- RenameIdents(nc2,new.cluster.ids)
nc2$CellType <- Idents(nc2)
FP<-DimPlot(nc2,reduction = "tsne", label = TRUE, pt.size = 0.5, repel = TRUE, group.by = "CellType")+NoLegend()
fname = paste(filepath,ProjectTitle,"_tsne_PC",pcnum,"+res",res,"_CellTypeLabels.tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(FP)
dev.off()


#Print Dendrogram
nc2<- BuildClusterTree(object = nc2)
PlotClusterTree(nc2)
fname = paste(filepath,ProjectTitle,"_",pcnum,"+res",res,"ClusterDendrogram.tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(PlotClusterTree(nc2))
dev.off()

nc2.tsne_markers <- FindAllMarkers(object = nc2, only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
nc2.tsne_markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

top10 <- nc2.tsne_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10HM <- DoHeatmap(nc2,group.by = "seurat_clusters") + NoLegend()
fname = paste(filepath,ProjectTitle,"_Top10MarkerHeatMap_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize+2000, height = imagesize+1000,units = "px",res = 300,pointsize = 14)
print(top10HM)
dev.off()



# write.csv(nc2.tsne_markers, "69hpf_tsne_markers.csv")

# # cluster1.markers <- FindMarkers(object = nc2, ident.1 = 0, thresh.use = 0.25,
# #                                 test.use = "roc", only.pos = TRUE)
# VlnPlot(nc2, features = c("dlx4a", "lect1"))
# VlnPlot(nc2, features = c("sox10", "lect1"), slot = "counts", log = TRUE)
# 
# FeaturePlot(nc2, features = c("sox10", "lect1", "dlx4a", "fabp7a",
#                                             "phox2bb", "foxc1a", "mitfa", "foxd3", "mycn"))
# 
# top100 <- nc2.tsne_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
# DoHeatmap(nc2, features = top10$gene) + NoLegend()
# 
#write.csv(top100, "69hpf_tsne_top100markers.csv")

# NC2 UMAP --------------------------------------------------------------------

nc2 <- RunUMAP(nc2,dims = 1:pcnum)
nc2_umap <- DimPlot(nc2,reduction = "umap",label = FALSE)
fname = paste(filepath,"69hpf_umap_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc2_umap)
dev.off()

nc2_umap <- DimPlot(nc2,reduction = "umap",label = TRUE,repel = TRUE, group.by = "CellType")
fname = paste(filepath,"69hpf_umap_PC",pcnum,"+res",res,"._Labeled.tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc2_umap)
dev.off()

FP<-DimPlot(nc2,reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE, group.by = "CellType")+NoLegend()
fname = paste(filepath,ProjectTitle,"_umap_PC",pcnum,"+res",res,"_CellTypeLabels.tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(FP)
dev.off()

#find all unique cluster biomarkers
#nc2.umap_markers <- FindAllMarkers(nc2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#nc2.umap_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#write.csv(nc2.umap_markers, "69hpf_umap_markers.csv")
# 
# top10<- nc2.markers %>% group_by(cluster) %>% top_n(n=pcnum, wt = avg_logFC)
# DoHeatmap(nc2,features = top10$gene) + NoLegend()

saveRDS(nc2, file = paste(filepath,ProjectTitle,"_PC",pcnum,"+res",res,".rds",sep = ""))
