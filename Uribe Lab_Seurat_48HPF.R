# Set up the 48hpf data for processing. -----------------------------------
library(dplyr)
library(Seurat)
library(reticulate)
library(pryr)
rm(list = ls())

nc1.data <- Read10X(data.dir = "C:/Users/agh8/Desktop/RU_10XGenSC42/Sox10Enriched-48h_analysis/outs/filtered_gene_bc_matrices/DanioGRCz10/")

#nc1.data <- Read10X(data.dir = "C:/Users/agh8/Desktop/RU_10XGenSC42/Sox10Enriched-48h_analysis/outs/Ensb_filtered_gene_bc_matrices/DanioGRCz10/")

ProjectTitle <- "48HPF"
nc1 <- CreateSeuratObject(counts =nc1.data, project = ProjectTitle, min.cells = 3, min.features = 200)
nc1

pcnum <- 20
res <- 1.2
imagesize<- 4000
filepath<- "C:/Users/agh8/Desktop/Rice University/Uribe lab/Transcriptome/R programming/Seurat Plots/"

# nc1 Quality Control ---------------------------------------------------------
# Here we will filter cells that above 2500 UMI and under 200 UMI
#Also will remove cells with >5% Mitochondrial genes (indicates dying or damaged cell)

#Add Mitocohondrial QC to Seurat object
nc1[["percent.mt"]] <- PercentageFeatureSet(nc1, pattern = "^MT-")

#Add Cell Cycle QC to Seurat Object
nc1[["CellCycle"]] <- PercentageFeatureSet(nc1)

#Visualize QC metrics prior to thresholding
VlnPlot(nc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1<- FeatureScatter(nc1, feature1= "nCount_RNA", feature2 = "percent.mt")
plot2<- FeatureScatter(nc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
nc1 <- subset(nc1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
plot3<-FeatureScatter(nc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot2,plot3))

# nc1 Normilization, scaling, and feature selection of Data -------------------

#normalizaing the data
nc1 <- NormalizeData(nc1, normalization.method = "LogNormalize", scale.factor = 10000)


#find variable features (eg. id most significant dimensions)
nc1 <- FindVariableFeatures(nc1, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nc1),10)

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(nc1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

#scaling data
all.genes <- rownames(nc1)
nc1 <- ScaleData(nc1, features = all.genes)

#Remove unwanted sources of variation
# nc1 <- ScaleData(nc1, vars.to.regress = "percent.mt") or refer to the Seurat V3 vignette on sctransform

#linear dimensional reduction (i.e. find the PCs)
nc1 <- RunPCA(nc1, features = VariableFeatures(object =nc1))
#Print the  top PC gene lists
print(nc1[["pca"]],dims = 1:pcnum, nfeatures = 5)

nc1qc.pcaVis <- VizDimLoadings(nc1,dims = 1:pcnum, reduction = "pca")
fname = paste(filepath,"/QC PLOTS/","48hpf_PCAViz_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(VizDimLoadings(nc1,dims = 1:pcnum, reduction = "pca"))
dev.off()


nc1qc.pcaHM %<a-%{ DimHeatmap(nc1, dims = 1:pcnum, cells = 100, balanced = TRUE, fast = TRUE)}
fname = paste(filepath,"/QC PLOTS/","48hpf_PCAHM_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize*4, height = imagesize*5,units = "px",res = 300,pointsize = 24)
nc1qc.pcaHM
p1.save <-recordPlot()
dev.off()

nc1 <- JackStraw(nc1, reduction = "pca", num.replicate = 100,dims = pcnum)
nc1 <- ScoreJackStraw(nc1, dims = 1:pcnum)
nc1qc.jk<-JackStrawPlot(nc1, dims = 1:pcnum)
fname = paste(filepath,"/QC PLOTS/","48hpf_JackStraw_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc1qc.jk)
dev.off()

nc1qc.EP<-ElbowPlot(nc1, ndims = pcnum)
fname = paste(filepath,"/QC PLOTS/","48hpf_ElbowPlot_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc1qc.EP)
dev.off()


# nc1 Clustering of Cells -----------------------------------------------------


nc1 <- FindNeighbors(nc1, dims = 1:pcnum)
nc1 <- FindClusters(nc1, resolution = res)

head(Idents(nc1),5)


# nc1 Run TSEN ----------------------------------------------------------------
nc1 <- RunTSNE(object = nc1, dims = 1:pcnum, do.fast = TRUE )
nc1_tsne <- TSNEPlot(nc1, label = FALSE)
fname = paste(filepath,"48hpf_tsne_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc1_tsne)
dev.off()

nc1_tsne <- TSNEPlot(nc1, label = TRUE)
fname = paste(filepath,"48hpf_tsne_PC",pcnum,"+res",res,"_Labled.tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(nc1_tsne)
dev.off()


#proposed cluster labels
new.cluster.ids <- c("cl0:Autonomic neuronal progenitor",
                     "cl1:Mesenchymal progenitor 1",
                     "cl2:Chondrogenic mesenchyme proliferative 1",
                     "cl3:Chondrogenic mesenchyme proliferative 2",
                     "cl4:Chondrogenic mesenchyme stem-like",
                     "cl5:Migratory/enteric neural crest"    ,        
                     "cl6:Mesenchymal differentiating 1"     ,
                     "cl7:CNS neurons"  ,
                     "cl8:Melanophore"  ,
                     "cl9:Chondrogenic migratory mesenchyme",
                     "cl10:Mesenchymal differentiating 2"         ,       
                     "cl11:Mesenchymal migratory",
                     "cl12:Chondrogenic mesenchyme proliferative 3",
                     "cl13:Neural progenitor",
                     "cl14:Fin bud" ,
                     "cl15:Peripherial glial progenitor" ,
                     "cl16:Otic epithelium"           ,              
                     "cl17:Sensory Neuronal progenitor",             
                     "cl18:Muscle")
names(new.cluster.ids) <- levels(nc1)
nc1 <- RenameIdents(nc1,new.cluster.ids)
nc1$CellType <- Idents(nc1)
FP<-DimPlot(nc1,reduction = "tsne", label = TRUE, pt.size = 0.5, group.by = "CellType", repel = TRUE)+NoLegend()
fname = paste(filepath,ProjectTitle,"_tsne_PC",pcnum,"+res",res,"_CellTypeLabels.tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(FP)
dev.off()


#Print Dendrogram
nc1<- BuildClusterTree(object = nc1)
PlotClusterTree(nc1)
fname = paste(filepath,ProjectTitle,"_",pcnum,"+res",res,"ClusterDendrogram.tiff",sep = "")
tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
print(PlotClusterTree(nc1))
dev.off()



#saveRDS(nc1, file = "~/Desktop/UribeSeurat_48hpf_scRNA.rds")

nc1.tsne_markers <- FindAllMarkers(object = nc1, only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
nc1.tsne_markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.csv(nc1.tsne_markers, "48hpf_tsne_markers.csv")

# # cluster1.markers <- FindMarkers(object = nc1, ident.1 = 0, thresh.use = 0.25,
# #                                 test.use = "roc", only.pos = TRUE)
# VlnPlot(nc1, features = c("dlx4a", "lect1"))
# VlnPlot(nc1, features = c("sox10", "lect1"), slot = "counts", log = TRUE)
# 
# FeaturePlot(nc1, features = c("sox10", "lect1", "dlx4a", "fabp7a",
#                                             "phox2bb", "foxc1a", "mitfa", "foxd3", "mycn"))
# 
top10 <- nc1.tsne_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10HM <- DoHeatmap(nc1, features = top10$gene) + NoLegend()
fname = paste(filepath,ProjectTitle,"_Top10MarkerHeatMap_PC",pcnum,"+res",res,".tiff",sep = "")
tiff(filename=fname,width = imagesize+2000, height = imagesize+1000,units = "px",res = 300,pointsize = 14)
print(top10HM)
dev.off()

# 
#write.csv(top100, "48hpf_tsne_top100markers.csv")
