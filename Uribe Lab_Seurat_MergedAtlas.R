#Intergrated Data Analysis

# Inputs for Data processing -----------------------------------
library(dplyr)
library(Seurat)
library(reticulate)
library(cowplot)
rm(list = ls())

pcnum<- 20
res <- 1.2
imagesize<- 4000
filepath<- "C:/Users/agh8/Desktop/Rice University/Uribe lab/Transcriptome/R programming/Seurat Plots/"

ProjectTitle <- "Merge"


nc1 <- readRDS(file = paste(filepath,"48HPF","_PC",20,"+res",1.2,".rds", sep = ""))
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
nc1$CellType <- paste("48HPF",nc1$CellType,sep = "")


nc2 <- readRDS(file = paste(filepath,"69HPF","_PC",20,"+res",1.2,".rds", sep = ""))
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
nc2$CellType <- paste("69HPF",nc2$CellType,sep = "")

nc3.anchors <- FindIntegrationAnchors(object.list = list(nc1,nc2), dims = 1:pcnum)
nc3 <- IntegrateData(anchorset = nc3.anchors, dims = 1:pcnum)


#scaling data?
nc3 <- ScaleData(nc3)

#linear dimensional reduction (i.e. find the PCs)
nc3 <- RunPCA(nc3, features = VariableFeatures(object =nc3), npcs = pcnum)
nc3 <- RunUMAP(nc3, reduction = "pca", dims = 1:pcnum)
nc3 <- RunTSNE(nc3, reduction = "pca", dims = 1:pcnum)
nc3 <- FindNeighbors(nc3, reduction = "pca", dims = 1:pcnum)
nc3 <- FindClusters(nc3, resolution = res)


DefaultAssay(nc3) <- "integrated"

plottype <- c("tsne","umap")
printpath<- "C:/Users/agh8/Desktop/Rice University/Uribe lab/Transcriptome/R programming/Seurat Plots/Merge2/"
saveRDS(nc3, file = paste(printpath,ProjectTitle,"_PC",pcnum,"+res",res,".rds",sep = ""))



plottype <- c("tsne","umap")
printpath<- "C:/Users/agh8/Desktop/Rice University/Uribe lab/Transcriptome/R programming/Seurat Plots/Merge2/"
nc3 <- readRDS(file = paste(printpath,ProjectTitle,"_PC",pcnum,"+res",res,".rds",sep = ""))


# Neural Sub out ---------------------------------------------------------
Idents(nc3) <-"CellType"
nc3.Neural <- subset(nc3, idents = c("48HPFcl5:Migratory/Enteric Neural Crest",
                                     "69HPFcl10:Neural progenitor",
                                     "69HPFcl14:Glia",
                                     "48HPFcl15:Glia",
                                     "48HPFcl13:Neural progenitors",
                                     "48HPFcl17:Sensory Neuronal progenitor",
                                     "69HPFcl5:Immature Neuron",
                                     "69HPFcl12:Enteric Neurons",
                                     "48HPFcl0:Autonomic Neuronal Progenitor",
                                     "48HPFcl7:Neural"))

Subset.type <- "NeuralSubclusters"
for (ii in 1:length(plottype)){
  
  #Original Cell IDents
  FP <- DimPlot(nc3.Neural, reduction = plottype[ii], group.by = "seurat_clusters", label = FALSE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"_OriginalCellIdents",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Native Cluster Idents
  FP <- DimPlot(nc3.Neural, reduction = plottype[ii], group.by = "orig.ident", label = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"_NativeClusterIdents",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Cell Type Idents
  FP <- DimPlot(nc3.Neural, reduction = plottype[ii], group.by = "CellType", label = TRUE, repel = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"CellTypeMapping",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Cell Type Idents
  FP <- DimPlot(nc3.Neural, reduction = plottype[ii], group.by = "CellType", label = FALSE, repel = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"CellTypeMapping_NoLabel",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
}

# Pigment Sub out for paper ---------------------------------------------------------
Idents(nc3) <-"CellType"
nc3.pigment <- subset(nc3, idents = c("69HPFcl4:Melanophores",
                                     "48HPFcl8:Melanophore",
                                     "69HPFcl18:Proliferative Melanophores",
                                     "69HPFcl13:Pigment Progenitor",
                                     "69HPFcl16:Iridiophore",
                                     "69HPFcl15:Xanthophore"))

Subset.type <- "PigmentSubclusters"
for (ii in 1:length(plottype)){
  
  #Original Cell IDents
  FP <- DimPlot(nc3.pigment, reduction = plottype[ii], group.by = "seurat_clusters", label = FALSE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"_OriginalCellIdents",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Native Cluster Idents
  FP <- DimPlot(nc3.pigment, reduction = plottype[ii], group.by = "orig.ident", label = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"_NativeClusterIdents",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Cell Type Idents
  FP <- DimPlot(nc3.pigment, reduction = plottype[ii], group.by = "CellType", label = TRUE, repel = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"CellTypeMapping",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Cell Type Idents
  FP <- DimPlot(nc3.pigment, reduction = plottype[ii], group.by = "CellType", label = FALSE, repel = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"CellTypeMapping_NoLabel",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
}
# Mesenchyme Sub out for paper ---------------------------------------------------------
Idents(nc3) <-"CellType"
nc3.mesenchyme <- subset(nc3, idents = c("69HPFcl8:Proliferative Chondrocytic mesenchyme",
                                     "48HPFcl6:Mesenchymal progenitor 2",
                                     "48HPFcl1:Mesenchymal progenitor 1",
                                     "69HPFcl2:Mesenchymal progenitor",
                                     "69HPFcl0:Differentiating chondrocytic mesenchyme",
                                     "48HPFcl11:Mesenchymal progenitor 4",
                                     "48HPFcl12:Mesenchymal progenitor 5",
                                     "48HPFcl10:Mesenchymal progenitor 3",
                                     "69HPFcl6:Differentiating chondrogenic mesenchyme",
                                     "69HPFcl20:Migratory chondrocytic mesenchyme",
                                     "48HPFcl4:Chondrocytic mesenchyme stem-like",
                                     "69HPFcl21:Unknown",
                                     "48HPFcl9:Chondrocytic migratory mesenchyme",
                                     "48HPFcl3:Chondrocytic mesenchyme proliferative 2",
                                     "48HPFcl2:Chondrocytic mesenchyme proliferative 1",
                                     "69HPFcl17:Mesenchymal progenitor",
                                     "69HPFcl22:Mesenchymal progenitor"))
printpath<- "C:/Users/agh8/Desktop/Rice University/Uribe lab/Transcriptome/R programming/Seurat Plots/Merge2/"
saveRDS(nc3, file = paste(printpath,ProjectTitle,"_PC",pcnum,"+res",res,"_MesenchymeSubset.rds",sep = ""))

nc3.mesenchyme <- NormalizeData(nc3.mesenchyme, normalization.method = "LogNormalize", scale.factor = 10000)
nc3.mesenchyme <- FindVariableFeatures(nc3.mesenchyme, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(nc3.mesenchyme)
nc3.mesenchyme <- ScaleData(nc3.mesenchyme, features = all.genes)
nc3.mesenchyme <- RunPCA(nc3.mesenchyme, features = VariableFeatures(object =nc3.mesenchyme))
nc3.mesenchyme <- JackStraw(nc3.mesenchyme, reduction = "pca", num.replicate = 100,dims = pcnum)
nc3.mesenchyme <- ScoreJackStraw(nc3.mesenchyme, dims = 1:pcnum)
nc3.mesenchyme <- FindNeighbors(nc3.mesenchyme, dims = 1:pcnum)
nc3.mesenchyme <- FindClusters(nc3.mesenchyme, resolution = res)
head(Idents(nc3.mesenchyme),5)
nc3.mesenchyme <- RunTSNE(object = nc3.mesenchyme, dims = 1:pcnum, do.fast = TRUE )
nc3.mesenchyme <- RunUMAP(nc3.mesenchyme,dims = 1:pcnum)

Subset.type <- "MesenchyeSubclusters"
for (ii in 1:length(plottype)){
  
  #Original Cell IDents
  FP <- DimPlot(nc3.mesenchyme, reduction = plottype[ii], group.by = "seurat_clusters", label = FALSE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"_OriginalCellIdents",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Native Cluster Idents
  FP <- DimPlot(nc3.mesenchyme, reduction = plottype[ii], group.by = "orig.ident", label = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"_NativeClusterIdents",".tiff",sep = "")
  tiff(filename=fname,width = imagesize, height = imagesize,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Cell Type Idents
  FP <- DimPlot(nc3.mesenchyme, reduction = plottype[ii], group.by = "CellType", label = TRUE, repel = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"CellTypeMapping",".tiff",sep = "")
  tiff(filename=fname,width = imagesize+1000, height = imagesize+1000,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
  
  #Cell Type Idents
  FP <- DimPlot(nc3.mesenchyme, reduction = plottype[ii], group.by = "CellType", label = FALSE, repel = TRUE)
  fname = paste(printpath,ProjectTitle,"_",Subset.type,"_",plottype[ii],"_PC",pcnum,"+res",res,"CellTypeMapping_NoLabel",".tiff",sep = "")
  tiff(filename=fname,width = imagesize+1000, height = imagesize+1000,units = "px",res = 300,pointsize = 12)
  print(FP)
  dev.off()
}
  
  #find all unique cluster biomarkers
  nc3.mesenchyme <- FindAllMarkers(nc3.mesenchyme, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  nc3.mesenchyme %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  figtitle <- paste(printpath,"Merge2_MesenchymeSubcluster_Genes",sep="")
  write.csv(nc3.mesenchyme,figtitle)
# Pigment Sub out ---------------------------------------------------------
  
  Idents(nc3)<-"CellType"
  nc3.pigment <- subset(nc3, idents = c("69HPFcl4:Melanophores","69HPFcl13: Pigment Progenitor","69HPFcl15:Xanthophores",
                                        "69HPFcl16:Iridiophore","69HPFcl18:Melanophores","48HPFcl8:Melanophores", 
                                        "48HPFcl5:Neural Crest",'48HPFcl5:Neural Crest'))
  
  
  
  #scaling data
  nc3.pigment <- ScaleData(nc3.pigment)
  
  #linear dimensional reduction (i.e. find the PCs)
  nc3.pigment <- RunPCA(nc3.pigment, features = VariableFeatures(object =nc3.pigment), npcs = pcnum)
  nc3.pigment <- RunUMAP(nc3.pigment, reduction = "pca", dims = 1:pcnum)
  nc3.pigment <- RunTSNE(nc3.pigment, reduction = "pca", dims = 1:pcnum)
  nc3.pigment <- FindNeighbors(nc3.pigment, reduction = "pca", dims = 1:pcnum)
  nc3.pigment <- FindClusters(nc3.pigment, resolution = res)
  
  for (ii in 1:length(plottype)){
    
    #Original Cell IDents
    FP <- DimPlot(nc3.pigment, reduction = plottype[ii], group.by = "orig.ident", label = FALSE)
    fname = paste(printpath,ProjectTitle,"_PigmentSubclusters_",plottype[ii],"_PC",pcnum,"+res",res,"_OriginalCellIdents",".tiff",sep = "")
    tiff(filename=fname,width = imagesize-2000, height = imagesize-2000,units = "px",res = 300,pointsize = 12)
    print(FP)
    dev.off()
    
    #Native Cluster Idents
    FP <- DimPlot(nc3.pigment, reduction = plottype[ii], label = TRUE)
    fname = paste(printpath,ProjectTitle,"_PigmentSubclusters_",plottype[ii],"_PC",pcnum,"+res",res,"_NativeClusterIdents",".tiff",sep = "")
    tiff(filename=fname,width = imagesize-2000, height = imagesize-2000,units = "px",res = 300,pointsize = 12)
    print(FP)
    dev.off()
    
    #Cell Type Idents
    FP <- DimPlot(nc3.pigment, reduction = plottype[ii], group.by = "CellType", label = TRUE, repel = TRUE) + NoLegend()
    fname = paste(printpath,ProjectTitle,"_PigmentSubclusters_",plottype[ii],"_PC",pcnum,"+res",res,"CellTypeMapping",".tiff",sep = "")
    tiff(filename=fname,width = imagesize+1000, height = imagesize+1000,units = "px",res = 300,pointsize = 12)
    print(FP)
    dev.off()
  }
    

# ENS Subcluster formation ------------------------------------------------

  
  filepath<- "C:/Users/agh8/Desktop/Rice University/Uribe lab/Transcriptome/R programming/Seurat Plots/ENSSubcluster_Analysis/"
  printpath_ENS <- paste(filepath,"ENSSubcluster_Analysis/",sep = "")
  ProjectTitle.ens<- "69hpf subcluster 5,12"
  
  nc1.ENS.sub <- readRDS(file = paste(filepath,ProjectTitle.ens,".rds", sep = ""))
  nc1.ENS.sub<-RunUMAP(object = nc1.ENS.sub, dims = 1:6 )
  
  
  
  nc1.ENS.sub <- BuildClusterTree(nc1.ENS.sub, verbose = FALSE, reorder = TRUE)
  data.tree <- Tool(object = nc1.ENS.sub, slot = "BuildClusterTree")
  fname = paste(printpath,ProjectTitle,"_","PC",pcnum,"+res",res,"_ClustTree_AllHox",".tiff",sep = "")
  tiff(filename=fname,width = 1900, height = 2100,units = "px",res = 300,pointsize = 12)
  plot.phylo(x = data.tree,type = "phy", 
             show.node.label = TRUE, 
             no.margin = TRUE,node.pos = 2,
             edge.width = 2, label.offset = 0.2 )
  dev.off()
  
  saveRDS(object = nc1.ENS.sub, file = paste(filepath,ProjectTitle.ens,"_V2.rds", sep = ""))
  