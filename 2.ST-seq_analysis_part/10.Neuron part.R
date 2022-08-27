# 10.1 Loading package

library(Seurat)
library(ggplot2)
library(ggsci)
library(RCurl)
library(viridis)

# 10.2 loading dataset and subset Neuron

data <- readRDS("CM_ST_dataset.rds") 
Neuron <- subset(data,subset = SingleR.labels == "Neuron")
Neuron <- SCTransform(Neuron,return.only.var.genes = FALSE,assay = "Spatial",verbose = FALSE)
Neuron <- RunPCA(object = Neuron,verbose = FALSE)
Neuron <- FindNeighbors(Neuron,dim=1:20)
Neuron <- FindClusters(Neuron,resolution = 0.8)
Neuron <- RunUMAP(Neuron,reduction="pca", dims = 1:20)
Neuron <- RunTSNE(Neuron,reduction="pca", dims = 1:20)

DimPlot(Neuron,reduction = "umap",pt.size = 1,label = T,label.size = 5)
DimPlot(Neuron,reduction = "tsne",pt.size = 1,label = T,label.size = 5)

new.cluster.ids <- c("HY","STR","STR","CTX-2","CTX-1","HB","CTX-1","CB",
                       "STR","OB","TH","CB","CTX-2","HY","MB","OB","STR",
                       "CTX-1","HY","HC","HB","HC","CTX-2","STR","CB","CB","HC","HY")
names(new.cluster.ids) <- levels(Neuron)
Neuron <- RenameIdents(Neuron,new.cluster.ids)
saveRDS(Neuron,"Neuron_anno.rds")
  
table(Neuron2@active.ident)
table(Neuron2$region_type)

# 10.3 Doheatmap visulization
                                                         
DoHeatmap(Neuron2,group.by = "region_type",
          c( "Dcx","Slc6a11",#"Omp","Shisa3","Cdhr1",# OB
              "Nrgn","Stx1a","Fezf2","Camk2n1",# CTX
              "Hpca","Itpka","Dsp","Lct","Spink8",#"Camkv",# HC
             "Drd2","Tal1",#"Pak1", # MB
             "Pcp2","Cbln1","Cbln3","Pvalb","Calb1",# CB
             "Slc18a3","Slc10a4","Calca", # HB
              "Prkcd","Rora","Rgs16","Tnnt1","Pcp4",# TH
             "Zcchc12","Ahi1",#"Trh",# HY
              "Plp1","Mbp"#,"Gpr6","Rgs9"# STR
              ))

DoHeatmap(Neuron,group.by = "integrated_snn_res.0.8",
          c( "Dcx","Slc6a11","Omp","Shisa3","Cdhr1" # OB
            "Nrgn","Stx1a","Fezf2","Camk2n1" # CTX
            "Hpca","Itpka","Dsp","Lct","Spink8","Camkv" # HC
            "Prkcd","Rora","Rgs16","Tnnt1","Pcp4" # TH
            "Plp1","Mbp","Gpr6","Drd2","Rgs9" # STR
            "Zcchc12","Ahi1","Trh" # HY
            "Pcp2","Cbln1","Cbln3","Pvalb","Calb1" # CB
            "Slc18a3","Slc10a4","Hoxb5","Calca" # HB
          ))

Neuron2$region_type <- Neuron2@active.ident
table(Neuron2$region_type)
Neuron2$region_type <- factor(Neuron2$region_type,levels = c("OB","CTX-1","CTX-2","HC","MB","CB","HB","TH","HY","STR"))

# 10.4 SpatialDimPlot

SpatialDimPlot(Neuron,crop = F,pt.size = 1.5,images = c("anterior1","posterior1"))

# 10.5 SpatialFeaturePlot

# 10.5.1 Olfactory bulb（OB）
SpatialFeaturePlot(Neuron,c("Dcx"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Slc6a11"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Omp"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Shisa3"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Cdhr1"),pt.size = 1.5,images = c( "anterior1","posterior1"))

# 10.5.2 Cortex (CTX-1,CTX-2)
SpatialFeaturePlot(Neuron,c("Nrgn"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Stx1a"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Fezf2"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Camk2n1"),pt.size = 1.5,images = c( "anterior1","posterior1"))

# 10.5.3 Hippocampus (HC)
SpatialFeaturePlot(Neuron,c("Hpca"),pt.size = 1.5,images = c( "anterior1","posterior1"))/ 
  SpatialFeaturePlot(Neuron,c("Itpka"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Dsp"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Lct"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Spink8"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Camkv"),pt.size = 1.5,images = c( "anterior1","posterior1"))

# 10.5.4 Thalamus（TH）
SpatialFeaturePlot(Neuron,c("Prkcd"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Rora"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Rgs16"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Tnnt1"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Pcp4"),pt.size = 1.5,images = c( "anterior1","posterior1"))

# 10.5.5 Striatum（STR）
SpatialFeaturePlot(Neuron,c("Plp1"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Mbp"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Gpr6"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Drd2"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Rgs9"),pt.size = 1.5,images = c( "anterior1","posterior1"))

# 10.5.6 Hypothalamus（HY）
SpatialFeaturePlot(Neuron,c("Zcchc12"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Ahi1"),pt.size = 1.5,images = c("anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Trh"),pt.size = 1.5,images = c("anterior1","posterior1"))

# 10.5.7 Cerebellum (CB)
SpatialFeaturePlot(Neuron,c("Pcp2"),pt.size = 1.5,images = c( "anterior1","posterior1")) /
  SpatialFeaturePlot(Neuron,c("Cbln1"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Cbln3"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Pvalb"),pt.size = 1.5,images = c( "anterior1","posterior1"))/
  SpatialFeaturePlot(Neuron,c("Calb1"),pt.size = 1.5,images = c( "anterior1","posterior1"))

# 10.5.8 Pons & medulla (MB)
SpatialFeaturePlot(Neuron,c("Slc18a3"),pt.size = 1.5,images = c("BM1_ST_P","BM2_ST_P","BA1_ST_P","BA2_ST_P")) /
  SpatialFeaturePlot(Neuron,c("Slc10a4"),pt.size = 1.5,images = c("BM1_ST_P","BM2_ST_P","BA1_ST_P","BA2_ST_P")) /
  SpatialFeaturePlot(Neuron,c("Hoxb5"),pt.size = 1.5,images = c("BM1_ST_P","BM2_ST_P","BA1_ST_P","BA2_ST_P")) /
  SpatialFeaturePlot(Neuron,c("Calca"),pt.size = 1.5,images = c("BM1_ST_P","BM2_ST_P","BA1_ST_P","BA2_ST_P")) 

# 10.6 DEGs and GO enrichment 

# 10.6.1 Extract different region of Neuron dataset

OB <- subset(Neuron,idents = "OB")  
CTX_1 <- subset(Neuron,idents = "CTX-1")
CTX_2 <- subset(Neuron,idents = "CTX-2")
HC <- subset(Neuron,idents = "HC")
MB <- subset(Neuron,idents = "MB")
CB <- subset(Neuron,idents = "CB") 
HB <- subset(Neuron,idents = "HB") 
TH <- subset(Neuron,idents = "TH")
HY <- subset(Neuron,idents = "HY") 
STR <- subset(Neuron,idents = "STR")

OB@active.ident <- OB$type
CTX_1@active.ident <- CTX_1$type
CTX_2@active.ident <- CTX_2$type
HC@active.ident <- HC$type
MB@active.ident <- MB$type
CB@active.ident <- CB$type
HB@active.ident <- HB$type
TH@active.ident <- TH$type
HY@active.ident <- HY$type
STR@active.ident <- STR$type

# 10.6.2 DEGs analysis of Model versus control (MVC)

OB_DEG_MVC <- FindMarkers(OB,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.threshold = 0.25) 
OB_DEG_MVC$change = ifelse(OB_DEG_MVC$p_val_adj < 0.05 & abs(OB_DEG_MVC$avg_log2FC) >= 0.25, 
                           ifelse(OB_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(OB_DEG_MVC$change)
# Down Stable     Up 
# 289   1467    424 

CTX_1_DEG_MVC <- FindMarkers(CTX_1,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.threshold = 0.25) 
CTX_1_DEG_MVC$change = ifelse(CTX_1_DEG_MVC$p_val_adj < 0.05 & abs(CTX_1_DEG_MVC$avg_log2FC) >= 0.25, 
                              ifelse(CTX_1_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CTX_1_DEG_MVC$change)
# Down Stable     Up 
# 93   1190    191

CTX_2_DEG_MVC <- FindMarkers(CTX_2,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.threshold = 0.25) 
CTX_2_DEG_MVC$change = ifelse(CTX_2_DEG_MVC$p_val_adj < 0.05 & abs(CTX_2_DEG_MVC$avg_log2FC) >= 0.25, 
                              ifelse(CTX_2_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CTX_2_DEG_MVC$change)
# Down Stable     Up 
# 72    820    157

HC_DEG_MVC <- FindMarkers(HC,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.threshold = 0.25) 
HC_DEG_MVC$change = ifelse(HC_DEG_MVC$p_val_adj < 0.05 & abs(HC_DEG_MVC$avg_log2FC) >= 0.25, 
                           ifelse(HC_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(HC_DEG_MVC$change)
# Down Stable     Up 
# 74   1353     96

MB_DEG_MVC <- FindMarkers(MB,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.threshold = 0.25) 
MB_DEG_MVC$change = ifelse(MB_DEG_MVC$p_val_adj < 0.05 & abs(MB_DEG_MVC$avg_log2FC) >= 0.25, 
                           ifelse(MB_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(MB_DEG_MVC$change)
# Down Stable     Up 
# 50   1322    180

CB_DEG_MVC <- FindMarkers(CB,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.threshold = 0.25) 
CB_DEG_MVC$change = ifelse(CB_DEG_MVC$p_val_adj < 0.05 & abs(CB_DEG_MVC$avg_log2FC) >= 0.25, 
                           ifelse(CB_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CB_DEG_MVC$change)
# Down Stable     Up 
# 65    906    186

HB_DEG_MVC <- FindMarkers(HB,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.threshold = 0.25) 
HB_DEG_MVC$change = ifelse(HB_DEG_MVC$p_val_adj < 0.05 & abs(HB_DEG_MVC$avg_log2FC) >= 0.25, 
                           ifelse(HB_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(HB_DEG_MVC$change)
# Down Stable     Up 
# 38    994    106

TH_DEG_MVC <- FindMarkers(TH,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.threshold = 0.25) 
TH_DEG_MVC$change = ifelse(TH_DEG_MVC$p_val_adj < 0.05 & abs(TH_DEG_MVC$avg_log2FC) >= 0.25, 
                           ifelse(TH_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(TH_DEG_MVC$change)
# Down Stable     Up 
# 36   1002    140 

HY_DEG_MVC <- FindMarkers(HY,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.HYreshold = 0.25) 
HY_DEG_MVC$change = ifelse(HY_DEG_MVC$p_val_adj < 0.05 & abs(HY_DEG_MVC$avg_log2FC) >= 0.25, 
                           ifelse(HY_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(HY_DEG_MVC$change)
# Down Stable     Up 
# 197   1084    504 

STR_DEG_MVC <- FindMarkers(STR,ident.1 = "Model",ident.2 = "Control",min.pct = 0.25, logfc.STRreshold = 0.25) 
STR_DEG_MVC$change = ifelse(STR_DEG_MVC$p_val_adj < 0.05 & abs(STR_DEG_MVC$avg_log2FC) >= 0.25, 
                            ifelse(STR_DEG_MVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(STR_DEG_MVC$change)
# Down Stable     Up 
# 79   1670    249 

# 10.6.3 Upsetplot of MVC DEGs visualization

OB_DEG_MVC <- row.names(OB_DEG_MVC[OB_DEG_MVC$change != "Stable",]) # n = 722
CTX_1_DEG_MVC <- row.names(CTX_1_DEG_MVC[CTX_1_DEG_MVC$change != "Stable",]) # n = 377
CTX_2_DEG_MVC <- row.names(CTX_2_DEG_MVC[CTX_2_DEG_MVC$change != "Stable",]) # n = 424
HC_DEG_MVC <- row.names(HC_DEG_MVC[HC_DEG_MVC$change != "Stable",])# n = 253
MB_DEG_MVC <- row.names(MB_DEG_MVC[MB_DEG_MVC$change != "Stable",])# n = 301
CB_DEG_MVC <- row.names(CB_DEG_MVC[CB_DEG_MVC$change != "Stable",])# n = 488
HB_DEG_MVC <- row.names(HB_DEG_MVC[HB_DEG_MVC$change != "Stable",])# n = 236
TH_DEG_MVC <- row.names(TH_DEG_MVC[TH_DEG_MVC$change != "Stable",])# n = 322
HY_DEG_MVC <- row.names(HY_DEG_MVC[HY_DEG_MVC$change != "Stable",])# n = 820
STR_DEG_MVC <- row.names(STR_DEG_MVC[STR_DEG_MVC$change != "Stable",])# n = 334

ALL_DEG <- union(OB_DEG_MVC,c(CTX_1_DEG_MVC,CTX_2_DEG_MVC,HC_DEG_MVC,MB_DEG_MVC,
                              CB_DEG_MVC,HB_DEG_MVC,TH_DEG_MVC,HY_DEG_MVC,STR_DEG_MVC))
ALL_DEG <- as.data.frame(ALL_DEG)
row.names(ALL_DEG) <- ALL_DEG$ALL_DEG
colnames(ALL_DEG) <- c("ID")
head(ALL_DEG)

OB_DEG_MVC <-  as.data.frame(OB_DEG_MVC)
CTX_1_DEG_MVC <- as.data.frame(CTX_1_DEG_MVC)
CTX_2_DEG_MVC <- as.data.frame(CTX_2_DEG_MVC)
HC_DEG_MVC <- as.data.frame(HC_DEG_MVC)
MB_DEG_MVC <- as.data.frame(MB_DEG_MVC)
CB_DEG_MVC <- as.data.frame(CB_DEG_MVC)
HB_DEG_MVC <- as.data.frame(HB_DEG_MVC)
TH_DEG_MVC <- as.data.frame(TH_DEG_MVC)
HY_DEG_MVC <- as.data.frame(HY_DEG_MVC)
STR_DEG_MVC <- as.data.frame(STR_DEG_MVC)

rownames(OB_DEG_MVC) <- OB_DEG_MVC$OB_DEG_MVC
rownames(CTX_1_DEG_MVC) <- CTX_1_DEG_MVC$CTX_1_DEG_MVC
rownames(CTX_2_DEG_MVC) <- CTX_2_DEG_MVC$CTX_2_DEG_MVC
rownames(HC_DEG_MVC) <- HC_DEG_MVC$HC_DEG_MVC
rownames(MB_DEG_MVC) <- MB_DEG_MVC$MB_DEG_MVC
rownames(CB_DEG_MVC) <- CB_DEG_MVC$CB_DEG_MVC
rownames(HB_DEG_MVC) <- HB_DEG_MVC$HB_DEG_MVC
rownames(TH_DEG_MVC) <- TH_DEG_MVC$TH_DEG_MVC
rownames(HY_DEG_MVC) <- HY_DEG_MVC$HY_DEG_MVC
rownames(STR_DEG_MVC) <- STR_DEG_MVC$STR_DEG_MVC

OB_DEG_MVC$ID <- OB_DEG_MVC$OB_DEG_MVC
CTX_1_DEG_MVC$ID <- CTX_1_DEG_MVC$CTX_1_DEG_MVC
CTX_2_DEG_MVC$ID <- CTX_2_DEG_MVC$CTX_2_DEG_MVC
HC_DEG_MVC$ID <- HC_DEG_MVC$HC_DEG_MVC
MB_DEG_MVC$ID <- MB_DEG_MVC$MB_DEG_MVC
CB_DEG_MVC$ID <- CB_DEG_MVC$CB_DEG_MVC
HB_DEG_MVC$ID <- HB_DEG_MVC$HB_DEG_MVC
TH_DEG_MVC$ID <- TH_DEG_MVC$TH_DEG_MVC
HY_DEG_MVC$ID <- HY_DEG_MVC$HY_DEG_MVC
STR_DEG_MVC$ID <- STR_DEG_MVC$STR_DEG_MVC

merge_DEG <- left_join(ALL_DEG,OB_DEG_MVC,by="ID") 
merge_DEG <- left_join(merge_DEG,CTX_1_DEG_MVC,by="ID")
merge_DEG <- left_join(merge_DEG,CTX_2_DEG_MVC,by="ID")
merge_DEG <- left_join(merge_DEG,HC_DEG_MVC,by="ID")
merge_DEG <- left_join(merge_DEG,MB_DEG_MVC,by="ID")
merge_DEG <- left_join(merge_DEG,CB_DEG_MVC,by="ID")
merge_DEG <- left_join(merge_DEG,HB_DEG_MVC,by="ID")
merge_DEG <- left_join(merge_DEG,TH_DEG_MVC,by="ID")
merge_DEG <- left_join(merge_DEG,HY_DEG_MVC,by="ID")
merge_DEG <- left_join(merge_DEG,STR_DEG_MVC,by="ID")

head(merge_DEG)
write.csv(merge_DEG,"Neuron_merge_DEG.csv")

merge_DEG_t <- merge_DEG
rownames(merge_DEG_t) <- merge_DEG_t$ID

head(merge_DEG_t)
merge_DEG_t[which(!is.na(merge_DEG_t),arr.ind = T)]<-1
merge_DEG_t[which(is.na(merge_DEG_t),arr.ind = T)]<-0
head(merge_DEG_t)
merge_DEG_t <- as.data.frame(lapply(merge_DEG_t,as.numeric))
merge_DEG_t <- merge_DEG_t[,-1]
str(merge_DEG_t)

colnames(merge_DEG_t) <- c("OB","CTX-1","CTX-2","HC","MB",
                           "CB","HB","TH","HY","STR")

upset(merge_DEG_t, nsets = 10,nintersects = 25,
      mb.ratio = c(0.6, 0.4),
      order.by = c("freq"),
      decreasing = c(TRUE,T), 
      sets.bar.color = cols,
      mainbar.y.label = "Intersection size of DEGs(Up)", sets.x.label = "DEG number",
      text.scale = 1.4) # 6*5

# 10.6.4 Subset MVC upregulated DEGS

OB_DEG_MVC_up <- row.names(OB_DEG_MVC[OB_DEG_MVC$change == "Up",]) # n = 430
CTX_1_DEG_MVC_up <- row.names(CTX_1_DEG_MVC[CTX_1_DEG_MVC$change == "Up",]) # n = 259
CTX_2_DEG_MVC_up <- row.names(CTX_2_DEG_MVC[CTX_2_DEG_MVC$change == "Up",]) # n = 289
HC_DEG_MVC_up <- row.names(HC_DEG_MVC[HC_DEG_MVC$change == "Up",])# n = 142
MB_DEG_MVC_up <- row.names(MB_DEG_MVC[MB_DEG_MVC$change == "Up",])# n = 232
CB_DEG_MVC_up <- row.names(CB_DEG_MVC[CB_DEG_MVC$change == "Up",])# n = 372
HB_DEG_MVC_up <- row.names(HB_DEG_MVC[HB_DEG_MVC$change == "Up",])# n = 134
TH_DEG_MVC_up <- row.names(TH_DEG_MVC[TH_DEG_MVC$change == "Up",])# n = 256
HY_DEG_MVC_up <- row.names(HY_DEG_MVC[HY_DEG_MVC$change == "Up",])# n = 626
STR_DEG_MVC_up <- row.names(STR_DEG_MVC[STR_DEG_MVC$change == "Up",])# n = 253

# 10.6.5 GO enrichment of DEGs (MVC_up)

OB_DEG_MVC_up_GO <-  enrichGO(gene = OB_DEG_MVC_up$OB_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CTX_1_DEG_MVC_up_GO <- enrichGO(gene = CTX_1_DEG_MVC_up$CTX_1_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CTX_2_DEG_MVC_up_GO <- enrichGO(gene = CTX_2_DEG_MVC_up$CTX_2_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
HC_DEG_MVC_up_GO <- enrichGO(gene = HC_DEG_MVC_up$HC_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
MB_DEG_MVC_up_GO <- enrichGO(gene = MB_DEG_MVC_up$MB_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CB_DEG_MVC_up_GO <- enrichGO(gene = CB_DEG_MVC_up$CB_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
HB_DEG_MVC_up_GO <- enrichGO(gene = HB_DEG_MVC_up$HB_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
TH_DEG_MVC_up_GO <- enrichGO(gene = TH_DEG_MVC_up$TH_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
HY_DEG_MVC_up_GO <- enrichGO(gene = HY_DEG_MVC_up$HY_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
STR_DEG_MVC_up_GO <- enrichGO(gene = STR_DEG_MVC_up$STR_DEG_MVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

(dotplot(OB_DEG_MVC_up_GO,title="OB region") + dotplot(CTX_1_DEG_MVC_up_GO,title="CTX_1 region")) /
  (dotplot(CTX_2_DEG_MVC_up_GO,title="CTX_2 region") + dotplot(HC_DEG_MVC_up_GO,title="HC region")) /
  (dotplot(MB_DEG_MVC_up_GO,title="MB region") + dotplot(CB_DEG_MVC_up_GO,title="CB region")) /
  (dotplot(HB_DEG_MVC_up_GO,title="HB region") + dotplot(TH_DEG_MVC_up_GO,title="TH region")) /
  (dotplot(HY_DEG_MVC_up_GO,title="HY region") + dotplot(STR_DEG_MVC_up_GO,title="STR region"))

OB_DEG_MVC_up_GO
OB_DEG_MVC_up_GO_3 <- OB_DEG_MVC_up_GO@result[c(1,2,4),c("Description","pvalue","Count")]
OB_DEG_MVC_up_GO_3$log10Pvalue <- -log10(OB_DEG_MVC_up_GO_3$pvalue)
OB_DEG_MVC_up_GO_3$subtype <- "OB"
OB_DEG_MVC_up_GO_3

CTX_1_DEG_MVC_up_GO
CTX_1_DEG_MVC_up_GO_3 <- CTX_1_DEG_MVC_up_GO@result[c(3:5),c("Description","pvalue","Count")]
CTX_1_DEG_MVC_up_GO_3$log10Pvalue <- -log10(CTX_1_DEG_MVC_up_GO_3$pvalue)
CTX_1_DEG_MVC_up_GO_3$subtype <- "CTX_1"
CTX_1_DEG_MVC_up_GO_3

CTX_2_DEG_MVC_up_GO
CTX_2_DEG_MVC_up_GO_3 <- CTX_2_DEG_MVC_up_GO@result[c(1,2,5),c("Description","pvalue","Count")]
CTX_2_DEG_MVC_up_GO_3$log10Pvalue <- -log10(CTX_2_DEG_MVC_up_GO_3$pvalue)
CTX_2_DEG_MVC_up_GO_3$subtype <- "CTX_2"
CTX_2_DEG_MVC_up_GO_3

HC_DEG_MVC_up_GO
HC_DEG_MVC_up_GO_3 <- HC_DEG_MVC_up_GO@result[c(6,8,10),c("Description","pvalue","Count")]
HC_DEG_MVC_up_GO_3$log10Pvalue <- -log10(HC_DEG_MVC_up_GO_3$pvalue)
HC_DEG_MVC_up_GO_3$subtype <- "HC"
HC_DEG_MVC_up_GO_3

MB_DEG_MVC_up_GO
MB_DEG_MVC_up_GO_3 <- MB_DEG_MVC_up_GO@result[c(2,8,10),c("Description","pvalue","Count")]
MB_DEG_MVC_up_GO_3$log10Pvalue <- -log10(MB_DEG_MVC_up_GO_3$pvalue)
MB_DEG_MVC_up_GO_3$subtype <- "MB"
MB_DEG_MVC_up_GO_3

CB_DEG_MVC_up_GO
CB_DEG_MVC_up_GO_3 <- CB_DEG_MVC_up_GO@result[c(1,6,10),c("Description","pvalue","Count")]
CB_DEG_MVC_up_GO_3$log10Pvalue <- -log10(CB_DEG_MVC_up_GO_3$pvalue)
CB_DEG_MVC_up_GO_3$subtype <- "CB"
CB_DEG_MVC_up_GO_3

HB_DEG_MVC_up_GO
HB_DEG_MVC_up_GO_3 <- HB_DEG_MVC_up_GO@result[c(3,6,9),c("Description","pvalue","Count")]
HB_DEG_MVC_up_GO_3$log10Pvalue <- -log10(HB_DEG_MVC_up_GO_3$pvalue)
HB_DEG_MVC_up_GO_3$subtype <- "HB"
HB_DEG_MVC_up_GO_3

TH_DEG_MVC_up_GO
TH_DEG_MVC_up_GO_3 <- TH_DEG_MVC_up_GO@result[c(1,4,3,9),c("Description","pvalue","Count")]
TH_DEG_MVC_up_GO_3$log10Pvalue <- -log10(TH_DEG_MVC_up_GO_3$pvalue)
TH_DEG_MVC_up_GO_3$subtype <- "TH"
TH_DEG_MVC_up_GO_3

HY_DEG_MVC_up_GO
HY_DEG_MVC_up_GO_3 <- HY_DEG_MVC_up_GO@result[c(1:3),c("Description","pvalue","Count")]
HY_DEG_MVC_up_GO_3$log10Pvalue <- -log10(HY_DEG_MVC_up_GO_3$pvalue)
HY_DEG_MVC_up_GO_3$subtype <- "HY"
HY_DEG_MVC_up_GO_3

STR_DEG_MVC_up_GO
STR_DEG_MVC_up_GO_3 <- STR_DEG_MVC_up_GO@result[c(1,7,10),c("Description","pvalue","Count")]
STR_DEG_MVC_up_GO_3$log10Pvalue <- -log10(STR_DEG_MVC_up_GO_3$pvalue)
STR_DEG_MVC_up_GO_3$subtype <- "STR"
STR_DEG_MVC_up_GO_3

Go_DEG_up <- rbind(OB_DEG_MVC_up_GO_3,CTX_1_DEG_MVC_up_GO_3,CTX_2_DEG_MVC_up_GO_3,
                   HC_DEG_MVC_up_GO_3,MB_DEG_MVC_up_GO_3,CB_DEG_MVC_up_GO_3,HB_DEG_MVC_up_GO_3,
                   TH_DEG_MVC_up_GO_3,HY_DEG_MVC_up_GO_3,STR_DEG_MVC_up_GO_3)
Go_DEG_up$item <- row.names(Go_DEG_up)
Go_DEG_up$subtype <- factor(Go_DEG_up$subtype,levels = c("OB","CTX_1","CTX_2","HC","MB",
                                                         "CB","HB","TH","HY","STR"))
Go_DEG_up$pathways <- paste0(Go_DEG_up$item," _ ",Go_DEG_up$Description)
ggbarplot(Go_DEG_up, x="pathways", y="log10Pvalue", fill = "subtype", 
          color = "subtype",
          palette =  cols,
          sort.val = "asc",
          sort.by.grodowns=TRUE, 
          x.text.angle=45, 
          xlab = NULL) 

# 10.6.6 Subset AVC upregulated DEGS

OB_DEG_AVC <- FindMarkers(OB,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.threshold = 0.25) 
OB_DEG_AVC$change = ifelse(OB_DEG_AVC$p_val_adj < 0.05 & abs(OB_DEG_AVC$avg_log2FC) >= 0.25, 
                           ifelse(OB_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(OB_DEG_AVC$change)
# Down Stable     Up 
# 244   1612    518

CTX_1_DEG_AVC <- FindMarkers(CTX_1,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.threshold = 0.25) 
CTX_1_DEG_AVC$change = ifelse(CTX_1_DEG_AVC$p_val_adj < 0.05 & abs(CTX_1_DEG_AVC$avg_log2FC) >= 0.25, 
                              ifelse(CTX_1_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CTX_1_DEG_AVC$change)
# Down Stable     Up 
# 84   1788    192

CTX_2_DEG_AVC <- FindMarkers(CTX_2,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.threshold = 0.25) 
CTX_2_DEG_AVC$change = ifelse(CTX_2_DEG_AVC$p_val_adj < 0.05 & abs(CTX_2_DEG_AVC$avg_log2FC) >= 0.25, 
                              ifelse(CTX_2_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CTX_2_DEG_AVC$change)
# Down Stable     Up 
# 82   1378    184

HC_DEG_AVC <- FindMarkers(HC,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.threshold = 0.25) 
HC_DEG_AVC$change = ifelse(HC_DEG_AVC$p_val_adj < 0.05 & abs(HC_DEG_AVC$avg_log2FC) >= 0.25, 
                           ifelse(HC_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(HC_DEG_AVC$change)
# Down Stable     Up 
# 224   1680    240

MB_DEG_AVC <- FindMarkers(MB,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.threshold = 0.25) 
MB_DEG_AVC$change = ifelse(MB_DEG_AVC$p_val_adj < 0.05 & abs(MB_DEG_AVC$avg_log2FC) >= 0.25, 
                           ifelse(MB_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(MB_DEG_AVC$change)
# Down Stable     Up 
# 27   1849    124

CB_DEG_AVC <- FindMarkers(CB,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.threshold = 0.25) 
CB_DEG_AVC$change = ifelse(CB_DEG_AVC$p_val_adj < 0.05 & abs(CB_DEG_AVC$avg_log2FC) >= 0.25, 
                           ifelse(CB_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(CB_DEG_AVC$change)
# Down Stable     Up 
# 127   1155    422

HB_DEG_AVC <- FindMarkers(HB,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.threshold = 0.25) 
HB_DEG_AVC$change = ifelse(HB_DEG_AVC$p_val_adj < 0.05 & abs(HB_DEG_AVC$avg_log2FC) >= 0.25, 
                           ifelse(HB_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(HB_DEG_AVC$change)
# Down Stable     Up 
# 49   1696    142 

TH_DEG_AVC <- FindMarkers(TH,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.threshold = 0.25) 
TH_DEG_AVC$change = ifelse(TH_DEG_AVC$p_val_adj < 0.05 & abs(TH_DEG_AVC$avg_log2FC) >= 0.25, 
                           ifelse(TH_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(TH_DEG_AVC$change)
# Down Stable     Up 
# 110   1227    303

HY_DEG_AVC <- FindMarkers(HY,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.HYreshold = 0.25) 
HY_DEG_AVC$change = ifelse(HY_DEG_AVC$p_val_adj < 0.05 & abs(HY_DEG_AVC$avg_log2FC) >= 0.25, 
                           ifelse(HY_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(HY_DEG_AVC$change)
# Down Stable     Up 
# 290   1595    522

STR_DEG_AVC <- FindMarkers(STR,ident.1 = "ART",ident.2 = "Control",min.pct = 0.1, logfc.STRreshold = 0.25) 
STR_DEG_AVC$change = ifelse(STR_DEG_AVC$p_val_adj < 0.05 & abs(STR_DEG_AVC$avg_log2FC) >= 0.25, 
                            ifelse(STR_DEG_AVC$avg_log2FC > 0.25 ,'Up','Down'),'Stable')
table(STR_DEG_AVC$change)
# Down Stable     Up 
# 25   1915     49

# 10.6.7 Upsetplot of MVC DEGs visualization

OB_DEG_AVC <- row.names(OB_DEG_AVC[OB_DEG_AVC$change != "Stable",]) # n = 762
CTX_1_DEG_AVC <- row.names(CTX_1_DEG_AVC[CTX_1_DEG_AVC$change != "Stable",]) # n = 276
CTX_2_DEG_AVC <- row.names(CTX_2_DEG_AVC[CTX_2_DEG_AVC$change != "Stable",]) # n = 266
HC_DEG_AVC <- row.names(HC_DEG_AVC[HC_DEG_AVC$change != "Stable",])# n = 464
MB_DEG_AVC <- row.names(MB_DEG_AVC[MB_DEG_AVC$change != "Stable",])# n = 151
CB_DEG_AVC <- row.names(CB_DEG_AVC[CB_DEG_AVC$change != "Stable",])# n = 549
HB_DEG_AVC <- row.names(HB_DEG_AVC[HB_DEG_AVC$change != "Stable",])# n = 191
TH_DEG_AVC <- row.names(TH_DEG_AVC[TH_DEG_AVC$change != "Stable",])# n = 413
HY_DEG_AVC <- row.names(HY_DEG_AVC[HY_DEG_AVC$change != "Stable",])# n = 812
STR_DEG_AVC <- row.names(STR_DEG_AVC[STR_DEG_AVC$change != "Stable",])# n = 74

ALL_DEG <- union(OB_DEG_AVC,c(CTX_1_DEG_AVC,CTX_2_DEG_AVC,HC_DEG_AVC,MB_DEG_AVC,
                              CB_DEG_AVC,HB_DEG_AVC,TH_DEG_AVC,HY_DEG_AVC,STR_DEG_AVC))
ALL_DEG <- as.data.frame(ALL_DEG)
row.names(ALL_DEG) <- ALL_DEG$ALL_DEG
colnames(ALL_DEG) <- c("ID")
head(ALL_DEG)

OB_DEG_AVC <-  as.data.frame(OB_DEG_AVC)
CTX_1_DEG_AVC <- as.data.frame(CTX_1_DEG_AVC)
CTX_2_DEG_AVC <- as.data.frame(CTX_2_DEG_AVC)
HC_DEG_AVC <- as.data.frame(HC_DEG_AVC)
MB_DEG_AVC <- as.data.frame(MB_DEG_AVC)
CB_DEG_AVC <- as.data.frame(CB_DEG_AVC)
HB_DEG_AVC <- as.data.frame(HB_DEG_AVC)
TH_DEG_AVC <- as.data.frame(TH_DEG_AVC)
HY_DEG_AVC <- as.data.frame(HY_DEG_AVC)
STR_DEG_AVC <- as.data.frame(STR_DEG_AVC)

rownames(OB_DEG_AVC) <- OB_DEG_AVC$OB_DEG_AVC
rownames(CTX_1_DEG_AVC) <- CTX_1_DEG_AVC$CTX_1_DEG_AVC
rownames(CTX_2_DEG_AVC) <- CTX_2_DEG_AVC$CTX_2_DEG_AVC
rownames(HC_DEG_AVC) <- HC_DEG_AVC$HC_DEG_AVC
rownames(MB_DEG_AVC) <- MB_DEG_AVC$MB_DEG_AVC
rownames(CB_DEG_AVC) <- CB_DEG_AVC$CB_DEG_AVC
rownames(HB_DEG_AVC) <- HB_DEG_AVC$HB_DEG_AVC
rownames(TH_DEG_AVC) <- TH_DEG_AVC$TH_DEG_AVC
rownames(HY_DEG_AVC) <- HY_DEG_AVC$HY_DEG_AVC
rownames(STR_DEG_AVC) <- STR_DEG_AVC$STR_DEG_AVC

OB_DEG_AVC$ID <- OB_DEG_AVC$OB_DEG_AVC
CTX_1_DEG_AVC$ID <- CTX_1_DEG_AVC$CTX_1_DEG_AVC
CTX_2_DEG_AVC$ID <- CTX_2_DEG_AVC$CTX_2_DEG_AVC
HC_DEG_AVC$ID <- HC_DEG_AVC$HC_DEG_AVC
MB_DEG_AVC$ID <- MB_DEG_AVC$MB_DEG_AVC
CB_DEG_AVC$ID <- CB_DEG_AVC$CB_DEG_AVC
HB_DEG_AVC$ID <- HB_DEG_AVC$HB_DEG_AVC
TH_DEG_AVC$ID <- TH_DEG_AVC$TH_DEG_AVC
HY_DEG_AVC$ID <- HY_DEG_AVC$HY_DEG_AVC
STR_DEG_AVC$ID <- STR_DEG_AVC$STR_DEG_AVC

merge_DEG <- left_join(ALL_DEG,OB_DEG_AVC,by="ID") 
merge_DEG <- left_join(merge_DEG,CTX_1_DEG_AVC,by="ID")
merge_DEG <- left_join(merge_DEG,CTX_2_DEG_AVC,by="ID")
merge_DEG <- left_join(merge_DEG,HC_DEG_AVC,by="ID")
merge_DEG <- left_join(merge_DEG,MB_DEG_AVC,by="ID")
merge_DEG <- left_join(merge_DEG,CB_DEG_AVC,by="ID")
merge_DEG <- left_join(merge_DEG,HB_DEG_AVC,by="ID")
merge_DEG <- left_join(merge_DEG,TH_DEG_AVC,by="ID")
merge_DEG <- left_join(merge_DEG,HY_DEG_AVC,by="ID")
merge_DEG <- left_join(merge_DEG,STR_DEG_AVC,by="ID")

head(merge_DEG)
write.csv(merge_DEG,"Neuron_merge_DEG_up_AVC.csv")

merge_DEG_t <- merge_DEG
rownames(merge_DEG_t) <- merge_DEG_t$ID

head(merge_DEG_t)
merge_DEG_t[which(!is.na(merge_DEG_t),arr.ind = T)]<-1
merge_DEG_t[which(is.na(merge_DEG_t),arr.ind = T)]<-0
head(merge_DEG_t)
merge_DEG_t <- as.data.frame(lapply(merge_DEG_t,as.numeric))
merge_DEG_t <- merge_DEG_t[,-1]
str(merge_DEG_t)

colnames(merge_DEG_t) <- c("OB","CTX-1","CTX-2","HC","MB",
                           "CB","HB","TH","HY","STR")
cols = c("#C03830","#28B0E6D6","#C8DBC8","#5AAA46","#F8C2D2","#F1B13F","#72BEB7","#FF4848","#B989B9","#3A68AE")

upset(merge_DEG_t, nsets = 10,nintersects = 25,
      mb.ratio = c(0.6, 0.4),
      order.by = c("freq"),
      decreasing = c(TRUE,T), 
      sets.bar.color = cols,
      mainbar.y.label = "Intersection size of DEGs(Up)", sets.x.label = "DEG number",
      text.scale = 1.4) # 6*5

OB_DEG_AVC_up <- row.names(OB_DEG_AVC[OB_DEG_AVC$change == "Up",]) # n = 518
CTX_1_DEG_AVC_up <- row.names(CTX_1_DEG_AVC[CTX_1_DEG_AVC$change == "Up",]) # n = 192
CTX_2_DEG_AVC_up <- row.names(CTX_2_DEG_AVC[CTX_2_DEG_AVC$change == "Up",]) # n = 184
HC_DEG_AVC_up <- row.names(HC_DEG_AVC[HC_DEG_AVC$change == "Up",])# n = 240
MB_DEG_AVC_up <- row.names(MB_DEG_AVC[MB_DEG_AVC$change == "Up",])# n = 124
CB_DEG_AVC_up <- row.names(CB_DEG_AVC[CB_DEG_AVC$change == "Up",])# n = 422
HB_DEG_AVC_up <- row.names(HB_DEG_AVC[HB_DEG_AVC$change == "Up",])# n = 142
TH_DEG_AVC_up <- row.names(TH_DEG_AVC[TH_DEG_AVC$change == "Up",])# n = 303
HY_DEG_AVC_up <- row.names(HY_DEG_AVC[HY_DEG_AVC$change == "Up",])# n = 522
STR_DEG_AVC_up <- row.names(STR_DEG_AVC[STR_DEG_AVC$change == "Up",])# n = 49

# 10.6.8 GO enrichment of DEG_AVC_up
OB_DEG_AVC_up_GO <-  enrichGO(gene = OB_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CTX_1_DEG_AVC_up_GO <- enrichGO(gene = CTX_1_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CTX_2_DEG_AVC_up_GO <- enrichGO(gene = CTX_2_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
HC_DEG_AVC_up_GO <- enrichGO(gene = HC_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
MB_DEG_AVC_up_GO <- enrichGO(gene = MB_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
CB_DEG_AVC_up_GO <- enrichGO(gene = CB_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
HB_DEG_AVC_up_GO <- enrichGO(gene = HB_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
TH_DEG_AVC_up_GO <- enrichGO(gene = TH_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
HY_DEG_AVC_up_GO <- enrichGO(gene = HY_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
STR_DEG_AVC_up_GO <- enrichGO(gene = STR_DEG_AVC_up, OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP", pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.05)

(dotplot(OB_DEG_AVC_up_GO,title="OB region") + dotplot(CTX_1_DEG_AVC_up_GO,title="CTX_1 region")) /
  (dotplot(CTX_2_DEG_AVC_up_GO,title="CTX_2 region") + dotplot(HC_DEG_AVC_up_GO,title="HC region")) /
  (dotplot(MB_DEG_AVC_up_GO,title="MB region") + dotplot(CB_DEG_AVC_up_GO,title="CB region")) /
  (dotplot(HB_DEG_AVC_up_GO,title="HB region") + dotplot(TH_DEG_AVC_up_GO,title="TH region")) /
  (dotplot(HY_DEG_AVC_up_GO,title="HY region") + dotplot(STR_DEG_AVC_up_GO,title="STR region"))

OB_DEG_AVC_up_GO
OB_DEG_AVC_up_GO_3 <- OB_DEG_AVC_up_GO@result[c(1,3,4),c("Description","pvalue","Count")]
OB_DEG_AVC_up_GO_3$log10Pvalue <- -log10(OB_DEG_AVC_up_GO_3$pvalue)
OB_DEG_AVC_up_GO_3$subtype <- "OB"
OB_DEG_AVC_up_GO_3

CTX_1_DEG_AVC_up_GO
CTX_1_DEG_AVC_up_GO_3 <- CTX_1_DEG_AVC_up_GO@result[c(1,2,5),c("Description","pvalue","Count")]
CTX_1_DEG_AVC_up_GO_3$log10Pvalue <- -log10(CTX_1_DEG_AVC_up_GO_3$pvalue)
CTX_1_DEG_AVC_up_GO_3$subtype <- "CTX_1"
CTX_1_DEG_AVC_up_GO_3

CTX_2_DEG_AVC_up_GO
CTX_2_DEG_AVC_up_GO_3 <- CTX_2_DEG_AVC_up_GO@result[c(2,3,5),c("Description","pvalue","Count")]
CTX_2_DEG_AVC_up_GO_3$log10Pvalue <- -log10(CTX_2_DEG_AVC_up_GO_3$pvalue)
CTX_2_DEG_AVC_up_GO_3$subtype <- "CTX_2"
CTX_2_DEG_AVC_up_GO_3

HC_DEG_AVC_up_GO
HC_DEG_AVC_up_GO_3 <- HC_DEG_AVC_up_GO@result[c(2,3,4),c("Description","pvalue","Count")]
HC_DEG_AVC_up_GO_3$log10Pvalue <- -log10(HC_DEG_AVC_up_GO_3$pvalue)
HC_DEG_AVC_up_GO_3$subtype <- "HC"
HC_DEG_AVC_up_GO_3

MB_DEG_AVC_up_GO
MB_DEG_AVC_up_GO_3 <- MB_DEG_AVC_up_GO@result[c(1,5,9),c("Description","pvalue","Count")]
MB_DEG_AVC_up_GO_3$log10Pvalue <- -log10(MB_DEG_AVC_up_GO_3$pvalue)
MB_DEG_AVC_up_GO_3$subtype <- "MB"
MB_DEG_AVC_up_GO_3

CB_DEG_AVC_up_GO
CB_DEG_AVC_up_GO_3 <- CB_DEG_AVC_up_GO@result[c(1,2,9),c("Description","pvalue","Count")]
CB_DEG_AVC_up_GO_3$log10Pvalue <- -log10(CB_DEG_AVC_up_GO_3$pvalue)
CB_DEG_AVC_up_GO_3$subtype <- "CB"
CB_DEG_AVC_up_GO_3

HB_DEG_AVC_up_GO
HB_DEG_AVC_up_GO_3 <- HB_DEG_AVC_up_GO@result[c(1,8,10),c("Description","pvalue","Count")]
HB_DEG_AVC_up_GO_3$log10Pvalue <- -log10(HB_DEG_AVC_up_GO_3$pvalue)
HB_DEG_AVC_up_GO_3$subtype <- "HB"
HB_DEG_AVC_up_GO_3

TH_DEG_AVC_up_GO
TH_DEG_AVC_up_GO_3 <- TH_DEG_AVC_up_GO@result[c(2,4,6),c("Description","pvalue","Count")]
TH_DEG_AVC_up_GO_3$log10Pvalue <- -log10(TH_DEG_AVC_up_GO_3$pvalue)
TH_DEG_AVC_up_GO_3$subtype <- "TH"
TH_DEG_AVC_up_GO_3

HY_DEG_AVC_up_GO
HY_DEG_AVC_up_GO_3 <- HY_DEG_AVC_up_GO@result[c(2,3,7),c("Description","pvalue","Count")]
HY_DEG_AVC_up_GO_3$log10Pvalue <- -log10(HY_DEG_AVC_up_GO_3$pvalue)
HY_DEG_AVC_up_GO_3$subtype <- "HY"
HY_DEG_AVC_up_GO_3

STR_DEG_AVC_up_GO
STR_DEG_AVC_up_GO_3 <- STR_DEG_AVC_up_GO@result[c(1,4,5),c("Description","pvalue","Count")]
STR_DEG_AVC_up_GO_3$log10Pvalue <- -log10(STR_DEG_AVC_up_GO_3$pvalue)
STR_DEG_AVC_up_GO_3$subtype <- "STR"
STR_DEG_AVC_up_GO_3

Go_DEG_up <- rbind(OB_DEG_AVC_up_GO_3,CTX_1_DEG_AVC_up_GO_3,CTX_2_DEG_AVC_up_GO_3,
                   HC_DEG_AVC_up_GO_3,MB_DEG_AVC_up_GO_3,CB_DEG_AVC_up_GO_3,HB_DEG_AVC_up_GO_3,
                   TH_DEG_AVC_up_GO_3,HY_DEG_AVC_up_GO_3,STR_DEG_AVC_up_GO_3)
Go_DEG_up$item <- row.names(Go_DEG_up)
Go_DEG_up$subtype <- factor(Go_DEG_up$subtype,levels = c("OB","CTX_1","CTX_2","HC","MB",
                                                         "CB","HB","TH","HY","STR"))
Go_DEG_up$pathways <- paste0(Go_DEG_up$item," _ ",Go_DEG_up$Description)
ggbarplot(Go_DEG_up, x="pathways", y="log10Pvalue", fill = "subtype", 
          color = "subtype",
          palette =  cols,
          sort.val = "asc",
          sort.by.grodowns=TRUE, 
          x.text.angle=45, 
          xlab = NULL) 