# 9. loading packages
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(ggsci)
library(VennDiagram)

# 9.1 load ST-seq data
data <- readRDS("CM_ST_dataset.rds") 
meta <- data@meta.data

Control_meta <- subset(meta,sample_type=="Control")
Model_meta <- subset(meta,sample_type=="Model")
ART_meta <- subset(meta,sample_type=="ART")

Control_data <- as.matrix(data@assays$SCT@data[,as.character(row.names(Control_meta))])
Model_data <- as.matrix(data@assays$SCT@data[,as.character(row.names(Model_meta))])
ART_data <- as.matrix(data@assays$SCT@data[,as.character(row.names(ART_meta))])

color <- c(
"Astro"="#00a988",
"CPC"="#3f5089",
"Endo"="#fc9572",
"Fibro"="#888eb5",
"Immune"="#8bd3c3",
"Neuron"="#e55130",
"Oligo"="#32bedb"
)

# 9.2 cell-cell communication in each group

cellchat <- createCellChat(object = Control_data, meta = Control_meta, group.by = "SingleR.labels")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 40) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ST_Control_cellchat.RDS")

cellchat <- createCellChat(object = Model_data, meta = Model_meta, group.by = "SingleR.labels")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 40) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ST_Model_cellchat.RDS")

cellchat <- createCellChat(object = ART_data, meta = ART_meta, group.by = "SingleR.labels")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 40) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat) # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) #Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat) #Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use=color)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use=color)
saveRDS(cellchat, "ST_ART_cellchat.RDS")

# 9.3 Merge cell-cell communication in ST-seq dataset
control<-readRDS( "ST_Control_cellchat.RDS")
model<-readRDS( "ST_Model_cellchat.RDS") 
art<-readRDS( "ST_ART_cellchat.RDS")

object.list <- list(Control=control,Model=model,ART=art)
cellchat <- mergeCellChat(object.list, add.names = compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
netVisual_heatmap(cellchat,color.use=color)
netVisual_heatmap(cellchat,comparison=c("Control","Model"),title.name = "Different number of interactions of MVC\nin ST-seq dataset",font.size=16,font.size.title = 18)+ggtitle("ModvsCon")
netVisual_heatmap(cellchat,comparison=c("Control","ART"),title.name = "Different number of interactions of AVC\nin ST-seq dataset",font.size=16,font.size.title = 18)+ggtitle("ARTvsCon")

# 9.4 Different enriched pathways in ST-seq datasets

ST_MvsC_Path=c("MHC-I","ICAM","GALECTIN","ITGAL-ITGB2","COMPLEMENT","CD45","LCK","CCL","CXCL","MIF")
ST_AvsC_Path=c("MHC-I","GALECTIN","NEGR","ICAM","CD45","SELL","ITGAL-ITGB2","COMPLEMENT","CD22","SELPLG","PTPRM","IFN-II","ANNEXIN","NPNT","PARs","CCL","MIF","CXCL","IGF","APP")
rankNet(ST_cellchat, mode = "comparison", stacked = T,do.stat =TRUE,comparison = c(1,2,3),color.use=c("#3B4992","#EE0000","#008B45" ),signaling =unique(c(ST_MvsC_Path,ST_AvcC_Path)))+coord_flip()

# 9.5 different pathway overlap across different groups in scRNA-seq and ST-seq datastes

venn.diagram(x=list(SC_MvsC=SC_MvsC_Path,SC_AvsC=SC_AvsC_Path,ST_MvsC=ST_MvsC_Path,ST_AvsC=ST_AvsC_Path),filename=NULL,fill=pal_jama()(4),cex = 1.5, cat.col = 'black', cat.cex = 1.5, cat.fontface = "bold", margin = 0.05,rotation.degree = 0)
grid.draw(plot)

# 9.6 using stLearn to perform LR communication in ST-seq datasets by Python

import stlearn as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

data = st.Read10X("BM2_ST_P/outs/")
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath="BM2_ST_P/outs/spatial/tissue_hires_image.png",
             library_id="BM2_ST_P", visium=True)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data) # NOTE: no log1p
spot_mixtures = pd.read_csv('BM2_P_meta.csv', index_col=0)
labels = spot_mixtures.loc[:,"SingleR.labels"].values.astype(str)

print(labels)
print(spot_mixtures)
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values)) # Check is in correct order
data.obs['cell_type'] = labels # Adding the dominant cell type labels per spot
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures # Adding the cell type scores
st.pl.cluster_plot(data, use_label='cell_type',fname="BM2_P_celltype.pdf")
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))
Manual_lr = np.array(['H2-D1_Cd8b1','H2-K1_Cd8a','H2-M3_Cd8b1','H2-Q4_Cd8a','H2-Q6_Cd8a','H2-Q7_Cd8b1','H2-T22_Cd8a','H2-T23_Klrd1','H2-T23_Klrc1','H2-T23_Cd8b1'])
st.tl.cci.run(data, Manual_lr,
                  min_spots = 1, #Filter out any LR pairs with no scores for less than min_spots
                  distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=48, # Number of CPUs for parallel. If None, detects & use all available.
                  )

lr_info = data.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
print('\n', lr_info)

st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')

data.write("BM2_P_backup.h5ad")
stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
#"Ccl3_Ccr5","Ccl4_Ccr5","Ccl5_Ccr5","Ccl7_Ccr2","Ccl8_Ccr5",'Ifng_Ifngr1', 'Ifng_Ifngr2','Vcam1_Itga4','Vcam1_Itgb1','Icam2_Itgam', 'Icam2_Itgb2','Icam2_Itgal','Icam1_Itgam', 'Icam1_Itgb2','Icam1_Itgal'
best_lr = ["Ccl3_Ccr5","Ccl4_Ccr5","Ccl5_Ccr5","Ccl7_Ccr2","Ccl8_Ccr5",'Ifng_Ifngr1', 'Ifng_Ifngr2','Vcam1_Itga4','Vcam1_Itgb1','Icam2_Itgam', 'Icam2_Itgb2','Icam2_Itgal','Icam1_Itgam', 'Icam1_Itgb2','Icam1_Itgal']

fig, axes = plt.subplots(ncols=2, figsize=(8,6))
st.pl.lr_result_plot(data, use_result='-log10(p_adjs)', use_lr=best_lr[14], show_color_bar=False, ax=axes[1])
st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=best_lr[14], show_color_bar=False, ax=axes[0])
axes[1].set_title(f'{best_lr[14]} -log10(p_adjs)')
axes[0].set_title(f'{best_lr[14]} lr_sig_scores')

Manual_lr = ['H2-D1_Cd8b1','H2-K1_Cd8a','H2-M3_Cd8b1','H2-Q4_Cd8a','H2-Q6_Cd8a','H2-Q7_Cd8b1','H2-T22_Cd8a','H2-T23_Klrd1','H2-T23_Klrc1','H2-T23_Cd8b1']
a=Manual_lr[3]

fig, axes = plt.subplots(ncols=2, figsize=(8,6))
st.pl.lr_result_plot(data, use_result='-log10(p_adjs)', use_lr=a, show_color_bar=True, ax=axes[1])
st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[0])
axes[1].set_title(f'{a} -log10(p_adjs)')
axes[0].set_title(f'{a} lr_sig_scores')





