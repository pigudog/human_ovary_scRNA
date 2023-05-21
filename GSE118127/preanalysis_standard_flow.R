library(rhdf5)
library(tidyverse)
library(patchwork)
library(stringr)
library(Matrix)
library(hdf5r)
suppressMessages(library(plotly))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(topGO))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(gplots))
suppressMessages(library(genefilter))
suppressMessages(library(scran))
library(DoubletFinder)
rm(list = ls())
rm_doublet <- function(name=NULL,input=NULL,dim.usage=30,auto="false") {
  ## get the path of the matric 
  ## The matrix passed in here is in the 10X format
  inpath <- list.files(path = input,pattern = name,full.names = T)
  if (auto=="true") {
    inpath <- paste0(inpath,"/","04.Matrix/")#这是另一种矩阵路径结构
  }
  
  ## before doubeltFinder, you need cluster, just choose the default parameters
  ## This section can be written more succinctly using the pipe character %>%
  EC <- Read10X_h5(inpath)
  EC <- CreateSeuratObject(EC, min.cells = 3, min.features = 100)
  EC <- SCTransform(EC)
  EC <- RunPCA(EC,verbose=F)
  EC <- RunUMAP(EC, dims = 1:dim.usage) 
  EC <- FindNeighbors(EC, dims = 1:dim.usage) %>% FindClusters(resolution = 0.3) 
  ##DoubletFinder去双胞的标准流程，封装成一个函数
  Find_doublet <- function(data){
    ### pK Identification (no ground-truth)
    sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = T)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    ### Homotypic Doublet Proportion Estimate
    DoubletRate = ncol(pbmc)*8*1e-6
    homotypic.prop <- modelHomotypic(data$seurat_clusters)   # 最好提供celltype
    nExp_poi <- round(DoubletRate*ncol(data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ### Run DoubletFinder with varying classification stringencies
    data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    return(data)
  }
  
  ##调用上面写好的函数，返回的是一个Seurat对象，meta.data信息里会有双胞信息，需要自己手动删除
  EC<-Find_doublet(EC)
  EC<-subset(EC,subset=doublet_info=="Singlet")
  EC@meta.data$library = name #顺便打上这个样本的label
  
  ##生成的Seurat对象有个问题，会在meta.data里多了很多pANN_开头的列，需要手动删除
  c <- grep("pANN_",colnames(EC@meta.data))
  EC@meta.data <- EC@meta.data[,-c]
  ##输出此样本的细胞数
  print(paste0(name," cells: ", length(EC@meta.data$orig.ident)," is read!"," Time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
  return(EC)
}

folders=list.files('../raw_data/','^GSM')
write.table(folders, file ="/home/data/t070421/GSE118127/scRNA/infile.txt", sep ="", row.names =FALSE, col.names =FALSE, quote =FALSE)

## set parameters
infile = "/home/data/t070421/GSE118127/scRNA/infile.txt"
input = "/home/data/t070421/GSE118127/raw_data/"
samples_name = c('1.1','1.2','1.3','1.4','1.5','1.6','1.7','1.8','3.13','3.14','3.15','3.16','3.17','3.18','3.5','3.6','01',
                 '10','11','12','13','145','02','03','04','05','6a','07','8a','8b','C1')
merge_seob <- function(infile,input,batch=NULL){
  filelist <- readLines(infile)
  seob_list <- lapply(filelist, rm_doublet, input=input)
  n <- length(seob_list)
  for(i in 1:n){
    #给细胞barcode加个前缀，防止合并后barcode重名
    seob_list[[i]] <- RenameCells(seob_list[[i]], add.cell.id = samples_name[i])   
  }
  names(seob_list) <- samples_name
  all <- merge(seob_list[[1]], seob_list[c(2:n)])
  if (!is.null(batch)) {
    all$batch <- batch
  }
  return(all)
}

ovary <- merge_seob(infile=infile,input=input)
table(ovary$library)
save(ovary,file = 'ovary_after_doublet.Rdata')


## Print the overview of cell counts per sample
table(ovary@meta.data$library)
## Infer the run (old or June) from the sample name and save it
indexJuneRun <- grep("-", ovary@meta.data$library, fixed = TRUE) 
origin.run <- orgin.exper <- rep("OldRun", times=length(ovary@active.ident))
origin.run[indexJuneRun] = "JuneRun"
#ovary <- AddMetaData(object = ovary, metadata = origin.run, col.name = "origin.run")
ovary@meta.data$origin.run <- origin.run 
table(ovary@meta.data$origin.run)
# JuneRun  OldRun 
# 44343   11758

## Create a meta.data field with patient info and sample info
ovary@meta.data$orig.ident <- ovary@meta.data$library
table(ovary@meta.data$orig.ident)
ovary_index <- ovary@meta.data$orig.ident
ovary_patient <- ovary@meta.data$orig.ident
ovary_sample <- ovary@meta.data$orig.ident
# 1.1
ovary.1.1 <- grep("1-1", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "1.1"
ovary_patient[ovary.1.1] = "P9"
ovary_sample[ovary.1.1] = "stroma"
# 1.2
ovary.1.1 <- grep("1-2", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "1.2"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "stroma"
# 1.3
ovary.1.1 <- grep("1-3", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "1.3"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "stroma"
# 1.4
ovary.1.1 <- grep("1-4", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "1.4"
ovary_patient[ovary.1.1] = "P9"
ovary_sample[ovary.1.1] = "stroma"
# 1.5
ovary.1.1 <- grep("1-5", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "1.5"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "Follicle 2-5mm"
# 1.6
ovary.1.1 <- grep("1-6", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "1.6"
ovary_patient[ovary.1.1] = "P9"
ovary_sample[ovary.1.1] = "stroma"
# 1.7
ovary.1.1 <- grep("1-7", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "1.7"
ovary_patient[ovary.1.1] = "P9"
ovary_sample[ovary.1.1] = "stroma"
# 1.8
ovary.1.1 <- grep("1-8", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "1.8"
ovary_patient[ovary.1.1] = "P9"
ovary_sample[ovary.1.1] = "stroma"
# 3-13
ovary.1.1 <- grep("3-13", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "3.13"
ovary_patient[ovary.1.1] = "P3"
ovary_sample[ovary.1.1] = "Follicle 1-2mm"
# 3-14
ovary.1.1 <- grep("3-14", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "3-14"
ovary_patient[ovary.1.1] = "P3"
ovary_sample[ovary.1.1] = "Follicle 1-2mm"
# 3-15
ovary.1.1 <- grep("3-15", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "3-15"
ovary_patient[ovary.1.1] = "P3"
ovary_sample[ovary.1.1] = "Follicle 1-2mm"
# 3-16
ovary.1.1 <- grep("3-16", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "3-16"
ovary_patient[ovary.1.1] = "P3"
ovary_sample[ovary.1.1] = "Follicle 2-5mm"
# 3-17
ovary.1.1 <- grep("3-17", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "3-17"
ovary_patient[ovary.1.1] = "P3"
ovary_sample[ovary.1.1] = "stroma"
# 3.18
ovary.1.1 <- grep("3-18", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "3-18"
ovary_patient[ovary.1.1] = "P2"
ovary_sample[ovary.1.1] = "stroma"
# 3.5
ovary.1.1 <- grep("3-5", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "3-5"
ovary_patient[ovary.1.1] = "P0"
ovary_sample[ovary.1.1] = "stroma"
# 3.6
ovary.1.1 <- grep("3-6", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "3-6"
ovary_patient[ovary.1.1] = "P0"
ovary_sample[ovary.1.1] = "stroma"
# 01
ovary.1.1 <- grep("sample1", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "01"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "stroma"
# 10
ovary.1.1 <- grep("sample10", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "10"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "Follicle 1-2mm"
# 11
ovary.1.1 <- grep("sample11", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "11"
ovary_patient[ovary.1.1] = "P9"
ovary_sample[ovary.1.1] = "stroma"
# 12
ovary.1.1 <- grep("sample12", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "12"
ovary_patient[ovary.1.1] = "P0"
ovary_sample[ovary.1.1] = "stroma"
# 13
ovary.1.1 <- grep("sample13", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "13"
ovary_patient[ovary.1.1] = "P0"
ovary_sample[ovary.1.1] = "stroma"
# 145
ovary.1.1 <- grep("sample145", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "145"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "Follicle 2-5mm"
# 02
ovary.1.1 <- grep("sample2", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "02"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "stroma"
# 03
ovary.1.1 <- grep("sample3", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "03"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "stroma"
# 04
ovary.1.1 <- grep("sample4", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "04"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "stroma"
# 05
ovary.1.1 <- grep("sample5", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "05"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "stroma"
# 6a
ovary.1.1 <- grep("sample6a", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "6a"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "stroma"
# 07
ovary.1.1 <- grep("sample7", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "07"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "Follicle 1-2mm"
# 8a
ovary.1.1 <- grep("sample8a", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "8a"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "Follicle 1-2mm"
# 8b
ovary.1.1 <- grep("sample8b", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "8b"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "Follicle 1-2mm"
# C1
ovary.1.1 <- grep("sampleC1", ovary@meta.data$orig.ident, fixed = TRUE) 
ovary_index[ovary.1.1] = "C1"
ovary_patient[ovary.1.1] = "P7"
ovary_sample[ovary.1.1] = "Follicle 1-2mm"
ovary@meta.data$ovary_index <- ovary_index
ovary@meta.data$ovary_patient <- ovary_patient
ovary@meta.data$ovary_sample <- ovary_sample
table(ovary@meta.data$ovary_index)
table(ovary@meta.data$ovary_patient)
table(ovary@meta.data$ovary_sample)
# > table(ovary@meta.data$ovary_index)
# 
# 01   02   03   04   05   07  1.1  1.2  1.3  1.4  1.5  1.6  1.7  1.8   10   11 
# 879  738  107  554  491  720 5437 1710 3598 2948 2561 2087 4349 2698  511 1224 
# 12   13  145 3-14 3-15 3-16 3-17 3-18  3-5  3-6 3.13   6a   8a   8b   C1 
# 843  434 1184 1285  818 5125 4533 5954  562  164  514 2546  528  759  240 
# > table(ovary@meta.data$ovary_patient)
# 
# P0    P2    P3    P7    P9 
# 2003  5954 12275 17126 18743 
# > table(ovary@meta.data$ovary_sample)
# 
# Follicle 1-2mm Follicle 2-5mm         stroma 
# 5375           8870          41856

#####################################################################
### Standard pre-processing workflow
# 1. Screen
##mitochondrion, ribosome、erythrocyte
ovary[["percent.mt"]] <- PercentageFeatureSet(ovary, pattern = "^MT-")
ovary[["percent.rb"]] <- PercentageFeatureSet(ovary, pattern = "^RP[SL]")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(ovary))
ovary[["percent.HB"]]<-PercentageFeatureSet(ovary, features=HB.genes)
### QC violin picture
# possible theme
theme.set2 = theme(axis.title.x=element_blank())
# the elements of picture
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
group = "ovary_sample"
# vlnplot before qc
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(ovary, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
ggsave("./preanalysis/vlnplot_before_qc.png", plot = violin, width = 12, height = 8)
### set QC strandards
minGene=100
maxGene=5000
maxUMI=20000
pctMT=10
pctHB=1
### vlnplot after QC
ovary <- subset(ovary, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB <pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(ovary, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)     
ggsave("./preanalysis/vlnplot_after_qc.png", plot = violin, width = 12, height = 8) 


cellinfo <- subset(ovary@meta.data, select = c("orig.ident", "percent.mt", "percent.rb", "percent.HB","origin.run",
                                                  "ovary_index","ovary_patient","ovary_sample"))
ovary <- CreateSeuratObject(ovary@assays$RNA@counts, meta.data = cellinfo)
save(ovary,file = 'ovary_after_screen.Rdata')

# 2. harmony for integrate
library(harmony)
library(tidyverse)
library(patchwork)
### standardize data with SCT v2!
ovary <- SCTransform(ovary, vst.flavor = "v2")
DefaultAssay(ovary)

### PCA
ovary <- RunPCA(ovary, npcs=50, verbose=FALSE)

ovary <- RunHarmony(ovary, 
                    group.by.vars=c("ovary_patient","origin.run"),
                    reduction = "pca", 
                    assay.use="SCT", 
                    max.iter.harmony = 20) 
# group.by.vars: The parameter is to set which group to integrate by
# max.iter.harmony: Set the number of iterations, default is 10. When running RunHarmony, results will indicate how many iterations have elapsed before convergence is complete.

elowplot <- ElbowPlot(ovary, ndims = 50)
ggsave("./preanalysis/elbow_plot.png", plot = elowplot, width = 12, height = 8) 

pc.num=1:30
ovary <- RunUMAP(ovary, reduction="harmony", dims=pc.num)
ovary <- FindNeighbors(ovary, reduction = "harmony",
                       dims = 1:30)
ovary<-FindClusters(ovary,resolution=0.5)
p <- DimPlot(ovary, group.by = "ovary_patient")
ggsave("./preanalysis/UMAP_ovary_patient_primality_0.5.png", p, width = 8, height = 6)
p <- DimPlot(ovary, group.by = "seurat_clusters")
ggsave("./preanalysis/UMAP_ovary_seurat_clusters_primality_0.5.png", p, width = 8, height = 6)
p <- DimPlot(ovary, group.by = "ovary_patient", split.by = "ovary_patient", ncol = 4)
ggsave("./preanalysis/UMAP_ovary_patient_spilit_primality_0.5.png", p, width = 18, height = 12)
save(ovary,file = 'ovary_after_harmony.Rdata')

###########################################################################################
# Find HVG and cluster
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(stringr)
rm(list=ls())
options(stringsAsFactors = F)
ovary.markers <- FindAllMarkers(object = ovary,assay = "SCT", only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
# wriet to decidual/HVGmarkers_res_0.5
write.csv(ovary.markers,file=paste0("./preanalysis",'/HVGmarkers_res_0.5.csv'))
save(ovary,ovary.markers,file = 'ovary_HVG.Rdata')

library(dplyr) 
p<-DimHeatmap(object = ovary,dims = 1:15)
ggsave(p,filename=paste0("preanalysis",'/ovary_dimheatmap.pdf'),width = 24,height = 18)
top10 <- ovary.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
p <- DoHeatmap(ovary,top10$gene,size=3)
ggsave(p,filename=paste0("preanalysis",'/top10_heatmap.png'),width = 24,height = 18)
# Specify genes  - Mac
genes_to_check = c( "FCN1","MS4A7","CD14")
# featureplot
p <- FeaturePlot(ovary, features = genes_to_check)
ggsave(p,filename=paste0("preanalysis",'/fearureplot_sepcify_mac.png'),width = 16,height = 10)
# All on Dotplot 
p <- DotPlot(ovary, features = genes_to_check) + coord_flip()
ggsave(p,filename=paste0("preanalysis",'/dotplot_sepcify_specific_MAC.pdf'),width = 16,height = 12)
# 13 macrophage


#################################################################################################
# cell annotation
# Specify genes  - Mac
genes_to_check = c("AMH","HSD17B1","SERPINE2","GSTA1","DCN","LUM","TAGLN","RGS5","VWF","CLDN5","CD53","CXCR4","PTPRC")
# featureplot
p <- FeaturePlot(ovary, features = genes_to_check)
ggsave(p,filename=paste0("preanalysis",'/fearureplot_sepcify.png'),width = 16,height = 10)
# All on Dotplot 
p <- DotPlot(ovary, features = genes_to_check) + coord_flip()
ggsave(p,filename=paste0("preanalysis",'/dotplot_sepcify_specific.png'),width = 16,height = 12)


#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:15,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(8,13),2]='Immune' 
celltype[celltype$ClusterID %in% c(6,10,11),2]='Granulosa'
celltype[celltype$ClusterID %in% c(0,1,5,7,9,12),2]='Theca&stroma'
celltype[celltype$ClusterID %in% c(3),2]='Smooth muscle'

celltype[celltype$ClusterID %in% c(2,4,14,15),2]='Endo'


head(celltype)
celltype 
table(celltype$celltype)

new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(ovary)
ovary <- RenameIdents(ovary, new.cluster.ids)
save(ovary,file = 'ovary_after_cluster.Rdata')
# nolegend
p<-DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(p,filename=paste0("preanalysis",'/umap_firstcluster.pdf'),width = 16,height = 12)
# legend
p<-DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(p,filename=paste0("preanalysis",'/umap_firstcluster_legend.pdf'),width = 16,height = 12)

sample_table <- as.data.frame(table(ovary@meta.data$ovary_sample,ovary@active.ident))
names(sample_table) <- c("Samples","celltype","CellNumber")
# color
library(ggsci)
library(Cairo)
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.title = element_text(lineheight=2, face="bold", hjust=2, size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20)
  )+labs(y="Percentage")+RotatedAxis()
ggsave(plot_sample,filename=paste0("preanalysis",'/percentage_ovary.pdf'),width = 6,height = 8)
#####################################################################################
# re_cluster
ovary_immune <- ovary[,ovary@active.ident %in% "Immune"]
sce <-ovary_immune
sce <- SCTransform(sce, vst.flavor = "v2")
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
sce <- RunUMAP(object = sce, dims = 1:30, do.fast = TRUE)
res.used <- 0.5
set.seed(1314)
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(object = sce, verbose = T, resolution = res.used)
p <- DimPlot(sce,reduction = "umap",label=T)
ggsave(p,filename=paste0("preanalysis",'/dimplot_immune.png'),width = 16,height = 10)
save(sce,file = 'sce_immune_before_cluster.Rdata')

# Specify genes  - Mac - NK - T -B - complement - Innate -Mast -dc
genes_to_check = c( "FCN1","MS4A7","CD14","NCAM1","GZMB","CD3D","CD3G","MS4A1","IGKC","C3","C1QA",
                    "CD68","IFI30","TPSB2","TPSAB1","CLEC9A")
# featureplot
p <- FeaturePlot(sce, features = genes_to_check)
ggsave(p,filename=paste0("preanalysis",'/fearureplot_sepcify_mac_immune.png'),width = 16,height = 10)
# All on Dotplot 
p <- DotPlot(sce, features = genes_to_check) + coord_flip()
ggsave(p,filename=paste0("preanalysis",'/dotplot_sepcify_specific_immune.png'),width = 16,height = 12)

#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:14,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(1,5,9,11),2]='Mac' 
celltype[celltype$ClusterID %in% c(6),2]='NK'
celltype[celltype$ClusterID %in% c(0,2,3,4,7,8,10,13),2]='T'
celltype[celltype$ClusterID %in% c(12,14),2]='B'

# celltype[celltype$ClusterID %in% c(3),2]='Complement'
# celltype[celltype$ClusterID %in% c(9),2]='Innate'

head(celltype)
celltype 
table(celltype$celltype)

new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
save(sce,file = 'sce_immune_after_cluster.Rdata')
# legend
p<-DimPlot(sce, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(p,filename=paste0("preanalysis",'/umap_firstcluster_legend_immune.pdf'),width = 16,height = 12)


# boxplot
# level0
library(ggpubr)
library(ggsignif)
# 设置可能用到的主题
theme.set2 = theme(axis.title.x=element_blank())
mac <- sce[,sce@active.ident %in% "Mac"]
group = "ovary_sample"
plot.featrures = c("CD68","ENPP2","S100A4")
my_comparisons <- list(c("Follicle 1-2mm","Follicle 2-5mm"))
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(mac, group.by=group, pt.size = 0,
                       features = plot.featrures[i])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    stat_compare_means(comparisons = my_comparisons) +
    ylim(-1, 5)+
    geom_boxplot()
  theme.set2 + NoLegend()}
p <- wrap_plots(plots = plots, nrow=2)     
ggsave(p,filename=paste0("preanalysis",'/umap_MAC_CD68_boxplot.pdf'),width = 16,height = 12)

sample_table <- as.data.frame(table(sce@meta.data$ovary_sample,sce@active.ident))
names(sample_table) <- c("Samples","celltype","CellNumber")
# color
library(ggsci)
library(Cairo)
color = c(pal_d3("category20")(20),
          pal_d3("category20b")(20),
          pal_d3("category20c")(20),
          pal_d3("category10")(10))

plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.title = element_text(lineheight=2, face="bold", hjust=2, size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20)
  )+labs(y="Percentage")+RotatedAxis()
ggsave(plot_sample,filename=paste0("preanalysis",'/percentage_immune.pdf'),width = 6,height = 8)





