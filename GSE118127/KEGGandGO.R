library(clusterProfiler)
# library(org.Mm.eg.db) ##load mouse
library(org.Hs.eg.db) ##load human
library(reticulate)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(tidyverse)
library(Matrix)
library( magrittr)
library(sp)
library(monocle3)
options(stringsAsFactors = FALSE)
rm(list=ls())
load(file = 'sce_immune_after_cluster.Rdata')
sce@meta.data$cell_type <- sce@active.ident

sce_MAC<- sce[,sce$cell_type %in% "Mac"]
# CTRL VS EM
# future::plan("multicore", workers = 8)
MAC_FOLLI <- FindMarkers(sce_MAC, ident.1 = "Follicle 2-5mm",ident.2 = "Follicle 1-2mm", 
                          group.by = 'ovary_sample',logfc.threshold = 0,min.pct = 0)
sig_dge.all <- subset(MAC_FOLLI, p_val_adj<0.05&abs(avg_log2FC)>1) 

ego_CC <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05) 

display_number = c(10, 10, 10)# These three number represent the number of pathway of BP、CC、MF
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

##integrate
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))
# color and plot
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))#这一步是必须的，为了让柱子按顺序显示，不至于很乱
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")#设定颜色

plotc <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
  geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
  scale_fill_manual(values = COLS) + ###颜色
  coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
  xlab("GO term") +
  ylab("Gene_Number") +
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

ggsave('./KEGGandGO/MAC_FOLLI_enrichGO.png', plotc, width = 12,height = 14)

# KEGG
genelist <- bitr(row.names(sig_dge.all), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',pvalueCutoff = 0.1)
p1 <- barplot(ekegg, showCategory=15)
p2 <- dotplot(ekegg, showCategory=15)
plotc = p1/p2
ggsave("./KEGGandGO/MAC_FOLLI_enrichKEGG.png", plot = plotc, width = 12, height = 14)
