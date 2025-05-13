#Extraction of MFAP5-positive cells
fibro<-readRDS("fibro.rds")
MF<-subset(fibro,MFAP5>1,slot="data")
MF<-NormalizeData(MF)
MF<-FindVariableFeatures(MF)
MF<-ScaleData(MF)
MF<- RunPCA(MF,features = VariableFeatures(object = MF)) 
MF<-RunHarmony(MF,reduction="pca", group.by.vars="orig.ident",reduction.save="harmony")
MF<- FindNeighbors(MF, dims = 1:5,reduction = "harmony")
MF<- FindClusters(MF, resolution = 0.2 )
Idents(MF)<-"RNA_snn_res.0.2"
MF<-RunUMAP(MF, dims = 1:5, reduction = "harmony", reduction.name = "umap")
DimPlot(MF)
p<-DimPlot(MF,group.by = "group",pt.size = 1)
ggsave("UMAP.tiff", plot = p, bg = "transparent", width =6.2, height =5, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)
markers<-FindAllMarkers(object = MF,test.use = "wilcox",only.pos = T,logfc.threshold = 0.25)
write.csv(markers,"marker.csv")

#Cells were annotated based on the expression of canonical marker genes and the results of GO and KEGG pathway enrichment analyses.

MF$celltype<-plyr::mapvalues(MF$RNA_snn_res.0.2,from = 0:4, to=c("Wnt-responsive fibroblasts",
                                                                 "Inflammatory fibroblasts",
                                                                 "Wnt-responsive fibroblasts",
                                                                 "TGFβ/BMP-responsive fibroblasts",
                                                                 "Pathological fibroblasts"))#待修改
new.cluster.ids <-c("Wnt-responsive fibroblasts",
                    "Inflammatory fibroblasts",
                    "Wnt-responsive fibroblasts",
                    "TGFβ/BMP-responsive fibroblasts",
                    "Pathological fibroblasts")
names(new.cluster.ids) <- levels(MF)
MF <- RenameIdents(MF, new.cluster.ids)
#Plot
p<-DimPlot(MF, reduction = "umap",pt.size = 0.1,group.by = "group")+scale_color_npg()+labs(title = "Cell cluster")+
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust =0.03))+
  theme(
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 14,face = "bold")
  )
ggsave("umap.tiff", plot = p, bg = "transparent", width =6.2, height =5, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)



#Cell proportion analysis
ratio<-table(MF$celltype,MF$group)

ratio<-as.data.frame(ratio)

write.csv(ratio,"ratio.csv")
ratio1<-ratio
for (i in 1:nrow(ratio)) {
  if (ratio$Var2[i] == "Control") {
    ratio1$Freq[i] <- ratio$Freq[i] / sum(ratio$Freq[ratio$Var2 == "Control"])
  } else {
    
    ratio1$Freq[i] <- ratio$Freq[i] / sum(ratio$Freq[ratio$Var2 == "IPF"])
  }
  
}
}

colors<- c("#DC0000FF","#4DBBD5FF","#00A087FF","#3C5488FF")
library(ggalluvial)

p<-ggplot(ratio1,aes(x=Var2,
                     stratum=Var1,
                     alluvium=Var1,
                     y=Freq,
                     fill = Var1,
                     label=Var1))+
  scale_y_continuous(
    label = scales::percent_format(), 
    expand = c(0.02, 0)  
  ) +
  scale_fill_manual(values = as.character(colors),breaks = names(colors)) +  
  geom_stratum(width = 0.4, colour = "grey40")+
  theme(
    axis.title = element_blank(),  
    axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1),  
    
    panel.background = element_blank(),  
    axis.line = element_line(),  
    axis.ticks = element_line(), 
    legend.position = "right",
    axis.text = element_text(size = 12) 
  )+
  theme(
    legend.position = "right",       
    legend.title = element_blank(),  
    legend.text = element_text(face="bold",size = 10),  
    legend.key = element_rect(fill = "white", color = "black")  
  )+
  geom_text(aes(label = paste0(round(Freq*100, 2),"%")), 
            position = position_stack(vjust = 0.5), 
            size =10, color = "black")
ggsave("proportion.tiff", plot = p, bg = "transparent", width =6.2, height =5, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)
#Pseudotime analysis
library(monocle3)
MF.cds <- as.cell_data_set(MF)
MF.cds <- cluster_cells(cds = MF.cds, reduction_method = "UMAP")
MF.cds <- learn_graph(MF.cds, use_partition = TRUE,learn_graph_control = list(minimal_branch_len=3))
MF.cds=order_cells(MF.cds)
#Note: Save image 
plot_cells(MF.cds,color_cells_by = "pseudotime",cell_size = 0.1,group_label_size = 5,trajectory_graph_segment_size = 2.5,label_branch_points = F,cell_stroke =10)+theme_dr(xlength = 0.2, 
                                                                                                                                                                        ylength = 0.2,
                                                                                                                                                                           arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+ theme(panel.grid = element_blank(),
                                                                                                                                                                                                                                               axis.title = element_text(face = 2,hjust =0.03))
MF.cds <- estimate_size_factors(MF.cds)

## Add gene names into CDS
MF.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(MF[["RNA"]])

tgen<-graph_test(MF.cds,neighbor_graph = "principal_graph",cores = 6)
library(dplyr)
tgen_sig<-tgen%>%top_n(n=10,morans_I)%>%pull(gene_short_name)%>%as.character()
#Note: Save image 
plot_genes_in_pseudotime(MF.cds[tgen_sig,],color_cells_by="celltype",min_expr=1,ncol = 2)+scale_color_manual(values = as.character(colors),breaks = names(colors))




