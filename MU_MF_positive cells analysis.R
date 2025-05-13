setwd("your dir")


#Extraction of MF-positive cells
BLM_fibro<-subset(GSEBLM_TIME,RNA_snn_res.0.1%in% c(0,1) )
BLM_fibro<-NormalizeData(BLM_fibro) %>%FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%  ScaleData() %>%RunPCA(npcs = 30,verbose = T)
BLM_fibro<-RunHarmony(BLM_fibro,reduction="pca", group.by.vars="orig.ident",reduction.save="harmony")
BLM_fibro<-FindNeighbors(BLM_fibro,reduction = "harmony",dims = 1:20)
BLM_fibro<-FindClusters(BLM_fibro,resolution = 0.3)
BLM_fibro<- RunUMAP(BLM_fibro,reduction = "harmony",dims=1:20)
MF<-subset(BLM_fibro,Mfap5>1)
MF<-NormalizeData(MF) %>%FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%  ScaleData() %>%RunPCA(npcs = 30,verbose = T)
MF<-RunHarmony(MF,reduction="pca", group.by.vars="orig.ident",reduction.save="harmony")
MF<-FindNeighbors(MF,reduction = "harmony",dims = 1:20)
MF<-FindClusters(MF,resolution = 0.3)
MF<- RunUMAP(MF,reduction = "harmony",dims=1:20)
markers<-FindAllMarkers(object = MF,test.use = "wilcox",only.pos = T,logfc.threshold = 0.25)
#Gene profiles
list<-c("Col13a1",
        "Inmt",
        "Itga8",
        "Npnt",
        "Pi16",
        "Dcn",
        "Ly6c1",
        "Scara5",
        "Cthrc1",
        "Spp1",
        "Postn",
        "Inhba",
        "Ube2c",
        "Birc5",
        "Cdc20",
        "Nusap1",
        "Rsad2",
        "Cxcl10",
        "Ifit3",
        "Ifit1"
        )
p<-DotPlot(MF, features = list, cols = c("lightgrey", "blue"), 
           dot.scale =20)+
  scale_size_continuous(trans = "sqrt", range = c(0.1, 10)) + 
  scale_color_viridis(option = "E", direction = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 20)) +
  theme(axis.text.y = element_text(size = 20))+
  theme(legend.title =element_text(size = 20))
ggsave("dot.tiff", plot = p, bg = "transparent", width =14.5, height =6, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)
#Cell annotation
MF$cellType<-plyr::mapvalues(MF$RNA_snn_res.0.3,from = 0:4, to=c("Alveolar fibroblasts",
                                                                 "Adventitial fibroblasts",
                                                                 "Pathological fibroblasts",
                                                                 "Proliferative fibroblasts",
                                                                 "Inflammatory fibroblast"
))


new.cluster.ids <- c("Alveolar fibroblasts",
                     "Adventitial fibroblasts",
                     "Pathological fibroblasts",
                     "Proliferative fibroblasts",
                     "Inflammatory fibroblast"
)
names(new.cluster.ids) <- levels(MF)
MF<- RenameIdents(MF, new.cluster.ids)

#Plot
p<-DimPlot(MF, reduction = "umap",pt.size = 0.1)+scale_color_simpsons()+labs(title = "Cell cluster")+
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust =0.03))+
  theme(
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 14,face = "bold")
  )

ggsave("UMAP_MU.tiff", plot = p, bg = "transparent", width =7.8, height =5, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)












#Define the fibrotic state of mice under BLM treatment
library(scales)
library(grid)
library(cowplot)
MF@meta.data$stage<-1

for (i in 1:nrow(MF@meta.data)){
  MF@meta.data$stage[i]<-ifelse(stringr::str_sub(MF@meta.data$group[i],4,5)=="0_","Normal state",
                                ifelse(as.numeric(stringr::str_sub(MF@meta.data$group[i],4,5))<8,"Inflammatory phase",
                                       ifelse(as.numeric(stringr::str_sub(MF@meta.data$group[i],4,5))<28,"Fibrotic phase",
                                              ifelse(as.numeric(stringr::str_sub(MF@meta.data$group[i],4,5))<55,"Recovery phase")
                                       )
                                       
                                )
  )
  print(label_percent()(i/nrow(MF@meta.data)))
  
}
#Plot
p1<-DimPlot(MF,reduction = "umap",cells = WhichCells(MF, expression = stage=="Normal state"),pt.size = 0.1)+
  scale_color_manual(values = rep("#2166AC", length(levels(MF))))+
  annotation_custom(
    grob = roundrectGrob(x=0.5, y=0.5, width=1, height=1,
                         r = unit(0, "snpc"),
                         gp = gpar(col = "black", fill = NA, lwd = 2))
  )+
  theme(
    axis.line = element_blank(),  
    axis.ticks = element_blank(),  
    axis.text = element_blank(),     
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  theme(legend.position = "none")

p2<-DimPlot(MF,reduction = "umap",cells = WhichCells(MF, expression = stage=="Inflammatory phase"),pt.size = 0.1)+
  scale_color_manual(values = rep("#80D5D9", length(levels(MF))))+
  annotation_custom(
    grob = roundrectGrob(x=0.5, y=0.5, width=1, height=1,
                         r = unit(0, "snpc"),
                         gp = gpar(col = "black", fill = NA, lwd = 2))
  )+
  theme(
    axis.line = element_blank(),  
    axis.ticks = element_blank(),  
    axis.text = element_blank(),    
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  theme(legend.position = "none") 

p3<-DimPlot(MF,reduction = "umap",cells = WhichCells(MF, expression = stage=="Fibrotic phase"),pt.size = 0.1)+
  scale_color_manual(values = rep("#F7B27A", length(levels(MF))))+
  annotation_custom(
    grob = roundrectGrob(x=0.5, y=0.5, width=1, height=1,
                         r = unit(0, "snpc"),
                         gp = gpar(col = "black", fill = NA, lwd = 2))
  )+
  theme(
    axis.line = element_blank(),   
    axis.ticks = element_blank(),  
    axis.text = element_blank(),    
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  theme(legend.position = "none") 

p4<-DimPlot(MF,reduction = "umap",cells = WhichCells(MF, expression = stage=="Recovery phase"),pt.size = 0.1)+
  scale_color_manual(values = rep("#B2182B", length(levels(MF))))+
  annotation_custom(
    grob = roundrectGrob(x=0.5, y=0.5, width=1, height=1,
                         r = unit(0, "snpc"),
                         gp = gpar(col = "black", fill = NA, lwd = 2))
  )+
  theme(
    axis.line = element_blank(),   
    axis.ticks = element_blank(),  
    axis.text = element_blank(),    
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  theme(legend.position = "none")

p<-plot_grid(p1, p2,p3, p4, nrow = 1)

ggsave("abc.tiff", plot = p, bg = "transparent", width = 17, height =4, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)

grad_colors <- c(
  "#F7B27A",  "#80D5D9", "#2166AC", "#B2182B"
)

p<-DimPlot(MF,pt.size = 0.1,reduction = "umap",group.by = "stage")+scale_color_manual(values = grad_colors)+
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust =0.03))+
  theme(
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 14,face = "bold")
  )
ggsave("a.tiff", plot = p, bg = "transparent", width = 7.35, height =5, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)

#Cell proportion analysis
ratio<-table(MF$cellType,MF$stage)
ratio<-as.data.frame(ratio)#Calculate proportions
write.csv(ratio,"ratio.csv")
ratio1<-read.csv("ratio.csv",header = T)

colors<- c("#709AE1FF" ,"#FED439FF" ,"#FD7446FF","#8A9197FF","#D2AF81FF"  )
library(ggalluvial)
ratio1$Var2 <- factor(ratio1$Var2, levels = c("Normal state", "Inflammatory phase", "Fibrotic phase","Recovery phase"))
p<-ggplot(ratio1,aes(x=Var2,
                     stratum=Var1,
                     alluvium=Var1,
                     y=Freq,
                     fill = Var1,
                     label=Var1))+
  scale_y_continuous(
    label = scales::percent_format(),  
    expand = c(0.01, 0)  
  ) +
  scale_fill_manual(values = as.character(colors),breaks = names(colors)) +
  geom_stratum(width = 0.3, colour = "grey40")+
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
            size =5, color = "black")+
  scale_x_discrete(position = "top")
ggsave("cycle.tiff", plot = p, bg = "transparent", width = 12, height =6, dpi = 300, device = "tiff", compression = "none",limitsize = FALSE)

#Cells were extracted according to their assigned states
time_cell1=MF[,MF@meta.data$stage %in% c("Normal state")]
time_cell2=MF[,MF@meta.data$stage %in% c("Inflammatory phase")]
time_cell3=MF[,MF@meta.data$stage %in% c("Fibrotic phase")]
time_cell4=MF[,MF@meta.data$stage %in% c("Recovery phase")]

p<-DimPlot(time_cellx,reduction = "umap",pt.size = 0.1)+scale_color_simpsons()+
  theme(
    axis.line = element_blank(),   
    axis.ticks = element_blank(),   
    axis.text = element_blank(),     
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  theme(legend.position = "none")
ggsave("x.tiff", plot = p, bg = "transparent", width = 5, height =5, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)

#Pseudotime analysis
library(monocle3)
MF.cds <- as.cell_data_set(MF)
MF.cds <- cluster_cells(cds = MF.cds, reduction_method = "UMAP")
MF.cds <- learn_graph(MF.cds)
MF.cds=order_cells(MF.cds,)
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
colors<- c("#FED439FF","#709AE1FF","#8A9197FF", "#D2AF81FF", "#FD7446FF")
plot_genes_in_pseudotime(MF.cds[tgen_sig,],color_cells_by="cellType",min_expr=0.5,ncol = 2)+scale_color_manual(values = as.character(colors),breaks = names(colors))

#Construct a pseudotime heatmap
##Please note that this step was performed using Seurat v4, and thus requires prior downgrading of the Seurat package to version 4
genes<-row.names(subset(tgen,morans_I>0.25))
pt.matrix <- exprs(MF.cds) [match(genes, rownames(rowData(MF.cds))), order(pseudotime(MF.cds))]
pt.matrix <- t(apply(pt.matrix, 1, function(x){smooth.spline(x, df=3)$y}))
pt.matrix <- t(apply(pt.matrix, 1, function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
mycol = rev(RColorBrewer::brewer.pal(11, "Spectral"))
library(ComplexHeatmap)
p = ComplexHeatmap::Heatmap(
  pt.matrix, name = "z-score", show_row_names = T, show_column_names = F,
  col = circlize::colorRamp2(seq(from=-2, to=2, length=11), mycol),
  row_names_gp = gpar(fontsize = 12), row_title_rot = 0, km = 4,
  cluster_rows = TRUE, cluster_row_slices = FALSE, cluster_columns = FALSE, use_raster = F
)















