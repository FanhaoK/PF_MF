#This analysis was based on publicly available data from GSE207851. 
#GSE207851 included 30 BLM-treated mice, which underwent single-cell sequencing at different time points.
library(stringr)
library(scales)
library(Seurat)
library(ggplot2)
library(clustree)
library(ggsci)
library(tidydr)
#Improt data
Day0_M2<-read.table("Day0_M2.txt",header = T)
rownames(Day0_M2)<-Day0_M2[,1]
Day0_M2<-Day0_M2[,-1]
Day0_M2<-CreateSeuratObject(counts = Day0_M2,project = 'Day0_M2',min.cells = 3, min.features = 200)

list.files()
a<-list.files()
GSEBLM_TIME<-Day0_M2
for (i in 2:30){
  new.data = read.table(a[i],header = T)
  rownames(new.data)<-new.data[,1]
  new.data.filter<-new.data[,-1]
  GSEBLM_TIME<-append(GSEBLM_TIME,CreateSeuratObject(counts = new.data.filter,project = strsplit(a[i], ".",fixed=T)[[1]][1],min.cells = 3, min.features = 200))
  names(GSEBLM_TIME)[[i]]<-strsplit(a[i], ".",fixed=T)[[1]][1]
  print(label_percent()(i/30))
}
names(GSEBLM_TIME)[[1]]<-strsplit(a[2], ".",fixed=T)[[1]][1]
#Quality control
a=1
for (a in 1:30) {
  GSEBLM_TIME[[a]][["percent.mt"]]<-PercentageFeatureSet(GSEBLM_TIME[[a]], pattern="^mt-")
  
  print(label_percent()(a/83))
  
}
GSEBLM_TIME<-merge(x=GSEBLM_TIME[[1]],y=GSEBLM_TIME[-1])
GSEBLM_TIME<-NormalizeData(GSEBLM_TIME) %>%FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%  ScaleData() %>%RunPCA(npcs = 30,verbose = T)
group<-str_sub(GSEBLM_TIME@meta.data$orig.ident,1,5)
GSEBLM_TIME$group<-group
GSEBLM_TIME<-IntegrateLayers(GSEBLM_TIME,method = HarmonyIntegration,orig.reduction = "pca", new.reduction="harmony",verbose=F)
GSEBLM_TIME[["RNA"]]<-JoinLayers(GSEBLM_TIME[["RNA"]])
HB.genes.mouse <- c("Hba-a1", "Hba-a2", "Hbb-bs")
HB_m <- match(HB.genes.mouse, rownames(GSEBLM_TIME@assays$RNA)) 
HB.genes.mouse <- rownames(GSEBLM_TIME@assays$RNA)[HB_m] 
HB.genes.mouse <- HB.genes[!is.na(HB.genes.mouse)] 
GSEBLM_TIME[["percent.HB"]]<-PercentageFeatureSet(GSEBLM_TIME, features=HB.genes.mouse) 
GSEBLM_TIME<- subset(GSEBLM_TIME, subset = nFeature_RNA > 500& nFeature_RNA < 2500 & percent.mt < 5 & percent.HB < 3 )
#Cell cycle proportion score
g2m_genes<-cc.genes$g2m.genes
g2m_genes<-CaseMatch(search = g2m_genes,match = rownames(GSEBLM_TIME))
s_genes<-cc.genes$s.genes
s_genes<-CaseMatch(search = s_genes,match = rownames(GSEBLM_TIME))
GSEBLM_TIME<-CellCycleScoring(object = GSEBLM_TIME,g2m.features = g2m_genes,s.features = s_genes,set.ident = T)
GSEBLM_TIME<- RunPCA(GSEBLM_TIME,npcs = 15) 
GSEBLM_TIME<-FindNeighbors(GSEBLM_TIME,reduction = "harmony",dims = 1:15)
GSEBLM_TIME<-FindClusters(GSEBLM_TIME,resolution = 0.1)
GSEBLM_TIME<- RunUMAP(GSEBLM_TIME, dims = 1:15, reduction = "harmony", reduction.name = "umap")
clustree(GSEBLM_TIME)
GSEBLM_TIME$RNA_snn_res.0.1
Idents(GSEBLM_TIME)<-"RNA_snn_res.0.1"
GSEBLM_TIME<- RunTSNE(GSEBLM_TIME, dims = 1:10, reduction = "harmony", reduction.name = "tsne")
#Findmarker
markers<-FindAllMarkers(object = GSEBLM_TIME,test.use = "wilcox",only.pos = T,logfc.threshold = 0.25)
write.csv(markers,"marker.csv")
#Annotation of cells
GSEBLM_TIME$cellType<-plyr::mapvalues(GSEBLM_TIME$RNA_snn_res.0.1,from = 0:7, to=c("Fibroblasts",
                                                                                   "Myofibroblasts",
                                                                                   "Mesothelial cells",
                                                                                   "Epithelial cells",
                                                                                   "Lymphatic endothelial cells",
                                                                                   "Mesenchymal cells",
                                                                                   "Smooth muscle cells",
                                                                                   "Vascular endothelial cells"
))

new.cluster.ids <- c("Fibroblasts",
                     "Myofibroblasts",
                     "Mesothelial cells",
                     "Epithelial cells",
                     "Lymphatic endothelial cells",
                     "Mesenchymal cells",
                     "Smooth muscle cells",
                     "Vascular endothelial cells"
)
names(new.cluster.ids) <- levels(GSEBLM_TIME)
GSEBLM_TIME<- RenameIdents(GSEBLM_TIME, new.cluster.ids)
#Plot
umap_postion<-GSEBLM_TIME@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% cbind(cellType = GSEBLM_TIME@meta.data$cellType)
celltypepos <- umap_postion %>%group_by(cellType)%>%summarise(umap_a = median(umap_1),umap_b= median(umap_2))
p<-DimPlot(GSEBLM_TIME,pt.size = 0.05)+scale_color_simpsons()+
theme_dr(xlength = 0.2, 
         ylength = 0.2,
         arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust =0.03))+
  theme(legend.text = element_text(size = 15,face = "bold"))+
  NoLegend() +
  geom_label_repel(aes(x = umap_a,y = umap_b,label = cellType,color = cellType), 
                   fontface = "bold",
                   data = celltypepos,
                   box.padding = 0.5,size=5)
p<-FeaturePlot(GSEBLM_TIME,features = c("Mfap5"),cols=c("lightgrey","#197EC0"),pt.size = 0.05)+
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust =0.03))+
  theme(legend.text = element_text(size = 15,face = "bold"))
ggsave("UMAP.tiff", plot = p, bg = "transparent", width =8.3, height =5, dpi = 300, device = "tiff", compression = "none",limitsize = FALSE)

#Isolate the fibroblast population
GSEBLM_TIME_FIBRO= GSEBLM_TIME[,GSEBLM_TIME@meta.data$RNA_snn_res.0.1 %in% c(0,1,6)]
GSEBLM_TIME_FIBRO<-NormalizeData(GSEBLM_TIME_FIBRO) %>%FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%  ScaleData() %>%RunPCA(npcs = 30,verbose = T)
GSEBLM_TIME_FIBRO<- RunPCA(GSEBLM_TIME_FIBRO,npcs = 15) 
GSEBLM_TIME_FIBRO<-FindNeighbors(GSEBLM_TIME_FIBRO,reduction = "harmony",dims = 1:15)
GSEBLM_TIME_FIBRO<-FindClusters(GSEBLM_TIME_FIBRO,resolution = 0.4)
GSEBLM_TIME_FIBRO<- RunUMAP(GSEBLM_TIME_FIBRO, dims = 1:15, reduction = "harmony", reduction.name = "umap")


