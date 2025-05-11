#This analysis was based on publicly available data from GSE136831. 
#The analysis incorporated 60 samples from 32 patients with idiopathic pulmonary fibrosis (IPF) as well as from 28 healthy donors.
setwd("your dir")
rm(list = ls())
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")

library(BiocManager)
BiocManager::install("multtest")

library(multtest)
BiocManager::install("Seurat")
BiocManager::install("rtracklayer")
library(Seurat)
library(mindr)
library(tidydr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(magrittr)
library(ggplot2)
library(ggsci)
library(Matrix)
#Import data

GSM4058900<-readRDS("GSM4058900_133C-a.dgecounts.rds")

GSM4058900DATA<-GSM4058900$umicount$exon$all   #Use umicounts
a<-list.files()
GSEDATA<-GSM4058900DATA
for (i in 2:83){
  new.data = readRDS(a[i])
  new.data.filter<-new.data$umicount$exon$all
  GSEDATA<-append(GSEDATA,new.data.filter)
  names(GSEDATA)[[i]]<-strsplit(a[i], ".",fixed=T)[[1]][1]
  print(label_percent()(i/83))
}
names(GSEDATA)[[1]]<-strsplit(a[1], ".",fixed=T)[[1]][1]



#Replace gene names(From ensembl to gene symbol )

GENID<-read.table(file="GSE136831_AllCells.GeneIDs.txt",header=F)
n=1
for(b in 1:83){
  
  
  n<-length(GSEDATA[[b]]@Dimnames[[1]])
  
  for(c in 1:n){
    GSEDATA[[b]]@Dimnames[[1]][c]<-GENID[grep(GSEDATA[[b]]@Dimnames[[1]][c],GENID[,1]),2]
    
    
  }
  print(label_percent()(b/83))
}#This process may take some time. Please be patient.

#Create seurat objects
library(scales)

Seurat_obj<-CreateSeuratObject(counts = GSEDATA[[1]],project = names(GSEDATA)[1],min.cells=3,min.features=200) 
for (x in 2:83) {
  Seurat_obj<- append(Seurat_obj,CreateSeuratObject(counts = GSEDATA[[names(GSEDATA)[x]]],project =names(GSEDATA)[x],min.cells=3,min.features=200))
  names(Seurat_obj)[[x]]<-names(GSEDATA)[x]
  print(label_percent()(x/83))
  
}
names(Seurat_obj)[[1]]<-names(GSEDATA)[1]

a=1
for (a in 1:83) {
  Seurat_obj[[a]][["percent.mt"]]<-PercentageFeatureSet(Seurat_obj[[a]], pattern="^MT-")
  
  print(label_percent()(a/83))
  
}
#Quality control screening

filter1<-function(a){
  subset(Seurat_obj[[a]],subset = nFeature_RNA>200 & nFeature_RNA <7500 & percent.mt <15)
}
Seurat_obj_filted<- lapply(1:83,filter1)
#Consolidate the data into a Seurat object
Seurat_obj_filted<-merge(x=Seurat_obj_filted[[1]],y=Seurat_obj_filted[-1])


#Data standardization

Seurat_obj_filted<-NormalizeData(Seurat_obj_filted) %>%FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%  ScaleData() %>%RunPCA(npcs = 30,verbose = T)


#Grouping according to metadata

library(stringr)


group<-as.numeric(str_sub(Normalized_harmony@meta.data$orig.ident,8,10))
group<-ifelse(group<928,"Control",ifelse(group<946,"COPD","IPF"))
Normalized_harmony$group<-group
Normalized_harmony$group

GSE159354@meta.data$sample<-stringr::str_sub(GSE159354$orig.ident,12,20)
Normalized_harmony$sample<-sample
Normalized_harmony$group

#HarmonyIntegration
final_harmony<-IntegrateLayers(Seurat_obj_filted,method = HarmonyIntegration,orig.reduction = "pca", new.reduction="harmony",verbose=F)
final_harmony[["RNA"]]<-JoinLayers(final_harmony[["RNA"]])

#Cell cycle proportion score
g2m_genes<-cc.genes$g2m.genes
g2m_genes<-CaseMatch(search = g2m_genes,match = rownames(final_harmony))
s_genes<-cc.genes$s.genes
s_genes<-CaseMatch(search = s_genes,match = rownames(final_harmony))
final_harmony<-CellCycleScoring(object = final_harmony,g2m.features = g2m_genes,s.features = s_genes,set.ident = T)
Normalized_harmonyOnelayer@meta.data%>%ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase), size=0.5)+theme_minimal()


#Umap plot
final_harmony<-FindNeighbors(final_harmony,reduction = "harmony",dims = 1:30)
final_harmony<-FindClusters(final_harmony,resolution = seq(from=0.1, to=1.0,by=0.1))
final_harmony<- RunUMAP(final_harmony, dims = 1:30, reduction = "harmony", reduction.name = "umap")
Idents(final_harmony)<-"RNA_snn_res.0.2"




p<-DimPlot(final_harmony,reduction = "umap",group.by = "RNA_snn_res.1")+
  theme_dr(xlength = 0.5, 
           ylength = 0.5,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust =0.03))

ggsave("umap.tiff", plot = p, bg = "transparent", width = , height =, dpi = 600, device = "tiff", compression = "none",limitsize = FALSE)

#find marker

markers<-FindAllMarkers(object = final_harmony,test.use = "wilcox",only.pos = T,logfc.threshold = 0.25)
write.csv(markers,"marker.csv")









