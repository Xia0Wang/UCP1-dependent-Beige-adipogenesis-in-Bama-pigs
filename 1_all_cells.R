
options(stringsAsFactors = F)
library(Seurat)
packageVersion("Seurat")
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(stringr) #GO
setwd("E:/BIIF/liutianxia/Adipocyte_scRNA_seq_8/3_结果/0")
getwd()

#list=ls()
#rm(list=ls()[c(2:28)])
#rm(list=ls()[1:10])
barplot(1:36,col = DiscretePalette(36, palette = "polychrome"),names.arg = c(1:36))
barplot(1:36,col = DiscretePalette(36, palette = "polychrome")[c(4:8,11:16,18:36)],names.arg = c(1:36))
#颜色
col_1 <- DiscretePalette(36, palette = "polychrome")[c(2:8,11:16,18:36)]
col_3 <- DiscretePalette(36, palette = "polychrome")[c(4,7:8,27,11:15,32,18,20,28:36)]

xincol <- c("#FAD38E","#EAA4FC","#DF90BF","#6ABFD4","#BB806D",
            "#4E9134","#F6A8A7","#A188FE","#94ABA9","#CFEF89")


#########step1######################
#read
yangpin = c("6609_JJ","6609_SZ","7103_JJ","7103_SZ")
zhang_List = lapply(yangpin,function(pro){
  folder=file.path("E:/BIIF/liutianxia/Adipocyte_scRNA_seq_7/1_juzhen/",pro)
  CreateSeuratObject(counts = Read10X(folder),
                     project = pro )
})

#########step2######################
#merge
duixiang_zhifang <- merge(x = zhang_List[[1]],
                        y = c(zhang_List[[2]],zhang_List[[3]],zhang_List[[4]]),
                        add.cell.ids = c("6609_JJ","6609_SZ","7103_JJ","7103_SZ"), 
                        project = "zhang")

#########step2######################
#mito genes
duixiang_zhifang[["percent.mt"]] <- PercentageFeatureSet(duixiang_zhifang,
                                                       features = c("ND1","ND2","COX1","COX2",
                                                                    "ATP8","ATP6","COX3","ND3",
                                                                    "ND4L","ND4","ND5","ND6"))
#RPgenes
RPgenes <- read.csv("RP核糖体蛋白.csv")
duixiang_zhifang[["percent.RP"]]<- PercentageFeatureSet(duixiang_zhifang,
                                                        features = RPgenes[c(1:84),]$x)
duixiang_zhifang[["percent.RN"]]<- PercentageFeatureSet(duixiang_zhifang,
                                                        features = RPgenes[c(85:86),]$x)

#plot
plot1 <- VlnPlot(duixiang_zhifang, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP","percent.RN"), ncol = 3)
ggsave(plot1, file='plot1-过滤前质控.pdf', width=15, height=8)
plot2 <- VlnPlot(duixiang_zhifang, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP","percent.RN"),pt.size = 0, ncol = 3)
ggsave(plot2, file='plot2-过滤前质控-没有点.pdf', width=15, height=8)

#remove
dim(duixiang_zhifang)
a <- duixiang_zhifang@assays[["RNA"]]@data@Dimnames[[1]]
c <- match(setdiff(a,RPgenes[c(85:86),]$x),a)
duixiang_zhifang1 <- duixiang_zhifang[c,]
dim(duixiang_zhifang1)

#Filtering
duixiang_zhifang2 <- subset(duixiang_zhifang1, subset = nFeature_RNA > 500 
                            & nFeature_RNA < 6000 
                            & nCount_RNA > 500 
                            & percent.mt < 10
                            & percent.RP < 20
                            & percent.RN < 5)
table(duixiang_zhifang2$orig.ident)

#plot
plot3 <- VlnPlot(duixiang_zhifang2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP"), ncol = 2)
ggsave(plot3, file='Plot3-过滤后质控.pdf', width=10, height=8)
plot4 <- VlnPlot(duixiang_zhifang2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP"),pt.size = 0, ncol = 2)
ggsave(plot4, file='Plot4-过滤后质控-没有点.pdf', width=10, height=8)

#########step4######################
g.list <- SplitObject(duixiang_zhifang2, split.by = "orig.ident")

#########step5######################
#
for(i in 1:length(g.list)) {
  g.list[[i]] <- NormalizeData(g.list[[i]],
                               normalization.method =  "LogNormalize",
                               scale.factor = 10000)
  g.list[[i]] <- FindVariableFeatures(g.list[[i]], 
                                      selection.method = "vst", 
                                      nfeatures = 2000)
}

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = g.list) 
AllBatch.anchors <- FindIntegrationAnchors(object.list = g.list,
                                           anchor.features = features,
                                           dims = 1:20,
                                           k.anchor = 5,
                                           k.filter = 200,
                                           k.score = 30,)
zhifang <- IntegrateData(anchorset = AllBatch.anchors, dims = 1:20)

DefaultAssay(zhifang) <- "integrated"

all.genes <- rownames(zhifang)
zhifang <- ScaleData(zhifang, features = all.genes)
zhifang <- RunPCA(object = zhifang, pc.genes = VariableFeatures(zhifang))
#rm(AllBatch.anchors)
#rm(zhang_List)
#rm(duixiang_zhifang1)
#rm(g.list)
#ElbowPlot
plot8 <- ElbowPlot(zhifang)
ElbowPlot(zhifang,ndims = 50, reduction = "pca") 
table(zhifang@meta.data$orig.ident)

#########step6######################
set.seed(123)
zhifang_umap <- zhifang
#umap
zhifang_umap <- RunUMAP(zhifang_umap, dims = 1:20)
#cluster
zhifang_umap <- FindNeighbors(zhifang_umap, dims = 1:20)
zhifang_umap <- FindClusters(zhifang_umap, resolution = 0.8)
#plot
plot13 <- DimPlot(zhifang_umap,reduction = "umap",label=F, group.by = "orig.ident",pt.size = 0.0001,
                  cols = DiscretePalette(26, palette = "alphabet2")[c(19:26)])
ggsave(plot13, file='Plot13-umap-1.pdf', width=6, height=5)
plot14 <- DimPlot(zhifang_umap,reduction = "umap",label=F, group.by = "orig.ident", split.by ='orig.ident',
                  cols = DiscretePalette(26, palette = "alphabet2")[c(19:26)])
ggsave(plot14, file='Plot14-umap-2.pdf', width=18, height=5)
plot15 <- DimPlot(zhifang_umap,reduction = "umap",label=T,
                  cols = col_1)
ggsave(plot15, file='Plot15-umap-3.pdf', width=6.2, height=5)
plot16 <- DimPlot(zhifang_umap,reduction = "umap",label=T,split.by ='orig.ident',
                  cols = col_1)
ggsave(plot16, file='Plot16-umap-4.pdf', width=18, height=5)

#Assay
DefaultAssay(zhifang_umap) <- "RNA"
DefaultAssay(zhifang_umap)

#plot
plot5 <- VlnPlot(zhifang_umap, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP"),
                 group.by = "seurat_clusters",
                 split.by = "seurat_clusters" ,ncol = 2)
ggsave(plot5, file='Plot5-过滤后质控.pdf', width=20, height=8)
plot6 <- VlnPlot(zhifang_umap, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP"),
                 group.by  = "seurat_clusters",
                 split.by = "seurat_clusters",pt.size = 0, ncol = 2)
ggsave(plot6, file='Plot6-过滤后质控-没有点.pdf', width=20, height=8)

#########step7######################
markers.to.plot2 <- c("PLIN1","ADIPOQ",#脂肪
                      "PDGFRA","ABCA6",#aspc
                      "CYYR1","VWF",#内皮
                      "PROX1","RELN",#淋巴内皮
                      "GUCY1A2","TRPC6",#周皮
                      "XKR4","SCN7A",#施旺细胞
                      "MYOCD","MYH11",#血管平滑肌
                      "NEB","TRDN",#骨骼肌
                      "CD163","MRC1",#巨噬细胞
                      "SKAP1","ITGA4",#T细胞
                      "TBXAS1","LYN",#DC
                      "MKI67","TOP2A",#KI67
                      "EEF1A1","TPT1")
#总体
dotplot <- DotPlot(zhifang_umap, 
                   features = markers.to.plot2, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪.pdf', width=10.5, height=7)

markers.to.plot3 <- c("PLIN1","ADIPOQ",#脂肪
                      "PDGFRA","ABCA6",#aspc
                      "ITGB4","BICD1", #干细胞
                      "CYYR1","VWF",#内皮
                      "PROX1","RELN",#淋巴内皮
                      "GPM6A","WWC1",#间皮
                      "GUCY1A2","TRPC6",#周皮
                      "MYOCD","MYH11",#血管平滑肌
                      "XKR4","SCN7A",#施旺细胞
                      "CD163","MRC1",#巨噬细胞
                      "SKAP1","ITGA4",#T细胞
                      "TBXAS1","LYN",#DC
                      "MKI67","TOP2A")#KI67
#总体
dotplot <- DotPlot(zhifang_umap, 
                   features = markers.to.plot3, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪_2.pdf', width=10.5, height=7)

#########step8######################
#为了判断double,还算一下每簇细胞的平均feature和counts数。
mean_Count <- aggregate(zhifang_umap@meta.data, nCount_RNA ~ seurat_clusters,mean)
mean_Feature <- aggregate(zhifang_umap@meta.data, nFeature_RNA ~ seurat_clusters,mean)

#脂肪： 1 3 5 15 16 24 27 29
#FAP: 0 2 5 7 8 11 22 23
#内皮：9 10 12 16 20 23 30
#周皮：13 21 27 30
#淋巴内皮：17 29
#施旺细胞：26
#平滑肌： 19 21
#SKC: 28 去掉
#巨噬细胞：6 22 24
#T：14
#DC：18
#4是低质量的那组,count数和feature数都低，过渡态Transitional
#通过上述，判断出double包括  5 16 21-24 27-30

guolv_1 <- subset(zhifang_umap, seurat_clusters %in% c(5,16,21:24,27:30), invert = T)
DefaultAssay(guolv_1) <- "integrated"
#设置随机数
#可选项，计算umap的分群
guolv_1 <- RunUMAP(guolv_1, dims = 1:20)
#（3）分群。
guolv_1 <- FindNeighbors(guolv_1, dims = 1:20)
guolv_1 <- FindClusters(guolv_1, resolution = 0.7)

#画图
plot13 <- DimPlot(guolv_1,reduction = "umap",label=F, group.by = "orig.ident",
                  cols = DiscretePalette(26, palette = "alphabet2")[c(19:26)])
ggsave(plot13, file='Plot17-umap-1-去除double之后.pdf', width=6, height=5)
plot14 <- DimPlot(guolv_1,reduction = "umap",label=F, group.by = "orig.ident", split.by ='orig.ident',
                  cols = DiscretePalette(26, palette = "alphabet2")[c(19:26)])
ggsave(plot14, file='Plot18-umap-2.pdf', width=18, height=5)
plot15 <- DimPlot(guolv_1,reduction = "umap",label=T,
                  cols = col_1)
ggsave(plot15, file='Plot19-umap-3.pdf', width=6.2, height=5)
plot16 <- DimPlot(guolv_1,reduction = "umap",label=T,split.by ='orig.ident',
                  cols = col_1)
ggsave(plot16, file='Plot20-umap-4.pdf', width=18, height=5)

mean_Count22 <- aggregate(guolv_1@meta.data, nCount_RNA ~ seurat_clusters,mean)
mean_Feature22 <- aggregate(guolv_1@meta.data, nFeature_RNA ~ seurat_clusters,mean)
#统计分群结果
table(guolv_1@meta.data$seurat_clusters)
# 在两个样本里各个cluster的分布情况
table(guolv_1@meta.data$seurat_clusters,guolv_1@meta.data$orig.ident)
#修改默认的Assay
DefaultAssay(guolv_1) <- "RNA"
DefaultAssay(guolv_1)

#########step9######################
markers.to.plot2 <- c("PLIN1","ADIPOQ",#脂肪
                      "PDGFRA","ABCA6",#aspc
                      "CYYR1","VWF",#内皮
                      "PROX1","RELN",#淋巴内皮
                      "GUCY1A2","TRPC6",#周皮
                      "XKR4","SCN7A",#施旺细胞
                      "MYOCD","MYH11",#血管平滑肌
                      "NEB","TRDN",#骨骼肌
                      "CD163","MRC1",#巨噬细胞
                      "SKAP1","ITGA4",#T细胞
                      "TBXAS1","LYN",#DC
                      "MKI67","TOP2A",#KI67
                      "EEF1A1","TPT1")
#总体
dotplot <- DotPlot(guolv_1, 
                   features = markers.to.plot2, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪_3_去除double之后.pdf', width=10.5, height=9)

#肌肉系列
markers.to.plot <- c("PDGFRA","DCN",#FAP
                     "PAX7","NCAM1",#肌卫星
                     "DIAPH3","MKI67",#ki67
                     "OBSCN","MYPN",#肌肉
                     "ADIPOQ", "PPARG",#脂肪
                     "PECAM1", "CYYR1",#内皮
                     "GUCY1A2","FHL5",#周皮
                     "PROX1","TBX1",#淋巴内皮
                     "COL28A1","NEGR1",#成纤维
                     "MBP","CDH19",#施旺细胞
                     "CD163","MRC1",#巨噬细胞
                     "AFF3","ITGB4",#B细胞
                     "SKAP1","PTPRC")#T细胞
#总体
dotplot <- DotPlot(guolv_1, 
                   features = markers.to.plot, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_肌肉_去除double之后.pdf', width=10.5, height=7)

#########step10######################
guolv_2 <- subset(guolv_1, seurat_clusters %in% c(21), invert = T)
DefaultAssay(guolv_2) <- "integrated"
#设置随机数
#可选项，计算umap的分群
guolv_2 <- RunUMAP(guolv_2, dims = 1:20)
#（3）分群。
guolv_2 <- FindNeighbors(guolv_2, dims = 1:20)
guolv_2 <- FindClusters(guolv_2, resolution = 0.7)

#画图
plot13 <- DimPlot(guolv_2,reduction = "umap",label=F, group.by = "orig.ident",
                  cols = DiscretePalette(26, palette = "alphabet2")[c(19:26)])
ggsave(plot13, file='Plot21-umap-1-再去除double之后.pdf', width=6, height=5)
plot14 <- DimPlot(guolv_2,reduction = "umap",label=F, group.by = "orig.ident", split.by ='orig.ident',
                  cols = DiscretePalette(26, palette = "alphabet2")[c(19:26)])
ggsave(plot14, file='Plot22-umap-2.pdf', width=18, height=5)
plot15 <- DimPlot(guolv_2,reduction = "umap",label=T,
                  cols = col_1)
ggsave(plot15, file='Plot23-umap-3.pdf', width=6.2, height=5)
plot16 <- DimPlot(guolv_2,reduction = "umap",label=T,split.by ='orig.ident',
                  cols = col_1)
ggsave(plot16, file='Plot24-umap-4.pdf', width=18, height=5)

mean_Count22 <- aggregate(guolv_2@meta.data, nCount_RNA ~ seurat_clusters,mean)
mean_Feature22 <- aggregate(guolv_2@meta.data, nFeature_RNA ~ seurat_clusters,mean)
#统计分群结果
table(guolv_2@meta.data$seurat_clusters)
# 在两个样本里各个cluster的分布情况
table(guolv_2@meta.data$seurat_clusters,guolv_2@meta.data$orig.ident)
#修改默认的Assay
DefaultAssay(guolv_2) <- "RNA"
DefaultAssay(guolv_2)

#########step11######################
markers.to.plot2 <- c("PLIN1","ADIPOQ",#脂肪
                      "PDGFRA","ABCA6",#aspc
                      "CYYR1","VWF",#内皮
                      "PROX1","RELN",#淋巴内皮
                      "GUCY1A2","TRPC6",#周皮
                      "XKR4","SCN7A",#施旺细胞
                      "MYOCD","MYH11",#血管平滑肌
                      "NEB","TRDN",#骨骼肌
                      "CD163","MRC1",#巨噬细胞
                      "SKAP1","ITGA4",#T细胞
                      "TBXAS1","LYN",#DC
                      "MKI67","TOP2A",#KI67
                      "EEF1A1","TPT1")
#总体
dotplot <- DotPlot(guolv_2, 
                   features = markers.to.plot2, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪_4_再去除double之后.pdf', width=10.5, height=9)

#########step12######################
#FindAllMarkers
ALL.markers <- FindAllMarkers(guolv_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ALL.markers,file="ALL.markers.csv")
ALL.markers <- read.csv(file="ALL.markers.csv")

#top
top10 <- ALL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="ALL.marker.top10.csv")
top50 <- ALL.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50,file="ALL.marker.top50.csv")

#########step13######################
#0 3 4 5 7 9前脂肪
#1 2 14 脂肪
#8 11 12 15 内皮
#16 淋巴内皮
#6 巨噬细胞
#13 19 T细胞
#10 周皮
#18 DCs
#17 血管平滑肌
#20 施旺细胞

new.cluster.ids <- c("ASPC", "Adipocyte","Adipocyte", "ASPC",  "ASPC", "ASPC",
                     "Macrophage", "ASPC","Endothelial","ASPC","Pericyte",
                     "Endothelial","Endothelial","T cell","Adipocyte","Endothelial",
                     "LEC","Smooth muscle","Dendritic cell","T cell","Schwann cell")
names(new.cluster.ids) <- levels(guolv_2)
levels(guolv_2)
guolv_2 <- RenameIdents(guolv_2, new.cluster.ids)
levels(guolv_2)
guolv_2@active.ident <- factor(guolv_2@active.ident, 
                               levels = c("Adipocyte","ASPC","Endothelial","Pericyte","Smooth muscle",
                                          "LEC","Schwann cell",
                                          "Macrophage","Dendritic cell","T cell"))
levels(guolv_2)
table(guolv_2@active.ident)

#加入一列细胞
guolv_2@meta.data$xibao <- guolv_2@active.ident
#修改顺序
#guolv_2@meta.data$orig.ident <- factor(guolv_2@meta.data$orig.ident,
#                                       levels = c("29907_FGG","29907_SZ","6609_JJ",
#                                                  "6609_SZ","7103_JJ","7103_SZ"))
guolv_2@meta.data$orig.ident <- factor(guolv_2@meta.data$orig.ident,
                                       levels = c("6609_JJ","7103_JJ",
                                                  "6609_SZ","7103_SZ"))

#不带legend的图
plot33 <- DimPlot(guolv_2, reduction = "umap", label = T, pt.size = 0.5,
                  cols = col_3) + NoLegend()
ggsave(plot33, file='plot33-umap-25-细胞命名.pdf', width=5.5, height=5)
plot34 <- DimPlot(guolv_2, reduction = "umap", pt.size = 0.5,
                  label = T, split.by ='orig.ident',
                  cols = col_3) + NoLegend()
ggsave(plot34, file='plot34-umap-26.pdf', width=16, height=5)
plot35 <- DimPlot(guolv_2, reduction = "umap", label = F,pt.size = 0.5,
                  cols = col_3)
ggsave(plot35, file='plot35-umap-27.pdf', width=6.8, height=5)
plot36 <- DimPlot(guolv_2, reduction = "umap", pt.size = 0.5,
                  label = F, split.by ='orig.ident',
                  cols = col_3)
ggsave(plot36, file='plot36-umap-28.pdf', width=17.2, height=5)

#新颜色
plot33 <- DimPlot(guolv_2, reduction = "umap", label = T, pt.size = 0.5,
                  cols = xincol) + NoLegend()
ggsave(plot33, file='plot1033-umap-25-细胞命名.pdf', width=5.5, height=5)
plot34 <- DimPlot(guolv_2, reduction = "umap", pt.size = 0.5,
                  label = T, split.by ='orig.ident',
                  cols = xincol) + NoLegend()
ggsave(plot34, file='plot1034-umap-26.pdf', width=16, height=5)
plot35 <- DimPlot(guolv_2, reduction = "umap", label = F,pt.size = 0.5,
                  cols = xincol)
ggsave(plot35, file='plot1035-umap-27.pdf', width=6.8, height=5)
plot36 <- DimPlot(guolv_2, reduction = "umap", pt.size = 0.5,
                  label = F, split.by ='orig.ident',
                  cols = xincol)
ggsave(plot36, file='plot1036-umap-28.pdf', width=17.2, height=5)

#########step14######################
table(guolv_2@meta.data$orig.ident)
guolv_2@meta.data$weizhi <- substring(guolv_2@meta.data$orig.ident,7,7)
table(guolv_2@meta.data$weizhi)
guolv_2@meta.data$weizhi[which(guolv_2@meta.data$weizhi == "J")] <- "sWAT"
guolv_2@meta.data$weizhi[which(guolv_2@meta.data$weizhi == "Z")] <- "pWAT"
guolv_2@meta.data$weizhi <- factor(guolv_2@meta.data$weizhi,
                                   levels = c("sWAT","pWAT"))

guolv_2@meta.data$chuli <- substring(guolv_2@meta.data$orig.ident,1,1)
table(guolv_2@meta.data$chuli)
guolv_2@meta.data$chuli[which(guolv_2@meta.data$chuli == "6")] <- "WT"
guolv_2@meta.data$chuli[which(guolv_2@meta.data$chuli == "7")] <- "KI"
guolv_2@meta.data$chuli <- factor(guolv_2@meta.data$chuli,
                                  levels = c("WT","KI"))
#分别是野生常温，野生冷刺激，KI冷刺激
#根据WT和Mu区分2种颜色来画图
plot37 <- DimPlot(guolv_2,reduction = "umap",label=F, split.by ='chuli',pt.size = 0.5,
                  cols = col_3)
ggsave(plot37, file='plot37-umap-29.pdf', width=11, height=5)
#根据2个位置，用2种颜色来画图
plot38 <- DimPlot(guolv_2,reduction = "umap",label=F, split.by ='weizhi',pt.size = 0.5,
                  cols = col_3)
ggsave(plot38, file='Plot38-umap-30.pdf', width=11, height=5)

#########step15######################
#FindAllMarkers
table(guolv_2@meta.data$xibao)
ALL.markers_cells <- FindAllMarkers(guolv_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ALL.markers_cells,file="ALL.markers细胞命名之后.csv")
#ALL.markers_cells <- read.csv(file="ALL.markers_cells.csv")

#top
top10 <- ALL.markers_cells %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="ALL.markers细胞命名之后_top10.csv")
top50 <- ALL.markers_cells %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50,file="ALL.markers细胞命名之后_top50.csv")

#########step16######################
markers.to.plot5 <- c("PLIN1","ADIPOQ",#脂肪
                      "PDGFRA","ABCA6",#fap
                      "CYYR1","VWF",#内皮
                      "GUCY1A2","TRPC6",#周皮
                      "MYOCD","MYH11",#血管平滑肌
                      "PROX1","RELN",#淋巴内皮
                      "XKR4","SCN7A",#施旺细胞
                      "CD163","MRC1",#巨噬细胞
                      "TBXAS1","LYN",#DC
                      "SKAP1","ITGA4"#T细胞
                      )
#总体
dotplot <- DotPlot(guolv_2, 
                   features = markers.to.plot5, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪_5_细胞命名之后.pdf', width=10, height=4.5)

#########step17######################
plot401 <- VlnPlot(guolv_2,
                   features = c("PDGFRA","ABCA6",#fap
                                "PLIN1","ADIPOQ",#脂肪
                                "EEF1A1","TPT1", #Transitional
                                "CYYR1","VWF",#内皮
                                "GUCY1A2","TRPC6",#周皮
                                "MYOCD","MYH11"#血管平滑肌
                   ),
                   cols = col_3,
                   #group.by = "WTMu", 
                   pt.size = 0, 
                   combine = T,
                   same.y.lims = F,
                   ncol = 1)
#plot401 <- CombinePlots(plots = plot400, ncol = 1)
ggsave(plot401, file='plot401-小提琴图1.pdf', width=6, height=29)

plot401 <- VlnPlot(guolv_2,
                   features = c("PROX1","RELN",#淋巴内皮
                                "CD163","MRC1",#巨噬细胞
                                "SKAP1","ITGA4",#T细胞
                                "TBXAS1","LYN",#DC
                                "XKR4","SCN7A",#施旺细胞
                                "NEB","TRDN"#骨骼肌
                   ),
                   cols = col_3,
                   #group.by = "WTMu", 
                   pt.size = 0, 
                   combine = T,
                   same.y.lims = F,
                   ncol = 1)
#plot401 <- CombinePlots(plots = plot400, ncol = 1)
ggsave(plot401, file='plot401-小提琴图2.pdf', width=6, height=29)

#########step18######################
cell.numbers <- table(guolv_2@meta.data$xibao,guolv_2@meta.data$orig.ident)
#横向合并
cell.numbers <- cbind(cell.numbers)
#纵向求和
qiuhe1 <- colSums(cell.numbers)
#纵向合并
cell.numbers <- rbind(cell.numbers,qiuhe1)
##野生和突变组的细胞数目不同，需要先标准化成，各10000个细胞中，有几个这种细胞
for (i in 1:ncol(cell.numbers)){
  cell.numbers[,i] <- cell.numbers[,i]/cell.numbers[nrow(cell.numbers),i]
}
#去掉求和那一列
cell.numbers <- cell.numbers[-nrow(cell.numbers),]
#write.csv(cell.numbers,"plot501-num-1.csv")
#转换数据框
cell.numbers <- as.data.frame(cell.numbers) 
cell.numbers$clusters <- row.names(cell.numbers)
#宽表变长表
cell.numbers <- reshape2::melt(cell.numbers, id.vars = "clusters")
#折线图中横坐标必须为离散型数值和数值型数值，否则会画不出来！！！！！
#cell.numbers$clusters <- as.numeric(cell.numbers$clusters)
cell.numbers$clusters <- factor(cell.numbers$clusters, levels = levels(guolv_2))
#转换为因子
cell.numbers$variable <- factor(cell.numbers$variable)
cell.numbers$value <- as.numeric(cell.numbers$value)

#画柱状图
plot408 <- ggplot(cell.numbers, aes(x = clusters, y = value, fill = variable) ) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "white", size=0, width = 0.7)+
  scale_fill_manual(values = DiscretePalette(36, palette = "polychrome")[c(8,11,15,33,31,32)])+
  theme_classic()+
  #geom_text(aes(label=value),position = position_dodge(0.8),size=1.5,vjust=-0.5)+ #增加柱状图上的数字
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, size = 6), #修改X轴文字
        legend.position = c(0.9,0.7))+ #修改图例位置
  labs(y = "Average fraction of all nuclei")+ #title = "muscle"
  scale_y_continuous(expand = c(0,0))#+
  #scale_y_continuous(expand = c(0,0),limits = c(0,4000)) #修改Y轴跨度
ggsave(plot408, file='plot501-num-1.pdf', width=5, height=4)

#########step19######################
cell.numbers <- table(guolv_2@active.ident)
#横向合并
cell.numbers <- cbind(cell.numbers)
#纵向求和
qiuhe1 <- colSums(cell.numbers)
#纵向合并
cell.numbers <- rbind(cell.numbers,qiuhe1)
##野生和突变组的细胞数目不同，需要先标准化成，各10000个细胞中，有几个这种细胞
for (i in 1:ncol(cell.numbers)){
  cell.numbers[,i] <- cell.numbers[,i]/cell.numbers[nrow(cell.numbers),i]
}
#去掉求和那一列
cell.numbers <- cell.numbers[-nrow(cell.numbers),]
#write.csv(cell.numbers,"plot502-num-总体的大群.csv")
#转换数据框
cell.numbers <- as.data.frame(cell.numbers) 
cell.numbers$clusters <- row.names(cell.numbers)
#宽表变长表
cell.numbers <- reshape2::melt(cell.numbers, id.vars = "clusters")
#折线图中横坐标必须为离散型数值和数值型数值，否则会画不出来！！！！！
#cell.numbers$clusters <- as.numeric(cell.numbers$clusters)
cell.numbers$clusters <- factor(cell.numbers$clusters, levels = levels(guolv_2))
#转换为因子
cell.numbers$variable <- factor(cell.numbers$variable)
cell.numbers$value <- as.numeric(cell.numbers$value)

#画柱状图
plot408 <- ggplot(cell.numbers, aes(x = clusters, y = value, fill = clusters) ) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "white", size=0, width = 0.7)+
  scale_fill_manual(values = col_3)+
  theme_classic()+
  #geom_text(aes(label=value),position = position_dodge(0.8),size=1.5,vjust=-0.5)+ #增加柱状图上的数字
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, size = 6), #修改X轴文字
  )+ #修改图例位置legend.position = c(1,1)
  labs(y = "Average fraction of all nuclei")+ #title = "muscle"
  scale_y_continuous(expand = c(0,0))#+
#scale_y_continuous(expand = c(0,0),limits = c(0,4000)) #修改Y轴跨度
ggsave(plot408, file='plot502-num-总体的大群.pdf', width=6, height=4)

#画柱状图
plot408 <- ggplot(cell.numbers, aes(x = clusters, y = value, fill = clusters) ) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "white", size=0, width = 0.7)+
  scale_fill_manual(values = xincol)+
  theme_classic()+
  #geom_text(aes(label=value),position = position_dodge(0.8),size=1.5,vjust=-0.5)+ #增加柱状图上的数字
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, size = 6), #修改X轴文字
  )+ #修改图例位置legend.position = c(1,1)
  labs(y = "Average fraction of all nuclei")+ #title = "muscle"
  scale_y_continuous(expand = c(0,0))#+
#scale_y_continuous(expand = c(0,0),limits = c(0,4000)) #修改Y轴跨度
ggsave(plot408, file='plot1502-num-总体的大群.pdf', width=6, height=4)

