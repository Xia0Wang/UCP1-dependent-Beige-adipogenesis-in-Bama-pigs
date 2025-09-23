#在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类
options(stringsAsFactors = F)
library(Seurat)
packageVersion("Seurat")
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(stringr) #GO用到
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

###先读取，再合并，再过滤，再拆分，再去批次效应。####
#########步骤1 读取######################
yangpin = c("6609_JJ","6609_SZ","7103_JJ","7103_SZ")
#或者直接读取这个文件夹里面的文件夹的名字。
#yangpin = list.files(path = "1_juzhen/")
#用函数，挨个创建对象。
zhang_List = lapply(yangpin,function(pro){
  folder=file.path("E:/BIIF/liutianxia/Adipocyte_scRNA_seq_7/1_juzhen/",pro)
  CreateSeuratObject(counts = Read10X(folder),
                     project = pro )
})

#########步骤2 合并######################
#读取之后，要合并
duixiang_zhifang <- merge(x = zhang_List[[1]],
                        y = c(zhang_List[[2]],zhang_List[[3]],zhang_List[[4]]),
                        add.cell.ids = c("6609_JJ","6609_SZ","7103_JJ","7103_SZ"), 
                        project = "zhang")

#########步骤3 过滤######################
#统计线粒体的基因。基因名字从jiyin中找来自线粒体的。
duixiang_zhifang[["percent.mt"]] <- PercentageFeatureSet(duixiang_zhifang,
                                                       features = c("ND1","ND2","COX1","COX2",
                                                                    "ATP8","ATP6","COX3","ND3",
                                                                    "ND4L","ND4","ND5","ND6"))
#读取核糖体蛋白相关基因，和核糖体RNA相关基因，这些基因是把RPL和RPS开头的基因都读进来。
RPgenes <- read.csv("RP核糖体蛋白.csv")
duixiang_zhifang[["percent.RP"]]<- PercentageFeatureSet(duixiang_zhifang,
                                                        features = RPgenes[c(1:84),]$x)
duixiang_zhifang[["percent.RN"]]<- PercentageFeatureSet(duixiang_zhifang,
                                                        features = RPgenes[c(85:86),]$x)

#先画
plot1 <- VlnPlot(duixiang_zhifang, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP","percent.RN"), ncol = 3)
ggsave(plot1, file='plot1-过滤前质控.pdf', width=15, height=8)
plot2 <- VlnPlot(duixiang_zhifang, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP","percent.RN"),pt.size = 0, ncol = 3)
ggsave(plot2, file='plot2-过滤前质控-没有点.pdf', width=15, height=8)

#去掉核糖体RNA的基因
#不去除核糖体蛋白
dim(duixiang_zhifang)
a <- duixiang_zhifang@assays[["RNA"]]@data@Dimnames[[1]]
c <- match(setdiff(a,RPgenes[c(85:86),]$x),a)
duixiang_zhifang1 <- duixiang_zhifang[c,]
dim(duixiang_zhifang1)

#按照三个指标过滤细胞
duixiang_zhifang2 <- subset(duixiang_zhifang1, subset = nFeature_RNA > 500 
                            & nFeature_RNA < 6000 
                            & nCount_RNA > 500 
                            & percent.mt < 10
                            & percent.RP < 20
                            & percent.RN < 5)
table(duixiang_zhifang2$orig.ident)

#作图
plot3 <- VlnPlot(duixiang_zhifang2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP"), ncol = 2)
ggsave(plot3, file='Plot3-过滤后质控.pdf', width=10, height=8)
plot4 <- VlnPlot(duixiang_zhifang2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP"),pt.size = 0, ncol = 2)
ggsave(plot4, file='Plot4-过滤后质控-没有点.pdf', width=10, height=8)

#############################步骤3.5.画质控信息##########################
plot101 <- VlnPlot(duixiang_zhifang2,
                   features = "nCount_RNA",
                   group.by = 'orig.ident',
                   split.by = 'orig.ident',
                   pt.size = 0,)
ggsave(plot101, file='Plot11-Mean-UMI.pdf', width=7, height=5)

mean1 <- aggregate(duixiang_zhifang2@meta.data, nCount_RNA ~ orig.ident,mean)
write.csv(mean1,"Plot11-Mean-UMI.csv")

plot101 <- VlnPlot(duixiang_zhifang2,
                   features = "nFeature_RNA",
                   group.by = 'orig.ident',
                   split.by = 'orig.ident',
                   pt.size = 0,)
ggsave(plot101, file='Plot12-Mean-gene.pdf', width=7, height=5)

mean2 <- aggregate(duixiang_zhifang2@meta.data, nFeature_RNA ~ orig.ident,mean)
write.csv(mean2,"Plot12-Mean-gene.csv")

num2 <- as.data.frame(table(duixiang_zhifang2@meta.data$orig.ident))
write.csv(num2,"Plot12-质控后每个样品cellnum.csv")

#########步骤4 拆分######################
g.list <- SplitObject(duixiang_zhifang2, split.by = "orig.ident")

#########步骤5 整合######################
#先分别标准化和找可变基因
for(i in 1:length(g.list)) {
  g.list[[i]] <- NormalizeData(g.list[[i]],
                               normalization.method =  "LogNormalize",
                               scale.factor = 10000)
  g.list[[i]] <- FindVariableFeatures(g.list[[i]], 
                                      selection.method = "vst", 
                                      nfeatures = 2000)
}

# select features that are repeatedly variable across datasets for integration
#默认是2000个
features <- SelectIntegrationFeatures(object.list = g.list) 
AllBatch.anchors <- FindIntegrationAnchors(object.list = g.list,
                                           anchor.features = features,
                                           dims = 1:20,
                                           k.anchor = 5,
                                           k.filter = 200,
                                           k.score = 30,)
zhifang <- IntegrateData(anchorset = AllBatch.anchors, dims = 1:20)
#设置数据
DefaultAssay(zhifang) <- "integrated"
#再标准化
all.genes <- rownames(zhifang)
zhifang <- ScaleData(zhifang, features = all.genes)
zhifang <- RunPCA(object = zhifang, pc.genes = VariableFeatures(zhifang))
#rm(AllBatch.anchors)
#rm(zhang_List)
#rm(duixiang_zhifang1)
#rm(g.list)
#拐点
plot8 <- ElbowPlot(zhifang)
ElbowPlot(zhifang,ndims = 50, reduction = "pca") 
table(zhifang@meta.data$orig.ident)

#################################步骤6.细胞分群############
#设置随机数
set.seed(123)
#提前保存一个
zhifang_umap <- zhifang
#可选项，计算umap的分群
zhifang_umap <- RunUMAP(zhifang_umap, dims = 1:20)
#（3）分群。
zhifang_umap <- FindNeighbors(zhifang_umap, dims = 1:20)
zhifang_umap <- FindClusters(zhifang_umap, resolution = 0.8)
#画图
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
#统计分群结果
table(zhifang_umap@meta.data$seurat_clusters)
# 在两个样本里各个cluster的分布情况
table(zhifang_umap@meta.data$seurat_clusters,zhifang_umap@meta.data$orig.ident)
#修改默认的Assay
DefaultAssay(zhifang_umap) <- "RNA"
DefaultAssay(zhifang_umap)

#作图
plot5 <- VlnPlot(zhifang_umap, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP"),
                 group.by = "seurat_clusters",
                 split.by = "seurat_clusters" ,ncol = 2)
ggsave(plot5, file='Plot5-过滤后质控.pdf', width=20, height=8)
plot6 <- VlnPlot(zhifang_umap, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.RP"),
                 group.by  = "seurat_clusters",
                 split.by = "seurat_clusters",pt.size = 0, ncol = 2)
ggsave(plot6, file='Plot6-过滤后质控-没有点.pdf', width=20, height=8)

##################################步骤7.marker基因点图########################
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

#################################步骤8.去除double之后细胞分群############
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

##################################步骤9.去double之后，画marker基因点图########################
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

######发现21簇还是double，去掉#######
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


##################################步骤9.去double之后，画marker基因点图########################
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

#################################步骤6.细胞的marker基因注释---############
#正负可以理解为上调、下调的意思；正marker表示在这个cluster中表达量高，而在其他的cluster中低；
#负marker表示在这个cluster中表达量低，而在其他的cluster中高。感觉用处不大。
#cluster0.markers <- FindMarkers(guolv_1_jw, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25)
  
#write.csv(cluster0.markers,file="cluster0.markers.csv")

#FindAllMarkers() 在细胞群之间比较，找到每个细胞群的正、负表达biomarker
#寻找每一cluster的marker
ALL.markers <- FindAllMarkers(guolv_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ALL.markers,file="ALL.markers.csv")
ALL.markers <- read.csv(file="ALL.markers.csv")

# 每个细胞群选前十个基因，所以一共有300个基因。这是官网的算法，有bug。顺序不对。
top10 <- ALL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="ALL.marker.top10.csv")
top50 <- ALL.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50,file="ALL.marker.top50.csv")

#########################步骤6.8.亚群重命名和改颜色-1###############
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


#################################步骤5.5. 根据WTMU和组织，画几个点图############
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

#########################步骤6.9 在细胞中，细胞的marker基因再注释--用于热图############
#正负可以理解为上调、下调的意思；正marker表示在这个cluster中表达量高，而在其他的cluster中低；
#负marker表示在这个cluster中表达量低，而在其他的cluster中高。感觉用处不大。
#cluster0.markers <- FindMarkers(guolv_1_jw, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(cluster0.markers,file="cluster0.markers.csv")
#FindAllMarkers() 在细胞群之间比较，找到每个细胞群的正、负表达biomarker
#寻找每一cluster的marker
table(guolv_2@meta.data$xibao)
ALL.markers_cells <- FindAllMarkers(guolv_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ALL.markers_cells,file="ALL.markers细胞命名之后.csv")
#ALL.markers_cells <- read.csv(file="ALL.markers_cells.csv")

# 每个细胞群选前50个基因
top10 <- ALL.markers_cells %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="ALL.markers细胞命名之后_top10.csv")
top50 <- ALL.markers_cells %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50,file="ALL.markers细胞命名之后_top50.csv")

##########################步骤7.细胞的marker做热图############
#slot = "scale.data",如果用counts,data。不好看。
#热图，默认slot = 'scale.data',所以有很多基因找不到。因为前面只用了2000个高可变基因做中心化。
#所以，我前面把所有的基因都做了中心化。生信技能书有代码能解决这个问题。对所有数据进行中心化。
#对所有基因进行中心化
guolv_3 <- guolv_2
all.genes <- rownames(guolv_3)
guolv_3 <- ScaleData(guolv_3, features = all.genes)
#每个细胞类型，选出来50个细胞。
#提取所有细胞的细胞类型，形成一个向量
suoyou_zhushi <- guolv_3@meta.data$xibao
#细胞类型的向量
xibao_names <- levels(guolv_3)
#每种细胞选50个，保存细胞在向量中的所在位置。
cluster50_all <- NULL
for (i in xibao_names) {
  cluster50_temp <- tail(head(which(suoyou_zhushi== i),50),50)
  cluster50_all <- union(cluster50_all, cluster50_temp)
}
#画图
plot300 <- DoHeatmap(guolv_3, 
                     features = top50$gene,
                     slot = "scale.data",
                     cells = cluster50_all,
                     group.colors = col_3,
                     draw.lines = F,
                     hjust = 0, #文字居中。
                     raster = F, #这个设置成F，图片就不会模糊
                     angle = 45)
#这个梯度，是设置的热图红色和蓝色的过渡。
plot300 <- plot300 + scale_fill_gradientn(colors = c("#00aedb","white","#d11141"))
ggsave(plot300, file='plot300-热图1.pdf', width=8, height=7)

##################################步骤9.细胞命名，画marker基因点图########################
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


##########步骤7.1细胞的marker做小提琴图###########
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



#########################步骤7.细胞数目统计和柱状图###############
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

#########################步骤7.1 命名之后，细胞数目统计和柱状图###############
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

