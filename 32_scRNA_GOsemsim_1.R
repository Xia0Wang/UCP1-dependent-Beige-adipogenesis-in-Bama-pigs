library(dplyr)
library(Seurat)
library(patchwork)
library(Rcpp)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(clusterProfiler)
library(org.Ss.eg.db)
library(stringr)
library(dplyr)
library(reshape2)

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/32_亚群")
adipoFap <- readRDS("AD+FAP_SeuratObject_Celltype.rds")

DefaultAssay(adipoFap)

table(adipoFap@meta.data$sample)
table(adipoFap@meta.data$Celltype)
table(adipoFap@meta.data$Celltype,adipoFap@meta.data$sample)


##########UMAP
adipoFap <- NormalizeData(adipoFap, normalization.method =  "LogNormalize",
                         scale.factor = 10000)
#查看:GetAssay(guolv_1,assay = "RNA")

#前2000个高变feature RNA，默认是2000个。
adipoFap <- FindVariableFeatures(adipoFap,
                                selection.method = "vst", 
                                nfeatures = 2000)

adipoFap <- ScaleData(adipoFap)
adipoFap <- RunPCA(object = adipoFap, pc.genes = VariableFeatures(adipoFap))

adipoFap <- RunUMAP(object = adipoFap, dims = 1:25)
adipoFap <- FindNeighbors(adipoFap, dims = 1:25)
#备注
adipoFap <- FindClusters(adipoFap, resolution = 0.3)


x <- 0.2
################################步骤5. DEGs 1#########################
#S11和 29906都是冷刺激后的
adipoFap_FGG <-  subset(adipoFap, sample == c("S11-FGG","S29906FGG"))
adipoFap_SZ <-  subset(adipoFap, sample == c("S11-SZ","S29906SZ"))


#依次计算DEGs--FGG
b <- 0
for (i in levels(adipoFap_FGG)[1:7]){
  b <- b + 1
  #wenjianm = paste("DEGs_in_cluster_",i,"_betweenWTMu",sep = "")
  DEGs_up <- FindMarkers(adipoFap_FGG, 
                         ident.1 = "S11-FGG", 
                         ident.2 = "S29906FGG",
                         group.by = "sample", 
                         subset.ident = i,
                         min.pct = x,
                         logfc.threshold = 0.5,
                         test.use = "wilcox",
                         only.pos = F
  ) %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0) 
  DEGs_up$genename <- row.names(DEGs_up)
  DEGs_up$genename2 <- substring(DEGs_up$genename,1,7)
  DEGs_up <- subset(DEGs_up, DEGs_up$genename2 != "ENSSSCG")
  if (nrow(DEGs_up) !=0){
    assign(paste("FGG_DEGs_",b,"_",i,"_UP",sep = ""), DEGs_up)
    write.csv(DEGs_up,file= paste("FGG_DEGs_",b," ",i,"_UP.csv",sep = "")) 
  }
  DEGs_down <- FindMarkers(adipoFap_FGG, 
                           ident.1 = "S11-FGG", 
                           ident.2 = "S29906FGG",
                           group.by = "sample", 
                           subset.ident = i,
                           min.pct = x,
                           logfc.threshold = 0.5, 
                           test.use = "wilcox",
                           only.pos = F
  ) %>%
    filter(p_val_adj < 0.05 & avg_log2FC < 0)
  DEGs_down$genename <- row.names(DEGs_down)
  DEGs_down$genename2 <- substring(DEGs_down$genename,1,7)
  DEGs_down <- subset(DEGs_down, DEGs_down$genename2 != "ENSSSCG")
  if (nrow(DEGs_down) !=0){
    assign(paste("FGG_DEGs_",b,"_",i,"_DOWN",sep = ""), DEGs_down)
    write.csv(DEGs_down,file= paste("FGG_DEGs_",b," ",i,"_DOWN.csv",sep = "")) 
  }
}



#依次计算DEGs--SZ
b <- 0
for (i in levels(adipoFap_SZ)[1:7]){
  b <- b + 1
  #wenjianm = paste("DEGs_in_cluster_",i,"_betweenWTMu",sep = "")
  DEGs_up <- FindMarkers(adipoFap_SZ, 
                         ident.1 = "S11-SZ", 
                         ident.2 = "S29906SZ",
                         group.by = "sample", 
                         subset.ident = i,
                         min.pct = x,
                         logfc.threshold = 0.5,
                         test.use = "wilcox",
                         only.pos = F
  ) %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0) 
  DEGs_up$genename <- row.names(DEGs_up)
  DEGs_up$genename2 <- substring(DEGs_up$genename,1,7)
  DEGs_up <- subset(DEGs_up, DEGs_up$genename2 != "ENSSSCG")
  if (nrow(DEGs_up) !=0){
    assign(paste("SZ_DEGs_",b,"_",i,"_UP",sep = ""), DEGs_up)
    write.csv(DEGs_up,file= paste("SZ_DEGs_",b," ",i,"_UP.csv",sep = "")) 
  }
  DEGs_down <- FindMarkers(adipoFap_SZ, 
                           ident.1 = "S11-SZ", 
                           ident.2 = "S29906SZ",
                           group.by = "sample", 
                           subset.ident = i,
                           min.pct = x,
                           logfc.threshold = 0.5, 
                           test.use = "wilcox",
                           only.pos = F
  ) %>%
    filter(p_val_adj < 0.05 & avg_log2FC < 0)
  DEGs_down$genename <- row.names(DEGs_down)
  DEGs_down$genename2 <- substring(DEGs_down$genename,1,7)
  DEGs_down <- subset(DEGs_down, DEGs_down$genename2 != "ENSSSCG")
  if (nrow(DEGs_down) !=0){
    assign(paste("SZ_DEGs_",b,"_",i,"_DOWN",sep = ""), DEGs_down)
    write.csv(DEGs_down,file= paste("SZ_DEGs_",b," ",i,"_DOWN.csv",sep = "")) 
  }
}


#####总体的 FGG########
DEGs_up <- FindMarkers(adipoFap_FGG, 
                       ident.1 = "S11-FGG", 
                       ident.2 = "S29906FGG",
                       group.by = "sample", 
                       subset.ident = c("LYA","TGA","LGA"),
                       min.pct = x,
                       logfc.threshold = 0.5,
                       test.use = "wilcox",
                       only.pos = F
) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) 
DEGs_up$genename <- row.names(DEGs_up)
DEGs_up$genename2 <- substring(DEGs_up$genename,1,7)
DEGs_up <- subset(DEGs_up, DEGs_up$genename2 != "ENSSSCG")
if (nrow(DEGs_up) !=0){
  write.csv(DEGs_up,file= paste("FGG_DEGs_","总体","_UP.csv",sep = "")) 
}

DEGs_down <- FindMarkers(adipoFap_FGG, 
                       ident.1 = "S11-FGG", 
                       ident.2 = "S29906FGG",
                       group.by = "sample", 
                       subset.ident = c("LYA","TGA","LGA"),
                       min.pct = x,
                       logfc.threshold = 0.5,
                       test.use = "wilcox",
                       only.pos = F
) %>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) 
DEGs_down$genename <- row.names(DEGs_down)
DEGs_down$genename2 <- substring(DEGs_down$genename,1,7)
DEGs_down <- subset(DEGs_down, DEGs_down$genename2 != "ENSSSCG")
if (nrow(DEGs_down) !=0){
  write.csv(DEGs_down,file= paste("FGG_DEGs_","总体","_DOWN.csv",sep = "")) 
}



#####总体的SZ########
DEGs_up <- FindMarkers(adipoFap_SZ, 
                       ident.1 = "S11-SZ", 
                       ident.2 = "S29906SZ",
                       group.by = "sample", 
                       subset.ident = c("LYA","TGA","LGA"),
                       min.pct = x,
                       logfc.threshold = 0.5,
                       test.use = "wilcox",
                       only.pos = F
) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) 
DEGs_up$genename <- row.names(DEGs_up)
DEGs_up$genename2 <- substring(DEGs_up$genename,1,7)
DEGs_up <- subset(DEGs_up, DEGs_up$genename2 != "ENSSSCG")
if (nrow(DEGs_up) !=0){
  assign(paste("SZ_DEGs_","zongti","_UP",sep = ""), DEGs_up)
  write.csv(DEGs_up,file= paste("SZ_DEGs_","总体","_UP.csv",sep = "")) 
}

DEGs_down <- FindMarkers(adipoFap_SZ, 
                         ident.1 = "S11-SZ", 
                         ident.2 = "S29906SZ",
                         group.by = "sample", 
                         subset.ident = c("LYA","TGA","LGA"),
                         min.pct = x,
                         logfc.threshold = 0.5,
                         test.use = "wilcox",
                         only.pos = F
) %>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) 
DEGs_down$genename <- row.names(DEGs_down)
DEGs_down$genename2 <- substring(DEGs_down$genename,1,7)
DEGs_down <- subset(DEGs_down, DEGs_down$genename2 != "ENSSSCG")
if (nrow(DEGs_down) !=0){
  assign(paste("SZ_DEGs_","zongti","_DOWN",sep = ""), DEGs_down)
  write.csv(DEGs_down,file= paste("SZ_DEGs_","总体","_DOWN.csv",sep = "")) 
}


library(ggvenn)
setwd("E:/BIIF/liujiali/zangzhu/3_结果/34_维恩")
#############维恩图###############
LTX_zongti_UP <- read.csv("pWAT_DEGs_1 Adipocyte_UP.csv")
LTX_zongti_UP$genename2 <- substring(LTX_zongti_UP$X,1,3)
LTX_zongti_UP <- subset(LTX_zongti_UP, LTX_zongti_UP$genename2 != "LOC")
LTX_zongti_DOWN <- read.csv("pWAT_DEGs_1 Adipocyte_DOWN.csv")
LTX_zongti_DOWN$genename2 <- substring(LTX_zongti_DOWN$X,1,3)
LTX_zongti_DOWN <- subset(LTX_zongti_DOWN, LTX_zongti_DOWN$genename2 != "LOC")
LTX_C4_UP <- read.csv("SZ_DEGs_4 C4_UP.csv")
LTX_C4_UP$genename2 <- substring(LTX_C4_UP$X,1,3)
LTX_C4_UP <- subset(LTX_C4_UP, LTX_C4_UP$genename2 != "LOC")
LTX_C4_DOWN <- read.csv("SZ_DEGs_4 C4_DOWN.csv")
LTX_C4_DOWN$genename2 <- substring(LTX_C4_DOWN$X,1,3)
LTX_C4_DOWN <- subset(LTX_C4_DOWN, LTX_C4_DOWN$genename2 != "LOC")


plot_venn <- ggvenn(data = list(LJL_zongti_UP = row.names(SZ_DEGs_zongti_UP),
                                LTX_zongti_UP = LTX_zongti_UP$X),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("维恩图_总体_up.pdf",sep = ""), width=3, height=3)
RP1 <- subset(LTX_zongti_UP,LTX_zongti_UP$X %in% row.names(SZ_DEGs_zongti_UP))
write.csv(RP1,"维恩图_总体_up.csv")


plot_venn <- ggvenn(data = list(LJL_zongti_DOWN = row.names(SZ_DEGs_zongti_DOWN),
                                LTX_zongti_DOWN = LTX_zongti_DOWN$X),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("维恩图_总体_DOWN.pdf",sep = ""), width=3, height=3)
RP2 <- subset(LTX_zongti_DOWN,LTX_zongti_DOWN$X %in% row.names(SZ_DEGs_zongti_DOWN))
write.csv(RP2,"维恩图_总体_DOWN.csv")



plot_venn <- ggvenn(data = list(LJL_TGA_UP = row.names(SZ_DEGs_5_TGA_UP),
                                LTX_C4_UP = LTX_C4_UP$X),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("维恩图_C4_up.pdf",sep = ""), width=3, height=3)
RP3 <- subset(LTX_C4_UP,LTX_C4_UP$X %in% row.names(SZ_DEGs_5_TGA_UP))
write.csv(RP3,"维恩图_C4_up.csv")


plot_venn <- ggvenn(data = list(LJL_TGA_DOWN = row.names(SZ_DEGs_5_TGA_DOWN),
                                LTX_C4_DOWN = LTX_C4_DOWN$X),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("维恩图_C4_DOWN.pdf",sep = ""), width=3, height=3)
RP4 <- subset(LTX_C4_DOWN,LTX_C4_DOWN$X %in% row.names(SZ_DEGs_5_TGA_DOWN))
write.csv(RP4,"维恩图_C4_DOWN.csv")

#自己造一个阶乘的函数
jiecheng <- function(m,n){
  r <- 1
  for (i in 1:n){
    r <- r*m
    m <- m-1
  }
  return(r)
}
#实验
jiecheng(5,3)

#算概率方法，具体见图片
#这是只约分右边
jieguo <- (choose(13349-b,a)*choose(13349,b)*jiecheng(b+c,b))/(choose(13349,a+b)*jiecheng(13349,b))
#这是左右都约分
jieguo2 <- (jiecheng(b+a,b)*choose(13349,b)*jiecheng(b+c,b))/(jiecheng(13349,b)*jiecheng(13349,b))
jieguo2

#所有的脂肪，UP
a <- 17
b <- 6
c <- 242
jieguo <- (choose(13349-b,a)*choose(13349,b)*jiecheng(b+c,b))/(choose(13349,a+b)*jiecheng(13349,b))
jieguo
#3.909624e-06

#所有的脂肪，down
a <- 19
b <- 9
c <- 244
jieguo <- (choose(13349-b,a)*choose(13349,b)*jiecheng(b+c,b))/(choose(13349,a+b)*jiecheng(13349,b))
jieguo
#1.892344e-09

#产热脂肪up
a <- 21
b <- 3
c <- 229
jieguo <- (choose(13349-b,a)*choose(13349,b)*jiecheng(b+c,b))/(choose(13349,a+b)*jiecheng(13349,b))
jieguo
#0.01049034

#产热脂肪 down
a <- 23
b <- 6
c <- 271
jieguo <- (choose(13349-b,a)*choose(13349,b)*jiecheng(b+c,b))/(choose(13349,a+b)*jiecheng(13349,b))
jieguo
#3.595102e-05



############看 IL1相关基因##########
#修改
markers.to.plot <- c("IL1RAP","IL1RN","IL1B","IL1A","IL1R1","IL1RAPL2","IL1R2")
#总体
dotplot <- DotPlot(adipoFap, 
                   features = markers.to.plot, 
                   group.by = "sample",
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  labs(y = NULL,x = NULL)+
  RotatedAxis()
ggsave(dotplot, file='plot300-点图_RNA.pdf', width=6, height=3)


#修改
markers.to.plot <- c("IL1RAP","IL1RN","IL1B","IL1A","IL1R1","IL1RAPL2","IL1R2")
#总体
dotplot <- DotPlot(adipoFap, 
                   features = markers.to.plot, 
                   #group.by = "sample",
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  labs(y = NULL,x = NULL)+
  RotatedAxis()
ggsave(dotplot, file='plot301-点图_RNA.pdf', width=6, height=3)


#############
plotMarker1 <- FeaturePlot(adipoFap, features = c("IL1RAP"),split.by ='sample')
ggsave(plotMarker1, file='Rplot295-Marker-IL1RAP.pdf', width=24, height=4)

plotMarker1 <- FeaturePlot(adipoFap, features = c("ARNTL"),split.by ='sample')
ggsave(plotMarker1, file='Rplot295-Marker-ARNTL.pdf', width=24, height=4)



