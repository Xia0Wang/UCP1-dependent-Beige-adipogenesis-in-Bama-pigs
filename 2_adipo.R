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
setwd("E:/BIIF/liutianxia/Adipocyte_scRNA_seq_8/3_结果/12.0_亚群_纯adipo")
getwd()

#########step1######################
#Adipo
levels(guolv_2)
yaqun_1_Adipo <- subset(x = guolv_2, xibao %in% c("Adipocyte"))
#integrated
DefaultAssay(yaqun_1_Adipo) <- "integrated"

set.seed(123)
yaqun_1_Adipo <- RunUMAP(object = yaqun_1_Adipo, dims = 1:30)
yaqun_1_Adipo <- FindNeighbors(yaqun_1_Adipo, dims = 1:30)
yaqun_1_Adipo <- FindClusters(yaqun_1_Adipo, resolution = 0.6)
DefaultAssay(yaqun_1_Adipo) <- "RNA"
table(yaqun_1_Adipo$orig.ident)

x <- 1
plot113 <- DimPlot(yaqun_1_Adipo,reduction = "umap",label=F, group.by = "orig.ident", pt.size = x,
                   cols = DiscretePalette(36, palette = "polychrome")[c(6,9,7,8,18,12)])
ggsave(plot113, file='Plot113-umap-1.pdf', width=5.5, height=4)
plot114 <- DimPlot(yaqun_1_Adipo,reduction = "umap",label=F, group.by = "orig.ident", split.by ='orig.ident', pt.size = x,
                   cols = DiscretePalette(36, palette = "polychrome")[c(6,9,7,8,18,12)])
ggsave(plot114, file='Plot114-umap-2.pdf', width=18, height=4)
plot115 <- DimPlot(yaqun_1_Adipo,reduction = "umap",label=T, pt.size = x,
                   cols = DiscretePalette(36, palette = "polychrome")[c(4:36)])
ggsave(plot115, file='Plot115-umap-3!!!.pdf', width=5, height=4)
plot116 <- DimPlot(yaqun_1_Adipo,reduction = "umap",label=T,split.by ='orig.ident', pt.size = x,
                   cols = DiscretePalette(36, palette = "polychrome")[c(4:36)])
ggsave(plot116, file='Plot116-umap-4.pdf', width=18, height=4)
plot117 <- DimPlot(yaqun_1_Adipo,reduction = "umap",label=T,split.by ='weizhi', pt.size = x,
                   cols = DiscretePalette(36, palette = "polychrome")[c(4:36)])
ggsave(plot117, file='Plot117-umap-5.pdf', width=8, height=4)
plot118 <- DimPlot(yaqun_1_Adipo,reduction = "umap",label=T,split.by ='chuli', pt.size = x,
                   cols = DiscretePalette(36, palette = "polychrome")[c(4:36)])
ggsave(plot118, file='Plot118-umap-6.pdf', width=12, height=4)

mean_Count <- aggregate(yaqun_1_Adipo@meta.data, nCount_RNA ~ integrated_snn_res.0.6,mean)
mean_Feature <- aggregate(yaqun_1_Adipo@meta.data, nFeature_RNA ~ integrated_snn_res.0.6,mean)

#########step2######################
cell.numbers <- table(yaqun_1_Adipo@meta.data$integrated_snn_res.0.6,yaqun_1_Adipo@meta.data$orig.ident)
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
#write.csv(cell.numbers,"7.4.1.csv")
#转换数据框
cell.numbers <- as.data.frame(cell.numbers) 
cell.numbers$clusters <- row.names(cell.numbers)
#宽表变长表
cell.numbers <- reshape2::melt(cell.numbers, id.vars = "clusters")
#折线图中横坐标必须为离散型数值和数值型数值，否则会画不出来！！！！！
#cell.numbers$clusters <- as.numeric(cell.numbers$clusters)
cell.numbers$clusters <- factor(cell.numbers$clusters, levels = levels(yaqun_1_Adipo))
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
ggsave(plot408, file='plot501-num-adipo-1.pdf', width=5, height=4)

#########step3######################
DefaultAssay(yaqun_1_Adipo) <- "RNA"
ALL.markers_Adipo <- FindAllMarkers(yaqun_1_Adipo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ALL.markers_Adipo,file="ALL.markers单独Adipo.csv")
#top
top10 <- ALL.markers_Adipo %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="ALL.markers单独Adipo_top10.csv")
top50 <- ALL.markers_Adipo %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50,file="ALL.markers单独Adipo_top50.csv")

#########step4######################
markers.to.plot1 <- c("ADIPOQ", "PPARG","PLIN1",
                     "COL1A2","COL1A1","COL28A1","NEGR1",
                     "IGF2R","ACLY","ACSL1","ACACA","ELOVL6",
                     "LIPA","LRP3","APOE","LPIN1","PDK4","PCK1","LIPE","VEGFA","ANGPTL4","PNPLA2",
                     "SYNE2","DHRS4","UCP3","PRDM16")
markers.to.plot2 <- c("PRDM16","CIDEA","PPARGC1A","IRF4","CS","ACE2","VEGFA","MAPK14","GLUL","ANGPT4","ANGPTL4",
                      "TGFBR3","THBS1","ANG","GRB10","ZFAND5","FOXO1","IRS1","SPRY1","F3",
                      "PLPP3","ABHD5","SOCS2","PTEN","LPL","DHRS4","FASN","IDH1","COL3A1","LUM")
markers.to.plot3 <- c("ACE2","EBF2","THRB","GATM","DHRS4","SGK2","IRF4","PRDM16",
                      "ABHD5","IGFBP5","ELOVL6","FASN",
                      "VEGFA","NRP1","FGFR2","MAPK14","GLUL","B4GALT1","F3","ANGPTL4",
                      "GK","PCK1","IDH1","PDK4","DDIT4","UGP2",
                      "EGR1","PTGIS","FAM162A","GRB10","PID1","SLC25A33","PIK3R1","FOXO1")
markers.to.plot4 <- c("PPARG", "ACACA","ELOVL6","ACLY","IGF2R",
                      "NNAT","LRP3","ABCD2","APOE","ABCG1","TREM2","CD36",
                      "HIF1A","NEDD9","GADD45G","LEP")
#总体
dotplot <- DotPlot(yaqun_1_Adipo, 
                   features = markers.to.plot1, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪亚群1.pdf', width=10, height=4)
dotplot <- DotPlot(yaqun_1_Adipo, 
                   features = markers.to.plot2, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪亚群2.pdf', width=10, height=4)
dotplot <- DotPlot(yaqun_1_Adipo, 
                   features = markers.to.plot3, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪亚群3.pdf', width=10, height=4)
dotplot <- DotPlot(yaqun_1_Adipo, 
                   features = markers.to.plot4, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪亚群4.pdf', width=10, height=4)

#########step5######################
#thermogenesis
markers.to.plot5 <- c("GRB10","ACE2","PRDM16","EBF2","CITED1",
                      "PPARGC1A","TMEM26","IRF4","TBX1","FOXO1",
                      "PDGFA","PPARG","ADIPOQ","HOXC8","HOXC9","PDGFRA","ABCA6")
len <- length(markers.to.plot5)
gene_cell_exp_all <- AverageExpression(yaqun_1_Adipo,
                                      features = markers.to.plot5,
                                      group.by = 'integrated_snn_res.0.6',
                                      assays = "RNA",
                                      slot = 'data') 
EXP_ALL2 <- gene_cell_exp_all[["RNA"]]
#z-score
RS_F0_15_2 = t(scale(t(EXP_ALL2),center = TRUE, scale = TRUE))
#限定上限，使表达量大于2的等于2
RS_F0_15_2[RS_F0_15_2> 2] <- 2 
#限定下限，使表达量小于-2的等于-2
RS_F0_15_2[RS_F0_15_2< -2] <-  -2 
#热图
library(pheatmap)
bb <- pheatmap(RS_F0_15_2,
               cluster_cols=F, #不对列进行聚类
               cluster_rows=F,
               show_rownames = T, #不显示行名
               border_color = "white",
               cellwidth = 15,
               cellheight = 15,
               color = colorRampPalette(colors =  c("#372ac1","#3060e9","#378afb","#59acfd","#aad5ff","#dcebfb",
                                                    "#FFFFFF","#fde3e4","#ffb5ba","#fd6f6e","#f52a27","#e12112","#d01c05"))(10000))
ggsave(bb, file='plot热图-产热marker-all.pdf', width=5, height=6)

#########step6######################
#产热相关
markers.to.plot6 <- c("ACE2","EBF2","THRB","GATM","DHRS4","SGK2","IRF4","PRDM16",
                      "CD36","LIPE","ABHD5","IGFBP5","ELOVL6","FASN","DGAT2","ACLY","ACSL1","ACACA",
                      "VEGFA","NRP1","FGFR2","MAPK14","GLUL","B4GALT1","F3","ANGPTL4",
                      "GK","PCK1","IDH1","PDK4","DDIT4","UGP2",
                      "EGR1","PTGIS","FAM162A","GRB10","PID1","SLC25A33","PIK3R1","FOXO1")
len <- length(markers.to.plot6)
gene_cell_exp_all <- AverageExpression(yaqun_1_Adipo,
                                       features = markers.to.plot6,
                                       group.by = 'integrated_snn_res.0.6',
                                       assays = "RNA",
                                       slot = 'data') 
EXP_ALL2 <- gene_cell_exp_all[["RNA"]]
#进行z-score
RS_F0_15_2 = t(scale(t(EXP_ALL2),center = TRUE, scale = TRUE))
#限定上限，使表达量大于2的等于2
RS_F0_15_2[RS_F0_15_2> 2] <- 2 
#限定下限，使表达量小于-2的等于-2
RS_F0_15_2[RS_F0_15_2< -2] <-  -2 
#热图
library(pheatmap)
bb <- pheatmap(RS_F0_15_2,
               cluster_cols=F, #不对列进行聚类
               cluster_rows=F,
               show_rownames = T, #不显示行名
               border_color = "white",
               cellwidth = 15,
               cellheight = 15,
               color = colorRampPalette(colors =  c("#372ac1","#3060e9","#378afb","#59acfd","#aad5ff","#dcebfb",
                                                    "#FFFFFF","#fde3e4","#ffb5ba","#fd6f6e","#f52a27","#e12112","#d01c05"))(10000))
ggsave(bb, file='plot热图-产热marker-all2.pdf', width=5, height=10)

#########step7######################
yaqun_1_Adipo_sWAT <-  subset(yaqun_1_Adipo, weizhi == "sWAT")
yaqun_1_Adipo_sWAT <-  subset(yaqun_1_Adipo, weizhi == "sWAT")
#依次计算DEGs--sWAT
b <- 0
for (i in levels(yaqun_1_Adipo_sWAT)){
  b <- b + 1
  #wenjianm = paste("DEGs_in_cluster_",i,"_betweenWTMu",sep = "")
  DEGs_up <- FindMarkers(yaqun_1_Adipo_sWAT, 
                         ident.1 = "KI", 
                         ident.2 = "WT",
                         group.by = "chuli", 
                         subset.ident = i,
                         min.pct = 0.3,
                         logfc.threshold = 0.5,
                         test.use = "wilcox",
                         only.pos = F
  ) %>%
    filter(p_val < 0.05 & avg_log2FC > 0) 
  if (nrow(DEGs_up) !=0){
    assign(paste("JJ_DEGs_",b,"_",i,"_UP",sep = ""), DEGs_up)
    write.csv(DEGs_up,file= paste("JJ_DEGs_",b," ",i,"_UP.csv",sep = "")) 
  }
  DEGs_down <- FindMarkers(yaqun_1_Adipo_sWAT, 
                           ident.1 = "KI", 
                           ident.2 = "WT",
                           group.by = "chuli", 
                           subset.ident = i,
                           min.pct = 0.3,
                           logfc.threshold = 0.5, 
                           test.use = "wilcox",
                           only.pos = F
  ) %>%
    filter(p_val < 0.05 & avg_log2FC < 0)
  if (nrow(DEGs_down) !=0){
    assign(paste("sWAT_DEGs_",b,"_",i,"_DOWN",sep = ""), DEGs_down)
    write.csv(DEGs_down,file= paste("JJ_DEGs_",b," ",i,"_DOWN.csv",sep = "")) 
  }
}

#依次计算DEGs--sWAT
b <- 0
for (i in levels(yaqun_1_Adipo_sWAT)){
  b <- b + 1
  #wenjianm = paste("DEGs_in_cluster_",i,"_betweenWTMu",sep = "")
  DEGs_up <- FindMarkers(yaqun_1_Adipo_sWAT, 
                         ident.1 = "KI", 
                         ident.2 = "WT",
                         group.by = "chuli", 
                         subset.ident = i,
                         min.pct = 0.3,
                         logfc.threshold = 0.5, 
                         test.use = "wilcox",
                         only.pos = F
  ) 
  #这块稍微调整了一下，是因为2簇没有DEGs。这样就可以继续循环了。
  if (nrow(DEGs_up) !=0){
    DEGs_up <- filter(DEGs_up,p_val < 0.05 & avg_log2FC > 0) 
    assign(paste("SZ_DEGs_",b,"_",i,"_UP",sep = ""), DEGs_up)
    write.csv(DEGs_up,file= paste("SZ_DEGs_",b," ",i,"_UP.csv",sep = "")) 
  }
  
  DEGs_down <- FindMarkers(yaqun_1_Adipo_sWAT, 
                           ident.1 = "KI", 
                           ident.2 = "WT",
                           group.by = "chuli", 
                           subset.ident = i,
                           min.pct = 0.3,
                           logfc.threshold = 0.5, 
                           test.use = "wilcox",
                           only.pos = F
  ) 
  if (nrow(DEGs_down) !=0){
    DEGs_down <- filter(DEGs_down,p_val < 0.05 & avg_log2FC < 0) 
    assign(paste("sWAT_DEGs_",b,"_",i,"_DOWN",sep = ""), DEGs_down)
    write.csv(DEGs_down,file= paste("SZ_DEGs_",b," ",i,"_DOWN.csv",sep = ""))
  }
}


yaqun_1_Adipo_sWAT <- subset(yaqun_1_Adipo, weizhi =="sWAT")
yaqun_1_Adipo_sWAT <- subset(yaqun_1_Adipo, weizhi =="sWAT")

#########step8######################
DefaultAssay(yaqun_1_Adipo_sWAT)
yaqun_1_Adipo_sWAT@meta.data$zhifangyaqun <- yaqun_1_Adipo_sWAT@active.ident
yaqun_1_Adipo_sWAT_WT <- subset(x = yaqun_1_Adipo_sWAT, chuli == "WT")
yaqun_1_Adipo_sWAT_KI <- subset(x = yaqun_1_Adipo_sWAT, chuli == "KI")
levels(yaqun_1_Adipo_sWAT)
#产热相关

markers.to.plot5 <- c("GRB10","ACE2","PRDM16","EBF2","CITED1","CD44",
                      "PPARGC1A","TMEM26","IRF4","TBX1","FOXO1",
                      "PDGFA","PPARG","ADIPOQ","HOXC8","HOXC9","PDGFRA","ABCA6")

len <- length(markers.to.plot5)

gene_cell_exp_WT <- AverageExpression(yaqun_1_Adipo_sWAT_WT,
                                      features = markers.to.plot5,
                                      group.by = 'zhifangyaqun',
                                      assays = "RNA",
                                      slot = 'data') 
gene_cell_exp_KI <- AverageExpression(yaqun_1_Adipo_sWAT_KI,
                                      features = markers.to.plot5,
                                      group.by = 'zhifangyaqun',
                                      assays = "RNA",
                                      slot = 'data') 

#宽变长
EXP_WT <- reshape2::melt(gene_cell_exp_WT[["RNA"]])
EXP_KI <- reshape2::melt(gene_cell_exp_KI[["RNA"]])
#组合
EXP_ALL <- cbind(EXP_WT,EXP_KI)
EXP_ALL <- EXP_ALL[,c(3,6)]
colnames(EXP_ALL) <- c("WT","KI")
#为了使，一个亚群里，各个样品都一起z-score，再变格式
EXP_ALL_1 <- EXP_ALL[1:len,]
EXP_ALL_2 <- EXP_ALL[(len+1):(len*2),]
EXP_ALL_3 <- EXP_ALL[(len*2+1):(len*3),]
EXP_ALL_4 <- EXP_ALL[(len*3+1):(len*4),]
EXP_ALL_5 <- EXP_ALL[(len*4+1):(len*5),]
EXP_ALL_6 <- EXP_ALL[(len*5+1):(len*6),]
EXP_ALL_7 <- EXP_ALL[(len*6+1):(len*7),]
EXP_ALL_8 <- EXP_ALL[(len*7+1):(len*8),]

EXP_ALL2 <- cbind(cbind(cbind(cbind(cbind(cbind(cbind(EXP_ALL_1,EXP_ALL_2),EXP_ALL_3),EXP_ALL_4),EXP_ALL_5),EXP_ALL_6),EXP_ALL_7),EXP_ALL_8)
row.names(EXP_ALL2) <- markers.to.plot5
#进行z-score
RS_F0_15_2 = t(scale(t(EXP_ALL2),center = TRUE, scale = TRUE))
#限定上限，使表达量大于2的等于2
RS_F0_15_2[RS_F0_15_2> 2] <- 2 
#限定下限，使表达量小于-2的等于-2
RS_F0_15_2[RS_F0_15_2< -2] <-  -2 
#热图
library(pheatmap)
bb <- pheatmap(RS_F0_15_2,
               cluster_cols=F, #不对列进行聚类
               cluster_rows=F,
               show_rownames = T, #不显示行名
               border_color = "white",
               cellwidth = 15,
               cellheight = 15,
               color = colorRampPalette(colors =  c("#372ac1","#3060e9","#378afb","#59acfd","#aad5ff","#dcebfb",
                                                    "#FFFFFF","#fde3e4","#ffb5ba","#fd6f6e","#f52a27","#e12112","#d01c05"))(10000))
ggsave(bb, file='RNA_seq_产热marker-sWAT.pdf', width=5, height=6)

#########step9######################
DefaultAssay(yaqun_1_Adipo_sWAT)
yaqun_1_Adipo_sWAT@meta.data$zhifangyaqun <- yaqun_1_Adipo_sWAT@active.ident
yaqun_1_Adipo_sWAT_WT <- subset(x = yaqun_1_Adipo_sWAT, chuli == "WT")
yaqun_1_Adipo_sWAT_KI <- subset(x = yaqun_1_Adipo_sWAT, chuli == "KI")
levels(yaqun_1_Adipo_sWAT)
#产热相关

markers.to.plot5 <- c("GRB10","ACE2","PRDM16","EBF2","CITED1","CD44",
                      "PPARGC1A","TMEM26","IRF4","TBX1","FOXO1",
                      "PDGFA","PPARG","ADIPOQ","HOXC8","HOXC9","PDGFRA","ABCA6")

len <- length(markers.to.plot5)

gene_cell_exp_WT <- AverageExpression(yaqun_1_Adipo_sWAT_WT,
                                      features = markers.to.plot5,
                                      group.by = 'zhifangyaqun',
                                      assays = "RNA",
                                      slot = 'data') 
gene_cell_exp_KI <- AverageExpression(yaqun_1_Adipo_sWAT_KI,
                                      features = markers.to.plot5,
                                      group.by = 'zhifangyaqun',
                                      assays = "RNA",
                                      slot = 'data') 

#宽变长
EXP_WT <- reshape2::melt(gene_cell_exp_WT[["RNA"]])
EXP_KI <- reshape2::melt(gene_cell_exp_KI[["RNA"]])
#组合
EXP_ALL <- cbind(EXP_WT,EXP_KI)
EXP_ALL <- EXP_ALL[,c(3,6)]
colnames(EXP_ALL) <- c("WT","KI")
#为了使，一个亚群里，各个样品都一起z-score，再变格式
EXP_ALL_1 <- EXP_ALL[1:len,]
EXP_ALL_2 <- EXP_ALL[(len+1):(len*2),]
EXP_ALL_3 <- EXP_ALL[(len*2+1):(len*3),]
EXP_ALL_4 <- EXP_ALL[(len*3+1):(len*4),]
EXP_ALL_5 <- EXP_ALL[(len*4+1):(len*5),]
EXP_ALL_6 <- EXP_ALL[(len*5+1):(len*6),]
EXP_ALL_7 <- EXP_ALL[(len*6+1):(len*7),]
EXP_ALL_8 <- EXP_ALL[(len*7+1):(len*8),]

EXP_ALL2 <- cbind(cbind(cbind(cbind(cbind(cbind(cbind(EXP_ALL_1,EXP_ALL_2),EXP_ALL_3),EXP_ALL_4),EXP_ALL_5),EXP_ALL_6),EXP_ALL_7),EXP_ALL_8)
row.names(EXP_ALL2) <- markers.to.plot5
#进行z-score
RS_F0_15_2 = t(scale(t(EXP_ALL2),center = TRUE, scale = TRUE))
#限定上限，使表达量大于2的等于2
RS_F0_15_2[RS_F0_15_2> 2] <- 2 
#限定下限，使表达量小于-2的等于-2
RS_F0_15_2[RS_F0_15_2< -2] <-  -2 
#热图
library(pheatmap)
bb <- pheatmap(RS_F0_15_2,
               cluster_cols=F, #不对列进行聚类
               cluster_rows=F,
               show_rownames = T, #不显示行名
               border_color = "white",
               cellwidth = 15,
               cellheight = 15,
               color = colorRampPalette(colors =  c("#372ac1","#3060e9","#378afb","#59acfd","#aad5ff","#dcebfb",
                                                    "#FFFFFF","#fde3e4","#ffb5ba","#fd6f6e","#f52a27","#e12112","#d01c05"))(10000))
ggsave(bb, file='RNA_seq_产热marker-sWAT.pdf', width=5, height=6)


setwd("E:/BIIF/liutianxia/Adipocyte_scRNA_seq_8/3_结果/12.1_亚群_纯adipo_整合之后")
#########step10######################
#通过分析，决定合并一些亚簇，然后重新命名
#成熟脂肪：7（c1.米脂）；2,3,5（c2.脂肪发生）0,1（c3，脂肪水解）；4,6（c4.脂肪产热）；

yaqun_1_Adipo_mm <- yaqun_1_Adipo
new.cluster.ids <- c("C3","C3","C2","C2", "C4","C2","C4","C1")
names(new.cluster.ids) <- levels(yaqun_1_Adipo_mm)
levels(yaqun_1_Adipo_mm)
yaqun_1_Adipo_mm <- RenameIdents(yaqun_1_Adipo_mm, new.cluster.ids)
levels(yaqun_1_Adipo_mm)
yaqun_1_Adipo_mm@active.ident <- factor(yaqun_1_Adipo_mm@active.ident,
                                        levels = c("C1","C2","C3","C4"))
levels(yaqun_1_Adipo_mm)
table(yaqun_1_Adipo_mm@active.ident)
yaqun_1_Adipo_mm@meta.data$yaqun <- yaqun_1_Adipo_mm@active.ident
DefaultAssay(yaqun_1_Adipo_mm) 
#设置颜色，跟之前对应，取代表性颜色。
col_3 <- DiscretePalette(36, palette = "polychrome")[c(4:36)]
#不带legend的图
y <- 1
plot33 <- DimPlot(yaqun_1_Adipo_mm, reduction = "umap", label = T, pt.size = y,
                  cols = col_3) + NoLegend()
ggsave(plot33, file='plot133-umap-25-合并之后.pdf', width=4.5, height=4)
plot34 <- DimPlot(yaqun_1_Adipo_mm, reduction = "umap", pt.size = y,
                  label = T, split.by ='orig.ident',
                  cols = col_3) + NoLegend()
ggsave(plot34, file='plot134-umap-26.pdf', width=16, height=4)
plot35 <- DimPlot(yaqun_1_Adipo_mm, reduction = "umap", label = F,pt.size = y,
                  cols = col_3)
ggsave(plot35, file='plot135-umap-27.pdf', width=5, height=4)
plot36 <- DimPlot(yaqun_1_Adipo_mm, reduction = "umap", pt.size = y,
                  label = F, split.by ='orig.ident',
                  cols = col_3)
ggsave(plot36, file='plot136-umap-28.pdf', width=16, height=4)

#########step11######################
cell.numbers <- table(yaqun_1_Adipo_mm@meta.data$yaqun,yaqun_1_Adipo_mm@meta.data$orig.ident)
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
write.csv(cell.numbers,"plot502-num-adipo命名之后-1.csv")
#转换数据框
cell.numbers <- as.data.frame(cell.numbers) 
cell.numbers$clusters <- row.names(cell.numbers)
#宽表变长表
cell.numbers <- reshape2::melt(cell.numbers, id.vars = "clusters")
#折线图中横坐标必须为离散型数值和数值型数值，否则会画不出来！！！！！
#cell.numbers$clusters <- as.numeric(cell.numbers$clusters)
cell.numbers$clusters <- factor(cell.numbers$clusters, levels = levels(yaqun_1_Adipo_mm))
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
        legend.position = c(0.9,0.8))+ #修改图例位置
  labs(y = "Average fraction of all nuclei")+ #title = "muscle"
  scale_y_continuous(expand = c(0,0))#+
#scale_y_continuous(expand = c(0,0),limits = c(0,4000)) #修改Y轴跨度
ggsave(plot408, file='plot502-num-adipo命名之后-1.pdf', width=5, height=4)

#########step12######################
cell.numbers <- table(yaqun_1_Adipo_mm@active.ident)
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
#write.csv(cell.numbers,"plot502-num-adipo命名之后-2-总体的.csv")
#转换数据框
cell.numbers <- as.data.frame(cell.numbers) 
cell.numbers$clusters <- row.names(cell.numbers)
#宽表变长表
cell.numbers <- reshape2::melt(cell.numbers, id.vars = "clusters")
#折线图中横坐标必须为离散型数值和数值型数值，否则会画不出来！！！！！
#cell.numbers$clusters <- as.numeric(cell.numbers$clusters)
cell.numbers$clusters <- factor(cell.numbers$clusters, levels = levels(yaqun_1_Adipo_mm))
#转换为因子
cell.numbers$variable <- factor(cell.numbers$variable)
cell.numbers$value <- as.numeric(cell.numbers$value)

#画柱状图
plot408 <- ggplot(cell.numbers, aes(x = clusters, y = value, fill = clusters) ) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), colour = "white", size=0, width = 0.7)+
  scale_fill_manual(values = DiscretePalette(36, palette = "polychrome")[c(4:20)])+
  theme_classic()+
  #geom_text(aes(label=value),position = position_dodge(0.8),size=1.5,vjust=-0.5)+ #增加柱状图上的数字
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, size = 6), #修改X轴文字
        )+ #修改图例位置legend.position = c(1,1)
  labs(y = "Average fraction of all nuclei")+ #title = "muscle"
  scale_y_continuous(expand = c(0,0))#+
#scale_y_continuous(expand = c(0,0),limits = c(0,4000)) #修改Y轴跨度
ggsave(plot408, file='plot502-num-adipo命名之后-2-总体的.pdf', width=4, height=4)

#########step13######################
#DefaultAssay(yaqun_1_Adipo_mm) <- "RNA"
ALL.markers_Adipo <- FindAllMarkers(yaqun_1_Adipo_mm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ALL.markers_Adipo,file="ALL.markers_单独Adipo_mm.csv")
#以后就不用每次都计算了，时间有点长，读取之前的计算结果即可。
#ALL.markers_Adipo <- read.csv(file="ALL.markers_单独Adipo_mm_RNA.csv")

# 每个细胞群选前十个基因，所以一共有300个基因。这是官网的算法，有bug。顺序不对。
top10 <- ALL.markers_Adipo %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="ALL.markers_单独Adipo_mm_top10.csv")
top50 <- ALL.markers_Adipo %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50,file="ALL.markers_单独Adipo_mm_top50.csv")

#########step14######################
#读取marker文件
#ALL.markers_Adipo <- read.csv(file="ALL.markers_单独Adipo_mm_integrated.csv")
xibao_names <- levels(yaqun_1_Adipo_mm)
b <- 0
for (i in xibao_names){
  b <- b + 1
  cluster_top <- head(subset(ALL.markers_Adipo, ALL.markers_Adipo$cluster == i),100)
  #提取差异表达基因SYMBOL的向量
  difgene <- cluster_top$gene
  #将提取的SYMBOL向量转换为ENTREZID(这是Entrez gene数据库通用的ID)
  difgeneID <- bitr(difgene, 
                    fromType="SYMBOL", 
                    toType="ENTREZID", 
                    OrgDb="org.Ss.eg.db")
  difgoALL <- enrichGO(gene = difgeneID$ENTREZID, 
                       OrgDb = org.Ss.eg.db, 
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pvalueCutoff = 0.1, #默认是0.05
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,#默认是0.2
                       minGSSize = 10, #加不加没区别
                       maxGSSize = 500, #加不加没区别
                       readable = TRUE,
                       pool = FALSE)
  if (nrow(difgoALL) !=0){
    GOplot <- dotplot(difgoALL,font.size = 15) +
      scale_color_continuous(low="#7bc043", high="#ee4035")+
      scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
      theme_classic()+
      theme(axis.text.x = element_text(
        angle = 45,
        vjust = 0.5, hjust = 0.5
      ))
    ggsave(GOplot, file= paste("细胞群_mm_GO_第",b,"簇","_",i,".pdf",sep = ""), width=5.8, height=4)
    write.csv(difgoALL@result,file= paste("细胞群_mm_GO_第",b,"簇","_",i,".csv",sep = ""),row.names = F)
    # write.csv(cluster_top,file= paste("GOplot_",i,"_gene.csv",sep = ""),row.names = F)
  }
}

#########step15######################
#按照cluster 和 细胞种类分别做一个。
#热图，默认slot = 'scale.data',所以有很多基因找不到。因为前面只用了2000个高可变基因做中心化。
#所以，我前面把所有的基因都做了中心化。生信技能书有代码能解决这个问题。对所有数据进行中心化。
#对所有基因进行中心化
all.genes <- rownames(yaqun_1_Adipo_mm)
yaqun_1_Adipo_mm <- ScaleData(yaqun_1_Adipo_mm, features = all.genes)

yaqun_1_Adipo_mm@meta.data$yl <- yaqun_1_Adipo_mm@active.ident
#每个细胞类型，选出来50个细胞。
#提取所有细胞的细胞类型，形成一个向量
suoyou_zhushi <- yaqun_1_Adipo_mm@meta.data$yl
#细胞类型的向量
xibao_names <- levels(yaqun_1_Adipo_mm)
#每种细胞选50个，保存细胞在向量中的所在位置。
cluster50_all <- NULL
for (i in xibao_names) {
  cluster50_temp <-head(which(suoyou_zhushi== i),500)
  cluster50_all <- union(cluster50_all, cluster50_temp)
}
#画图
plot300 <- DoHeatmap(yaqun_1_Adipo_mm, 
                     features = top10$gene,
                     slot = "scale.data",
                     cells = cluster50_all,
                     group.colors = col_3,
                     draw.lines = F,
                     #hjust = 0.5, #文字居中。
                     raster = T #这个设置成F，图片就不会模糊
)#angle = 0
#这个梯度，是设置的热图红色和蓝色的过渡。
plot300 <- plot300 + scale_fill_gradientn(colors = c("#00aedb","white","#d11141"))
ggsave(plot300, file='plot300-热图10.pdf', width=8, height=8)

#########step16######################
#产热相关
markers.to.plot5 <- c("PPARG","LPL","ADD1","LEP","LIFR", 
                      "DGAT1", "FGFR2", "ACLY", "SCD", "ACSS2", 
                      "ANXA1", "PHLDB2", "ELOVL6", "ACACA", "ACSL1",
                      "LPIN1", "PLIN2", "ELL2", "ZBTB20", "WNT5B", 
                      "TGFB2", "EFEMP1")
len <- length(markers.to.plot5)
gene_cell_exp_all <- AverageExpression(yaqun_1_Adipo_mm,
                                       features = markers.to.plot5,
                                       group.by = 'yaqun',
                                       assays = "RNA",
                                       slot = 'data') 
EXP_ALL2 <- gene_cell_exp_all[["RNA"]]
#进行z-score
RS_F0_15_2 = t(scale(t(EXP_ALL2),center = TRUE, scale = TRUE))
#限定上限，使表达量大于2的等于2
RS_F0_15_2[RS_F0_15_2> 2] <- 2 
#限定下限，使表达量小于-2的等于-2
RS_F0_15_2[RS_F0_15_2< -2] <-  -2 
#热图
library(pheatmap)
bb <- pheatmap(RS_F0_15_2,
               cluster_cols=F, #不对列进行聚类
               cluster_rows=F,
               show_rownames = T, #不显示行名
               border_color = "white",
               cellwidth = 15,
               cellheight = 15,
               color = colorRampPalette(colors =  c("#372ac1","#3060e9","#378afb","#59acfd","#aad5ff","#dcebfb",
                                                    "#FFFFFF","#fde3e4","#ffb5ba","#fd6f6e","#f52a27","#e12112","#d01c05"))(10000))
ggsave(bb, file='plot热图-产热marker-合并之后.pdf', width=5, height=6)

#########step17######################
#产热相关
markers.to.plot6 <- c("ABCA6", "PDGFRA", "HOXC9", "DCLK1", "CD44",
                      "DGAT1", "DGAT2", "FGFR2", "TGFB2",
                      "LEP", "ADD1", "ACSS2", "PDGFA", "GATM", 
                      "FOXO1", "NRP1", "IRF4", "GLUL", "CD81",
                      "DHRS4", "PDK4", "ALDH1A1", "TMEM26", "PRDM16", 
                      "F3", "MAPK14", "ACE2", "ACSL1", "B4GALT1", 
                      "SGK2", "PPARG",  "VEGFA", "ABHD3", 
                      "IDH1", "ABHD5", "CD36", "LIPE", "PHLDB2")
len <- length(markers.to.plot6)
gene_cell_exp_all <- AverageExpression(yaqun_1_Adipo_mm,
                                       features = markers.to.plot6,
                                       group.by = 'yaqun',
                                       assays = "RNA",
                                       slot = 'data') 
EXP_ALL2 <- gene_cell_exp_all[["RNA"]]
#进行z-score
RS_F0_15_2 = t(scale(t(EXP_ALL2),center = TRUE, scale = TRUE))
#限定上限，使表达量大于2的等于2
RS_F0_15_2[RS_F0_15_2> 2] <- 2 
#限定下限，使表达量小于-2的等于-2
RS_F0_15_2[RS_F0_15_2< -2] <-  -2 
#热图
library(pheatmap)
bb <- pheatmap(RS_F0_15_2,
               cluster_cols=F, #不对列进行聚类
               cluster_rows=F,
               show_rownames = T, #不显示行名
               border_color = "white",
               cellwidth = 20,
               cellheight = 10,
               color = colorRampPalette(colors =  c("#372ac1","#3060e9","#378afb","#59acfd","#aad5ff","#dcebfb",
                                                    "#FFFFFF","#fde3e4","#ffb5ba","#fd6f6e","#f52a27","#e12112","#d01c05"))(10000))
ggsave(bb, file='plot热图-产热marker-合并之后2.pdf', width=5, height=10)

#########step18######################
#产热相关
markers.to.plot6 <- c("EBF2","ABCA6", "PDGFRA", "HOXC9", "DCLK1", "CD44",
                      "DGAT1", "DGAT2", "FGFR2", "TGFB2",
                      "LEP", "ADD1", "ACSS2", "PDGFA", "GATM", 
                      "FOXO1", "NRP1", "IRF4", "GLUL", "CD81",
                      "DHRS4", "PDK4", "ALDH1A1", "TMEM26", "PRDM16", 
                      "F3", "MAPK14", "ACE2", "ACSL1", "B4GALT1", 
                      "SGK2", "PPARG",  "VEGFA", "ABHD3", 
                      "IDH1", "ABHD5", "CD36", "LIPE", "PHLDB2")
len <- length(markers.to.plot6)
gene_cell_exp_all <- AverageExpression(yaqun_1_Adipo_mm,
                                       features = markers.to.plot6,
                                       group.by = 'yaqun',
                                       assays = "RNA",
                                       slot = 'data') 
EXP_ALL2 <- gene_cell_exp_all[["RNA"]]
#进行z-score
RS_F0_15_2 = t(scale(t(EXP_ALL2),center = TRUE, scale = TRUE))
#限定上限，使表达量大于2的等于2
RS_F0_15_2[RS_F0_15_2> 2] <- 2 
#限定下限，使表达量小于-2的等于-2
RS_F0_15_2[RS_F0_15_2< -2] <-  -2 
#热图
library(pheatmap)
bb <- pheatmap(RS_F0_15_2,
               cluster_cols=F, #不对列进行聚类
               cluster_rows=F,
               show_rownames = T, #不显示行名
               border_color = "white",
               cellwidth = 20,
               cellheight = 10,
               color = colorRampPalette(colors =  c("#372ac1","#3060e9","#378afb","#59acfd","#aad5ff","#dcebfb",
                                                    "#FFFFFF","#fde3e4","#ffb5ba","#fd6f6e","#f52a27","#e12112","#d01c05"))(10000))
ggsave(bb, file='plot热图-产热marker-合并之后3.pdf', width=5, height=10)

#########step19######################
markers.to.plot10 <-  c("ABCA6", "PDGFRA", "HOXC9", "DCLK1", "CD44",
                       "DGAT1", "DGAT2", "FGFR2", "TGFB2",
                       "LEP", "ADD1", "ACSS2", "PDGFA", "GATM", 
                       "FOXO1", "NRP1", "IRF4", "GLUL", "CD81",
                       "DHRS4", "PDK4", "ALDH1A1", "TMEM26", "PRDM16", 
                       "F3", "MAPK14", "ACE2", "ACSL1", "B4GALT1", 
                       "SGK2", "PPARG",  "VEGFA", "ABHD3", 
                       "IDH1", "ABHD5", "CD36", "LIPE", "PHLDB2")
#总体
dotplot <- DotPlot(yaqun_1_Adipo_mm, 
                   features = markers.to.plot10, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪亚群_整合之后.pdf', width=14, height=3)

markers.to.plot10 <-  c("IRF4")
#总体
dotplot <- DotPlot(yaqun_1_Adipo_mm, 
                   features = markers.to.plot10, 
                   cols = c("#00ff85","#e90052"),
                   dot.scale = 8) +
  RotatedAxis()
ggsave(dotplot, file='plot-点图marker_脂肪亚群_整合之后2.pdf', width=14, height=3)

#########step20######################
yaqun_1_Adipo_mm_sWAT <-  subset(yaqun_1_Adipo_mm, weizhi == "sWAT")
yaqun_1_Adipo_mm_sWAT <-  subset(yaqun_1_Adipo_mm, weizhi == "sWAT")

#依次计算DEGs--sWAT
b <- 0
for (i in levels(yaqun_1_Adipo_mm_sWAT)){
  b <- b + 1
  #wenjianm = paste("DEGs_in_cluster_",i,"_betweenWTMu",sep = "")
  DEGs_up <- FindMarkers(yaqun_1_Adipo_mm_sWAT, 
                         ident.1 = "KI", 
                         ident.2 = "WT",
                         group.by = "chuli", 
                         subset.ident = i,
                         min.pct = 0.3,
                         logfc.threshold = 0.5,
                         test.use = "wilcox",
                         only.pos = F
  ) %>%
    filter(p_val < 0.05 & avg_log2FC > 0) 
  if (nrow(DEGs_up) !=0){
    assign(paste("JJ_DEGs_",b,"_",i,"_UP",sep = ""), DEGs_up)
    write.csv(DEGs_up,file= paste("JJ_DEGs_",b," ",i,"_UP.csv",sep = "")) 
  }
  DEGs_down <- FindMarkers(yaqun_1_Adipo_mm_sWAT, 
                           ident.1 = "KI", 
                           ident.2 = "WT",
                           group.by = "chuli", 
                           subset.ident = i,
                           min.pct = 0.3,
                           logfc.threshold = 0.5, 
                           test.use = "wilcox",
                           only.pos = F
  ) %>%
    filter(p_val < 0.05 & avg_log2FC < 0)
  if (nrow(DEGs_down) !=0){
    assign(paste("sWAT_DEGs_",b,"_",i,"_DOWN",sep = ""), DEGs_down)
    write.csv(DEGs_down,file= paste("JJ_DEGs_",b," ",i,"_DOWN.csv",sep = "")) 
  }
}

#依次计算DEGs--sWAT
b <- 0
for (i in levels(yaqun_1_Adipo_mm_sWAT)){
  b <- b + 1
  #wenjianm = paste("DEGs_in_cluster_",i,"_betweenWTMu",sep = "")
  DEGs_up <- FindMarkers(yaqun_1_Adipo_mm_sWAT, 
                         ident.1 = "KI", 
                         ident.2 = "WT",
                         group.by = "chuli", 
                         subset.ident = i,
                         min.pct = 0.3,
                         logfc.threshold = 0.5, 
                         test.use = "wilcox",
                         only.pos = F
  ) 
  #这块稍微调整了一下，是因为2簇没有DEGs。这样就可以继续循环了。
  if (nrow(DEGs_up) !=0){
    DEGs_up <- filter(DEGs_up,p_val < 0.05 & avg_log2FC > 0) 
    assign(paste("SZ_DEGs_",b,"_",i,"_UP",sep = ""), DEGs_up)
    write.csv(DEGs_up,file= paste("SZ_DEGs_",b," ",i,"_UP.csv",sep = "")) 
  }
  
  DEGs_down <- FindMarkers(yaqun_1_Adipo_mm_sWAT, 
                           ident.1 = "KI", 
                           ident.2 = "WT",
                           group.by = "chuli", 
                           subset.ident = i,
                           min.pct = 0.3,
                           logfc.threshold = 0.5, 
                           test.use = "wilcox",
                           only.pos = F
  ) 
  if (nrow(DEGs_down) !=0){
    DEGs_down <- filter(DEGs_down,p_val < 0.05 & avg_log2FC < 0) 
    assign(paste("sWAT_DEGs_",b,"_",i,"_DOWN",sep = ""), DEGs_down)
    write.csv(DEGs_down,file= paste("SZ_DEGs_",b," ",i,"_DOWN.csv",sep = "")) 
  }
}

#########step21######################
#依次计算DEGs--swat
b <- 0
for (i in levels(yaqun_1_Adipo_mm_sWAT)){
  b <- b + 1
  #wenjianm = paste("DEGs_in_cluster_",i,"_betweenWTMu",sep = "")
  DEGs_up <- FindMarkers(yaqun_1_Adipo_mm_sWAT, 
                         ident.1 = "KI", 
                         ident.2 = "WT",
                         group.by = "chuli", 
                         subset.ident = i,
                         min.pct = 0.3,
                         logfc.threshold = 0.5, 
                         test.use = "wilcox",
                         only.pos = F
  ) %>%
    filter(p_val < 0.05 & avg_log2FC > 0) 
  #提取差异表达基因SYMBOL的向量
  DEGs_genes <- rownames(DEGs_up)
  #将提取的SYMBOL向量转换为ENTREZID(这是Entrez gene数据库通用的ID)
  difgeneID <- bitr(DEGs_genes, 
                    fromType="SYMBOL", 
                    toType="ENTREZID", 
                    OrgDb="org.Ss.eg.db")
  difgoALL <- enrichGO(gene = difgeneID$ENTREZID, 
                       OrgDb = org.Ss.eg.db, 
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pvalueCutoff = 0.05, #默认是0.05
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,#默认是0.2
                       minGSSize = 10, #加不加没区别
                       maxGSSize = 500, #加不加没区别
                       readable = TRUE,
                       pool = FALSE)
  if (nrow(difgoALL) !=0){
    GOplot <- dotplot(difgoALL,font.size = 15) +
      scale_color_continuous(low="#7bc043", high="#ee4035")+
      scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
      theme_classic()+
      theme(axis.text.x = element_text(
        angle = 45,
        vjust = 0.5, hjust = 0.5
      ))
    ggsave(GOplot, file= paste("swat_GOplot_",b," ",i,"_UP.pdf",sep = ""), width=5.8, height=4)
    write.csv(difgoALL@result,file= paste("swat_GOplot_",b," ",i,"_UP.csv",sep = ""),row.names = F) 
  }
  
  DEGs_down <- FindMarkers(yaqun_1_Adipo_mm_sWAT, 
                           ident.1 = "KI", 
                           ident.2 = "WT",
                           group.by = "chuli", 
                           subset.ident = i,
                           min.pct = 0.3,
                           logfc.threshold = 0.5, 
                           test.use = "wilcox",
                           only.pos = F
  ) %>%
    filter(p_val < 0.05 & avg_log2FC < 0)
  #提取差异表达基因SYMBOL的向量
  DEGs_genes <- rownames(DEGs_down)
  #将提取的SYMBOL向量转换为ENTREZID(这是Entrez gene数据库通用的ID)
  difgeneID <- bitr(DEGs_genes, 
                    fromType="SYMBOL", 
                    toType="ENTREZID", 
                    OrgDb="org.Ss.eg.db")
  difgoALL <- enrichGO(gene = difgeneID$ENTREZID, 
                       OrgDb = org.Ss.eg.db, 
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pvalueCutoff = 0.05, #默认是0.05
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,#默认是0.2
                       minGSSize = 10, #加不加没区别
                       maxGSSize = 500, #加不加没区别
                       readable = TRUE,
                       pool = FALSE)
  if (nrow(difgoALL) !=0){
    GOplot <- dotplot(difgoALL,font.size = 15) +
      scale_color_continuous(low="#7bc043", high="#ee4035")+
      scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
      theme_classic()+
      theme(axis.text.x = element_text(
        angle = 45,
        vjust = 0.5, hjust = 0.5
      ))
    ggsave(GOplot, file= paste("swat_GOplot_",b," ",i,"_DOWN.pdf",sep = ""), width=5.8, height=4)
    write.csv(difgoALL@result,file= paste("swat_GOplot_",b," ",i,"_DOWN.csv",sep = ""),row.names = F) 
  }
  #这个命令可以循环生成变量名
  #assign(wenjianm,marker)
  #这个命令可以循环生成文件名
  #write.csv(marker,file=paste("DEGs_in_cluster_",i,"_betweenWTMu.csv",sep = ""))
}

#依次计算DEGs--sWAT
#备注 SMCs骨骼肌，肾周KI里没有，所以要跳过。
b <- 0
for (i in levels(yaqun_1_Adipo_mm_sWAT)){
  b <- b + 1
  #wenjianm = paste("DEGs_in_cluster_",i,"_betweenWTMu",sep = "")
  DEGs_up <- FindMarkers(yaqun_1_Adipo_mm_sWAT, 
                         ident.1 = "KI", 
                         ident.2 = "WT",
                         group.by = "chuli", 
                         subset.ident = i,
                         min.pct = 0.3,
                         logfc.threshold = 0.5, 
                         test.use = "wilcox",
                         only.pos = F
  ) %>%
    filter(p_val < 0.05 & avg_log2FC > 0) 
  #提取差异表达基因SYMBOL的向量
  DEGs_genes <- rownames(DEGs_up)
  #将提取的SYMBOL向量转换为ENTREZID(这是Entrez gene数据库通用的ID)
  difgeneID <- bitr(DEGs_genes, 
                    fromType="SYMBOL", 
                    toType="ENTREZID", 
                    OrgDb="org.Ss.eg.db")
  difgoALL <- enrichGO(gene = difgeneID$ENTREZID, 
                       OrgDb = org.Ss.eg.db, 
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pvalueCutoff = 0.05, #默认是0.05
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,#默认是0.2
                       minGSSize = 10, #加不加没区别
                       maxGSSize = 500, #加不加没区别
                       readable = TRUE,
                       pool = FALSE)
  if (nrow(difgoALL) !=0){
    GOplot <- dotplot(difgoALL,font.size = 15) +
      scale_color_continuous(low="#7bc043", high="#ee4035")+
      scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
      theme_classic()+
      theme(axis.text.x = element_text(
        angle = 45,
        vjust = 0.5, hjust = 0.5
      ))
    ggsave(GOplot, file= paste("sWAT_GOplot_",b," ",i,"_UP.pdf",sep = ""), width=5.8, height=4)
    write.csv(difgoALL@result,file= paste("sWAT_GOplot_",b," ",i,"_UP.csv",sep = ""),row.names = F) 
  }
  
  DEGs_down <- FindMarkers(yaqun_1_Adipo_mm_sWAT, 
                           ident.1 = "KI", 
                           ident.2 = "WT",
                           group.by = "chuli", 
                           subset.ident = i,
                           min.pct = 0.3,
                           logfc.threshold = 0.5, 
                           test.use = "wilcox",
                           only.pos = F
  ) %>%
    filter(p_val < 0.05 & avg_log2FC < 0)
  #提取差异表达基因SYMBOL的向量
  DEGs_genes <- rownames(DEGs_down)
  #将提取的SYMBOL向量转换为ENTREZID(这是Entrez gene数据库通用的ID)
  difgeneID <- bitr(DEGs_genes, 
                    fromType="SYMBOL", 
                    toType="ENTREZID", 
                    OrgDb="org.Ss.eg.db")
  difgoALL <- enrichGO(gene = difgeneID$ENTREZID, 
                       OrgDb = org.Ss.eg.db, 
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pvalueCutoff = 0.05, #默认是0.05
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.2,#默认是0.2
                       minGSSize = 10, #加不加没区别
                       maxGSSize = 500, #加不加没区别
                       readable = TRUE,
                       pool = FALSE)
  if (nrow(difgoALL) !=0){
    GOplot <- dotplot(difgoALL,font.size = 15) +
      scale_color_continuous(low="#7bc043", high="#ee4035")+
      scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
      theme_classic()+
      theme(axis.text.x = element_text(
        angle = 45,
        vjust = 0.5, hjust = 0.5
      ))
    ggsave(GOplot, file= paste("sWAT_GOplot_",b," ",i,"_DOWN.pdf",sep = ""), width=5.8, height=4)
    write.csv(difgoALL@result,file= paste("sWAT_GOplot_",b," ",i,"_DOWN.csv",sep = ""),row.names = F) 
  }
  #这个命令可以循环生成变量名
  #assign(wenjianm,marker)
  #这个命令可以循环生成文件名
  #write.csv(marker,file=paste("DEGs_in_cluster_",i,"_betweenWTMu.csv",sep = ""))
}

#########step22######################
plot401 <- VlnPlot(yaqun_1_Adipo_mm,
                   features = c("GLUL","CD81","DHRS4","ACE2","ZFAND5",
                                "ACSL1","EBF2","PCK1","IRF4"
                   ),
                   idents = "C1",
                   cols = col_3,
                   group.by = "orig.ident", 
                   pt.size = 0, 
                   combine = T,
                   same.y.lims = F,
                   ncol = 1)
ggsave(plot401, file='plot401-小提琴图-adipo-C1.pdf', width=6, height=22.5)

plot401 <- VlnPlot(yaqun_1_Adipo_mm,
                   features = c("ACSL1","GRB10","ABHD5","ACE2","MCTP2",
                                "ZFAND5","EBF2","VEGFA","IDH1","DHRS4",
                                "THRB","F3","MAPK14","B4GALT1","IRF4"
                   ),
                   idents = "C4",
                   cols = col_3,
                   group.by = "orig.ident", 
                   pt.size = 0, 
                   combine = T,
                   same.y.lims = F,
                   ncol = 1)
ggsave(plot401, file='plot402-小提琴图-adipo-C4.pdf', width=6, height=37.5)

plot401 <- VlnPlot(yaqun_1_Adipo_mm,
                   features = c("WNT5B", "ADIG", "ITGB5", "ELOVL7", 
                                 "ARHGAP21", "LYPLAL1", "GPR180"
                   ),
                   idents = "C2",
                   cols = col_3,
                   group.by = "orig.ident", 
                   pt.size = 0, 
                   combine = T,
                   same.y.lims = F,
                   ncol = 1)
ggsave(plot401, file='plot403-小提琴图-adipo-C2.pdf', width=6, height=17.5)

plot401 <- VlnPlot(yaqun_1_Adipo_mm,
                   features = c("ACSL1", "ABHD5", "ADAMTS20", "ANGPTL4", "AUTS2",
                                "B4GALT1", "DHRS4", "FOXO1", "GRB10", "IDH1",
                                "MGLL", "PCK1", "ZFAND5"
                   ),
                   idents = "C3",
                   cols = col_3,
                   group.by = "orig.ident", 
                   pt.size = 0, 
                   combine = T,
                   same.y.lims = F,
                   ncol = 1)
ggsave(plot401, file='plot404-小提琴图-adipo-C3.pdf', width=6, height=32.5)



