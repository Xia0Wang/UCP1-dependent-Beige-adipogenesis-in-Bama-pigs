library(ggplot2)
library(Seurat)
library(velocyto.R)
library(SeuratWrappers) 
library(hdf5r)

setwd("E:/BIIF/liutianxia/Adipocyte_scRNA_seq_8/3_结果/28_RNA速率 - 脂肪_整合之后")
getwd()

#########step1######################
#ldat <- read.loom.matrices(file = "~/possorted_genome_bam_LYJXT.loom")
ldat1 <- ReadVelocity(file = "E:/BIIF/liutianxia/Adipocyte_scRNA_seq_8/3_结果/26_RNA速率 - 脂肪/6609_JJ.loom")
#修改ldat1中的细胞名字。
colnames(ldat1$spliced) <- gsub("possorted_genome_bam_LYJXT:","6609_JJ_",colnames(ldat1$spliced))
colnames(ldat1$spliced) <- paste(colnames(ldat1$spliced),"-1",sep = "")
colnames(ldat1$unspliced) <- colnames(ldat1$spliced)
colnames(ldat1$ambiguous) <- colnames(ldat1$spliced)
#转换
bm1 <- as.Seurat(x = ldat1)

ldat2 <- ReadVelocity(file = "E:/BIIF/liutianxia/Adipocyte_scRNA_seq_8/3_结果/26_RNA速率 - 脂肪/6609_SZ.loom")
#修改ldat2中的细胞名字。
colnames(ldat2$spliced) <- gsub("possorted_genome_bam_V1E3I:","6609_SZ_",colnames(ldat2$spliced))
colnames(ldat2$spliced) <- paste(colnames(ldat2$spliced),"-1",sep = "")
colnames(ldat2$unspliced) <- colnames(ldat2$spliced)
colnames(ldat2$ambiguous) <- colnames(ldat2$spliced)
#转换
bm2 <- as.Seurat(x = ldat2)

ldat3 <- ReadVelocity(file = "E:/BIIF/liutianxia/Adipocyte_scRNA_seq_8/3_结果/26_RNA速率 - 脂肪/7103_JJ.loom")
#修改ldat3中的细胞名字。
colnames(ldat3$spliced) <- gsub("possorted_genome_bam_X2464:","7103_JJ_",colnames(ldat3$spliced))
colnames(ldat3$spliced) <- paste(colnames(ldat3$spliced),"-1",sep = "")
colnames(ldat3$unspliced) <- colnames(ldat3$spliced)
colnames(ldat3$ambiguous) <- colnames(ldat3$spliced)
#转换
bm3 <- as.Seurat(x = ldat3)

ldat4 <- ReadVelocity(file = "E:/BIIF/liutianxia/Adipocyte_scRNA_seq_8/3_结果/26_RNA速率 - 脂肪/7103_SZ.loom")
#修改ldat4中的细胞名字。
colnames(ldat4$spliced) <- gsub("possorted_genome_bam_8YF5V:","7103_SZ_",colnames(ldat4$spliced))
colnames(ldat4$spliced) <- paste(colnames(ldat4$spliced),"-1",sep = "")
colnames(ldat4$unspliced) <- colnames(ldat4$spliced)
colnames(ldat4$ambiguous) <- colnames(ldat4$spliced)
#转换
bm4 <- as.Seurat(x = ldat4)

#########step2######################
#4个样品分别做
yaqun_1_Adipo_tem <- subset(yaqun_1_Adipo_mm, orig.ident == "6609_JJ")
#选取
bm_all_no1 <- bm1
#取细胞的子集
maichi <- match(row.names(yaqun_1_Adipo_tem@meta.data),row.names(bm_all_no1@meta.data))
bm_all_no1@meta.data$mchi <- "A"
bm_all_no1@meta.data$mchi[maichi] <- "B"
bm_all_no1_2 <- subset(bm_all_no1, mchi == "B")
#加一个assay，seurat包带的一个标准化的功能
bm_all_no1_2 <- SCTransform(object = bm_all_no1_2, assay = "spliced")
bm_all_no1_2 <- RunPCA(object = bm_all_no1_2, verbose = FALSE)
bm_all_no1_2 <- FindNeighbors(object = bm_all_no1_2, dims = 1:20)
bm_all_no1_2 <- FindClusters(object = bm_all_no1_2)
bm_all_no1_2 <- RunUMAP(object = bm_all_no1_2, dims = 1:20)

#向bm_all_no1_2中填入UMAP信息。
int.embed <- Embeddings(yaqun_1_Adipo_tem, reduction = "umap")
bm_all_no1_2@reductions[["umap"]]@cell.embeddings <- int.embed
#向bm_all_no1_2中填入cluster信息。
maichi2 <- match(row.names(bm_all_no1_2@meta.data),row.names(yaqun_1_Adipo_tem@meta.data))
bm_all_no1_2@meta.data$yaqun <- yaqun_1_Adipo_tem@meta.data$yaqun[maichi2]
#把yaqun定为 active.ident
Idents(bm_all_no1_2) <- "yaqun"
#速率
bm_all_no1_2 <- RunVelocity(object = bm_all_no1_2, deltaT = 1, kCells = 25, fit.quantile = 0.02)

#画图
#ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm_all_no1_2)))
col_3 <- DiscretePalette(36, palette = "polychrome")[c(4:36)]
ident.colors <- col_3
names(x = ident.colors) <- levels(x = bm_all_no1_2)
cell.colors <- ident.colors[Idents(object = bm_all_no1_2)]
names(x = cell.colors) <- colnames(x = bm_all_no1_2)

pdf("plot_velocity_1_6609_JJ_2.pdf",width = 5,height = 5)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm_all_no1_2, reduction = "umap"), 
                               vel = Tool(object = bm_all_no1_2, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.9), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0)
dev.off()

yaqun_1_Adipo_tem <- subset(yaqun_1_Adipo_mm, orig.ident == "6609_SZ")
#选取
bm_all_no2 <- bm2
#取细胞的子集
maichi <- match(row.names(yaqun_1_Adipo_tem@meta.data),row.names(bm_all_no2@meta.data))
bm_all_no2@meta.data$mchi <- "A"
bm_all_no2@meta.data$mchi[maichi] <- "B"
bm_all_no2_2 <- subset(bm_all_no2, mchi == "B")
#加一个assay，seurat包带的一个标准化的功能
bm_all_no2_2 <- SCTransform(object = bm_all_no2_2, assay = "spliced")
bm_all_no2_2 <- RunPCA(object = bm_all_no2_2, verbose = FALSE)
bm_all_no2_2 <- FindNeighbors(object = bm_all_no2_2, dims = 1:20)
bm_all_no2_2 <- FindClusters(object = bm_all_no2_2)
bm_all_no2_2 <- RunUMAP(object = bm_all_no2_2, dims = 1:20)

#向bm_all_no2_2中填入UMAP信息。
int.embed <- Embeddings(yaqun_1_Adipo_tem, reduction = "umap")
bm_all_no2_2@reductions[["umap"]]@cell.embeddings <- int.embed
#向bm_all_no2_2中填入cluster信息。
maichi2 <- match(row.names(bm_all_no2_2@meta.data),row.names(yaqun_1_Adipo_tem@meta.data))
bm_all_no2_2@meta.data$yaqun <- yaqun_1_Adipo_tem@meta.data$yaqun[maichi2]
#把yaqun定为 active.ident
Idents(bm_all_no2_2) <- "yaqun"
#速率
bm_all_no2_2 <- RunVelocity(object = bm_all_no2_2, deltaT = 1, kCells = 25, fit.quantile = 0.02)

#画图
#ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm_all_no2_2)))
col_3 <- DiscretePalette(36, palette = "polychrome")[c(4:36)]
ident.colors <- col_3
names(x = ident.colors) <- levels(x = bm_all_no2_2)
cell.colors <- ident.colors[Idents(object = bm_all_no2_2)]
names(x = cell.colors) <- colnames(x = bm_all_no2_2)

pdf("plot_velocity_2_6609_SZ_2.pdf",width = 5,height = 5)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm_all_no2_2, reduction = "umap"), 
                               vel = Tool(object = bm_all_no2_2, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.9), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0)
dev.off()

yaqun_1_Adipo_tem <- subset(yaqun_1_Adipo_mm, orig.ident == "7103_JJ")
#选取
bm_all_no3 <- bm3
#取细胞的子集
maichi <- match(row.names(yaqun_1_Adipo_tem@meta.data),row.names(bm_all_no3@meta.data))
bm_all_no3@meta.data$mchi <- "A"
bm_all_no3@meta.data$mchi[maichi] <- "B"
bm_all_no3_2 <- subset(bm_all_no3, mchi == "B")
#加一个assay，seurat包带的一个标准化的功能
bm_all_no3_2 <- SCTransform(object = bm_all_no3_2, assay = "spliced")
bm_all_no3_2 <- RunPCA(object = bm_all_no3_2, verbose = FALSE)
bm_all_no3_2 <- FindNeighbors(object = bm_all_no3_2, dims = 1:20)
bm_all_no3_2 <- FindClusters(object = bm_all_no3_2)
bm_all_no3_2 <- RunUMAP(object = bm_all_no3_2, dims = 1:20)

#向bm_all_no3_2中填入UMAP信息。
int.embed <- Embeddings(yaqun_1_Adipo_tem, reduction = "umap")
bm_all_no3_2@reductions[["umap"]]@cell.embeddings <- int.embed
#向bm_all_no3_2中填入cluster信息。
maichi2 <- match(row.names(bm_all_no3_2@meta.data),row.names(yaqun_1_Adipo_tem@meta.data))
bm_all_no3_2@meta.data$yaqun <- yaqun_1_Adipo_tem@meta.data$yaqun[maichi2]
#把yaqun定为 active.ident
Idents(bm_all_no3_2) <- "yaqun"
#速率
bm_all_no3_2 <- RunVelocity(object = bm_all_no3_2, deltaT = 1, kCells = 25, fit.quantile = 0.02)

#画图
#ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm_all_no3_2)))
col_3 <- DiscretePalette(36, palette = "polychrome")[c(4:36)]
ident.colors <- col_3
names(x = ident.colors) <- levels(x = bm_all_no3_2)
cell.colors <- ident.colors[Idents(object = bm_all_no3_2)]
names(x = cell.colors) <- colnames(x = bm_all_no3_2)

pdf("plot_velocity_3_7103_JJ_2.pdf",width = 5,height = 5)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm_all_no3_2, reduction = "umap"), 
                               vel = Tool(object = bm_all_no3_2, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.9), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0)
dev.off()

yaqun_1_Adipo_tem <- subset(yaqun_1_Adipo_mm, orig.ident == "7103_SZ")
#选取
bm_all_no4 <- bm4
#取细胞的子集
maichi <- match(row.names(yaqun_1_Adipo_tem@meta.data),row.names(bm_all_no4@meta.data))
bm_all_no4@meta.data$mchi <- "A"
bm_all_no4@meta.data$mchi[maichi] <- "B"
bm_all_no4_2 <- subset(bm_all_no4, mchi == "B")
#加一个assay，seurat包带的一个标准化的功能
bm_all_no4_2 <- SCTransform(object = bm_all_no4_2, assay = "spliced")
bm_all_no4_2 <- RunPCA(object = bm_all_no4_2, verbose = FALSE)
bm_all_no4_2 <- FindNeighbors(object = bm_all_no4_2, dims = 1:20)
bm_all_no4_2 <- FindClusters(object = bm_all_no4_2)
bm_all_no4_2 <- RunUMAP(object = bm_all_no4_2, dims = 1:20)

#向bm_all_no4_2中填入UMAP信息。
int.embed <- Embeddings(yaqun_1_Adipo_tem, reduction = "umap")
bm_all_no4_2@reductions[["umap"]]@cell.embeddings <- int.embed
#向bm_all_no4_2中填入cluster信息。
maichi2 <- match(row.names(bm_all_no4_2@meta.data),row.names(yaqun_1_Adipo_tem@meta.data))
bm_all_no4_2@meta.data$yaqun <- yaqun_1_Adipo_tem@meta.data$yaqun[maichi2]
#把yaqun定为 active.ident
Idents(bm_all_no4_2) <- "yaqun"
#速率
bm_all_no4_2 <- RunVelocity(object = bm_all_no4_2, deltaT = 1, kCells = 25, fit.quantile = 0.02)

#画图
#ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm_all_no4_2)))
col_3 <- DiscretePalette(36, palette = "polychrome")[c(4:36)]
ident.colors <- col_3
names(x = ident.colors) <- levels(x = bm_all_no4_2)
cell.colors <- ident.colors[Idents(object = bm_all_no4_2)]
names(x = cell.colors) <- colnames(x = bm_all_no4_2)

pdf("plot_velocity_4_7103_SZ_2.pdf",width = 5,height = 5)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm_all_no4_2, reduction = "umap"), 
                               vel = Tool(object = bm_all_no4_2, slot = "RunVelocity"), 
                               n = 200, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.9), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0)
dev.off()




