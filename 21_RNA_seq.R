#一共需要两个文件，一个是countData，一个是colData，然后design是包含在colData中的condition中的。
#用右边的Impoer Dataset功能导入gene_count文件，然后heading选择yes，row names选择first colomu。
#自己创造一个表格为coldata，按照要求的格式，只需要有condition就行,下面是untreat treat。

#读入库
library("DESeq2")
library("ggplot2")

setwd("E:/BIIF/liutianxia/lilan_RNAseq")

coldata <- read.table("coldata.txt",header = TRUE)

gene_count <- read.csv(file = "gene_count.csv", header = TRUE,row.names = 1)

###过滤
#看看这行有多少大于10
he <- rowSums(gene_count > 0)
gene_count$he <- he
gene_count <- subset(gene_count, he > 2)
gene_count <- gene_count[,-9]

#创建矩阵dds
dds <- DESeqDataSetFromMatrix(countData = gene_count,
                              colData = coldata,
                              design = ~condition)

#下面这个命令可以查看组别，如何分组，以及level。默认是按照英文字母的先后顺序的。tr在un前面。
dds$condition 
#如果组别不满意，需要用下面的命令来修改对照组和处理组。修改完之后就别再点上一条命令了，否则又会修改回去。
dds$condition <- relevel(dds$condition, ref = "WT")
dds$condition 

dds <- DESeq(dds) ##一步到位，不需要太多步骤
res <- results(dds) ##得到结果
res #查看结果

#输出csv格式的文件
#write.csv(res,file = "res.csv")

#resdata是把原始数据和res数据 merge在一个矩阵里
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)

#indexOf(str,str2)，str是被搜索的字符串，str2是搜索的字符串，返回搜索的字符串在被搜索的字符串中第一次出现的位置，找不到返回0
indexOf = function(str,str2){
  cd=nchar(str);
  cd2=nchar(str2);
  if(cd==0||cd2==0){
    return(0);
  }
  for(i in 1:cd){
    t=substr(str,i,i);
    for(j in 1:cd2){
      if(t==substr(str2,j,j)&&j==1){
        if(cd2==1){
          return(i);
        }else{
          c=TRUE;
          for(k in 1:(cd2-1)){
            if(substr(str,i+k,i+k)!=substr(str2,j+k,j+k)){
              c=FALSE;
              break;
            }
          }
          if(c==TRUE){
            return(i);
          }
        }
      }else{
        break;
      }
    }
  }
  return(0);
}


#截取|竖线之后的字符
for (i in 1:length(resdata$Row.names)){
  resdata$Row.names[i] <- substring(resdata$Row.names[i],indexOf(resdata$Row.names[i],"|")+1,50)
}

#大概看看结果,这个能看多少上调，多少下调
summary(res,0.01)
#这个命令可以看有多少P值小于0.01的基因
sum(res$padj < 0.001, na.rm=TRUE)

# 取padj小于0.05的数据，得到行
table(res$padj<0.001) 

# 按照padj的大小将resdata重新排列，并且输出文件
resdata <- resdata[order(resdata$padj),] 
write.csv(resdata,file = "resdata.csv")

# 筛选（获取padj小于0.01，并且L2FC大于0.4或者小于-0.4的差异表达基因） 并赋值给resdataDEG
#resdataDEG <- subset(resdata,padj < 0.05 & (log2FoldChange >1 | log2FoldChange < -1))
#write.csv(resdataDEG,file = "resdataDEG.csv")

resdataDEG <- subset(resdata,pvalue < 0.05)
write.csv(resdataDEG,file = "resdataDEG_pvalue.csv")

#MAplot图 p值小于0.001的为蓝色点
#plotMA1 <- plotMA(res, alpha = 0.05, ylim=c(-3,3))
#ggsave(plotMA1, file = "plot1-plotMA.pdf", width=6, height=6 )
#有差异时，it is important to set blind=FALSE for downstream analysis
#做PCA需要对表达量标准化，大于30样本用vst函数，因为快。低于30样本建议使用rlog
guiyihua <- rlog(dds, blind= FALSE)

#PCA图,方法1.
plotPCA(guiyihua, intgroup="condition")

#PCA图，方法4
library(ggplot2)
library(ggrepel)
pcaData <- plotPCA(guiyihua, intgroup="condition", returnData = TRUE)
plot2 <- ggplot(data=pcaData, aes(x=PC1, y=PC2, color=condition)) +
  geom_text_repel(aes(PC1,PC2,label = name),show.legend = FALSE, box.padding = 0.5) + 
  #这里show.legend或者show_guide=F可以吧图例的a去掉。box.padding = 0.5是增加文字和圆圈的距离
  geom_point(alpha=1, size = 3) +    #alpha=0.8为透明度
  xlab(paste0("PC1: 41% variance")) +  #这俩数是从方法1获得的
  ylab(paste0("PC2: 25% variance")) +  #这俩数是从方法1获得的
  scale_colour_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33" )) + 
  scale_alpha_discrete(range = c(0.5,1)) +
  coord_fixed() +  
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(plot2, file = "plot2-PCA.pdf", width=6, height=6)


#热图
library(pheatmap)
#输入RPKM矩阵
resdata <- read.csv(file = "resdata.csv", header = TRUE)
#输入想要的基因
adipogene <- read.csv(file = "adipogene_ncbi.csv", header = TRUE)
#提取
tiqu <- merge(adipogene, resdata,by= "Row.names",sort=FALSE)
#行名
row.names(tiqu) <- tiqu$Row.names
#保留相应的列
tiqu <- tiqu[,9:16]
#进行z-score，,center = TRUE, scale = TRUE 这个可加可不加
tiqu2=t(scale(t(tiqu),center = TRUE, scale = TRUE))
tiqu2[tiqu2>2]=2 #限定上限，使表达量大于2的等于2
tiqu2[tiqu2< -2]= -2 #限定下限，使表达量小于-2的等于-2
#热图

aaa <- pheatmap(tiqu2,
                cluster_cols=F, #不对列进行聚类
                cluster_rows=F,
                show_rownames = T, #不显示行名
                color = colorRampPalette(colors = c("#1e7574","#1ea09f","#74d1d0","#84eeed",
                                                    "#ffc100","#ff9a00","#ff7400","#ff4d00"))(10000))
ggsave(aaa, file='adipo-heat.pdf', width=5, height=8)


aaa <- pheatmap(tiqu2,
                cluster_cols=F, #不对列进行聚类
                show_rownames = T, #不显示行名
                color = colorRampPalette(colors = c("#1e7574","#1ea09f","#74d1d0","#84eeed",
                                                    "#ffc100","#ff9a00","#ff7400","#ff4d00"))(10000))
ggsave(aaa, file='adipo-heat2.pdf', width=5, height=8)

  
  
