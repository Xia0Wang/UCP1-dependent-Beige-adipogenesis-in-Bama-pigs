
#options(connectionObserver = NULL)
#library(topGO)
#library(pathview)
#install.packages("readxl")
library("readxl")
library(clusterProfiler)
library(org.Ss.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(stringr)
library("ggplot2")
library("GOSemSim")
setwd("E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes")

#使用GOsemsim证明两个GO集合的相关性强。来证明KI猪和藏猪的相似性
#参考张勇的文章qu-et-al-2021-the-evolution-of-ancestral-and-species-specific-adaptations-in-snowfinches-at-the-qinghai-tibet-plateau
#随机，然后线性拟合。然后看看实验组在这些随机的什么位置，为了证明实验组的相似性数值大。

#数字备注：
#DEGs数： 丽兰：上调1941，下调1902
#DEGs数： 田侠：上调676，下调599
#上调基因的GO数：丽兰，176
#下调基因的GO数：丽兰，508
#上调基因的GO数：田侠，215
#下调基因的GO数：田侠，136

###############第一步，随机生成基因名###########
#备注，使用了全部的基因，包括LOC的。来进行随机。
Genes <- read.csv("gene_count_ncbi.csv")
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
for (i in 1:length(Genes$gene_id)){
  Genes$gene_id[i] <- substring(Genes$gene_id[i],indexOf(Genes$gene_id[i],"|")+1,50)
}

setwd("E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes/Lilan_up_random_genes")
for (i in c(1:100)){
  lilan_up <- sample(Genes$gene_id,size = 1941)
  lilan_up <- as.data.frame(lilan_up)
  write.csv(lilan_up, paste("Lilan_up_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes/Tianxia_up_random_genes")
for (i in c(1:100)){
  Tianxia_up <- sample(Genes$gene_id,size = 676)
  Tianxia_up <- as.data.frame(Tianxia_up)
  write.csv(Tianxia_up, paste("Tianxia_up_random_genes_",i,".csv",sep = ""),
            row.names = F)}
#Tianxia_up <- sample(Genes$gene_id,size = 676)
#write.csv(Tianxia_up, "1.csv")

setwd("E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes/Lilan_down_random_genes")
for (i in c(1:100)){
  Lilan_down <- sample(Genes$gene_id,size = 1902)
  Lilan_down <- as.data.frame(Lilan_down)
  write.csv(Lilan_down, paste("Lilan_down_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes/Tianxia_down_random_genes")
for (i in c(1:100)){
  Tianxia_down <- sample(Genes$gene_id,size = 599)
  Tianxia_down <- as.data.frame(Tianxia_down)
  write.csv(Tianxia_down, paste("Tianxia_down_random_genes_",i,".csv",sep = ""),
            row.names = F)}
#Tianxia_down <- sample(Genes$gene_id,size = 599)
#write.csv(Tianxia_up, "2.csv")

###############第二步，手动到网站里去计算GO###############
#https://metascape.org/gp/index.html#/main/step1
#一共400组基因。

############第三部，计算,上调##########
setwd("E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes")
HsGO <- godata('org.Hs.eg.db', ont="BP")
path1 <- "E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes/Lilan_up_GO/"
path2 <- "E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes/Tianxia_up_GO/"
file_names_lilan <- list.files(path1)
file_names_tianxia <- list.files(path2)
#ss_score_max <- vector()
#ss_score_avg <- vector()
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:100)){
  GO1 <- read_excel(paste0(path1,file_names_lilan[k]), sheet = 2)
  #截取|竖线之后的字符
  for (i in 1:length(GO1$GroupID)){
    GO1$GroupID[i] <- substring(GO1$GroupID[i],indexOf(GO1$GroupID[i],"_")+1,50)
  }
  GO1 <- subset(GO1, GroupID == "Member" & Category == "GO Biological Processes")
  
  GO2 <- read_excel(paste0(path2,file_names_tianxia[k]), sheet = 2)
  #截取|竖线之后的字符
  for (i in 1:length(GO2$GroupID)){
    GO2$GroupID[i] <- substring(GO2$GroupID[i],indexOf(GO2$GroupID[i],"_")+1,50)
  }
  GO2 <- subset(GO2, GroupID == "Member" & Category == "GO Biological Processes")
  
  go1 <- GO1$Term
  go2 <- GO2$Term
  #One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
  #score1 <- mgoSim(go1, go2, semData=ssGO, measure="Wang", combine="max")
  #ss_score_max[k] <- score1
  #score2 <- mgoSim(go1, go2, semData=ssGO, measure="Wang", combine="avg")
  #ss_score_avg[k] <- score2
  score3 <- mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
  ss_score_rcmax[k] <- score3
  score4 <- mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
  ss_score_BMA[k] <- score4
}
#write.csv(ss_score_max[1:100],"ss_score_max_random_gene_UP.csv")
#write.csv(ss_score_avg[1:100],"ss_score_avg_random_gene_UP.csv")
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_random_gene_UP.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_random_gene_UP.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_random_gene_UP.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#fba6ab")+
  geom_vline(xintercept  = 0.718, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_random_gene_UP.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_random_gene_UP.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#fba6ab")+
  geom_vline(xintercept  = 0.679, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_random_gene_UP.pdf', width=6, height=3)

############第4部，计算,下调##########
ssGO <- godata('org.Ss.eg.db', ont="BP")
path1 <- "E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes/Lilan_down_GO/"
path2 <- "E:/BIIF/liutianxia/lilan_RNAseq/ss值随机_genes/Tianxia_down_GO/"
file_names_lilan <- list.files(path1)
file_names_tianxia <- list.files(path2)
#ss_score_max <- vector()
#ss_score_avg <- vector()
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:100)){
  GO1 <- read_excel(paste0(path1,file_names_lilan[k]), sheet = 2)
  #截取|竖线之后的字符
  for (i in 1:length(GO1$GroupID)){
    GO1$GroupID[i] <- substring(GO1$GroupID[i],indexOf(GO1$GroupID[i],"_")+1,50)
  }
  GO1 <- subset(GO1, GroupID == "Member" & Category == "GO Biological Processes")
  
  GO2 <- read_excel(paste0(path2,file_names_tianxia[k]), sheet = 2)
  #截取|竖线之后的字符
  for (i in 1:length(GO2$GroupID)){
    GO2$GroupID[i] <- substring(GO2$GroupID[i],indexOf(GO2$GroupID[i],"_")+1,50)
  }
  GO2 <- subset(GO2, GroupID == "Member" & Category == "GO Biological Processes")
  
  go1 <- GO1$Term
  go2 <- GO2$Term
  #One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
  #score1 <- mgoSim(go1, go2, semData=ssGO, measure="Wang", combine="max")
  #ss_score_max[k] <- score1
  #score2 <- mgoSim(go1, go2, semData=ssGO, measure="Wang", combine="avg")
  #ss_score_avg[k] <- score2
  score3 <- mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
  ss_score_rcmax[k] <- score3
  score4 <- mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
  ss_score_BMA[k] <- score4
}
#write.csv(ss_score_max[1:100],"ss_score_max_random_gene_down.csv")
#write.csv(ss_score_avg[1:100],"ss_score_avg_random_gene_down.csv")
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_random_gene_down.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_random_gene_down.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_random_gene_down.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#6bd2db")+
  geom_vline(xintercept  = 0.703, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_random_gene_down.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_random_gene_down.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#6bd2db")+
  geom_vline(xintercept  = 0.482, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_random_gene_down.pdf', width=6, height=3)


##########100次循环，从GO的随机选择出发###########
#上调基因的GO数：丽兰，176
#下调基因的GO数：丽兰，508
#上调基因的GO数：田侠，215
#下调基因的GO数：田侠，136
#读取总共的BP。
AllGO_BP <- unique(HsGO@geneAnno$GO)

########上调的#############
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:500)){
  Lilan_up <- sample(AllGO_BP,size = 176)
  Tianxia_up <- sample(AllGO_BP,size = 215)
  go1 <- Lilan_up
  go2 <- Tianxia_up
  #One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
  #score1 <- mgoSim(go1, go2, semData=ssGO, measure="Wang", combine="max")
  #ss_score_max[k] <- score1
  #score2 <- mgoSim(go1, go2, semData=ssGO, measure="Wang", combine="avg")
  #ss_score_avg[k] <- score2
  score3 <- mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
  ss_score_rcmax[k] <- score3
  score4 <- mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
  ss_score_BMA[k] <- score4
}
#write.csv(ss_score_max[1:100],"ss_score_max_random_gene_UP.csv")
#write.csv(ss_score_avg[1:100],"ss_score_avg_random_gene_UP.csv")
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_random_GOBP_UP.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_random_GOBP_UP.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_random_GOBP_UP.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#e9a0b9")+
  geom_vline(xintercept  = 0.718, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_random_GOBP_UP.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_random_GOBP_UP.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#e9a0b9")+
  geom_vline(xintercept  = 0.679, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_random_GOBP_UP.pdf', width=6, height=3)

########下调的#############
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:500)){
  Lilan_down <- sample(AllGO_BP,size = 508)
  Tianxia_down <- sample(AllGO_BP,size = 136)
  go1 <- Lilan_down
  go2 <- Tianxia_down
  #One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
  #score1 <- mgoSim(go1, go2, semData=ssGO, measure="Wang", combine="max")
  #ss_score_max[k] <- score1
  #score2 <- mgoSim(go1, go2, semData=ssGO, measure="Wang", combine="avg")
  #ss_score_avg[k] <- score2
  score3 <- mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
  ss_score_rcmax[k] <- score3
  score4 <- mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
  ss_score_BMA[k] <- score4
}
#write.csv(ss_score_max[1:100],"ss_score_max_random_gene_down.csv")
#write.csv(ss_score_avg[1:100],"ss_score_avg_random_gene_down.csv")
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_random_GOBP_down.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_random_GOBP_down.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_random_GOBP_down.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#72C8E5")+
  geom_vline(xintercept  = 0.703, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_random_GOBP_down.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_random_GOBP_down.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#72C8E5")+
  geom_vline(xintercept  = 0.482, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_random_GOBP_down.pdf', width=6, height=3)








