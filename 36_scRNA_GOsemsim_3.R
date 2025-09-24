
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
setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因")

#使用GOsemsim证明两个GO集合的相关性强。来证明KI猪和藏猪的相似性
#参考张勇的文章qu-et-al-2021-the-evolution-of-ancestral-and-species-specific-adaptations-in-snowfinches-at-the-qinghai-tibet-plateau
#随机，然后线性拟合。然后看看实验组在这些随机的什么位置，为了证明实验组的相似性数值大。

#数字备注：
#整个adipo
#DEGs数： 嘉莉：上调23，下调28
#DEGs数： 田侠：上调248，下调253
#TGA
#DEGs数： 嘉莉：上调24，下调232
#DEGs数： 田侠：上调29，下调277

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

#去掉LOC基因
Genes$id2 <- substring(Genes$gene_id,1,3)
Genes <- subset(Genes, Genes$id2 != "LOC")

#整个adipo
#DEGs数： 嘉莉：上调23，下调28
#DEGs数： 田侠：上调248，下调253
#TGA
#DEGs数： 嘉莉：上调24，下调29
#DEGs数： 田侠：上调232，下调277

dir.create("Adipo_jiali_up")
dir.create("Adipo_jiali_down")
dir.create("Adipo_tianxia_up")
dir.create("Adipo_tianxia_down")
dir.create("TGA_jiali_up")
dir.create("TGA_jiali_down")
dir.create("TGA_tianxia_up")
dir.create("TGA_tianxia_down")

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Adipo_jiali_up")
for (i in c(1:110)){
  RandomGenes <- sample(Genes$gene_id,size = 23)
  RandomGenes <- as.data.frame(RandomGenes)
  write.csv(RandomGenes, paste("Adipo_jiali_up_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Adipo_jiali_down")
for (i in c(1:110)){
  RandomGenes <- sample(Genes$gene_id,size = 28)
  RandomGenes <- as.data.frame(RandomGenes)
  write.csv(RandomGenes, paste("Adipo_jiali_down_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Adipo_tianxia_up")
for (i in c(1:110)){
  RandomGenes <- sample(Genes$gene_id,size = 248)
  RandomGenes <- as.data.frame(RandomGenes)
  write.csv(RandomGenes, paste("Adipo_tianxia_up_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Adipo_tianxia_down")
for (i in c(1:110)){
  RandomGenes <- sample(Genes$gene_id,size = 253)
  RandomGenes <- as.data.frame(RandomGenes)
  write.csv(RandomGenes, paste("Adipo_tianxia_down_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/随机基因_1/TGA_jiali_up")
for (i in c(1:112)){
  RandomGenes <- sample(Genes$gene_id,size = 24)
  RandomGenes <- as.data.frame(RandomGenes)
  write.csv(RandomGenes, paste("TGA_jiali_up_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/TGA_jiali_down")
for (i in c(1:110)){
  RandomGenes <- sample(Genes$gene_id,size = 29)
  RandomGenes <- as.data.frame(RandomGenes)
  write.csv(RandomGenes, paste("TGA_jiali_down_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/TGA_tianxia_up")
for (i in c(1:110)){
  RandomGenes <- sample(Genes$gene_id,size = 232)
  RandomGenes <- as.data.frame(RandomGenes)
  write.csv(RandomGenes, paste("TGA_tianxia_up_random_genes_",i,".csv",sep = ""),
            row.names = F)}

setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/TGA_tianxia_down")
for (i in c(1:110)){
  RandomGenes <- sample(Genes$gene_id,size = 277)
  RandomGenes <- as.data.frame(RandomGenes)
  write.csv(RandomGenes, paste("TGA_tianxia_down_random_genes_",i,".csv",sep = ""),
            row.names = F)}

###############第二步，手动到网站里去计算GO###############
#https://metascape.org/gp/index.html#/main/step1
#一共800组基因。


########删除没GO上的######
path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/Adipo_jiali_up/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/Adipo_tianxia_up/"
file_names_Jiali <- list.files(path1)
file_names_tianxia <- list.files(path2)
for (k in c(1:100)){
  nnn <- length(excel_sheets(paste0(path1,file_names_Jiali[k])))
  if (nnn == 1){
    file.remove(paste0(path1,file_names_Jiali[k]))
  }
}

path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/Adipo_jiali_down/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/Adipo_tianxia_down/"
file_names_Jiali <- list.files(path1)
file_names_tianxia <- list.files(path2)
for (k in c(1:100)){
  nnn <- length(excel_sheets(paste0(path1,file_names_Jiali[k])))
  if (nnn == 1){
    file.remove(paste0(path1,file_names_Jiali[k]))
  }
}

path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/TGA_jiali_down/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/TGA_tianxia_down/"
file_names_Jiali <- list.files(path1)
file_names_tianxia <- list.files(path2)
for (k in c(1:100)){
  nnn <- length(excel_sheets(paste0(path1,file_names_Jiali[k])))
  if (nnn == 1){
    file.remove(paste0(path1,file_names_Jiali[k]))
  }
}

path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/TGA_jiali_up/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/TGA_tianxia_up/"
file_names_Jiali <- list.files(path1)
file_names_tianxia <- list.files(path2)
for (k in c(1:100)){
  nnn <- length(excel_sheets(paste0(path1,file_names_Jiali[k])))
  if (nnn == 1){
    file.remove(paste0(path1,file_names_Jiali[k]))
  }
}



############第三部，计算,上调##########
setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term")
HsGO <- godata('org.Hs.eg.db', ont="BP")
path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/Adipo_jiali_up/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/Adipo_tianxia_up/"
file_names_Jiali <- list.files(path1)
file_names_tianxia <- list.files(path2)
#ss_score_max <- vector()
#ss_score_avg <- vector()
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:100)){
  GO1 <- read_excel(paste0(path1,file_names_Jiali[k]), sheet = 2)
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
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_Adipo_random_gene_UP.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_Adipo_random_gene_UP.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_Adipo_random_gene_UP.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#fba6ab")+
  geom_vline(xintercept  = 0.811, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_Adipo_random_gene_UP.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_Adipo_random_gene_UP.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#fba6ab")+
  geom_vline(xintercept  = 0.51, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_Adipo_random_gene_UP.pdf', width=6, height=3)

############第4部，计算,下调##########
path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/Adipo_jiali_down/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/Adipo_tianxia_down/"
file_names_Jiali <- list.files(path1)
file_names_tianxia <- list.files(path2)
#ss_score_max <- vector()
#ss_score_avg <- vector()
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:100)){
  GO1 <- read_excel(paste0(path1,file_names_Jiali[k]), sheet = 2)
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
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_Adipo_random_gene_down.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_Adipo_random_gene_down.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_Adipo_random_gene_down.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#6bd2db")+
  geom_vline(xintercept  = 0.633, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_Adipo_random_gene_down.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_Adipo_random_gene_down.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#6bd2db")+
  geom_vline(xintercept  = 0.485, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_Adipo_random_gene_down.pdf', width=6, height=3)


############第5部，计算,上调##########
path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/TGA_jiali_up/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/TGA_tianxia_up/"
file_names_Jiali <- list.files(path1)
file_names_tianxia <- list.files(path2)
#ss_score_max <- vector()
#ss_score_avg <- vector()
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:100)){
  GO1 <- read_excel(paste0(path1,file_names_Jiali[k]), sheet = 2)
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
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_TGA_random_gene_UP.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_TGA_random_gene_UP.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_TGA_random_gene_UP.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#fba6ab")+
  geom_vline(xintercept  = 0.816, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_TGA_random_gene_UP.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_TGA_random_gene_UP.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#fba6ab")+
  geom_vline(xintercept  = 0.488, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_TGA_random_gene_UP.pdf', width=6, height=3)

############第6部，计算,下调##########
path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/TGA_jiali_down/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term/TGA_tianxia_down/"
file_names_Jiali <- list.files(path1)
file_names_tianxia <- list.files(path2)
#ss_score_max <- vector()
#ss_score_avg <- vector()
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:100)){
  GO1 <- read_excel(paste0(path1,file_names_Jiali[k]), sheet = 2)
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
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_TGA_random_gene_down.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_TGA_random_gene_down.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_TGA_random_gene_down.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#6bd2db")+
  geom_vline(xintercept  = 0.671, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_TGA_random_gene_down.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_TGA_random_gene_down.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#6bd2db")+
  geom_vline(xintercept  = 0.534, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_TGA_random_gene_down.pdf', width=6, height=3)


setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/36随机取基因/Random_GO_term(直接从随机的GOterm出发)")
##########100次循环，从GO的随机选择出发###########
#Adipo
#上调基因的GO数：嘉莉，21
#下调基因的GO数：嘉莉，50
#上调基因的GO数：田侠，197
#下调基因的GO数：田侠，104

#TGA
#上调基因的GO数：嘉莉，26
#下调基因的GO数：嘉莉，82
#上调基因的GO数：田侠，252
#下调基因的GO数：田侠，181

#读取总共的BP。
HsGO <- godata('org.Hs.eg.db', ont="BP")
AllGO_BP <- unique(HsGO@geneAnno$GO)

###########Adipo############
########上调的#############
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:500)){
  Jiali_up <- sample(AllGO_BP,size = 21)
  Tianxia_up <- sample(AllGO_BP,size = 197)
  go1 <- Jiali_up
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
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_Adipo_random_GOBP_UP.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_Adipo_random_GOBP_UP.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_Adipo_random_GOBP_UP.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#e9a0b9")+
  geom_vline(xintercept  = 0.811, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_Adipo_random_GOBP_UP.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_Adipo_random_GOBP_UP.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#e9a0b9")+
  geom_vline(xintercept  = 0.51, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_Adipo_random_GOBP_UP.pdf', width=6, height=3)

########下调的#############
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:500)){
  Jiali_down <- sample(AllGO_BP,size = 50)
  Tianxia_down <- sample(AllGO_BP,size = 104)
  go1 <- Jiali_down
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
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_Adipo_random_GOBP_down.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_Adipo_random_GOBP_down.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_Adipo_random_GOBP_down.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#72C8E5")+
  geom_vline(xintercept  = 0.633, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_Adipo_random_GOBP_down.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_Adipo_random_GOBP_down.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#72C8E5")+
  geom_vline(xintercept  = 0.485, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_Adipo_random_GOBP_down.pdf', width=6, height=3)


###########TGA############
########上调的#############
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:500)){
  Jiali_up <- sample(AllGO_BP,size = 26)
  Tianxia_up <- sample(AllGO_BP,size = 252)
  go1 <- Jiali_up
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
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_TGA_random_GOBP_UP.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_TGA_random_GOBP_UP.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_TGA_random_GOBP_UP.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#e9a0b9")+
  geom_vline(xintercept  = 0.816, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_TGA_random_GOBP_UP.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_TGA_random_GOBP_UP.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#e9a0b9")+
  geom_vline(xintercept  = 0.488, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_TGA_random_GOBP_UP.pdf', width=6, height=3)

########下调的#############
ss_score_rcmax <- vector()
ss_score_BMA <- vector()
for (k in c(1:500)){
  Jiali_down <- sample(AllGO_BP,size = 82)
  Tianxia_down <- sample(AllGO_BP,size = 181)
  go1 <- Jiali_down
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
write.csv(ss_score_rcmax[1:100],"ss_score_rcmax_TGA_random_GOBP_down.csv")
write.csv(ss_score_BMA[1:100],"ss_score_BMA_TGA_random_GOBP_down.csv")

#密度曲线
ss_score_rcmax <- read.csv("ss_score_rcmax_TGA_random_GOBP_down.csv")
p_rcmax <- ggplot(ss_score_rcmax, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#72C8E5")+
  geom_vline(xintercept  = 0.671, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_rcmax, file='ss_score_rcmax_TGA_random_GOBP_down.pdf', width=6, height=3)

ss_score_BMA <- read.csv("ss_score_BMA_TGA_random_GOBP_down.csv")
p_BMA <- ggplot(ss_score_BMA, aes(x= x))+ 
  geom_density(color ="#666666", fill ="#72C8E5")+
  geom_vline(xintercept  = 0.534, color ="#666666", linetype="dashed")+
  theme_classic()
ggsave(p_BMA, file='ss_score_BMA_TGA_random_GOBP_down.pdf', width=6, height=3)


