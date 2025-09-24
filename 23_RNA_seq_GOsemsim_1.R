#org.Ss.eg.db包用不了的时候，先运行这个。
#options(connectionObserver = NULL)
#library(topGO)
#library(pathview)

library(clusterProfiler)
library(org.Ss.eg.db)
library(stringr)
library("ggplot2")
library("GOSemSim")

setwd("E:/BIIF/liutianxia/lilan_RNAseq/综合比较丽兰和田侠的")
DEG_lilan <- read.csv("resdataDEG_lilan.csv")
DEG_tianxia <- read.csv("resdataDEG_tianxia.csv")

############上调，下调分开做##########
resdataDEG_go_shangtiao_lilan <- DEG_lilan[DEG_lilan$log2FoldChange > 0,]
difgene_shangtiao <- resdataDEG_go_shangtiao_lilan$Row.names
write.csv(difgene_shangtiao,"resdataDEG_lilan_up.csv")
difgeneID_shangtiao_lilan <- bitr(difgene_shangtiao, 
                  fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Ss.eg.db")

resdataDEG_go_xiatiao_lilan <- DEG_lilan[DEG_lilan$log2FoldChange < 0,]
difgene_xiatiao <- resdataDEG_go_xiatiao_lilan$Row.names
write.csv(difgene_xiatiao,"resdataDEG_lilan_down.csv")
difgeneID_xiatiao_lilan <- bitr(difgene_xiatiao, 
                            fromType="SYMBOL", 
                            toType="ENTREZID", 
                            OrgDb="org.Ss.eg.db")

#######GO用网页工具做#####
###############GO##############
#上调
#pvalueCutoff和qvalueCutoff两个参数反复调整，卡出来比较好的go term。
difgoALL_shangtiao_lilan <- enrichGO(gene = difgeneID_shangtiao_lilan$ENTREZID, 
                     OrgDb = org.Ss.eg.db, 
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pvalueCutoff = 0.05, #默认是0.05。其实这个卡的是p adj值。
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.2, #默认是0.2
                     readable = TRUE,
                     pool = FALSE)
write.csv(difgoALL_shangtiao_lilan@result,file= "difgoALL_shangtiao_lilan.csv",row.names = F) 
GOplot <- dotplot(difgoALL_shangtiao_lilan,font.size = 15)+
  scale_color_continuous(low="#f37736", high="#7bc043")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme_classic()+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
ggsave(GOplot,filename = "GO_shangtiao_lilan.pdf" ,width = 5,height = 5)

#下调
#pvalueCutoff和qvalueCutoff两个参数反复调整，卡出来比较好的go term。
difgoALL_xiatiao_lilan <- enrichGO(gene = difgeneID_xiatiao_lilan$ENTREZID, 
                     OrgDb = org.Ss.eg.db, 
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pvalueCutoff = 0.05, #默认是0.05。其实这个卡的是p adj值。
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.2, #默认是0.2
                     readable = TRUE,
                     pool = FALSE)
write.csv(difgoALL_xiatiao_lilan@result,file= "difgoALL_xiatiao_lilan.csv",row.names = F) 
GOplot <- dotplot(difgoALL_xiatiao_lilan,font.size = 15)+
  scale_color_continuous(low="#f37736", high="#7bc043")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme_classic()+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
ggsave(GOplot,filename = "GO_xiatiao_lilan.pdf" ,width = 5,height = 5)



############上调，下调分开做##########
resdataDEG_go_shangtiao_tianxia <- DEG_tianxia[DEG_tianxia$log2FoldChange > 0,]
difgene_shangtiao <- resdataDEG_go_shangtiao_tianxia$Row.names
write.csv(difgene_shangtiao,"resdataDEG_tianxia_up.csv")
difgeneID_shangtiao_tianxia <- bitr(difgene_shangtiao, 
                            fromType="SYMBOL", 
                            toType="ENTREZID", 
                            OrgDb="org.Ss.eg.db")

resdataDEG_go_xiatiao_tianxia <- DEG_tianxia[DEG_tianxia$log2FoldChange < 0,]
difgene_xiatiao <- resdataDEG_go_xiatiao_tianxia$Row.names
write.csv(difgene_xiatiao,"resdataDEG_tianxia_down.csv")
difgeneID_xiatiao_tianxia <- bitr(difgene_xiatiao, 
                          fromType="SYMBOL", 
                          toType="ENTREZID", 
                          OrgDb="org.Ss.eg.db")


#######GO用网页工具做#####
###############GO##############
#上调
#pvalueCutoff和qvalueCutoff两个参数反复调整，卡出来比较好的go term。
difgoALL_shangtiao_tianxia <- enrichGO(gene = difgeneID_shangtiao_tianxia$ENTREZID, 
                                     OrgDb = org.Ss.eg.db, 
                                     keyType = "ENTREZID",
                                     ont = "ALL",
                                     pvalueCutoff = 0.05, #默认是0.05。其实这个卡的是p adj值。
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2, #默认是0.2
                                     readable = TRUE,
                                     pool = FALSE)
write.csv(difgoALL_shangtiao_tianxia@result,file= "difgoALL_shangtiao_tianxia.csv",row.names = F) 
GOplot <- dotplot(difgoALL_shangtiao_tianxia,font.size = 15)+
  scale_color_continuous(low="#f37736", high="#7bc043")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme_classic()+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
ggsave(GOplot,filename = "GO_shangtiao_tianxia.pdf" ,width = 5,height = 5)

#下调
#pvalueCutoff和qvalueCutoff两个参数反复调整，卡出来比较好的go term。
difgoALL_xiatiao_tianxia <- enrichGO(gene = difgeneID_xiatiao_tianxia$ENTREZID, 
                                   OrgDb = org.Ss.eg.db, 
                                   keyType = "ENTREZID",
                                   ont = "ALL",
                                   pvalueCutoff = 0.05, #默认是0.05。其实这个卡的是p adj值。
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.2, #默认是0.2
                                   readable = TRUE,
                                   pool = FALSE)
write.csv(difgoALL_xiatiao_tianxia@result,file= "difgoALL_xiatiao_tianxia.csv",row.names = F) 
GOplot <- dotplot(difgoALL_xiatiao_tianxia,font.size = 15)+
  scale_color_continuous(low="#f37736", high="#7bc043")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme_classic()+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
ggsave(GOplot,filename = "GO_xiatiao_tianxia.pdf" ,width = 5,height = 5)




##########GOsemsim计算语义分析############
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

HsGO <- godata('org.Hs.eg.db', ont="BP")
#SsGO <- godata('org.Ss.eg.db', ont="BP")
#########上调############
lll <- read.csv("GO_up_lilan.csv")
#截取|竖线之后的字符
for (i in 1:length(lll$GroupID)){
  lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
}
lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")

ttt <- read.csv("GO_up_tianxia.csv")
#截取|竖线之后的字符
for (i in 1:length(ttt$GroupID)){
  ttt$GroupID[i] <- substring(ttt$GroupID[i],indexOf(ttt$GroupID[i],"_")+1,50)
}
ttt <- subset(ttt, GroupID == "Member" & Category == "GO Biological Processes")
go1 <- lll$Term
go2 <- ttt$Term

#备注，这3个都可以
#One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="avg")
#0.155
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
#0.679
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
#0.718
intersect(go1,go2)


#########下调############
lll <- read.csv("GO_down_lilan.csv")
#截取|竖线之后的字符
for (i in 1:length(lll$GroupID)){
  lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
}
lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")

ttt <- read.csv("GO_down_tianxia.csv")
#截取|竖线之后的字符
for (i in 1:length(ttt$GroupID)){
  ttt$GroupID[i] <- substring(ttt$GroupID[i],indexOf(ttt$GroupID[i],"_")+1,50)
}
ttt <- subset(ttt, GroupID == "Member" & Category == "GO Biological Processes")
go1 <- lll$Term
go2 <- ttt$Term

#One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="avg")
#0.1,不行
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
#0.482，不行
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
#0.703
intersect(go1,go2)


###########维恩######
library(ggvenn)
plot_venn <- ggvenn(data = list(lilan_up_Gene = resdataDEG_go_shangtiao_lilan$Row.names,
                                tianxia_up_Gene = resdataDEG_go_shangtiao_tianxia$Row.names),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("RNA_seq_维恩图_up_Gene_共有的.pdf",sep = ""), width=3, height=3)
gene <- intersect(resdataDEG_go_shangtiao_lilan$Row.names, resdataDEG_go_shangtiao_tianxia$Row.names)
write.csv(gene,"RNA_seq_维恩图_up_Gene_共有的.csv")

plot_venn <- ggvenn(data = list(lilan_up_Gene = resdataDEG_go_xiatiao_lilan$Row.names,
                                tianxia_up_Gene = resdataDEG_go_xiatiao_tianxia$Row.names),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("RNA_seq_维恩图_down_Gene_共有的.pdf",sep = ""), width=3, height=3)
gene <- intersect(resdataDEG_go_xiatiao_lilan$Row.names, resdataDEG_go_xiatiao_tianxia$Row.names)
write.csv(gene,"RNA_seq_维恩图_down_Gene_共有的.csv")

#########GO上调############
lll <- read.csv("GO_up_lilan.csv")
#截取|竖线之后的字符
for (i in 1:length(lll$GroupID)){
  lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
}
lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")

ttt <- read.csv("GO_up_tianxia.csv")
#截取|竖线之后的字符
for (i in 1:length(ttt$GroupID)){
  ttt$GroupID[i] <- substring(ttt$GroupID[i],indexOf(ttt$GroupID[i],"_")+1,50)
}
ttt <- subset(ttt, GroupID == "Member" & Category == "GO Biological Processes")
go1 <- lll$Term
go2 <- ttt$Term
plot_venn <- ggvenn(data = list(lilan_up_GO = go1,
                                tianxia_up_GO = go2),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("RNA_seq_维恩图_up_GO_共有的.pdf",sep = ""), width=3, height=3)
RP1 <- subset(lll,lll$Term %in% intersect(go1,go2))
write.csv(RP1,"RNA_seq_维恩图_up_GO_共有的.csv")

#########GO下调############
lll <- read.csv("GO_down_lilan.csv")
#截取|竖线之后的字符
for (i in 1:length(lll$GroupID)){
  lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
}
lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")

ttt <- read.csv("GO_down_tianxia.csv")
#截取|竖线之后的字符
for (i in 1:length(ttt$GroupID)){
  ttt$GroupID[i] <- substring(ttt$GroupID[i],indexOf(ttt$GroupID[i],"_")+1,50)
}
ttt <- subset(ttt, GroupID == "Member" & Category == "GO Biological Processes")
go1 <- lll$Term
go2 <- ttt$Term
plot_venn <- ggvenn(data = list(lilan_down_GO = go1,
                                tianxia_down_GO = go2),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("RNA_seq_维恩图_down_GO_共有的.pdf",sep = ""), width=3, height=3)
RP1 <- subset(lll,lll$Term %in% intersect(go1,go2))
write.csv(RP1,"RNA_seq_维恩图_down_GO_共有的.csv")



