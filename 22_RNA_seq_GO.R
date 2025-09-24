#org.Ss.eg.db包用不了的时候，先运行这个。
#options(connectionObserver = NULL)
#library(topGO)
#library(pathview)

library(clusterProfiler)
library(org.Ss.eg.db)
library(stringr)
library(KEGG.db)

#DEGs
resdataDEG <- subset(resdata,pvalue < 0.05)
#write.csv(resdataDEG,file = "resdataDEG.csv")

############上调，下调分开做##########
resdataDEG_go_shangtiao <- resdataDEG [resdataDEG$log2FoldChange > 0,]
difgene_shangtiao <- resdataDEG_go_shangtiao$Row.names
difgeneID_shangtiao <- bitr(difgene_shangtiao, 
                  fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Ss.eg.db")

resdataDEG_go_xiatiao <- resdataDEG [resdataDEG$log2FoldChange < 0,]
difgene_xiatiao <- resdataDEG_go_xiatiao$Row.names
difgeneID_xiatiao <- bitr(difgene_xiatiao, 
                            fromType="SYMBOL", 
                            toType="ENTREZID", 
                            OrgDb="org.Ss.eg.db")

#上调
#pvalueCutoff和qvalueCutoff两个参数反复调整，卡出来比较好的go term。
difgoALL <- enrichGO(gene = difgeneID_shangtiao$ENTREZID, 
                     OrgDb = org.Ss.eg.db, 
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pvalueCutoff = 0.05, #默认是0.05。其实这个卡的是p adj值。
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.2, #默认是0.2
                     readable = TRUE,
                     pool = FALSE)
write.csv(difgoALL@result,file= "difgoALL_shangtiao.csv",row.names = F) 
GOplot <- dotplot(difgoALL,font.size = 15)+
  scale_color_continuous(low="#f37736", high="#7bc043")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme_classic()+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
ggsave(GOplot,filename = "GO_shangtiao.pdf" ,width = 5,height = 5)

#下调
#pvalueCutoff和qvalueCutoff两个参数反复调整，卡出来比较好的go term。
difgoALL <- enrichGO(gene = difgeneID_xiatiao$ENTREZID, 
                     OrgDb = org.Ss.eg.db, 
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pvalueCutoff = 0.05, #默认是0.05。其实这个卡的是p adj值。
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.2, #默认是0.2
                     readable = TRUE,
                     pool = FALSE)
write.csv(difgoALL@result,file= "difgoALL_xiatiao.csv",row.names = F) 
GOplot <- dotplot(difgoALL,font.size = 15)+
  scale_color_continuous(low="#f37736", high="#7bc043")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme_classic()+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
ggsave(GOplot,filename = "GO_xiatiao.pdf" ,width = 5,height = 5)




#KEGG是日本主导的一个项目对gene和genome进行了非常详细的注释
#KEGG 需要包KEGG.db
difkk <- enrichKEGG(gene = difgeneID_shangtiao$ENTREZID,
                    organism ="ssc", #猪
                    keyType = "kegg",
                    pvalueCutoff = 0.01,#默认是0.05
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 500,
                    qvalueCutoff = 0.05,#默认是0.2
                    #readable = TRUE,
                    #use_internal_data =TRUE
)
write.csv(difkk@result,"KEGG_上调.csv")
GOplot <- dotplot(difkk,font.size = 15)+
  scale_color_continuous(low="#f37736", high="#7bc043")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme_classic()+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
ggsave(GOplot,filename = "KEGG_上调.pdf" ,width = 5,height = 5)


difkk <- enrichKEGG(gene = difgeneID_xiatiao$ENTREZID,
                    organism ="ssc", #猪
                    keyType = "kegg",
                    pvalueCutoff = 0.01,#默认是0.05
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 500,
                    qvalueCutoff = 0.05,#默认是0.2
                    #readable = TRUE,
                    #use_internal_data =TRUE
)
write.csv(difkk@result,"KEGG_下调.csv")
GOplot <- dotplot(difkk,font.size = 15)+
  scale_color_continuous(low="#f37736", high="#7bc043")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme_classic()+
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
ggsave(GOplot,filename = "KEGG_下调.pdf" ,width = 5,height = 5)




