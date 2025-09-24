library("GOSemSim")
library("readxl")
library("ggplot2")
setwd("E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/35语义分析GOsemsim")

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

path1 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/35语义分析GOsemsim/嘉莉metascape/"
path2 <- "E:/BIIF/liutianxia/liujiali_藏猪单细胞/3_结果/35语义分析GOsemsim/田侠metascape/"

jiali_Adipo_up_GO <- read_excel(paste0(path1,"SZ_DEGs_总体_UP_GO.xlsx"), sheet = 2)
jiali_Adipo_down_GO <- read_excel(paste0(path1,"SZ_DEGs_总体_DOWN_GO.xlsx"), sheet = 2)
jiali_TGA_up_GO <- read_excel(paste0(path1,"SZ_DEGs_5 TGA_UP_GO.xlsx"), sheet = 2)
jiali_TGA_down_GO <- read_excel(paste0(path1,"SZ_DEGs_5 TGA_DOWN_GO.xlsx"), sheet = 2)

Tianxia_Adipo_up_GO <- read_excel(paste0(path2,"pWAT_DEGs_1 Adipocyte_UP_GO.xlsx"), sheet = 2)
Tianxia_Adipo_down_GO <- read_excel(paste0(path2,"pWAT_DEGs_1 Adipocyte_DOWN_GO.xlsx"), sheet = 2)
Tianxia_TGA_up_GO <- read_excel(paste0(path2,"SZ_DEGs_4 C4_UP_GO.xlsx"), sheet = 2)
Tianxia_TGA_down_GO <- read_excel(paste0(path2,"SZ_DEGs_4 C4_DOWN_GO.xlsx"), sheet = 2)



HsGO <- godata('org.Hs.eg.db', ont="BP")
#SsGO <- godata('org.Ss.eg.db', ont="BP")

#########Adipo 上调############
lll <- jiali_Adipo_up_GO
#截取|竖线之后的字符
for (i in 1:length(lll$GroupID)){
  lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
}
lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")

ttt <- Tianxia_Adipo_up_GO
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
#0.179
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
#0.51
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
#0.811
intersect(go1,go2)

#########Adipo 下调############
lll <- jiali_Adipo_down_GO
#截取|竖线之后的字符
for (i in 1:length(lll$GroupID)){
  lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
}
lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")

ttt <- Tianxia_Adipo_down_GO
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
#0.108
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
#0.485
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
#0.633
intersect(go1,go2)

#########TGA 上调############
lll <- jiali_TGA_up_GO
#截取|竖线之后的字符
for (i in 1:length(lll$GroupID)){
  lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
}
lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")

ttt <- Tianxia_TGA_up_GO
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
#0.148
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
#0.488
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
#0.816
intersect(go1,go2)


#########TGA 下调############
lll <- jiali_TGA_down_GO
#截取|竖线之后的字符
for (i in 1:length(lll$GroupID)){
  lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
}
lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")

ttt <- Tianxia_TGA_down_GO
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
#0.11
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="BMA")
#0.534
mgoSim(go1, go2, semData=HsGO, measure="Wang", combine="rcmax")
#0.671
intersect(go1,go2)


########维恩图###############

jiali_Adipo_up_GO <- read_excel(paste0(path1,"SZ_DEGs_总体_UP_GO.xlsx"), sheet = 2)
jiali_Adipo_down_GO <- read_excel(paste0(path1,"SZ_DEGs_总体_DOWN_GO.xlsx"), sheet = 2)
jiali_TGA_up_GO <- read_excel(paste0(path1,"SZ_DEGs_5 TGA_UP_GO.xlsx"), sheet = 2)
jiali_TGA_down_GO <- read_excel(paste0(path1,"SZ_DEGs_5 TGA_DOWN_GO.xlsx"), sheet = 2)

Tianxia_Adipo_up_GO <- read_excel(paste0(path2,"pWAT_DEGs_1 Adipocyte_UP_GO.xlsx"), sheet = 2)
Tianxia_Adipo_down_GO <- read_excel(paste0(path2,"pWAT_DEGs_1 Adipocyte_DOWN_GO.xlsx"), sheet = 2)
Tianxia_TGA_up_GO <- read_excel(paste0(path2,"SZ_DEGs_4 C4_UP_GO.xlsx"), sheet = 2)
Tianxia_TGA_down_GO <- read_excel(paste0(path2,"SZ_DEGs_4 C4_DOWN_GO.xlsx"), sheet = 2)

#########这个代码非常高级###########
#循环，把每个Term 都只留下 GO BP的 term
for (k in c("jiali_Adipo_up_GO","jiali_Adipo_down_GO","jiali_TGA_up_GO","jiali_TGA_down_GO",
           "Tianxia_Adipo_up_GO","Tianxia_Adipo_down_GO","Tianxia_TGA_up_GO","Tianxia_TGA_down_GO")){
  lll <- get(k) #获得这个字符串所包含的变量
  #截取|竖线之后的字符
  for (i in 1:length(lll$GroupID)){
    lll$GroupID[i] <- substring(lll$GroupID[i],indexOf(lll$GroupID[i],"_")+1,50)
  }
  lll <- subset(lll, GroupID == "Member" & Category == "GO Biological Processes")
  assign(k,lll) # 把变量重新赋值给一个变量名
}

#1
plot_venn <- ggvenn(data = list(jiali_Adipo_up_GO = jiali_Adipo_up_GO$Term,
                                Tianxia_Adipo_up_GO = Tianxia_Adipo_up_GO$Term),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("jiali_Adipo_up_GO和Tianxia_Adipo_up_GO.pdf",sep = ""), width=3, height=3)
RP1 <- intersect(jiali_Adipo_up_GO$Description,Tianxia_Adipo_up_GO$Description)
write.csv(RP1,"jiali_Adipo_up_GO和Tianxia_Adipo_up_GO.csv")

#2
plot_venn <- ggvenn(data = list(jiali_Adipo_down_GO = jiali_Adipo_down_GO$Term,
                                Tianxia_Adipo_down_GO = Tianxia_Adipo_down_GO$Term),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("jiali_Adipo_down_GO和Tianxia_Adipo_down_GO.pdf",sep = ""), width=3, height=3)
RP1 <- intersect(jiali_Adipo_down_GO$Description,Tianxia_Adipo_down_GO$Description)
write.csv(RP1,"jiali_Adipo_down_GO和Tianxia_Adipo_down_GO.csv")

#3
plot_venn <- ggvenn(data = list(jiali_TGA_up_GO = jiali_TGA_up_GO$Term,
                                Tianxia_TGA_up_GO = Tianxia_TGA_up_GO$Term),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("jiali_TGA_up_GO和Tianxia_TGA_up_GO.pdf",sep = ""), width=3, height=3)
RP1 <- intersect(jiali_TGA_up_GO$Description,Tianxia_TGA_up_GO$Description)
write.csv(RP1,"jiali_TGA_up_GO和Tianxia_TGA_up_GO.csv")

#4
plot_venn <- ggvenn(data = list(jiali_TGA_down_GO = jiali_TGA_down_GO$Term,
                                Tianxia_TGA_down_GO = Tianxia_TGA_down_GO$Term),
                    show_percentage = F,
                    stroke_color = "black",
                    fill_alpha = 0.7,
                    stroke_size = 0.8,
                    set_name_size = 3,
                    text_size = 4,
                    fill_color = c("#46667a","#6f58a9"),
                    set_name_color = c("#46667a","#6f58a9"))
ggsave(plot_venn, file = paste("jiali_TGA_down_GO和Tianxia_TGA_down_GO.pdf",sep = ""), width=3, height=3)
RP1 <- intersect(jiali_TGA_down_GO$Description,Tianxia_TGA_down_GO$Description)
write.csv(RP1,"jiali_TGA_down_GO和Tianxia_TGA_down_GO.csv")




