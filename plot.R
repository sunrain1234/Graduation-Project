profile_text<-read.table('D:/clinical/TCGA-HNSC/TCGA-HNSC-mRNA/Cisplatin-AUC/mRNA-AUC-Cisplatin_new.csv',header = T, sep = ',')
library(ggplot2)
library(reshape2)
data_m<-melt(profile_text)
head(data_m)
p <- ggplot(data_m, aes(x=variable, y=value),color=variable) + 
  geom_boxplot(aes(fill=factor(variable))) + 
  geom_point(position=position_jitter(width=.2, height=0))+
  labs(title = "mRNA-Cisplatin-HNSC-AUC", x = "Classifier", y = "AUC")+
  theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = seq(0.0,1.0,0.1))
p

AUC<-read.csv('D:/clinical/TCGA-BLCA/TCGA-BLCA-mRNA/Cisplatin-AUC/mRNA_AUC_Cisplatin.csv',header = T)
typeof(AUC$DummmyClassifer)
typeof(AUC$RBF.SVM)
as.numeric(AUC$DummmyClassifer)
typeof(AUC$DummmyClassifer)
as.numeric(AUC$RBF.SVM)
is.numeric(AUC$DummmyClassifer)
is.numeric(AUC$RBF.SVM)
wilcox.test(AUC$DummmyClassifier,AUC$RBF.SVM,alternative = "less",exact=F)
typeof(AUC)
new <- c(76,152,243,82,240,220,205,38,243,44,190,100)
wilcox.test(old,new,paired = TRUE)
new
typeof(new)

miRNA<-read.csv('D:/clinical/TCGA-BLCA/TCGA-BLCA-miRNA/Merge_matrix.csv',row.names = 1,header = T)
sum(miRNA[1,]==0)

GOTERM_BP_FAT<-read_xlsx('D:/clinical/TCGA-BLCA/TCGA-BLCA-miRNA/Carboplatin-AUC/GOTERM_MF_FAT.xlsx')
GOTERM_BP_FAT_1 <- GOTERM_BP_FAT[which(GOTERM_BP_FAT$PValue<0.05),]

GOTERM_BP_FAT_1$PValue<--log(GOTERM_BP_FAT_1$PValue,2)
#GOTERM_BP_FAT_2 <- GOTERM_BP_FAT_1[order(GOTERM_BP_FAT_1[,12]),]
GOTERM_BP_FAT_2 <- arrange(GOTERM_BP_FAT_1,-PValue)
GOTERM_BP_FAT_2 <- GOTERM_BP_FAT_2[1:10,]
#barplot(rev(GOTERM_BP_FAT_2$Benjamini),horiz=T,xlim=c(-4,1),axes=F,col=rep(brewer.pal(9,'YlOrRd'),each=15))  
#text(seq(from=0.7,length.out=135,by=1.2),x=-2,label=rev(GOTERM_BP_FAT_2$Term))  
#axis(3,c(0,2.5,5,7.5,10),c('0','2.5','5.0','7.5','10.0'))  
#barplot(GOTERM_BP_FAT_2$Benjamini,beside = TRUE,col = 'azure4',axes = FALSE,horiz = TRUE)
#axis(1)
#axis(2,labels = GOTERM_BP_FAT_2$Term,tick = FALSE)

#ggplot(GOTERM_BP_FAT_2,aes(Term,fill=Benjamini))+
#  geom_bar(stat="identity",position="dodge")+
#  theme_wsj()+
#  scale_fill_wsj("rgby", "")+
#  theme(axis.ticks.length=unit(0.5,'cm'))+
#  guides(fill=guide_legend(title=NULL))+
#  ggtitle("The Financial Performance of Five Giant")+
#  theme(axis.title = element_blank(),legend.position='none')+
#  facet_grid(.~Benjamini)+
#  coord_flip()
GOTERM_BP_FAT_2$Term <- substr(GOTERM_BP_FAT_2$Term,12,100)
GOTERM_BP_FAT_2$Term = factor(GOTERM_BP_FAT_2$Term,levels=GOTERM_BP_FAT_2$Term)
ggplot(data = GOTERM_BP_FAT_2,mapping = aes(x=Term,y=PValue))+
   geom_bar(stat = 'identity',width = 0.2)+
   coord_flip()+
  labs(title = 'GOTERM_MF_FAT', y = "-log2(p-value)")

Cisplatin_BLCA_mRNA_sorted_again<-read.csv('D:/clinical/TCGA-BLCA/TCGA-BLCA-mRNA/Cisplatin-BLCA-mRNA-sorted-again.csv',header=T)
pheatmap(Cisplatin_BLCA_mRNA_sorted_again)

#对差异基因做显著性分析
mRNA_Cisplatin_BLCA <- read.csv('D:/clinical/TCGA-BLCA/TCGA-BLCA-mRNA/Cisplatin-BLCA-mRNA-sorted-again.csv',header = T)
sensitive <- mRNA_Cisplatin_BLCA[,2:38]
insensitive <- mRNA_Cisplatin_BLCA[,39:58]
p_value = NULL
for (i in 1:length(sensitive[,1])){
  p_value <- union(p_value,t.test(sensitive[i,],insensitive[i,])$p.value)
  
}

#画基因聚类热图
Cisplatin_BLCA_mRNA_sorted_again <- read.csv('D:/clinical/TCGA-BRCA/TCGA-BRCA-mRNA/Paclitaxel-BRCA-mRNA-sorted-again.csv',header = T,row.names = 1)
response <- rep('sensitive',53)
non_response <-rep('insensitive',7)
a<-data.frame(response=factor(c(response,non_response)))

rownames(a)<-colnames(Cisplatin_BLCA_mRNA_sorted_again)
#b<-as.data.frame(Cisplatin_BLCA_mRNA_sorted_again)
#rownames(b)<-Cisplatin_BLCA_mRNA_sorted_again$gene
pheatmap(as.data.frame(Cisplatin_BLCA_mRNA_sorted_again),cluster_rows = FALSE,cluster_cols = TRUE,
          annotation_col = a,annotation_names_col = F,scale = 'row',clustering_method = "complete")


result<-read.csv('D:/Bio/result.csv',header = F,row.names = 1)
result_1<-result[1:12634,1:179]
result_2<-result[1:12634,180:286]
rows<-NULL
for(i in 1:length(result_1[,1])){
  pvalue<-t.test(result_1[i,],result_2[i,])$p.value
  if(pvalue<0.01){
    rows<-union(rows,i)
  }
}
new_result<-result[rows,]
write.csv(new_result,'D:/Bio/new_result_2.csv',row.names = F)
new_data<-read.csv('D:/Bio/new_result_2.csv',header = T)

new_data<-t(new_data)
new_data<-data.frame(new_data)
sen<-rep(1,179)
insen<-rep(2,286-179)
total<-data.frame(response=factor(c(sen,insen)))
rownames(total)<-rownames(new_data)

tsne<-Rtsne(new_data,dim=2,perplexity = 50)
final_data<-cbind(tsne$Y,total)
final_data$response<-factor(final_data$response)
ggplot(data=final_data,mapping = aes(x=`1`,y=`2`,colour=response))+geom_point(size=3)

pheatmap(final_data[,1:2],cluster_rows = T,cluster_cols = F,annotation_row = total,annotation_names_row = F,clustering_method = "complete")
pheatmap(new_data,cluster_rows = T,cluster_cols = F,annotation_row = total,annotation_names_row = F,clustering_method = "complete")
