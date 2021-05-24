library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(regexPipes)
library(dplyr)
library(stringr)
library(DESeq2)
library(ggplot2)

#####下載TCGA_BRCA資料
query <- GDCquery(project = "TCGA-BRCA", 
                  legacy = FALSE, 
                  experimental.strategy = "RNA-Seq", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
data <- GDCprepare(query, 
                   save = TRUE, 
                   directory =  "D:/user/Documents/TNBC/GDCdata",   
                   save.filename = "BRCA.RData")
count_matrix <- assay(data)
dataExp <- as.data.frame(count_matrix)

#####下載臨床資料
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
clinical_2<-colData(data)

#####找出Ensembl ID對應的Gene ID
testEID <- rownames(dataExp)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = testEID, mart= mart)
write.csv(gene_IDs,file ="gene_IDs")

#####找出三陰性乳癌患者的樣品
test2 <-as.data.frame(t(dataExp))
cancer <-subset(test2,substring(rownames(test2),14,15)=="01")
cancer <-mutate(cancer,caseID=rownames(cancer),case=str_sub(rownames(cancer),1,12))
tnbcid<-as.data.frame(TNBC_case$bcr_patient_barcode) 
test3=dplyr::inner_join(cancer,tnbcid,by="case")

#####將ensembl ID轉成gene ID
dataExp <-mutate(dataExp,ensembl_gene_id=rownames(dataExp))
dataExp$ensembl_gene_id<-as.character(dataExp$ensembl_gene_id)
dataExp <- dataExp %>%
  inner_join(gene_IDs,by="ensembl_gene_id") %>%
  select(-ensembl_gene_id) %>%
  select(symbol,everything())%>%
  mutate(rowmean=rowMeans(.[grep("TCGA",names(.))]))%>%
  arrange(desc(rowmean))%>%
  distinct(symbol,.keep_all = T)%>%
  select(-rowmean)%>%
  tibble::column_to_rownames(colnames(.)[1])

#####找出TCGA-BRCA正常樣本
metadata <- data.frame(TCGA_id=colnames(dataExp))
table(substring(metadata$TCGA_id,14,15))
tdataEXP <-as.data.frame(t(dataExp))
normal <-subset(tdataEXP,substring(rownames(tdataEXP),14,15)=="11")
normal <-as.data.frame(t(normal))
write.csv(normal,file = "normal_Exp.csv")

#####將正常與三陰性行樣本merge另存成新的表達矩陣
cancer <-mutate(cancer,gene_id=rownames(cancer))
normal <- mutate(normal,gene_id=rownames(normal))
BRCAexp <- inner_join(cancer,normal,by="gene_id")
rownames(BRCAexp) <-BRCAexp$gene_id
BRCAexp <-dplyr::select(BRCAexp,-gene_id)
write.csv(BRCAexp,file ="BRCAexp.csv")


######DEseq2
dds <- DESeqDataSetFromMatrix(countData = mycounts,
                              colData = coldata,
                              design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <-DESeq(dds)
save(dds,file="BRCA_DEseq2(1)_dds.Rdata")
res<-results(dds,contrast = c("condition","cancer","normal"))
res<-as.data.frame(res)
res <-cbind(rownames(res),res)
colnames(res)[1] <-"gene_id"
res$change <-as.factor(ifelse(res$padj<0.05 & abs(res$log2FoldChange)>1,ifelse(res$log2FoldChange>1,"UP","DOWN"),"NOT"))
resSig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1),]
write.csv(res,file = "BRCA_DEseq2(1).csv")
write.csv(resSig,file = "BRCA_DEseq2(1)_DEG.csv")

BRCA_DEG_gene_id <-rownames(resSig)
BRCA_DEG_gene_id <-as.data.frame(BRCA_DEG_gene_id)
write.csv(BRCA_DEG_gene_id,file = "BRCA_DEG(1)_gene_id.csv")


#####提取有DEG意義的基因表達矩陣
BRCAexp_vst <-cbind(rownames(BRCAexp_vst),BRCAexp_vst)
colnames(BRCAexp_vst)[1]<-"gene_id"
BRCAexp_F<-dplyr::inner_join(BRCAexp_vst,resSig,by="gene_id")

#####畫DEG火山圖
ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(alpha=0.8, size = 1) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) + 
  geom_hline(yintercept=2,linetype=4) + 
  geom_vline(xintercept=c(-1,1) ,linetype=4 ) +
  scale_color_manual(name = "", values = c("red", "green", "black"), limits = c("UP", "DOWN", "NOT"))