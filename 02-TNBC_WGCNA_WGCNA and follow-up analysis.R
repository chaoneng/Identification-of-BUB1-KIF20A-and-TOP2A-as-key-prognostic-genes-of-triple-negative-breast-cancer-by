library(dplyr)
library(metaMA)
library(WGCNA)
library(flashClust)
library(stringr)
#####
count_matrix <- assay(data)
dataExp <- as.data.frame(count_matrix)

#####將TCGA-BRCA的FPKM數據的ensembl ID轉成gene ID
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

#####將TCGA-BRCA的FPKM數據進行整理
TNBC_FPKM<-as.data.frame(t(dataExp))
dataExp <-TNBC_FPKM
dataExp <-mutate(dataExp,gene_id=rownames(dataExp))
dataExp<-dplyr::inner_join(dataExp,BRCA_GENE_id,by="gene_id")
rownames(dataExp)<-dataExp$gene_id
dataExp <-dplyr::select(dataExp,-gene_id)
dataExp <-as.data.frame(t(dataExp))
dataExp <-mutate(dataExp,tcga_id=rownames(dataExp))
dataExp <-dplyr::inner_join(dataExp,BRCA_TCGA_id,by="tcga_id")
rownames(dataExp) <-dataExp$tcga_id
dataExp <-dplyr::select(dataExp,-tcga_id)
write.csv(TNBC_FPKM,file = "BRCA.dataEXP_FPKM.CSV")
dataExp<-as.data.frame(t(BRCA_dataEXP_FPKM))

######WGCNA樣本聚類離群值檢查
gsg = goodSamplesGenes(dataExp, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(dataExp), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
clust = cutreeStatic(sampleTree, cutHeight = 35000, minSize = 10)
keepSamples = (clust==1)
dataExp = dataExp[keepSamples, ]
sampleTree = hclust(dist(dataExp), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="")
nGenes = ncol(dataExp)
nSamples = nrow(dataExp)

#####WGCNA軟閥值篩選
powers= c(seq(1,10,by=1), seq(from =12, to=20, by=2))
sft = pickSoftThreshold(dataExp, powerVector=powers, verbose =5,networkType="unsigned")
par(mfrow = c(1,2))
cex1=0.8
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1,
     col="red")
abline(h=0.80, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean
Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
power = sft$powerEstimate

#####WGCNA畫出K值方圖
k <- softConnectivity(dataExp,power=5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main=paste("Check Scale free topology\n","power=",5))

#####sample traitdata
sample <- ifelse(substring(rownames(dataExp),14,15)!="11","cancer","normal")
traitdata<-data.frame(row.names = rownames(dataExp),condition=sample)
datT1 <-with(traitdata,outer(condition,levels(condition), `==`)*1)
rownames(datT1) <-rownames(traitdata)
colnames(datT1) <-paste(levels(traitdata$condition))
datatrait<-as.data.frame(datT1)

######建構共表達矩陣WGCNA
adjacency=adjacency(dataExp, power=5, type="unsigned")
TOM= TOMsimilarity(adjacency, TOMType="unsigned")
dissTOM=1-TOM
geneTree= flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels=
       FALSE, hang=0.04)
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2,
                           pamRespectsDendro= FALSE, minClusterSize= 100)
dynamicColors= labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE,
                    hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")
MEList= moduleEigengenes(dataExp, colors= dynamicColors)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
MEDissThres=0.25
abline(h=MEDissThres, col="red")
merge= mergeCloseModules(dataExp, dynamicColors, cutHeight=0.25, verbose =3)
mergedColors= merge$colors
mergedMEs= merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut","Merged dynamic"), dendroLabels= FALSE,
                    hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file= "SamplesAndColors.RData")
MEs0 = moduleEigengenes(dataExp, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleGeneCor=cor(MEs,dataExp)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#####特徵性狀關聯表&圖
moduleTraitCor = cor(MEs, datatrait, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix =paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitPvalue)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(datatrait), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix =
                 textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))

#####導出關鍵模塊基因
module="blue"
gene=colnames(dataExp)
inmodule=(moduleColors==module)
modgene=as.data.frame(gene[inmodule])
write.csv(modgene,file = "bluemod_gene.csv")

#####篩選hub gene
ADJ1=abs(cor(dataExp,use="p"))^power 
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
datKME=signedKME(dataExp, MEs, outputColumnName="MM.")
GS1=as.data.frame(cor(dataExp,datatrait[,1],use = "p"))
FilterGenes_spe = ((abs(GS1) > 0.2) & (abs(datKME$MM.blue)>0.8))
trait_hubGenes_spe <- as.data.frame(colnames(dataExp)[FilterGenes_spe])
write.csv(trait_hubGenes_spe,file = "allmod_hub_gene.csv")

#####篩選在blue模塊內的 hub gene
modgene1<-modgene
names(modgene1)="gene"
blue_hubgene <-trait_hubGenes_spe
names(blue_hubgene)="gene"
blue_hubgene <-dplyr::inner_join(blue_hubgene,modgene1,by="gene")
write.csv(blue_hubgene,file = "blue_hub_gene.csv")

#####GO & KEGG分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
testmodgene = bitr(blue_hub_gene$gene,fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene=testmodgene$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
go_ALL <-as.data.frame(ego)
dotplot(ego, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

kk <-enrichKEGG(gene = testmodgene$ENTREZID,
                organism = 'hsa',
                pvalueCutoff = 1
                )
kegg <-as.data.frame(kk)
barplot(kk,showCategory=9,title="Enrichment KEGG")
save(kk, ego, file= "GO&KEGG.RData")

#####檢查關鍵模塊內基因的表達差異狀況
names(modgene)<-"gene_id"
pp<-modgene
pp<-dplyr::inner_join(pp,`BRCA_DEseq2(1)_DEG`,by="gene_id")
up_pp=(pp$change=="UP")
UP_modgene=as.data.frame(pp$gene_id[up_pp])
d_pp=(pp$change=="DOWN")
down_modgene=as.data.frame(pp$gene_id[d_pp])
write.csv(UP_modgene,file = "up_modgene.csv")
write.csv(down_modgene,file = "down_modgene.csv")


#####生存分析
library(RTCGA)
library(RTCGA.clinical)
clin<-survivalTCGA(BRCA.clinical)
library(stringr)
library(survival)
library(survminer)
tnbcTID <-as.data.frame(substring(tnbc_tcga_id$`colnames(TNBC_Exp)`,1,12))
tnbctt<-as.data.frame(substring(colnames(TNBC_Exp),1,16))
clin <-dplyr::inner_join(clin,TNBC_case,by="bcr_patient_barcode")
BRCA_FPKM <-as.data.frame(t(BRCA_dataEXP_FPKM))
BRCA_FPKM <-dplyr::mutate(BRCA_FPKM,bcr_patient_barcode = substring(rownames(BRCA_FPKM), 1, 16))
BRCA_FPKM <-dplyr::select(BRCA_FPKM,bcr_patient_barcode,PBK,TOP2A,CDCA8,ASPM,CCNA2,KIF20A,BUB1,AURKB,CDK1,CCNB2)
TNBC_SUR<-dplyr::inner_join(tnbctt,BRCA_FPKM,by="bcr_patient_barcode")
TNBC_SUR$bcr_patient_barcode<-substring(TNBC_SUR$bcr_patient_barcode,1,12)
TNBC_SURT<-dplyr::inner_join(clin,TNBC_SUR,by="bcr_patient_barcode")

mysurv <-Surv(TNBC_SURT$times,TNBC_SURT$patient.vital_status)
lgrp<-apply(TNBC_SURT[,4:13],2,function(v1){
  group=ifelse(v1>median(v1),'high','low')
  kmfit2<-survfit(mysurv~group,data=TNBC_SURT)
  #plot(kmfit2)
  data.survdiff=survdiff(mysurv~group)
  pv= 1 - pchisq(data.survdiff$chisq,length(data.survdiff)-1)
})

group <-ifelse(TNBC_SURT$TOP2A > median(TNBC_SURT$TOP2A) 
               & TNBC_SURT$KIF20A > median(TNBC_SURT$KIF20A)
               & TNBC_SURT$CDCA8 > median(TNBC_SURT$CDCA8)
               & TNBC_SURT$CCNB2 > median(TNBC_SURT$CCNB2)
               ,'high','low')
sfit<-survfit(Surv(times,patient.vital_status)~group,data = TNBC_SURT)
ggsurvplot(sfit,conf.int = F,pval = TRUE,risk.table =T,
           ,xlab="Time(days)",title="TOP2A,KIF20A,CDCA8,CCNB2")


#####GO & KEGG(TOP10)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
testmodgeneX = bitr(top10$X1,fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
egoX <- enrichGO(gene=testmodgeneX$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
go_ALLX <-as.data.frame(egoX)
dotplot(egoX, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
