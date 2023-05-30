
#Installing the annotation for K. nitens
library(devtools)
install_github("fran-romero-campero/AlgaeFUN/packages/txdb_packages/TxDb.Knitens.Phycocosm")
library(TxDb.Knitens.Phycocosm)
install.packages("./org.Knitens.eg.db", repos = NULL, type = "source")
library(org.Knitens.eg.db)

txdb <- TxDb.Knitens.Phycocosm

#Read peak file obtained after peak calling
library(ChIPseeker)
library(clusterProfiler)
library(ChIPpeakAnno)

peaks <- readPeakFile(peakfile = "IPvsMock.bed", header=F) 

#Annotate peaks
peakAnno <- annotatePeak(peak=peaks, tssRegion=c(-1000, 100),
                         TxDb=txdb, annoDb ="org.Knitens.eg.db", ignoreOverlap = F)

Anno= as.data.frame(peakAnno)
colnames(Anno)[19] ="geneID" 
sum(Anno$width)/103921766 #% of genome covered by H3K27me3

genes <- getBioRegion(TxDb = txdb,
                         by = "gene")
df_genes <- data.frame(seqnames=seqnames(genes),
                       start=start(genes) - 1,
                       end=end(genes),
                       names=c(rep(".", length(genes))),
                       scores=c(rep(".", length(genes))),
                       strands=strand(genes),genes$gene_id)
length(unique(Anno$geneID))/length(unique(df_genes$genes.gene_id))  #% of genes covered by H3K27me3

#Create Figure 10A
peakLengths= end(peaks) - start(peaks)+1

hist(peakLengths, breaks = 15, col = "#4981BF", 
     border = "white", xlab = "Peak size (bp)", 
     ylab = "Number of peaks", main=NULL)
med <- median(peakLengths)
abline(v = med, col = "red", lwd = 2)

#Create Figure 10B
  #Read the table obtained during RNA-seq data mining
merged=read.table("C:/Users/jacob/OneDrive/Escritorio/TFG/klebso/03no_expresado_genoma/Studies_merged.txt")

MY_NO_Expressed_genome = merged[merged$low_light_1 <= 0 & merged$low_light_2 <= 0, ]
MY_NO_Expressed_genes= MY_NO_Expressed_genome$geneID
matrix_cistrome= matrix[matrix$geneID%in%Anno$geneID, ] #matrix file from RNA_seq data mining
cistrome_noexpressedgenes= list(matrix_cistrome$geneID, MY_NO_Expressed_genes)
library(ggVennDiagram)
ggVennDiagram(cistrome_noexpressedgenes, label_alpha = 0, category.names = c("H3K27me3 target genes","Unexpressed genes")
) + ggplot2::scale_fill_gradient(low = "#F4FAFE", 
                                 high = "#4981BF")
#Testing the statistical significance of the enrichment of repressed genes in H3K27me3
library("GeneOverlap")
go.obj <- newGeneOverlap(matrix_cistrome$geneID,
                         MY_NO_Expressed_genes,
                         genome.size=17290)
go.obj <- testGeneOverlap(go.obj)
print(go.obj) #p-value

### Get genebody.bed 
#This .bed file is to compute the Matrix with deepTools
gr <- getBioRegion(TxDb = txdb,
                         by = "gene", type = "body")

df <- data.frame(seqnames=seqnames(gr),
                 starts=start(gr)-1,
                 ends=end(gr),
                 names=c(rep(".", length(gr))),
                 scores=c(rep(".", length(gr))),
                 strands=strand(gr))

write.table(df, file="genebody.bed", quote=F, sep="\t", row.names=F, col.names=F)


