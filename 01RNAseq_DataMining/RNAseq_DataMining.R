### Download and readthe three available RNAseq studies
#Our study
matrix = read.table(gzfile("GSE198330_gene_expression.tsv.gz"), header = T)

#OPDA study
matrix_opda = read.table(gzfile("GSE141417_Normalized_counts.tab.gz"), header = T)

#Transcription heat factor study
heat2022 <- as.data.frame(read_xlsx("GSE178776_3_raw_reads.xlsx"))

###The format of gene names is different among experiments, we should make it all the same:
#
heat2022$...1 <- gsub("gene-KFL_","kfl", heat2022$...1) 
heat2022$...1 <- gsub("kfl_","kfl", heat2022$...1) 

heat2022$...1 <- sub("(.{8})(.*)", "\\1_\\2", heat2022$...1 )

#
opda_genes <- rownames(matrix_opda)
matrix_opda$geneID = opda_genes
matrix_opda$geneID <- gsub("_v1.1","", matrix_opda$geneID)


###Merge all the studies in one table and save it
merged_2 = merge(matrix, heat2022, by=1)
merged_3 = merge(merged_2, matrix_opda, by ="geneID")

write.table(merged_3, "Studies_merged.txt")

### Get subgroups of unexpressed genes for each study 
merged_all = read.table("Studies_merged.txt")
MY_NO_Expressed_genome = cbind(merged_all[1],merged_all[2:5])
MY_NO_Expressed_genome = MY_NO_Expressed_genome[merged_all$low_light_1 <= 0 & merged_all$low_light_2 <= 0 & 
                                  merged_all$high_light_1 <= 0  & merged_all$high_light_2 <= 0, ]

HEAT_NO_Expressed_genome = cbind(merged_all[1],merged_all[6:11])
HEAT_NO_Expressed_genome = HEAT_NO_Expressed_genome[merged_all$c1 <= 0 & merged_all$c2 <= 0 & 
                                                      merged_all$c2 <= 0 & merged_all$c3 <= 0
                                                    & merged_all$hs1 <= 0 & merged_all$hs2 <= 0, ]

OPDA_NO_Expressed_genome = cbind(merged_all[1],merged_all[12:17])
OPDA_NO_Expressed_genome = OPDA_NO_Expressed_genome[merged_all$OPDA1_kleb.htseq <= 0 &
                                                      merged_all$OPDA2_kleb.htseq <= 0 & merged_all$OPDA4_kleb.htseq <= 0 & 
                                                      merged_all$Mock1_kleb.htseq <= 0  & merged_all$Mock2_kleb.htseq <= 0
                                                    & merged_all$Mock3_kleb.htseq <= 0, ]
   
   
### Intersect the subgroups with a venn diagram
library(ggVennDiagram)
library(ggplot2)

x=list(MY_NO_Expressed_genome$geneID, HEAT_NO_Expressed_genome$geneID, OPDA_NO_Expressed_genome$geneID)
ggVennDiagram(x, label_alpha = 0, color= "black", lwd =1, lty=1, category.names = 
                c("A","B", 
                  "C"), stroke_size = 1) + ggplot2::scale_fill_gradient(low = "#F4FAFE", 
                                                                        high = "#4981BF") + theme(legend.position = "right")
