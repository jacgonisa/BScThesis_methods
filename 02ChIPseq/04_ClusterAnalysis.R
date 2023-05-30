#Continue from 03, with the output of sorting clusters

regions <- read.table("3sorted_clusters.txt", header = F)
cluster1 <- regions[regions$V13 == "C1", ]
cluster2 <- regions[regions$V13 == "C2", ]

#Analyze mean FPKM in cluster 1
colnames(cluster1)[1]="seqnames"
colnames(cluster1)[2]="start"
colnames(cluster1)[3]="end"
cluster1_df_genes= merge(cluster1, df_genes, by=c("start","end","seqnames") )

cluster1_df_genes_tx =merge(cluster1_df_genes, df_transcript, by=c("start","end","seqnames") )
colnames(cluster1_df_genes_tx)[21]= "transcriptId"
cluster1_df_genes_tx_Annot = merge(cluster1_df_genes_tx, AnnoWithExpression, by="transcriptId" )
mean(cluster1_df_genes_tx_Annot$low_light_1)
sd(cluster1_df_genes_tx_Annot$low_light_1)/sqrt(length(cluster1_df_genes_tx_Annot$low_light_1))

nrow(cluster1)
nrow(cluster1_df_genes_tx)

write.table(cluster1_df_genes_tx$genes.gene_id,row.names = F, "cluster1.txt")

#Analyze mean FPKM in cluster 2
cluster2=read.table("../06ChIP_klebs/NOSKIPcluster2.txt", header = T)
nrow(cluster2)

colnames(cluster2)[1]="seqnames"
colnames(cluster2)[2]="start"
colnames(cluster2)[3]="end"
cluster2_df_genes= merge(cluster2, df_genes, by=c("start","end","seqnames") )

cluster2_df_genes_tx =merge(cluster2_df_genes, df_transcript, by=c("start","end","seqnames") )
colnames(cluster2_df_genes_tx)[21]= "transcriptId"
cluster2_df_genes_tx_Annot = merge(cluster2_df_genes_tx, AnnoWithExpression, by="transcriptId" )
mean(cluster2_df_genes_tx_Annot$low_light_1)
sd(cluster2_df_genes_tx_Annot$low_light_1)/sqrt(length(cluster2_df_genes_tx_Annot$low_light_1))
mean(matrix$low_light_1)
sd(matrix$low_light_1)/sqrt(length(matrix$low_light_1))

write.table(cluster2_df_genes_tx$genes.gene_id,row.names = F, "Cluster2NOSKIP.txt")
