regions <- read.table("H3K27me3.txt", header = F)

CLUSTER1_noskip=read.table("../06ChIP_klebs/NOSKIPcluster1.txt", header = T)


colnames(CLUSTER1_noskip)[1]="seqnames"
colnames(CLUSTER1_noskip)[2]="start"
colnames(CLUSTER1_noskip)[3]="end"
CLUSTER1_noskip_df_genes= merge(CLUSTER1_noskip, df_genes, by=c("start","end","seqnames") )

CLUSTER1_noskip_df_genes_tx =merge(CLUSTER1_noskip_df_genes, df_transcript, by=c("start","end","seqnames") )
colnames(CLUSTER1_noskip_df_genes_tx)[21]= "transcriptId"
CLUSTER1_noskip_df_genes_tx_Annot = merge(CLUSTER1_noskip_df_genes_tx, AnnoWithExpression, by="transcriptId" )
mean(CLUSTER1_noskip_df_genes_tx_Annot$low_light_1)
sd(CLUSTER1_noskip_df_genes_tx_Annot$low_light_1)/sqrt(length(CLUSTER1_noskip_df_genes_tx_Annot$low_light_1))

nrow(CLUSTER1_noskip)
nrow(CLUSTER1_noskip_df_genes_tx)

write.table(CLUSTER1_noskip_df_genes_tx$genes.gene_id,row.names = F, "Cluster1NOSKIP.txt")

##
CLUSTER2_noskip=read.table("../06ChIP_klebs/NOSKIPcluster2.txt", header = T)
nrow(CLUSTER2_noskip)
getwd()


colnames(CLUSTER2_noskip)[1]="seqnames"
colnames(CLUSTER2_noskip)[2]="start"
colnames(CLUSTER2_noskip)[3]="end"
CLUSTER2_noskip_df_genes= merge(CLUSTER2_noskip, df_genes, by=c("start","end","seqnames") )

CLUSTER2_noskip_df_genes_tx =merge(CLUSTER2_noskip_df_genes, df_transcript, by=c("start","end","seqnames") )
colnames(CLUSTER2_noskip_df_genes_tx)[21]= "transcriptId"
CLUSTER2_noskip_df_genes_tx_Annot = merge(CLUSTER2_noskip_df_genes_tx, AnnoWithExpression, by="transcriptId" )
mean(CLUSTER2_noskip_df_genes_tx_Annot$low_light_1)
sd(CLUSTER2_noskip_df_genes_tx_Annot$low_light_1)/sqrt(length(CLUSTER2_noskip_df_genes_tx_Annot$low_light_1))
mean(matrix$low_light_1)
sd(matrix$low_light_1)/sqrt(length(matrix$low_light_1))

write.table(CLUSTER2_noskip_df_genes_tx$genes.gene_id,row.names = F, "Cluster2NOSKIP.txt")
