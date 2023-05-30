
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' IPvsINPUT.bed > IPvsINPUT.bedgraph
sort -k1,1 -k2,2n IPvsINPUT.bedgraph > IPvsINPUT_sorted.bedgraph
bedGraphToBigWig IPvsINPUT_sorted.bedgraph myChrom.sizes IPvsINPUT.bw


computeMatrix  --scale-regions  -b 3000 -a 3000 -R genebody.bed -S IPvsINPUT.bw -o matrix.gz  --outFileSortedRegions regions_H3K27me3.bed  --regionBodyLength 5000  -o matrix.gz --outFileNameMatrix matrix.tab --outFileSortedRegions matrix_sorted_regions.bed
plotHeatmap -m ../matrix.gz --outFileName HEATMAP.png --colorMap RdBu  --heatmapWidth 10 --heatmapHeight 10 --plotTitle "H3K27me3 gene body enrichment" --kmeans 3
