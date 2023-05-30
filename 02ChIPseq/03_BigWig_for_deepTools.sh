cd results
# Get chromosomes lengths for analysis
samtools faidx ../genome/knitens.fa
cut -f1,2 ../genome/knitens.fa.fai > 
cp ../genome/chromosome_lengths.txt .

bedClip IPvsINPUT.bed chromosome_lengths.txt bedClipped.bed #Filtering peaks out of the chromosomes range
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' bedClipped.bed > IPvsINPUT.bedgraph # Convert BED to BEDGRAPH
sort -k1,1 -k2,2n IPvsINPUT.bedgraph > IPvsINPUT_sorted.bedgraph #Sort BEDGRAPH
bedGraphToBigWig IPvsINPUT_sorted.bedgraph chromosome_lengths.txt IPvsINPUT.bw #Convert BEDGRAPH to BIGWIG


computeMatrix  --scale-regions  -b 3000 -a 3000 -R genebody.bed -S IPvsINPUT.bw -o matrix.gz  --outFileSortedRegions regions_H3K27me3.bed  --regionBodyLength 5000  -o matrix.gz --outFileNameMatrix matrix.tab --outFileSortedRegions matrix_sorted_regions.bed
plotHeatmap -m ../matrix.gz --outFileName HEATMAP.png --colorMap RdBu  --heatmapWidth 10 --heatmapHeight 10 --plotTitle "H3K27me3 gene body enrichment" --kmeans 3
