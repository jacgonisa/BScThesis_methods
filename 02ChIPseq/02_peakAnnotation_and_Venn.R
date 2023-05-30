
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
peakLengths= end(peaks) - start(peaks)+1

#Crea
hist(peakLengths, breaks = 15, col = "#4981BF", 
     border = "white", xlab = "Peak size (bp)", 
     ylab = "Number of peaks", main=NULL)
### Get genebody.bed 



### Venn regulome
