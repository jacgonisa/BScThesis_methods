### Generate working space
mkdir genome  results samples 
cp $GENOME genome/knitens.fa  #copy the reference genome into the respective directory
cd samples
mkdir IP INPUT Mock
cp K27-IP-Klebs.fastq.gz IP/K27-IP-Klebs.fastq.gz
cp K27-INPUT-Klebs.fastq.gz INPUT/K27-INPUT-Klebs.fastq.gz
cp K27-Mock-Klebs.fastq.gz Mock/K27-Mock-Klebs.fastq.gz

### Workspace created


### Create index
cd ../genome
bowtie2-build knitens.fa index


### Processing samples
# Mapping
cd ../samples/IP
fastqc K27-IP-Klebs.fastq.gz
bowtie2 -x ../../../genome/index -U K27-IP-Klebs.fastq.gz-S IP.sam 

cd ../INPUT
fastqc K27-INPUT-Klebs.fastq.gz
bowtie2 -x ../../../genome/index -U K27-INPUT-Klebs.fastq.gz-S INPUT.sam 

# Convert SAM to BAM and create index
samtools sort -o INPUT.bam INPUT.sam
rm INPUT.sam
samtools index INPUT.bam

cd ../IP
samtools sort -o IP.bam IP.sam
rm IP.sam
samtools index IP.bam

# Convert BAM to BED and peak calling
bedtools bamtobed -i IP.bam > IP.bed
cd ../INPUT
bedtools bamtobed -i INPUT.bam > INPUT.bed
cd ../../results
epic2 --treatment ../samples/IP/IP.bed --control ../samples/INPUT/INPUT.bed > IPvsINPUT.bed
