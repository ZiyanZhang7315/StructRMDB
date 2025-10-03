suppressMessages(library(GenomicFeatures))
suppressMessages(library(BSgenome))
suppressMessages(library(stringr))
suppressMessages(library(Biostrings))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(AnnotationHub))
suppressMessages(library(GenomicRanges))
suppressMessages(library(RMAnno))
library(readr) 
library(dplyr)

source("matureExtract2.R")
##############################################################
site1 <-read.csv("Species.bed",sep = "", head = F) ##This is an example site information input
site1$V1 <- str_remove(site1$V1, "chr")
site1$V1 <- str_replace(site1$V1, "M", "MT")
site1 <- GRanges(seqnames = site1$V1,
                 strand = site1$V6,
                 ranges = IRanges(start = site1$V3, width = 1))
txdb_1 <- makeTxDbFromGFF("Species.annotation.gtf")  ##This is an example gtf input
genome_1 <- readDNAStringSet("genome.fa")  ##This is an example genome input
seqlevels(txdb_1) <- str_remove(seqlevels(txdb_1),"chr")
genome_1 <- (genome_1)[1:25]
seqlevels(genome_1) <- str_remove(seqlevels(genome_1),"chr")
seqlevels(genome_1) <- sub(".*\\s", "", seqlevels(genome_1))
seqlevels(genome_1)[seqlevels(genome_1)=="MT"] <- "M"
site1<- site1[seqnames(site1) %in% names(genome_1)]
seqlevels(site1) <-  seqlevels(site1)[1:24] 
list1 <- matureExtract(site1, txdb_1, genome_1)

Grange1<- list1[[1]]
exon_seq1<- list1[[2]]
list1_paste<- paste0(seqnames(Grange1), "_", Grange1$site, "_", strand(Grange1))
list1_paste<- data.frame(list1_paste)
colnames(list1_paste) <- "location"

# read information
Species<- read.csv("Species.bed",sep = "", head = F)
Species$V1 <- str_remove(Species$V1,"chr")
Species$V1 <- str_replace(Species$V1, "M", "Mt")
Species<-data.frame(Species$V4,paste0(Species$V1, "_", Species$V3, "_", Species$V6))
colnames(Species_info) <- c("ID","location")
merged1 <- left_join(list1_paste, Species_info, by = "location")
Grange1$ID<- merged1$ID

###write the exon information and sequence
write_rds(Grange1,"Species_exon")
writeXStringSet(exon_seq1,"Species_exon_seq")