#########################This file is used to extrace mature RNA sequence using the exon sequence
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

Species<- read_rds("Species_exon")   #####This is an example exon input
Species_exon_seq <-readDNAStringSet("Species_exon_seq")  #####This is an example exon sequence input
Species_exon <- GRanges(seqnames = Species$tx_name,IRanges(start=as.integer(Species$siteInMature),end =as.integer(Species$siteInMature),width = 1), strand= strand(Species), ID =Species$ID)
strand(Species_exon) <- "+" 
sum(letterFrequency(getSeq(Species_exon_seq,Species_exon),"A")) #check
index<- which(letterFrequency(getSeq(Species_exon_seq,Species_exon),"A") ==1 )
not_Species<- Species_exon[-index]
Species_exon<- Species_exon[index]
write_rds(not_Species,"notA_Species.rds")

location<-data.frame(seq_start=1,seq_end=width(Species_exon_seq),tx_id =names(Species_exon_seq)) 
Species_exon_frame<- data.frame(Species_exon)
names(Species_exon_frame)[1] <- "tx_id"
merged_df <- merge(Species_exon_frame, location, by = "tx_id")
merged_df <- merged_df[order(merged_df$ID),]
merged_df$m6A_location <- pmin(201,merged_df$start)
merged_df$start <- pmax(merged_df$seq_start,merged_df$start-200)
merged_df$end <- pmin(merged_df$seq_end,merged_df$end+200)
merged_df$width <- merged_df$end-merged_df$start+1

#only pick the longest one
merged_m6A <- merged_df %>%
  group_by(ID) %>%
  filter(width == max(width)) %>%
  slice(1) %>%
  ungroup()

merged_m6A$ID<- paste(merged_m6A$ID,"m",sep="")
m6A_location<- as.vector(merged_m6A[,9])
names(merged_m6A)[1] <- "seqnames"
Grange_m6A<- GRanges(data.frame(merged_m6A[,1:6]))
#out put transcriptom name
write.table(merged_m6A[,c(1,6)],"mat_Species_trans_name.txt", sep = "\t",row.names = FALSE)
extract<- getSeq(Species_exon_seq,Grange_m6A)
sum(letterFrequency(subseq(extract,start =m6A_location$m6A_location,end =m6A_location$m6A_location),"A")) == length(extract) #check
subseq(extract,start =m6A_location$m6A_location,end =m6A_location$m6A_location) <- "M"
names(extract)<- Grange_m6A$ID

indices_with_N <- grep("N", extract)
indices_with_N
m6A_sequence <- extract

####write out the sequence
m6A_sequence <- RNAStringSet(gsub("T", "U", m6A_sequence))
writeXStringSet(m6A_sequence, "mat_Species_sequence_modification.fa")
m6A_sequence <- RNAStringSet(gsub("M", "A", m6A_sequence))
writeXStringSet(m6A_sequence, "mat_Species_sequence.fa")
writeXStringSet(m6A_sequence_withN, "mat_Species_sequence_withn.fa")
merged_m6A_location<- merged_m6A[,c(6,9)]
write.table(merged_m6A_location,"mat_Species__modification_location.txt", sep = "\t",row.names = FALSE)


