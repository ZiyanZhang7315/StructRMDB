#####This script is used to extract modification site surrounding sequences for species
library(BSgenome)  
library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)
library(AnnotationHub)
library(GenomicRanges)

####load site information and genome
Species<- read.csv("Species.uniqID.anno.score.conserved.guide.resort.col25.bed",sep = "",head =F) ####This is an example bed file
genome <-readDNAStringSet("primary_assembly.genome.fa")  ####This is an example genome file
modification<- GRanges(seqnames=Species$V1, range = IRanges(start = Species$V3, end = Species$V3), strand = Species$V6,ID = Species$V4)
gtf_file <- ExtractGTFFile("Species.annotation.gtf") ####This is an example GTF file
txdb <- makeTxDbFromGFF("Species.annotation.gtf")  ####This is an example GTF file
modification<- subset(modification,seqnames(modification) %in% intersect(seqlevels(modification),gtf_file$seqlevels))

seqlevels(modification)<- seqlevels(modification)[1:21]
modification<- rmAnno(modification, gtf_file)

seqlevels(modification) <- str_remove(seqlevels(modification),"chr")
seqlevels(txdb) <- str_remove(seqlevels(txdb),"chr")
transcript<- transcripts(txdb) 
genome <- (genome)[1:22]
seqlevels(genome)<- seqlevels(genome)[1:22]
seqlevels(genome) <- str_remove(seqlevels(genome),"chr")
seqlevels(genome) <- sub(".*\\s", "", seqlevels(genome))
seqlevels(genome)[seqlevels(genome)=="MT"] <- "M"

getSeq(genome,modification)
A_index<- which(letterFrequency(getSeq(genome,modification),"A") ==1)
modification<- modification[A_index]
write_rds(modification,"Species_basic_info.rds")
########################################
modification[which(is.na(str_detect(modification$Region, "Exon")))]$Region <- "Intergenic"
pre_rna <-  modification[modification$Region != "Intergenic"]
mapped_transcripts <- mapToTranscripts(pre_rna, txdb) 
tx_name<- data.frame(transcript$tx_id,transcript$tx_name)
names(tx_name) <- c("tx_id","tx_name")
transcript_sequences <- getSeq(genome, transcript)
names(transcript_sequences) <- transcript$tx_id
#replace Xhit to ID
xHits<- mapped_transcripts$xHits
modification_ID<- data.frame(pre_rna$ID)
modification_Hits<- modification_ID[xHits,]
mapped_transcripts$xHits<- modification_Hits

######extract sequence
strand(mapped_transcripts) <- "+" 
#site<- getSeq(transcript_sequences,mapped_transcripts)
location<-data.frame(seq_start=1,seq_end=width(transcript_sequences),tx_id =names(transcript_sequences)) 
mapped_transcripts_frame<- data.frame(mapped_transcripts)
names(mapped_transcripts_frame)[7] <- "tx_id"
merged_df <- merge(mapped_transcripts_frame, location, by = "tx_id")
merged_df <- merge(merged_df,tx_name, by= "tx_id")
merged_df <- merged_df[order(merged_df$xHits), ]
merged_df$modification_location <- pmin(201,merged_df$start)
merged_df$start <- pmax(merged_df$seq_start,merged_df$start-200)
merged_df$end <- pmin(merged_df$seq_end,merged_df$end+200)
merged_df$width <- merged_df$end-merged_df$start+1

# keep only the longest one
merged_modification <- merged_df %>%
  group_by(xHits) %>%
  filter(width == max(width)) %>%
  slice(1) %>%
  ungroup()
merged_modification$xHits<- paste(merged_modification$xHits,"p",sep="")
modification_location<- as.vector(merged_modification[,11])
Grange_modification<- GRanges(data.frame(merged_modification[,1:7]))

######extract sequence
extract<- getSeq(transcript_sequences,Grange_modification)
sum(letterFrequency(subseq(extract,start =modification_location$modification_location,end =modification_location$modification_location),"A")) == length(extract)
subseq(extract,start =modification_location$modification_location,end =modification_location$modification_location) <- "M"
names(extract)<- Grange_modification$xHits
indices_with_N <- grep("N", extract)
modification_sequence<- extract

#output sequence
modification_sequence <- RNAStringSet(gsub("T", "U", modification_sequence))
writeXStringSet(modification_sequence, "pre_Species_sequence_modification.fa")
modification_sequence <- RNAStringSet(gsub("M", "A", modification_sequence))
writeXStringSet(modification_sequence, "pre_Species_sequence.fa")
writeXStringSet(modification_sequence_withN, "Species_sequence_withn.fa")

merged_modification<- merged_modification[merged_modification$xHits %in% names(modification_sequence),]
merged_modification_location<- merged_modification[,c(7,11)]
write.table(merged_modification_location,"Species_modification_location.txt", sep = "\t",row.names = FALSE)
#output transcript name
write.table(merged_modification[,c(7,10)],"Species_trans_name.txt", sep = "\t",row.names = FALSE)

