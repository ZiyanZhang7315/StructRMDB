library(RNAsmc)
library(parallel)

############################################################## RNAsmc scoring
folder1_path <- "CTm"
folder2_path <- "CTa"
dir<- getwd()
dir1<- paste(dir,"/CTm/",sep = "")
dir2<- paste(dir,"/CTa/",sep = "")

compare_files <- function(file1, file2) {
  data1 <- read.csv(paste(dir1, file1, sep = ""), sep = "", row.names = NULL)
  data2 <- read.csv(paste(dir2, file2, sep = ""), sep = "", row.names = NULL)
  
  result <- strCompare(data1, data2, randomTime = 2)
  smcscore <- result[[1]]
  p_value_smc <- result[[2]]
  
  return(c(smcscore, p_value_smc))
}

files_in_folder1 <- list.files(path = folder1_path)
files_in_folder2 <- list.files(path = folder2_path)

mc <- getOption("mc.cores", 5)
results <- t(mclapply(1:length(files_in_folder1), function(i) {
  compare_files(files_in_folder1[i], files_in_folder2[i])
},mc.cores = mc)) 


smcscore <- sapply(results, function(sublist) sublist[[1]])
p_value_smc <- sapply(results, function(sublist) sublist[[2]])

SMC <- data.frame(site = files_in_folder1,smcscore,p_value_smc)
SMC$site <- substring(SMC$site, 4, nchar(SMC$site))
SMC$site <- gsub("\\s", "", SMC$site)

############################################################## These code is used for combine the result
setwd("result")
similarity_score <- readLines("similarity_score.txt")
relative_score <- readLines("relative_score.txt")
forester_distance <- readLines("forester_distance.txt")
names<- readLines("site_name.txt")
forester <- data.frame(site = names,similarity_score,relative_score,forester_distance)


MFE_m6A<- read.csv("MFE_structure_m6A_R.txt",sep = " ",head = FALSE)
MFE<- read.csv("MFE_structure_R.txt",sep = " ",head = FALSE)
MFE_m6A<- MFE_m6A[-2]
MFE<- MFE[-2]
colnames(MFE_m6A) <- c("MFE_m6A", "site")
colnames(MFE) <- c("MFE", "site")
MFE_C<- merge(MFE_m6A,MFE,by = "site")

data_df <- merge(merge(MFE_C, forester, by = "site", all = TRUE), SMC, by = "site", all = TRUE)
write.table(data_df,file = "data.txt", sep = "\t", row.names = TRUE)
