setwd("foldm")
system("sed -i '1s/.* //' *")
system("sed -i '1s/^/>/' *")

setwd("../folda")
system("sed -i '1s/.* //' *")
system("sed -i '1s/^/>/' *")

setwd("..")

folder1_path <- "foldm"
folder2_path <- "folda"
dir<- getwd()
dir1<- paste(dir,"/foldm/",sep = "")
dir2<- paste(dir,"/folda/",sep = "")

files_in_folder1 <- list.files(path = folder1_path)
files_in_folder2 <- list.files(path = folder2_path)


total_site<- read.csv("Species_modification_location.txt", sep="")  ###This is an example input
total_site$ID <- paste(total_site$xHits,".fold",sep = "")
files_in_folder <- substr(files_in_folder1, 4, nchar(files_in_folder1))
site_201<- total_site[total_site$ID %in% files_in_folder & total_site$m6A_location == 201,]
site_no201<- total_site[total_site$ID %in% files_in_folder & total_site$m6A_location != 201,]


CTm_201<- paste("CTm",site_201$ID,sep = "")
CTa_201<- paste("CTa",site_201$ID,sep = "")

write.table(CTm_201, file = "foldm/CTm_201.txt", sep = "\t", col.names = FALSE, row.names = FALSE,quote = FALSE)
write.table(CTa_201, file = "folda/CTa_201.txt", sep = "\t", col.names = FALSE, row.names = FALSE,quote = FALSE)

#foldm
setwd("foldm")
command1 <- "cat CTm_201.txt | awk -F'.' '{print $1}' > CTm1_201.txt"
system(command1, intern = TRUE)

command2 <- "cat CTm1_201.txt | parallel -j 5 'RNAplot --pre \"201 201 12 RED omark 201 cmark 201\" < {}.fold'"
system(command2, intern = TRUE)

command_no201 <- paste('RNAplot --pre "',site_no201$m6A_location ,sep = "",collapse = NULL)
command_no201 <- paste(command_no201,site_no201$m6A_location)
command_no201 <- paste(command_no201,'12 RED omark')
command_no201 <- paste(command_no201,site_no201$m6A_location)
command_no201 <- paste(command_no201,'cmark')
command_no201 <- paste(command_no201,site_no201$m6A_location)
command_no201 <- paste(command_no201,'" <')
command_no201 <- paste(command_no201,"CTm")
command_no201 <- paste0(command_no201,site_no201$ID)
lapply(command_no201, function(cmd) system(cmd, intern = TRUE))

system("mkdir ../photo_m")
system("mv *.ps ../photo_m")


#folda
setwd("../folda")
command1 <- "cat CTa_201.txt | awk -F'.' '{print $1}' > CTa1_201.txt"
system(command1, intern = TRUE)

command2 <- "cat CTa1_201.txt | parallel -j 5 'RNAplot --pre \"201 201 12 RED omark 201 cmark 201\" < {}.fold'"
system(command2, intern = TRUE)

command_no201 <- paste('RNAplot --pre "',site_no201$m6A_location ,sep = "",collapse = NULL)
command_no201 <- paste(command_no201,site_no201$m6A_location)
command_no201 <- paste(command_no201,'12 RED omark')
command_no201 <- paste(command_no201,site_no201$m6A_location)
command_no201 <- paste(command_no201,'cmark')
command_no201 <- paste(command_no201,site_no201$m6A_location)
command_no201 <- paste(command_no201,'" <')
command_no201 <- paste(command_no201,"CTa")
command_no201 <- paste0(command_no201,site_no201$ID)
lapply(command_no201, function(cmd) system(cmd, intern = TRUE))

system("mkdir ../photo_a")
system("mv *.ps ../photo_a")



