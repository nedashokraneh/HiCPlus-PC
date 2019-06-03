num_interaction_matrix <- matrix(nrow = 23, ncol = 23)


setwd("10kb_resolution_interchromosomal/")
list_dir <- list.dirs(recursive = F)
for (dir in list_dir){
  chr_pair = unlist(strsplit(substr(dir,3,30), split="_"))
  first_chr = sub("chr", "", chr_pair[1])
  second_chr = gsub("chr", "", chr_pair[2])
  if (first_chr == "X"){
    first_chr = 23
    second_chr = as.integer(second_chr)
  } else if (second_chr == "X"){
    second_chr = 23
    first_chr = as.integer(first_chr)
  } else{
    first_chr = as.integer(first_chr)
    second_chr = as.integer(second_chr)
  }
  interaction_file = list.files(paste(dir, "/MAPQG0/", sep = ""), pattern = "\\.RAWobserved$")
  print(interaction_file)
  interaction_file = paste(dir, "/MAPQG0/", interaction_file, sep = "")
  print(interaction_file)
  interaction_data <- read.table(interaction_file)
  total_num_interactions = sum(interaction_data$V3)
  num_interaction_matrix[first_chr,second_chr] = total_num_interactions
}

setwd("../10kb_resolution_intrachromosomal/")
list_dir <- list.dirs(recursive = F)
for (dir in list_dir){
  chr_num = gsub("chr", "", substr(dir,3,20))
  if (chr_num == "X"){
    chr_num = 23
  }
  else {
    chr_num = as.integer(chr_num)
  }
  
  interaction_file = list.files(paste(dir, "/MAPQG0/", sep = ""), pattern = "\\.RAWobserved$")
  print(interaction_file)
  interaction_file = paste(dir, "/MAPQG0/", interaction_file, sep = "")
  print(interaction_file)
  interaction_data <- read.table(interaction_file)
  total_num_interactions = sum(interaction_data$V3)
  num_interaction_matrix[chr_num,chr_num] = total_num_interactions
}



