setwd("supplementary-data/")
chr_pair_interactions <- read.table("chr_pair_interactions.csv")
total_num_reads <- sum(chr_pair_interactions$V3)

intra_chromosomal_data = list()
intra_chromosomal_data_names = list()
setwd("../10kb_resolution_intrachromosomal/")
list_dir <- paste(list.dirs(".", recursive = F),"/MAPQG0/", sep = "")
for (d in 1:length(list_dir)){
  chr_file <- list.files(list_dir[d], pattern = "\\.RAWobserved$")
  chr_dir <- paste(list_dir[d], chr_file, sep = "")
  intra_chromosomal_data[[d]] <- read.table(chr_dir)
  intra_chromosomal_data_names[[d]] <- sub("\\..*", "", chr_file)
}

rm(list_dir, d, chr_file, chr_dir)

total_intra_reads = sum(chr_pair_interactions[chr_pair_interactions$V1 == chr_pair_interactions$V2, "V3"])
total_enter_reads = total_num_reads - total_intra_reads

prob_vector = vector()
for (d in 1:length(intra_chromosomal_data)){
  prob_vector = c(prob_vector, as.integer(intra_chromosomal_data[[d]]$V3))
}
prob_vector = c(prob_vector, as.integer(total_enter_reads))

down_sampled_vector = rmultinom(1, total_num_reads/16, prob_vector)
down_sampled_vector = as.vector(down_sampled_vector)

setwd("../10kb_resolution_intrachromosomal_down/")
cnt = 0
for (d in 1:length(intra_chromosomal_data)){
  down_sampled_intra_chromosomal_data = data.frame(intra_chromosomal_data[[d]]$V1,
                                                        intra_chromosomal_data[[d]]$V2,
                                                        down_sampled_vector[(cnt+1):(cnt+nrow(intra_chromosomal_data[[d]]))])
  file_name = paste(intra_chromosomal_data_names[[d]], "_down.RAWobserved", sep = "")
  write.table(down_sampled_intra_chromosomal_data, file = file_name, row.names = F, col.names = F)
  cnt = cnt + nrow(intra_chromosomal_data[[d]])
}


