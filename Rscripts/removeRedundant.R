library(tidyverse)
library("optparse")
library("purrr")
library(msa)

#Options parser
option_list = list(
  make_option(c("-f", "--fire-input"),
              type="character",
              dest = "fire_input",
              help="Merged fire kmers motifs"),
  make_option(c("-o", "--output-name"),
              type="character",
              dest = "output",
              help="Name of output, full path from .")
)

arguments <- parse_args(OptionParser(option_list = option_list))

# setwd("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/a549_slam//motif_analysis/all_transcripts/")
# arguments <- data.frame(fire_input = "fire.dir/halflife_highstab.allkmer.signif.motifs",
#                         output = "fire.dir/lowstab.allkmer.fireMotifs")

my_msa <- function (x, my_method) {
  seq_msa <- msa(BStringSet(x),method = my_method,
                 type = "rna")
  consensus <- msaConsensusSequence(seq_msa)
  return(consensus)
}

inputs <- str_split(arguments$fire_input, pattern = ",")[[1]]

#Only one motif file
if (length(inputs) ==1) {
    fireMotifs <- read.delim(inputs, header = F)
  } else {
    fireMotifs <- lapply(inputs, read.delim, header = F) %>%
    purrr::reduce(full_join)
}

#Check for redundancy, cluster similar motifs, keep "consensus"

motif_vector <- fireMotifs[,2]

self_match <- sapply(motif_vector, function(i) motif_vector[str_detect(motif_vector, i)])

if (class(self_match) != "list") {
  #No match
  write.table(fireMotifs, arguments$output, 
              col.names = F, row.names = F, quote = F, sep = "\t")
} else {
  
final <- self_match[lapply(self_match, length) == 1]

several <- self_match[lapply(self_match,length) >1]

if ( all(lapply(several,  length) ==  2)       ) {
  #No need to cluster, keep final, it keeps the longest
  final_motifs <- final %>% do.call(c, .)
  Fire_final <- fireMotifs[fireMotifs$V2 %in% final_motifs,]
  write.table(Fire_final, arguments$output, 
              col.names = F, row.names = F, quote = F, sep = "\t")
  
} else {
  #Cluster them
clustered <- lapply(several, 
                    function(x) unique(unlist(several[sapply(several, 
                                                             function(y) any(x %in% y))]))) %>%
  unique()
msa_consensus <- lapply(clustered, function(x) my_msa(x, my_method = "ClustalOmega") )
list_several <- do.call(c, msa_consensus) %>% 
  str_remove_all("[^AUCG]") 

to_filter_in_final <- clustered %>% do.call(c,.)

final <- final[!(final %in% to_filter_in_final)]

final_motifs <- final %>% do.call(c, .) %>% c(., list_several) %>%
  unique()

Fire_final <- fireMotifs[fireMotifs$V2 %in% final_motifs,]

write.table(Fire_final, arguments$output, 
            col.names = F, row.names = F, quote = F, sep = "\t")
}

}