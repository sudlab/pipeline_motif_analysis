library(tidyverse)
library("optparse")

#Options parser
option_list = list(
  make_option(c("-f", "--fire-input"),
              type="character",
              dest = "fire_path",
              help="full path directory of fire motif analysis output,
                    .signif.motifs.rep")
)

arguments <- parse_args(OptionParser(option_list = option_list))

#arguments <- data.frame(fire_path = "/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/hepg2_motif_analysis/one_transcript/fire_halflife.txt_FIRE/RNA/fire_halflife.txt.signif.motifs.rep")

outdir <- str_remove(arguments$fire_path, ".signif.motifs.rep")

title <- str_split(outdir, ".txt")

outdir <- paste(title[[1]][1], ".txt", title[[1]][2], title[[1]][3], sep = "")


#Processing FIRE output
signif.motifs.rep <- read.delim2(arguments$fire_path,
                                 head = F,
                                 #colClasses = c("character", rep("numeric", 7)),
                                 )
#Have to do it like that cause it doesn't work with colClasses, I don't know why...
signif.motifs.rep$V5 <- as.numeric(signif.motifs.rep$V5)
signif.motifs.rep$V6 <- as.numeric(signif.motifs.rep$V6)
signif.motifs.rep$V7 <- as.numeric(signif.motifs.rep$V7)


#col V7 is proportion of genes in that bin that have that motif
#So selecting the 2 bin with the highest proportion
signif.motifs.rep %>%
  group_by(V1) %>%
  slice_max(V7, n = 2) -> highest.motif


#Dividing motif in high stab bin or low stab bin
#Taking top and bottom 10 bins (corresponds to more or
#less 20% of top and bottom genes)
lowbin = round(max(signif.motifs.rep$V2) * 0.2)
topbin = round(max(signif.motifs.rep$V2) * 0.8)

#High stab
  highest.motif %>% filter(V2 >= topbin) %>%
    filter(V7 >= 0.25) %>%
  select(V1) %>% distinct() %>% pull(V1)  %>%
  str_remove_all("\\.")   %>%
    str_replace_all("T", "U") %>%
    as.data.frame() %>%
  write.table(paste0(outdir, "_highstab.signif.motifs"),
              col.names = F, row.names = F, quote = F)


#Low stab
highest.motif %>% filter(V2 <= lowbin) %>%
  filter(V7 >= 0.25) %>%
  select(V1) %>% distinct() %>% pull(V1)  %>%
  str_remove_all("\\.")  %>%
    str_replace_all("T", "U") %>%
    as.data.frame() %>%
  write.table(paste0(outdir, "_lowstab.signif.motifs"),
              col.names = F, row.names = F, quote = F)
