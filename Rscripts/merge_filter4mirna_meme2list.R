library("universalmotif")
library("optparse")
library("stringr")
library("tidyverse")

#Options parser
option_list = list(
  make_option(c("-i", "--input-meme"),
              type="character",
              dest = "input",
              help="high/lowstab input file in meme format"),
  make_option(c("-m", "--mirna-seeds"),
              type="character",
              dest = "miRNA_seeds",
              help="miRNA seeds target motifs in meme format"),
  make_option(c("-t", "--tomtom-file"),
              type="character",
              dest = "tomtom",
              help="tomtom results of highstab vs lowstab"),
  make_option(c("-f", "--fire-file"),
              type="character",
              dest = "fire",
              help="fire motifs"),
  make_option(c("-l", "--filter-fire"),
              type="character",
              dest = "filter_fire",
              help="fire motifs"),
  make_option(c("-o", "--output-file"),
              type="character",
              dest = "output",
              help="Name output file with path from ."),
  make_option(c("-r", "--relevant-miRNA"),
              type="character",
              dest = "relevant",
              help="File with names of relevant miRNAs for a cell line if available")
)

arguments <- parse_args(OptionParser(option_list = option_list))

#setwd("/mnt/sharc/shared/sudlab1/General/projects/SLAMseq_CHO_Mikayla/cho_slam/slam_picr_newUTRS/slam_highest_exons/motif_analysis/one_transcript")
#setwd("/mnt/sharc/shared/sudlab1/General/projects/SLAMseq_CHO_Mikayla/cho_slam/slam_picr_newUTRS/slam_highest_exons/motif_analysis/new_lasso/all_transcripts")
# setwd("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/hepg2_tripseq/motif_analysis/random/one_transcript")
# arguments <- data.frame(fire = "fire.dir/lowstab.allkmer.fireMotifs",
#                         input = "final_motifs/lowstab_merge_homer_streme.meme",
#                         miRNA_seeds = "/mnt/sharc/shared/sudlab1/General/mirror/meme_motif_db/motif_databases/MIRBASE/22/Homo_sapiens_hsa.seeds.meme",
#                         tomtom = "final_motifs/merge_homer_streme_highstab_vs_lowstab.tomtom",
#                         filter_fire = "fire.dir/lowstab_in_lowstab.list",
#                         output = "final_motifs/lowstab_final_motifs.list",
#                         relevant = "/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/miRmine/relevant_miRNA_a549_hepg2.tsv")

# arguments <- data.frame(fire = "/shared/sudlab1/General/projects/SLAMseq_CHO_Mikayla/cho_slam/slam_picr_newUTRS/slam_highest_exons/motif_analysis/one_transcript/fire.dir/highstab.allkmer.fireMotifs",
#                         input = "final_motifs/highstab_merge_homer_streme.meme",
#                         miRNA_seeds = "/shared/sudlab1/General/mirror/meme_motif_db/motif_databases/MIRBASE/22/Cricetulus_griseus_cgr.seeds.meme",
#                         tomtom = "final_motifs/merge_homer_streme_highstab_vs_lowstab.tomtom",
#                         filter_fire = "fire.dir/highstab_in_lowstab.listt",
#                         output = "final_motifs/highstab_final_motifs.list")





#Modify read_meme function so it accepts spaces in the file 
#Otherwise it throws an error
#raw_lines <- readLines(con <- file("final_motifs/highstab_final_motifs.meme"))
read_meme2 <- function (file, skip = 0, readsites = FALSE, readsites.meta = FALSE) 
{
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file), 1, 
                                 FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(skip = args$skip), 1, 
                                FALSE, TYPE_NUM)
  logi_check <- check_fun_params(list(readsites = args$readsites, 
                                      readsites.meta = args$readsites.meta), numeric(), logical(), 
                                 TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) 
    stop(all_checks_collapse(all_checks))
  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) 
    raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]
  raw_lines <- raw_lines[!grepl("^#", raw_lines)]
  raw_lines <- raw_lines[!grepl("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", 
                                raw_lines)]
  raw_lines <- raw_lines[!grepl("------------", raw_lines)]
  check_meme_version(raw_lines)
  alph <- get_meme_alph(raw_lines)
  alph <- "RNA"
  alph.len <- get_meme_alph_len(alph)
  alph.split <- switch(alph, DNA = DNA_BASES, RNA = RNA_BASES, 
                       AA = AA_STANDARD2, safeExplode(alph))
  strands <- raw_lines[grepl("^strands:", raw_lines)]
  if (length(strands) > 0) {
    strands <- strsplit(strands, "\\s+")[[1]][-1]
  }
  else {
    message("Could not find strand info, assuming +.")
    strands <- "+"
  }
  if (all(c("+", "-") %in% strands)) {
    strands <- "+-"
  }
  bkg.start <- grep("^Background letter frequencies", raw_lines)
  if (length(bkg.start)) {
    bkg.offset <- 1 
    bkg <- raw_lines[bkg.start + bkg.offset]
    bkg <- strsplit(bkg, "\\s+")[[1]]
    bkg <- bkg[!grepl("^$", bkg)] #adding this to remove space
    bkg <- as.numeric(bkg[seq_len(length(bkg))%%2 == 0])
    while (length(bkg) < alph.len) {
      bkg.offset <- bkg.offset + 1
      bkg.tmp <- raw_lines[bkg.start + bkg.offset]
      bkg.tmp <- strsplit(bkg.tmp, "\\s+")[[1]]
      bkg.tmp <- as.numeric(bkg.tmp[seq_along(bkg.tmp)%%2 == 
                                      0])
      bkg <- c(bkg, bkg.tmp)
    }
    if (anyNA(bkg)) 
      stop("Could not parse background frequencies, check that they match alphabet")
  }
  else {
    message("Could not find background, assuming uniform frequencies.")
    bkg <- rep(1/length(alph.split), length(alph.split))
  }
  motif_meta <- grep("^letter-probability matrix:", raw_lines)
  motif_names_i <- grep("^MOTIF ", raw_lines)
  motif_names <- lapply(raw_lines[motif_names_i], function(x) {
    x <- strsplit(x, "\\s+")[[1]]
    if (x[1] == "") 
      x[3]
    else x[2]
  })
  motif_altnames <- lapply(raw_lines[motif_names_i], function(x) {
    x <- strsplit(x, "\\s+")[[1]]
    if (x[1] == "") 
      x[4]
    else x[3]
  })
  motif_starts <- motif_meta + 1
  motif_stops <- sapply(raw_lines[motif_meta], function(x) strsplit(x, 
                                                                    "\\s+")[[1]][6])
  motif_stops <- motif_meta + as.numeric(motif_stops)
  motif_meta <- lapply(raw_lines[motif_meta], function(x) {
    x <- strsplit(x, "\\s+")[[1]]
    c(nsites = x[8], eval = x[10])
  })
  motif_list <- mapply(function(x, y) {
    z <- trimws(raw_lines[x:y])
    z <- sapply(z, function(x) strsplit(x, "\\s+")[[1]])
    if (nrow(z) != alph.len) 
      stop("Alphabet length does not match motif length")
    z <- z[order(alph.split, method = "radix"), ]
    as.numeric(z)
  }, motif_starts, motif_stops, SIMPLIFY = FALSE)
  motif_list <- mapply(function(x, y, z, x2) {
    mot <- universalmotif_cpp(name = x, type = "PPM", altname = x2, 
                              nsites = as.numeric(y[1]), eval = as.numeric(y[2]), 
                              bkg = bkg, alphabet = alph, strand = strands, extrainfo = c(eval.string = unname(y[2])), 
                              motif = t(matrix(z, ncol = alph.len, byrow = TRUE)))
    validObject_universalmotif(mot)
    mot
  }, motif_names, motif_meta, motif_list, motif_altnames, 
  SIMPLIFY = FALSE)
  if (length(motif_list) == 1) 
    motif_list <- motif_list[[1]]
  if (readsites) {
    if (is.list(motif_list)) 
      mot.names <- vapply(motif_list, function(x) x@name, 
                          character(1))
    else mot.names <- motif_list@name
    block.starts <- grep("in BLOCKS format", raw_lines)
    if (length(block.starts) == 0) {
      warning("could not find BLOCKS formatted motifs in MEME file")
      motif_list <- list(motifs = motif_list, sites = NULL)
    }
    else {
      block.len <- vapply(block.starts, function(x) strsplit(raw_lines[x + 
                                                                         1], "seqs=")[[1]][2], character(1))
      block.len <- as.numeric(block.len)
      block.starts <- block.starts + 2
      block.stops <- block.starts + block.len - 1
      blocks <- mapply(function(x, y) read.table(text = raw_lines[x:y], 
                                                 stringsAsFactors = FALSE), block.starts, block.stops, 
                       SIMPLIFY = FALSE)
      sites <- lapply(blocks, function(x) x$V4)
      site.names <- lapply(blocks, function(x) x$V1)
      sites <- switch(alph, DNA = lapply(sites, DNAStringSet), 
                      RNA = lapply(sites, RNAStringSet), AA = lapply(sites, 
                                                                     AAStringSet), lapply(sites, BStringSet))
      sites <- mapply(function(x, y) {
        names(x) <- y
        x
      }, sites, site.names, SIMPLIFY = FALSE)
      names(sites) <- mot.names
      if (length(sites) == 1) 
        sites <- sites[[1]]
      if (is.list(sites) && is.list(motif_list)) 
        if (length(sites) != length(motif_list)) 
          sites <- sites[seq_len(length(motif_list))]
      motif_list <- list(motifs = motif_list, sites = sites)
    }
    if (readsites.meta) {
      site.starts <- grep("sites sorted by position p-value", 
                          raw_lines)
      site.stops <- grep("block diagrams$", raw_lines)
      if (length(site.starts) == 0 || length(site.stops) == 
          0) {
        warning("Could not find site P-values in MEME file")
      }
      else {
        site.starts <- site.starts + 2
        site.stops <- site.stops - 1
        site.tables <- mapply(function(x, y) {
          z <- raw_lines[x:y]
          lapply(z, function(x) strsplit(x, "\\s+")[[1]])
        }, site.starts, site.stops, SIMPLIFY = FALSE)
        col.seqname <- 1
        if (all(grepl("Strand", raw_lines[site.starts - 
                                          1]))) {
          col.pos <- 3
          col.pval <- 4
          col.seq <- 6
        }
        else {
          col.pos <- 2
          col.pval <- 3
          col.seq <- 5
        }
        site.tables <- lapply(site.tables, function(x) {
          z1 <- vapply(x, function(x) x[col.seqname], 
                       character(1))
          z2 <- vapply(x, function(x) x[col.pos], character(1))
          z3 <- vapply(x, function(x) x[col.pval], character(1))
          z4 <- vapply(x, function(x) x[col.seq], character(1))
          data.frame(Sequence = z1, Position = as.numeric(z2), 
                     Pvalue = as.numeric(z3), Site = z4, stringsAsFactors = FALSE)
        })
        names(site.tables) <- mot.names
        if (length(site.tables) == 1) 
          site.tables <- site.tables[[1]]
        motif_list <- list(motifs = motif_list$motifs, 
                           sites = motif_list$sites, sites.meta = site.tables)
        summ.start <- grep("Combined block diagrams: non-overlapping sites", 
                           raw_lines)
        if (length(summ.start) == 0) {
          warning("Could not find combined P-values in MEME file")
        }
        else {
          summ.start <- summ.start + 2
          summ.raw <- raw_lines[summ.start:(length(raw_lines) - 
                                              2)]
          need.fix <- grep("\\", summ.raw, fixed = TRUE)
          if (length(need.fix) > 0) {
            need.fix2 <- need.fix + 1
            summ.raw[need.fix] <- vapply(summ.raw[need.fix], 
                                         function(x) strsplit(x, "\\", fixed = TRUE)[[1]], 
                                         character(1))
            summ.raw[need.fix2] <- gsub("\\s+", "", 
                                        summ.raw[need.fix2])
            summ.raw[need.fix] <- mapply(function(x, 
                                                  y) paste0(x, y), summ.raw[need.fix], summ.raw[need.fix2])
            summ.raw <- summ.raw[-need.fix2]
          }
          summ.tab <- read.table(text = summ.raw, stringsAsFactors = FALSE)
          colnames(summ.tab) <- c("Sequence", "Combined.Pvalue", 
                                  "Diagram")
          motif_list <- list(motifs = motif_list$motifs, 
                             sites = motif_list$sites, sites.meta = motif_list$sites.meta, 
                             sites.meta.combined = summ.tab)
        }
      }
    }
  }
  else if (readsites.meta) {
    warning("'readsites.meta' is not valid if 'readsites = FALSE'")
  }
  motif_list
}

environment(read_meme2) <- asNamespace('universalmotif')
assignInNamespace("read_meme", read_meme2, ns = "universalmotif")


#Now reading Streme/Homer motif file 
motif_file <- read_meme2(arguments$input)

#Extract sequences
if (length(motif_file) == 0) {
  motifs_sequences <- data.frame()
 } else if (length(motif_file) == 1) {
  motifs_sequences <- data.frame(name = motif_file@name, 
                                 sequence = get_matches(motif_file@motif, motif_score(motif_file@motif)[2]) )
 } else {
   motifs_sequences <- sapply(motif_file, function(x) cbind(x@name, get_matches(x@motif, motif_score(x@motif)[2])))
 }

#Sometimes it return as a list, sometimes as an array, so have to check
if (length(motif_file) == 0 | length(motif_file) == 1) {
  #Do nothing
} else if (is.list(motifs_sequences)) {
  motifs_sequences <- Reduce(rbind, motifs_sequences)   
  motifs_sequences <- data.frame(name = motifs_sequences[,1],
                                 sequence = motifs_sequences[,2])
  } else { 
   motifs_sequences <- data.frame(name = motifs_sequences[1,],
                               sequence = motifs_sequences[2,]) 
   }

#Remove similar motifs in highstab between lowstab
if  (grepl(x = arguments$input, pattern = "highstab", fixed = TRUE) ) {
  tomtom <- try(read.table(arguments$tomtom, fill = T, header = T))
  #Check if tomtom file is empty
  if (class(tomtom) == "try-error" | any(is.na(tomtom))) {
    to_remove <- vector()
    filtered_sequences <- motifs_sequences
  } else {
    tomtom <- (tomtom[1:(nrow(tomtom)-4),])
    to_remove <- unique(tomtom$Query_ID)
    tomtom_remove <- motifs_sequences[motifs_sequences[, 1] %in% to_remove,]
    filtered_sequences <- motifs_sequences[!(motifs_sequences[, 1] %in% to_remove),]
  }
} else {
  filtered_sequences <- motifs_sequences
}

#Read fire input
fire <- try(read.delim(arguments$fire, header = F))

if (class(fire) != "try-error") {
  fire_filter <- try(read.delim(arguments$filter_fire, header = F))
  if (class(fire_filter) != "try-error") {
    #Remove Highstab in lowstab for fire
    fire <- fire[!(fire$V2 %in% fire_filter$V2),]
  } else if (class(fire_filter) == "try-error" ){
    fire_filter <- data.frame()  
    }
  #Add to homer and streme  
  colnames(fire) <- c("name","sequence")
  filtered_sequences <- rbind(filtered_sequences, fire)
} else {
  fire <- data.frame()
  fire_filter <- data.frame()  
}

if (nrow(filtered_sequences) == 0) {
  write.table("The pipeline didn't output any motifs.", 
              paste0(arguments$output, ".log") ,
              sep = "\n",
              col.names = F,
              row.names = F,
              quote = F)
  file.create(arguments$output)
  
} else {

#Checking for miRNA seed targets
seeds_motifs <- read_meme2(arguments$miRNA_seeds)  
seeds_sequences <- sapply(seeds_motifs, function(x) cbind(x@name, get_matches(x@motif, motif_score(x@motif)[2]))) 
seeds_sequences <- data.frame(name = seeds_sequences[1,],
                              sequence = seeds_sequences[2,])


#######Filter for relevant miRNAs for A549 and Hepg2 if necessary######
relevant <- try(read.table(arguments$relevant)) 
if (class(relevant) != "try-error") {
 #intersect(relevant[,1], seeds_sequences[,1]) %>% length()
#missing 17 entries
seeds_sequences <- try(seeds_sequences %>% 
  dplyr::filter(seeds_sequences[,1] %in% relevant[,1])) 
}

#######################################################

match_sirna <- lapply(seeds_sequences[,2], function(x) filtered_sequences[str_detect(filtered_sequences[,2], x),]  )
names(match_sirna) <- paste(seeds_sequences[,1] , seeds_sequences[,2], sep = ":")
match_sirna <-  match_sirna[lapply(match_sirna,nrow) >0]
match_sirna_table <- do.call (rbind , match_sirna)
#matched_sequences <- sapply(match_sirna, function(x) x$sequence) %>% unique()

if  (grepl(x = arguments$input, pattern = "highstab", fixed = TRUE) ) {
  filtered_sequences <- filtered_sequences[!(filtered_sequences[, 2] %in% match_sirna_table$sequence),]
}

#Check for most common polyA signal(s)
poly <- length(str_which(filtered_sequences[,2], "AUAAA"))

if (poly != 0 && grepl(x = arguments$input, pattern = "highstab", fixed = TRUE)) {
  polyA <- str_subset(filtered_sequences[,2], "AUAAA")
  filtered_sequences <- filtered_sequences[!(filtered_sequences[, 2] %in% polyA),]
}



#Log file printing
if (grepl(x = arguments$input, pattern = "highstab", fixed = TRUE) ) {
  
  log <- c(paste0("Number of motifs from Streme and Homer: ", length(motif_file)),
           paste0("Total number of motifs sequences from Streme and Homer: ", nrow(motifs_sequences)),
           paste0("Number of motifs sequences from Fire: (if 0, or empty or lowstab filtered) ", nrow(fire)),
           paste0("Number of motifs similar to lowstab motifs for Streme and Homer : ", length(to_remove)),
           paste0("Number of motifs similar to lowstab motifs for Fire : ", nrow(fire_filter)),
           paste0("Number of sequences matching miRNA seed targets: ", length(match_sirna_table$sequence)),
           paste0("Number of sequences with polyA signal AUAAA: ", poly),
           paste0("Final number of sequences: ", nrow(filtered_sequences)),
           "See .matching.mirna.seeds file to find motifs targeted by miRNAs")
  write.table(filtered_sequences ,file = arguments$output, sep = "\t", col.names = T, quote = F,
              row.names = F)
  write.table(log, 
              paste0(arguments$output, ".log") ,
              sep = "\n",
              col.names = F,
              row.names = F,
              quote = F)
} else {
  log <- c(paste0("Number of motifs from Streme and Homer: ", length(motif_file)),
           paste0("Total number of motifs sequences from Streme and Homer: ", nrow(motifs_sequences)),
           paste0("Number of motifs sequences from Fire: ", nrow(fire)),
           paste0("Total number of motifs sequences: ", nrow(filtered_sequences)),
           paste0("Number of sequences matching miRNA seed targets: ", length(match_sirna_table$sequence)),
           paste0("Number of sequences with polyA signal AUAAA: ", poly),
           "See .matching.mirna.seeds file to find motifs targeted by miRNAs",
           "Lowstab motifs have not been filetred")
  write.table(filtered_sequences,file = arguments$output, sep = "\t", col.names = T, quote = F,
              row.names = F)
  write.table(log, 
              paste0(arguments$output, ".log") ,
              sep = "\n",
              col.names = F,
              row.names = F,
              quote = F)
}

write.table(match_sirna_table, 
            str_replace(arguments$output, ".list", ".matching.mirna.seeds"),
            sep = "\t",
            row.names = T,
            quote = F)

}