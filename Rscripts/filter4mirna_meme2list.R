library("universalmotif")
library("optparse")
library("stringr")

#Options parser
option_list = list(
  make_option(c("-i", "--input-meme"),
              type="character",
              dest = "input",
              help="highstab input file in meme format"),
  make_option(c("-m", "--mirna-seeds"),
              type="character",
              dest = "miRNA_seeds",
              help="miRNA seeds target sequences"),
  make_option(c("-l", "--lax-score"),
              action="store_true",
              default=FALSE,
              dest = "score",
              help="If option added, icscore is used to select sequences from the 
              motif PWM, otherwise only 1 sequence is selected from it.
              See universalmotif vignette for details.")
)

arguments <- parse_args(OptionParser(option_list = option_list))

# setwd("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/hepg2_slam/motif_analysis/one_transcript/")
# arguments <- data.frame(score = FALSE,
#                         input = "final_motifs/highstab_final_motifs.meme",
#                         miRNA_seeds = "/mnt/sharc/shared/sudlab1/General/mirror/meme_motif_db/motif_databases/MIRBASE/22/Homo_sapiens_hsa.seeds.sequences")

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


#Now reading file 
motif_file <- read_meme2(arguments$input)

seeds_sequences <- (read.table(arguments$miRNA_seeds, header = F))
seeds_list <- list()
for(i in 1:nrow(seeds_sequences)) {        
  seeds_list[[i]] <- seeds_sequences[ i, 1]
}

motif_sequences <- character()
for (i in 1:length(motif_file)) {
  if (arguments$score) {
    #icscore threshold used
    motifs_i <- get_matches(motif_file[[i]]@motif, motif_file[[i]]@icscore)
    motif_sequences <- cbind(motif_sequences, motifs_i)
  } else {
    #Only 1
    max_score <- motif_score(motif_file[[i]]@motif)[2]
    motifs_i <- get_matches(motif_file[[i]]@motif, max_score) 
    motif_sequences <- cbind(motif_sequences, motifs_i)
  }
}

matches <- lapply(seeds_list, FUN = function(x) motif_sequences[str_which(motif_sequences, x)])

matches <- matches[lapply(matches,length) >0]

matches <- Reduce(c, matches) %>% Reduce(c, .) %>% unique() 


filtered_sequences <- motif_sequences[!(motif_sequences %in% matches)]
poly <- length(str_which(filtered_sequences, "AAUAAA"))
filtered_sequences <- filtered_sequences[-(str_which(filtered_sequences, "AAUAAA"))]

log <- c(paste0("Total number of motifs: ", length(motif_file)), 
         paste0("Total number of motifs sequences: ", length(motif_sequences)),
         paste0("Number of sequences matching miRNA seed targets: ", length(matches)),
         paste0("Number of sequences with polyA signal AUAAA: ", length(poly)),
         "Matching motifs to miRNA seeds:")


out <- str_replace(arguments$input, ".meme", ".list")

write.table(filtered_sequences ,file = out, sep = "\n", col.names = F, quote = F,
            row.names = F)

write.table(c(log,matches), 
            paste0(out, ".log") ,
            sep = "\n",
            col.names = F,
            row.names = F,
            quote = F)



