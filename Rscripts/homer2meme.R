library("optparse")
library(stringr)

###
#This script has been adapted from https://gist.github.com/rtraborn/e395776b965398c54c4d
###


#Options parser
option_list = list(
  make_option(c("-i", "--homer-input"),
              type="character",
              dest = "homer_path",
              help="full path directory of homer motif analysis output, 
                    homerMotifs.all.motifs"),
  make_option(c("-b", "--background-input"),
              type="character",
              dest = "background_path",
              help="full path directory of background.fasta.bg"),
  make_option(c("-t", "--FDR-threshold"),
              type="numeric",
              dest = "fdr_threshold",
              help="FDR threshold value, keep only motifs below this threshold")
  
)

arguments <- parse_args(OptionParser(option_list = option_list))

# arguments <- data.frame(homer_path = "/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/a549_motif_analysis/one_transcript/halflife_lowstab_homer.dir/homerMotifs.all.motifs",
#                         background_path = "/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/hepg2_motif_analysis/one_transcript/background.fasta.bg",
#                         fdr_threshold = 0.1
#                         )
inFile = arguments$homer_path
homer.bg = arguments$background_path

#converts a Homer .motif PWM/PFM into MEME format
motif2meme <- function(inFile, homer.bg) {
    library(tools)
    #establishing the number of distinct motifs in the file and parsing it accordingly
    stopifnot(is.character(inFile))
    outFile <- paste(inFile,"meme",sep=".")
    thisFile <- file(outFile)
    fileName <- file_path_sans_ext(inFile)
    #reading the input file
    motif.file <- scan(file=inFile,character(0), sep="\n",quote=NULL)
    motif.file <- str_replace_all(motif.file, "T", "U")
    motif.index <- grep(pattern="^>",motif.file)
    n.motifs <- length(motif.index)
    total.len <-  length(motif.file)
    bg <- scan(file=homer.bg,
         character(0), sep="\n",quote=NULL)
    bg <- bg[4:7]
    #print(n.motifs)
    sink(thisFile,append=TRUE)
    cat("MEME version 5\n\n",file=thisFile,append=TRUE)
    cat("ALPHABET= ACGU\n\n",file=thisFile,append=TRUE)
    cat("strands: + \n\n")
    cat("Background letter frequencies\n", file=thisFile,append=TRUE)
    cat(bg, "\n\n",file=thisFile,append=TRUE)
    for (i in 1:(n.motifs-1)) {
        #print(i)
        (motif.index[i]+1) -> index.start
        ((motif.index[i+1])-1) -> index.end
        motif.file[index.start:index.end] -> this_motif
        strsplit(this_motif,split="\t") -> motif_split
        #print(head(motif_split))
        length(this_motif) -> motif_row
        array(NA,c(motif_row,4)) -> motif_array
        for (n in 1:motif_row) {
        as.numeric(motif_split[[n]]) -> motif_array[n,]
        }
        motif_array <- as.data.frame(motif_array)
        motif.file[motif.index[i]] -> header.string
        strsplit(header.string,split="[\t]") -> prob.string
        strsplit(prob.string[[1]][6],"[,]") -> prob.string3
        prob.string3[[1]][4] -> prob.fdr
        str_remove_all(prob.fdr, "FDR:") -> prob.fdr
        as.numeric(prob.fdr) -> prob.fdr
        if (prob.fdr > 0.1) {next}
        prob.string3[[1]][3] -> this.p.val
        str_replace(this.p.val,":", "= ") -> this.p.val
        prob.string[[1]][2] -> name.string
        motif_name <- paste("HOMER",i,sep="-")
        cat("MOTIF",name.string,motif_name,"\n",file=thisFile, append=TRUE)
        cat("letter-probability matrix: ",file=thisFile, append=TRUE)
        cat("alength= 4 w=", motif_row, file=thisFile, append=TRUE)
        cat(" nsites= 20 ",file=thisFile, append=TRUE)
        cat(this.p.val,"\n",file=thisFile, append=TRUE)
        write.table(motif_array,file=thisFile,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t")
        cat("\n",file=thisFile, append=TRUE)
    }
    #Same for last motif
   (motif.index[n.motifs]+1) -> index.start
   total.len -> index.end
   motif.file[index.start:index.end] -> this_motif
   strsplit(this_motif,split="\t") -> motif_split
   #print(head(motif_split))
   length(this_motif) -> motif_row
   array(NA,c(motif_row,4)) -> motif_array
   for (n in 1:motif_row) {
     as.numeric(motif_split[[n]]) -> motif_array[n,]
   }
   motif_array <- as.data.frame(motif_array)
   motif.file[motif.index[n.motifs]] -> header.string
   strsplit(header.string,split="[\t]") -> prob.string
   strsplit(prob.string[[1]][6],"[,]") -> prob.string3
   prob.string3[[1]][4] -> prob.fdr
   str_remove_all(prob.fdr, "FDR:") -> prob.fdr
   as.numeric(prob.fdr) -> prob.fdr
   if (prob.fdr < 0.1) {
   prob.string3[[1]][3] -> this.p.val
   str_replace(this.p.val,":", "= ") -> this.p.val
   prob.string[[1]][2] -> name.string
   motif_name <- paste("HOMER",n.motifs,sep="-")
   cat("MOTIF",name.string,motif_name,"\n",file=thisFile, append=TRUE)
   cat("letter-probability matrix: ",file=thisFile, append=TRUE)
   cat("alength= 4 w=", motif_row, file=thisFile, append=TRUE)
   cat(" nsites= 20 ",file=thisFile, append=TRUE)
   cat(this.p.val,"\n",file=thisFile, append=TRUE)
   write.table(motif_array,file=thisFile,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t")
   cat("\n",file=thisFile, append=TRUE)
   }
   sink()
   close(thisFile) 
   print("matrix has been converted to MEME")
}
 
motif2meme(inFile = arguments$homer_path,
           homer.bg = arguments$background_path)
