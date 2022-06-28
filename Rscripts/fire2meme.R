library("optparse")
library(stringr)

#Options parser
option_list = list(
  make_option(c("-i", "--fire-input"),
              type="character",
              dest = "fire_path",
              help="full path directory of fire extracted motif file)"),
  make_option(c("-b", "--background-input"),
              type="character",
              dest = "background_path",
              help="full path directory of background.fasta.bg")
  
)

arguments <- parse_args(OptionParser(option_list = option_list))

# arguments <- data.frame(fire_path = "fire_residuals.txt_FIRE/RNA/fire_residuals_highstab.signif.motifs.temp",
#                         background_path = "fire.fasta.bg"
#                         )
# inFile = arguments$fire_path
# fire.bg = arguments$background_path

#converts a Homer .motif PWM/PFM into MEME format
motif2meme <- function(inFile, fire.bg) {
    library(tools)
    #establishing the number of distinct motifs in the file and parsing it accordingly
    stopifnot(is.character(inFile))
    outFile <- str_remove(inFile, ".temp")
    outFile <- paste(outFile,"meme",sep=".")
    thisFile <- file(outFile)
    fileName <- file_path_sans_ext(inFile)
    #reading the input file
    motif.file <- scan(file=inFile,character(0), sep="\n",quote=NULL)
    if (length(motif.file) == 0) {
      file.create(outFile)
    } else {
    motif.file <- str_replace_all(motif.file, "T", "U")
    motif.index <- grep(pattern="^MEME version 5.",motif.file)
    n.motifs <- length(motif.index)
    total.len <-  length(motif.file)
    bg <- scan(file=fire.bg,
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
        (motif.index[i]+7) -> index.start
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
        motif.file[motif.index[i]+5] -> header.string
        strsplit(header.string, split = " ") -> header_split
        motif.file[motif.index[i]+6] -> prob.string
        header_split[[1]][2] -> name.string
        motif_name <- paste("FIRE",i,sep="-")
        cat("MOTIF",name.string,motif_name,"\n",file=thisFile, append=TRUE)
        cat(prob.string,"\n",file=thisFile, append=TRUE)
        write.table(motif_array,file=thisFile,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t")
        cat("\n",file=thisFile, append=TRUE)
    }
    #Same for last motif
   (motif.index[n.motifs]+7) -> index.start
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
   motif.file[motif.index[n.motifs]+5] -> header.string
   strsplit(header.string, split = " ") -> header_split
   motif.file[motif.index[n.motifs]+6] -> prob.string
   header_split[[1]][2] -> name.string
   motif_name <- paste("FIRE",n.motifs,sep="-")
   cat("MOTIF",name.string,motif_name,"\n",file=thisFile, append=TRUE)
   cat(prob.string,"\n",file=thisFile, append=TRUE)
   write.table(motif_array,file=thisFile,append=TRUE,col.names=FALSE,row.names=FALSE,sep="\t")
   cat("\n",file=thisFile, append=TRUE)
   sink()
   close(thisFile) 
   print("matrix has been converted to MEME")

   } #end else
}
motif2meme(inFile = arguments$fire_path,
           fire.bg = arguments$background_path)
