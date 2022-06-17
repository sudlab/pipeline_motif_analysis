
library(glmnet)
library(tidyverse)
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(caret)

#This R script Generate bed files for motif analysis, i will determine sampling
# size for testing and training, it is determine by doing the lasso_CodonFreq
# R markdown.


###INPUTS###
############

#Percent of testing for best model outcomes after running lasso_Codonfreq
i = 0.85

#Output directory
outdir <- "/mnt/sharc/shared/sudlab1/General/projects/SLAMseq_CHO_Mikayla/cho_slam/slam_picr/lasso/"
#Fasta file of CDS sequence
cds_fasta <- "/mnt/sharc/shared/sudlab1/General/mirror/ensembl/criGri_PICR/Cricetulus_griseus_picr.CriGri-PICR.cds.all.rename.fa.gz"
#Table of transcript to gene ids
transcript_to_gene <- "/mnt/sharc/shared/sudlab1/General/mirror/ensembl/criGri_PICR/Cricetulus_griseus_picr.CriGri-PICR.105.ena.tsv"
#Table of half-life for transcripts from slamdunk
halflife_table <- "/mnt/sharc/shared/sudlab1/General/projects/SLAMseq_CHO_Mikayla/cho_slam/slam_picr/halflife_2rep_aboveCPM/halflife_filtered2.tsv"
#Seed number to reproduce results
seed_num <- 44

###CDSs###
##########
cds_Cfreq <- readDNAStringSet(cds_fasta,
                         format = "fasta")

cds_Cfreq_cf <- oligonucleotideFrequency(cds_Cfreq, width = 3, step = 3)
cds_Cfreq_id <- cds_Cfreq@ranges@NAMES
row.names(cds_Cfreq_cf) <- cds_Cfreq_id
cds_lengths = data.frame(transcript_id = cds_Cfreq_id, cds_length = rowSums(cds_Cfreq_cf))

tx2gene <- read.table(transcript_to_gene, header=TRUE,
                      fill = T) %>%
  select(c(gene_stable_id, transcript_stable_id))

###Stab data for both cell lines###
###################################

stability <- read.delim(halflife_table, header = T) %>%
  select(chr, start, end, transcript_id, length, strand, halflife) %>%
  inner_join(cds_lengths) %>%
  inner_join(tx2gene %>% distinct(), by = c("transcript_id" = "transcript_stable_id")) %>%
  group_by(gene_stable_id) %>%
  slice_max(order_by=cds_length, n=1, with_ties = FALSE) %>% ungroup()


#prep codon freq
selected_cf <- cds_Cfreq_cf[stability$transcript_id,]

#Normalize
normed_cf <- sweep(selected_cf, 1, rowSums(selected_cf), "/")
# Add the CDS length of the condon frequencies
normed_cf <- cbind(normed_cf, rowSums(selected_cf))
colnames(normed_cf)[ncol(normed_cf)] <- "cds_length"

final_data <- merge(stability %>%
                      select(transcript_id, halflife), normed_cf, by.x = "transcript_id", by.y = 0)
#Normalize half-lives
final_data$halflife <- log(final_data$halflife)

#Scale
final_data[,3:ncol(final_data)] <- scale(final_data[,3:ncol(final_data)])

###LASSO###
###########

  #Splitting data
  set.seed(seed_num)
  index <- createDataPartition(final_data$halflife, p= i,
                               list = F, times = 1)
  train_data<-   final_data[index,]
  test_data <- final_data[-index,]
  
  #Checking that training and testing don't share same transcripts
  shared_transcript <- test_data[(test_data$transcript_id %in% train_data$transcript_id),]

  #Remove them from the testing, go to training
  test_data <- test_data[!(test_data$transcript_id %in% train_data$transcript_id),]
  #Check again to be sure
  if ((nrow(test_data[(test_data$transcript_id %in% train_data$transcript_id),]) == 0) == F) {
    print("Training and testing set contain same genes")
    break}

  #Add removed rows to train
  train_data <- rbind(train_data, shared_transcript)
  #Check no duplicate, if FALSE something is not ok
  if ((nrow(distinct(train_data)) == nrow(train_data)) == F) {
    print("Duplicates in traning set")
    break}

  #Check split %
  if ((nrow(train_data) + nrow(test_data) == nrow(final_data)) == F) {
    print("Traning and testing data aren't correct")
    break}
  how_many <- round((nrow(test_data)/(nrow(final_data)) *100), 0)


  #Creating 5 grouping folds for k-fold cv
  groud_folds <- groupKFold(train_data$transcript_id, k=5)
  ctrl_kfold <- trainControl(method = 'cv', index = groud_folds,
                             savePredictions = "final")

  # Create vector of lambda
  lambda_vector <- 10^seq(5, -5, length=500)

  #Formula for Lasso, (just so transcript_ids aren't in it)
  my_x_vars <- paste0(colnames(train_data)[3:ncol(train_data)], collapse = "+")
  formula <- as.formula(paste('halflife ~ ', my_x_vars))

  #Traning
  set.seed(seed_num)
  model1 <- train(formula,
                  data = train_data,
                  method= "glmnet",
                  tuneGrid=expand.grid(alpha=1, lambda = lambda_vector),
                  trControl=ctrl_kfold)

  training_perf_model1 <- getTrainPerf(model1)[,1:2]

  #Predict
  pred_model1 <- predict(model1,
                         newdata = select(test_data, - c(halflife, transcript_id)))
pred_model1 <- cbind(test_data[,1:2], pred_model1)


# Get residuals and generate bed files -----------------------------------------------------------

#For training and testing set of model
testing_residuals <- pred_model1
testing_residuals$Residuals <- testing_residuals$halflife - testing_residuals$pred_model1
colnames(testing_residuals)[3] <- "Predictions"
training_residuals <- predict(model1,
                              newdata = select(train_data, - c(halflife, transcript_id)))
training_residuals <- cbind(train_data[,1:2], training_residuals)
training_residuals$Residuals <- training_residuals$halflife - training_residuals$training_residuals
colnames(training_residuals)[3] <- "Predictions"

final_residual_table <- rbind(testing_residuals, training_residuals)
colnames(final_residual_table)[2] <- "LogHalflife"
final_residual_table$halflife <- exp(final_residual_table$LogHalflife)

final_residual_table<- final_residual_table %>%
  inner_join(stability %>% select(chr, start, end, transcript_id, strand, length))

final_residual_table %>%
  arrange(halflife) %>%
  write.table(paste0(outdir,"stab_residual_one_transcript.tsv"),quote = F, sep = "\t", row.names = F)

#Bed files
number = round(nrow(final_residual_table)*0.20)

#Half-life ranking
final_residual_table %>% arrange(halflife) %>% head(n = number) %>%
  bind_cols(rep(".", n = nrow(final_residual_table))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/one_transcript/","halflife_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
final_residual_table %>% arrange(halflife) %>% tail(n = number) %>%
  bind_cols(rep(".", n = nrow(final_residual_table))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/one_transcript/","halflife_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Regressed out ranking
final_residual_table %>% arrange(Residuals) %>% head(n = number) %>%
  bind_cols(rep(".", n = nrow(final_residual_table))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/one_transcript/","residual_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
final_residual_table %>% arrange(Residuals) %>% tail(n = number) %>%
  bind_cols(rep(".", n = nrow(final_residual_table))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/one_transcript/","residual_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Background
final_residual_table %>% arrange(halflife) %>%
  bind_cols(rep(".", n = nrow(final_residual_table))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/one_transcript/","background.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Fire INPUTS
final_residual_table %>% select(transcript_id, halflife) %>%
  write.table(paste0(outdir,"fire_halflife.txt"), row.names = F, quote = F, sep = "\t")
final_residual_table %>% select(transcript_id, Residuals) %>%
  write.table(paste0(outdir,"fire_residuals.txt"), row.names = F, quote = F, sep = "\t")

final_residual_table %>%
  bind_cols(rep(".", n = nrow(final_residual_table))) %>%
  filter(length < 10000 & length > 6) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  arrange(transcript_id, start) %>%
  write.table(paste0(outdir,"fire.bed"), row.names = F, quote = F, sep = "\t", col.names = F)

# Get Residuals  for all transcripts -----------------------------------------------------------

#Let's predict for all
#Mean transcripts that have several half-life per transcripts (the ones that
#have several entries in the bed, causeexons in 3'UTR, this is like 50 of them,
#so nothing compared to the 16500 I have)
stability_all <-  read.delim(halflife_table, header = T) %>%
  select(chr, start, end, transcript_id, length, strand, halflife) %>%
  group_by(transcript_id) %>%
  transform(halflife, mean(halflife) )


#prep codon freq
selected_cf_all <- cds_Cfreq_cf[stability_all$transcript_id,]

#Normalize
normed_cf_all <- sweep(selected_cf_all, 1, rowSums(selected_cf_all), "/")
# Add the CDS length of the condon frequencies
normed_cf_all <- cbind(normed_cf_all, rowSums(selected_cf_all))
colnames(normed_cf_all)[ncol(normed_cf_all)] <- "cds_length"

final_stab_all <- merge(stability_all %>%
                          select(transcript_id, halflife), normed_cf_all, by.x = "transcript_id", by.y = 0)

#Normalize half-lives
final_stab_all$halflife <- log(final_stab_all$halflife)

#Scale
final_stab_all[,3:ncol(final_stab_all)] <- scale(final_stab_all[,3:ncol(final_stab_all)])

pred_bestModel_all <- predict(model1,
                              newdata = select(final_stab_all, - c(halflife, transcript_id)))


pred_bestModel_all <- cbind(final_stab_all[,1:2],
                            pred_bestModel_all)
head(pred_bestModel_all)

# Residuals
all_residuals <- pred_bestModel_all
all_residuals$Residuals <- all_residuals$halflife - all_residuals$pred_bestModel_all
colnames(all_residuals)[3] <- "Predictions"
colnames(all_residuals)[2] <- "LogHalflife"
all_residuals$halflife <- exp(all_residuals$LogHalflife)

comparison_testing_all_only_one_transcript <- merge(all_residuals, final_residual_table, by = 1)
cor.test(comparison_testing_all_only_one_transcript$Predictions.x, comparison_testing_all_only_one_transcript$Predictions.y)
#0.9997261

sd(final_residual_table$Residuals) / mean((final_residual_table$Residuals) )
sd(all_residuals$Residuals) / mean(all_residuals$Residuals)

all_residuals <- all_residuals %>% inner_join(stability_all %>% select(chr, start, end, transcript_id, strand, length))
all_residuals %>% arrange(halflife) %>%
  write.table(paste0(outdir,"/all_transcripts/stab_residuals_all.tsv"),quote = F, sep = "\t", row.names = F)


# Generate bed files ------------------------------------------------------

number = round(nrow(all_residuals)*0.20)

#Half-life ranking
all_residuals %>% arrange(halflife) %>% head(n = number) %>%
  bind_cols(rep(".", n = nrow(all_residuals))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/all_transcripts/","halflife_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
all_residuals %>% arrange(halflife) %>% tail(n = number) %>%
  bind_cols(rep(".", n = nrow(all_residuals))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/all_transcripts/","halflife_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")


#Regressed out ranking
all_residuals %>% arrange(Residuals) %>% head(n = number) %>%
  bind_cols(rep(".", n = nrow(all_residuals))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/all_transcripts/","residual_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
all_residuals %>% arrange(Residuals) %>% tail(n = number) %>%
  bind_cols(rep(".", n = nrow(all_residuals))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/all_transcripts/","residual_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Background
all_residuals %>% arrange(halflife) %>%
  bind_cols(rep(".", n = nrow(all_residuals))) %>%
  select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(outdir, "/all_transcripts/","background.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Fire INPUTS
all_residuals %>% select(transcript_id, halflife) %>%
    write.table(paste0(outdir, "/all_transcripts/","fire_halflife.txt"), row.names = F, quote = F, sep = "\t")
all_residuals %>% select(transcript_id, Residuals) %>%
    write.table(paste0(outdir, "/all_transcripts/","fire_residuals.txt"), row.names = F, quote = F, sep = "\t")

all_residuals %>%
    bind_cols(rep(".", n = nrow(all_residuals))) %>%
    filter(length < 10000 & length > 6) %>%
    select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
    arrange(transcript_id, start) %>%
    write.table(paste0(outdir, "/all_transcripts/","fire.bed"), row.names = F, quote = F, sep = "\t", col.names = F)
