
library(glmnet)
library(tidyverse)
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(caret)
library(reshape)


###INPUTS###
############
i = 0.7 #gives sampling best model 

outdir <- "/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/a549_slam/motif_analysis/lasso/"
cds_fasta <- "/mnt/sharc/shared/sudlab1/General/projects/example_project/cds_hg38_noalt/Homo_sapiens.GRCh38.cds.all.fa.gz"
transcript_to_gene <- "/mnt/sharc/shared/sudlab1/General/mirror/ensembl/hg38_ensembl93/Homo_sapiens.GRCh38.93.IDs.tsv"
seed_num <- 44
###CDSs###
##########
hg38 <- readDNAStringSet(cds_fasta,
                         format = "fasta")

hg38_cf <- oligonucleotideFrequency(hg38, width = 3, step = 3)
hg38_id <- hg38@ranges@NAMES
row.names(hg38_cf) <- hg38_id
cds_lengths = data.frame(transcript_id = hg38_id, cds_length = rowSums(hg38_cf))

tx2gene <- read.table(transcript_to_gene, header=TRUE,
                      fill = T) %>%
  dplyr::select(c(gene_stable_id, transcript_stable_id))

###Stab data for both cell lines###
###################################

stability_a549 <- read.delim("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/a549_slam/slam/halflife_2rep_aboveCPM//halflife_filtered.tsv",
                             header = T) %>%
  dplyr::select(chr, start, end, transcript_id, length, strand, halflife) %>%
  inner_join(cds_lengths) %>%
  inner_join(tx2gene %>% distinct(), 
             by = c("transcript_id" = "transcript_stable_id")) %>% 
  group_by(gene_stable_id) %>%
  slice_max(order_by=cds_length, n=1, with_ties = FALSE)  %>% ungroup()

stability_hep <- read.delim("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/hepg2_slam/slam_downsamplingbg//halflife_2rep_aboveCPM//halflife_filtered.tsv",
                            header = T) %>%
  dplyr::select(chr, start, end, transcript_id, length, strand, halflife) %>%
  inner_join(cds_lengths) %>%
  inner_join(tx2gene %>% distinct(), 
             by = c("transcript_id" = "transcript_stable_id")) %>% 
  group_by(gene_stable_id) %>%
  slice_max(order_by=cds_length, n=1, with_ties = FALSE)  %>% ungroup()

###Creating merged table### 
###########################

stabilities <- bind_rows(a549 = stability_a549, hep = stability_hep, .id = "cell_line") %>% 
  dplyr::select(cell_line, transcript_id, halflife) 

#define one-hot encoding function
dummy <- dummyVars("~ cell_line", data=stabilities)

#perform one-hot encoding on data frame
final_stab <- data.frame(transcript_id = stabilities$transcript_id,
                         halflife = stabilities$halflife,
                         predict(dummy, newdata=stabilities))

#prep codon freq
selected_cf <- hg38_cf[final_stab$transcript_id,]

#Normalize
normed_cf <- sweep(selected_cf, 1, rowSums(selected_cf), "/")
# Add the CDS length of the condon frequencies
normed_cf <- cbind(normed_cf, rowSums(selected_cf))
colnames(normed_cf)[ncol(normed_cf)] <- "cds_length"

final_data <- merge(final_stab, normed_cf, by.x = "transcript_id", by.y = 0)


#Normalize half-lives
final_data$halflife <- log(final_data$halflife)


###Interaction terms###
#######################

inter_a549 <- final_data[,"cell_linea549"] * final_data[,5:ncol(final_data)]
colnames(inter_a549) <- paste0(colnames(inter_a549), "_ia549")
inter_hep <- final_data[,"cell_linehep"] * final_data[,5:ncol(final_data)]
colnames(inter_hep) <- paste0(colnames(inter_hep), "_ihep")

final_data <- cbind(final_data, inter_a549, inter_hep)

#Scale
final_data[,5:ncol(final_data)] <- scale(final_data[,5:ncol(final_data)])


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
  
  round((nrow(test_data)/(nrow(final_data)) *100), 0)
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
                         newdata = dplyr::select(test_data, - c(halflife, transcript_id)))
  pred_model1 <- cbind(test_data[,1:4], pred_model1)
  

###Residuals and bed files###
#############################

# Get Residuals  -----------------------------------------------------------

testing_residuals <- pred_model1
testing_residuals$Residuals <- testing_residuals$halflife - testing_residuals$pred_model1 
colnames(testing_residuals)[5] <- "Predictions"


cor.test(x = testing_residuals$halflife, y = testing_residuals$Residuals)

testing_residuals %>% dplyr::select(transcript_id, cell_linea549, halflife, Residuals) %>%
  #unite(c(transcript_id, cell_linea549),col = "trans_cell" ,sep = "_") %>%
  ggplot(aes(halflife, Residuals,col = as.factor(cell_linea549)))  +
  geom_point(size = 0.5) + theme_light()  + 
  stat_smooth(method="lm", color = "black") +
  xlab("Predicted half-life") + ylab("Observed half-life")
ggsave(paste0(outdir,"corr_grap", ".jpeg"))

training_residuals <- predict(model1, 
                              newdata = dplyr::select(train_data, - c(halflife, transcript_id)))
training_residuals <- cbind(train_data[,1:4], training_residuals)
training_residuals$Residuals <- training_residuals$halflife - training_residuals$training_residuals 
colnames(training_residuals)[5] <- "Predictions"

final_residual_table <- rbind(testing_residuals, training_residuals)
colnames(final_residual_table)[2] <- "LogHalflife"
final_residual_table$halflife <- exp(final_residual_table$LogHalflife)

residual_table_a549 <- final_residual_table %>%
  filter(cell_linea549 == 1) %>%
  inner_join(y = (stability_a549 %>% 
                    dplyr::select(chr, start, end, transcript_id, length, strand, length)), 
             by = "transcript_id")

residual_table_hep <- final_residual_table %>%
  filter(cell_linehep == 1) %>%
  inner_join(y = (stability_hep %>% 
                    dplyr::select(chr, start, end, transcript_id, length, strand, length)), 
             by = "transcript_id")

sd(final_residual_table$Residuals) / mean((final_residual_table$Residuals) )

residual_table_a549 %>% arrange(halflife) %>% 
  dplyr::select(- c(cell_linea549, cell_linehep)) %>%
  write.table(paste0(outdir,"stab_residual_a549.tsv"),quote = F, sep = "\t", row.names = F)

residual_table_hep %>% arrange(halflife) %>%
  dplyr::select(- c(cell_linea549, cell_linehep)) %>%
  write.table(paste0(outdir,"stab_residual_hep.tsv"),quote = F, sep = "\t", row.names = F)

residual_table_hep <- read_delim(paste0(outdir,"stab_residual_hep.tsv"))
residual_table_a549 <- read_delim(paste0(outdir,"stab_residual_a549.tsv"))

bind_rows(list(a549 = residual_table_a549, hepg2 = residual_table_hep),.id = "cell_line") %>% 
  select(transcript_id, LogHalflife, Residuals, cell_line) %>%
ggplot(aes(LogHalflife, Residuals,col = as.factor(cell_line)))  +
  geom_point(size = 0.5, alpha = 0.5, shape = 3) 
+ theme_light()  + 
  stat_smooth(method="lm", color = "black") +
  xlab("Predicted half-life") + ylab("Observed half-life")
ggsave(paste0(outdir,"corr_grap", ".jpeg"))


# Generate bed files ------------------------------------------------------

number_a549 = round(nrow(residual_table_a549)*0.20)
number_hep = round(nrow(residual_table_hep)*0.20)
out_a549 <- "/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/a549_slam/motif_analysis/"
dir.create(paste0(out_a549,"one_transcript"))
out_hep <- "/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/hepg2_slam/motif_analysis/"
dir.create(paste0(out_hep,"one_transcript"))
#residual_table_a549 %>% arrange(halflife) %>% tail(n = number) %>% head()
#residual_table_a549 %>% arrange(halflife) %>% head(n = number) %>% tail()

#Half-life ranking
residual_table_a549 %>% arrange(halflife) %>% head(n = number_a549) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/one_transcript/","halflife_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
residual_table_a549 %>% arrange(halflife) %>% tail(n = number_a549) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/one_transcript/","halflife_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

residual_table_hep %>% arrange(halflife) %>% head(n = number_hep) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/one_transcript/","halflife_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
residual_table_hep %>% arrange(halflife) %>% tail(n = number_hep) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/one_transcript/","halflife_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Regressed out ranking 
residual_table_a549 %>% arrange(Residuals) %>% head(n = number_a549) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/one_transcript/","residual_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
residual_table_a549 %>% arrange(Residuals) %>% tail(n = number_a549) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/one_transcript/","residual_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

residual_table_hep %>% arrange(Residuals) %>% head(n = number_hep) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/one_transcript/","residual_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
residual_table_hep %>% arrange(Residuals) %>% tail(n = number_hep) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/one_transcript/","residual_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Background
residual_table_a549 %>% arrange(halflife) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/one_transcript/","background.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

residual_table_hep %>% arrange(halflife) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/one_transcript/","background.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#FIRE
residual_table_a549 %>% dplyr::select(transcript_id, halflife) %>%
  write.table(paste0(out_a549, "/one_transcript/","fire_halflife.txt"), row.names = F, quote = F, sep = "\t")
residual_table_a549 %>% dplyr::select(transcript_id, Residuals) %>%
  write.table(paste0(out_a549, "/one_transcript/","fire_residual.txt"), row.names = F, quote = F, sep = "\t")
residual_table_a549 %>%
  bind_cols(rep(".", n = nrow(residual_table_a549))) %>%
  filter(length < 10000 & length > 6) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  arrange(transcript_id, start) %>%
  write.table(paste0(out_a549,"/one_transcript/","fire.bed"), row.names = F, quote = F, sep = "\t", col.names = F)

residual_table_hep %>% dplyr::select(transcript_id, halflife) %>%
  write.table(paste0(out_hep,"/one_transcript/","fire_halflife.txt"), row.names = F, quote = F, sep = "\t")
residual_table_hep %>% dplyr::select(transcript_id, Residuals) %>%
  write.table(paste0(out_hep,"/one_transcript/","fire_residual.txt"), row.names = F, quote = F, sep = "\t")
residual_table_hep %>%
  bind_cols(rep(".", n = nrow(residual_table_hep))) %>%
  filter(length < 10000 & length > 6) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  arrange(transcript_id, start) %>%
  write.table(paste0(out_hep,"/one_transcript/","fire.bed"), row.names = F, quote = F, sep = "\t", col.names = F)


# Get Residuals  for all transcripts -----------------------------------------------------------

#Let's predict for all
#Getting rid of transcripts with double 3'UTR entries
stability_a549_all <- read.delim("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/a549_slam/slam/halflife_2rep_aboveCPM//halflife_filtered.tsv",
                                 header = T) %>%
  dplyr::select(chr, start, end, transcript_id, length, strand, halflife) %>%
  distinct()

stability_hep_all <- read.delim("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/hepg2_slam/slam_downsamplingbg//halflife_2rep_aboveCPM//halflife_filtered.tsv",
                                header = T) %>% 
  dplyr::select(chr, start, end, transcript_id, length, strand, halflife) %>% 
  distinct()


#Spearman R
jpeg(paste0(outdir, "spearman_a549_hep_all", ".jpeg"))
merge(stability_a549_all, stability_hep_all, by = c(1:6), 
      suffixes = c("_a549", "_hep")) %>% 
  ggscatter(x = "halflife_hep", y = "halflife_a549", 
            add = "reg.line", conf.int = TRUE, size = 0.5,
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Half-life A549", ylab = "Half-life HepG2")
dev.off()

# stability_a549_all %>%
#   filter(transcript_id %in% c("ENST00000644137","ENST00000644652",
#                               "ENST00000370128", "ENST00000368935",
#                               "ENST00000400479", "ENST00000262160"))   %>% 
#   bind_cols(rep(".", n = nrow(stability_a549_all))) %>%
#   dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
#   arrange(transcript_id, start) %>%
#   write.table("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/bed_testing/test_sorted.bed", row.names = F, quote = F, sep = "\t", col.names = F)
# 
# stability_a549_all %>%
#   filter(transcript_id %in% c("ENST00000644137","ENST00000644652",
#                               "ENST00000370128", "ENST00000368935",
#                               "ENST00000400479", "ENST00000262160"))   %>% 
#   bind_cols(rep(".", n = nrow(stability_a549_all))) %>%
#   dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
#   #arrange(transcript_id, start) %>%
#   write.table("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/bed_testing/test.bed", row.names = F, quote = F, sep = "\t", col.names = F)

# 
# 
# stability_a549_all %>% filter(duplicated(stability_a549_all$transcript_id)) %>% arrange(transcript_id)
# stability_a549_all %>% filter(transcript_id %in% c("ENST00000644137","ENST00000644652"))
# 
# duplicated(stability_a549_all$transcript_id)

###Creating merged table### 
###########################

stabilities_all <- bind_rows(a549 = stability_a549_all, hep = stability_hep_all, .id = "cell_line") %>% 
  dplyr::select(cell_line, transcript_id, halflife) 

#define one-hot encoding function
dummy <- dummyVars("~ cell_line", data=stabilities_all)

#perform one-hot encoding on data frame
final_stab_all <- data.frame(transcript_id = stabilities_all$transcript_id,
                             halflife = stabilities_all$halflife,
                             predict(dummy, newdata=stabilities_all))

#prep codon freq
selected_cf_all <- hg38_cf[final_stab_all$transcript_id,]

#Normalize
normed_cf_all <- sweep(selected_cf_all, 1, rowSums(selected_cf_all), "/")
# Add the CDS length of the condon frequencies
normed_cf_all <- cbind(normed_cf_all, rowSums(selected_cf_all))
colnames(normed_cf_all)[ncol(normed_cf_all)] <- "cds_length"

final_stab_all <- merge(final_stab_all, normed_cf_all, by.x = "transcript_id", by.y = 0)
#Normalize half-lives
final_stab_all$halflife <- log(final_stab_all$halflife)


###Interaction terms###
#######################

inter_a549 <- final_stab_all[,"cell_linea549"] * final_stab_all[,5:ncol(final_stab_all)]
colnames(inter_a549) <- paste0(colnames(inter_a549), "_ia549")
inter_hep <- final_stab_all[,"cell_linehep"] * final_stab_all[,5:ncol(final_stab_all)]
colnames(inter_hep) <- paste0(colnames(inter_hep), "_ihep")

final_stab_all <- cbind(final_stab_all, inter_a549, inter_hep)

#Scale
final_stab_all[,5:ncol(final_stab_all)] <- scale(final_stab_all[,5:ncol(final_stab_all)])

#Predict
pred_bestModel_all <- predict(model1, 
                              newdata = dplyr::select(final_stab_all, - c(halflife, transcript_id)))


pred_bestModel_all <- cbind(final_stab_all[,1:4],
                            pred_bestModel_all)

#Residuals
pred_bestModel_all$Residuals <- pred_bestModel_all$halflife - pred_bestModel_all$pred_bestModel_all 
colnames(pred_bestModel_all)[5] <- "Predictions"
colnames(pred_bestModel_all)[2] <- "LogHalflife"
pred_bestModel_all$halflife <- exp(pred_bestModel_all$LogHalflife)

comparison_testing_all_only_one_transcript <- merge(pred_bestModel_all, final_residual_table, by = 1)
cor.test(comparison_testing_all_only_one_transcript$Predictions.x, comparison_testing_all_only_one_transcript$Predictions.y)
#0.777
jpeg(paste0(outdir, "pearson_one_transcript_all_transcript", ".jpeg"))
comparison_testing_all_only_one_transcript %>% 
  ggscatter(x = "Predictions.x", y = "Predictions.y", 
            add = "reg.line", conf.int = TRUE, size = 0.5,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Predictions one transcript", ylab = "Predictions all transcripts")
dev.off()


residual_table_a549_all <- pred_bestModel_all %>%
  filter(cell_linea549 == 1) %>%
  inner_join(y = (stability_a549_all %>% 
                    dplyr::select(chr, start, end, transcript_id, length, strand, length)), 
             by = "transcript_id")



residual_table_hep_all <- pred_bestModel_all %>%
  filter(cell_linehep == 1) %>%
  inner_join(y = (stability_hep_all %>% 
                    dplyr::select(chr, start, end, transcript_id, length, strand, length)), 
             by = "transcript_id")

residual_table_a549_all %>% arrange(halflife) %>% 
  dplyr::select(- c(cell_linea549, cell_linehep)) %>%
  write.table(paste0(outdir,"stab_residual_a549_all.tsv"),quote = F, sep = "\t", row.names = F)

residual_table_hep_all %>% arrange(halflife) %>%
  dplyr::select(- c(cell_linea549, cell_linehep)) %>%
  write.table(paste0(outdir,"stab_residual_hep_all.tsv"),quote = F, sep = "\t", row.names = F)

residual_table_hep_all <- read_delim(paste0(outdir,"stab_residual_hep_all.tsv"))
residual_table_a549_all <- read_delim(paste0(outdir,"stab_residual_a549_all.tsv"))
#Bed files

number_a549 = round(nrow(residual_table_a549_all)*0.20)
number_hep = round(nrow(residual_table_hep_all)*0.20)
dir.create(paste0(out_a549,"all_transcripts"))
dir.create(paste0(out_hep,"all_transcripts"))


#Half-life ranking
residual_table_a549_all %>% arrange(halflife) %>% head(n = number_a549) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/all_transcripts/","halflife_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
residual_table_a549_all %>% arrange(halflife) %>% tail(n = number_a549) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/all_transcripts/","halflife_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

residual_table_hep_all %>% arrange(halflife) %>% head(n = number_hep) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/all_transcripts/","halflife_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
residual_table_hep_all %>% arrange(halflife) %>% tail(n = number_hep) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep_all))) %>% 
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/all_transcripts/","halflife_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Regressed out ranking 
residual_table_a549_all %>% arrange(Residuals) %>% head(n = number_a549) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/all_transcripts/","residual_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
residual_table_a549_all %>% arrange(Residuals) %>% tail(n = number_a549) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/all_transcripts/","residual_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

residual_table_hep_all %>% arrange(Residuals) %>% head(n = number_hep) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/all_transcripts/","residual_lowstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")
residual_table_hep_all %>% arrange(Residuals) %>% tail(n = number_hep) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/all_transcripts/","residual_highstab.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#Background
residual_table_a549_all %>% arrange(halflife) %>%
  bind_cols(rep(".", n = nrow(residual_table_a549_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_a549, "/all_transcripts/","background.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

residual_table_hep_all %>% arrange(halflife) %>%
  bind_cols(rep(".", n = nrow(residual_table_hep_all))) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  write.table(paste0(out_hep, "/all_transcripts/","background.bed"), col.names = F, quote = F, row.names = F, sep = "\t")

#FIRE
residual_table_a549_all <- residual_table_a549_all[!duplicated(residual_table_a549_all$transcript_id), ] 
residual_table_a549_all %>% dplyr::select(transcript_id, halflife) %>% distinct() %>% 
  write.table(paste0(out_a549, "/all_transcripts/","fire_halflife.txt"), row.names = F, quote = F, sep = "\t")
residual_table_a549_all %>% dplyr::select(transcript_id, Residuals) %>%  distinct() %>%
  write.table(paste0(out_a549, "/all_transcripts/","fire_residual.txt"), row.names = F, quote = F, sep = "\t")
residual_table_a549_all %>% distinct() %>%
  bind_cols(rep(".", n = nrow(residual_table_a549_all))) %>%
  filter(length < 10000 & length > 6) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  arrange(transcript_id, start) %>%
  write.table(paste0(out_a549,"/all_transcripts/","fire.bed"), row.names = F, quote = F, sep = "\t", col.names = F)

residual_table_hep_all <- residual_table_hep_all[!duplicated(residual_table_hep_all$transcript_id), ] 
residual_table_hep_all %>% dplyr::select(transcript_id, halflife) %>%  distinct() %>%
  write.table(paste0(out_hep,"/all_transcripts/","fire_halflife.txt"), row.names = F, quote = F, sep = "\t")
residual_table_hep_all %>% dplyr::select(transcript_id, Residuals) %>%  distinct() %>%
  write.table(paste0(out_hep,"/all_transcripts/","fire_residual.txt"), row.names = F, quote = F, sep = "\t")
residual_table_hep_all %>% distinct() %>%
  bind_cols(rep(".", n = nrow(residual_table_hep_all))) %>%
  filter(length < 10000 & length > 6) %>%
  dplyr::select(chr, start, end, transcript_id, paste0("...", ncol(.)),strand) %>%
  arrange(transcript_id, start) %>%
  write.table(paste0(out_hep,"/all_transcripts/","fire.bed"), row.names = F, quote = F, sep = "\t", col.names = F)



