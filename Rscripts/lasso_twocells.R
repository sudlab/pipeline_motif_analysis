#Final correct script, with one or all transcript bed generator

library(glmnet)
library(tidyverse)
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(caret)


###INPUTS###
############
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

stability_a549 <- read.delim("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/a549_slam/slam/halflife_2rep_aboveCPM/halflife_filtered.tsv",
                             header = T) %>%
  dplyr::select(chr, start, end, transcript_id, length, strand, halflife) %>%
  inner_join(cds_lengths) %>%
  inner_join(tx2gene %>% distinct(), 
             by = c("transcript_id" = "transcript_stable_id")) %>% 
  group_by(gene_stable_id) %>%
  slice_max(order_by=cds_length, n=1, with_ties = FALSE)  %>% ungroup()

stability_hep <- read.delim("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/hepg2_slam/slam_downsamplingbg/halflife_2rep_aboveCPM/halflife_filtered.tsv",
                            header = T) %>%
  dplyr::select(chr, start, end, transcript_id, length, strand, halflife) %>%
  inner_join(cds_lengths) %>%
  inner_join(tx2gene %>% distinct(), 
             by = c("transcript_id" = "transcript_stable_id")) %>% 
  group_by(gene_stable_id) %>%
  slice_max(order_by=cds_length, n=1, with_ties = FALSE)  %>% ungroup()

merge(stability_a549, stability_hep, by = c(1:6), 
      suffixes = c("_a549", "_hep")) %>% nrow()

#Mean h-l
mean(stability_a549$halflife)
median(stability_a549$halflife)
mean(stability_hep$halflife)
median(stability_hep$halflife)

#Spearman R
jpeg(paste0(outdir, "pearson_a549_hep", ".jpeg"))
merge(stability_a549, stability_hep, by = c(1:6), 
      suffixes = c("_a549", "_hep")) %>% 
  ggscatter(x = "halflife_hep", y = "halflife_a549", 
            add = "reg.line", conf.int = TRUE, size = 0.5,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Half-life A549", ylab = "Half-life HepG2") +
  theme_classic(base_size = 16) 
dev.off()

#H-l density plot 
jpeg(paste0(outdir, "density_a549_hep", ".jpeg"))
merge(stability_a549, stability_hep, by = c(1:6), 
      suffixes = c("_A549", "_HepG2")) %>%
  dplyr::select(transcript_id, halflife_A549, halflife_HepG2) %>%
  melt() %>%
  ggplot(aes(x=value, fill=variable)) +
  geom_density(alpha=.25) +
  theme_classic(base_size = 16) +
  labs(fill = "Cell line", x = "Half-life")
dev.off()

#Number shared
merge(stability_a549, stability_hep, by = c(1:6), 
      suffixes = c("_A549", "_HepG2")) %>%
  dplyr::select(transcript_id, halflife_A549, halflife_HepG2) %>% nrow()

###Creating merged table### 
###########################

stabilities <- bind_rows(a549 = stability_a549, hep = stability_hep, .id = "cell_line") %>% 
  dplyr::select(cell_line, transcript_id, halflife) 

#Check duplicates 
if (stabilities %>% filter(cell_line == 1) %>%
  group_by(transcript_id) %>%
  filter(n()>1) %>% nrow() != 0 ) {
  print("Duplicates in stabilities data")
}

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

sampling = c(0.55, 0.6, 0.65, 0.7, 0.75)
summary_models <- matrix(ncol=13, nrow= 0)

for (i in sampling) {
  
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
  #summary(model1)
  how_many <- round((nrow(test_data)/(nrow(final_data)) *100), 0)
  df_k_folds_stats <- as.data.frame(model1$resample)
  df_k_folds_stats %>% 
    ggplot(aes(RMSE,Rsquared)) + 
    geom_point() + geom_text(aes(label=Resample),hjust=0, vjust=0)
  ggsave(paste0(outdir,"stats_training_folds_", how_many, ".jpeg"))
  write.table(df_k_folds_stats, paste0(outdir,"performances_traning_folds_", how_many,".tsv"),
              quote = F, sep = "\t ", col.names = F)
  
  # Optimal tuning parameters
  lambda_kfold = model1$bestTune$lambda
  #Saving lamdba plot
  jpeg(paste0(outdir,"lambdaplot_model_", how_many, ".jpeg"))
  plot(log(model1$results$lambda), model1$results$RMSE,
       xlab = "log(lambda)", ylab="RMSE",
       xlim=c(-6, -1))
  abline(v=log(lambda_kfold), col="coral")
  dev.off()
  
  #Predict
  pred_model1 <- predict(model1, 
                         newdata = dplyr::select(test_data, - c(halflife, transcript_id)))
  pred_model1 <- cbind(test_data[,1:4], pred_model1)
  model1_perf <- data.frame(RMSE=RMSE(pred_model1$pred_model1, pred_model1$halflife),
                            Rsquared=R2(pred_model1$pred_model1, pred_model1$halflife))
  
  #Create summary table
  get_stats_best_model_stats <- df_k_folds_stats[which.min(df_k_folds_stats$RMSE),]
  summary_model1 <- cbind(model_sampling = paste0("model_", how_many),
                          training_perf_model1, 
                          Lambda = lambda_kfold,
                          RMSE_bestModel = get_stats_best_model_stats$RMSE,
                          Rsquared_bestModel = get_stats_best_model_stats$Rsquared,
                          model1_perf, 
                          PercentTraining = nrow(train_data)/(nrow(final_data)) *100,
                          PercentTesting = nrow(test_data)/(nrow(final_data)) *100,
                          Number_testing = nrow(test_data),
                          Number_training = nrow(train_data),
                          TraningRsquared_minus_TestRsuared = (get_stats_best_model_stats$RMSE)-(model1_perf[1,"Rsquared"]))
  colnames(summary_models) <- colnames(summary_model1)
  summary_models <- rbind(summary_models, summary_model1)
}


write.table(summary_models, paste0(outdir,"performances_summary", ".tsv"),
            quote = F, sep = "\t ", row.names = F)


###Determine best sampling###
#############################

summary_models %>%
  ggplot(aes(Number_training, Rsquared)) + 
  geom_point() + geom_text(aes(label=model_sampling),hjust=0, vjust=0)
ggsave(paste0(outdir,"Ntraining_vs_Rsquared", ".jpeg"))

summary_models %>%
  ggplot(aes(TraningRsquared_minus_TestRsuared,Rsquared)) + 
  geom_point() + geom_text(aes(label=model_sampling),hjust=0, vjust=0)
ggsave(paste0(outdir,"DiffRsquared_vs_Rsquared", ".jpeg"))

summary_models %>%
  ggplot(aes(RMSE,Rsquared)) + 
  geom_point() + geom_text(aes(label=model_sampling),hjust=0, vjust=0) +
  xlab("RMSE") + ylab("Rsquared")
ggsave(paste0(outdir,"Rsquared_vs_RMSE", ".jpeg"))


