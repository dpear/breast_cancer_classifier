library("tidyverse")
# in case you cannot install tidyverse package, please uncomment the packages below
# library("readr")
# library("tibble")
# library("dplyr")
# library("tidyr")
# library("ggplot2")
library("magrittr")
library("DESeq")
library("e1071")
library("ROCR")

# load data ---------------------------------------------------------------
# load the expression level data. TPM and readCoutns, the first row shows the sample ID
tpm <- read_tsv("data/tpm_96_nodup.tsv", col_names = T)
readcounts <- read_tsv("data/readcounts_96_nodup.tsv", col_names = T)

# load sample meta data file
# the column "sampleID" matchs the colnames of variable tpm and readcounts
# the column "recurStatus" shows the recurrence status of the patient, R = recurrence, N = non-recurrence 
patient_info <- read_csv("data/patient_info.csv", col_names = T) # recurrence status file

# generate sample_recurStatus dataframe with the same sample order in tpm and reacounts
sampel_recurStatus <- patient_info %>% arrange(match(sampleID, colnames(readcounts)[-1])) %>%
  select(sampleID, recurStatus)

# load the pre-selected marker gene list
preselectedList <- readLines("data/preselectedList") # breast cancer biomarkers

# differential expression (DE) analysis of the biomarker genes ---------------
# we use DESeq to calculate the DE p-value of each biomarker gene, 
# then rank the marker genes by their p-values

cds <- newCountDataSet(countData = readcounts %>% filter(gene_id %in% preselectedList) %>% column_to_rownames("gene_id"), 
                       conditions = sampel_recurStatus$recurStatus)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds, "R", "N")
rank_marker_gene <- res %>% filter(id %in% preselectedList) %>%
  arrange(pval) %>% .$id

# calculate and plot Fig5C in PNAS paper ----------------------------------
# number of gene in biomarker gene sets selected by their rank in rank_marker_gene,
# top 30, 60, 90, ... 720, 750
geneNum <- seq(30, 750, 30) 

# calculate average AUC of SVM based on each subset of biomarker genes 
# need several minutes to complete the loop

averageAUC <- rep(0, length(geneNum)) 
for(i in 1:length(geneNum)){
  temp.tpm <- tpm %>% filter(gene_id %in% rank_marker_gene[1:geneNum[i]]) %>%
    column_to_rownames("gene_id") %>% t() %>% as.data.frame() %>%
    mutate(status = as.factor(ifelse(sampel_recurStatus$recurStatus == "R", 1, 0)))
  
  idx_pos <- which(temp.tpm$status == 1)
  idx_neg <- which(temp.tpm$status == 0)
  
  set.seed(123) # seed for random number, for reproducible results
  
  # random sampling of trainging and testing dataset, 100 times average auc
  auc <- sapply(1:100, FUN = function(x){
    # trainning set 20 pos vs 20 neg, test set 8 pos vs 48 neg
    idx_testset <- c(sample(idx_pos, 8),sample(idx_neg, 48))
    
    testset <- temp.tpm[idx_testset, ]
    trainset <- temp.tpm[-idx_testset, ]
    
    idx_zero_train <- which(colSums(trainset[, -ncol(trainset)]) == 0)
    if(length(idx_zero_train) > 0){
      trainset <- trainset[, -idx_zero_train]
      testset <- testset[, -idx_zero_train]
    }

    # use svm function in e1071 package to train SVM model,
    # and then use ROCR function to calculate ROC curve data and auc
    svm.model <- svm(status ~ ., data = trainset, cost = 5, kernel = 'sigmoid')
    svm.pred <- predict(svm.model, testset[, -ncol(testset)], decision.values = T)
    pred <- prediction(attr(svm.pred, "decision.values"), testset$status)
    perf_auc <- as.numeric(performance(pred, "auc")@y.values)

    return(perf_auc)
  })
  
  averageAUC[i] <- mean(auc)
}

# plot the number of biomarker against average AUC plot (Fig 5C)

data.frame(geneNum = geneNum, averageAUC = averageAUC) %>%
  ggplot(aes(x = geneNum, y = averageAUC)) +
  geom_point(shape=16, color = 'red',size=1.2) + 
  scale_y_continuous(limits = 0:1) +
  theme(plot.background = element_rect(fill = 'white', colour = 'grey',linetype = 'solid'),
        panel.background = element_rect(fill = 'white',color = 'grey',size=0.75,linetype='solid'),
        panel.grid.major = element_line(colour = "grey",size=0.2),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.75)) +
  labs(x = 'Number of genes', y = 'Average AUC',fill = '') + 
  theme(axis.text.x = element_text(size=8,hjust = 1),axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=11),axis.title.x = element_text(size=11),
        legend.text = element_text(size = 8)) 
ggsave('averageAUC.pdf', width = 4, height = 4)


# select top 215 biomarker genes to plot Fig 5D ---------------------------

temp.tpm <- tpm %>% filter(gene_id %in% rank_marker_gene[1:215]) %>%
  column_to_rownames("gene_id") %>% t() %>% as.data.frame() %>%
  mutate(status = as.factor(ifelse(sampel_recurStatus$recurStatus == "R", 1, 0)))

idx_pos <- which(temp.tpm$status == 1)
idx_neg <- which(temp.tpm$status == 0)

set.seed(123)
# random sampling of trainging and testing dataset, 3 times

roc_data <- c()
auc <- c()
for(i in 1:3){
  idx_testset <- c(sample(idx_pos, 8),sample(idx_neg, 48))
  
  testset <- temp.tpm[idx_testset, ]
  trainset <- temp.tpm[-idx_testset, ]
  
  idx_zero_train <- which(colSums(trainset[, -ncol(trainset)]) == 0)
  if(length(idx_zero_train) > 0){
    trainset <- trainset[, -idx_zero_train]
    testset <- testset[, -idx_zero_train]
  }
  
  # use svm function in e1071 package to train SVM model,
  # and then use ROCR function to calculate ROC curve data and auc
  svm.model <- svm(status ~ ., data = trainset, cost = 5, kernel = 'sigmoid')
  svm.pred <- predict(svm.model, testset[, -ncol(testset)], decision.values = T)
  pred <- prediction(attr(svm.pred, "decision.values"), testset$status)
  perf_roc <- performance(pred, "tpr", "fpr")
  roc_data <- rbind(roc_data, 
                    data.frame(runtime = as.factor(i), 
                               x = perf_roc@x.values[[1]], 
                               y = perf_roc@y.values[[1]]))
  perf_auc <- performance(pred, "auc")
  auc <- c(auc, as.numeric(perf_auc@y.values))
  
}

# plot roc curves from randomly generated roc_data

ggplot(roc_data, aes(x = x, y = y, group = runtime, col = runtime)) + geom_line() +
  labs(x = '1 - Specificity', y = 'Sensitivity',fill = '') +
  theme(plot.background = element_rect(fill = 'white', colour = 'grey',linetype = 'solid'),
        panel.background = element_rect(fill = 'white',color = 'grey',size=0.75,linetype='solid'),
        panel.grid.major = element_line(colour = "grey",size=0.2),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.75)) +
  theme(axis.text.x = element_text(size=8,hjust = 1),axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=11),axis.title.x = element_text(size=11),
        legend.position = "none") +
  scale_y_continuous(limits = 0:1) +
  scale_x_continuous(limits = 0:1) +
  annotate('text',label = paste0('Average AUC = ',format(round(mean(auc), 2), nsmall = 3)), 
           x = 0.65, y=0.2, size = 3)

ggsave('ROCplot.pdf', width = 4, height = 4)
column_to_rownames()
