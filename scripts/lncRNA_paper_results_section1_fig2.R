###########################################################################
###########################################################################
###                                                                     ###
###      Code to accompany paper Results section 1: Bone comparisons    ###     
### 								                                                    ###
###########################################################################
###########################################################################


library(DESeq2)
library(pheatmap)
library(data.table)
library(tidyverse)
library(ggplot2)
library(plotrix)

#1) Perform differential expression analysis between bone fragments and bone biopsy 
CountData_time_series <- read.csv("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/diff_exp_analysis/gene_count_matrix.csv", row.names = 1)
ColData_time_series <- read.csv("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/diff_exp_analysis/phenotypes.csv", row.names = 1)
reorder_idx_time_series <- match(rownames(ColData_time_series),colnames(CountData_time_series))
# Reorder the columns of the count data
reordered_CountData_time_series <- CountData_time_series[ , reorder_idx_time_series]
all(rownames(ColData_time_series) %in% colnames(CountData_time_series))
#We also see a warning message, where our condition was converted to a factor these messages are ok to see!https://angus.readthedocs.io/en/2019/diff-ex-and-viz.html 
dds_time_series <- DESeqDataSetFromMatrix(countData =  reordered_CountData_time_series,
                                          colData =  ColData_time_series,
                                          design = ~ condition)
dds_time_series <- estimateSizeFactors(dds_time_series)
# Extract the normalized counts
normalized_counts_time_series <- counts(dds_time_series, normalize = TRUE)
# Transform the normalized counts 
vsd_time_series <- vst(dds_time_series, blind = TRUE)
# Extract the matrix of transformed counts
vsd_mat_time_series <- assay(vsd_time_series)
# Compute the correlation values between samples
vsd_cor_time_series <- cor(vsd_mat_time_series, vsd_mat_time_series)
# Plot the PCA of PC1 and PC2
plotPCA(vsd_time_series, intgroup= "condition")
dds <- DESeq(dds_time_series,minReplicatesForReplace=Inf)
# Plot dispersions
plotDispEsts(dds)
# Subset the results to only return the significant genes with p-adjusted values less than 0.05
IH_lit <- results(dds, 
                  contrast = c("condition","IH", "lit"), 
                  alpha = 0.05)
summary(IH_lit)
plotMA(IH_lit, ylim=c(-2,2))
IH_lit_sig <- subset(IH_lit, padj < 0.05)
IH_lit_sig_order <- IH_lit_sig[order(IH_lit_sig$log2FoldChange),]
top_1000_lit <- head(IH_lit_sig_order, n=1000)
top_1000_IH <- tail(IH_lit_sig_order, n=1000)
# Extract the 1000 differentially expressed genes to identify GO terms with Panther db
foo <- data.frame(do.call('rbind', strsplit(as.character(row.names(top_1000_IH)),'|',fixed=TRUE)))
top_1000_IH_id <- as.data.frame(foo$X2)
#write.table(top_1000_IH_id, "/scratch/aa9gj/top_1000_IH_DESEQ", quote=F)

foo <- data.frame(do.call('rbind', strsplit(as.character(row.names(top_1000_lit)),'|',fixed=TRUE)))
top_1000_lit_id <- as.data.frame(foo$X2)
#write.table(top_1000_lit_id, "/scratch/aa9gj/top_1000_lit_DESEQ", quote=F, row.names = F)

#2) Perform correaltion analysis (Figure 2A)
all_17 <- fread("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/ftp_download/all17.abund")
all_58 <- fread("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/ftp_download/all58.abund")
all_17 <- filter(all_17, all_17$`Gene Name` != "-")
all_58 <- filter(all_58, all_58$`Gene Name` != "-")
both <- full_join(all_17, all_58, by = "Gene Name")
both$tpm_use.x <- log2(both$TPM.x+1)
both$tpm_use.y <- log2(both$TPM.y+1)
both$tpm_use.x <- replace_na(both$tpm_use.x, 0)
both$tpm_use.y <- replace_na(both$tpm_use.y, 0)
cor(both$tpm_use.x, both$tpm_use.y)
cor.test(both$tpm_use.x, both$tpm_use.y)
ggplot(both, aes(x=log2(TPM.x)+1, y=log2(TPM.y)+1)) + geom_point() + geom_abline(col="blue") + xlab("This study (log2(TPM)+1)") + ylab("Farr et al. (log2(TPM)+1)") + theme_classic()

#3) Perform bone RNA comparison with genes from the amazing osteocyte paper
setwd("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/diff_exp_analysis/amazing_osteocyte")
counts_files <- read.delim("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/diff_exp_analysis/files", header = FALSE)
setDT(counts_files)
IH_list <- as.data.frame(counts_files[V1 %like% "IH"])
lit_list <-as.data.frame(counts_files[V1 %like% "lit"]) 
counts_data_IH <- list()
for (i in seq_along(IH_list$V1)) {
  counts_data_IH[[i]] <- read.delim(file = IH_list$V1[i], header=F)
}

counts_data_lit <- list()
for (i in seq_along(IH_list$V1)) {
  counts_data_lit[[i]] <- read.delim(file = lit_list$V1[i], header=F)
}

counts_data_IH
counts_data_lit
##Problems with RUNX2 and PHEX (grepped AS1 with them) so a quick fix is this
# counts_data_IH[[13]]$V2 <- gsub("RUNX2-AS1", "no", counts_data_IH[[13]]$V2)
# counts_data_IH[[13]] <- filter(counts_data_IH[[13]], V2 != "no")
# counts_data_lit[[13]]$V2 <- gsub("RUNX2-AS1", "no", counts_data_lit[[13]]$V2)
# counts_data_lit[[13]] <- filter(counts_data_lit[[13]], V2 != "no")
# 
# counts_data_IH[[12]]$V2 <- gsub("PHEX-AS1", "no", counts_data_IH[[12]]$V2)
# counts_data_IH[[12]] <- filter(counts_data_IH[[12]], V2 != "no")
# counts_data_lit[[12]]$V2 <- gsub("PHEX-AS1", "no", counts_data_lit[[12]]$V2)
# counts_data_lit[[12]] <- filter(counts_data_lit[[12]], V2 != "no")

##do the analysis function
prepare_genes <- function(x,y) {
  x$event <- rep('in_house', length(x$V1))
  x$mean <- mean(x$V9)
  x$sd <- sd(x$V9)
  x$se <- std.error(x$V9)
  y$event <- rep('lit', length(y$V1))
  y$mean <- mean(y$V9)
  y$sd <- sd(y$V9)
  y$se <- std.error(y$V9)
  gene <- rbind(x, y)
  return(gene)
}

prepare_gene_list <- list()
for (i in seq_along(IH_list$V1)) {
  prepare_gene_list[[i]] <- prepare_genes(counts_data_IH[[i]], counts_data_lit[[i]])
}

file_names <- gsub("_IH", "",IH_list$V1)
file_names <- as.data.frame(file_names)
split_df_counts<-function(list){
  for (i in 1:length(list)){
    assign(file_names$file_names[i], list[[i]], envir = .GlobalEnv)
  }
}
split_df_counts(prepare_gene_list)

## late osteoblast/osteocyte markers enrichment figure (Figure 2B)
BGLAP$V2 <- gsub("PMF1-BGLAP", "BGLAP", BGLAP$V2)
test <- rbind(SOST, PHEX, MEPE,DMP1, PDPN, RUNX2, SP7, ALPL, BGLAP, CSNK2A1, DSTN, MMP14, CAPG, HYOU1, FGF23)
plot <- ggplot(test, aes(as.factor(V2), log2(mean+1), fill=event))
plot <- plot + geom_bar(stat = "identity", position = 'dodge')
plot
unique(test$V2)
test$V2 <- factor(test$V2,levels = c("FGF23", "PHEX", "SOST","PDPN", "MEPE","DMP1", "SP7","RUNX2", "ALPL", "CSNK2A1", "DSTN", "MMP14", "HYOU1","CAPG", "BGLAP"))
plot <- ggplot(test, aes(x=as.factor(V2), y=log2(mean+1), fill=event)) +
  geom_bar(position=position_dodge(), stat = "identity", colour='black') + 
  geom_errorbar(aes(ymin=log2(mean-se+1), ymax=log2(mean+se+1)), width=.2, position=position_dodge(.9)) + 
  xlab("Osteoblast/Osteocyte representative genes") + ylab("Average expression (log2(TPM+1)") + theme_classic()+theme(axis.text.x = element_text(angle = 45))
plot

# Bonemarrow enrichment figure

setwd("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/diff_exp_analysis/bonemarrow_genes/")
counts_files <- read.delim("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/diff_exp_analysis/bonemarrow_files.txt", header = FALSE)
setDT(counts_files)
IH_list <- as.data.frame(counts_files[V1 %like% "IH"])
lit_list <-as.data.frame(counts_files[V1 %like% "lit"]) 
counts_data_IH <- list()
for (i in seq_along(IH_list$V1)) {
  counts_data_IH[[i]] <- read.delim(file = IH_list$V1[i], header=F)
}

counts_data_lit <- list()
for (i in seq_along(IH_list$V1)) {
  counts_data_lit[[i]] <- read.delim(file = lit_list$V1[i], header=F)
}

counts_data_IH
counts_data_lit

##do the analysis function
prepare_genes <- function(x,y) {
  x$event <- rep('in_house', length(x$V1))
  x$mean <- mean(x$V9)
  x$sd <- sd(x$V9)
  x$se <- std.error(x$V9)
  y$event <- rep('lit', length(y$V1))
  y$mean <- mean(y$V9)
  y$sd <- sd(y$V9)
  y$se <- std.error(y$V9)
  gene <- rbind(x, y)
  return(gene)
}

prepare_gene_list <- list()
for (i in seq_along(IH_list$V1)) {
  prepare_gene_list[[i]] <- prepare_genes(counts_data_IH[[i]], counts_data_lit[[i]])
}

file_names <- gsub("_IH", "",IH_list$V1)
file_names <- as.data.frame(file_names)
split_df_counts<-function(list){
  for (i in 1:length(list)){
    assign(file_names$file_names[i], list[[i]], envir = .GlobalEnv)
  }
}
split_df_counts(prepare_gene_list)

# Figure for paper (Figure 2C)
test <- rbind(AZU1, CTSG, DEFA1B,DEFA1,DEFA4, ELANE, GYPA, HBB, HBD, MPO, OR10Z1, RHAG)
plot <- ggplot(test, aes(as.factor(V2), log2(mean+1), fill=event))
plot <- plot + geom_bar(stat = "identity", position = 'dodge')
plot
unique(test$V2)
test$V2 <- factor(test$V2,levels = c("OR10Z1","RHAG","GYPA","AZU1", "CTSG","DEFA4", "MPO", "ELANE", "HBD","DEFA1B","DEFA1","HBB"))
plot <- ggplot(test, aes(x=as.factor(V2), y=log2(mean+1), fill=event)) +
  geom_bar(position=position_dodge(), stat = "identity", colour='black') + 
  geom_errorbar(aes(ymin=log2(mean-se+1), ymax=log2(mean+se+1)), width=.2, position=position_dodge(.9)) + 
  xlab("Bone marrow representative genes") + ylab("Average expression (log2(TPM+1)")+theme_classic() + theme(axis.text.x = element_text(angle = 45))
plot

# Ratio analysis + Figure 2D-E
all_17 <- fread("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/ftp_download/all17.abund")
all_58 <- fread("/scratch/aa9gj/wrk_backup/lncRNA_paper/bone_RNA_comparison/ftp_download/all58.abund")
all_17 <- filter(all_17, all_17$`Gene Name` != "-")
all_58 <- filter(all_58, all_58$`Gene Name` != "-")
all_17 <- filter(all_17, all_17$TPM != 0)
all_58 <- filter(all_58, all_58$TPM != 0)
both <- inner_join(all_17, all_58, by = "Gene Name")
both$tpm_use.x <- both$TPM.x
both$tpm_use.y <- both$TPM.y
both$ratio <- both$tpm_use.x / both$tpm_use.y
head(both)
both <- both[order(both$ratio, decreasing = T),]

osteocyte_transcriptome <- fread("/scratch/aa9gj/osteocyte_transcriptome.csv", header = T)
head(osteocyte_transcriptome)
osteocyte_transcriptome <- filter(osteocyte_transcriptome, Biotype == "protein_coding")
osteocyte_transcriptome <- filter(osteocyte_transcriptome, Human_GeneSymbol != "NA")#18309 pc genes
osteocyte_transcriptome <- filter(osteocyte_transcriptome, LFC != "NA")#12608 pc genes with LFC value
full_table <- inner_join(osteocyte_transcriptome, both, by = c("Human_GeneSymbol" = "Gene Name"))#12187 genes
colnames(full_table)
full_table <- full_table[, c(1,7,42,46:48,51,78:80)]
osteocyte_sig <- filter(full_table, Osteocyte_transcriptome_signature == TRUE)
ggplot(osteocyte_sig, aes(ratio))+geom_histogram()
cor.test(log2(osteocyte_sig$ratio), osteocyte_sig$LFC)
par(mfrow=c(1,1))
ggplot(osteocyte_sig, aes(LFC, log2(ratio)))+geom_point()
mean(osteocyte_sig$ratio)
median(osteocyte_sig$ratio)
subset(osteocyte_sig, ratio >= 1)
osteocyte_sig$tissue <- rep("osteocyte", length(osteocyte_sig$Ensembl_ID))

osteocyte_nosig <- filter(full_table, Osteocyte_transcriptome_signature == FALSE)
cor.test(log2(osteocyte_nosig$ratio), osteocyte_nosig$LFC)
mean(osteocyte_nosig$ratio)
median(osteocyte_nosig$ratio)

blood <- filter(full_table, Ayturk_blood == TRUE)
cor.test(log2(blood$ratio), blood$LFC)
mean(blood$ratio)
median(blood$ratio)

marrow <- filter(full_table, Ayturk_marrow == TRUE)
cor.test(log2(marrow$ratio), marrow$LFC)
mean(marrow$ratio)
median(marrow$ratio)
subset(marrow, ratio < 1)

marrow$tissue <- rep('marrow', length(marrow$Ensembl_ID))

muscle <- filter(full_table, Ayturk_muscle == TRUE)
cor.test(log2(full_table$ratio), full_table$LFC)
mean(muscle$ratio)
median(muscle$ratio)

full_table
#ggplot(test, aes(x=Human_GeneSymbol, y=ratio, color = as.factor(Osteocyte_transcriptome_signature))) + geom_point()

test <- full_table

test$Osteocyte_transcriptome_signature<- ifelse(test$Osteocyte_transcriptome_signature==TRUE, 1 , 0)

ggplot(full_table, aes(log2(ratio), color=as.factor(Osteocyte_transcriptome_signature)))+geom_density()+
  theme_classic()+ labs(color = "Osteocyte signature") + xlab("Log2(Ratio of Expression)")+geom_vline(xintercept = 0)

test_ks <- ks.test(full_table$ratio[full_table$Osteocyte_transcriptome_signature == TRUE],
                   full_table$ratio[full_table$Osteocyte_transcriptome_signature == FALSE])
test_ks
wilcox.test(ratio ~ Osteocyte_transcriptome_signature, data = full_table)

ggplot(full_table, aes(log2(ratio), color=as.factor(Ayturk_marrow)))+geom_density()+theme_classic()+
  labs(color = "Bone marrow enriched")+ xlab("Log2(Ratio of Expression)")+geom_vline(xintercept = 0)
wilcox.test(ratio ~ Ayturk_marrow, data = full_table)


ggplot(full_table, aes(row.names(full_table),log2(ratio),color=as.factor(Ayturk_marrow)))+geom_point()
ggplot(full_table, aes(row.names(full_table),log2(ratio),color=as.factor(Osteocyte_transcriptome_signature)))+geom_point()
