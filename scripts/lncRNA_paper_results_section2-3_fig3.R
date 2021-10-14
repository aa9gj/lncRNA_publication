###########################################################################
###########################################################################
###                                                                     ###
###      Allelic Imbalance R analysis for 17 bone fragment samples      ###     
### 								                                                    ###
###########################################################################
###########################################################################
# Load packages and set working directory
library(data.table)
library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(liftOver)
setwd("/scratch/aa9gj/wrk_backup/lncRNA_paper/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/test_wasp/WASP_analysis/test_rename")
# Read in all 17 samples ASEreadcounter output from command line and split them to environment
counts_files <- read.delim("files", header = FALSE) #files is file names of interest
counts_data <- list()
for (i in seq_along(counts_files$V1)) {
  counts_data[[i]] <- read.delim(file = counts_files$V1[i], header=TRUE)
}
#split them to their respective files using appropriate sample names
file_names <- gsub("wasp_counter", "",counts_files$V1)
file_names <- as.data.frame(file_names)
split_df_counts<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "WASP_counts"), list[[i]], envir = .GlobalEnv)
  }
}
split_df_counts(counts_data)

# Perform basic histogram of the ratio of ref (or alt) reads to 
# total counts and filter sites accordingly
# In this example I decided to keep only sites with 
# ratios between 0.1 and 0.9 and did the same thing for all 17 samples
hist(a102_WASP_counts$refCount/a102_WASP_counts$totalCount)
hist(a117_WASP_counts$refCount/a117_WASP_counts$totalCount)

# Function to filter potentially bogus sites and perform binomial test

ASE_analysis <- function(x) {
  x<-filter(x, totalCount > 20)
  x<-filter(x, (refCount/totalCount) > .1,(refCount/totalCount) < .9)
  for (i in 1:length(x$contig)) {
    temp = binom.test(x$altCount[i], x$totalCount[i], p = 0.5)
    x$binom_p[i] = temp$p.value
    x$binom_success[i]= temp$estimate
    x$fdr[i] <- p.adjust(x$binom_p[i], method = "fdr", n=length(x$binom_p))
  }
  x_filterd<-filter(x, fdr < 0.05)
  return(x_filterd)
}

# Perform on all 17 samples

my_ASE_results <- list()
for (i in seq_along(counts_files$V1)) {
  my_ASE_results[[i]] <- ASE_analysis(counts_data[[i]])
  my_ASE_results[[i]] <- makeGRangesFromDataFrame(my_ASE_results[[i]], keep.extra.columns = TRUE, start.field = "position", end.field = "position", seqnames.field = "contig")
}

# Split results into their respective files

split_df_ASE<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "ASE_final_grange"), list[[i]], envir = .GlobalEnv)
  }
}
split_df_ASE(my_ASE_results)


# Read in Morris et al. BMD GWAS lead SNPs and lift them 
# over from Grch37 to Grch38 and expand to 400kb

morris_lead_snps <- fread("/scratch/aa9gj/wrk_backup/lncRNA_paper/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/test_wasp/WASP_analysis/test_rename/morris_gwas.txt")
ch <- import.chain("/scratch/aa9gj/wrk_backup/hg19ToHg38.over.chain")
colnames(morris_lead_snps)[3] <- "seqnames"
morris_lead_snps$seqnames <- gsub("23", "X", morris_lead_snps$seqnames)
cur <- makeGRangesFromDataFrame(morris_lead_snps, keep.extra.columns = TRUE, start.field = "BP", end.field = "BP")
seqlevelsStyle(cur) = "UCSC"  # necessary
cur38 = liftOver(cur, ch)
class(cur38)
cur38 = unlist(cur38)
genome(cur38) = "hg38"
cur38 = new("gwaswloc", cur38)
cur38 <- as.data.table(cur38)
cur38$start<- cur38$start-200000
cur38$end<- cur38$end+200000
morris_grange <- makeGRangesFromDataFrame(cur38, keep.extra.columns = TRUE, start.field = "start", end.field = "end")

# Create a GENCODE v34 database

database <- rtracklayer::import("/scratch/aa9gj/wrk_backup/genomes/gencode_v34_grch38/gencode.v34.primary_assembly.annotation.gtf")
database <- as.data.frame(database)
database_grange<-makeGRangesFromDataFrame(database, keep.extra.columns = TRUE, start.field = "start", end.field = "end")

# Read in stringtie 2nd output files from command line and separate to known genes

gtf_list <- read.delim("/scratch/aa9gj/wrk_backup/lncRNA_paper/human_ASE/alignment_RF/merged_all_bam/gtf_list", header = FALSE)
gtf_files <- list()
for (i in seq_along(gtf_list$V1)) {
  gtf_files[[i]] <- as.data.frame(rtracklayer::import(gtf_list$V1[i]))
  gtf_files[[i]] <- filter(gtf_files[[i]], TPM > 0.1)
  gtf_files[[i]] <- filter(gtf_files[[i]], ref_gene_name != is.na(gtf_files[[i]]$ref_gene_name))
  gtf_files[[i]] <- makeGRangesFromDataFrame(gtf_files[[i]], keep.extra.columns = TRUE, start.field = "start", end.field = "end")
}

# count the total number of known lncRNAs in all the samples 
gtf_list <- read.delim("/scratch/aa9gj/wrk_backup/lncRNA_paper/human_ASE/alignment_RF/merged_all_bam/gtf_list", header = FALSE)
gtf_files_total <- list()
setwd("/scratch/aa9gj/wrk_backup/lncRNA_paper/human_ASE/alignment_RF/merged_all_bam/")
for (i in seq_along(gtf_list$V1)) {
  gtf_files_total[[i]] <- as.data.frame(rtracklayer::import(gtf_list$V1[i]))
  gtf_files_total[[i]] <- filter(gtf_files_total[[i]], TPM > 0.1)
  gtf_files_total[[i]] <- filter(gtf_files_total[[i]], ref_gene_name != is.na(gtf_files_total[[i]]$ref_gene_name))
}
df_known <- do.call(rbind, gtf_files_total)
df_known_genes <- as.data.frame(unique(df_known$ref_gene_name))
df_known_genes <- inner_join(df_known_genes, database, by=c("unique(df_known$ref_gene_name)" = "gene_name"))
df_known_genes <- filter(df_known_genes, gene_type == "lncRNA")
df_known_genes <- filter(df_known_genes, type == "gene")
#write.table(df_known_genes, "/scratch/aa9gj/known_lnc_gtf", row.names = F, quote = F)
length(unique(df_known_genes$`unique(df_known$ref_gene_name)`))# 6612

df_known_genes <- df_known_genes %>% distinct(`unique(df_known$ref_gene_name)`, width, .keep_all = T)
mean(df_known_genes$width)
median(df_known_genes$width)
df_known_genes[order(df_known_genes$width, decreasing = T),]
df_known_genes <- inner_join(df_known_genes, df_known, by=c("unique(df_known$ref_gene_name)" = "ref_gene_name"))
mean(as.numeric(df_known_genes$TPM))
median(as.numeric(df_known_genes$TPM))
# Read in stringtie output for all genes since the one above was for only known genes

gtf_all_files <- list()
for (i in seq_along(gtf_list$V1)) {
  gtf_all_files[[i]] <- as.data.frame(rtracklayer::import(gtf_list$V1[i]))
  gtf_all_files[[i]] <- makeGRangesFromDataFrame(gtf_all_files[[i]], keep.extra.columns = TRUE, start.field = "start", end.field = "end")
}

# Split known granges objects to their respective files

split_df_gtf<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "known_grange"),list[[i]], envir = .GlobalEnv)
  }
}
split_df_gtf(gtf_files)

# Split the "all_gtf" granges objects to their respective files

split_df_all_gtf<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "all_grange"),list[[i]], envir = .GlobalEnv)
  }
}
split_df_all_gtf(gtf_all_files)

# Read in stringtie 2nd output files from command line and separate 
# to unknown lncRNAs with info from CPAT

cpat_list <- read.delim("/scratch/aa9gj/wrk_backup/lncRNA_paper/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/test_wasp/WASP_analysis/test_rename/cpat_list", header = FALSE)
setwd("/scratch/aa9gj/wrk_backup/lncRNA_paper/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/test_wasp/WASP_analysis/test_rename")
cpat_files <- list()
for (i in seq_along(cpat_list$V1)) {
  cpat_files[[i]] <- fread(cpat_list$V1[i])
  cpat_files[[i]] <- filter(cpat_files[[i]], coding_prob < 0.364)
  cpat_files[[i]] <- inner_join(cpat_files[[i]], as.data.frame(gtf_all_files[[i]]), by = c("V1" = "transcript_id"))
  cpat_files[[i]] <- filter(cpat_files[[i]], type == "transcript")
  cpat_files[[i]] <- makeGRangesFromDataFrame(cpat_files[[i]], keep.extra.columns = TRUE, start.field = "start", end.field = "end")
}
# How many novel lncRNAs (make sure the above command doesn't have makegranges to caluclate this)
cpat_list <- read.delim("/scratch/aa9gj/wrk_backup/lncRNA_paper/human_ASE/alignment_RF/merged_all_bam/WASP_corrected_ASE/test_wasp/WASP_analysis/test_rename/cpat_list", header = FALSE)
cpat_files_total <- list()
for (i in seq_along(cpat_list$V1)) {
  cpat_files_total[[i]] <- fread(cpat_list$V1[i])
  cpat_files_total[[i]] <- filter(cpat_files_total[[i]], coding_prob < 0.364)
  cpat_files_total[[i]] <- inner_join(cpat_files_total[[i]], as.data.frame(gtf_all_files[[i]]), by = c("V1" = "transcript_id"))
  cpat_files_total[[i]] <- filter(cpat_files_total[[i]], type == "transcript")
  cpat_files_total[[i]]$sample <- rep(cpat_list$V1[i], length(cpat_files_total[[i]]$V1))
}

df<-do.call(rbind, cpat_files_total)
length(unique(df$V1)) #2440
df <- df %>% distinct(V1, width, .keep_all = TRUE)
#write.table(df, "/scratch/aa9gj/unknown_lncRNAs_gtf", row.names = F, quote = F)
mean(as.numeric(df$TPM))#123 tpm
median(as.numeric(df$TPM))#1.91 tpm 

mean(df$width)
median(df$width)
min(df$width)
max(df$width)
# Split unknown lncRNAs granges objects to their respective files

split_df_unknown_granges<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "unknown_lnc_grange"),list[[i]], envir = .GlobalEnv)
  }
}
split_df_unknown_granges(cpat_files)

# Perform overlap analysis for GWAS (known and unknown lncRNAs in the 17 samples that fall within 400kb of GWAS loci)
# unknown lncRNAs that are within +- 200kb from 1103 bmd gwas lead loci using unknown grange objects

morris_overlap_unknown <- function(x) {
  hits_unknown_lnc  <- findOverlaps(query = x, subject = morris_grange)
  olap_unknown_lnc  <- pintersect(x[queryHits(hits_unknown_lnc)],
                                  morris_grange[subjectHits(hits_unknown_lnc)])
  return(olap_unknown_lnc)
}

morris_overlap_unknown_results <- list()
for (i in seq_along(counts_files$V1)) {
  morris_overlap_unknown_results[[i]] <- morris_overlap_unknown(cpat_files[[i]])
}

# Split files
split_df_unknown_overlap<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "morris_unknown_overlap"), list[[i]], envir = .GlobalEnv)
  }
}
split_df_unknown_overlap(morris_overlap_unknown_results)

# Known lncRNAs that are within +- 200kb from 1103 bmd gwas lead loci using known grange objects

morris_overlap_known <- function(x) {
  hits_known_lnc  <- findOverlaps(query = x, subject = morris_grange)
  olap_known_lnc  <- pintersect(x[queryHits(hits_known_lnc)],
                                morris_grange[subjectHits(hits_known_lnc)])
  olap_known_lnc <- inner_join(as.data.frame(olap_known_lnc), as.data.frame(database_grange), by = c("ref_gene_name"="gene_name"))
  olap_known_lnc <- filter(olap_known_lnc, gene_type == "lncRNA")
  olap_known_lnc <- filter(olap_known_lnc, type.y == "gene")
  return(olap_known_lnc)
}

morris_overlap_known_results <- list()
for (i in seq_along(counts_files$V1)) {
  morris_overlap_known_results[[i]] <- morris_overlap_known(gtf_files[[i]])
  morris_overlap_known_results[[i]] <- makeGRangesFromDataFrame(morris_overlap_known_results[[i]], keep.extra.columns = TRUE, start.field = "start.x", end.field = "end.x", seqnames.field = "seqnames.x")
}

# Split files
split_df_known_overlap<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "morris_known_overlap") ,list[[i]], envir = .GlobalEnv)
  }
}
split_df_known_overlap(morris_overlap_known_results)

# Now we need unknown and known lncRNAs that are implicated by ASE only

ASE_only_unknown <- function(x,y) {
  hits_ase <- findOverlaps(query = x[[i]], subject = y[[i]])
  olap_ase  <- pintersect(x[[i]][queryHits(hits_ase)],
                          y[[i]][subjectHits(hits_ase)])
  return(olap_ase)
}

ASE_only_unknown_results <- list()
for (i in seq_along(counts_files$V1)) {
  ASE_only_unknown_results[[i]] <- ASE_only_unknown(cpat_files,my_ASE_results)
}

#Split it
split_df_ASE_only_unknown<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "ASE_unknown") ,list[[i]], envir = .GlobalEnv)
  }
}
split_df_ASE_only_unknown(ASE_only_unknown_results)


ASE_only_known <- function(x,y) {
  hits_ase <- findOverlaps(query = x[[i]], subject = y[[i]])
  olap_ase  <- pintersect(x[[i]][queryHits(hits_ase)],
                          y[[i]][subjectHits(hits_ase)])
  olap_ase <- inner_join(as.data.frame(olap_ase), as.data.frame(database_grange), by = c("ref_gene_name"="gene_name"))
  olap_ase <- filter(olap_ase, gene_type == "lncRNA")
  olap_ase <- filter(olap_ase, type.y == "gene")
  return(olap_ase)
}

ASE_only_known_results <- list()
for (i in seq_along(counts_files$V1)) {
  ASE_only_known_results[[i]] <- ASE_only_known(gtf_files,my_ASE_results)
  ASE_only_known_results[[i]] <- makeGRangesFromDataFrame(ASE_only_known_results[[i]], keep.extra.columns = TRUE, start.field = "start.x", end.field = "end.x", seqnames.field = "seqnames.x")
}

##split it
split_df_ASE_only_known<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "ASE_known") ,list[[i]], envir = .GlobalEnv)
  }
}
split_df_ASE_only_known(ASE_only_known_results)

# overlap between GWAS and ASE for unknown lncRNAs

GWAS_ASE_unknown <- function(x,y) {
  hits_ase <- findOverlaps(query = x[[i]], subject = y[[i]])
  olap_ase  <- pintersect(x[[i]][queryHits(hits_ase)],
                          y[[i]][subjectHits(hits_ase)])
  return(olap_ase)
}

GWAS_ASE_unknown_results <- list()
for (i in seq_along(counts_files$V1)) {
  GWAS_ASE_unknown_results[[i]] <- GWAS_ASE_unknown(morris_overlap_unknown_results,my_ASE_results)
}
##split it
split_df_GWAS_ASE_unknown<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "GWAS_ASE_unknown") ,list[[i]], envir = .GlobalEnv)
  }
}
split_df_GWAS_ASE_unknown(GWAS_ASE_unknown_results)

###now for known lncRNAs

GWAS_ASE_known <- function(x,y) {
  hits_ase <- findOverlaps(query = x[[i]], subject = y[[i]])
  olap_ase <- pintersect(x[[i]][queryHits(hits_ase)],
                         y[[i]][subjectHits(hits_ase)])
  return(olap_ase)
}

GWAS_ASE_known_results <- list()
for (i in seq_along(counts_files$V1)) {
  GWAS_ASE_known_results[[i]] <- GWAS_ASE_known(morris_overlap_known_results,my_ASE_results)
}
##split it
split_df_GWAS_ASE_known<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "GWAS_ASE_known") ,list[[i]], envir = .GlobalEnv)
  }
}
split_df_GWAS_ASE_known(GWAS_ASE_known_results)

GWAS_ASE_unknown <- function(x,y) {
  hits_ase <- findOverlaps(query = x[[i]], subject = y[[i]])
  olap_ase  <- pintersect(x[[i]][queryHits(hits_ase)],
                          y[[i]][subjectHits(hits_ase)])
  return(olap_ase)
}

GWAS_ASE_unknown_results <- list()
for (i in seq_along(counts_files$V1)) {
  GWAS_ASE_unknown_results[[i]] <- GWAS_ASE_unknown(morris_overlap_unknown_results,my_ASE_results)
}
##split it
split_df_GWAS_ASE_unknown<-function(list){
  for (i in 1:length(list)){
    assign(paste0(file_names$file_names[i], "GWAS_ASE_unknown") ,list[[i]], envir = .GlobalEnv)
  }
}
split_df_GWAS_ASE_unknown(GWAS_ASE_unknown_results)

# Now let's make dataframes for information in all data

ASE_only_all_known_lncRNAs <- do.call("c", ASE_only_known_results)
ASE_only_all_unknown_lncRNAs <- do.call("c", ASE_only_unknown_results)
GWAS_ASE_all_known_lncRNAs <- do.call("c", GWAS_ASE_known_results)   
GWAS_ASE_all_unknown_lncRNAs <- do.call("c", GWAS_ASE_unknown_results)
GWAS_only_all_known_lncRNAs <- do.call("c", morris_overlap_known_results) 
GWAS_only_all_unknown_lncRNAs <- do.call("c", morris_overlap_unknown_results)

length(unique(ASE_only_all_known_lncRNAs$ref_gene_name))#47
length(unique(ASE_only_all_unknown_lncRNAs$gene_id))#127
length(unique(GWAS_ASE_all_known_lncRNAs$ref_gene_name))#8
as.data.frame(unique(GWAS_ASE_all_known_lncRNAs$ref_gene_name))
length(unique(GWAS_ASE_all_unknown_lncRNAs$gene_id))#19
as.data.frame(unique(GWAS_ASE_all_unknown_lncRNAs$gene_id))
length(unique(GWAS_only_all_known_lncRNAs$ref_gene_name))#1215
length(unique(GWAS_only_all_unknown_lncRNAs$gene_id))#281