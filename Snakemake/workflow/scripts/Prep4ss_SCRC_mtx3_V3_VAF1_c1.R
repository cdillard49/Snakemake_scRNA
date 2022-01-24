library("tidyverse")
library("data.table")

args <- commandArgs(trailingOnly = TRUE)

#setwd("~/Documents/Papers/neuroblastoma/")
#read the ori, summarize and filter it
sample_id <- args[1]
pattern <-'.vaf-m3.matrix.tsv'
ori_name <- paste0(sample_id, pattern)
ori = fread(ori_name)
ori$SUM = rowSums(ori[,2:ncol(ori)], na.rm = T)
ori$nonZeroCOUNT = rowSums(ori[,2:(ncol(ori)-1)]!=0, na.rm = T)

ori_vaf2 = ori %>% filter(SUM>=1)
ori_vaf2_c2 = ori_vaf2 %>% filter(nonZeroCOUNT>0)
ori_vaf2_c2_summary = ori_vaf2_c2 %>% select(SNV,SUM,nonZeroCOUNT)
to_ss = ori_vaf2_c2 %>% select(1) %>% distinct() %>% separate(1, c("chr","pos","ref","alt"), sep = "[:_>]")

#write outputs
fwrite(to_ss, paste0(sample_id, "_2ss_R3.txt"), sep = '\t')
fwrite(ori_vaf2_c2, paste0(sample_id, "_VAF1_c1_mtx_R3.txt"), sep = '\t')
fwrite(ori_vaf2_c2_summary, paste0(sample_id, "_VAF1_c1_sum_R3.txt"), sep = '\t')
