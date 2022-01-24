library("tidyverse")
library("data.table")

args <- commandArgs(trailingOnly = TRUE)

#setwd("~/Documents/Papers/neuroblastoma/")
#read the ori, summarize and filter it
sample_id <- args[1]
pattern <-'.vaf-m5.matrix.tsv'
ori_name <- paste0(sample_id, pattern)
ori = fread(ori_name)
ori$SUM = rowSums(ori[,2:ncol(ori)], na.rm = T)
ori$nonZeroCOUNT = rowSums(ori[,2:(ncol(ori)-1)]!=0, na.rm = T)

ori_vaf_non0 = ori %>% filter(SUM>0)
ori_vaf_above1 = ori %>% filter(SUM>1)
ori_summary = ori_vaf_non0 %>% select(SNV,SUM,nonZeroCOUNT)
ori_summary_above1 = ori_vaf_above1 %>% select(SNV,SUM,nonZeroCOUNT)
to_ss = ori_vaf_non0 %>% select(1) %>% distinct() %>% separate(1, c("chr","pos","ref","alt"), sep = "[:_>]")

#write ss-ready inputs
fwrite(to_ss, paste0(sample_id, "_2ss.txt"), sep = '\t')
fwrite(ori_vaf_non0, paste0(sample_id, "_mtx_non0.txt"), sep = '\t')
#fwrite(ori_vaf_above1,paste0(sample_id, "_mtx_above1.txt"), sep = '\t')
fwrite(ori_summary, paste0(sample_id, "_non0_summary.txt"), sep = '\t')
