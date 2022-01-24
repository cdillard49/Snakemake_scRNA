#!/bin/bash
#For RNA-Seq, preprocess the RNA-Seq data 
#splits reads into exon segments (getting rid of Ns but maintaining grouping information) 
#hard-clip any sequences overhanging into the intronic regions



"filename={input.ind_bam_list} while read line; do {gatk AddOrReplaceReadGroups \ I=${{i}}'.bam' \ O=${{i}}'_rg.bam' \ RGLB=lib1 \ RGPL=illumina \ RGPU=unit1 \ RGSM=20 ; samtools index -@ 10 ${{i}}'_rg.bam' ; gatk SplitNCigarReads \ -R {input.FA} \ -I ${{i}}'_rg.bam' \ -O ${{i}}'_split_rg.bam' \ -RF MappingQualityReadFilter ; samtools index -@ 10 ${{i}}'_split_rg.bam' ; gatk --java-options '-Xmx4g' HaplotypeCaller \ -R {input.FA} \ -I ${{i}}'_split_rg.bam' \ -O ${{i}}'_gatk.vcf' ; bcftools view ${{i}}'_gatk.vcf' -Oz -o ${{i}}'_gatk.vcf.gz' ; bcftools index ${{i}}'_gatk.vcf.gz' ; bcftools filter \ -i 'TYPE='snp' && MIN(INFO/DP)>=3 && QUAL>=100 && MQ==60' \ -r1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y ${{i}}'_gatk.vcf.gz' -Ov -o ${{i}}'_gatk_filt.vcf' ; 
rm ${{i}}'_rg.bam'; 
rm ${{i}}'_rg.bam.bai';
rm ${{i}}'_split_rg.bam';
rm ${{i}}'_split_rg.bai';
rm ${{i}}'_split_rg.bam.bai';
rm ${{i}}'_gatk.vcf.idx';
rm ${{i}}'_gatk.vcf.gz';
rm ${{i}}'_gatk.vcf.gz.csi';"

