#!/bin/bash
#For RNA-Seq, preprocess the RNA-Seq data 
#splits reads into exon segments (getting rid of Ns but maintaining grouping information) 
#hard-clip any sequences overhanging into the intronic regions
ml gatk
ml samtools
ml bcftools

list=${1}
while read line
do
{
gatk AddOrReplaceReadGroups \
       I=$line".bam" \
       O=$line"_rg.bam" \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

samtools index -@ 10 $line"_rg.bam"

gatk SplitNCigarReads \
	-R /data/lab/Homo_sapiens/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	-I $line"_rg.bam" \
	-O $line"_split_rg.bam" \
	-RF MappingQualityReadFilter

samtools index -@ 10 $line"_split_rg.bam"

gatk --java-options "-Xmx4g" HaplotypeCaller \
	-R /data/lab/Homo_sapiens/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	-I $line"_split_rg.bam" \
	-O $line"_gatk.vcf"

bcftools view $line"_gatk.vcf" -Oz -o $line"_gatk.vcf.gz"

bcftools index $line"_gatk.vcf.gz"

bcftools filter \
	-i 'TYPE="snp" && MIN(INFO/DP)>=3 && QUAL>=100 && MQ==60' \
	-r1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $line"_gatk.vcf.gz" -Ov -o $line"_gatk_filt.vcf"

rm $line"_rg.bam"
rm $line"_rg.bam.bai"
rm $line"_split_rg.bam"
rm $line"_split_rg.bai"
rm $line"_split_rg.bam.bai"
rm $line"_gatk.vcf.idx"
rm $line"_gatk.vcf.gz"
rm $line"_gatk.vcf.gz.csi"
}
done < ${list}
echo "DONE"
