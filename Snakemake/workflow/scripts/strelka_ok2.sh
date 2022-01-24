#!/usr/bin/bash
# How to run the script:


# Germline: sh strelka.sh Germline list /path/to/fastaref dir_name_prefix
# list is a list of paths to the bam file
fa=/data/lab/Homo_sapiens/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

ml strelka
ml bcftools
if [[ "$1" == "Germline" ]]
then
	while read l
	do
	{
		configureStrelkaGermlineWorkflow.py --bam ${l}".bam" --referenceFasta ${fa} --rna --runDir ${3}"_"${l}
		python ${3}"_"${l}/runWorkflow.py -m local &> log_${l}
		mv ${3}"_"${l}"/results/variants/variants.vcf.gz" ${3}"_"${l}"/results/variants/"${l}"_strelka.vcf.gz"
		mv ${3}"_"${l}"/results/variants/"${l}"_strelka.vcf.gz" .
		bcftools index ${l}"_strelka.vcf.gz"
                bcftools filter \
                        -i 'TYPE="snp" && FILTER="PASS"' \
                        -r1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y ${l}"_strelka.vcf.gz" -Ov -o ${l}"_strelka_filt.vcf"
	}
	done<${2}
else
	echo "Strelka does not have this option"
fi
