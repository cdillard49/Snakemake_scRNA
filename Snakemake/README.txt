Steps for Snakemake:
--------------------

1. activate conda environment where the snakemake program is downloaded
    
2. add SRR to the config.yaml file under `SRR`

3. from the Snakemake dir, run `snakemake --cores [N] --use-conda --conda-frontend conda --verbose -U STARsolo`
    #This will generate the fastqs
    #The -U [RULENAME] means run snakemake until the specified rule. 
    #This could be a different name than in this readme so check the workflow/rules/preprossesing.smk for exact name.
    
  3a. If needed update the trimmomatic list with the SRRs, and run trimmomatic script to trim the R1 to the approriate barcode length.
  
  3b. Move trimmed fastq back into the fastq dir with the same name (important for the script!)
  
4. from the Snakemake dir, run `snakemake --cores [N] --use-conda --conda-frontend conda --verbose -U samtoolsprelude`
    #This will generate the bams
    
  4a. cat the bam files together based on the biosample (i.e `cat SRRfoo.bam SRRbar.bam SRRspam.bam > SRRfoobarspam.bam`) 
  4b. 2a. add the SRR's biosample into the config.yaml file under `Sample`

5. from Snakemake dir, run `snakemake --cores [N] --use-conda --conda-frontend conda --verbose -U GATK_ind`
  
