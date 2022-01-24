rule sort:
    input:
#        barcode_pat_gatk= make_barcodes_gatk,
#        barcode_pat_strelka= make_barcodes_strelka
#        gatk_filt= expand("fastqs/gatk_output/{sample}_wasp_Aligned.sortedByCoord_vW_filt.{barcode}_gatk_filt.vcf", sample=config["Biosample"], barcode=barcodes), 
#        strelka_filt= expand("fastqs/strelka/{sample}_wasp_Aligned.sortedByCoord_vW_filt.{barcode}_strelka_filt.vcf", sample=config["Biosample"], barcode=barcodes)
        strelka_filt_dir= directory("fastqs/strelka"), 
        gatk_filt_dir= directory("fastqs/gatk_output")

    output:
        sorted_txt= expand("fastqs/ind_bams/{sample}_wasp_Aligned.sortedByCoord_vW_filt.{barcode}_{ext_sorted}", sample=config["Biosample"], barcode=barcodes, ext_sorted=["strelka_sorted.txt", "gatk_sorted.txt"])
#        sorted_strelka= "fastqs/ind_bams/{sample}_sample/{barcode}_strelka_sorted.txt", 
#        sorted_gatk= "fastqs/ind_bams/{sample}_sample/{barcode}_gatk_sorted.txt"

    shell:
        """for i in {config[Biosample]}; do mkdir -p {config[Path2output]}; done; mv {config[Path2gatk]}*_gatk_filt.vcf {config[Path2samples]}; mv {config[Path2strelka]}*_strelka_filt.vcf {config[Path2samples]}; cd {config[Path2samples]}; for i in {config[Path2samples]}*_filt.vcf; do grep -v '#' ${{i}} | awk '{{print $1":"$2"_"$4">"$5 }}' | sort > ${{i%%_filt.vcf}}_sorted.txt; done; for i in{config[Biosample]}; do mv ${{i}}_wasp_Aligned.sortedByCoord_vW_filt. ;done"""

rule common rows:
    input:
        sorted_txt= expand("fastqs/ind_bams/{sample}_wasp_Aligned.sortedByCoord_vW_filt.{barcode}_{ext_sorted}", sample=config["Biosample"], barcode=barcodes, ext_sorted=["strelka_sorted.txt", "gatk_sorted.txt"])
#        sorted_strelka= "fastqs/ind_bams/{sample}_sample/{barcode}_strelka_sorted.txt", 
#        sorted_gatk= "fastqs/ind_bams/{sample}_sample/{barcode}_gatk_sorted.txt"
    output:
        comm_file= expand("fastqs/ind_bams/{sample}_wasp_Aligned.sortedByCoord_vW_filt.{barcode}_comm.txt",sample=config["Biosample"], barcode=barcodes)
#    log:
#        "2> log/common_row.log"
    shell:
        "cd {config[Path2samples]}; for i in {config[Path2samples]}*_gatk_sorted.txt; do comm -12 ${{i}} ${{i%%_gatk_sorted.txt}}_strelka_sorted.txt > ${{i%%_gatk_sorted.txt}}_comm.txt; done"

rule combine:
    input:
        comm_file= expand("fastqs/ind_bams/{sample}_wasp_Aligned.sortedByCoord_vW_filt.{barcode}_comm.txt", sample=config["Biosample"], barcode=barcodes)
    output:
        comb_file= expand("fastqs/ind_bams/{sample}_comb.txt", sample=config["Biosample"])
#    log:
#        "2> log/combine.log"
    shell:
        """cd {config[Path2samples]}; for i in {config[Path2samples]}*_comm.txt; do while read line ; do echo ${{line}} ${{i}} >> ${{i%%_*}}_comb.txt; done < ${{i}}; done"""


rule ss_prep:
    input:
        comb_file= expand("fastqs/ind_bams/{sample}_comb.txt", sample=config["Biosample"])
    output:
        to_ss_file= expand("fastqs/ind_bams/{sample}_2ss.txt", sample=config["Biosample"])
#    log:
#        "2> log/ss_prep.log"
    shell:
        """cd {config[Path2samples]}; for i in {config[Path2samples]}*_comb.txt; do awk '{{print $1}}' ${{i}} | awk -F'[:_>]' '{{print $1 "\t" $2 "\t" $3 "\t" $4}}' | awk '!seen[$0]++' > ${{i%%_comb.txt}}_2ss.txt; rm {config[Path2samples]}*_sorted.txt; rm {config[Path2samples]}*_comm.txt; done"""

  ##download the 2ss file and upload to seatle seq
  
##############################################