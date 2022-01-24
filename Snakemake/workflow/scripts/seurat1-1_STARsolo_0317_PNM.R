# Seurat1-1.R
# LAST UPDATED ON DEC.22 2020
# This script is to generate the Seurat original datset to visualization the feature distribution for downstream analysis

# load package

print('loading required packages: data.table, tidyverse, Seurat, biomaRt')
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(optparse))
# manages command line arguments
option_list <- list(
  make_option(c("-s", "--samplelist"), action = "store", type ="character", 
              default = NULL, help= "File (without header) that lists all samples (only sample names)")
  #make_option(c("-p", "--pattern"), action="store", type = "character", default="_wide_counts.tsv",
              #help="provide the common suffix to the sample files containg GE values. DEFAULT: _wide_counts.tsv")
)
opt <- parse_args(OptionParser(option_list=option_list, description = "-s option is necessary!!!!", usage = "usage: Rscript seurat1-1.R -s <sample_list>"))
opt$samplelist <- "list"

if (is.null(opt$samplelist)) {
  stop("Please make sure you have provided the sample list")
}

if (packageVersion("Seurat") < "3.0.0") {
  stop(paste0("You have Seurat version", packageVersion("Seurat"), "installed. Please make sure you have Seurat version > 3.0.0"))
}


# sample list
sample_list <<- fread(opt$samplelist, header = F)
sample_id <<- sample_list$V1
# check whether sample files exists or not
# for (i in 1:nrow(sample_list)){
#   if (file.exists(paste0(sample_id[i], opt$pattern) == F))
#     print(paste0('provided ', sample_id[i],' DOES NOT exist!'))
#   else { print(paste0('Processing: ', sample_id[i]))}
# }

for (i in 1:nrow(sample_list)){
  #import featureCount gene_expression matrix
    #Gene_matrix <- fread(paste0(sample_id[i], '_wide_counts.tsv'))
    #Gene_matrix <- spread(Gene_matrix, cell, count)
    # ensemble id change to gene_name
    data <- Read10X(data.dir = sample_id[i])
    Gene_matrix = as.data.frame(data)
    Gene_matrix$gene <- row.names(data)
    Gene_matrix <- Gene_matrix[,c(ncol(Gene_matrix),1:(ncol(Gene_matrix)-1))]
    # # row.names(Gene_matrix) <- NULL
    # library(biomaRt)
    # ##listEnsemblArchives()
    # # use human sapines ensembl verson100 data to annotate (should change by different dataset, keep the lastest version)
    # ensemble <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "ensembl", version =100)
    # query <- getBM(attributes =
    #                  c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position"),
    #                filters = "external_gene_name", values = Gene_matrix$gene, mart = ensemble)
    # colnames(query)[2] <- "gene"
    # print("a")
    # Gene_matrix_name <- inner_join(query,Gene_matrix, by = "gene")
    # Gene_matrix_name <- Gene_matrix_name[,-c(1,3,4,5)]
    # # deduplicate file
    # index <- duplicated(Gene_matrix_name$external_gene_name)
    # Gene_matrix_name <- Gene_matrix_name[!index,]
    # rownames(Gene_matrix_name) <- Gene_matrix_name$external_gene_name
    # Gene_matrix_name <- Gene_matrix_name[,-1]
    # print("b")
    #Create Seurat object with specific paramteres (command-line arguments)
    Gene_Seurat <- CreateSeuratObject(counts = Gene_matrix, project = gsub("*_wasp_Solo.out/Gene/filtered","",sample_id[i]), min.cells = 300, min.features =2000)
    #calculate mito QC metrics
    Gene_Seurat[["percent.mt"]] <- PercentageFeatureSet(Gene_Seurat, pattern = "^MT-")
    print(head(Gene_Seurat[["percent.mt"]]))
    
    # visualization for the feature distribution
    v <- VlnPlot(object = Gene_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    png(paste0(gsub("*_wasp_Solo.out/Gene/filtered","",sample_id[i]),'_','feature_distribution_vlnplot.png'), width = 1000, height = 400)
    print(v)
    dev.off()
    
    # density plot showing proportion of transcripts expressing genes and mito genes per cell vs UMIs per cell
    plot1 <- FeatureScatter(object = Gene_Seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_point()+ scale_fill_viridis_c() + stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") + geom_hex(bins = 70)  + theme_bw()
    plot2 <- FeatureScatter(object = Gene_Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ geom_point() + stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +scale_fill_viridis_c() + geom_hex(bins = 70)  + theme_bw()
    png(paste0(gsub("*_wasp_Solo.out/Gene/filtered","",sample_id[i]), '_feature_distribution.png'), width = 1000, height = 400)
    print(plot1 + plot2)
    dev.off()
    
    # saving original Seurat dataset
    saveRDS(Gene_Seurat,paste0(gsub("*_wasp_Solo.out/Gene/filtered","",sample_id[i]),'_beforefilter.rds'))
}

print('DONE')
