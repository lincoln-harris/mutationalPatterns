# install package (LEAVE COMMENTED OUT)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MutationalPatterns")

# load libraries
library(MutationalPatterns)
library(BSgenome)

# bullshit to install hg38 (LEAVE COMMENTED OUT)
#if (!require("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# load hg38
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# list some sample vcfs
vcf_files <- list.files(path = './laud_samples_sub', pattern = ".vcf", full.names = TRUE)

# get cell names
cell_names <- stringr::str_remove_all(vcf_files, ".vcf")
cell_names <- stringr::str_remove_all(cell_names, "./laud_samples_sub/")

# load vcfs, as Grange objs
vcfs <- read_vcfs_as_granges(vcf_files, cell_names, ref_genome, check_alleles = TRUE)
print(' ')
print('finished loading vcfs')
print(' ')

# get type of each mut
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

# plot, and save
pdf("mut_type_laud_gFilter.pdf", width=4, height=4)
plot_spectrum(type_occurrences)
dev.off()

