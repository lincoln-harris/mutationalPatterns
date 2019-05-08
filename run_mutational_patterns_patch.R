# load libraries
setwd('/home/rstudio/mutationalPatterns/')

library(devtools)
devtools::install_github("UMCUGenetics/MutationalPatterns", ref = "lincoln-harris", force = TRUE)

library(MutationalPatterns)
library(BSgenome)

# load hg38
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# list some sample vcfs
vcf_files <- list.files(path = './samples_filter_all', pattern = ".vcf", full.names = TRUE)

# get cell names
cell_names <- stringr::str_remove_all(vcf_files, ".vcf")
cell_names <- stringr::str_remove_all(cell_names, "./samples_filter_all/")

# load vcfs, as Grange objs
vcfs <- read_vcfs_as_granges(vcf_files, cell_names, ref_genome, check_alleles = TRUE)
print(' ')
print('finished loading vcfs')
print(' ')

# get type of each mut
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

# lets just save this mut_type_occurances dataframe
write.csv(type_occurrences, 'type_occurences_df.csv')

