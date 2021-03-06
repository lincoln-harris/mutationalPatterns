# load libraries
library(MutationalPatterns)
library(BSgenome)

# load hg38
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# list some sample vcfs
vcf_files <- list.files(path = './laud_samples', pattern = ".vcf", full.names = TRUE)

# get cell names
cell_names <- stringr::str_remove_all(vcf_files, ".vcf")
cell_names <- stringr::str_remove_all(cell_names, "./laud_samples/")

# load vcfs, as Grange objs
vcfs <- read_vcfs_as_granges(vcf_files, cell_names, ref_genome)

# get type of each mut
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

# plot, and save
pdf("myplot2.pdf", width=4, height=4)
p1 <- plot_spectrum(type_occurrences)
dev.off()
