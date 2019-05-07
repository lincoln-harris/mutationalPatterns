my_read_vcfs_as_granges <- function (vcf_files, sample_names, genome, group = "auto+sex", 
          check_alleles = TRUE) 
{
  library(MutationalPatterns)
  library(BSgenome)
  
  if (length(vcf_files) != length(sample_names)) 
    stop("Please provide the same number of sample names as VCF files")
  ref_genome <- base::get(genome)
  ref_organism <- GenomeInfoDb::organism(ref_genome)
  ref_style <- seqlevelsStyle(ref_genome)
  genome_name <- genome(ref_genome)[[1]]
  if (!(class(ref_genome) == "BSgenome")) 
    stop("Please provide the name of a BSgenome object.")
  num_cores = detectCores()
  if (!(.Platform$OS.type == "windows" || is.na(num_cores))) 
    num_cores <- detectCores()
  else num_cores = 1
  original_warn_state = getOption("warn")
  options(warn = -1)
  warnings <- NULL
  if (!check_alleles) {
    warnings <- c(warnings, paste("check_alleles is set to FALSE.  Make sure your", 
                                  "input VCF does not contain any positions with", 
                                  "insertions, deletions or multiple alternative", 
                                  "alleles, as these positions cannot be analysed", 
                                  "with MutationalPatterns and cause obscure", "errors."))
  }
  vcf_list <- mclapply(seq_along(vcf_files), function(index) {
    file <- vcf_files[index]
    print(file)
    vcf <- rowRanges(readVcf(file, genome_name))
    seqlevelsStyle(vcf) <- ref_style[1]
    groups <- c()
    if (group != "none") {
      if (group == "auto+sex") {
        groups <- c(extractSeqlevelsByGroup(species = ref_organism, 
                                            style = ref_style, group = "auto"), extractSeqlevelsByGroup(species = ref_organism, 
                                                                                                        style = ref_style, group = "sex"))
        groups_names <- names(groups)
        if (!is.null(groups_names)) {
          unique_names <- unique(groups_names)
          groups <- llply(unique_names, function(x) groups[groups_names == 
                                                             x])
          groups <- llply(groups, unlist, recursive = FALSE)
          groups <- unique(as.vector(groups[[1]]))
        }
      }
      else {
        groups <- extractSeqlevelsByGroup(species = ref_organism, 
                                          style = ref_style, group = group)
        groups <- unique(as.vector(t(groups)))
      }
      groups <- intersect(groups, seqlevels(vcf))
      vcf <- keepSeqlevels(vcf, groups, pruning.mode = "tidy")
    }
    if (check_alleles) {
      rem <- which(all(!(!is.na(match(vcf$ALT, DNA_BASES)) & 
                           !is.na(match(vcf$REF, DNA_BASES)) & (lengths(vcf$ALT) == 
                                                                  1))))
      if (length(rem) > 0) {
        vcf = vcf[-rem]
        warnings <- c(warnings, paste(length(rem), "position(s) with indels and/or multiple", 
                                      "alternative alleles are excluded from", paste(sample_names[[index]], 
                                                                                     ".", sep = "")))
      }
    }
    return(list(vcf, warnings))
  }, mc.cores = num_cores)
  options(warn = original_warn_state)
  vcf_list <- lapply(vcf_list, function(item) {
    if (class(item) == "try-error") 
      stop(item)
    if (!is.null(item[[2]])) 
      for (i in item[[2]]) warning(i)
    return(item[[1]])
  })
  vcf_list <- GRangesList(vcf_list)
  names(vcf_list) <- sample_names
  return(vcf_list)
}