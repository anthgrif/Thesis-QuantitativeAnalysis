# function to convert relative abundance to quantitative abundance
rel_to_quan <- function(physeq, counts) {
  for (i in 1:nsamples(physeq)) {
    otu_table(physeq)[,i] = get_taxa(physeq, sample_names(physeq)[i]) * counts[i] /100
  }
  return(otu_table(physeq))
}