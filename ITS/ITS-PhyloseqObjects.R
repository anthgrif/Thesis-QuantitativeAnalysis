library("phyloseq")
packageVersion("phyloseq")

library("biomformat")
packageVersion("biomformat")

biom_data <- import_biom(BIOMfilename = "feature-table-with-taxa.biom", 
                         treefilename = "midpoint-rooted-tree.nwk")
mapping_file <- import_qiime_sample_data(mapfilename = "ITS-metadata-withcounts.txt")

physeq <- merge_phyloseq(biom_data, mapping_file)
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")

samples_to_keep = c("T1R1", "T1R10", "T1R4", "T1R7", "T2R1", "T2R10", "T2R4", "T2R5", "T2R7","T2R9", "T310", "T3R1", "T3R4", "T3R5", "T3R7", "T3R9", "T410", "T4R1", "T4R4", "T4R5", "T4R7", "T4R9", "T5R1", "T5R10", "T5R4", "T5R5", "T5R7", "T5R9")
physeq.28s <- subset_samples(physeq, SampleID %in% samples_to_keep)

# exclude R10 samples because they don't have counts data
physeq.23s <- subset_samples(physeq.28s, SampleID != "T1R10" & SampleID != "T2R10" & SampleID != "T310" & SampleID != "T410" & SampleID != "T5R10")

# function to convert relative abundance to quantitative abundance
rel_to_quan <- function(physeq, counts) {
  for (i in 1:nsamples(physeq)) {
    otu_table(physeq)[,i] = get_taxa(physeq, sample_names(physeq)[i]) * counts[i] /100
  }
  return(otu_table(physeq))
}
