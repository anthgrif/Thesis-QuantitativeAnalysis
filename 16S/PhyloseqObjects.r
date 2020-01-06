library("phyloseq")
packageVersion("phyloseq")

library("biomformat")
packageVersion("biomformat")

biom_data <- import_biom(BIOMfilename = "table-with-taxa.biom", 
                         treefilename = "tree.nwk")
mapping_file <- import_qiime_sample_data(mapfilename = "16s-metadata-with-counts.tsv")
physeq.a <- merge_phyloseq(biom_data, mapping_file)

colnames(tax_table(physeq.a))= c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")

whole.samples <- c("T1R1","T1R4","T1R5","T1R7","T1R9","T1R10","T2R1","T2R4","T2R5","T2R7","T2R9","T2R10","T3R1","T3R4","T3R5","T3R7","T3R9","T3R10","T4R1","T4R4","T4R5","T4R7","T4R9","T4R10","T5R1","T5R4","T5R5","T5R7","T5R9","T5R10")
live.samples <- c("T2R1L","T2R4L","T2R5L","T2R7L","T2R9L","T3R1L","T3R4L","T3R5L","T3R7L","T3R9L","T3R10L","T4R1L","T4R4L","T4R5L","T4R7L","T4R9L","T4R10L","T5R1L","T5R4L","T5R5L","T5R7L","T5R9L","T5R10L")
dead.samples <- c("T2R1D","T2R4D","T2R5D","T2R7D","T2R9D","T2R10D","T3R1D","T3R4D","T3R5D","T3R7D","T3R9D","T3R10D","T4R1D","T4R4D","T4R5D","T4R7D","T4R9D","T4R10D","T5R1D","T5R4D","T5R5D","T5R7D","T5R9D","T5R10D")
live.dead.samples <- c("T2R1L","T2R4L","T2R5L","T2R7L","T2R9L","T3R1L","T3R4L","T3R5L","T3R7L","T3R9L","T3R10L","T4R1L","T4R4L","T4R5L","T4R7L","T4R9L","T4R10L","T5R1L","T5R4L","T5R5L","T5R7L","T5R9L","T5R10L","T2R1D","T2R4D","T2R5D","T2R7D","T2R9D","T2R10D","T3R1D","T3R4D","T3R5D","T3R7D","T3R9D","T3R10D","T4R1D","T4R4D","T4R5D","T4R7D","T4R9D","T4R10D","T5R1D","T5R4D","T5R5D","T5R7D","T5R9D","T5R10D")

physeq.whole <- subset_samples(physeq.a, SampleID %in% whole.samples)
physeq.live <- subset_samples(physeq.a, SampleID %in% live.samples)
physeq.dead <- subset_samples(physeq.a, SampleID %in% dead.samples)
physeq.lnd <- subset_samples(physeq.a, SampleID %in% live.dead.samples)

physeq.whole.percent <- transform_sample_counts(physeq.whole, function(x) 100 * x/sum(x))
physeq.w.percent.gyp <- subset_samples(physeq.whole.percent, Material == "Gypsum")
physeq.w.percent.mdf <- subset_samples(physeq.whole.percent, Material == "MDF")

physeq.l.percent <- transform_sample_counts(physeq.live, function(x) 100 * x/sum(x))
physeq.l.percent.gyp <- subset_samples(physeq.l.percent, Material == "Gypsum")
physeq.l.percent.mdf <- subset_samples(physeq.l.percent, Material == "MDF")

physeq.d.percent <- transform_sample_counts(physeq.dead, function(x) 100 * x/sum(x))  
physeq.d.percent.gyp <- subset_samples(physeq.d.percent, Material == "Gypsum")
physeq.d.percent.mdf <- subset_samples(physeq.d.percent, Material == "MDF")

# get counts
#sample_data(physeq.whole.percent)[,9]
count.whole <- as.data.frame(sample_data(physeq.whole.percent))$Count
count.live <- as.data.frame(sample_data(physeq.l.percent))$Count
count.dead <- as.data.frame(sample_data(physeq.d.percent))$Count

# function to convert relative abundance to quantitative abundance

rel_to_quan <- function(physeq, counts) {
    for (i in 1:nsamples(physeq)) {
        otu_table(physeq)[,i] = get_taxa(physeq, sample_names(physeq)[i]) * counts[i] /100
    }
    return(otu_table(physeq))
}

physeq.w.quan <- physeq.whole.percent
# replace relative otu table with a new quantitative one by using rel_to_quan function
otu_table(physeq.w.quan) <- rel_to_quan(physeq.w.quan, count.whole)

physeq.w.quan.gyp <- subset_samples(physeq.w.quan, Material == "Gypsum")
physeq.w.quan.mdf <- subset_samples(physeq.w.quan, Material == "MDF")

physeq.l.quan <- physeq.l.percent
# replace relative otu table with a new quantitative one by using rel_to_quan function
otu_table(physeq.l.quan) <- rel_to_quan(physeq.l.quan, count.live)

physeq.l.quan.gyp <- subset_samples(physeq.l.quan, Material == "Gypsum")
physeq.l.quan.mdf <- subset_samples(physeq.l.quan, Material == "MDF")

physeq.d.quan <- physeq.d.percent
# replace relative otu table with a new quantitative one by using rel_to_quan function
otu_table(physeq.d.quan) <- rel_to_quan(physeq.d.quan, count.dead)

physeq.d.quan.gyp <- subset_samples(physeq.d.quan, Material == "Gypsum")
physeq.d.quan.mdf <- subset_samples(physeq.d.quan, Material == "MDF")


