library(tidyverse)
library(grid)
library(gridExtra)

# some data preprocessing
input_directory <- "/Users/kevinxchan/Documents/UBC/YEAR4/MICB449/data/marker_gene_data/"
opt <- data.frame(input_directory, prefix, stringsAsFactors = FALSE)
opt$input_directory <- input_directory

file_names <- list.files(path = opt$input_directory)

for (file in file_names) {
  # if the merged dataset doesn't exist, create it
  if (!exists("acc_dat")) {
    acc_dat <- read.table(file.path(opt$input_directory, file), header = TRUE, sep = "\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("acc_dat")) {
    temp_dataset <- read.table(file.path(opt$input_directory, file), header=TRUE, sep="\t")
    acc_dat <- rbind(acc_dat, temp_dataset)
    rm(temp_dataset)
  }
}

acc_header <- c("MarkerGene", "DatasetID", "ReferenceID", "Software", "ReferenceTaxonomy", "MeanAlignmentPercentage", "TotalReadsAligned", "NumReadsAlignedOver80", "PercReadsAlignedOver80", "TaxDistance", "CumDistance")
names(acc_dat) <- acc_header
acc_dat$Software <- sub("graphmap", "GraphMap", acc_dat$Software)
head(acc_dat)

graphmap_dat <- acc_dat %>%
  filter(Software == "GraphMap" & TaxDistance >= 0)
minimap_dat <- acc_dat %>%
  filter(Software == "minimap2" & TaxDistance >= 0)

## overall stats including mean and median
summary(graphmap_dat$TaxDistance)
summary(minimap_dat$TaxDistance)
summary(graphmap_dat$CumDistance)
summary(minimap_dat$CumDistance)

# marker gene specific mean, median, IQR for taxa distance and contiguous taxa distance
markers <- c("pyrG", "recA", "ribosomal_L10P", "rpoB", "rps8")

for (marker in markers) {
  if (!exists("tax_distance_summary_graphmap")) {
    tax_distance_summary_graphmap <- graphmap_dat %>%
      filter(MarkerGene == marker) %>%
      summarise(mean_td = mean(TaxDistance),
                median_td = median(TaxDistance),
                iqr_td = IQR(TaxDistance),
                max_td = max(TaxDistance),
                min_td = min(TaxDistance),
                mean_ctd = mean(CumDistance),
                median_ctd = median(CumDistance),
                iqr_ctd = IQR(CumDistance),
                max_ctd = max(CumDistance),
                min_ctd = min(CumDistance)) %>%
      mutate(marker = marker)
    tax_distance_summary_graphmap <- select(tax_distance_summary_graphmap, marker, mean_td:min_ctd)
  } else {
    temp <- graphmap_dat %>%
      filter(MarkerGene == marker) %>%
      summarise(mean_td = mean(TaxDistance),
                median_td = median(TaxDistance),
                iqr_td = IQR(TaxDistance),
                max_td = max(TaxDistance),
                min_td = min(TaxDistance),
                mean_ctd = mean(CumDistance),
                median_ctd = median(CumDistance),
                iqr_ctd = IQR(CumDistance),
                max_ctd = max(CumDistance),
                min_ctd = min(CumDistance)) %>%
      mutate(marker = marker)
    temp <- select(temp, marker, mean_td:min_ctd)
    tax_distance_summary_graphmap <- rbind(tax_distance_summary_graphmap, temp)
  }
}

for (marker in markers) {
  if (!exists("tax_distance_summary_minimap")) {
    tax_distance_summary_minimap <- minimap_dat %>%
      filter(MarkerGene == marker) %>%
      summarise(mean_td = mean(TaxDistance),
                median_td = median(TaxDistance),
                iqr_td = IQR(TaxDistance),
                max_td = max(TaxDistance),
                min_td = min(TaxDistance),
                mean_ctd = mean(CumDistance),
                median_ctd = median(CumDistance),
                iqr_ctd = IQR(CumDistance),
                max_ctd = max(CumDistance),
                min_ctd = min(CumDistance)) %>%
      mutate(marker = marker)
    tax_distance_summary_minimap <- select(tax_distance_summary_minimap, marker, mean_td:min_ctd)
  } else {
    temp <- minimap_dat %>%
      filter(MarkerGene == marker) %>%
      summarise(mean_td = mean(TaxDistance),
                median_td = median(TaxDistance),
                iqr_td = IQR(TaxDistance),
                max_td = max(TaxDistance),
                min_td = min(TaxDistance),
                mean_ctd = mean(CumDistance),
                median_ctd = median(CumDistance),
                iqr_ctd = IQR(CumDistance),
                max_ctd = max(CumDistance),
                min_ctd = min(CumDistance)) %>%
      mutate(marker = marker)
    temp <- select(temp, marker, mean_td:min_ctd)
    tax_distance_summary_minimap <- rbind(tax_distance_summary_minimap, temp)
  }
}

fig3_dat <- group_by(acc_dat, Software, DatasetID, MarkerGene) %>% 
  summarise_at(vars(TotalReadsAligned), funs(sum)) %>% 
  rename(ReadPool = TotalReadsAligned)

harm_dist_dat <- fig3_dat %>% 
  merge(acc_dat, by = c("Software", "DatasetID", "MarkerGene")) %>% 
  mutate(Proportion = TotalReadsAligned/ReadPool) %>%
  mutate(CumTaxDistance = TaxDistance*Proportion) %>% 
  group_by(DatasetID, Software, MarkerGene) %>% 
  summarise_at(vars(CumTaxDistance), funs(sum))

for (marker in markers) {
  if (!exists("ctd_minimap")) {
    ctd_minimap <- harm_dist_dat %>%
      filter(MarkerGene == marker & Software == "minimap2") %>%
      ungroup() %>%
      summarise(mean_ctd = mean(CumTaxDistance),
                median_ctd = median(CumTaxDistance),
                iqr_ctd = IQR(CumTaxDistance),
                max_ctd = max(CumTaxDistance),
                min_ctd = min(CumTaxDistance)) %>%
      mutate(marker = marker)
    ctd_minimap <- select(ctd_minimap, marker, mean_ctd:min_ctd)
  } else {
    temp <- harm_dist_dat %>%
      filter(MarkerGene == marker & Software == "minimap2") %>%
      ungroup() %>%
      summarise(mean_ctd = mean(CumTaxDistance),
                median_ctd = median(CumTaxDistance),
                iqr_ctd = IQR(CumTaxDistance),
                max_ctd = max(CumTaxDistance),
                min_ctd = min(CumTaxDistance)) %>%
      mutate(marker = marker)
    temp <- select(temp, marker, mean_ctd:min_ctd)
    ctd_minimap <- rbind(ctd_minimap, temp)
  }
}

for (marker in markers) {
  if (!exists("ctd_graphmap")) {
    ctd_graphmap <- harm_dist_dat %>%
      filter(MarkerGene == marker & Software == "GraphMap") %>%
      ungroup() %>%
      summarise(mean_ctd = mean(CumTaxDistance),
                median_ctd = median(CumTaxDistance),
                iqr_ctd = IQR(CumTaxDistance),
                max_ctd = max(CumTaxDistance),
                min_ctd = min(CumTaxDistance)) %>%
      mutate(marker = marker)
    ctd_graphmap <- select(ctd_graphmap, marker, mean_ctd:min_ctd)
  } else {
    temp <- harm_dist_dat %>%
      filter(MarkerGene == marker & Software == "GraphMap") %>%
      ungroup() %>%
      summarise(mean_ctd = mean(CumTaxDistance),
                median_ctd = median(CumTaxDistance),
                iqr_ctd = IQR(CumTaxDistance),
                max_ctd = max(CumTaxDistance),
                min_ctd = min(CumTaxDistance)) %>%
      mutate(marker = marker)
    temp <- select(temp, marker, mean_ctd:min_ctd)
    ctd_graphmap <- rbind(ctd_graphmap, temp)
  }
}

write_tsv(tax_distance_summary_graphmap, path = "./td_graphmap.txt")
write_tsv(tax_distance_summary_minimap, path = "./td_minimap2.txt")
write_tsv(ctd_graphmap, path = "./cumtd_graphmap.txt")
write_tsv(ctd_minimap, path = "./cumtd_minimap2.txt")

