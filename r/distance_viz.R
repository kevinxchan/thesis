library(ggplot2)
library(RColorBrewer)
# library(Rmisc)
library(dplyr)
library(optparse)

option_list = list(make_option(c("-i", "--input_directory"), type="character", default=NULL,
                               help="A directory containing tab-separated value files output by lca.py", metavar="character"),
                   make_option(c("-p", "--prefix"), type="character", default=".",
                               help="Prefix for the output files. [current working directory]", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_directory)){
  print_help(opt_parser)
  stop("--input_directory must be provided!", call.=FALSE)
}

##
# For debugging
##
# prefix <- "viz"
# input_directory <- "/Users/kevinxchan/Documents/UBC/YEAR4/MICB449/data/marker_gene_data/with_lengths/"
# opt <- data.frame(input_directory, prefix, stringsAsFactors = FALSE)

spec_out <- file.path(opt$prefix, "Proportion_TaxDistance.png")
pdist_out <- file.path(opt$prefix, "WeightedDistance_MarkerGene.png")
cdist_out <- file.path(opt$prefix, "WeightedDistance_Dataset.png")

# start building up acc_dat by scanning through files in directory
file_names <- list.files(path = opt$input_directory)

for (file in file_names) {
  
  # if the merged dataset doesn't exist, create it
  if (!exists("acc_dat")) {
    foo <- file.path(opt$input_directory, file)
    write(foo, stdout())
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

##
# Figure 1: Rank-exclusion specificity in relation to optimal placement distance
##
sum_software <- acc_dat %>%
  group_by(MarkerGene, Software) %>%
  summarise(TotalBySoftware = sum(TotalReadsAligned, na.rm = T))

fig_1 <- acc_dat %>%
  filter(MeanAlignmentPercentage >= 20) %>%
  filter(TaxDistance >= 0) %>%
  select("MarkerGene", "Software", "MeanAlignmentPercentage", "TotalReadsAligned", "TaxDistance") %>%
  group_by(MarkerGene, Software, TaxDistance) %>%
  summarise(Total = sum(TotalReadsAligned))

tmp <- merge(sum_software, fig_1, by = c("MarkerGene", "Software"))
fig_1 <- mutate(tmp, Proportion = Total / TotalBySoftware * 100)

filter(fig_1, TaxDistance >= 4) %>% 
  filter(Total > 0)
pd <- position_dodge(width = 0.75)

fig_1 %>%
  mutate(MarkerGene = ifelse(MarkerGene == "ribosomal_L10P", "RPL10", as.character(MarkerGene))) %>%
  
  ggplot(aes(x=TaxDistance, y=Proportion, fill=MarkerGene)) +
  geom_bar(stat="identity", position=pd, colour="black", width = 0.75) +
  facet_wrap(~Software) +
  scale_fill_brewer(palette = "PuOr") +
  scale_y_continuous(breaks=seq(0,100,10)) +
  xlab("Distance from Optimal Rank") +
  ylab("Percentage of Queries") +
  labs(fill = "Marker Gene") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(filename = spec_out, width = 10, height = 6, dpi = 400)

##
# Figure 2: The taxonomic-hazard-weighted classification scores
##

# Determine the number of reads that were aligned by each tool for each dataset and marker gene
fig3_dat <- group_by(acc_dat, Software, DatasetID, MarkerGene) %>% 
  summarise_at(vars(TotalReadsAligned), funs(sum)) %>% 
  rename(ReadPool = TotalReadsAligned)

harm_dist_dat <- fig3_dat %>% 
  merge(acc_dat, by = c("Software", "DatasetID", "MarkerGene")) %>% 
  mutate(Proportion = TotalReadsAligned/ReadPool) %>%
  mutate(CumTaxDistance = TaxDistance*Proportion) %>% 
  group_by(DatasetID, Software, MarkerGene) %>% 
  summarise_at(vars(CumTaxDistance), funs(sum))

ggplot(harm_dist_dat, aes(x=DatasetID, y=CumTaxDistance)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "PuOr") +
  facet_wrap(~Software) +
  ylab("Cumulative Taxonomic Distance") +
  xlab("Dataset ID") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = cdist_out, width = 8, height = 5, dpi = 400)

harm_dist_dat %>%
  mutate(MarkerGene = ifelse(MarkerGene == "ribosomal_L10P", "RPL10", as.character(MarkerGene))) %>%
  
  ggplot(aes(x=MarkerGene, y=CumTaxDistance)) +
  geom_jitter(aes(fill=DatasetID), position=position_jitter(0.1), colour="black",
              shape=21, size=3, alpha=2/3) +
  scale_fill_brewer(palette = "PuOr") +
  facet_wrap(~Software) +
  ylab("Cumulative Taxonomic Distance") +
  xlab("Marker Gene") +
  labs(fill = "Dataset ID") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = pdist_out, width = 8, height = 5, dpi = 400)

stat <- harm_dist_dat %>%
  filter(!is.na(CumTaxDistance)) %>% 
  group_by(Software) %>% 
  summarise_at(vars(CumTaxDistance), funs(mean, median, var))

stat

##
# Welch Two Sample t-test between the two trials
##
cat("Trials evaluated: ")
unlist(unique(harm_dist_dat$Trial)[1:2])
t.test(filter(harm_dist_dat, Trial == unique(harm_dist_dat$Trial)[1])$PlaceDist,
       filter(harm_dist_dat, Trial == unique(harm_dist_dat$Trial)[2])$PlaceDist)
