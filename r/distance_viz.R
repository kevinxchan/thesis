library(ggplot2)
library(RColorBrewer)
# library(Rmisc)
library(dplyr)
library(optparse)

#########
# TODO: #
#########
# 1. after lca.py finished, remake the data frame and rename the columns appropriately
# 2. figure 1 = real taxonomic distance vs proportion
# 3. figure 2 = cumulative taxonomic distance vs marker gene. when comparing multiple markers, need to write
#    more code to parse through a directory full of these dataframes and merge them
# 4. instead of dividing by all totals, group by query dataset too and then divide by this total count.

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
prefix <- "viz"
input_directory <- "/Users/kevinxchan/Documents/UBC/YEAR4/MICB449/data/marker_gene_data/"
opt <- data.frame(input_directory, prefix, stringsAsFactors = FALSE)

spec_out <- file.path(opt$prefix, "Proportion_TaxDistance.png")
pdist_out <- file.path(opt$prefix, "WeightedDistance_MarkerGene.png")

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

acc_header <- c("MarkerGene", "DatasetID", "ReferenceID", "Software", "ReferenceTaxonomy", "TotalReadsAligned", "NumReadsAlignedOver80", "PercReadsAlignedOver80", "TaxDistance", "CumDistance")
names(acc_dat) <- acc_header
acc_dat$Software <- sub("^(PN:)(.*)", "\\2", acc_dat$Software)

##
# Figure 1: Rank-exclusion specificity in relation to optimal placement distance
##
sum_software <- acc_dat %>%
  group_by(MarkerGene, Software) %>%
  summarise(TotalBySoftware = sum(TotalReadsAligned, na.rm = T))

fig_1 <- acc_dat %>%
  select("MarkerGene", "Software", "TotalReadsAligned", "TaxDistance") %>%
  group_by(MarkerGene, Software, TaxDistance) %>%
  summarise(Total = sum(TotalReadsAligned)) %>%
  filter(TaxDistance >= 0)

tmp <- merge(sum_software, fig_1, by = c("MarkerGene", "Software"))
fig_1 <- mutate(tmp, Proportion = Total / TotalBySoftware * 100)

filter(fig_1, TaxDistance >= 4) %>% 
  filter(Total > 0)
pd <- position_dodge(width = 0.75)

ggplot(fig_1, aes(x=TaxDistance, y=Proportion, fill=MarkerGene)) +
  geom_bar(stat="identity", position=pd, colour="black", width = 0.75) +
  facet_wrap(~Software) +
  scale_fill_brewer(palette = "PuOr") +
  scale_y_continuous(breaks=seq(0,100,10)) +
  xlab("Distance from Optimal Rank") +
  ylab("Percentage of Queries") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

ggsave(filename = spec_out, width = 10, height = 6, dpi = 400)

##
# Figure 2: The taxonomic-hazard-weighted classification scores
##
to_plot <- mutate(to_plot, Marker = "recA")

harm_dist_dat <- acc_dat %>% 
  filter(Proportion > 0) %>% 
  merge(build_params, by.x = "Marker", by.y = "name") %>% 
  mutate(Marker = reorder(Marker, ref_sequences)) %>% 
  mutate(Rank = reorder(Rank, Position)) %>%
  mutate(PlaceDist = Distance*Proportion) %>% 
  group_by(Trial, Software, Marker, Rank, ref_sequences) %>% 
  summarise_at(vars(PlaceDist), funs(sum))

# TODO: break down by dataset
# sum_software <- acc_dat %>%
#   group_by(MarkerGene, DatasetID, Software) %>%
#   summarise(TotalBySoftware = sum(TotalReadsAligned, na.rm = T))
# 
# head(acc_dat)
# 
# fig_1 <- acc_dat %>%
#   select("MarkerGene", "DatasetID", "Software", "TotalReadsAligned", "TaxDistance") %>%
#   group_by(MarkerGene, DatasetID, Software, TaxDistance) %>%
#   summarise(Total = sum(TotalReadsAligned)) %>%
#   filter(TaxDistance >= 0)
# 
# tmp <- fig_1 %>%
#   group_by(MarkerGene, DatasetID, Software) %>%
#   summarise(TotalBySoftwareDataset = sum(Total))
# 
# tmp <- merge(tmp, fig_1, by = c("MarkerGene", "DatasetID", "Software"))
# fig_1 <- mutate(tmp, Proportion = Total / TotalBySoftwareDataset * 100)

ggplot(to_plot, aes(x=Marker, y=Distance)) +
  geom_point(aes(fill=Marker), colour="black",
             shape=21, size=3, alpha=2/3) +
  scale_fill_brewer(palette = "PuOr") +
  facet_wrap(~Software) +
  ylab("Cumulative Taxonomic Distance") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = pdist_out, width = 8, height = 5, dpi = 400)

harm_dist_dat %>%
  group_by(Trial) %>% 
  summarise_at(vars(PlaceDist), funs(mean, median))

##
# Welch Two Sample t-test between the two trials
##
cat("Trials evaluated: ")
unlist(unique(harm_dist_dat$Trial)[1:2])
t.test(filter(harm_dist_dat, Trial == unique(harm_dist_dat$Trial)[1])$PlaceDist,
       filter(harm_dist_dat, Trial == unique(harm_dist_dat$Trial)[2])$PlaceDist)
