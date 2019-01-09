library(ggplot2)
library(RColorBrewer)
library(Rmisc)
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

option_list = list(make_option(c("-i", "--input_table"), type="character", default=NULL,
                               help="Tab-separated value file output by lca.py", metavar="character"),
                   make_option(c("-p", "--prefix"), type="character", default=".",
                               help="Prefix for the output files.", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_table)){
  print_help(opt_parser)
  stop("--input_table must be provided!", call.=FALSE)
}

##
# For debugging
##
prefix <- "."
input_table_minimap2 <- "/Users/kevinxchan/Documents/UBC/YEAR4/MICB449/data/lca_scores/recA/recA_minimap2.txt"
input_table_graphmap <- "/Users/kevinxchan/Documents/UBC/YEAR4/MICB449/data/lca_scores/recA/recA_graphmap.txt"
opt <- data.frame(input_table_minimap2, input_table_graphmap, prefix, stringsAsFactors = FALSE)

spec_out <- paste(opt$prefix, "lca.png", sep='_')
sens_out <- paste(opt$prefix, "Classification_sensitivity.png", sep='_')
pdist_out <- paste(opt$prefix, "Classification_WeightedDistance.png", sep='_')

minimap2 <- read.table(opt$input_table_minimap2,
                      sep="\t", header=TRUE)
graphmap <- read.table(opt$input_table_graphmap,
                       sep="\t", header=TRUE)
acc_dat <- rbind(minimap2, graphmap)

lca_header <- c("Dataset", "Reference ID", "Software", "Parameters", "Reference Taxonomy", "Total", "# Aligned Reads >= 80%", "% Reads Aligned >= 80%", "Distance")
names(acc_dat) <- lca_header
acc_dat <- select(acc_dat, -Parameters)
acc_dat$Software <- sub("^(PN:)(.*)", "\\2", acc_dat$Software)

##
# Figure 1: Rank-exclusion specificity in relation to optimal placement distance
##
sum_software <- acc_dat %>%
  group_by(Dataset, Software) %>%
  summarise(TotalBySoftware = sum(Total, na.rm = T))
tmp <- merge(acc_dat, sum_software, by = c("Dataset", "Software"))
fig_1 <- tmp %>%
  select("Software", "Total", "Distance", "TotalBySoftware") %>%
  group_by(Dataset, Software, Distance) %>%
  summarise(Total = sum(Total)) %>%
  filter(Distance >= 0)
to_plot <- merge(to_plot, sum_software, by = "Software")
to_plot <- mutate(to_plot, Proportion = Total / TotalBySoftware * 100)

filter(to_plot, Distance >= 4) %>% 
  filter(Total > 0)
pd <- position_dodge(width = 0.75)

ggplot(to_plot, aes(x=Distance, y=Proportion)) +
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
